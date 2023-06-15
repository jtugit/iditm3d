static char help[] = "Inductive-Dynamic-Ionosphere-Thermosphere (IDIT) Model. \n";

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

#include "param.h"
#include "param_def.h"
#include "funcdef.h"

#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc,char **argv)
{
    Vec            X;
    Field          ***xx, ***xn;
    DM             da;
    AppCtx         params;               /* user-defined work context */

    PetscMPIInt    nprocs, rank; //, i;
    PetscInt       nprocx, nprocy, nprocz=1;
    PetscInt       xs, ys, zs, xm, ym, zm;
    char           buff[150];
    string         fname;
    time_t         start_t, end_t;
    struct         tm *now;
    ofstream       logfstr;
    int            aa;

    PetscInitialize(&argc,&argv,(char*)0,help);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    double allProcessesTime, start_time, end_time;
    start_t = time(NULL);
    start_time = (double)start_t;
    MPI_Reduce(&start_time, &allProcessesTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    start_time = allProcessesTime;

    if (!rank) cout << "Inductive-Dynamic-Ionosphere-Thermosphere (IDIT) Model!" << endl << endl;

    //get working directory
    if(getcwd(buff, sizeof(buff)) == NULL) {
        perror("getcwd() error");
        return -1;
    }
    size_t pos = string(buff).find("/src");
    params.workdir = string(buff).substr(0, pos);

    /* Input run setting parameters */
    if (input_param(&params) < 0) exit(-1);
    double startUT = params.sec;

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --------
    Create distributed array (DMDA) to manage parallel grid and vectors PETSC
    decides number of processes for a1 and a2 but specified as 1 for a3
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    DMDACreate3d(MPI_COMM_WORLD,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_PERIODIC,
           DMDA_STENCIL_STAR,a1,a2,a3,PETSC_DECIDE,PETSC_DECIDE,1,a4,2,NULL,NULL,NULL,&da);

    DMSetFromOptions(da);
    DMSetUp(da);

    DMDAGetInfo(da, NULL, NULL, NULL, NULL, &nprocx, &nprocy, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    if (!rank) output_param(argv, (int)nprocx, (int)nprocy, (int)nprocz, &params);

    if(a1/nprocx < 5 || a2/nprocy < 5) {
        cout<<"Too few grids on r or theta coordinate a process owns!"
            <<" Specify different number of MPI processes"<<endl;
            return -1;
    }
    if(a1 % nprocx !=0) {
        cout<<"Number of processes on r = "<<nprocx
            <<" must be a factor of grid number on r-axis = "<<a1<<endl;
            return -1;
    }
    if(a2 % nprocy !=0) {
        cout<<"Number of processes on theta = "<<nprocy
            <<" must be a factor of grid number on theta-axis = "<<a2<<endl;
            return -1;
    }

    /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Extract global vectors from DMDA;
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    DMCreateGlobalVector(da, &X);
    VecDuplicate(X, &params.U);
    VecDuplicate(X, &params.Xn);

    DMDAGetCorners(da, &xs ,&ys, &zs, &xm, &ym, &zm);
    aa=zm*ym*xm;

    /* set up grids and related geometric parameters. Input initial solution */
    if (initialize(da, X, &params) < 0) exit(-1);

    SNES snes;
    SNESCreate(MPI_COMM_WORLD, &snes);
    SNESSetType(snes, SNESNEWTONTR);
    SNESSetFromOptions(snes);
    SNESSetDM(snes, da);
    SNESSetFunction(snes, NULL, formfunctions, &params);

    KSP ksp;
    SNESGetKSP(snes, &ksp);
    KSPSetType(ksp, KSPFGMRES);
    KSPSetFromOptions(ksp);

    PC pc;
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCJACOBI);
    PCSetFromOptions(pc);

    Mat A;
    DMSetMatrixPreallocateOnly(da, PETSC_FALSE);
    DMSetMatType(da, MATMPIAIJ);
    DMDASetBlockFills(da, dfill, ofill);
    DMCreateMatrix(da, &A);

    SNESSetJacobian(snes, A, A, jacobian, &params);

    if (!rank) cout <<endl<< "Start time advancing ..." <<endl;

    params.ntot = params.ntot + params.npre;
    /*************** start time advancing *************************************/
    for (params.ndt = params.npre+1; params.ndt < params.ntot+1; params.ndt++) {
        params.sec += dt*t0;
        if (params.sec >= 86400) update_timedate(&params);

        top_bc_vel(&params, ys, ym, zs, zm);
        parameters(da, X, &params);

        SNESSetSolution(snes, X);   //set initial guess of the solution
        SNESSolve(snes, PETSC_NULL, X);   //iterative solver to find the solution

        DMDAVecGetArray(da, X, &xx);
        check_positivity(da, xx);

/*-----------------------------------------------------------------------------------------
 *----- smooth O+, H+, He+ ion velocities and temperatures 
 *----- multidomensional Shapiro filter is used to conduct smoothing (Falissard, JCP 2013)
 -----------------------------------------------------------------------------------------*/
        //smooth_multi_dim_X(da, X);

        //copy solution to previous time step vector
        DMDAVecGetArray(da, params.Xn, &xn);
        copy(&xx[zs][ys][xs], &xx[zs][ys][xs]+aa, &xn[zs][ys][xs]);

        // output in parallel to a hdf5 file at chosen time steps
        if (params.ndt % params.nout ==0 || params.ndt==params.ntot) output_solution(da, xx, &params);
        DMDAVecRestoreArray(da, X, &xx);
        DMDAVecRestoreArray(da, params.Xn, &xn);

        if (params.ndt % params.lognum ==0 || params.ndt==params.ntot) {
            end_t=time(NULL);
            end_time=(double)end_t;
            MPI_Reduce(&end_time, &allProcessesTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

            PetscPrintf(MPI_COMM_WORLD,"Step %D: Run time used %g (s), Phyiscal time simulated %12.5f (s)\n",
               params.ndt, allProcessesTime-start_time, params.sec-startUT);
        }
    }
    /*************** done time advancing *************************************/

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Free work space.
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    array_deallocate(xm, ym, zm);

    VecDestroy(&params.U);
    VecDestroy(&params.Xn);
    VecDestroy(&X);

    SNESDestroy(&snes);
    MatDestroy(&A);
    DMDestroy(&da);

    if (!rank) {
        end_t=time(NULL);
        now=localtime(&end_t);

        cout <<"Program ended: "
             <<now->tm_year+1900<<"-"<<now->tm_mon+1<<"-"<<now->tm_mday
             <<" " <<now->tm_hour<<":"<<now->tm_min<< ":"<<now->tm_sec<<endl;

        logfstr.open(params.outpdir + "/iditm3d.log", fstream::app);

        logfstr <<"Program ended: "
                <<now->tm_year+1900<<"-"<<now->tm_mon+1<<"-"<<now->tm_mday
                <<" " <<now->tm_hour<<":"<<now->tm_min<< ":"<<now->tm_sec<<endl;
        logfstr.close();
    }

    PetscFinalize();
    PetscFunctionReturn(0);
}
