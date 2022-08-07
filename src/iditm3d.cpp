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
    Vec            X, Xn, Xn1;
    Field          ***xx;
    DM             da;
    AppCtx         params;               /* user-defined work context */

    PetscMPIInt    nprocs, rank; //, i;
    int            nprocx, nprocy, nprocz=1;
    PetscInt       xs, ys, zs, xm, ym, zm;
    char           buff[150];
    string         fname;
    time_t         start_t, end_t;
    struct         tm *now;
    ofstream       logfstr;

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

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --------
    Create distributed array (DMDA) to manage parallel grid and vectors PETSC
    decides number of processes for a1 and a2 but specified as 1 for a3
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    DMDACreate3d(MPI_COMM_WORLD,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_PERIODIC,
           DMDA_STENCIL_STAR,a1,a2,a3,PETSC_DECIDE,PETSC_DECIDE,1,a4,2,NULL,NULL,NULL,&da);

    DMSetFromOptions(da);
    DMSetUp(da);

    DMDAGetInfo(da, NULL, NULL, NULL, NULL, &nprocx, &nprocy, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    if (!rank) output_param(argv, nprocx, nprocy, nprocz, &params);

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
    VecDuplicate(X, &Xn);
    VecDuplicate(X, &Xn1);
    VecDuplicate(X, &params.U);
    VecDuplicate(X, &params.V);
    VecDuplicate(X, &params.W);
    VecDuplicate(X, &params.Z);

    DMDAGetCorners(da, &xs ,&ys, &zs, &xm, &ym, &zm);

    /* set up grids and related geometric parameters. Input initial solution */
    if (initialize(da, Xn1, &params) < 0) exit(-1);

    /* using halftime step to get first step solution at tn + dt/2 */
    forward_scheme(da, Xn, Xn1, &params);

    if (!rank) cout <<endl<< "Start time advancing ..." <<endl;

    params.ntot = params.ntot + params.npre;
    /*************** start time advancing *************************************/
    for (params.ndt = params.npre+1; params.ndt < params.ntot+1; params.ndt++) {
        params.sec += dt;
        if (params.sec >= 86400) update_timedate(&params);

        //advance the first half time step for explicit part
        imex_leap_frog(da, X, Xn, Xn1, &params);

        // output in parallel to a hdf5 file at chosen time steps
        if (params.ndt % params.nout ==0 || params.ndt==params.ntot) {
            DMDAVecGetArray(da, X, &xx);
            output_solution(da, xx, &params);
            DMDAVecRestoreArray(da, X, &xx);
        }

        if (params.ndt % params.lognum ==0 || params.ndt==params.ntot) {
            end_t=time(NULL);
            end_time=(double)end_t;
            MPI_Reduce(&end_time, &allProcessesTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

            PetscPrintf(MPI_COMM_WORLD,"Step %D: Run time used %g (s), Phyiscal time simulated %12.5f (s)\n",
               params.ndt, allProcessesTime-start_time, (double)(params.ndt-params.npre)*dt);
        }
    }
    /*************** done time advancing *************************************/

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Free work space.
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    array_deallocate(xm, ym, zm);

    VecDestroy(&params.U);
    VecDestroy(&params.V);
    VecDestroy(&params.W);
    VecDestroy(&params.Z);
    VecDestroy(&Xn1);
    VecDestroy(&Xn);
    VecDestroy(&X);

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
