/***********************************************************************
!    This function evaluates convection velocities for given locations 
!    (specified by mlat and mlon).  The program also provides 
!    electrical potential. Weimer's model is used to calculate electric
!    field potential at h=110 km.
!
!    Author: Jiannan Tu
!    Date:   March 18, 2021
!***********************************************************************/
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
using namespace std;

#include "param.h"

#define tag 1

extern "C" void convection_veloc(int *, int *, int *, int *, int *, int *, 
        float *, float *, float *, float *, float *, float *, float *, float *, 
        float *, int *, float*, float *, float *, float *, float *, float *);

/* Top boundary BC for ion velocity evaluated at rC[Nr], thetaC[j], phi[k] */
void top_bc_vel(AppCtx *params, PetscInt ys, PetscInt ym, PetscInt zs, PetscInt zm)
{
    int     yr, mn, day, hr, imin, sec;
    int     j, k, UseAL, zk, yj;

    float  Bximf, Byimf, Bzimf, SWDen, vgsex, vgsey, vgsez, ALindex, gMLT;
    float  Vr, Vth, Vph, epot, mlat, mlong;

    yr=params->iyr;
    mn=params->mon;
    day=params->idate;
    hr=int(params->sec/3600.0);
    imin=int((params->sec-hr*3600)/60.0);
    sec=(int)(params->sec)-hr*3600-imin*60;

    Bximf=(float)(params->Bximf);
    Byimf=(float)(params->Byimf);
    Bzimf=(float)(params->Bzimf);
    SWDen=(float)(params->SWDen);
    vgsex=(float)(params->vgsex);
    vgsey=(float)(params->vgsey);
    vgsez=(float)(params->vgsez);
    UseAL=params->UseAL;
    ALindex=(float)(params->ALindex);

    for(k=zs; k<zs+zm; k++) {
        zk=k-zs;

        mlong=(float)(phi[k]/rad);

        for(j=ys; j<ys+ym; j++) {
            yj=j-ys;

            mlat=(float)(90.0-theta[j]/rad);

            if (abs(mlat) < 50.0) {
                vt[zk][yj]=-1.0e6;
                vp[zk][yj]=-1.0e6;
            }
            else {
                convection_veloc(&yr, &mn, &day, &hr, &imin, &sec, &mlat, &mlong,
                  &Bximf, &Byimf, &Bzimf, &SWDen, &vgsex, &vgsey, &vgsez, &UseAL,
                  &ALindex, &gMLT, &Vr, &Vth, &Vph, &epot);

                vt[zk][yj]=(double)Vth*params->rurb3/v0;
                vp[zk][yj]=(double)Vph*params->rurb3/v0;
            }
        }
    }
}

void print_top_bc_vel(AppCtx *params, PetscInt xsm, PetscInt ys, PetscInt ym, 
        PetscInt zs, PetscInt zm)
{
    int     yr, mn, day, hr, imin, sec;
    int     j, k, UseAL;

    float  Bximf, Byimf, Bzimf, SWDen, vgsex, vgsey, vgsez, ALindex, gMLT;
    float  Vr, Vth, Vph, epot, mlat, mlong;

    PetscMPIInt        rank, nproc, file_free=0;
    //MPI_Status status;

    fstream  fstr;

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);

    yr=params->iyr;
    mn=params->mon;
    day=params->idate;
    hr=int(params->sec/3600.0);
    imin=int((params->sec-hr*3600)/60.0);
    sec=(int)(params->sec)-hr*3600-imin*60;

    Bximf=(float)(params->Bximf);
    Byimf=(float)(params->Byimf);
    Bzimf=(float)(params->Bzimf);
    SWDen=(float)(params->SWDen);
    vgsex=(float)(params->vgsex);
    vgsey=(float)(params->vgsey);
    vgsez=(float)(params->vgsez);
    UseAL=params->UseAL;
    ALindex=(float)(params->ALindex);

    file_free=0;

    /* the first process does printing first */
    if (rank==0 || xsm == a1) {
        file_free=1;

        fstr.open("conveloc.dat", fstream::out);

        fstr<<"Date: "<<yr<<" "<<mn<<" "<<day<<" UT: "<<hr<<" "<<imin<<" "<<sec<<endl;
        fstr<<"IMF-Bx: "<<Bximf<<" IMF-By: "<<Byimf<<" IMF-Bz: "<<Bzimf<<endl;
        fstr<<"SW density: "<<SWDen<<" SW GSE-Vx: "<<vgsex<<" SW GSE-Vy: "<<vgsey
            <<" SW GSE-Vz: "<<vgsez<<endl;
        if (UseAL) fstr<<"Use AL: Yes"<<endl;
        else fstr<<"Use AL: No"<<endl;
        fstr<<"AL index: "<<ALindex<<endl;
        fstr<<setw(11)<<"MLAT"<<setw(10)<<"MLT"<<setw(12)<<"MLONG"<<setw(15)<<"poten (kV)"
            <<setw(8)<<"Vr"<<setw(12)<<"Vth"<<setw(12)<<"Vph"<<endl;
    }
    /*else if (rank > 0 && xsm == a1) {
        MPI_Recv(&file_free,1,MPI_INT,rank-1,tag,MPI_COMM_WORLD,&status);

        fstr.open("conveloc.dat",fstream::app);
        fstr.seekp(0,fstream::end);
    }*/

    if(file_free) {
        for(k=zs; k<zs+zm; k++) {
            mlong=(float)(phi[k]/rad);
            for(j=ys; j<ys+ym; j++) {
                mlat=(float)(90.0-theta[j]/rad);

                if (abs(mlat) < 50.0) continue;

                convection_veloc(&yr, &mn, &day, &hr, &imin, &sec, &mlat, &mlong,
                  &Bximf, &Byimf, &Bzimf, &SWDen, &vgsex, &vgsey, &vgsez, &UseAL,
                  &ALindex, &gMLT, &Vr, &Vth, &Vph, &epot);

                vt[k-zs][j-ys]=(double)Vth*params->rurb3;
                vp[k-zs][j-ys]=(double)Vph*params->rurb3;

                fstr<<scientific<<setw(12)<<setprecision(4)<<mlat;
                fstr<<scientific<<setw(12)<<setprecision(4)<<gMLT;
                fstr<<scientific<<setw(12)<<setprecision(4)<<mlong;
                fstr<<scientific<<setw(12)<<setprecision(4)<<epot;
                fstr<<scientific<<setw(12)<<setprecision(4)<<Vr*params->rurb3;
                fstr<<scientific<<setw(12)<<setprecision(4)<<vt[k-zs][j-ys];
                fstr<<scientific<<setw(12)<<setprecision(4)<<vp[k-zs][j-ys]<<endl;
            }
        }

        fstr.close();
    }

    if(rank != nproc-1 && xsm == a1) 
        MPI_Send(&file_free,1,MPI_INT,rank+1,tag,MPI_COMM_WORLD);

}
