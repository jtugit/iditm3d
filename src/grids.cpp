#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
using namespace std;

#include "param.h"

#define tag1 1

int grids(DM da, AppCtx *params)
{
    int        i, j, k;
    PetscInt   xs, ys, zs, xm, ym, zm;
    string     fname, ch1, tempstr;
    PetscMPIInt rank, nproc;

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);

    DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);

    dr=(params->ru-params->rb)/(double)Nr;
    dth = 180.0 / (double)Nth;
    dph= 360.0 / (double)a3;

/* --------- start run from the scratch: first set up grid points */
    /* r_{i} */
    for (i = 0; i < a1; i++) rr[i] = dr/2.0+params->rb + dr*(double)i; 

    /* theta_{j} in deg*/
    theta[0] = dth/2.0;
    for (j = 1; j < Nth; j++) theta[j] = theta[0] + dth*double(j);
    theta[Nth]=theta[Nthm];  //ghost grid on the other side of southern pole

    /* phi_{k} in deg */
    phi[0] = dph/2.0;
    for (k = 1; k <= Np; k++) phi[k] = phi[0] + dph*double(k);

    /* convert to radians */
    dth=dth*rad;
    dph=dph*rad;

    for (j = 0; j < a2; j++) theta[j]=theta[j]*rad;
    for (k = 0; k < a3; k++) phi[k]=phi[k]*rad;

/***************************************************************************
* rh[0]    rr[0]    rh[1]    rr[1]    rh[2]             rh[Nr-1] rr[Nr-1]  rh[Nr]
*   |--------x--------|--------X--------|--------X---...---|--------X--------|--------X
* -1/2               1/2               3/2               Nr-3/2            Nr-1/2
* im=0               im=1              im=2             im=Nr-1            im=Nr
*           i=0               i=1               i=2 ............. i=Nr-1            i=Nr (ghost)
* rh[i] is the left interface of the cell i centered at rr[i] (geometric center rc[i])
* rh[i+1] is the right interface of the cell i
****************************************************************************/
    rh[0] = rr[0] - dr/2.0;
    for (i = 1; i < a1; i++) {
        rh[i]=0.5*(rr[i]+rr[i-1]);
    }
    rh[a1]=rh[Nr]+dr; //extra rh for calculating ghost rC[Nr], zh[Nr]

    //rC[i] is geometric center of the cell i=0, 1, ..., Nrm (Nrm=Nr-1)
    for (i = 0; i < a1; i++) {
        rC[i] = 0.75*(pow(rh[i+1], 4.0) - pow(rh[i], 4.0))/(pow(rh[i+1], 3.0) - pow(rh[i], 3.0));
        rfavg[i]= 2.0/3.0*(pow(rh[i+1], 3.0) - pow(rh[i], 3.0))/(rh[i+1]*rh[i+1] - rh[i]*rh[i]);
        zh[i] = (rC[i]-Re)*1.0e-3;

        rh_d3[i] = 1.0/3.0*(rh[i+1]*rh[i+1]*rh[i+1] - rh[i]*rh[i]*rh[i]);
        rh_d2[i] = 0.5*(rh[i+1]*rh[i+1] - rh[i]*rh[i]);

        rh2[i]=rh[i]*rh[i];
        rfavg_dth[i]=rfavg[i]*dth;
    }
    //geometric center of the ghost cells
    double rhm=rh[0]-dr, rhup=rh[a1]+dr;
    rCm1 = 0.75*(pow(rh[0], 4.0)-pow(rhm, 4.0))/(pow(rh[0], 3.0)-pow(rhm, 3.0));
    rfavgm1 = 2.0/3.0*(pow(rh[0], 3.0)-pow(rhm, 3.0))/(rh[0]*rh[0]-rhm*rhm);
    rC[a1]=0.75*(pow(rhup,4.0)-pow(rh[a1],4.0))/(pow(rhup,3.0)-pow(rh[a1],3.0));    

/***************************************************************************
*thetah[0]=0 theta[0]  theath[1]   theta[1]  thetah[2]...theta[Nthm] thetah[Nth]=pi
*      |----------x-----------|----------X----------| ........ X----------|---------X
*                                                                                   | ghost grid
*                                                                theta[Nth]=theta[Nthm] kc=(k+a3/3) mod a3
* thetah[j] is the left interface of the grid theta[j] (geometric venter thetac[j])
* thetah[j+1] is the right interface of the same grid
****************************************************************************/
    double *costh_h = new double[a2+1];
    double *costh_hd = new double[a2];

    thetah[0] = 0.0;
    sinth_h[0]=sin(thetah[0]); sinth[0]=sin(theta[0]); costh_h[0]=cos(thetah[0]);
    for (j = 1; j <= Nthm; j++) {
        thetah[j]=0.5*(theta[j]+theta[j-1]);
        sinth_h[j]=sin(thetah[j]);
        sinth[j]=sin(theta[j]);

        costh_h[j]=cos(thetah[j]);
    }
    thetah[Nth]=pi; sinth_h[Nth]=0.0; costh_h[Nth]=-1.0;

    for (j = 0; j <= Nthm; j++) {
        costh_hd[j] = costh_h[j] - costh_h[j+1];
        thetaC[j] = (thetah[j]*costh_h[j]-thetah[j+1]*costh_h[j+1]+sinth_h[j+1]-sinth_h[j])/costh_hd[j];
        cotth[j] = cos(theta[j])/sin(theta[j]);
    }
    thetaC[Nth]=thetaC[Nthm]; //geometric center of the ghost cell on the other side of the pole

/***************************************************************************
*phih[0]=0 phi[0]  phih[1]   phi[1]  phih[2]...phi[Np-1] phih[Np]=2pi
*   |--------X--------|--------X--------| ........ X--------|---------X
*                                                                     | ghost grid
*                                                              phi[Nth]=phi[0]
* thetah[j] is the left interface of the grid theta[j] (geometric venter thetac[j])
* thetah[j+1] is the right interface of the same grid.
* Note phi[0] and phi[Np] is the same point
****************************************************************************/
    phih[0]=0.0;
    for (k = 1; k <= Np; k++) phih[k] = 0.5*(phi[k]+phi[k-1]);
    phih[Np+1] = 2.0*pi;

    /*----------------------------------------------------------------------------------*/
    for (i = xs; i< xs+xm; i++) {
        /* normalized gravitational acceleration at radial distance r */
        gr[i-xs]=ge*pow(Re/(rC[i]), 2.0);
    }

    for (j = 0; j < Nth; j++) {
        for (i = 0; i < Nr; i++) {
            rCsinC[j][i]=rC[i]*sin(thetaC[j]);
            rfavg_costh[j][i]=rfavg[i]*costh_hd[j];
            rfavg_costh_dth_dph[j][i]=rfavg_costh[j][i]*dph/dth;
            rfavg_sinth_dph[j][i]=rfavg[i]*sinth_h[j]*dph;
            rh_costh[j][i]=rh[i]*costh_hd[j];
            rh_costh_dth_dph[j][i]=rh_costh[j][i]*dph/dth;
            dAtheta_dV[j][i]=(sinth_h[j+1]-sinth_h[j])/(rfavg[i]*cotth[j]);
        }
    }

/*------------ output grids information --------------------------------------*/
    if (!rank) {
        fname = params->outpdir;
        fname = fname + "/grids.dat";
        fstream fstr(fname, fstream::out);

        fstr<<"# "<<setw(5)<<a1<<setw(5)<<a2<<setw(5)<<a3<<endl;
        fstr<<"  i        rC               r              rh            zh (km)"<<endl;
        for (i = 0; i < a1; i++) {
            fstr << setw(4) << i;
            fstr <<scientific<<setw(16)<<setprecision(8)<<rC[i];
            fstr <<scientific<<setw(16)<<setprecision(8)<<rr[i];
            fstr <<scientific<<setw(16)<<setprecision(8)<<rh[i];
            fstr <<scientific<<setw(16)<<setprecision(8)<<zh[i] << endl;
        }
        fstr<<"  j    thetaC (deg)     theta (deg)    thetah (deg)"<<endl;
        for (j = 0; j < a2; j++) {
            fstr << setw(4) << j;
            fstr <<scientific<<setw(16)<<setprecision(8)<<thetaC[j]*deg;
            fstr <<scientific<<setw(16)<<setprecision(8)<<theta[j]*deg;
            fstr <<scientific<<setw(16)<<setprecision(8)<<thetah[j]*deg << endl;
        }

        /* mean center = geometric center for phi when mesh along phi is uniform */
        fstr<<"  k      phi (deg)      phih (deg)"<<endl;
        for (k = 0; k < a3; k++) {
            fstr << setw(4) << k;
            fstr <<scientific<<setw(16)<<setprecision(8)<<phi[k]*deg;
            fstr <<scientific<<setw(16)<<setprecision(8)<<phih[k]*deg << endl;
        }

        fstr.close();
    }

    delete[] costh_h;
    delete[] costh_hd;

    return 0;
}
