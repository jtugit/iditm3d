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
    int        i, j, k, xs, ys, zs, xm, ym, zm, xi, yj, zk;
    uint32_t   yj_t;
    string     fname, ch1, tempstr;
    PetscMPIInt rank, nproc;

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);

    DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);

    //boundary radial distances normalized
    params->rb=(1.0e3*params->rb+Re)/r0;
    params->ru=(1.0e3*params->ru+Re)/r0;

    params->rurb3=pow(params->rb/params->ru, 3.0);

    dr=(params->ru-params->rb)/(double)Nr;
    dth = 180.0 / (double)Nthm;
    dph= 360.0 / (double)a3;

/* --------- start run from the scratch: first set up grid points */
    /* r_{i} */
    for (i = 0; i < a1; i++) rr[i] = params->rb + dr*(double)i;

    /* theta_{j} in deg*/
    theta[0]=dth/2.0;   //dummy element
    theta[1]=dth/2.0;
    for (j = 2; j < Nth; j++) theta[j] = theta[j-1] + dth;
    theta[Nth]=theta[Nthm];  //dummy element

    /* phi_{k} in deg */
    phi[0] = dph/2.0;
    for (k = 1; k < a3; k++) phi[k] = phi[0] + dph*double(k);

    /* convert to radians */
    dth=dth*rad;
    dph=dph*rad;

    for (j = 0; j < a2; j++) theta[j]=theta[j]*rad;
    for (k = 0; k < a3; k++) phi[k]=phi[k]*rad;

/*----------- transform matrices -------------------------------------------*/
    double sinphi, cosphi;
    vector<double> sintheta, costheta;

    for (j = ys; j < ys+ym; j++) {
        sintheta.push_back(sin(theta[j])); costheta.push_back(cos(theta[j]));
    }

    for (k = zs; k < zs+zm; k++) {
        zk=k-zs;

        sinphi=sin(phi[k]); cosphi=cos(phi[k]);

        for (j = ys; j < ys+ym; j++) {
            yj=j-ys; yj_t=(uint32_t)yj;

            J11[zk][yj] = sintheta[yj_t]*cosphi;
            J12[zk][yj] = sintheta[yj_t]*sinphi;

            Jiv11[zk][yj] = sintheta[yj_t]*cosphi;
            Jiv21[zk][yj] = sintheta[yj_t]*sinphi;

            K11[zk][yj] = sintheta[yj_t]*cosphi;
            K12[zk][yj] = costheta[yj_t]*cosphi;
            K21[zk][yj] = sintheta[yj_t]*sinphi;
            K22[zk][yj] = costheta[yj_t]*sinphi;

            for (i = xs; i < xs+xm; i++) {
                xi=i-xs;

                J21[zk][yj][xi] = costheta[yj_t]*cosphi/rr[i];
                J22[zk][yj][xi] = costheta[yj_t]*sinphi/rr[i];

                J31[zk][yj][xi] = -sinphi/(rr[i]*sintheta[yj_t]);
                J32[zk][yj][xi] =  cosphi/(rr[i]*sintheta[yj_t]);

                Jiv12[zk][yj][xi] = rr[i]*costheta[yj_t]*cosphi;
                Jiv13[zk][yj][xi] =-rr[i]*sintheta[yj_t]*sinphi;

                Jiv22[zk][yj][xi] = rr[i]*costheta[yj_t]*sinphi;
                Jiv23[zk][yj][xi] = rr[i]*sintheta[yj_t]*cosphi;
            }
        }

        K13[zk]=-sinphi; K23[zk]= cosphi;
    }

    for (j = ys; j < ys+ym; j++) {
        yj=j-ys; yj_t=(uint32_t)yj;

        J13[yj] = costheta[yj_t]; Jiv31[yj] = costheta[yj_t];

        K31[yj] = costheta[yj_t]; K32[yj] =-sintheta[yj_t];

        for (i = xs; i <xs+xm; i++) {
            xi=i-xs;

            J23[yj][xi] = sintheta[yj_t]/rr[i];
            Jiv32[yj][xi] = -rr[i]*sintheta[yj_t];
        }
    }
/*----------- end calculation of transform matrices -----------------------------------*/

/***************************************************************************************
* rh[0]    rr[0]    rh[1]    rr[1]    rh[2]             rh[Nr-1] rr[Nr-1]  rh[Nr]
*   |--------x--------|--------X--------|--------X---...---|--------X--------|--------X
* -1/2               1/2               3/2               Nr-3/2            Nr-1/2
* im=0               im=1              im=2             im=Nr-1            im=Nr
*           i=0               i=1               i=2 ............. i=Nr-1            i=Nr
* rh[i] is the left interface of the cell i centered at rr[i] (geometric center rc[i])
* rh[i+1] is the right interface of the cell i
***************************************************************************************/
    rh[0] = rr[0] - dr/2.0;
    zh[0]=(rr[0]*r0-Re)*1.0e-3;
    for (i = 1; i < a1; i++) {
        rh[i]=0.5*(rr[i]+rr[i-1]);
        zh[i] = (rr[i]*r0-Re)*1.0e-3;
    }
    rh[a1]=rh[Nr]+dr; //extra rh for calculating rC[Nr], zh[Nr]

/**************************************************************************************
*thetah[1]=0 theta[1]  theath[2]   theta[2]  thetah[3]...theta[Nthm] thetah[Nth]=pi
*      |----------x-----------|----------X----------| ........ X----------|
* thetah[j] is the left interface of the grid theta[j] (geometric center thetac[j])
* thetah[j+1] is the right interface of the same grid. Internal grids are in [1:Nth-1]
**************************************************************************************/
    thetah[0]=0.0;   //dummy element
    thetah[1]=0.0;
    for (j = 2; j < Nth; j++) {
        thetah[j]=0.5*(theta[j]+theta[j-1]);
    }
    thetah[Nth]=pi;

    for (j = ys; j < ys+ym; j++) {
        yj=j-ys; yj_t=(uint32_t)yj;

        for (i = xs; i < xs+xm; i++) {
            r2sintheta[yj][i-xs]=rr[i]*rr[i]*sintheta[yj_t];
            cot_div_r[yj][i-xs]=cos(theta[j])/(rr[i]*sintheta[yj_t]);
            rsin[yj][i-xs]=rr[i]*sintheta[yj_t];
        }
    }

/***************************************************************************
*phih[0]=0 phi[0]  phih[1]   phi[1]  phih[2]...phi[Np-1] phih[Np]=2pi
*   |--------X--------|--------X--------| ........ X--------|---------X
*                                                                     | ghost grid
*                                                              phi[Nth]=phi[0]
* thetah[j] is the left interface of the grid theta[j] (geometric venter thetac[j])
* thetah[j+1] is the right interface of the same grid.
* Note phih[0] and phih[a3] is the same point
****************************************************************************/
    phih[0]=0.0;
    for (k = 1; k < a3; k++) phih[k] = 0.5*(phi[k]+phi[k-1]);
    phih[a3] = 2.0*pi;

    /*----------------------------------------------------------------------------------*/
    for (i = xs; i< xs+xm; i++) {
        /* normalized gravitational acceleration at radial distance r */
        gr[i-xs]=gen*pow(Re/(rr[i]*r0), 2.0);
    }

/*------------ output grids information --------------------------------------*/
    if (!rank) {
        fname = params->outpdir;
        fname = fname + "/grids.dat";
        fstream fstr(fname, fstream::out);

        fstr<<"# "<<setw(5)<<a1<<setw(5)<<a2<<setw(5)<<a3<<endl;
        fstr<<"  i        r              rh            zh (km)"<<endl;
        for (i = 0; i < a1; i++) {
            fstr << setw(4) << i;
            fstr <<scientific<<setw(16)<<setprecision(8)<<rr[i]*r0;
            fstr <<scientific<<setw(16)<<setprecision(8)<<rh[i]*r0;
            fstr <<scientific<<setw(16)<<setprecision(8)<<zh[i] << endl;
        }
        fstr<<"  j    theta (deg)    thetah (deg)"<<endl;
        for (j = 0; j < a2; j++) {
            fstr << setw(4) << j;
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

    return 0;
}
