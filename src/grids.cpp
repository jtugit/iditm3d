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
    PetscInt   i, j, k;
    PetscInt   xs, ys, zs, xm, ym, zm;
    string     fname, ch1, tempstr;
    PetscMPIInt rank, nproc;

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);

    DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);

    dr=(params->ru-params->rb)/(double)Nr;
    dth = 180.0 / (double)(Nth);
    dph= 360.0 / (double)a3;

/* --------- start run from the scratch: first set up grid points */
    /* r_{i} */
    for (i = 0; i < a1; i++) rr[i] = dr/2.0+params->rb + dr*(double)i; 

    /* theta_{j} in deg*/
    theta[0]=dth/2.0;   //dummy element
    theta[1]=dth/2.0;
    for (j = 2; j < Nth; j++) theta[j] = theta[1] + dth*double(j);
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
    for (k = zs; k < zs+zm; k++) {
        for (j = ys; j < ys+ym; j++) {
            Jmat.J11.push_back(sin(theta[j])*cos(phi[k]));
            Jmat.J12.push_back(sin(theta[j])*sin(phi[k]));
        }
    }

    for (j = ys; j < ys+ym; j++) Jmat.J13.push_back(cos(theta[j]));

    for (k = zs; k < zs+zm; k++) {
        for (j = ys; j < ys+ym; j++) {
            for (i = xs; i < xs+xm; i++) {
                Jmat.J21.push_back(cos(theta[j])*cos(phi[k])/rr[i]);
                Jmat.J22.push_back(cos(theta[j])*sin(phi[k])/rr[i]);

                Jmat.J31.push_back(-sin(phi[k])/(rr[i]*sin(theta[j])));
                Jmat.J32.push_back( cos(phi[k])/(rr[i]*sin(theta[j])));
            }
        }
    }

    for (j = ys; j < ys+ym; j++) {
        for (i = xs; i < xs+xm; i++) {
            Jmat.J23.push_back(-sin(theta[j])/rr[i]);
        }
    }

    for (k = zs; k < zs+zm; k++) {
        for (j = ys; j < ys+ym; j++) {
            Jinv.Jiv11.push_back(sin(theta[j]*cos(phi[k])));
            Jinv.Jiv21.push_back(sin(theta[j]*sin(phi[k])));
        }
    }

    for (k = zs; k < zs+zm; k++) {
        for (j = ys; j < ys+ym; j++) {
            for (i = xs; i < xs+xm; i++) {
                Jinv.Jiv12.push_back( rr[i]*cos(theta[j])*cos(phi[k]));
                Jinv.Jiv13.push_back(-rr[i]*sin(theta[j])*sin(phi[k]));

                Jinv.Jiv22.push_back( rr[i]*cos(theta[j])*sin(phi[k]));
                Jinv.Jiv23.push_back( rr[i]*sin(theta[j])*cos(phi[k]));
            }
        }
    }

    for (j = ys; j < ys+ym; j++) Jinv.Jiv31.push_back(cos(theta[j]));

    for (j = ys; j < ys+ym; j++) {
        for (i = xs; i < xs+xm; i++) {
            Jinv.Jiv32.push_back(-rr[i]*sin(theta[j]));
        }
    }

    for (k = zs; k < zs+zm; k++) {
        for (j = ys; j < ys+ym; j++) {
            Kmat.K11.push_back(sin(theta[j])*cos(phi[k]));
            Kmat.K12.push_back(cos(theta[j])*cos(phi[k]));
            Kmat.K21.push_back(sin(theta[j])*sin(phi[k]));
            Kmat.K22.push_back(cos(theta[j])*sin(phi[k]));
        }
    }

    for (k = zs; k < zs+zm; k++) {
        Kmat.K13.push_back(-sin(phi[k]));
        Kmat.K23.push_back( cos(phi[k]));
    }

    for (j = ys; j < ys+ym; j++) {
        Kmat.K31.push_back( cos(theta[j]));
        Kmat.K32.push_back(-sin(theta[j]));
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
    for (i = 1; i < a1; i++) {
        rh[i]=0.5*(rr[i]+rr[i-1]);
        zh[i] = (rr[i]-Re)*1.0e-3;
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
        if (j < 1) j = 1;
        if (j > Nthm) j = Nthm;

        for (i = xs; i < xs+xm; i++) {
            r2sintheta[j-ys][i-xs]=rr[i]*rr[i]*sin(theta[j]);
            cot_div_r[j-ys][i-xs]=cos(theta[j])/(rr[i]*sin(theta[j]));
            rsin[j-ym][i-xs]=rr[i]*sin(theta[j])*dph;
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
            fstr <<scientific<<setw(16)<<setprecision(8)<<rr[i];
            fstr <<scientific<<setw(16)<<setprecision(8)<<rh[i];
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
