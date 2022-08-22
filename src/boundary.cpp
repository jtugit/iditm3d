/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <cmath>
using namespace std;

#include "funcdef.h"
#include "check_positivity.h"

void boundary_bc(Field ***xx, int xs, int xm, int ys, int ym, int zs, int zm, AppCtx *params)
{
    int     i, j, k, kc, s, yj, zk, Nrm2=Nr-2;
    double  y1, y2, sgn, rho, rrb, rrt, ne, Nn;

    rrb=(rC[0]-rC[1])/(rC[2]-rC[1]);
    rrt=(rC[Nr]-rC[Nrm])/(rC[Nrm]-rC[Nrm2]);

    double rrbh=(rh[0]-rh[1])/(rh[2]-rh[1]);
    double rrbi=(rr[0]-rr[1])/(rr[2]-rr[1]);

    double rrth=(rh[Nr]-rh[Nrm])/(rh[Nrm]-rh[Nrm2]);
    double rrti=(rr[Nr]-rr[Nrm])/(rr[Nrm]-rr[Nrm2]);

   if (xs+xm == a1) top_bc_vel(params, ys, ym, zs, zm);

    for ( k = zs; k < zs+zm; k++) {
        zk = k-zs;

        for (j =  ys; j < ys+ym; j++) {
            yj = j-ys;

/*---------------------------------------------------------------*/
/*------- lower boundary conditions -----------------------------*/
/*---------------------------------------------------------------*/
            if (xs == 0) {
                //density of ion and neutral species
                ne = 0.0; Nn = 0.0;
                for (s=0; s< sl; s++) {
                    y1=log(xx[k][j][1].fx[s]);
                    y2=log(xx[k][j][2].fx[s]);
                    xx[k][j][0].fx[s]= exp(y1 + rrb*(y2-y1));
                    ne += xx[k][j][0].fx[s];

                    y1=log(xx[k][j][1].fx[12+s]);
                    y2=log(xx[k][j][2].fx[12+s]);
                    xx[k][j][0].fx[12+s]= exp(y1 + rrb*(y2-y1));
                    Nn += xx[k][j][0].fx[12+s];
                }

                for (s=7; s <= 9; s++) {
                    //(rho ui_r), (rho ui_theta)
                    y1= xx[k][j][1].fx[s];
                    y2= xx[k][j][2].fx[s];
                    xx[k][j][0].fx[s]= y1 + rrb*(y2-y1);

                    //rhon un_r), (rhon un_theta)
                    y1= xx[k][j][1].fx[12+s];
                    y2= xx[k][j][2].fx[12+s];
                    xx[k][j][0].fx[12+s]= y1 + rrb*(y2-y1);
                }

                // pi = ni*kb*Ti
                y1= log(xx[k][j][1].fx[10]);
                y2= log(xx[k][j][2].fx[10]);
                xx[k][j][0].fx[10] = exp(y1 + rrb*(y2-y1));

                // pn = Nn*kb*Tn
                y1= log(xx[k][j][1].fx[22]);
                y2= log(xx[k][j][2].fx[22]);
                xx[k][j][0].fx[22] = exp(y1 + rrb*(y2-y1));

                if (xx[k][j][0].fx[11]/(kb*ne) < xx[k][j][0].fx[22]/(kb*Nn))
                    xx[k][j][0].fx[11] = ne/Nn*xx[k][j][0].fx[22];

                // pe = ni*kb*Te
                xx[k][j][0].fx[11]=xx[k][j][0].fx[10];

                //Br_{i=-1/2,j,k},  Btheta_{i=0,j-1/2,k} and Bphi_{i=0,j,k-1/2}
                y1=xx[k][j][1].fx[23];
                y2=xx[k][j][2].fx[23];
                xx[k][j][0].fx[23] = y1+rrbh*(y2-y1);

                for (s = 24; s < 26; s++) {
                    y1=xx[k][j][1].fx[s];
                    y2=xx[k][j][2].fx[s];
                    xx[k][j][0].fx[s] = y1+rrbi*(y2-y1);
                }
            }


/*---------------------------------------------------------------*/
/*------- upper boundary conditions -----------------------------*/
/*---------------------------------------------------------------*/
            if (xs+xm == a1) {
                //density of ion and neutral species
                for (s=0; s<sl; s++) {
                    y1=log(xx[k][j][Nrm].fx[s]);
                    y2=log(xx[k][j][Nrm2].fx[s]);
                    xx[k][j][Nr].fx[s] = exp(y1+rrt*(y1-y2));

                    y1=log(xx[k][j][Nrm].fx[12+s]);
                    y2=log(xx[k][j][Nrm2].fx[12+s]);
                    xx[k][j][Nr].fx[12+s] = exp(y1+rrt*(y1-y2));
                }

                //(rho ui_r)
                y1=xx[k][j][Nrm].fx[7];
                y2=xx[k][j][Nrm2].fx[7];
                xx[k][j][Nr].fx[7] = y1+rrt*(y1-y2);

                //Br_{Nr-1/2,j,k}
                y1=xx[k][j][Nrm].fx[23];
                y2=xx[k][j][Nrm2].fx[23];
                xx[k][j][Nr].fx[23] = y1+rrth*(y1-y2);

                if (vt[zk][yj] > -1.e5) {
                    rho=0.0;
                    for (s=0; s<sl; s++) rho += xx[k][j][Nr].fx[s]*ms[s];

                    // (rho ui_theta)                      
                    xx[k][j][Nr].fx[8] = rho*vt[zk][yj];

                    // Btheta_{Nr,j-1/2,k}
                    if (j<Nth/2) sgn=-1.0; else sgn=1.0;
                    xx[k][j][Nr].fx[24] = sgn*vt[zk][yj]*sqrt(rho*mu0);
                }
                else {
                    y1=xx[k][j][Nrm].fx[8];
                    y2=xx[k][j][Nrm2].fx[8];
                    xx[k][j][Nr].fx[8] = y1+rrt*(y1-y2);

                    y1=xx[k][j][Nrm].fx[24];
                    y2=xx[k][j][Nrm2].fx[24];
                    xx[k][j][Nr].fx[24] = y1+rrti*(y1-y2);
                }

                if (vp[zk][yj] > -1.e5) {
                    // (rho ui_phi)
                    rho=0.0;
                    for (s=0; s<sl; s++) rho += xx[k][j][Nr].fx[s]*ms[s];

                    xx[k][j][Nr].fx[9] = rho*vp[zk][yj];

                    // Bphi_{Nr,j,k-1/2}
                    if (j<Nth/2) sgn=-1.0; else sgn=1.0;
                    xx[k][j][Nr].fx[25] = sgn*vp[zk][yj]*sqrt(rho*mu0);
                }
                else {
                    y1=xx[k][j][Nrm].fx[9];
                    y2=xx[k][j][Nrm2].fx[9];
                    xx[k][j][Nr].fx[9] = y1+rrt*(y1-y2);

                    y1=xx[k][j][Nrm].fx[25];
                    y2=xx[k][j][Nrm2].fx[25];
                    xx[k][j][Nr].fx[25] = y1+rrti*(y1-y2);
                }

                // pi = ni*kb*Ti and pe = ne*kb*Te
                for (s = 10; s <= 11; s++) {
                    y1= log(xx[k][j][Nrm].fx[s]);
                    y2= log(xx[k][j][Nrm2].fx[s]);
                    xx[k][j][Nr].fx[s]= exp(y1 + rrt*(y1-y2));
                }

                //rhon un_r), (rhon un_theta), (rhon un_phi)
                for (s = 19; s <= 21; s++) {
                    xx[k][j][Nr].fx[s]= xx[k][j][Nrm].fx[s];
                    /*y1= xx[k][j][Nrm].fx[s];
                    y2= xx[k][j][Nrm2].fx[s];
                    xx[k][j][Nr].fx[s]= y1 + rrt*(y1-y2);*/
                }

                y1= log(xx[k][j][Nrm].fx[22]);
                y2= log(xx[k][j][Nrm2].fx[22]);
                xx[k][j][Nr].fx[22]= exp(y1 + rrt*(y1-y2));
            }
        }
    }

/*-----------------------------------------------------------------------------------
 * set values on the grids across the southern pole (that is j=Nth at the other side)
 *---------------------------------------------------------------------------------*/
    if (ys+ym == a2) {
        for ( k = zs; k < zs+zm; k++) {
            kc = (k+a3/2) % a3;

            for (i = xs; i < xs+xm; i++) {
                for (s = 0; s < sl; s++) {
                    xx[k][Nth][i].fx[s] = xx[kc][Nthm][i].fx[s];
                    xx[k][Nth][i].fx[12+s] = xx[kc][Nthm][i].fx[12+s];
                }

                xx[k][Nth][i].fx[7] = xx[kc][Nthm][i].fx[7];
                xx[k][Nth][i].fx[8] =-xx[kc][Nthm][i].fx[8];
                for (s = 9; s < 12; s++) xx[k][Nth][i].fx[s] = xx[kc][Nthm][i].fx[s];

                xx[k][Nth][i].fx[19] = xx[kc][Nthm][i].fx[19];
                xx[k][Nth][i].fx[20] =-xx[kc][Nthm][i].fx[20];
                xx[k][Nth][i].fx[21] = xx[kc][Nthm][i].fx[21];
                xx[k][Nth][i].fx[22] = xx[kc][Nthm][i].fx[22];

                xx[k][Nth][i].fx[23] = xx[kc][Nthm][i].fx[23];
                xx[k][Nth][i].fx[25] = xx[kc][Nthm][i].fx[25];
            }
        }
    }
}

void boundary_V_T(Field ***xx, int xs, int xm, int ys, int ym, int zs, int zm)
{
    int     i, j, k, kc, s, yj, zk, Nrm2=Nr-2;
    double  y1, y2, rho, rrb, rrt;

    rrb=(rC[0]-rC[1])/(rC[2]-rC[1]);
    rrt=(rC[Nr]-rC[Nrm])/(rC[Nrm]-rC[Nrm2]);

    for ( k = zs; k < zs+zm; k++) {
        zk = k-zs;

        for (j =  ys; j < ys+ym; j++) {
            yj = j-ys;

/*---------------------------------------------------------------*/
/*------- lower boundary conditions -----------------------------*/
/*---------------------------------------------------------------*/
            if (xs == 0) {
                for (s=7; s <= 9; s++) {
                    //(rho ui_r), (rho ui_theta)
                    y1= xx[k][j][1].fx[s];
                    y2= xx[k][j][2].fx[s];
                    xx[k][j][0].fx[s]= y1 + rrb*(y2-y1);

                    //rhon un_r), (rhon un_theta)
                    y1= xx[k][j][1].fx[12+s];
                    y2= xx[k][j][2].fx[12+s];
                    xx[k][j][0].fx[12+s]= y1 + rrb*(y2-y1);
                }

                // pi = ni*kb*Ti and pe = ne*kb*Te
                for (s=10; s <= 11; s++) {
                    y1= log(xx[k][j][1].fx[s]);
                    y2= log(xx[k][j][2].fx[s]);
                    xx[k][j][0].fx[s] = exp(y1 + rrb*(y2-y1));
                }

                // pn = Nn*kb*Tn
                y1= log(xx[k][j][1].fx[22]);
                y2= log(xx[k][j][2].fx[22]);
                xx[k][j][0].fx[22] = exp(y1 + rrb*(y2-y1));
            }

/*---------------------------------------------------------------*/
/*------- upper boundary conditions -----------------------------*/
/*---------------------------------------------------------------*/
            if (xs+xm == a1) {
                //(rho ui_r)
                y1=xx[k][j][Nrm].fx[7];
                y2=xx[k][j][Nrm2].fx[7];
                xx[k][j][Nr].fx[7] = y1+rrt*(y1-y2);

                if (vt[zk][yj] > -1.e5) {
                    rho=0.0;
                    for (s=0; s<sl; s++) rho += xx[k][j][Nr].fx[s]*ms[s];

                    // (rho ui_theta)                      
                    xx[k][j][Nr].fx[8] = rho*vt[zk][yj];
                }
                else {
                    y1=xx[k][j][Nrm].fx[8];
                    y2=xx[k][j][Nrm2].fx[8];
                    xx[k][j][Nr].fx[8] = y1+rrt*(y1-y2);
                }

                if (vp[zk][yj] > -1.e5) {
                    // (rho ui_phi)
                    rho=0.0;
                    for (s=0; s<sl; s++) rho += xx[k][j][Nr].fx[s]*ms[s];

                    xx[k][j][Nr].fx[9] = rho*vp[zk][yj];
                }
                else {
                    y1=xx[k][j][Nrm].fx[9];
                    y2=xx[k][j][Nrm2].fx[9];
                    xx[k][j][Nr].fx[9] = y1+rrt*(y1-y2);
                }

                // pi = ni*kb*Ti and pe = ne*kb*Te
                for (s = 10; s <= 11; s++) {
                    y1= log(xx[k][j][Nrm].fx[s]);
                    y2= log(xx[k][j][Nrm2].fx[s]);
                    xx[k][j][Nr].fx[s]= exp(y1 + rrt*(y1-y2));
                }

                //rhon un_r), (rhon un_theta), (rhon un_phi)
                for (s = 19; s <= 21; s++) {
                    y1= xx[k][j][Nrm].fx[s];
                    y2= xx[k][j][Nrm2].fx[s];
                    xx[k][j][Nr].fx[s]= y1 + rrt*(y1-y2);
                }

                y1= log(xx[k][j][Nrm].fx[22]);
                y2= log(xx[k][j][Nrm2].fx[22]);
                xx[k][j][Nr].fx[22]= exp(y1 + rrt*(y1-y2));
            }
        }
    }

/*-----------------------------------------------------------------------------------
 * set values on the grids across the southern pole (that is j=Nth at the other side)
 *---------------------------------------------------------------------------------*/
    if (ys+ym == a2) {
        for ( k = zs; k < zs+zm; k++) {
            kc = (k+a3/2) % a3;

            for (i = xs; i < xs+xm; i++) {
                xx[k][Nth][i].fx[7] = xx[kc][Nthm][i].fx[7];
                xx[k][Nth][i].fx[8] =-xx[kc][Nthm][i].fx[8];
                for (s = 9; s < 12; s++) xx[k][Nth][i].fx[s] = xx[kc][Nthm][i].fx[s];

                xx[k][Nth][i].fx[19] = xx[kc][Nthm][i].fx[19];
                xx[k][Nth][i].fx[20] =-xx[kc][Nthm][i].fx[20];
                xx[k][Nth][i].fx[21] = xx[kc][Nthm][i].fx[21];
                xx[k][Nth][i].fx[22] = xx[kc][Nthm][i].fx[22];
            }
        }
    }
}
