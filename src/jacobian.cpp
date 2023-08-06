#include <cmath>
#include <iostream>

using namespace std;

#include <petscts.h>
#include <petscdm.h>
#include <petscdmda.h>

#include "param.h"
#include "funcdef.h"
#include "ele_cooling_rate.h"
#include "functions.h"

int jacobian(SNES snes, Vec X, Mat Jac, Mat Jpre, void *ctx)
{
    DM         da;
    Vec        localX, localXn, localU;
    Field      ***xx, ***xn, ***uu;
    AppCtx     *params = (AppCtx*)ctx;
    int        i, j, k, ir, nv, xs, xm, ys, ym, zs, zm, xi, yj, zk, s;
    MatStencil row, *col;
    PetscErrorCode ierr=0;
    double     *vals;
    const int  nzer=a4;

    double     temp;
    const double dt_quart=0.25*dt, dx = 1.0e-8;

    double Qeuv, Qphoto, Ps[14], Ls[14];

    SNESGetDM(snes, &da);

    DMDAGetCorners(da, &xs ,&ys, &zs, &xm, &ym, &zm);

    DMGetLocalVector(da, &localX);
    DMGlobalToLocalBegin(da, X, INSERT_VALUES,localX);
    DMGlobalToLocalEnd(da, X, INSERT_VALUES, localX);
    DMDAVecGetArray(da, localX, &xx);

    DMGetLocalVector(da, &localXn);
    DMGlobalToLocalBegin(da,params->Xn, INSERT_VALUES, localXn);
    DMGlobalToLocalEnd(da,params->Xn, INSERT_VALUES, localXn);
    DMDAVecGetArray(da, localXn, &xn);

    DMGetLocalVector(da, &localU);
    DMGlobalToLocalBegin(da,params->U,INSERT_VALUES,localU);
    DMGlobalToLocalEnd(da,params->U,INSERT_VALUES,localU);
    DMDAVecGetArray(da, localU, &uu);

    vals=new double[nzer];
    col=new MatStencil[nzer];

    for (k = zs; k < zs+zm; k++) {
      row.k=k; zk=(k-zs);

      for (j = ys; j < ys+ym; j++) {
        row.j=j;

        if (j == 0 || j == Nth) {
            for (i = xs; i < xs+xm; i++) {
                row.i = i;

                for (ir = 0; ir < a4; ir++) {
                    row.c=ir; nv=0;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ir; vals[nv] = 1.0;
                    nv++;

                    MatSetValuesStencil(Jac,1,&row,nv,col,vals,INSERT_VALUES);
                }
            }
            continue;
        }

        yj=(j-ys);

        for (i = xs; i < xs+xm; i++) {
            row.i=i;

            if (i == 0 || i == Nr) {
                for (ir = 0; ir < a4; ir++) {
                    row.c=ir; nv=0;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ir; vals[nv] = 1.0;
                    nv++;

                    MatSetValuesStencil(Jac,1,&row,nv,col,vals,INSERT_VALUES);
                }
                continue;
            }

            xi=(i-xs);

            //production and loss rates
            prod_loss_rates(xn, uu, i, j, k, zk, yj, xi, Qeuv, Qphoto, Ps, Ls);

            for (ir = 0; ir < a4; ir++) {
                row.c=ir; nv=0;

                if (ir < sl) {
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ir;
                    vals[nv]=1.0+dt*Ls[ir]; nv++;
                }

                if ((ir >= 7 && ir < 20) || (ir > 26 && ir < 31) || ir > 33) 
                    temp=functions(xx, xn, uu, i, j, k, xi, yj, zk, ir);

                if (ir == 7) {
                    for (int ss = 7; ss <= 15; ss++) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ss;

                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]+dx;
                        vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]-dx;

                        nv++;
                    }

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=27;

                    xx[k][j][i].fx[27]=xx[k][j][i].fx[27]+dx;
                    vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                    xx[k][j][i].fx[27]=xx[k][j][i].fx[27]-dx;

                    nv++;
                }
                else if (ir == 8) {
                    for (int ss = 7; ss <= 15; ss++) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ss;

                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]+dx;
                        vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]-dx;

                        nv++;
                    }

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=28;
                    xx[k][j][i].fx[28]=xx[k][j][i].fx[28]+dx;
                    vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                    xx[k][j][i].fx[28]=xx[k][j][i].fx[28]-dx;

                    nv++;
                }
                else if (ir == 9) {
                    for (int ss = 7; ss <= 15; ss++) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ss;

                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]+dx;
                        vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]-dx;

                        nv++;
                    }

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=29;

                    xx[k][j][i].fx[29]=xx[k][j][i].fx[29]+dx;
                    vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                    xx[k][j][i].fx[29]=xx[k][j][i].fx[29]-dx;

                    nv++;
                }
                else if (ir == 10) {
                    for (int ss = 7; ss <= 15; ss++) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ss;

                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]+dx;
                        vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]-dx;

                        nv++;
                    }

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=27;

                    xx[k][j][i].fx[27]=xx[k][j][i].fx[27]+dx;
                    vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                    xx[k][j][i].fx[27]=xx[k][j][i].fx[27]-dx;

                    nv++;
                }
                else if (ir == 11) {
                    for (int ss = 7; ss <= 15; ss++) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ss;

                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]+dx;
                        vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]-dx;

                        nv++;
                    }

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=28;

                    xx[k][j][i].fx[28]=xx[k][j][i].fx[28]+dx;
                    vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                    xx[k][j][i].fx[28]=xx[k][j][i].fx[28]-dx;

                    nv++;
                }
                else if (ir == 12) {
                    for (int ss = 7; ss <= 15; ss++) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ss;

                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]+dx;
                        vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]-dx;

                        nv++;
                    }

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=29;

                    xx[k][j][i].fx[29]=xx[k][j][i].fx[29]+dx;
                    vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                    xx[k][j][i].fx[29]=xx[k][j][i].fx[29]-dx;

                    nv++;
                }
                else if (ir == 13) {
                    for (int ss = 7; ss <= 15; ss++) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ss;

                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]+dx;
                        vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]-dx;

                        nv++;
                    }

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=27;

                    xx[k][j][i].fx[27]=xx[k][j][i].fx[27]+dx;
                    vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                    xx[k][j][i].fx[27]=xx[k][j][i].fx[27]-dx;

                    nv++;
                }
                else if (ir == 14) {
                    for (int ss = 7; ss <= 15; ss++) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ss;

                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]+dx;
                        vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]-dx;

                        nv++;
                    }

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=28;

                    xx[k][j][i].fx[28]=xx[k][j][i].fx[28]+dx;
                    vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                    xx[k][j][i].fx[28]=xx[k][j][i].fx[28]-dx;

                    nv++;
                }
                else if (ir == 15) {
                    for (int ss = 7; ss <= 15; ss++) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ss;

                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]+dx;
                        vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]-dx;

                        nv++;
                    }

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=29;

                    xx[k][j][i].fx[29]=xx[k][j][i].fx[29]+dx;
                    vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                    xx[k][j][i].fx[29]=xx[k][j][i].fx[29]-dx;

                    nv++;
                }
                else if (ir == 16) {  //O+ temperature equation
                    for (int ss = 7; ss <= 19; ss++) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ss;

                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]+dx;
                        vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]-dx;

                        nv++;
                    }

                    for (int ss = 27; ss <= 30; ss++) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ss;

                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]+dx;
                        vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]-dx;

                        nv++;
                    }
                }
                else if (ir == 17) {  //H+ temperature equation
                    for (int ss = 7; ss <= 19; ss++) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ss;

                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]+dx;
                        vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]-dx;

                        nv++;
                    }

                    for (int ss = 27; ss <= 30; ss++) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ss;

                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]+dx;
                        vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]-dx;

                        nv++;
                    }
                }
                else if (ir == 18) {  //He+ temperature equation
                    for (int ss = 7; ss <= 19; ss++) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ss;

                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]+dx;
                        vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]-dx;

                        nv++;
                    }

                    for (int ss = 27; ss <= 30; ss++) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ss;

                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]+dx;
                        vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]-dx;

                        nv++;
                    }
                }
                else if (ir == 19) {  //electron temperature equation
                    for (int ss = 7; ss <= 19; ss++) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ss;

                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]+dx;
                        vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]-dx;

                        nv++;
                    }

                    for (int ss = 27; ss <= 30; ss++) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ss;

                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]+dx;
                        vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]-dx;

                        nv++;
                    }
                }
                else if (ir >= 20 && ir <= 26) {
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ir;
                    vals[nv]=1.0; nv++;
                }
                else if (ir == 27) { //for u_{n,r} equation
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //dF/dU_{O+,r}
                    xx[k][j][i].fx[7]=xx[k][j][i].fx[7]+dx;
                    vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                    xx[k][j][i].fx[7]=xx[k][j][i].fx[7]-dx;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10;  //dF/dU_{H+,r}
                    xx[k][j][i].fx[10]=xx[k][j][i].fx[10]+dx;
                    vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                    xx[k][j][i].fx[10]=xx[k][j][i].fx[10]-dx;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13;  //dF/dU_{He+,r}
                    xx[k][j][i].fx[13]=xx[k][j][i].fx[13]+dx;
                    vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                    xx[k][j][i].fx[13]=xx[k][j][i].fx[13]-dx;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=27;  //dF/dU_{n,r}
                    xx[k][j][i].fx[27]=xx[k][j][i].fx[27]+dx;
                    vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                    xx[k][j][i].fx[27]=xx[k][j][i].fx[27]-dx;
                    nv++;
                }
                else if (ir == 28) { //for u_{n,theta} equation
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //dF/dU_{O+,r}
                    xx[k][j][i].fx[8]=xx[k][j][i].fx[8]+dx;
                    vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                    xx[k][j][i].fx[8]=xx[k][j][i].fx[8]-dx;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11;  //dF/dU_{H+,r}
                    xx[k][j][i].fx[11]=xx[k][j][i].fx[11]+dx;
                    vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                    xx[k][j][i].fx[11]=xx[k][j][i].fx[11]-dx;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14;  //dF/dU_{He+,r}
                    xx[k][j][i].fx[14]=xx[k][j][i].fx[14]+dx;
                    vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                    xx[k][j][i].fx[14]=xx[k][j][i].fx[14]-dx;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=28;  //dF/dU_{n,r}
                    xx[k][j][i].fx[28]=xx[k][j][i].fx[28]+dx;
                    vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                    xx[k][j][i].fx[28]=xx[k][j][i].fx[28]-dx;
                    nv++;
                }
                else if (ir == 29) { //for u_{n,phi} equation
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //dF/dU_{O+,r}
                    xx[k][j][i].fx[9]=xx[k][j][i].fx[9]+dx;
                    vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                    xx[k][j][i].fx[9]=xx[k][j][i].fx[9]-dx;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12;  //dF/dU_{H+,r}
                    xx[k][j][i].fx[12]=xx[k][j][i].fx[12]+dx;
                    vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                    xx[k][j][i].fx[12]=xx[k][j][i].fx[12]-dx;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15;  //dF/dU_{He+,r}
                    xx[k][j][i].fx[15]=xx[k][j][i].fx[15]+dx;
                    vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                    xx[k][j][i].fx[15]=xx[k][j][i].fx[15]-dx;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=29;  //dF/dU_{n,r}
                    xx[k][j][i].fx[29]=xx[k][j][i].fx[29]+dx;
                    vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                    xx[k][j][i].fx[29]=xx[k][j][i].fx[29]-dx;
                    nv++;
                }
                else if (ir == 30) { //for T_n} equation
                    for (int ss = 7; ss <= 19; ss++) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ss;

                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]+dx;
                        vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]-dx;

                        nv++;
                    }

                    for (int ss = 27; ss <= 30; ss++) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ss;

                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]+dx;
                        vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]-dx;

                        nv++;
                    }
                }
                else if (ir == 31) { //for Br equation
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=31;  //dF/dB^r_{i,j,k}
                    vals[nv]= 1.0;
                    nv++;

                    col[nv].k=k; col[nv].j=j-1; col[nv].i=i; col[nv].c=36;  //dF/dB^phi_{i,j-1,k}
                    vals[nv]=-dt_quart/dth;
                    nv++;

                    col[nv].k=k; col[nv].j=j+1; col[nv].i=i; col[nv].c=36;  //dF/dB^phi_{i,j+1,k}
                    vals[nv]= dt_quart/dth;
                    nv++;

                    col[nv].k=k-1; col[nv].j=j; col[nv].i=i; col[nv].c=35;  //dF/dB^theta_{i,j,k-1}
                    vals[nv]= dt_quart/dph;
                    nv++;

                    col[nv].k=k+1; col[nv].j=j; col[nv].i=i; col[nv].c=35;  //dF/dB^theta_{i,j,k+1}
                    vals[nv]=-dt_quart/dph;
                    nv++;
                }
                else if (ir == 32) { //for Btheta equation
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=32;  //dF/dB^theta_{i,j,k}
                    vals[nv]= 1.0;
                    nv++;

                    col[nv].k=k-1; col[nv].j=j; col[nv].i=i; col[nv].c=34;  //dF/dB^r_{i,j,k-1}
                    vals[nv]=-dt_quart/dph;
                    nv++;

                    col[nv].k=k+1; col[nv].j=j; col[nv].i=i; col[nv].c=34;  //dF/dB^r_{i,j,k+1}
                    vals[nv]= dt_quart/dph;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i-1; col[nv].c=36;  //dF/dB^phi_{i-1,j,k}
                    vals[nv]= dt_quart/dr;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i+1; col[nv].c=36;  //dF/dB^phi_{i+1,j,k}
                    vals[nv]=-dt_quart/dr;
                    nv++;
                }
                else if (ir == 33) { //for Bphi equation
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=33;  //dF/dB^phi_{i,j,k}
                    vals[nv]= 1.0;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i-1; col[nv].c=35;  //dF/dB^theta_{i-1,j,k}
                    vals[nv]=-dt_quart/dr;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i+1; col[nv].c=35;  //dF/dB^theta_{i+1,j,k}
                    vals[nv]= dt_quart/dr;
                    nv++;

                    col[nv].k=k; col[nv].j=j-1; col[nv].i=i; col[nv].c=34;  //dF/dB^r_{i,j-1,k}
                    vals[nv]= dt_quart/dth;
                    nv++;

                    col[nv].k=k; col[nv].j=j+1; col[nv].i=i; col[nv].c=34;  //dF/dB^r_{i,j+1,k}
                    vals[nv]=-dt_quart/dth;
                    nv++;
                }
                else if (ir == 34) { //for E^r equation
                    for (int ss = 7; ss <= 15; ss++) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ss;

                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]+dx;
                        vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]-dx;

                        nv++;
                    }

                    for (int ss = 31; ss <= 33; ss++) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ss;

                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]+dx;
                        vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]-dx;

                        nv++;
                    }

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=34;  //dF/dE^r_{i,j,k}
                    vals[nv]=1.0;
                    nv++;
                }
                else if (ir == 35) { //for E^theta equation
                    for (int ss = 7; ss <= 15; ss++) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ss;

                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]+dx;
                        vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]-dx;

                        nv++;
                    }

                    for (int ss = 31; ss <= 33; ss++) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ss;

                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]+dx;
                        vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]-dx;

                        nv++;
                    }

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=35;  //dF/dE^theta_{i,j,k}
                    vals[nv]=1.0;
                    nv++;
                }
                else if (ir == 36) { //for E^phi equation
                    for (int ss = 7; ss <= 15; ss++) {
                        if (ss == 9 or ss == 12) continue;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ss;

                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]+dx;
                        vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]-dx;

                        nv++;
                    }

                    for (int ss = 31; ss <= 32; ss++) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ss;

                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]+dx;
                        vals[nv]=(functions(xx, xn, uu, i, j, k, xi, yj, zk, ir)-temp)/dx;
                        xx[k][j][i].fx[ss]=xx[k][j][i].fx[ss]-dx;

                        nv++;
                    }

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=36;  //dF/dE^phi_{i,j,k}
                    vals[nv]=1.0;
                    nv++;
                }

                for (s = 0; s < nv; s++) {
                    if (isnan(vals[s]) || isinf(vals[s])) {
                        cout<<"Jacobian value is Nan or inf at ("<<row.i<<", "<<row.j<<", "<<row.k<<", "<<row.c
                        <<col[nv].i<<", "<<col[nv].j<<", "<<col[nv].k<<", "<<col[nv].c<<")"<<endl;
                        exit(-1);
                    }
                }
                MatSetValuesStencil(Jac,1,&row,nv,col,vals,INSERT_VALUES);
            }
        }
      }
    }

    MatAssemblyBegin(Jac, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Jac, MAT_FINAL_ASSEMBLY);
    if (Jpre != Jac) {
        MatAssemblyBegin(Jpre, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(Jpre, MAT_FINAL_ASSEMBLY);
    }

    DMDAVecRestoreArray(da, localX, &xx);
    DMDAVecRestoreArray(da, localU, &uu);
    DMDAVecRestoreArray(da, localXn, &xn);

    DMRestoreLocalVector(da, &localX);
    DMRestoreLocalVector(da, &localU);
    DMRestoreLocalVector(da, &localXn);

    delete[] vals;
    delete[] col;

    return ierr;
}
