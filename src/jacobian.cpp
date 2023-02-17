#include <cmath>
#include <iostream>

using namespace std;

#include <petscts.h>
#include <petscdm.h>
#include <petscdmda.h>

#include "param.h"
#include "funcdef.h"
#include "ele_cooling_rate.h"

int jacobian(TS ts, double ftime, Vec X, Vec Xdt, double a, Mat Jac, Mat Jpre, void *ctx)
{
    DM         da;
    Vec        localX, localU;
    Field      ***xx, ***uu;
    AppCtx     *params = (AppCtx*)ctx;
    int        i, j, k, ir, nv, xs, xm, ys, ym, zs, zm, xi, yj, zk, s;
    MatStencil row, *col;
    PetscErrorCode ierr=0;
    double     *vals;
    const int  nzer=a4;

    double     temp, ne, rhon, Nn, ns[7], nsmore, nsmore_norm; //nn[7], 
    double     Bx, By, Bz, Br, Bt, Bp, uex, uey, uez, uir, uit, uip, rhos[7];

    const double four3rd=4.0/3.0; //two3rd=2.0/3.0, dT=1.0e-8;

    double temp1, uir_unr, uit_unt, uip_unp, rhossum_nusq[3], neme_nsms_nues; //, Ce, Cep, Te, Tn,
    double Bdve[3], neu_friction_coef[3], neu_tem_exchange_coef[4];

    params->ftime=ftime;

    TSGetDM(ts, &da);

    DMDAGetCorners(da, &xs ,&ys, &zs, &xm, &ym, &zm);

    DMGetLocalVector(da, &localX);
    DMGlobalToLocalBegin(da, X, INSERT_VALUES,localX);
    DMGlobalToLocalEnd(da, X, INSERT_VALUES, localX);
    DMDAVecGetArray(da, localX, &xx);

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

            for (s = 0; s < sl; s++) {
                ns[s]=uu[k][j][i].fx[s+27]; rhos[s]=ns[s]*ams[s];
            }

            ne=uu[k][j][i].fx[6]; Nn=uu[k][j][i].fx[10]; rhon=uu[k][j][i].fx[34];

            nsmore=(ns[0]+ns[3]+ns[4]+ns[5]+ns[6]);
            uir=(nsmore*xx[k][j][i].fx[7]+ns[1]*xx[k][j][i].fx[10]+ns[2]*xx[k][j][i].fx[13])/ne;
            uit=(nsmore*xx[k][j][i].fx[8]+ns[1]*xx[k][j][i].fx[11]+ns[2]*xx[k][j][i].fx[14])/ne;
            uip=(nsmore*xx[k][j][i].fx[9]+ns[1]*xx[k][j][i].fx[12]+ns[2]*xx[k][j][i].fx[15])/ne;

            nsmore_norm=nsmore/ne-1.0;

            //magnetic field components in terms of special spherical components
            Bx=( Jiv11[zk][yj]*xx[k][j][i].fx[31]+Jiv12[zk][yj][xi]*xx[k][j][i].fx[32]
                +Jiv13[zk][yj][xi]*xx[k][j][i].fx[33])/r2sintheta[yj][xi];
            By=( Jiv21[zk][yj]*xx[k][j][i].fx[31]+Jiv22[zk][yj][xi]*xx[k][j][i].fx[32]
                +Jiv23[zk][yj][xi] *xx[k][j][i].fx[33])/r2sintheta[yj][xi];
            Bz=(Jiv31[yj]*xx[k][j][i].fx[31]+Jiv32[yj][xi]*xx[k][j][i].fx[32])/r2sintheta[yj][xi];

            //uex, uey, uez in terms of uer, and ue_theta, and ue_phi
            uex=K11[zk][yj]*uir+K12[zk][yj]*uit+K13[zk]*uip;
            uey=K21[zk][yj]*uir+K22[zk][yj]*uit+K23[zk]*uip;
            uez=K31[yj]*uir+K32[yj]*uit;

            Br=uu[k][j][i].fx[0]; Bt=uu[k][j][i].fx[1]; Bp=uu[k][j][i].fx[2];

            rhossum_nusq[0]=( rhos[0]*nust[zk][yj][xi][11]+rhos[3]*nust[zk][yj][xi][30]
                             +rhos[4]*nust[zk][yj][xi][32]+rhos[5]*nust[zk][yj][xi][34]
                             +rhos[6]*nust[zk][yj][xi][36])/rhon;
            rhossum_nusq[1]=rhos[1]*nust[zk][yj][xi][19]/rhon;
            rhossum_nusq[2]=rhos[2]*nust[zk][yj][xi][27]/rhon;

            for (ir = 0; ir < a4; ir++) {
                row.c=ir; nv=0;

                if (ir < sl) {
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ir;
                    vals[nv]=a; nv++;
                }
                else if (ir == 7) {
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;
                    vals[nv]= a+nust[zk][yj][xi][11]+nust[zk][yj][xi][7]+nust[zk][yj][xi][8];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;
                    vals[nv]= qms[0]*nsmore_norm*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;
                    vals[nv]=-qms[0]*nsmore_norm*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10;
                    vals[nv]=-nust[zk][yj][xi][7];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11;
                    vals[nv]= qms[0]*ns[1]/ne*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12;
                    vals[nv]=-qms[0]*ns[1]/ne*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13;
                    vals[nv]=-nust[zk][yj][xi][8];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14;
                    vals[nv]= qms[0]*ns[2]/ne*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15;
                    vals[nv]=-qms[0]*ns[2]/ne*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=27;
                    vals[nv]=-nust[zk][yj][xi][11];
                    nv++;
                }
                else if (ir == 8) {
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;
                    vals[nv]=-qms[0]*nsmore_norm*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;
                    vals[nv]= a+nust[zk][yj][xi][11]+nust[zk][yj][xi][7]+nust[zk][yj][xi][8];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;
                    vals[nv]= qms[0]*nsmore_norm*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10;
                    vals[nv]=-qms[0]*ns[1]/ne*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11;
                    vals[nv]=-nust[zk][yj][xi][7];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12;
                    vals[nv]= qms[0]*ns[1]/ne*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13;
                    vals[nv]=-qms[0]*ns[2]/ne*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14;
                    vals[nv]=-nust[zk][yj][xi][8];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15;
                    vals[nv]= qms[0]*ns[2]/ne*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=28;
                    vals[nv]=-nust[zk][yj][xi][11];
                    nv++;
                }
                else if (ir == 9) {
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;
                    vals[nv]= qms[0]*nsmore_norm*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;
                    vals[nv]=-qms[0]*nsmore_norm*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;
                    vals[nv]= a+nust[zk][yj][xi][11]+nust[zk][yj][xi][7]+nust[zk][yj][xi][8];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10;
                    vals[nv]= qms[0]*ns[1]/ne*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11;
                    vals[nv]=-qms[0]*ns[1]/ne*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12;
                    vals[nv]=-nust[zk][yj][xi][7];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13;
                    vals[nv]= qms[0]*ns[2]/ne*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14;
                    vals[nv]=-qms[0]*ns[2]/ne*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15;
                    vals[nv]=-nust[zk][yj][xi][8];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=29;
                    vals[nv]=-nust[zk][yj][xi][11];
                    nv++;
                }
                else if (ir == 10) {
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //u_{O+,r}
                    vals[nv]=-nust[zk][yj][xi][14];
                    nv++;

                    temp=qms[1]*nsmore/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //u_{O+,theta}
                    vals[nv]= temp*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //u_{O+,phi}
                    vals[nv]=-temp*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10;  //u_{H+,r}
                    vals[nv]=a+nust[zk][yj][xi][19]+nust[zk][yj][xi][14]+nust[zk][yj][xi][15];
                    nv++;

                    temp=qms[1]*(ns[1]/ne-1.0);
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11; //u_{H+,theta}
                    vals[nv]= temp*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12; //u_{H+,phi}
                    vals[nv]=-temp*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13;
                    vals[nv]=-nust[zk][yj][xi][15];
                    nv++;

                    temp=qms[1]*ns[2]/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14;
                    vals[nv]= temp*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15;
                    vals[nv]=-temp*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=27;
                    vals[nv]=-nust[zk][yj][xi][19];
                    nv++;
                }
                else if (ir == 11) {
                    temp=qms[1]*nsmore/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //u_{O+,r}
                    vals[nv]=-temp*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //u_{O+,theta}
                    vals[nv]=-nust[zk][yj][xi][14];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //u_{O+,phi}
                    vals[nv]= temp*Br;
                    nv++;

                    temp=qms[1]*(ns[1]/ne-1.0);
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10; //u_{H+,r}
                    vals[nv]=-temp*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11;  //u_{H+,theta}
                    vals[nv]=a+nust[zk][yj][xi][19]+nust[zk][yj][xi][14]+nust[zk][yj][xi][15];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12; //u_{H+,phi}
                    vals[nv]= temp*Br;
                    nv++;

                    temp=qms[1]*ns[2]/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13;
                    vals[nv]=-temp*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14;
                    vals[nv]=-nust[zk][yj][xi][15];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15;
                    vals[nv]= temp*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=28;
                    vals[nv]=-nust[zk][yj][xi][19];
                    nv++;
                }
                else if (ir == 12) {
                    temp=qms[1]*nsmore/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //u_{O+,r}
                    vals[nv]= temp*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //u_{O+,theta}
                    vals[nv]=-temp*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //u_{O+,phi}
                    vals[nv]=-nust[zk][yj][xi][14];
                    nv++;

                    temp=qms[1]*(ns[1]/ne-1.0);
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10; //u_{H+,r}
                    vals[nv]= temp*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11; //u_{H+,theta}
                    vals[nv]=-temp*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12;  //u_{H+,phi}
                    vals[nv]=a+nust[zk][yj][xi][19]+nust[zk][yj][xi][14]+nust[zk][yj][xi][15];
                    nv++;

                    temp=qms[1]*ns[2]/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13;
                    vals[nv]= temp*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14;
                    vals[nv]=-temp*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15;
                    vals[nv]=-nust[zk][yj][xi][15];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=29;
                    vals[nv]=-nust[zk][yj][xi][19];
                    nv++;
                }
                else if (ir == 13) {
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //u_{O+,r}
                    vals[nv]=-nust[zk][yj][xi][22];
                    nv++;

                    temp=qms[2]*nsmore/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //u_{O+,theta}
                    vals[nv]= temp*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //u_{O+,phi}
                    vals[nv]=-temp*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10; //u_{H+,r}
                    vals[nv]=-nust[zk][yj][xi][23];
                    nv++;

                    temp=qms[2]*ns[1]/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11; //u_{H+,theta}
                    vals[nv]= temp*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12; //u_{H+,phi}
                    vals[nv]=-temp*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13;  //u_{He+,r}
                    vals[nv]=a+nust[zk][yj][xi][27]+nust[zk][yj][xi][22]+nust[zk][yj][xi][23];
                    nv++;

                    temp=qms[2]*(ns[2]/ne-1.0);
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14; //u_{He+,theta}
                    vals[nv]= temp*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15; //u_{He+,phi}
                    vals[nv]=-temp*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=27;
                    vals[nv]=-nust[zk][yj][xi][27];
                    nv++;
                }
                else if (ir == 14) {
                    temp=qms[2]*nsmore/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //u_{O+,r}
                    vals[nv]=-temp*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //u_{O+,theta}
                    vals[nv]=-nust[zk][yj][xi][22];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //u_{O+,phi}
                    vals[nv]= temp*Br;
                    nv++;

                    temp=qms[2]*ns[1]/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10; //u_{H+,r}
                    vals[nv]=-temp*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11; //u_{H+,theta}
                    vals[nv]=-nust[zk][yj][xi][23];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12; //u_{H+,phi}
                    vals[nv]= temp*Br;
                    nv++;

                    temp=qms[2]*(ns[2]/ne-1.0);
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13; //u_{He+,theta}
                    vals[nv]=-temp*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14;  //u_{He+,theta}
                    vals[nv]=a+nust[zk][yj][xi][27]+nust[zk][yj][xi][22]+nust[zk][yj][xi][23];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15; //u_{He+,phi}
                    vals[nv]= temp*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=28;
                    vals[nv]=-nust[zk][yj][xi][27];
                    nv++;
                }
                else if (ir == 15) {
                    temp=qms[2]*nsmore/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //u_{O+,r}
                    vals[nv]= temp*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //u_{O+,theta}
                    vals[nv]=-temp*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //u_{O+,phi}
                    vals[nv]=-nust[zk][yj][xi][22];
                    nv++;

                    temp=qms[2]*ns[1]/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10; //u_{H+,r}
                    vals[nv]= temp*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11; //u_{H+,theta}
                    vals[nv]=-temp*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12; //u_{H+,phi}
                    vals[nv]=-nust[zk][yj][xi][23];
                    nv++;

                    temp=qms[2]*(ns[2]/ne-1.0);
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13; //u_{He+,theta}
                    vals[nv]= temp*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14; //u_{He+,theta}
                    vals[nv]=-temp*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15;  //u_{He+,phi}
                    vals[nv]=a+nust[zk][yj][xi][27]+nust[zk][yj][xi][22]+nust[zk][yj][xi][23];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=29;
                    vals[nv]=-nust[zk][yj][xi][27];
                    nv++;
                }
                else if (ir == 16) {  //O+ temperature equation
                    temp=four3rd*ams[0];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //U_{O+,r}
                    vals[nv]=-temp*( nust[zk][yj][xi][13]*(xx[k][j][i].fx[7]-xx[k][j][i].fx[27])
                                    +ams[1]*nust[zk][yj][xi][9]*(xx[k][j][i].fx[7]-xx[k][j][i].fx[10])
                                    +ams[2]*nust[zk][yj][xi][10]*(xx[k][j][i].fx[7]-xx[k][j][i].fx[13]));
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //U_{O+,theta}
                    vals[nv]=-temp*( nust[zk][yj][xi][13]*(xx[k][j][i].fx[8]-xx[k][j][i].fx[28])
                                    +ams[1]*nust[zk][yj][xi][9]*(xx[k][j][i].fx[8]-xx[k][j][i].fx[11])
                                    +ams[2]*nust[zk][yj][xi][10]*(xx[k][j][i].fx[8]-xx[k][j][i].fx[14]));
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //U_{O+,phi}
                    vals[nv]=-temp*( nust[zk][yj][xi][13]*(xx[k][j][i].fx[9]-xx[k][j][i].fx[29])
                                    +ams[1]*nust[zk][yj][xi][9]*(xx[k][j][i].fx[9]-xx[k][j][i].fx[12])
                                    +ams[2]*nust[zk][yj][xi][10]*(xx[k][j][i].fx[9]-xx[k][j][i].fx[15]));
                    nv++;

                    temp1=temp*ams[1]*nust[zk][yj][xi][9];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10;  //U_{H+,r}
                    vals[nv]= temp1*(xx[k][j][i].fx[7]-xx[k][j][i].fx[10]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11;  //U_{H+,theta}
                    vals[nv]= temp1*(xx[k][j][i].fx[8]-xx[k][j][i].fx[11]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12;  //U_{He+,phi}
                    vals[nv]= temp1*(xx[k][j][i].fx[9]-xx[k][j][i].fx[12]);
                    nv++;

                    temp1=temp*ams[2]*nust[zk][yj][xi][10];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13;  //U_{He+,r}
                    vals[nv]= temp1*(xx[k][j][i].fx[7]-xx[k][j][i].fx[13]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14;  //U_{He+,theta}
                    vals[nv]= temp1*(xx[k][j][i].fx[8]-xx[k][j][i].fx[14]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15;  //U_{He+,phi}
                    vals[nv]= temp1*(xx[k][j][i].fx[9]-xx[k][j][i].fx[15]);
                    nv++;

                    temp=2.0*ams[0];
                    neme_nsms_nues=2.0*ne*ame/rhos[0]*nust[zk][yj][xi][0];

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=16;  //T_{O+}
                    vals[nv]=a+temp*(nust[zk][yj][xi][12]+nust[zk][yj][xi][9]+nust[zk][yj][xi][10])
                              +neme_nsms_nues;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=17;  //T_{H+}
                    vals[nv]=-temp*nust[zk][yj][xi][9];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=18;  //T_{He+}
                    vals[nv]=-temp*nust[zk][yj][xi][10];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=19;  //T_{e}
                    vals[nv]=-neme_nsms_nues;
                    nv++;

                    temp=four3rd*ams[0]*nust[zk][yj][xi][13];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=27;  //U_{n,r}
                    vals[nv]= temp*(xx[k][j][i].fx[7]-xx[k][j][i].fx[27]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=28;  //U_{n,theta}
                    vals[nv]= temp*(xx[k][j][i].fx[8]-xx[k][j][i].fx[28]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=29;  //U_{n,phi}
                    vals[nv]= temp*(xx[k][j][i].fx[9]-xx[k][j][i].fx[29]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=30;  //T_{n}
                    vals[nv]=-2.0*ams[0]*nust[zk][yj][xi][12];
                    nv++;
                }
                else if (ir == 17) {  //H+ temperature equation
                    temp=four3rd*ams[1]*nust[zk][yj][xi][18];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //dF/dU_{O+,r}
                    vals[nv]= temp*(xx[k][j][i].fx[10]-xx[k][j][i].fx[7]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //dF/dU_{O+,theta}
                    vals[nv]= temp*(xx[k][j][i].fx[11]-xx[k][j][i].fx[8]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //dF/dU_{O+,phi}
                    vals[nv]= temp*(xx[k][j][i].fx[12]-xx[k][j][i].fx[9]);
                    nv++;

                    temp=four3rd*ams[1];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10;  //dF/dU_{H+,r}
                    vals[nv]=-temp*( nust[zk][yj][xi][21]*(xx[k][j][i].fx[10]-xx[k][j][i].fx[27])
                                    +nust[zk][yj][xi][18]*(xx[k][j][i].fx[10]-xx[k][j][i].fx[7])
                                    +ams[2]*nust[zk][yj][xi][17]*(xx[k][j][i].fx[10]-xx[k][j][i].fx[13]));
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11;  //dF/dU_{H+,theta}
                    vals[nv]=-temp*( nust[zk][yj][xi][21]*(xx[k][j][i].fx[11]-xx[k][j][i].fx[28])
                                    +nust[zk][yj][xi][18]*(xx[k][j][i].fx[11]-xx[k][j][i].fx[8])
                                    +ams[2]*nust[zk][yj][xi][17]*(xx[k][j][i].fx[11]-xx[k][j][i].fx[14]));
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12;  //dF/dU_{H+,phi}
                    vals[nv]=-temp*( nust[zk][yj][xi][21]*(xx[k][j][i].fx[12]-xx[k][j][i].fx[29])
                                    +nust[zk][yj][xi][18]*(xx[k][j][i].fx[12]-xx[k][j][i].fx[9])
                                    +ams[2]*nust[zk][yj][xi][17]*(xx[k][j][i].fx[12]-xx[k][j][i].fx[15]));
                    nv++;

                    temp=four3rd*ams[2]*nust[zk][yj][xi][17];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13;  //dF/dU_{He+,r}
                    vals[nv]= temp*(xx[k][j][i].fx[10]-xx[k][j][i].fx[13]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14;  //dF/dU_{He+,theta}
                    vals[nv]= temp*(xx[k][j][i].fx[11]-xx[k][j][i].fx[14]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15;  //dF/dU_{He+,r}
                    vals[nv]= temp*(xx[k][j][i].fx[12]-xx[k][j][i].fx[15]);
                    nv++;

                    temp=2.0*ams[1];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=16;  //dF/dT_{O+}
                    vals[nv]=-temp*nust[zk][yj][xi][16];
                    nv++;

                    neme_nsms_nues=2.0*ne*ame/rhos[1]*nust[zk][yj][xi][1];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=17;  //dF/dT_{H+}
                    vals[nv]=a+temp*(nust[zk][yj][xi][20]+nust[zk][yj][xi][16]+nust[zk][yj][xi][17])
                              +neme_nsms_nues;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=18;  //dF/dT_{He+}
                    vals[nv]=-temp*nust[zk][yj][xi][17];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=19;  //dF/dT_{e}
                    vals[nv]=-neme_nsms_nues;
                    nv++;

                    temp=four3rd*ams[1]*nust[zk][yj][xi][21];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=27;  //dF/dU_{n,r}
                    vals[nv]= temp*(xx[k][j][i].fx[10]-xx[k][j][i].fx[27]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=28;  //dF/dU_{n,theta}
                    vals[nv]= temp*(xx[k][j][i].fx[11]-xx[k][j][i].fx[28]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=29;  //dF/dU_{n,phi}
                    vals[nv]= temp*(xx[k][j][i].fx[12]-xx[k][j][i].fx[29]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=30;  //dF/dT_{n}
                    vals[nv]=-2.0*ams[1]*nust[zk][yj][xi][20];
                    nv++;
                }
                else if (ir == 18) {  //He+ temperature equation
                    temp=four3rd*ams[2]*nust[zk][yj][xi][26];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //dF/dU_{O+,r}
                    vals[nv]= temp*(xx[k][j][i].fx[13]-xx[k][j][i].fx[7]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //dF/dU_{O+,theta}
                    vals[nv]= temp*(xx[k][j][i].fx[14]-xx[k][j][i].fx[8]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //dF/dU_{O+,phi}
                    vals[nv]= temp*(xx[k][j][i].fx[15]-xx[k][j][i].fx[9]);
                    nv++;

                    temp=four3rd*ams[1]*nust[zk][yj][xi][25];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10;  //dF/dU_{H+,r}
                    vals[nv]= temp*(xx[k][j][i].fx[13]-xx[k][j][i].fx[10]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11;  //dF/dU_{H+,theta}
                    vals[nv]= temp*(xx[k][j][i].fx[14]-xx[k][j][i].fx[11]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12;  //dF/dU_{H+,phi}
                    vals[nv]= temp*(xx[k][j][i].fx[15]-xx[k][j][i].fx[12]);
                    nv++;

                    temp=four3rd*ams[2];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13;  //dF/dU_{He+,r}
                    vals[nv]=-temp*( nust[zk][yj][xi][29]*(xx[k][j][i].fx[13]-xx[k][j][i].fx[27])
                                    +nust[zk][yj][xi][26]*(xx[k][j][i].fx[13]-xx[k][j][i].fx[7])
                                    +ams[1]*nust[zk][yj][xi][25]*(xx[k][j][i].fx[13]-xx[k][j][i].fx[10]));
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14;  //dF/dU_{He+,theta}
                    vals[nv]=-temp*( nust[zk][yj][xi][29]*(xx[k][j][i].fx[14]-xx[k][j][i].fx[28])
                                    +nust[zk][yj][xi][26]*(xx[k][j][i].fx[14]-xx[k][j][i].fx[8])
                                    +ams[1]*nust[zk][yj][xi][25]*(xx[k][j][i].fx[14]-xx[k][j][i].fx[11]));
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15;  //dF/dU_{He+,r}
                    vals[nv]=-temp*( nust[zk][yj][xi][29]*(xx[k][j][i].fx[15]-xx[k][j][i].fx[29])
                                    +nust[zk][yj][xi][26]*(xx[k][j][i].fx[15]-xx[k][j][i].fx[9])
                                    +ams[1]*nust[zk][yj][xi][25]*(xx[k][j][i].fx[15]-xx[k][j][i].fx[12]));
                    nv++;

                    temp=2.0*ams[2];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=16;  //dF/dT_{O+}
                    vals[nv]=-temp*nust[zk][yj][xi][24];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=17;  //dF/dT_{H+}
                    vals[nv]=-temp*nust[zk][yj][xi][25];
                    nv++;

                    neme_nsms_nues=2.0*ne*ame/rhos[2]*nust[zk][yj][xi][2];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=18;  //dF/dT_{He+}
                    vals[nv]=a+temp*(nust[zk][yj][xi][28]+nust[zk][yj][xi][24]+nust[zk][yj][xi][25])
                              +neme_nsms_nues;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=19;  //dF/dT_{e}
                    vals[nv]=-neme_nsms_nues;
                    nv++;

                    temp=four3rd*ams[2]*nust[zk][yj][xi][29];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=27;  //dF/dU_{n,r}
                    vals[nv]= temp*(xx[k][j][i].fx[13]-xx[k][j][i].fx[27]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=28;  //dF/dU_{n,theta}
                    vals[nv]= temp*(xx[k][j][i].fx[14]-xx[k][j][i].fx[28]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=29;  //dF/dU_{n,phi}
                    vals[nv]= temp*(xx[k][j][i].fx[15]-xx[k][j][i].fx[29]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=30;  //dF/dT_{n}
                    vals[nv]=-2.0*ams[2]*nust[zk][yj][xi][28];
                    nv++;
                }
                else if (ir == 19) {  //electron temperature equation
                    uir_unr=uir-xx[k][j][i].fx[27]; uit_unt=uit-xx[k][j][i].fx[28];
                    uip_unp=uip-xx[k][j][i].fx[29];

                    temp1=four3rd*ame*nust[zk][yj][xi][5];
                    temp=temp1*nsmore/ne;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //dF/dU_{O+,r}
                    vals[nv]=-temp*uir_unr;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //dF/dU_{O+,theta}
                    vals[nv]=-temp*uit_unt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //dF/dU_{O+,phi}
                    vals[nv]=-temp*uip_unp;
                    nv++;

                    temp=temp1*ns[1]/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10;  //dF/dU_{H+,r}
                    vals[nv]=-temp*uir_unr;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11;  //dF/dU_{H+,theta}
                    vals[nv]=-temp*uit_unt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12;  //dF/dU_{H+,phi}
                    vals[nv]=-temp*uip_unp;
                    nv++;

                    temp=temp1*ns[2]/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13;  //dF/dU_{He+,r}
                    vals[nv]=-temp*uir_unr;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14;  //dF/dU_{He+,theta}
                    vals[nv]=-temp*uit_unt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15;  //dF/dU_{He+,r}
                    vals[nv]=-temp*uip_unp;
                    nv++;

                    temp=2.0*ame;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=16;  //dF/dT_{O+}
                    vals[nv]=-temp*nust[zk][yj][xi][4];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=17;  //dF/dT_{H+}
                    vals[nv]=-temp*nust[zk][yj][xi][1]/ams[1];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=18;  //dF/dT_{He+}
                    vals[nv]=-temp*nust[zk][yj][xi][2]/ams[2];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=19;  //dF/dT_{e}
                    vals[nv]= a+temp*( nust[zk][yj][xi][4]+nust[zk][yj][xi][1]/ams[1]
                                      +nust[zk][yj][xi][2]/ams[2]+nust[zk][yj][xi][6]);
                    nv++;

                    temp=four3rd*ame*nust[zk][yj][xi][5];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=27;  //dF/dU_{n,r}
                    vals[nv]= temp*uir_unr;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=28;  //dF/dU_{n,theta}
                    vals[nv]= temp*uit_unt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=29;  //dF/dU_{n,phi}
                    vals[nv]= temp*uip_unp;
                    nv++; 

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=30;  //dF/dT_{n}
                    vals[nv]=-2.0*ame*nust[zk][yj][xi][6];
                    nv++;
                }
                else if (ir >= 20 && ir <= 26) {
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ir;
                    vals[nv]=a; nv++;
                }
                else if (ir == 27) { //for u_{n,r} equation
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //dF/dU_{O+,r}
                    vals[nv]=-rhossum_nusq[0];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10;  //dF/dU_{H+,r}
                    vals[nv]=-rhossum_nusq[1];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13;  //dF/dU_{He+,r}
                    vals[nv]=-rhossum_nusq[2];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=27;  //dF/dU_{n,r}
                    vals[nv]= a+rhossum_nusq[0]+rhossum_nusq[1]+rhossum_nusq[2];
                    nv++;
                }
                else if (ir == 28) { //for u_{n,theta} equation
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //dF/dU_{O+,r}
                    vals[nv]=-rhossum_nusq[0];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11;  //dF/dU_{H+,r}
                    vals[nv]=-rhossum_nusq[1];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14;  //dF/dU_{He+,r}
                    vals[nv]=-rhossum_nusq[2];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=28;  //dF/dU_{n,r}
                    vals[nv]= a+rhossum_nusq[0]+rhossum_nusq[1]+rhossum_nusq[2];
                    nv++;
                }
                else if (ir == 29) { //for u_{n,phi} equation
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //dF/dU_{O+,r}
                    vals[nv]=-rhossum_nusq[0];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12;  //dF/dU_{H+,r}
                    vals[nv]=-rhossum_nusq[1];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15;  //dF/dU_{He+,r}
                    vals[nv]=-rhossum_nusq[2];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=29;  //dF/dU_{n,r}
                    vals[nv]= a+rhossum_nusq[0]+rhossum_nusq[1]+rhossum_nusq[2];
                    nv++;
                }
                else if (ir == 30) { //for T_n} equation
                    neu_friction_coef[0]
                        = four3rd/Nn*( rhos[0]*ams[0]*nust[zk][yj][xi][12]+rhos[3]*ams[3]*nust[zk][yj][xi][31]
                                      +rhos[4]*ams[4]*nust[zk][yj][xi][33]+rhos[5]*ams[5]*nust[zk][yj][xi][35]
                                      +rhos[6]*ams[6]*nust[zk][yj][xi][37]);

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //dF/dU_{O+,r}
                    vals[nv]=-neu_friction_coef[0]*(xx[k][j][i].fx[7]-xx[k][j][i].fx[27]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //dF/dU_{O+,theta}
                    vals[nv]=-neu_friction_coef[0]*(xx[k][j][i].fx[8]-xx[k][j][i].fx[28]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //dF/dU_{O+,phi}
                    vals[nv]=-neu_friction_coef[0]*(xx[k][j][i].fx[9]-xx[k][j][i].fx[29]);
                    nv++;

                    neu_friction_coef[1]= four3rd/Nn*rhos[1]*ams[1]*nust[zk][yj][xi][20];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10;  //dF/dU_{H+,r}
                    vals[nv]=-neu_friction_coef[1]*(xx[k][j][i].fx[10]-xx[k][j][i].fx[27]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11;  //dF/dU_{H+,theta}
                    vals[nv]=-neu_friction_coef[1]*(xx[k][j][i].fx[11]-xx[k][j][i].fx[28]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12;  //dF/dU_{H+,phi}
                    vals[nv]=-neu_friction_coef[1]*(xx[k][j][i].fx[12]-xx[k][j][i].fx[29]);
                    nv++;

                    neu_friction_coef[2]=four3rd/Nn*rhos[2]*ams[2]*nust[zk][yj][xi][28];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13;  //dF/dU_{He+,r}
                    vals[nv]=-neu_friction_coef[2]*(xx[k][j][i].fx[13]-xx[k][j][i].fx[27]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14;  //dF/dU_{He+,theta}
                    vals[nv]=-neu_friction_coef[2]*(xx[k][j][i].fx[14]-xx[k][j][i].fx[28]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15;  //dF/dU_{He+,phi}
                    vals[nv]=-neu_friction_coef[2]*(xx[k][j][i].fx[15]-xx[k][j][i].fx[29]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=27;  //dF/dU_{n,r}
                    vals[nv]= neu_friction_coef[0]*(xx[k][j][i].fx[7]-xx[k][j][i].fx[27])
                             +neu_friction_coef[1]*(xx[k][j][i].fx[10]-xx[k][j][i].fx[27])
                             +neu_friction_coef[2]*(xx[k][j][i].fx[13]-xx[k][j][i].fx[27]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=28;  //dF/dU_{n,theta}
                    vals[nv]= neu_friction_coef[0]*(xx[k][j][i].fx[8]-xx[k][j][i].fx[28])
                             +neu_friction_coef[1]*(xx[k][j][i].fx[11]-xx[k][j][i].fx[28])
                             +neu_friction_coef[2]*(xx[k][j][i].fx[14]-xx[k][j][i].fx[28]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=29;  //dF/dU_{n,phi}
                    vals[nv]= neu_friction_coef[0]*(xx[k][j][i].fx[9]-xx[k][j][i].fx[29])
                             +neu_friction_coef[1]*(xx[k][j][i].fx[12]-xx[k][j][i].fx[29])
                             +neu_friction_coef[2]*(xx[k][j][i].fx[15]-xx[k][j][i].fx[29]);
                    nv++;

                    neu_tem_exchange_coef[0]
                        = 2.0/Nn*( rhos[0]*nust[zk][yj][xi][12]+rhos[3]*nust[zk][yj][xi][31]
                                  +rhos[4]*nust[zk][yj][xi][33]+rhos[5]*nust[zk][yj][xi][35]
                                  +rhos[6]*nust[zk][yj][xi][37]);
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=16;  //dF/dT_{O+}}
                    vals[nv]=-neu_tem_exchange_coef[0];
                    nv++;

                    neu_tem_exchange_coef[1]= 2.0/Nn*rhos[1]*nust[zk][yj][xi][20];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=17;  //dF/dT_{H+}
                    vals[nv]=-neu_tem_exchange_coef[1];
                    nv++;

                    neu_tem_exchange_coef[2]= 2.0/Nn*rhos[2]*nust[zk][yj][xi][28];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=18;  //dF/dT_{He+}
                    vals[nv]=-neu_tem_exchange_coef[2];
                    nv++;

                    neu_tem_exchange_coef[3]= 2.0*ne*ame/Nn*nust[zk][yj][xi][6];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=19;  //dF/dT_{e}
                    vals[nv]=-neu_tem_exchange_coef[3];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=30;  //dF/dT_{n}
                    vals[nv]= a+neu_tem_exchange_coef[0]+neu_tem_exchange_coef[1]+neu_tem_exchange_coef[2]
                               +neu_tem_exchange_coef[3];
                    nv++;
               }
                else if (ir == 31) { //for Br equation
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=31;  //dF/dB^r_{i,j,k}
                    vals[nv]= a;
                    nv++;

                    col[nv].k=k; col[nv].j=j-1; col[nv].i=i; col[nv].c=36;  //dF/dB^phi_{i,j-1,k}
                    vals[nv]=-0.5/dth;
                    nv++;

                    col[nv].k=k; col[nv].j=j+1; col[nv].i=i; col[nv].c=36;  //dF/dB^phi_{i,j+1,k}
                    vals[nv]= 0.5/dth;
                    nv++;

                    col[nv].k=k-1; col[nv].j=j; col[nv].i=i; col[nv].c=35;  //dF/dB^theta_{i,j,k-1}
                    vals[nv]= 0.5/dph;
                    nv++;

                    col[nv].k=k+1; col[nv].j=j; col[nv].i=i; col[nv].c=35;  //dF/dB^theta_{i,j,k+1}
                    vals[nv]=-0.5/dph;
                    nv++;
                }
                else if (ir == 32) { //for Btheta equation
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=32;  //dF/dB^theta_{i,j,k}
                    vals[nv]= a;
                    nv++;

                    col[nv].k=k-1; col[nv].j=j; col[nv].i=i; col[nv].c=34;  //dF/dB^r_{i,j,k-1}
                    vals[nv]=-0.5/dph;
                    nv++;

                    col[nv].k=k+1; col[nv].j=j; col[nv].i=i; col[nv].c=34;  //dF/dB^r_{i,j,k+1}
                    vals[nv]= 0.5/dph;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i-1; col[nv].c=36;  //dF/dB^phi_{i-1,j,k}
                    vals[nv]= 0.5/dr;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i+1; col[nv].c=36;  //dF/dB^phi_{i+1,j,k}
                    vals[nv]=-0.5/dr;
                    nv++;
                }
                else if (ir == 33) { //for Bphi equation
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=33;  //dF/dB^phi_{i,j,k}
                    vals[nv]= a;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i-1; col[nv].c=35;  //dF/dB^theta_{i-1,j,k}
                    vals[nv]=-0.5/dr;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i+1; col[nv].c=35;  //dF/dB^theta_{i+1,j,k}
                    vals[nv]= 0.5/dr;
                    nv++;

                    col[nv].k=k; col[nv].j=j-1; col[nv].i=i; col[nv].c=34;  //dF/dB^r_{i,j-1,k}
                    vals[nv]= 0.5/dth;
                    nv++;

                    col[nv].k=k; col[nv].j=j+1; col[nv].i=i; col[nv].c=34;  //dF/dB^r_{i,j+1,k}
                    vals[nv]=-0.5/dth;
                    nv++;
                }
                else if (ir == 34) { //for E^r equation
                    temp=nsmore/ne;
                    Bdve[0]=-( Jiv11[zk][yj]*(By*K31[yj] - Bz*K21[zk][yj])
                              +Jiv21[zk][yj]*(Bz*K11[zk][yj] - Bx*K31[yj])
                              +Jiv31[yj]*(Bx*K21[zk][yj] - By*K11[zk][yj]));
                    Bdve[1]=-( Jiv11[zk][yj]*(By*K32[yj] - Bz*K22[zk][yj])
                              +Jiv21[zk][yj]*(Bz*K12[zk][yj] - Bx*K32[yj])
                              +Jiv31[yj]*(Bx*K22[zk][yj] - By*K12[zk][yj]));
                    Bdve[2]=-( (Jiv21[zk][yj]*K13[zk] - Jiv11[zk][yj]*K23[zk])*Bz
                              +Jiv31[yj]*(Bx*K23[zk] - By*K13[zk]));

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //dF/du_{O+,r,i,j,k}
                    vals[nv]=Bdve[0]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //dF/du_{O+,theta,i,j,k}
                    vals[nv]=Bdve[1]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //dF/du_{O+,phi,i,j,k}
                    vals[nv]=Bdve[2]*temp;
                    nv++;

                    temp=ns[1]/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10;  //dF/du_{H+,r,i,j,k}
                    vals[nv]=Bdve[0]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11;  //dF/du_{H+,theta,i,j,k}
                    vals[nv]=Bdve[1]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12;  //dF/du_{H+,phi,i,j,k}
                    vals[nv]=Bdve[2]*temp;
                    nv++;

                    temp=ns[2]/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13;  //dF/du_{He+,r,i,j,k}
                    vals[nv]=Bdve[0]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14;  //dF/du_{He+,theta,i,j,k}
                    vals[nv]=Bdve[1]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15;  //dF/du_{He+,phi,i,j,k}
                    vals[nv]=Bdve[2]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=31;  //dF/dB^r_{i,j,k}
                    vals[nv]=-( Jiv11[zk][yj]*(Jiv21[zk][yj]*uez - Jiv31[yj]*uey)
                               +Jiv21[zk][yj]*(Jiv31[yj]*uex - Jiv11[zk][yj]*uez)
                               +Jiv31[yj]*(Jiv11[zk][yj]*uey - Jiv21[zk][yj]*uex))/r2sintheta[yj][xi];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=32;  //dF/dB^theta_{i,j,k}
                    vals[nv]=-( Jiv11[zk][yj]*(Jiv22[zk][yj][xi]*uez - Jiv32[yj][xi] *uey)
                               +Jiv21[zk][yj]*(Jiv32[yj][xi]*uex - Jiv12[zk][yj][xi]*uez)
                               +Jiv31[yj]*(Jiv12[zk][yj][xi]*uey - Jiv22[zk][yj][xi]*uex))/r2sintheta[yj][xi];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=33;  //dF/dB^phi_{i,j,k}
                    vals[nv]=-( (Jiv11[zk][yj]*Jiv23[zk][yj][xi]-Jiv21[zk][yj]*Jiv13[zk][yj][xi])*uez
                               +Jiv31[yj]*(Jiv13[zk][yj][xi]*uey - Jiv23[zk][yj][xi]*uex))/r2sintheta[yj][xi];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=34;  //dF/dE^r_{i,j,k}
                    vals[nv]=1.0;
                    nv++;
                }
                else if (ir == 35) { //for E^theta equation
                    temp=nsmore/ne;
                    Bdve[0]=-( Jiv12[zk][yj][xi]*(By*K31[yj] - Bz*K21[zk][yj])
                              +Jiv22[zk][yj][xi]*(Bz*K11[zk][yj] - Bx*K31[yj])
                              +Jiv32[yj][xi]*(Bx*K21[zk][yj] - By*K11[zk][yj]));
                    Bdve[1]=-( Jiv12[zk][yj][xi]*(By*K32[yj] - Bz*K22[zk][yj])
                              +Jiv22[zk][yj][xi]*(Bz*K12[zk][yj] - Bx*K32[yj])
                              +Jiv32[yj][xi]*(Bx*K22[zk][yj] - By*K12[zk][yj]));
                    Bdve[2]=-( (Jiv22[zk][yj][xi]*K13[zk] - Jiv12[zk][yj][xi]*K23[zk])*Bz
                              +Jiv32[yj][xi]*(Bx*K23[zk] - By*K13[zk]));

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //dF/du_{O+,r,i,j,k}
                    vals[nv]=Bdve[0]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //dF/du_{O+,theta,i,j,k}
                    vals[nv]=Bdve[1]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //dF/du_{O+,phi,i,j,k}
                    vals[nv]=Bdve[2]*temp;
                    nv++;

                    temp=ns[1]/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10;  //dF/du_{H+,r,i,j,k}
                    vals[nv]=Bdve[0]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11;  //dF/du_{H+,theta,i,j,k}
                    vals[nv]=Bdve[1]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12;  //dF/du_{H+,phi,i,j,k}
                    vals[nv]=Bdve[2]*temp;
                    nv++;

                    temp=ns[2]/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13;  //dF/du_{He+,r,i,j,k}
                    vals[nv]=Bdve[0]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14;  //dF/du_{He+,theta,i,j,k}
                    vals[nv]=Bdve[1]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15;  //dF/du_{He+,phi,i,j,k}
                    vals[nv]=Bdve[2]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=31;  //dF/dB^r_{i,j,k}
                    vals[nv]=-( Jiv12[zk][yj][xi]*(Jiv21[zk][yj]*uez - Jiv31[yj]*uey)
                               +Jiv22[zk][yj][xi]*(Jiv31[yj]*uex - Jiv11[zk][yj]*uez)
                               +Jiv32[yj][xi]*(Jiv11[zk][yj]*uey - Jiv21[zk][yj]*uex))/r2sintheta[yj][xi];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=32;  //dF/dB^theta_{i,j,k}
                    vals[nv]=-( Jiv12[zk][yj][xi]*(Jiv22[zk][yj][xi]*uez - Jiv32[yj][xi] *uey)
                               +Jiv22[zk][yj][xi]*(Jiv32[yj][xi]*uex - Jiv12[zk][yj][xi]*uez)
                               +Jiv32[yj][xi]*(Jiv12[zk][yj][xi]*uey - Jiv22[zk][yj][xi]*uex))/r2sintheta[yj][xi];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=33;  //dF/dB^phi_{i,j,k}
                    vals[nv]=-( (Jiv12[zk][yj][xi]*Jiv23[zk][yj][xi] - Jiv22[zk][yj][xi]*Jiv13[zk][yj][xi])*uez
                               +Jiv32[yj][xi] *(Jiv13[zk][yj][xi]*uey - Jiv23[zk][yj][xi]*uex))/r2sintheta[yj][xi];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=35;  //dF/dE^theta_{i,j,k}
                    vals[nv]=1.0;
                    nv++;
                }
                else if (ir == 36) { //for E^phi equation
                    temp=nsmore/ne;
                    Bdve[0]=-( Jiv13[zk][yj][xi]*(By*K31[yj] - Bz*K21[zk][yj])
                              +Jiv23[zk][yj][xi]*(Bz*K11[zk][yj] - Bx*K31[yj]));
                    Bdve[1]=-( Jiv13[zk][yj][xi]*(By*K32[yj] - Bz*K22[zk][yj])
                              +Jiv23[zk][yj][xi]*(Bz*K12[zk][yj] - Bx*K32[yj]));
                    Bdve[2]= (Jiv13[zk][yj][xi]*K23[zk] - Jiv23[zk][yj][xi]*K13[zk])*Bz;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //dF/du_{O+,r,i,j,k}
                    vals[nv]=Bdve[0]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //dF/du_{O+,theta,i,j,k}
                    vals[nv]=Bdve[1]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //dF/du_{O+,phi,i,j,k}
                    vals[nv]=Bdve[2]*temp;
                    nv++;

                    temp=ns[1]/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10;  //dF/du_{H+,r,i,j,k}
                    vals[nv]=Bdve[0]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11;  //dF/du_{H+,theta,i,j,k}
                    vals[nv]=Bdve[1]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12;  //dF/du_{H+,phi,i,j,k}
                    vals[nv]=Bdve[2]*temp;
                    nv++;

                    temp=ns[2]/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13;  //dF/du_{He+,r,i,j,k}
                    vals[nv]=Bdve[0]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14;  //dF/du_{He+,theta,i,j,k}
                    vals[nv]=Bdve[1]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15;  //dF/du_{He+,phi,i,j,k}
                    vals[nv]=Bdve[2]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=31;  //dF/dB^r_{i,j,k}
                    vals[nv]=-( Jiv13[zk][yj][xi]*(Jiv21[zk][yj]*uez - Jiv31[yj]*uey)
                               +Jiv23[zk][yj][xi]*(Jiv31[yj]*uex - Jiv11[zk][yj]*uez))/r2sintheta[yj][xi];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=32;  //dF/dB^theta_{i,j,k}
                    vals[nv]=-( Jiv13[zk][yj][xi]*(Jiv22[zk][yj][xi]*uez - Jiv32[yj][xi]*uey)
                               +Jiv23[zk][yj][xi]*(Jiv32[yj][xi]*uex - Jiv12[zk][yj][xi]*uez))/r2sintheta[yj][xi];
                    nv++;

                    //col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=33;  //dF/dB^phi_{i,j,k}
                    //vals[nv]=-(Jiv13[zk][yj][xi]*Jiv23[zk][yj][xi]-Jiv23[zk][yj][xi]*Jiv13[zk][yj][xi])
                    //          *uez/r2sintheta[yj][xi];
                    //nv++;

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

    DMRestoreLocalVector(da, &localX);
    DMRestoreLocalVector(da, &localU);

    delete[] vals;
    delete[] col;

    return ierr;
}
