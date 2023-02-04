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
    int        i, j, k, ir, nv, xs, xm, ys, ym, zs, zm, xi, yj, zk;
    MatStencil row, *col;
    PetscErrorCode ierr=0;
    double     *vals;
    const int  nzer=a4;

    double     temp, ne, rhon, Nn, ns[7], nsmore, nsmore_norm; //nn[7], 
    double     Bx, By, Bz, Br, Bt, Bp, uex, uey, uez;

    const double four3rd=4.0/3.0; //two3rd=2.0/3.0, dT=1.0e-8;

    int    s, s14;
    double sum_nusq[7], sum_nueq, sum_nues_div_ms, sum_nueq_div_mq, nusq;
    double sum_nuHiOi, sum_nuHeiOi, mq_nusq_msmq[7], mi_nust_mims[2], nust_mims[2];
    double nust_msmt, mt_nust_msmt, uir, uit, uip, rhos[7];
    double temp1, uir_unr, uit_unt, uip_unp, rhossum_nusq[3]; //, Ce, Cep, Te, Tn,
    double sum_nusq_msmq[7], neme_nsms_nues, me_nsms_nues;
    //double Lorentz_r_Oi_Oi, Lorentz_t_Oi_Oi, Lorentz_p_Oi_Oi;
    //double Lorentz_r_Hi_Oi, Lorentz_t_Hi_Oi, Lorentz_p_Hi_Oi;
    //double Lorentz_r_Hei_Oi, Lorentz_t_Hei_Oi, Lorentz_p_Hei_Oi;
    //double unr, unt, unp, usr_uir[3], ust_uit[3], usp_uip[3], fn, us_un_sq[3];
    //double ExdnsOi, EydnsOi, EzdnsOi, ExdnsHi, EydnsHi, EzdnsHi, ExdnsHei, EydnsHei, EzdnsHei;
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

            //ne=0.0; Nn=0.0; rhon=0.0; 
            sum_nueq_div_mq=0.0; sum_nueq=0.0;
            for (s = 0; s < sl; s++) {
                ns[s]=uu[k][j][i].fx[s+27]; rhos[s]=ns[s]*ams[s];
                //nn[s]=exp(xx[k][j][i].fx[20+s]); Nn =+nn[s]; rhon += nn[s]*ams[s];

                sum_nueq += nust[zk][yj][xi][7+s];
                sum_nueq_div_mq += nust[zk][yj][xi][7+s]/ams[s];

                s14=14*s; sum_nusq[s]=0.0; sum_nusq_msmq[s]=0.0; mq_nusq_msmq[s]=0.0;
                for (int t = 0; t < sl; t++) {
                    nusq=nust[zk][yj][xi][21+s14+t];
                    sum_nusq[s] += nusq;

                    temp=nusq/(ams[s]+ams[t]);
                    sum_nusq_msmq[s] += temp; mq_nusq_msmq[s] += ams[t]*temp;
                }
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

            //usr_uir[0]=xx[k][j][i].fx[7]-uir;  ust_uit[0]=xx[k][j][i].fx[8]-uit;  usp_uip[0]=xx[k][j][i].fx[9]-uip;
            //usr_uir[1]=xx[k][j][i].fx[10]-uir; ust_uit[1]=xx[k][j][i].fx[11]-uit; usp_uip[1]=xx[k][j][i].fx[12]-uip;
            //usr_uir[2]=xx[k][j][i].fx[13]-uir; ust_uit[2]=xx[k][j][i].fx[14]-uit; usp_uip[2]=xx[k][j][i].fx[15]-uip;

            /*ExdnsOi=( By*(Kmat.K31[yj]*usr_uir[0]+Kmat.K32[yj]*ust_uit[0]-uez/ne)
                     -Bz*(Kmat.K21[kj]*usr_uir[0]+Kmat.K22[kj]*ust_uit[0]+Kmat.K23[zk]*usp_uip[0])-uey/ne)/ner2sintheta;
            EydnsOi=( Bz*(Kmat.K11[kj]*usr_uir[0]+Kmat.K12[kj]*ust_uit[0]+Kmat.K13[zk]*usp_uip[0]-uex/ne)
                     -Bx*(Kmat.K31[yj]*usr_uir[0]+Kmat.K32[yj]*ust_uit[0])-uez/ne)/ner2sintheta;
            EzdnsOi=( Bx*(Kmat.K21[kj]*usr_uir[0]+Kmat.K22[kj]*ust_uit[0]+Kmat.K23[zk]*usp_uip[0]-uey/ne)
                     -By*(Kmat.K11[kj]*usr_uir[0]+Kmat.K12[kj]*ust_uit[0]+Kmat.K13[zk]*usp_uip[0])-uex/ne)/ner2sintheta;

            ExdnsHi=( By*(Kmat.K31[yj]*usr_uir[1]+Kmat.K32[yj]*ust_uit[1]-uez/ne)
                     -Bz*(Kmat.K21[kj]*usr_uir[1]+Kmat.K22[kj]*ust_uit[1]+Kmat.K23[zk]*usp_uip[1])-uey/ne)/ner2sintheta;
            EydnsHi=( Bz*(Kmat.K11[kj]*usr_uir[1]+Kmat.K12[kj]*ust_uit[1]+Kmat.K13[zk]*usp_uip[1]-uex/ne)
                     -Bx*(Kmat.K31[yj]*usr_uir[1]+Kmat.K32[yj]*ust_uit[1])-uez/ne)/ner2sintheta;
            EzdnsHi=( Bx*(Kmat.K21[kj]*usr_uir[1]+Kmat.K22[kj]*ust_uit[1]+Kmat.K23[zk]*usp_uip[1]-uey/ne)
                     -By*(Kmat.K11[kj]*usr_uir[1]+Kmat.K12[kj]*ust_uit[1]+Kmat.K13[zk]*usp_uip[1])-uex/ne)/ner2sintheta;

            ExdnsHei=( By*(Kmat.K31[yj]*usr_uir[2]+Kmat.K32[yj]*ust_uit[2]-uez/ne)
                      -Bz*(Kmat.K21[kj]*usr_uir[2]+Kmat.K22[kj]*ust_uit[2]+Kmat.K23[zk]*usp_uip[2])-uey/ne)/ner2sintheta;
            EydnsHei=( Bz*(Kmat.K11[kj]*usr_uir[2]+Kmat.K12[kj]*ust_uit[2]+Kmat.K13[zk]*usp_uip[2]-uex/ne)
                      -Bx*(Kmat.K31[yj]*usr_uir[2]+Kmat.K32[yj]*ust_uit[2])-uez/ne)/ner2sintheta;
            EzdnsHei=( Bx*(Kmat.K21[kj]*usr_uir[2]+Kmat.K22[kj]*ust_uit[2]+Kmat.K23[zk]*usp_uip[2]-uey/ne)
                      -By*(Kmat.K11[kj]*usr_uir[2]+Kmat.K12[kj]*ust_uit[2]+Kmat.K13[zk]*usp_uip[2])-uex/ne)/ner2sintheta;*/

            /*Br=(Kmat.K11[kj]*Bx+Kmat.K21[kj]*By+Kmat.K31[yj]*Bz)/r2sintheta[yj][xi]+uu[k][j][i].fx[0];
            Bt=(Kmat.K12[kj]*Bx+Kmat.K22[kj]*By+Kmat.K32[yj]*Bz)/r2sintheta[yj][xi]+uu[k][j][i].fx[1];
            Bp=(Kmat.K13[zk]*Bx+Kmat.K23[zk]*By)/r2sintheta[yj][xi]+uu[k][j][i].fx[2];*/

            Br=uu[k][j][i].fx[0]; Bt=uu[k][j][i].fx[1]; Bp=uu[k][j][i].fx[2];

            sum_nuHiOi= nust[zk][yj][xi][29]+nust[zk][yj][xi][31]+nust[zk][yj][xi][32]
                       +nust[zk][yj][xi][33]+nust[zk][yj][xi][34];
            sum_nuHeiOi= nust[zk][yj][xi][43]+nust[zk][yj][xi][45]+nust[zk][yj][xi][46]
                        +nust[zk][yj][xi][47]+nust[zk][yj][xi][48];

            /*Lorentz_r_Oi_Oi = ((xx[k][j][i].fx[8]-uit)*Bp-(xx[k][j][i].fx[9]-uip)*Bt)/ne;
            Lorentz_t_Oi_Oi = ((xx[k][j][i].fx[9]-uip)*Br-(xx[k][j][i].fx[7]-uir)*Bp)/ne;
            Lorentz_p_Oi_Oi = ((xx[k][j][i].fx[7]-uit)*Bt-(xx[k][j][i].fx[8]-uit)*Br)/ne;

            Lorentz_r_Hi_Oi = ((xx[k][j][i].fx[11]-uit)*Bp-(xx[k][j][i].fx[12]-uip)*Bt)/ne;
            Lorentz_t_Hi_Oi = ((xx[k][j][i].fx[12]-uip)*Br-(xx[k][j][i].fx[10]-uir)*Bp)/ne;
            Lorentz_p_Hi_Oi = ((xx[k][j][i].fx[10]-uit)*Bt-(xx[k][j][i].fx[11]-uit)*Br)/ne;

            Lorentz_r_Hei_Oi = ((xx[k][j][i].fx[14]-uit)*Bp-(xx[k][j][i].fx[15]-uip)*Bt)/ne;
            Lorentz_t_Hei_Oi = ((xx[k][j][i].fx[15]-uip)*Br-(xx[k][j][i].fx[13]-uir)*Bp)/ne;
            Lorentz_p_Hei_Oi = ((xx[k][j][i].fx[13]-uit)*Bt-(xx[k][j][i].fx[14]-uit)*Br)/ne;*/

            //unr=xx[k][j][i].fx[27]; unt=xx[k][j][i].fx[28]; unp=xx[k][j][i].fx[29];

            rhossum_nusq[0]=( rhos[0]*sum_nusq[0]+rhos[3]*sum_nusq[3]+rhos[4]*sum_nusq[4]
                             +rhos[5]*sum_nusq[5]+rhos[6]*sum_nusq[6])/rhon;
            rhossum_nusq[1]=rhos[1]*sum_nusq[1]/rhon;
            rhossum_nusq[2]=rhos[2]*sum_nusq[2]/rhon;

            for (ir = 0; ir < a4; ir++) {
                row.c=ir; nv=0;

                if (ir < sl) {
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ir;
                    vals[nv]=a; nv++;
                }
                else if (ir == 7) {
                    /*temp = qms[0]*Lorentz_r_Oi_Oi;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=0;  //dF/d(xx[k][j][i].fx[0])
                    vals[nv]= ns[0]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=1;  //dF/d(xx[k][j][i].fx[1])
                    vals[nv]= qms[0]*ns[1]*Lorentz_r_Hi_Oi;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=2;  //dF/d(xx[k][j][i].fx[2])
                    vals[nv]= qms[0]*ns[2]*Lorentz_r_Hei_Oi;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=3;
                    vals[nv]= ns[3]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=4;
                    vals[nv]= ns[4]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=5;
                    vals[nv]= ns[5]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=6;
                    vals[nv]= ns[6]*temp;
                    nv++;*/

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;
                    vals[nv]= a+sum_nusq[0]+nust[zk][yj][xi][15]+nust[zk][yj][xi][16];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;
                    vals[nv]= qms[0]*nsmore_norm*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;
                    vals[nv]=-qms[0]*nsmore_norm*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10;
                    vals[nv]=-nust[zk][yj][xi][15];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11;
                    vals[nv]= qms[0]*ns[1]/ne*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12;
                    vals[nv]=-qms[0]*ns[1]/ne*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13;
                    vals[nv]=-nust[zk][yj][xi][16];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14;
                    vals[nv]= qms[0]*ns[2]/ne*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15;
                    vals[nv]=-qms[0]*ns[2]/ne*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=27;
                    vals[nv]=-sum_nusq[0];
                    nv++;
                }
                else if (ir == 8) {
                    /*temp = qms[0]*Lorentz_t_Oi_Oi;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=0;
                    vals[nv]= ns[0]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=1;
                    vals[nv]= qms[0]*ns[1]*Lorentz_t_Hi_Oi;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=2;
                    vals[nv]= qms[0]*ns[2]*Lorentz_t_Hei_Oi;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=3;
                    vals[nv]= ns[3]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=4;
                    vals[nv]= ns[4]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=5;
                    vals[nv]= ns[5]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=6;
                    vals[nv]= ns[6]*temp;
                    nv++;*/

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;
                    vals[nv]=-qms[0]*nsmore_norm*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;
                    vals[nv]= a+sum_nusq[0]+nust[zk][yj][xi][15]+nust[zk][yj][xi][16];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;
                    vals[nv]= qms[0]*nsmore_norm*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10;
                    vals[nv]=-qms[0]*ns[1]/ne*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11;
                    vals[nv]=-nust[zk][yj][xi][15];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12;
                    vals[nv]= qms[0]*ns[1]/ne*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13;
                    vals[nv]=-qms[0]*ns[2]/ne*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14;
                    vals[nv]=-nust[zk][yj][xi][16];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15;
                    vals[nv]= qms[0]*ns[2]/ne*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=28;
                    vals[nv]=-sum_nusq[0];
                    nv++;
                }
                else if (ir == 9) {
                    /*temp = qms[0]*Lorentz_p_Oi_Oi;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=0;
                    vals[nv]= ns[0]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=1;
                    vals[nv]= qms[0]*ns[1]*Lorentz_p_Hi_Oi;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=2;
                    vals[nv]= qms[0]*ns[2]*Lorentz_p_Hei_Oi;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=3;
                    vals[nv]= ns[3]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=4;
                    vals[nv]= ns[4]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=5;
                    vals[nv]= ns[5]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=6;
                    vals[nv]= ns[6]*temp;
                    nv++;*/

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;
                    vals[nv]= qms[0]*nsmore_norm*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;
                    vals[nv]=-qms[0]*nsmore_norm*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;
                    vals[nv]= a+sum_nusq[0]+nust[zk][yj][xi][15]+nust[zk][yj][xi][16];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10;
                    vals[nv]= qms[0]*ns[1]/ne*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11;
                    vals[nv]=-qms[0]*ns[1]/ne*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12;
                    vals[nv]=-nust[zk][yj][xi][15];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13;
                    vals[nv]= qms[0]*ns[2]/ne*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14;
                    vals[nv]=-qms[0]*ns[2]/ne*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15;
                    vals[nv]=-nust[zk][yj][xi][16];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=29;
                    vals[nv]=-sum_nusq[0];
                    nv++;
                }
                else if (ir == 10) {
                    /*temp = qms[1]*Lorentz_r_Oi_Oi;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=0;
                    vals[nv]= ns[0]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=1;
                    vals[nv]= qms[1]*ns[1]*Lorentz_r_Hi_Oi;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=2;
                    vals[nv]= qms[1]*ns[2]*Lorentz_r_Hei_Oi;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=3;
                    vals[nv]= ns[3]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=4;
                    vals[nv]= ns[4]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=5;
                    vals[nv]= ns[5]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=6;
                    vals[nv]= ns[6]*temp;
                    nv++;*/

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //u_{O+,r}
                    vals[nv]=-sum_nuHiOi;
                    nv++;

                    temp=qms[1]*nsmore/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //u_{O+,theta}
                    vals[nv]= temp*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //u_{O+,phi}
                    vals[nv]=-temp*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10;  //u_{H+,r}
                    vals[nv]=a+sum_nusq[1]+sum_nuHiOi+nust[zk][yj][xi][30];
                    nv++;

                    temp=qms[1]*(ns[1]/ne-1.0);
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11; //u_{H+,theta}
                    vals[nv]= temp*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12; //u_{H+,phi}
                    vals[nv]=-temp*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13;
                    vals[nv]=-nust[zk][yj][xi][30];
                    nv++;

                    temp=qms[1]*ns[2]/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14;
                    vals[nv]= temp*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15;
                    vals[nv]=-temp*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=27;
                    vals[nv]=-sum_nusq[1];
                    nv++;
                }
                else if (ir == 11) {
                    /*temp = qms[1]*Lorentz_t_Oi_Oi;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=0;
                    vals[nv]= ns[0]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=1;
                    vals[nv]= qms[1]*ns[1]*Lorentz_t_Hi_Oi;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=2;
                    vals[nv]= qms[1]*ns[2]*Lorentz_t_Hei_Oi;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=3;
                    vals[nv]= ns[3]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=4;
                    vals[nv]= ns[4]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=5;
                    vals[nv]= ns[5]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=6;
                    vals[nv]= ns[6]*temp;
                    nv++;*/

                    temp=qms[1]*nsmore/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //u_{O+,r}
                    vals[nv]=-temp*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //u_{O+,theta}
                    vals[nv]=-sum_nuHiOi;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //u_{O+,phi}
                    vals[nv]= temp*Br;
                    nv++;

                    temp=qms[1]*(ns[1]/ne-1.0);
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10; //u_{H+,r}
                    vals[nv]=-temp*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11;  //u_{H+,theta}
                    vals[nv]=a+sum_nusq[1]+sum_nuHiOi+nust[zk][yj][xi][30];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12; //u_{H+,phi}
                    vals[nv]= temp*Br;
                    nv++;

                    temp=qms[1]*ns[2]/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13;
                    vals[nv]=-temp*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14;
                    vals[nv]=-nust[zk][yj][xi][30];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15;
                    vals[nv]= temp*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=28;
                    vals[nv]=-sum_nusq[1];
                    nv++;
                }
                else if (ir == 12) {
                    /*temp = qms[1]*Lorentz_p_Oi_Oi;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=0;
                    vals[nv]= ns[0]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=1;
                    vals[nv]= qms[1]*ns[1]*Lorentz_p_Hi_Oi;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=2;
                    vals[nv]= qms[1]*ns[2]*Lorentz_p_Hei_Oi;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=3;
                    vals[nv]= ns[3]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=4;
                    vals[nv]= ns[4]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=5;
                    vals[nv]= ns[5]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=6;
                    vals[nv]= ns[6]*temp;
                    nv++;*/

                    temp=qms[1]*nsmore/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //u_{O+,r}
                    vals[nv]= temp*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //u_{O+,theta}
                    vals[nv]=-temp*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //u_{O+,phi}
                    vals[nv]=-sum_nuHiOi;
                    nv++;

                    temp=qms[1]*(ns[1]/ne-1.0);
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10; //u_{H+,r}
                    vals[nv]= temp*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11; //u_{H+,theta}
                    vals[nv]=-temp*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12;  //u_{H+,phi}
                    vals[nv]=a+sum_nusq[1]+sum_nuHiOi+nust[zk][yj][xi][30];
                    nv++;

                    temp=qms[1]*ns[2]/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13;
                    vals[nv]= temp*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14;
                    vals[nv]=-temp*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15;
                    vals[nv]=-nust[zk][yj][xi][30];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=29;
                    vals[nv]=-sum_nusq[1];
                    nv++;
                }
                else if (ir == 13) {
                    /*temp = qms[2]*Lorentz_r_Oi_Oi;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=0;
                    vals[nv]= ns[0]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=1;
                    vals[nv]= qms[2]*ns[1]*Lorentz_r_Hi_Oi;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=2;
                    vals[nv]= qms[2]*ns[2]*Lorentz_r_Hei_Oi;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=3;
                    vals[nv]= ns[3]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=4;
                    vals[nv]= ns[4]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=5;
                    vals[nv]= ns[5]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=6;
                    vals[nv]= ns[6]*temp;
                    nv++;*/

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //u_{O+,r}
                    vals[nv]=-sum_nuHeiOi;
                    nv++;

                    temp=qms[2]*nsmore/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //u_{O+,theta}
                    vals[nv]= temp*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //u_{O+,phi}
                    vals[nv]=-temp*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10; //u_{H+,r}
                    vals[nv]=-nust[zk][yj][xi][44];
                    nv++;

                    temp=qms[2]*ns[1]/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11; //u_{H+,theta}
                    vals[nv]= temp*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12; //u_{H+,phi}
                    vals[nv]=-temp*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13;  //u_{He+,r}
                    vals[nv]=a+sum_nusq[2]+sum_nuHeiOi+nust[zk][yj][xi][44];
                    nv++;

                    temp=qms[2]*(ns[2]/ne-1.0);
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14; //u_{He+,theta}
                    vals[nv]= temp*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15; //u_{He+,phi}
                    vals[nv]=-temp*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=27;
                    vals[nv]=-sum_nusq[2];
                    nv++;
                }
                else if (ir == 14) {
                    /*temp = qms[2]*Lorentz_t_Oi_Oi;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=0;
                    vals[nv]= ns[0]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=1;
                    vals[nv]= qms[2]*ns[1]*Lorentz_t_Hi_Oi;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=2;
                    vals[nv]= qms[2]*ns[2]*Lorentz_t_Hei_Oi;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=3;
                    vals[nv]= ns[3]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=4;
                    vals[nv]= ns[4]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=5;
                    vals[nv]= ns[5]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=6;
                    vals[nv]= ns[6]*temp;
                    nv++;*/

                    temp=qms[2]*nsmore/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //u_{O+,r}
                    vals[nv]=-temp*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //u_{O+,theta}
                    vals[nv]=-sum_nuHeiOi;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //u_{O+,phi}
                    vals[nv]= temp*Br;
                    nv++;

                    temp=qms[2]*ns[1]/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10; //u_{H+,r}
                    vals[nv]=-temp*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11; //u_{H+,theta}
                    vals[nv]=-nust[zk][yj][xi][44];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12; //u_{H+,phi}
                    vals[nv]= temp*Br;
                    nv++;

                    temp=qms[2]*(ns[2]/ne-1.0);
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13; //u_{He+,theta}
                    vals[nv]=-temp*Bp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14;  //u_{He+,theta}
                    vals[nv]=a+sum_nusq[2]+sum_nuHeiOi+nust[zk][yj][xi][44];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15; //u_{He+,phi}
                    vals[nv]= temp*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=28;
                    vals[nv]=-sum_nusq[2];
                    nv++;
                }
                else if (ir == 15) {
                    /*temp = qms[2]*Lorentz_p_Oi_Oi;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=0;
                    vals[nv]= ns[0]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=1;
                    vals[nv]= qms[2]*ns[1]*Lorentz_p_Hi_Oi;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=2;
                    vals[nv]= qms[2]*ns[2]*Lorentz_p_Hei_Oi;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=3;
                    vals[nv]= ns[3]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=4;
                    vals[nv]= ns[4]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=5;
                    vals[nv]= ns[5]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=6;
                    vals[nv]= ns[6]*temp;
                    nv++;*/

                    temp=qms[2]*nsmore/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //u_{O+,r}
                    vals[nv]= temp*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //u_{O+,theta}
                    vals[nv]=-temp*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //u_{O+,phi}
                    vals[nv]=-sum_nuHeiOi;
                    nv++;

                    temp=qms[2]*ns[1]/ne;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10; //u_{H+,r}
                    vals[nv]= temp*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11; //u_{H+,theta}
                    vals[nv]=-temp*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12; //u_{H+,phi}
                    vals[nv]=-nust[zk][yj][xi][44];
                    nv++;

                    temp=qms[2]*(ns[2]/ne-1.0);
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=13; //u_{He+,theta}
                    vals[nv]= temp*Bt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14; //u_{He+,theta}
                    vals[nv]=-temp*Br;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15;  //u_{He+,phi}
                    vals[nv]=a+sum_nusq[2]+sum_nuHeiOi+nust[zk][yj][xi][44];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=29;
                    vals[nv]=-sum_nusq[2];
                    nv++;
                }
                else if (ir == 16) {  //O+ temperature equation
                    me_nsms_nues=2.0*ame/(ns[0]*ams[0])*nust[zk][yj][xi][0];

                   /* temp=me_nsms_nues*(xx[k][j][i].fx[16]-xx[k][j][i].fx[19]);
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=0;  //U_{O+,r}
                    vals[nv]= ns[0]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=1;  //U_{O+,r}
                    vals[nv]= ns[1]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=2;  //U_{O+,r}
                    vals[nv]= ns[2]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=3;  //U_{O+,r}
                    vals[nv]= ns[3]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=4;  //U_{O+,r}
                    vals[nv]= ns[4]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=5;  //U_{O+,r}
                    vals[nv]= ns[5]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=6;  //U_{O+,r}
                    vals[nv]= ns[6]*temp;
                    nv++;*/

                    mi_nust_mims[0]=ams[1]*nust[zk][yj][xi][15]/(ams[0]+ams[1]);
                    mi_nust_mims[1]=ams[2]*nust[zk][yj][xi][16]/(ams[0]+ams[2]);

                    temp=four3rd*ams[0];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //U_{O+,r}
                    vals[nv]=-temp*( mq_nusq_msmq[0]*(xx[k][j][i].fx[7]-xx[k][j][i].fx[27])
                                    +mi_nust_mims[0]*(xx[k][j][i].fx[7]-xx[k][j][i].fx[10])
                                    +mi_nust_mims[1]*(xx[k][j][i].fx[7]-xx[k][j][i].fx[13]));
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //U_{O+,theta}
                    vals[nv]=-temp*( mq_nusq_msmq[0]*(xx[k][j][i].fx[8]-xx[k][j][i].fx[28])
                                    +mi_nust_mims[0]*(xx[k][j][i].fx[8]-xx[k][j][i].fx[11])
                                    +mi_nust_mims[1]*(xx[k][j][i].fx[8]-xx[k][j][i].fx[14]));
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //U_{O+,phi}
                    vals[nv]=-temp*( mq_nusq_msmq[0]*(xx[k][j][i].fx[9]-xx[k][j][i].fx[29])
                                    +mi_nust_mims[0]*(xx[k][j][i].fx[9]-xx[k][j][i].fx[12])
                                    +mi_nust_mims[1]*(xx[k][j][i].fx[9]-xx[k][j][i].fx[15]));
                    nv++;

                    temp1=temp*mi_nust_mims[0];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10;  //U_{H+,r}
                    vals[nv]= temp1*(xx[k][j][i].fx[7]-xx[k][j][i].fx[10]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11;  //U_{H+,theta}
                    vals[nv]= temp1*(xx[k][j][i].fx[8]-xx[k][j][i].fx[11]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12;  //U_{He+,phi}
                    vals[nv]= temp1*(xx[k][j][i].fx[9]-xx[k][j][i].fx[12]);
                    nv++;

                    temp1=temp*mi_nust_mims[1];
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
                    neme_nsms_nues=ne*me_nsms_nues;
                    nust_mims[0]=nust[zk][yj][xi][15]/(ams[0]+ams[1]);
                    nust_mims[1]=nust[zk][yj][xi][16]/(ams[0]+ams[2]);

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=16;  //T_{O+}
                    vals[nv]=a+temp*(sum_nusq_msmq[0]+nust_mims[0]+nust_mims[1])+neme_nsms_nues;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=17;  //T_{H+}
                    vals[nv]=-temp*nust_mims[0];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=18;  //T_{He+}
                    vals[nv]=-temp*nust_mims[1];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=19;  //T_{e}
                    vals[nv]=-neme_nsms_nues;
                    nv++;

                    temp=four3rd*ams[0]*mq_nusq_msmq[0];
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
                    vals[nv]=-2.0*ams[0]*sum_nusq_msmq[0];
                    nv++;
                }
                else if (ir == 17) {  //H+ temperature equation
                    me_nsms_nues=2.0*ame/(ns[1]*ams[1])*nust[zk][yj][xi][1];

                    /*temp=me_nsms_nues*(xx[k][j][i].fx[17]-xx[k][j][i].fx[19]);
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=0;  //U_{O+,r}
                    vals[nv]= ns[0]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=1;  //U_{O+,r}
                    vals[nv]= ns[1]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=2;  //U_{O+,r}
                    vals[nv]= ns[2]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=3;  //U_{O+,r}
                    vals[nv]= ns[3]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=4;  //U_{O+,r}
                    vals[nv]= ns[4]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=5;  //U_{O+,r}
                    vals[nv]= ns[5]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=6;  //U_{O+,r}
                    vals[nv]= ns[6]*temp;
                    nv++;*/

                    mt_nust_msmt= ams[0]*nust[zk][yj][xi][29]/(ams[1]+ams[0])
                                 +ams[3]*nust[zk][yj][xi][31]/(ams[1]+ams[3])
                                 +ams[4]*nust[zk][yj][xi][32]/(ams[1]+ams[4])
                                 +ams[5]*nust[zk][yj][xi][33]/(ams[1]+ams[5])
                                 +ams[6]*nust[zk][yj][xi][34]/(ams[1]+ams[6]);

                    temp=four3rd*ams[1]*mt_nust_msmt;
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
                    mi_nust_mims[1]=ams[2]/(ams[1]+ams[2])*nust[zk][yj][xi][30];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10;  //dF/dU_{H+,r}
                    vals[nv]=-temp*( mq_nusq_msmq[1]*(xx[k][j][i].fx[10]-xx[k][j][i].fx[27])
                                    +mt_nust_msmt*(xx[k][j][i].fx[10]-xx[k][j][i].fx[7])
                                    +mi_nust_mims[1]*(xx[k][j][i].fx[10]-xx[k][j][i].fx[13]));
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11;  //dF/dU_{H+,theta}
                    vals[nv]=-temp*( mq_nusq_msmq[1]*(xx[k][j][i].fx[11]-xx[k][j][i].fx[28])
                                    +mt_nust_msmt*(xx[k][j][i].fx[11]-xx[k][j][i].fx[8])
                                    +mi_nust_mims[1]*(xx[k][j][i].fx[11]-xx[k][j][i].fx[14]));
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12;  //dF/dU_{H+,phi}
                    vals[nv]=-temp*( mq_nusq_msmq[1]*(xx[k][j][i].fx[12]-xx[k][j][i].fx[29])
                                    +mt_nust_msmt*(xx[k][j][i].fx[12]-xx[k][j][i].fx[9])
                                    +mi_nust_mims[1]*(xx[k][j][i].fx[12]-xx[k][j][i].fx[15]));
                    nv++;

                    temp=four3rd*ams[1]*mi_nust_mims[1];
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
                    nust_msmt= nust[zk][yj][xi][29]/(ams[1]+ams[0])+nust[zk][yj][xi][31]/(ams[1]+ams[3])
                              +nust[zk][yj][xi][32]/(ams[1]+ams[4])+nust[zk][yj][xi][33]/(ams[1]+ams[5])
                              +nust[zk][yj][xi][34]/(ams[1]+ams[6]);
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=16;  //dF/dT_{O+}
                    vals[nv]=-temp*nust_msmt;
                    nv++;

                    neme_nsms_nues=ne*me_nsms_nues;
                    nust_mims[1]=nust[zk][yj][xi][30]/(ams[1]+ams[2]);
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=17;  //dF/dT_{H+}
                    vals[nv]=a+temp*(sum_nusq_msmq[1]+nust_msmt+nust_mims[1])+neme_nsms_nues;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=18;  //dF/dT_{He+}
                    vals[nv]=-temp*nust_mims[1];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=19;  //dF/dT_{e}
                    vals[nv]=-neme_nsms_nues;
                    nv++;

                    temp=four3rd*ams[1]*mq_nusq_msmq[1];
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
                    vals[nv]=-2.0*ams[1]*sum_nusq_msmq[1];
                    nv++;
                }
                else if (ir == 18) {  //He+ temperature equation
                    me_nsms_nues=2.0*ame/(ns[2]*ams[2])*nust[zk][yj][xi][2];
                    
                    /*temp=me_nsms_nues*(xx[k][j][i].fx[18]-xx[k][j][i].fx[19]);
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=0;  //U_{O+,r}
                    vals[nv]= ns[0]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=1;  //U_{O+,r}
                    vals[nv]= ns[1]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=2;  //U_{O+,r}
                    vals[nv]= ns[2]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=3;  //U_{O+,r}
                    vals[nv]= ns[3]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=4;  //U_{O+,r}
                    vals[nv]= ns[4]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=5;  //U_{O+,r}
                    vals[nv]= ns[5]*temp;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=6;  //U_{O+,r}
                    vals[nv]= ns[6]*temp;
                    nv++;*/

                    mt_nust_msmt= ams[0]*nust[zk][yj][xi][43]/(ams[2]+ams[0])
                                 +ams[3]*nust[zk][yj][xi][45]/(ams[2]+ams[3])
                                 +ams[4]*nust[zk][yj][xi][46]/(ams[2]+ams[4])
                                 +ams[5]*nust[zk][yj][xi][47]/(ams[2]+ams[5])
                                 +ams[6]*nust[zk][yj][xi][48]/(ams[2]+ams[6]);

                    temp=four3rd*ams[2]*mt_nust_msmt;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //dF/dU_{O+,r}
                    vals[nv]= temp*(xx[k][j][i].fx[13]-xx[k][j][i].fx[7]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //dF/dU_{O+,theta}
                    vals[nv]= temp*(xx[k][j][i].fx[14]-xx[k][j][i].fx[8]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //dF/dU_{O+,phi}
                    vals[nv]= temp*(xx[k][j][i].fx[15]-xx[k][j][i].fx[9]);
                    nv++;

                    mi_nust_mims[0]=ams[1]*nust[zk][yj][xi][44]/(ams[2]+ams[1]);
                    temp=four3rd*ams[2]*mi_nust_mims[0];
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
                    vals[nv]=-temp*( mq_nusq_msmq[2]*(xx[k][j][i].fx[13]-xx[k][j][i].fx[27])
                                    +mt_nust_msmt*(xx[k][j][i].fx[13]-xx[k][j][i].fx[7])
                                    +mi_nust_mims[0]*(xx[k][j][i].fx[13]-xx[k][j][i].fx[10]));
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=14;  //dF/dU_{He+,theta}
                    vals[nv]=-temp*( mq_nusq_msmq[2]*(xx[k][j][i].fx[14]-xx[k][j][i].fx[28])
                                    +mt_nust_msmt*(xx[k][j][i].fx[14]-xx[k][j][i].fx[8])
                                    +mi_nust_mims[0]*(xx[k][j][i].fx[14]-xx[k][j][i].fx[11]));
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=15;  //dF/dU_{He+,r}
                    vals[nv]=-temp*( mq_nusq_msmq[2]*(xx[k][j][i].fx[15]-xx[k][j][i].fx[29])
                                    +mt_nust_msmt*(xx[k][j][i].fx[15]-xx[k][j][i].fx[9])
                                    +mi_nust_mims[0]*(xx[k][j][i].fx[15]-xx[k][j][i].fx[12]));
                    nv++;

                    temp=2.0*ams[2];
                    nust_msmt= nust[zk][yj][xi][43]/(ams[2]+ams[0])+nust[zk][yj][xi][45]/(ams[2]+ams[3])
                              +nust[zk][yj][xi][46]/(ams[2]+ams[4])+nust[zk][yj][xi][47]/(ams[2]+ams[5])
                              +nust[zk][yj][xi][48]/(ams[2]+ams[6]);
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=16;  //dF/dT_{O+}
                    vals[nv]=-temp*nust_msmt;
                    nv++;

                    nust_mims[1]=nust[zk][yj][xi][44]/(ams[2]+ams[1]);
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=17;  //dF/dT_{H+}
                    vals[nv]=-temp*nust_mims[1];
                    nv++;

                    neme_nsms_nues=ne*me_nsms_nues;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=18;  //dF/dT_{He+}
                    vals[nv]=a+temp*(sum_nusq_msmq[2]+nust_msmt+nust_mims[1])+neme_nsms_nues;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=19;  //dF/dT_{e}
                    vals[nv]=-neme_nsms_nues;
                    nv++;

                    temp=four3rd*ams[2]*mq_nusq_msmq[2];
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
                    vals[nv]=-2.0*ams[2]*sum_nusq_msmq[2];
                    nv++;
                }
                else if (ir == 19) {  //electron temperature equation
                    sum_nues_div_ms= nust[zk][yj][xi][0]/ams[0]+nust[zk][yj][xi][3]/ams[3]
                                    +nust[zk][yj][xi][4]/ams[4]+nust[zk][yj][xi][5]/ams[5]
                                    +nust[zk][yj][xi][6]/ams[6];

                    uir_unr=uir-xx[k][j][i].fx[27];
                    uit_unt=uit-xx[k][j][i].fx[28];
                    uip_unp=uip-xx[k][j][i].fx[29];

                    /*Te=xx[k][j][i].fx[19]; Tn=xx[k][j][i].fx[30];

                    Ce=ele_cooling_rate(xx, Te, Tn, ne, i, j, k);

                    temp1=four3rd*ame*sum_nueq/ne;
                    temp2=two3rd/(ne*ne)*Ce;

                    temp= temp1*( uir_unr*(xx[k][j][i].fx[7]-uir)+uit_unt*(xx[k][j][i].fx[8]-uit)
                                 +uip_unp*(xx[k][j][i].fx[9]-uip)) + temp2;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=0;  //dF/d(xx[k][j][i].fx[0])
                    vals[nv]=-temp*ns[0];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=1;  //dF/d(xx[k][j][i].fx[1])
                    vals[nv]=-(temp1*( uir_unr*(xx[k][j][i].fx[10]-uir)+uit_unt*(xx[k][j][i].fx[11]-uit)
                                      +uip_unp*(xx[k][j][i].fx[12]-uip)) - temp2)*ns[1];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=2;
                    vals[nv]=-(temp1*( uir_unr*(xx[k][j][i].fx[13]-uir)+uit_unt*(xx[k][j][i].fx[14]-uit)
                                      +uip_unp*(xx[k][j][i].fx[15]-uip)) - temp2)*ns[2];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=3;
                    vals[nv]=-temp*ns[3];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=4;
                    vals[nv]=-temp*ns[4];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=5;
                    vals[nv]=-temp*ns[5];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=6;
                    vals[nv]=-temp*ns[6];
                    nv++; */

                    temp1=four3rd*ame*sum_nueq;
                    temp=temp1*nsmore/ne;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //dF/dU_{O+,r}
                    vals[nv]=-temp*uir_unr;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //dF/dU_{O+,theta}
                    vals[nv]=-temp1*uit_unt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //dF/dU_{O+,phi}
                    vals[nv]= temp1*uip_unp;
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
                    vals[nv]=-temp*sum_nues_div_ms;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=17;  //dF/dT_{H+}
                    vals[nv]=-temp*nust[zk][yj][xi][1]/ams[1];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=18;  //dF/dT_{He+}
                    vals[nv]=-temp*nust[zk][yj][xi][2]/ams[2];
                    nv++;

                    //Te=xx[k][j][i].fx[19]+dT;
                    //Cep=ele_cooling_rate(xx, Te, Tn, ne, i, j, k);

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=19;  //dF/dT_{e}
                    vals[nv]= a+temp*( sum_nues_div_ms+nust[zk][yj][xi][1]/ams[1]+nust[zk][yj][xi][2]/ams[2]
                                      +sum_nueq_div_mq); //+two3rd*(Cep-Ce)/(ne*dT);
                    nv++;

                    temp=four3rd*ame*sum_nueq;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=27;  //dF/dU_{n,r}
                    vals[nv]= temp*uir_unr;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=28;  //dF/dU_{n,theta}
                    vals[nv]= temp*uit_unt;
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=29;  //dF/dU_{n,phi}
                    vals[nv]= temp*uip_unp;
                    nv++; 

                    //Tn=xx[k][j][i].fx[30]+dT;
                    //Cep=ele_cooling_rate(xx, Te, Tn, ne, i, j, k);

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=30;  //dF/dT_{n}
                    //vals[nv]=two3rd*(Cep-Ce)/(ne*dT)-2.0*ame*sum_nueq_div_mq;
                    vals[nv]=-2.0*ame*sum_nueq_div_mq;
                    nv++;
                }
                else if (ir >= 20 && ir <= 26) {
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ir;
                    vals[nv]=a; nv++;
                }
                else if (ir == 27) { //for u_{n,r} equation
                    /*temp=ams[0]/rhon*(xx[k][j][i].fx[27]-xx[k][j][i].fx[7]);

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=0;  //dF/d(xx[k][j][i].fx[0])
                    vals[nv]= temp*sum_nusq[0]*ns[0];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=1;  //dF/d(xx[k][j][i].fx[1])
                    vals[nv]= ams[1]*sum_nusq[1]/rhon*(xx[k][j][i].fx[27]-xx[k][j][i].fx[10])*ns[1];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=2;
                    vals[nv]= ams[2]*sum_nusq[2]/rhon*(xx[k][j][i].fx[27]-xx[k][j][i].fx[13])*ns[2];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=3;
                    vals[nv]= temp*sum_nusq[3]*ns[3];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=4;
                    vals[nv]= temp*sum_nusq[4]*ns[4];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=5;
                    vals[nv]= temp*sum_nusq[5]*ns[5];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=6;
                    vals[nv]= temp*sum_nusq[6]*ns[6];
                    nv++;*/
                    
                    /*col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=20;
                    vals[nv]= temp1*nn[0];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=21;
                    vals[nv]= temp1*nn[1];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=22;
                    vals[nv]= temp1*nn[2];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=23;
                    vals[nv]= temp1*nn[3];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=24;
                    vals[nv]= temp1*nn[4];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=25;
                    vals[nv]= temp1*nn[5];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=26;
                    vals[nv]= temp1*nn[6];
                    nv++;*/

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
                    /*temp=ams[0]/rhon*(xx[k][j][i].fx[28]-xx[k][j][i].fx[8]);

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=0;  //dF/d(xx[k][j][i].fx[0])
                    vals[nv]= temp*sum_nusq[0]*ns[0];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=1;  //dF/d(xx[k][j][i].fx[1])
                    vals[nv]= ams[1]*sum_nusq[1]/rhon*(xx[k][j][i].fx[28]-xx[k][j][i].fx[11])*ns[1];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=2;
                    vals[nv]= ams[2]*sum_nusq[2]/rhon*(xx[k][j][i].fx[28]-xx[k][j][i].fx[14])*ns[1];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=3;
                    vals[nv]= temp*sum_nusq[3]*ns[3];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=4;
                    vals[nv]= temp*sum_nusq[4]*ns[4];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=5;
                    vals[nv]= temp*sum_nusq[5]*ns[5];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=6;
                    vals[nv]= temp*sum_nusq[6]*ns[6];
                    nv++;*/
                    
                    /*col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=20;
                    vals[nv]= temp1*nn[0];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=21;
                    vals[nv]= temp1*nn[1];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=22;
                    vals[nv]= temp1*nn[2];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=23;
                    vals[nv]= temp1*nn[3];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=24;
                    vals[nv]= temp1*nn[4];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=25;
                    vals[nv]= temp1*nn[5];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=26;
                    vals[nv]= temp1*nn[6];
                    nv++;*/

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
                    /*temp=ams[0]/rhon*(xx[k][j][i].fx[29]-xx[k][j][i].fx[9]);

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=0;  //dF/d(xx[k][j][i].fx[0])
                    vals[nv]= temp*sum_nusq[0]*ns[0];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=1;  //dF/d(xx[k][j][i].fx[1])
                    vals[nv]= ams[1]*sum_nusq[1]/rhon*(xx[k][j][i].fx[29]-xx[k][j][i].fx[12])*ns[1];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=2;
                    vals[nv]= ams[2]*sum_nusq[2]/rhon*(xx[k][j][i].fx[29]-xx[k][j][i].fx[15])*ns[1];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=3;
                    vals[nv]= temp*sum_nusq[3]*ns[3];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=4;
                    vals[nv]= temp*sum_nusq[4]*ns[4];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=5;
                    vals[nv]= temp*sum_nusq[5]*ns[5];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=6;
                    vals[nv]= temp*sum_nusq[6]*ns[6];
                    nv++;*/

                    /*col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=20;
                    vals[nv]= temp1*nn[0];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=21;
                    vals[nv]= temp1*nn[1];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=22;
                    vals[nv]= temp1*nn[2];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=23;
                    vals[nv]= temp1*nn[3];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=24;
                    vals[nv]= temp1*nn[4];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=25;
                    vals[nv]= temp1*nn[5];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=26;
                    vals[nv]= temp1*nn[6];
                    nv++;*/

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
                        = four3rd/Nn*( rhos[0]*ams[0]*sum_nusq_msmq[0]+rhos[3]*ams[3]*sum_nusq_msmq[3]
                                     +rhos[4]*ams[4]*sum_nusq_msmq[4]+rhos[5]*ams[5]*sum_nusq_msmq[5]
                                     +rhos[6]*ams[6]*sum_nusq_msmq[6]);

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;  //dF/dU_{O+,r}
                    vals[nv]=-neu_friction_coef[0]*(xx[k][j][i].fx[7]-xx[k][j][i].fx[27]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;  //dF/dU_{O+,theta}
                    vals[nv]=-neu_friction_coef[0]*(xx[k][j][i].fx[8]-xx[k][j][i].fx[28]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;  //dF/dU_{O+,phi}
                    vals[nv]=-neu_friction_coef[0]*(xx[k][j][i].fx[9]-xx[k][j][i].fx[29]);
                    nv++;

                    neu_friction_coef[1]= four3rd/Nn*rhos[1]*ams[1]*sum_nusq_msmq[1];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10;  //dF/dU_{H+,r}
                    vals[nv]=-neu_friction_coef[1]*(xx[k][j][i].fx[10]-xx[k][j][i].fx[27]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11;  //dF/dU_{H+,theta}
                    vals[nv]=-neu_friction_coef[1]*(xx[k][j][i].fx[11]-xx[k][j][i].fx[28]);
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=12;  //dF/dU_{H+,phi}
                    vals[nv]=-neu_friction_coef[1]*(xx[k][j][i].fx[12]-xx[k][j][i].fx[29]);
                    nv++;

                    neu_friction_coef[2]=four3rd/Nn*rhos[2]*ams[2]*sum_nusq_msmq[2];
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
                        = 2.0/Nn*( rhos[0]*sum_nusq_msmq[0]+rhos[3]*sum_nusq_msmq[3]
                                  +rhos[4]*sum_nusq_msmq[4]+rhos[5]*sum_nusq_msmq[5]
                                  +rhos[6]*sum_nusq_msmq[6]);
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=16;  //dF/dT_{O+}}
                    vals[nv]=-neu_tem_exchange_coef[0];
                    nv++;

                    neu_tem_exchange_coef[1]= 2.0/Nn*ns[1]*ams[1]*sum_nusq_msmq[1];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=17;  //dF/dT_{H+}
                    vals[nv]=-neu_tem_exchange_coef[1];
                    nv++;

                    neu_tem_exchange_coef[2]= 2.0/Nn*ns[2]*ams[2]*sum_nusq_msmq[2];
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=18;  //dF/dT_{He+}
                    vals[nv]=-neu_tem_exchange_coef[2];
                    nv++;

                    neu_tem_exchange_coef[3]= 2.0*ne*ame/Nn*sum_nueq_div_mq;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=19;  //dF/dT_{e}
                    vals[nv]=-neu_tem_exchange_coef[3];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=30;  //dF/dT_{n}
                    vals[nv]= a+neu_tem_exchange_coef[0]+neu_tem_exchange_coef[1]+neu_tem_exchange_coef[2]
                               +neu_tem_exchange_coef[3];
                    nv++;

                    /*us_un_sq[0]= (xx[k][j][i].fx[7]-unr)*(xx[k][j][i].fx[7]-unr)
                                +(xx[k][j][i].fx[8]-unt)*(xx[k][j][i].fx[8]-unt)
                                +(xx[k][j][i].fx[9]-unp)*(xx[k][j][i].fx[9]-unp);
                    us_un_sq[1]= (xx[k][j][i].fx[10]-unr)*(xx[k][j][i].fx[10]-unr)
                                +(xx[k][j][i].fx[11]-unt)*(xx[k][j][i].fx[11]-unt)
                                +(xx[k][j][i].fx[12]-unp)*(xx[k][j][i].fx[12]-unp);
                    us_un_sq[2]= (xx[k][j][i].fx[13]-unr)*(xx[k][j][i].fx[13]-unr)
                                +(xx[k][j][i].fx[14]-unt)*(xx[k][j][i].fx[14]-unt)
                                +(xx[k][j][i].fx[15]-unp)*(xx[k][j][i].fx[15]-unp);*/

                    /*Tn=xx[k][j][i].fx[30];
                    fn= 0.5*( neu_friction_coef[0]*us_un_sq[0]+neu_friction_coef[1]*us_un_sq[1]
                             +neu_friction_coef[2]*us_un_sq[2])
                       -neu_tem_exchange_coef[0]*(Tn-xx[k][j][i].fx[16])
                       -neu_tem_exchange_coef[1]*(Tn-xx[k][j][i].fx[17])
                       -neu_tem_exchange_coef[2]*(Tn-xx[k][j][i].fx[18])
                       -neu_tem_exchange_coef[3]*(Tn-xx[k][j][i].fx[19]);

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=0;  //dF/dn_{O+}
                    vals[nv]=ns[0]/Nn*ams[0]*(2.0*(Tn-xx[k][j][i].fx[16])-two3rd*ams[0]*us_un_sq[0])*sum_nusq_msmq[0];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=1;  //dF/dn_{H+}
                    vals[nv]= ns[1]/Nn*ams[1]*(2.0*(Tn-xx[k][j][i].fx[17])-two3rd*ams[1]*us_un_sq[1])*sum_nusq_msmq[1];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=2;  //dF/dn_{He+}
                    vals[nv]= ns[2]/Nn*ams[2]*(2.0*(Tn-xx[k][j][i].fx[18])-two3rd*ams[2]*us_un_sq[2])*sum_nusq_msmq[2];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=3;  //dF/dn_{O2+}
                    vals[nv]= ns[3]/Nn*ams[3]*(2.0*(Tn-xx[k][j][i].fx[16])-two3rd*ams[3]*us_un_sq[0])*sum_nusq_msmq[3];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=4;  //dF/dn_{N2+}
                    vals[nv]= ns[4]/Nn*ams[4]*(2.0*(Tn-xx[k][j][i].fx[16])-two3rd*ams[4]*us_un_sq[0])*sum_nusq_msmq[4];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=5;  //dF/dn_{NO+}
                    vals[nv]= ns[5]/Nn*ams[5]*(2.0*(Tn-xx[k][j][i].fx[16])-two3rd*ams[5]*us_un_sq[0])*sum_nusq_msmq[5];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=6;  //dF/dn_{N+}
                    vals[nv]= ns[6]/Nn*ams[6]*(2.0*(Tn-xx[k][j][i].fx[16])-two3rd*ams[6]*us_un_sq[0])*sum_nusq_msmq[6];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=20;  //dF/dn_{O}
                    vals[nv]= fn*nn[0];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=21;  //dF/dn_{H}
                    vals[nv]= fn*nn[1];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=22;  //dF/dn_{He}
                    vals[nv]= fn*nn[2];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=23;  //dF/dn_{O2}
                    vals[nv]= fn*nn[3];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=24;  //dF/dn_{N2}
                    vals[nv]= fn*nn[4];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=25;  //dF/dn_{NO}
                    vals[nv]= fn*nn[5];
                    nv++;
 
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=26;  //dF/dn_{N}
                    vals[nv]= fn*nn[6];
                    nv++;*/
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
                    /*temp=-(Jinv.Jiv11[kj]*ExdnsOi+Jinv.Jiv21[kj]*EydnsOi+Jinv.Jiv31[yj]*EzdnsOi);

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=0;  //dF/dn_{O+,i,j,k}
                    vals[nv]=temp*ns[0];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=1;  //dF/dn_{H+,i,j,k}
                    vals[nv]=-(Jinv.Jiv11[kj]*ExdnsHi+Jinv.Jiv21[kj]*EydnsHi+Jinv.Jiv31[yj]*EzdnsHi)*ns[1];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=2;  //dF/dn_{H+,i,j,k}
                    vals[nv]=-(Jinv.Jiv11[kj]*ExdnsHei+Jinv.Jiv21[kj]*EydnsHei+Jinv.Jiv31[yj]*EzdnsHei)*ns[2];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=3;  //dF/dn_{O2+,i,j,k}
                    vals[nv]=temp*ns[3];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=4;  //dF/dn_{N2+,i,j,k}
                    vals[nv]=temp*ns[4];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=5;  //dF/dn_{NO+,i,j,k}
                    vals[nv]=temp*ns[5];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=6;  //dF/dn_{N+,i,j,k}
                    vals[nv]=temp*ns[6];
                    nv++;*/

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
                    /*temp=-(Jinv.Jiv12[kji]*ExdnsOi+Jinv.Jiv22[kji]*EydnsOi+Jinv.Jiv32[ji]*EzdnsOi);

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=0;  //dF/dn_{O+,i,j,k}
                    vals[nv]=temp*ns[0];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=1;  //dF/dn_{H+,i,j,k}
                    vals[nv]=-(Jinv.Jiv12[kji]*ExdnsHi+Jinv.Jiv22[kji]*EydnsHi+Jinv.Jiv32[ji]*EzdnsHi)*ns[1];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=2;  //dF/dn_{H+,i,j,k}
                    vals[nv]=-(Jinv.Jiv12[kji]*ExdnsHei+Jinv.Jiv22[kji]*EydnsHei+Jinv.Jiv32[ji]*EzdnsHei)*ns[2];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=3;  //dF/dn_{O2+,i,j,k}
                    vals[nv]=temp*ns[3];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=4;  //dF/dn_{N2+,i,j,k}
                    vals[nv]=temp*ns[4];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=5;  //dF/dn_{NO+,i,j,k}
                    vals[nv]=temp*ns[5];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=6;  //dF/dn_{N+,i,j,k}
                    vals[nv]=temp*ns[6];
                    nv++;*/

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
                    /*temp=-(Jinv.Jiv13[kji]*ExdnsOi+Jinv.Jiv23[kji]*EydnsOi);

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=0;  //dF/dn_{O+,i,j,k}
                    vals[nv]=temp*ns[0];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=1;  //dF/dn_{H+,i,j,k}
                    vals[nv]=-(Jinv.Jiv13[kji]*ExdnsHi+Jinv.Jiv23[kji]*EydnsHi)*ns[1];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=2;  //dF/dn_{H+,i,j,k}
                    vals[nv]=-(Jinv.Jiv13[kji]*ExdnsHei+Jinv.Jiv23[kji]*EydnsHei)*ns[2];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=3;  //dF/dn_{O2+,i,j,k}
                    vals[nv]=temp*ns[3];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=4;  //dF/dn_{N2+,i,j,k}
                    vals[nv]=temp*ns[4];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=5;  //dF/dn_{NO+,i,j,k}
                    vals[nv]=temp*ns[5];
                    nv++;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=6;  //dF/dn_{N+,i,j,k}
                    vals[nv]=temp*ns[6];
                    nv++;*/

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

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=33;  //dF/dB^phi_{i,j,k}
                    vals[nv]=-(Jiv13[zk][yj][xi]*Jiv23[zk][yj][xi]-Jiv23[zk][yj][xi]*Jiv13[zk][yj][xi])
                              *uez/r2sintheta[yj][xi];
                    nv++;

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
