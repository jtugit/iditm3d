#include <cmath>
#include <iostream>

using namespace std;

#include <petscts.h>
#include <petscdm.h>
#include <petscdmda.h>

#include "param.h"
#include "funcdef.h"

int jacobian(TS ts, double ftime, Vec X, Vec Xdt, double a, Mat Jac, Mat Jpre, void *ctx)
{
    DM         da;
    Vec        localX, localU, localZ;
    Field      ***xx, ***uu, ***zz;
    AppCtx     *params = (AppCtx*)ctx;
    PetscInt   i, j, k, ir, nv;
    PetscInt   xs, xm, ys, ym, zs, zm;
    MatStencil row, *col;
    PetscErrorCode ierr=0;
    double     *vals;
    const PetscInt  nzer=a4;

    PetscInt   xi, yj, zk;
    double     temp, ne, rhon, Nn;

    const double two3rd=2.0/3.0, four3rd=4.0/3.0, one6th=1.0/6.0;
    double dr2=dr*dr;

    int    im, ip, jm, jp, km, kp, s, t, kc, kcm, jcm, kcp, jcp;
    double rhos, rhoi, sum_nues, sum_nueq, sum_rhonusq, sum_nues_div_ms, sum_nueq_div_mq;
    double rhos_nusq_msmq, rhos_nusqmq_msmq, rhos_nusqms_msmq;
    double rhoe_sum_nueq; //ne = ni

    double a_imjk[2], a_imjk_multi_a_imjk_subtr, a_ipjk[2], a_ipjk_multi_a_ipjk_subtr;
    double b_ijmk[2], b_ijmk_multi_b_ijmk_subtr, b_ijpk[2], b_ijpk_multi_b_ijpk_subtr;
    double c_ijkm[2], c_ijkm_multi_c_ijkm_subtr, c_ijkp[2], c_ijkp_multi_c_ijkp_subtr;

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

    DMGetLocalVector(da, &localZ);
    DMGlobalToLocalBegin(da,params->Z,INSERT_VALUES,localZ);
    DMGlobalToLocalEnd(da,params->Z,INSERT_VALUES,localZ);
    DMDAVecGetArray(da, localZ, &zz);

    vals=new double[nzer];
    col=new MatStencil[nzer];

    for (k = zs; k < zs+zm; k++) {
      row.k=k; zk=k-zs; 

      km = k-1; kp = k+1;

      for (j = ys; j < ys+ym; j++) {
        row.j=j; yj=j-ys;

        if (j == 0) {
            for (i = xs; i< xs+xm; i++) {
                row.i = i;

                for (ir = 0; ir < nvar; ir++) {
                    row.c=ir; nv=0;

                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ir; vals[nv]=1.0;
                    nv++;

                    MatSetValuesStencil(Jac,1,&row,nv,col,vals,INSERT_VALUES);
                }
            }
            continue;
        }
        else if (j == Nth) {
            for (i = xs; i< xs+xm; i++) {
                row.i = i;

                for (ir = 0; ir < nvar; ir++) {
                    row.c=ir; nv=0;

                    if (ir != 24) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ir; vals[nv]=1.0;
                        nv++;
                    }
                    else {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ir; vals[nv]=a;
                        nv++;
                    }

                    MatSetValuesStencil(Jac,1,&row,nv,col,vals,INSERT_VALUES);
                }
            }
            continue;
        }

        jm =j-1; jp = j+1;

        if (j == 1) {kcm = (k+a3/2) % a3; jcm=1;}
        else {kcm = k; jcm = j;}

        if (j < Nthm) {kcp = k; jcp = j;}
        else {kcp = (k+a3/2) % a3; jcp = Nthm;}

        for (i = xs; i < xs+xm; i++) {
            row.i=i; xi=i-xs;

            if (i == 0 || i == Nr) {
                for (ir = 0; ir < nvar; ir++) {
                    row.c=ir; nv=0;
                    col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ir; vals[nv] = 1.0;
                    nv++;

                    MatSetValuesStencil(Jac,1,&row,nv,col,vals,INSERT_VALUES);
                }
            }
            else {
                im = i-1; ip = i+1;

                a_imjk[0]=uu[k][j][i].fx[3];
                a_imjk[1]=uu[k][j][i].fx[4];
                a_imjk_multi_a_imjk_subtr = a_imjk[0]*a_imjk[1]/(a_imjk[0]-a_imjk[1]);
                a_ipjk[0]=uu[k][j][ip].fx[3];
                a_ipjk[1]=uu[k][j][ip].fx[4];
                a_ipjk_multi_a_ipjk_subtr = a_ipjk[0]*a_ipjk[1]/(a_ipjk[0]-a_ipjk[1]);

                b_ijmk[0]=uu[k][j][i].fx[5];
                b_ijmk[1]=uu[k][j][i].fx[6];
                b_ijmk_multi_b_ijmk_subtr = b_ijmk[0]*b_ijmk[1]/(b_ijmk[0]-b_ijmk[1]);
                if (jp < Nth) {
                    b_ijpk[0]=uu[k][jp][i].fx[5];
                    b_ijpk[1]=uu[k][jp][i].fx[6];
                }
                else {
                    kc=(k+a3/2) % a3;
                    b_ijpk[0]=uu[kc][Nthm][i].fx[5];
                    b_ijpk[1]=uu[kc][Nthm][i].fx[6];
                }
                b_ijpk_multi_b_ijpk_subtr = b_ijpk[0]*b_ijpk[1]/(b_ijpk[0]-b_ijpk[1]);

                c_ijkm[0]=uu[k][j][i].fx[15];
                c_ijkm[1]=uu[k][j][i].fx[16];
                c_ijkm_multi_c_ijkm_subtr = c_ijkm[0]*c_ijkm[1]/(c_ijkm[0]-c_ijkm[1]);
                c_ijkp[0]=uu[kp][j][i].fx[15];
                c_ijkp[1]=uu[kp][j][i].fx[16];
                c_ijkp_multi_c_ijkp_subtr = c_ijkp[0]*c_ijkp[1]/(c_ijkp[0]-c_ijkp[1]);

                rhoi=0.0; rhon=0.0;
                sum_nues=0.0; sum_nueq=0.0; sum_rhonusq=0.0; sum_nues_div_ms=0.0;
                rhos_nusq_msmq=0.0; rhos_nusqmq_msmq=0.0; sum_nueq_div_mq=0.0;
                rhos_nusqms_msmq=0.0; ne=0.0; Nn=0.0;

                for (s = 0; s < sl; s++) {
                    sum_nues += nust[zk][yj][xi][s];
                    sum_nueq += nust[zk][yj][xi][s+7];

                    ne += xx[k][j][i].fx[s];
                    rhos= xx[k][j][i].fx[s]*ms[s];
                    rhoi += rhos;

                    Nn += xx[k][j][i].fx[12+s];
                    rhon += xx[k][j][i].fx[12+s]*ms[s];

                    for (t = 0; t< sm; t++) {
                        temp=rhos*nust[zk][yj][xi][21+14*s+t];
                        sum_rhonusq += temp;

                        rhos_nusq_msmq += temp/(ms[s]+ms[t]);
                        rhos_nusqmq_msmq += temp*ms[t]/(ms[s]+ms[t]);
                        rhos_nusqms_msmq += temp*ms[s]/(ms[s]+ms[t]);
                    }

                    sum_nues_div_ms =+ nust[zk][yj][xi][s]/ms[s];
                    sum_nueq_div_mq =+ nust[zk][yj][xi][7+s]/ms[s];
                }

                rhoe_sum_nueq=ne*me*sum_nueq;

                for (ir = 0; ir < a4; ir++) {
                    row.c=ir; nv=0;

                    if (ir < sl) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ir;
                        vals[nv]=a + Ls[zk][yj][xi][ir]
                                   - (rh2[ip]*a_ipjk_multi_a_ipjk_subtr + rh2[i]*a_imjk_multi_a_imjk_subtr)/rh_d3[i]
                                   - (sinth_h[jp]*b_ijpk_multi_b_ijpk_subtr + sinth_h[j]*b_ijmk_multi_b_ijmk_subtr)
                                     /rfavg_costh[j][i]
                                   - (c_ijkp_multi_c_ijkp_subtr + c_ijkm_multi_c_ijkm_subtr)/rfavg_costh_dth_dph[j][i];
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=im; col[nv].c=ir;
                        vals[nv]=rh2[i]*a_imjk_multi_a_imjk_subtr/rh_d3[i];
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=ip; col[nv].c=ir;
                        vals[nv]=rh2[ip]*a_ipjk_multi_a_ipjk_subtr/rh_d3[i];
                        nv++;

                        col[nv].k=k; col[nv].j=jm; col[nv].i=i; col[nv].c=ir;
                        vals[nv]=sinth_h[j]*b_ijmk_multi_b_ijmk_subtr/rfavg_costh[j][i];
                        nv++;

                        col[nv].k=k; col[nv].j=jp; col[nv].i=i; col[nv].c=ir;
                        vals[nv]=sinth_h[jp]*b_ijpk_multi_b_ijpk_subtr/rfavg_costh[j][i];
                        nv++;

                        col[nv].k=km; col[nv].j=j; col[nv].i=i; col[nv].c=ir;
                        vals[nv]=c_ijkm_multi_c_ijkm_subtr/rfavg_costh_dth_dph[j][i];
                        nv++;

                        col[nv].k=kp; col[nv].j=j; col[nv].i=i; col[nv].c=ir;
                        vals[nv]=c_ijkp_multi_c_ijkp_subtr/rfavg_costh_dth_dph[j][i];
                        nv++;
                    }
                    else if (ir >= 7 && ir <= 9) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=im; col[nv].c=ir;
                        vals[nv]= a_imjk_multi_a_imjk_subtr*rh2[i]/rh_d3[i];
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ir;
                        vals[nv]= a + sum_rhonusq/rhoi
                                    - (rh2[ip]*a_ipjk_multi_a_ipjk_subtr + rh2[i]*a_imjk_multi_a_imjk_subtr)/rh_d3[i]
                                    - (sinth_h[jp]*b_ijpk_multi_b_ijpk_subtr + sinth_h[j]*b_ijmk_multi_b_ijmk_subtr)
                                      /rfavg_costh[j][i]
                                    - (c_ijkp_multi_c_ijkp_subtr + c_ijkm_multi_c_ijkm_subtr)/rfavg_costh_dth_dph[j][i];
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=ip; col[nv].c=ir;
                        vals[nv]= a_ipjk_multi_a_ipjk_subtr*rh2[ip]/rh_d3[i];
                        nv++;

                        col[nv].k=k; col[nv].j=jm; col[nv].i=i; col[nv].c=ir;
                        vals[nv]= sinth_h[j]*b_ijmk_multi_b_ijmk_subtr/rfavg_costh[j][i];
                        nv++;

                        col[nv].k=k; col[nv].j=jp; col[nv].i=i; col[nv].c=ir;
                        vals[nv]=sinth_h[jp]*b_ijpk_multi_b_ijpk_subtr/rfavg_costh[j][i];
                        nv++;

                        col[nv].k=km; col[nv].j=j; col[nv].i=i; col[nv].c=ir;
                        vals[nv]=c_ijkm_multi_c_ijkm_subtr/rfavg_costh_dth_dph[j][i];
                        nv++;

                        col[nv].k=kp; col[nv].j=j; col[nv].i=i; col[nv].c=ir;
                        vals[nv]=c_ijkp_multi_c_ijkp_subtr/rfavg_costh_dth_dph[j][i];
                        nv++;

                        if (ir == 7) {
                            col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=19;
                            vals[nv]= -sum_rhonusq/rhon;
                            nv++;
                        }
                        else if (ir == 8) {
                            col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=20;
                            vals[nv]= -sum_rhonusq/rhon;
                            nv++;
                        }
                        else {
                            col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=21;
                            vals[nv]= -sum_rhonusq/rhon;
                            nv++;
                        }
                    }
                    else if (ir == 10) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=im; col[nv].c=ir;
                        vals[nv] = rh2[i]*a_imjk_multi_a_imjk_subtr/rh_d3[i];
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ir;
                        vals[nv]=a + 2.0*(me*sum_nues_div_ms+rhos_nusq_msmq)/ne
                                   - (rh2[ip]*a_ipjk_multi_a_ipjk_subtr + rh2[i]*a_imjk_multi_a_imjk_subtr)/rh_d3[i]
                                   - (sinth_h[jp]*b_ijpk_multi_b_ijpk_subtr + sinth_h[j]*b_ijmk_multi_b_ijmk_subtr)
                                     /rfavg_costh[j][i]
                                   - (c_ijkp_multi_c_ijkp_subtr + c_ijkm_multi_c_ijkm_subtr)/rfavg_costh_dth_dph[j][i];
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=ip; col[nv].c=ir;
                        vals[nv] = rh2[ip]*a_ipjk_multi_a_ipjk_subtr/rh_d3[i];
                        nv++;

                        col[nv].k=k; col[nv].j=jm; col[nv].i=i; col[nv].c=ir;
                        vals[nv] = sinth_h[j]*b_ijmk_multi_b_ijmk_subtr/rfavg_costh[j][i];
                        nv++;

                        col[nv].k=k; col[nv].j=jp; col[nv].i=i; col[nv].c=ir;
                        vals[nv] = sinth_h[jp]*b_ijpk_multi_b_ijpk_subtr/rfavg_costh[j][i];
                        nv++;

                        col[nv].k=km; col[nv].j=j; col[nv].i=i; col[nv].c=ir;
                        vals[nv] = c_ijkm_multi_c_ijkm_subtr/rfavg_costh_dth_dph[j][i];
                        nv++;

                        col[nv].k=kp; col[nv].j=j; col[nv].i=i; col[nv].c=ir;
                        vals[nv] = c_ijkp_multi_c_ijkp_subtr/rfavg_costh_dth_dph[j][i];
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=11;
                        vals[nv]=-2.0*me*sum_nues_div_ms;
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=22;
                        vals[nv]=-2.0*rhos_nusq_msmq/uu[k][j][i].fx[18];
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;
                        vals[nv]=four3rd*rhos_nusqmq_msmq*(xx[k][j][i].fx[7]/rhoi-xx[k][j][i].fx[19]/rhon);
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;
                        vals[nv]=four3rd*rhos_nusqmq_msmq*(xx[k][j][i].fx[7]/rhoi-xx[k][j][i].fx[19]/rhon);
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;
                        vals[nv]=four3rd*rhos_nusqmq_msmq*(xx[k][j][i].fx[7]/rhoi-xx[k][j][i].fx[19]/rhon);
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=19;
                        vals[nv]=-four3rd*rhos_nusqmq_msmq*(xx[k][j][i].fx[7]/rhoi-xx[k][j][i].fx[19]/rhon);
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=20;
                        vals[nv]=-four3rd*rhos_nusqmq_msmq*(xx[k][j][i].fx[7]/rhoi-xx[k][j][i].fx[19]/rhon);
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=21;
                        vals[nv]=-four3rd*rhos_nusqmq_msmq*(xx[k][j][i].fx[7]/rhoi-xx[k][j][i].fx[19]/rhon);
                        nv++;
                    }
                    else if (ir == 11) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=im; col[nv].c=ir;
                        vals[nv] = rh2[i]*a_imjk_multi_a_imjk_subtr/rh_d3[i]
                                  +one6th*(zz[k][j][ip].fx[24]-xx[k][j][ip].fx[24])/(kb*dr2*uu[k][j][im].fx[17])
                                  +two3rd*zz[k][j][i].fx[24]*(1.0/(rfavg[i]*dr)-1.0/dr2)/(kb*uu[k][j][im].fx[17]);
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ir;
                        vals[nv]=a + 2.0*me*(sum_nues_div_ms+sum_nueq_div_mq)
                                   +four3rd*zz[k][j][i].fx[24]/(kb*uu[k][j][i].fx[17])
                                           *( 1.0/dr2+1.0/(rfavg_dth[i]*rfavg_dth[i])
                                             +1.0/(rfavg_sinth_dph[j][i]*rfavg_sinth_dph[j][i]))
                                   -(rh2[ip]*a_ipjk_multi_a_ipjk_subtr + rh2[i]*a_imjk_multi_a_imjk_subtr)/rh_d3[i]
                                   -(sinth_h[jp]*b_ijpk_multi_b_ijpk_subtr + sinth_h[j]*b_ijmk_multi_b_ijmk_subtr)
                                    /rfavg_costh[j][i]
                                   -(c_ijkp_multi_c_ijkp_subtr + c_ijkm_multi_c_ijkm_subtr)/rfavg_costh_dth_dph[j][i];
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=ip; col[nv].c=ir;
                        vals[nv] = rh2[ip]*a_ipjk_multi_a_ipjk_subtr/rh_d3[i]
                                  -one6th*(zz[k][j][ip].fx[24]-xx[k][j][ip].fx[24])/(kb*dr2*uu[k][j][ip].fx[17])
                                  -two3rd*zz[k][j][i].fx[24]*(1.0/(rfavg[i]*dr)+1.0/dr2)/(kb*uu[k][j][ip].fx[17]);
                        nv++;

                        col[nv].k=k; col[nv].j=jm; col[nv].i=i; col[nv].c=ir;
                        vals[nv] = sinth_h[j]*b_ijmk_multi_b_ijmk_subtr/rfavg_costh[j][i]
                                  +one6th*(zz[kcp][jcp][i].fx[24]-zz[kcm][jcm][i].fx[24])
                                          /(kb*rfavg_dth[i]*rfavg_dth[i]*uu[kcm][jcm][i].fx[17])
                                  +two3rd*zz[k][j][i].fx[24]/(kb*rfavg_dth[i]*uu[kcm][jcm][i].fx[17])
                                                            *(cotth[j]/rfavg[i]-1.0/rfavg_dth[i]);
                        nv++;

                        col[nv].k=k; col[nv].j=jp; col[nv].i=i; col[nv].c=ir;
                        vals[nv] = sinth_h[jp]*b_ijpk_multi_b_ijpk_subtr/rfavg_costh[j][i]
                                  -one6th*(zz[kcp][jcp][i].fx[24]-zz[kcm][jcm][i].fx[24])
                                          /(kb*rfavg_dth[i]*rfavg_dth[i]*uu[kcp][jcp][i].fx[17])
                                  -two3rd*zz[k][j][i].fx[24]/(kb*rfavg_dth[i]*uu[kcp][jcp][i].fx[17])
                                                            *(cotth[j]/rfavg[i]-1.0/rfavg_dth[i]);
                        nv++;

                        col[nv].k=km; col[nv].j=j; col[nv].i=i; col[nv].c=ir;
                        vals[nv] = c_ijkm_multi_c_ijkm_subtr/rfavg_costh_dth_dph[j][i]
                                  +(one6th*(zz[kp][j][i].fx[24]-zz[km][j][i].fx[24])-two3rd*zz[k][j][i].fx[17])
                                   /(kb*rfavg_sinth_dph[j][i]*rfavg_sinth_dph[j][i]*uu[km][j][i].fx[17]);
                        nv++;

                        col[nv].k=kp; col[nv].j=j; col[nv].i=i; col[nv].c=ir;
                        vals[nv] = c_ijkp_multi_c_ijkp_subtr/rfavg_costh_dth_dph[j][i]
                                  -(one6th*(zz[kp][j][i].fx[24]-zz[km][j][i].fx[24])+two3rd*zz[k][j][i].fx[17])
                                   /(kb*rfavg_sinth_dph[j][i]*rfavg_sinth_dph[j][i]*uu[kp][j][i].fx[17]);
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10;
                        vals[nv]=-2.0*me*sum_nues_div_ms;
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=22;
                        vals[nv]=-2.0*me*sum_nueq_div_mq*uu[k][j][i].fx[17]/uu[k][j][i].fx[18];
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;
                        vals[nv]=four3rd*rhoe_sum_nueq*(xx[k][j][i].fx[7]/rhoi-xx[k][j][i].fx[19]/rhon);
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;
                        vals[nv]=four3rd*rhoe_sum_nueq*(xx[k][j][i].fx[7]/rhoi-xx[k][j][i].fx[19]/rhon);
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;
                        vals[nv]=four3rd*rhoe_sum_nueq*(xx[k][j][i].fx[7]/rhoi-xx[k][j][i].fx[19]/rhon);
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=19;
                        vals[nv]=-four3rd*rhoe_sum_nueq*(xx[k][j][i].fx[7]/rhoi-xx[k][j][i].fx[19]/rhon);
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=20;
                        vals[nv]=-four3rd*rhoe_sum_nueq*(xx[k][j][i].fx[7]/rhoi-xx[k][j][i].fx[19]/rhon);
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=21;
                        vals[nv]=-four3rd*rhoe_sum_nueq*(xx[k][j][i].fx[7]/rhoi-xx[k][j][i].fx[19]/rhon);
                        nv++;
                    }
                    else if (ir >= 12 && ir <= 18) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ir;
                        vals[nv]=a+Ls[zk][yj][xi][s-5];
                        nv++;
                    }
                    else if (ir == 19) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ir;
                        vals[nv]=a + (rhoe_sum_nueq + sum_rhonusq)/rhon;
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;
                        vals[nv]= - (rhoe_sum_nueq + sum_rhonusq)/rhoi;
                        nv++;
                    }
                    else if (ir == 20) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ir;
                        vals[nv]=a + (rhoe_sum_nueq + sum_rhonusq)/rhon;
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;
                        vals[nv]= - (rhoe_sum_nueq + sum_rhonusq)/rhoi;
                        nv++;
                    }
                    else if (ir == 21) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ir;
                        vals[nv]=a + (rhoe_sum_nueq + sum_rhonusq)/rhon;
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;
                        vals[nv]= - (rhoe_sum_nueq + sum_rhonusq)/rhoi;
                        nv++;
                    }
                    else if (ir == 22) {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ir;
                        vals[nv]=a + 2.0*rhos_nusq_msmq/uu[k][j][i].fx[18];
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=10;
                        vals[nv]= - 2.0*rhos_nusq_msmq/uu[k][j][i].fx[17];
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=7;
                        vals[nv]=four3rd*rhos_nusqms_msmq*(xx[k][j][i].fx[7]/rhoi-xx[k][j][i].fx[19]/rhon);
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=8;
                        vals[nv]=four3rd*rhos_nusqms_msmq*(xx[k][j][i].fx[7]/rhoi-xx[k][j][i].fx[19]/rhon);
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=9;
                        vals[nv]=four3rd*rhos_nusqms_msmq*(xx[k][j][i].fx[7]/rhoi-xx[k][j][i].fx[19]/rhon);
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=19;
                        vals[nv]=-four3rd*rhos_nusqms_msmq*(xx[k][j][i].fx[7]/rhoi-xx[k][j][i].fx[19]/rhon);
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=20;
                        vals[nv]=-four3rd*rhos_nusqms_msmq*(xx[k][j][i].fx[7]/rhoi-xx[k][j][i].fx[19]/rhon);
                        nv++;

                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=21;
                        vals[nv]=-four3rd*rhos_nusqms_msmq*(xx[k][j][i].fx[7]/rhoi-xx[k][j][i].fx[19]/rhon);
                        nv++;
                    }
                    else {
                        col[nv].k=k; col[nv].j=j; col[nv].i=i; col[nv].c=ir;
                        vals[nv]=1.0;
                        nv++;
                    }

                    MatSetValuesStencil(Jac,1,&row,nv,col,vals,INSERT_VALUES);
                }
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
    DMDAVecRestoreArray(da, localZ, &zz);

    DMRestoreLocalVector(da, &localX);
    DMRestoreLocalVector(da, &localU);
    DMRestoreLocalVector(da, &localZ);

    delete[] vals;
    delete[] col;

    return ierr;
}
