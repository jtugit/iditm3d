#include "funcdef.h"
#include "boundary.h"

#include <fstream>
#include <iostream>

using namespace std;

#include "electric_field.h"
#include "ele_cooling_rate.h"
#include "operators.h"

int stifffunction(TS ts, double ftime, Vec X, Vec Xdt, Vec F, void* ctx)
{
    int    s, t, jcm, jcp, s3; //kcm, kcp,
    double Bx, By, Bz, Br, Bt, Bp, Er, Rt, Ep, Epartx, Eparty, Epartz, uex, uey, uez;
    double uer, uet, uep, usr[3], ust[3], usp[3], uir, uit, uip, unr, unt, unp, ns[3], rhos[7], Ts[3];
    double ne, Te, rhoe, Nn, rhon, Tn; //ne = ni
    double nusn[7], mtnust_msmt, sum_mqnusq_msmq, nusq_msmq[7], nuen, sum_nueq_div_mq, nusttemp, nuis;
    double mqnusq_msmq[3], nust_msmt, rhosnusn_rhon[3];
    double sum_nues, sum_nueq, sum_rhonusq, sum_nues_div_ms;
    double rhos_nusq_msmq[7], rhos_nusqmq_msmq, rhos_nusqms_msmq[7], temp;
    double rhoe_sum_nueq, uiminusun_sq[3], uOi_uHi_sq, uHi_uHei_sq, uOi_uHei_sq, ueminusui_sq, ueminusun_sq;
    double Qefric, Qifric, Qnfric, me_kb, neme;
    double nuis_ms, nuin[3];
    const double two3rd=2.0/3.0, one6th=1.0/6.0, kb_div_e=kb/e;

    Vec    localX, localU;
    AppCtx *params = (AppCtx*)ctx;
    Field  ***xx, ***dxdt, ***uu, ***ff;
    DM     da;
    int    xs, xm, ys, ym, zs, zm, zk, yj, xi, i, j, k, im, ip, jm, jp, km, kp, kji, kj, ji;

    params->ftime = ftime;
    TSGetDM(ts, &da);

    DMGetLocalVector(da, &localX);
    DMGlobalToLocalBegin(da, X, INSERT_VALUES,localX);
    DMGlobalToLocalEnd(da, X, INSERT_VALUES, localX);
    DMDAVecGetArray(da, localX, &xx);

    DMGetLocalVector(da, &localU);
    DMGlobalToLocalBegin(da, params->U, INSERT_VALUES,localU);
    DMGlobalToLocalEnd(da, params->U, INSERT_VALUES, localU);
    DMDAVecGetArray(da, localU, &uu);

    DMDAVecGetArray(da, Xdt, &dxdt);
    DMDAVecGetArray(da, F, &ff);

    DMDAGetCorners(da, &xs ,&ys, &zs, &xm, &ym, &zm);

    double dr2=dr*dr, deni, ene, ene_rhos;

    int rank, s14;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (k = zs; k < zs+zm; k++) {
        zk=k-zs;

        km = k-1; kp = k+1;

        for (j = ys; j < ys+ym; j++) {
            if (j == 0 || j == Nthm) {
                for (i = xs; i < xs+xm; i++) {
                    for (s = 0; s < a4; s++) ff[k][j][i].fx[s] = xx[k][j][i].fx[s];
                }
                continue;
            }

            yj=j-ys; jm = j-1; jp = j+1; kj=k*ym+j;

            for (i = xs; i < xs+xm; i++) {
                if (i == 0 || i == Nr) {
                    for (s = 0; s< a4; s++) ff[k][j][i].fx[s] = xx[k][j][i].fx[s];
                    continue;
                }

                xi=i-xs; im = i-1; ip = i+1; ji=j*xm+i; kji=k*ym*xm+j*ym+i;

                ne=0.0; Nn=0.0; rhon=0.0; uer=0.0; uet=0.0; uep=0.0;

                for (s = 0; s < 7; s++) {
                    nusn[s]=0.0; nusq_msmq[s]=0.0;
                }
                nuen=0.0; sum_nueq_div_mq=0.0;

                for (s = 0; s < sl; s++) {
                    ff[k][j][i].fx[s] = dxdt[k][j][i].fx[s];

                    ff[k][j][i].fx[20+s]=dxdt[k][j][i].fx[20+s];

                    sum_nues += nust[zk][yj][xi][s];
                    sum_nueq += nust[zk][yj][xi][s+7];

                    ne += xx[k][j][i].fx[s];
                    rhos[s]= xx[k][j][i].fx[s]*ms[s];

                    Nn += xx[k][j][i].fx[20+s];
                    rhon += xx[k][j][i].fx[20+s]*ms[s];

                    s14=14*s;
                    for (t = 0; t< sm; t++) {
                        nusttemp=nust[zk][yj][xi][21+s14+t];
                        nusn[s] += nusttemp;
                        mqnusq_msmq[s] += ms[t]*nusttemp/(ms[s]+ms[t]);
                        nusq_msmq[s] += nusttemp/(ms[s]+ms[t]);
                    }

                    sum_nueq_div_mq =+ nust[zk][yj][xi][7+s]/ms[s];

                    nuen += nust[zk][yj][xi][s+7];

                    if (s < 3) {
                        s3=s*3;
                        usr[s]=xx[k][j][i].fx[7+s3];
                        ust[s]=xx[k][j][i].fx[8+s3];
                        usp[s]=xx[k][j][i].fx[9+s3];

                        uer += xx[k][j][i].fx[s]*usr[s];
                        uet += xx[k][j][i].fx[s]*ust[s];
                        uep += xx[k][j][i].fx[s]*usp[s];
                    }
                }

                //ion bulk velocity
                deni=xx[k][j][i].fx[3]+xx[k][j][i].fx[4]+xx[k][j][i].fx[5]+xx[k][j][i].fx[6];
                uir=(uer+deni*xx[k][j][i].fx[7])/ne;
                uit=(uet+deni*xx[k][j][i].fx[8])/ne;
                uip=(uep+deni*xx[k][j][i].fx[9])/ne;

                //magnetic field in normal spherical components plus background magnetic field
                Br=uu[k][j][i].fx[0]; Bt=uu[k][j][i].fx[1]; Bp=uu[k][j][i].fx[2];

                ene=e*ne;

                //ion momentum equation of O+
                ff[k][j][i].fx[7] = dxdt[k][j][i].fx[7] + qms[0]*((uit-ust[0])*Bp-(uip-usp[0])*Bt)
                                   +nusn[0]*(usr[0]-unr) + nust[zk][yj][xi][15]*(usr[0]-usr[1])
                                   +nust[zk][yj][xi][16]*(usr[0]-usr[2]);
                ff[k][j][i].fx[8] = dxdt[k][j][i].fx[8] + qms[0]*((uip-usp[0])*Br-(uir-usr[0])*Bp)
                                   +nusn[0]*(ust[0]-unt) + nust[zk][yj][xi][15]*(ust[0]-ust[1])
                                   +nust[zk][yj][xi][16]*(ust[0]-ust[2]);
                ff[k][j][i].fx[9] = dxdt[k][j][i].fx[9] + qms[0]*((uir-usr[0])*Bt-(uit-ust[0])*Br)
                                   +nusn[0]*(usp[0]-unp) + nust[zk][yj][xi][15]*(usp[0]-usp[1])
                                   +nust[zk][yj][xi][16]*(usp[0]-usp[2]);

                //ion momentum equation of H+
                nuis= nust[zk][yj][xi][29]+nust[zk][yj][xi][31]+nust[zk][yj][xi][32]+nust[zk][yj][xi][33]
                     +nust[zk][yj][xi][34];

                ff[k][j][i].fx[10]= dxdt[k][j][i].fx[10] + qms[1]*((uit-ust[1])*Bp-(uip-usp[1])*Bt)
                                   +nusn[1]*(usr[1]-unr) + nuis*(usr[1]-usr[0]) + nust[zk][yj][xi][30]*(usr[1]-usr[2]);
                ff[k][j][i].fx[11]= dxdt[k][j][i].fx[11] + qms[1]*((uip-usp[1])*Br-(uir-usr[1])*Bp)
                                   +nusn[1]*(ust[1]-unt) + nuis*(ust[1]-ust[0]) + nust[zk][yj][xi][30]*(ust[1]-ust[2]);
                ff[k][j][i].fx[12]= dxdt[k][j][i].fx[12] + qms[1]*((uir-usr[1])*Bt-(uit-ust[1])*Br)
                                   +nusn[1]*(usp[1]-unp) + nuis*(usp[1]-usp[0]) + nust[zk][yj][xi][30]*(usp[1]-usp[2]);

                //ion momentum equation of He+
                nuis= nust[zk][yj][xi][43]+nust[zk][yj][xi][45]+nust[zk][yj][xi][46]+nust[zk][yj][xi][47]
                     +nust[zk][yj][xi][48];

                ff[k][j][i].fx[13]= dxdt[k][j][i].fx[13] + qms[2]*((uit-ust[2])*Bp-(uip-usp[2])*Bt)
                                   +nusn[2]*(usr[2]-unr) + nuis*(usr[2]-usr[0]) + nust[zk][yj][xi][44]*(usr[2]-usr[1]);
                ff[k][j][i].fx[14]= dxdt[k][j][i].fx[14] + qms[2]*((uip-usp[2])*Br-(uir-usr[2])*Bp)
                                   +nusn[2]*(ust[2]-unt) + nuis*(ust[2]-ust[0]) + nust[zk][yj][xi][44]*(ust[2]-ust[1]);
                ff[k][j][i].fx[15]= dxdt[k][j][i].fx[15] + qms[2]*((uir-usr[2])*Bt-(uit-ust[2])*Br)
                                   +nusn[2]*(usp[2]-unp) + nuis*(usp[2]-usp[0]) + nust[zk][yj][xi][44]*(usp[2]-usp[1]);

                //----- O+ temperature equation
                Tn=xx[k][j][i].fx[30]; Te=xx[k][j][i].fx[19];
                Ts[0]=xx[k][j][i].fx[16]; Ts[1]=xx[k][j][i].fx[17]; Ts[2]=xx[k][j][i].fx[18];
                unr=xx[k][j][i].fx[27]; unt=xx[k][j][i].fx[28]; unp=xx[k][j][i].fx[29];
                neme=ne*me;

                uiminusun_sq[0]=(usr[0]-unr)*(usr[0]-unr)+(ust[0]-unt)*(ust[0]-unt)+(usp[0]-unp)*(usp[0]-unp);
                uOi_uHi_sq=(usr[0]-usr[1])*(usr[0]-usr[1])+(ust[0]-ust[1])*(ust[0]-ust[1])+(usp[0]-usp[1])*(usp[0]-usp[1]);
                uOi_uHei_sq=(usr[0]-usr[2])*(usr[0]-usr[2])+(ust[0]-ust[2])*(ust[0]-ust[2])+(usp[0]-usp[2])*(usp[0]-usp[2]);

                ff[k][j][i].fx[16] = dxdt[k][j][i].fx[16]
                                    -two3rd*( mqnusq_msmq[0]*uiminusun_sq[0]
                                             +ms[1]*nust[zk][yj][xi][15]/(ms[0]+ms[1])*uOi_uHi_sq
                                             +ms[2]*nust[zk][yj][xi][16]/(ms[0]+ms[2])*uOi_uHei_sq)/ams[0]
                                    +2.0*ms[0]*( nusq_msmq[0]*(Ts[0]-Tn) + nust[zk][yj][xi][15]/(ms[0]+ms[1])*(Ts[0]-Ts[1])
                                                +nust[zk][yj][xi][16]/(ms[0]+ms[2])*(Ts[0]-Ts[2])
                                                +neme/(ns[0]*ms[0])*nust[zk][yj][xi][0]*(Ts[0]-Te));

                //----- H+ temperature equation
                uiminusun_sq[1]=(usr[1]-unr)*(usr[1]-unr)+(ust[1]-unt)*(ust[1]-unt)+(usp[1]-unp)*(usp[1]-unp);
                mtnust_msmt= ms[0]*nust[zk][yj][xi][29]/(ms[1]+ms[0])+ms[3]*nust[zk][yj][xi][31]/(ms[1]+ms[3])
                            +ms[4]*nust[zk][yj][xi][32]/(ms[1]+ms[4])+ms[5]*nust[zk][yj][xi][33]/(ms[1]+ms[5])
                            +ms[6]*nust[zk][yj][xi][34]/(ms[1]+ms[6]);
                nust_msmt= nust[zk][yj][xi][29]/(ms[1]+ms[0])+nust[zk][yj][xi][31]/(ms[1]+ms[3])
                          +nust[zk][yj][xi][32]/(ms[1]+ms[4])+nust[zk][yj][xi][33]/(ms[1]+ms[5])
                          +nust[zk][yj][xi][34]/(ms[1]+ms[6]);

                uHi_uHei_sq=(usr[1]-usr[2])*(usr[1]-usr[2])+(ust[1]-ust[2])*(ust[1]-ust[2])+(usp[1]-usp[2])*(usp[1]-usp[2]);
                ff[k][j][i].fx[17] = dxdt[k][j][i].fx[17]
                                    -two3rd*( mqnusq_msmq[1]*uiminusun_sq[1]+mtnust_msmt*uOi_uHi_sq
                                             +ms[2]*nust[zk][yj][xi][30]/(ms[1]+ms[2])*uHi_uHei_sq)/ams[1]
                                    +2.0*ms[1]*( nusq_msmq[1]*(Ts[1]-Tn)+nust_msmt*(Ts[1]-Ts[0])
                                                +nust[zk][yj][xi][30]/(ms[1]+ms[2])*(Ts[1]-Ts[2])
                                                +neme/(ns[1]*ms[1])*nust[zk][yj][xi][28]*(Ts[1]-Te));

                //----- He+ temperature equation
                uiminusun_sq[2]=(usr[2]-unr)*(usr[2]-unr)+(ust[2]-unt)*(ust[2]-unt)+(usp[2]-unp)*(usp[2]-unp);
                mtnust_msmt= ms[0]*nust[zk][yj][xi][43]/(ms[2]+ms[0])+ms[3]*nust[zk][yj][xi][45]/(ms[2]+ms[3])
                            +ms[4]*nust[zk][yj][xi][46]/(ms[2]+ms[4])+ms[5]*nust[zk][yj][xi][47]/(ms[2]+ms[5])
                            +ms[6]*nust[zk][yj][xi][48]/(ms[2]+ms[6]);
                nust_msmt= nust[zk][yj][xi][43]/(ms[2]+ms[0])+nust[zk][yj][xi][45]/(ms[2]+ms[3])
                          +nust[zk][yj][xi][46]/(ms[2]+ms[4])+nust[zk][yj][xi][47]/(ms[2]+ms[5])
                          +nust[zk][yj][xi][48]/(ms[2]+ms[6]);

                ff[k][j][i].fx[18] = dxdt[k][j][i].fx[18]
                                    -two3rd*( mqnusq_msmq[2]*uiminusun_sq[2]+mtnust_msmt*uOi_uHei_sq
                                             +ms[1]*nust[zk][yj][xi][44]/(ms[2]+ms[1])*uHi_uHei_sq)/ams[1]
                                    +2.0*ms[2]*( nusq_msmq[2]*(Ts[2]-Tn)+nust_msmt*(Ts[2]-Ts[0])
                                                +nust[zk][yj][xi][44]/(ms[2]+ms[1])*(Ts[2]-Ts[1])
                                                +neme/(ns[2]*ms[2])*nust[zk][yj][xi][42]*(Ts[2]-Te));

                //----- electron temperature equation
                nuis_ms= nust[zk][yj][xi][0]/ms[0]+nust[zk][yj][xi][3]/ms[3]+nust[zk][yj][xi][4]/ms[4]
                        +nust[zk][yj][xi][5]/ms[5]+nust[zk][yj][xi][6]/ms[6];
                nuis= nust[zk][yj][xi][0]+nust[zk][yj][xi][3]+nust[zk][yj][xi][4]+nust[zk][yj][xi][5]
                     +nust[zk][yj][xi][6];

                ff[k][j][i].fx[19] = dxdt[k][j][i].fx[19]
                                    +2.0*me*( nuis_ms*(Te-Ts[0])+nust[zk][yj][xi][1]/ms[1]*(Te-Ts[1])
                                             +nust[zk][yj][xi][2]/ms[2]*(Te-Ts[2])+sum_nueq_div_mq*(Te-Tn))
                                    -two3rd*me_kb*((uer-unr)*(uer-unr)+(uer-unr)*(uet-unt)+(uep-unp)*(uep-unp));

                //----- neutral momentum equation
                rhosnusn_rhon[0]=(rhos[0]*nuin[0]+rhos[3]*nuin[3]+rhos[4]*nuin[4]+rhos[5]*nuin[5]+rhos[6]*nuin[6])/rhon;
                rhosnusn_rhon[1]=rhos[1]*nuin[1]/rhon; rhosnusn_rhon[2]=rhos[2]*nuin[2]/rhon;

                ff[k][j][i].fx[27] = dxdt[k][j][i].fx[27]+rhosnusn_rhon[0]*(unr-usr[0])+rhosnusn_rhon[1]*(unr-usr[1])
                                                         +rhosnusn_rhon[2]*(unr-usr[2]);
                ff[k][j][i].fx[28] = dxdt[k][j][i].fx[28]+rhosnusn_rhon[0]*(unt-ust[0])+rhosnusn_rhon[1]*(unt-ust[1])
                                                         +rhosnusn_rhon[2]*(unt-ust[2]);
                ff[k][j][i].fx[29] = dxdt[k][j][i].fx[29]+rhosnusn_rhon[0]*(unp-usp[0])+rhosnusn_rhon[1]*(unp-usp[1])
                                                         +rhosnusn_rhon[2]*(unp-usp[2]);

                //neutral temperature equation
                nuis_ms = rhos_nusqms_msmq[0]+rhos_nusqms_msmq[3]+rhos_nusqms_msmq[4]+rhos_nusqms_msmq[5]
                         +rhos_nusqms_msmq[6];
                nuis = rhos_nusq_msmq[0]+rhos_nusq_msmq[3]+rhos_nusq_msmq[4]+rhos_nusq_msmq[5]+rhos_nusq_msmq[6];
                ff[k][j][i].fx[30] = dxdt[k][j][i].fx[30] 
                                    -two3rd/(Nn*kb)*( nuis_ms*uiminusun_sq[0]+rhos_nusqms_msmq[1]*uiminusun_sq[1]
                                                     +rhos_nusqms_msmq[2]*uiminusun_sq[2])
                                    +2.0/Nn*(nuis*(Tn-Ts[0])+rhos_nusq_msmq[1]*(Tn-Ts[1])+rhos_nusq_msmq[2]*(Tn-Ts[2]));

                //----------- B^r equation
                ff[k][j][i].fx[31]=0.5*( (xx[kp][j][i].fx[35]-xx[km][j][i].fx[35])/dph
                                        -(xx[k][jp][i].fx[36]-xx[k][jm][i].fx[36])/dth);

                //----------- B^theta equation
                ff[k][j][i].fx[32]=0.5*( (xx[k][j][ip].fx[36]-xx[k][j][im].fx[36])/dr
                                        -(xx[kp][j][i].fx[34]-xx[km][j][i].fx[34])/dph);

                //----------- B^phi equation
                ff[k][j][i].fx[33]=0.5*( (xx[k][jp][i].fx[35]-xx[k][jm][i].fx[35])/dth
                                        -(xx[k][j][ip].fx[34]-xx[k][j][im].fx[34])/dr);

                //magnetic field components in terms of special spherical components
                Bx=( Jinv.Jiv11[kj]*xx[k][j][i].fx[31]+Jinv.Jiv12[kji]*xx[k][j][i].fx[32]
                    +Jinv.Jiv13[kji]*xx[k][j][i].fx[33])/r2sintheta[yj][xi];
                By=( Jinv.Jiv21[kj]*xx[k][j][i].fx[31]+Jinv.Jiv22[kji]*xx[k][j][i].fx[32]
                    +Jinv.Jiv23[ji] *xx[k][j][i].fx[33])/r2sintheta[yj][xi];
                Bz=(Jinv.Jiv31[yj]*xx[k][j][i].fx[31]+Jinv.Jiv32[ji] *xx[k][j][i].fx[32])/r2sintheta[yj][xi];

                //uex, uey, uez in terms of uer, and ue_theta, and ue_phi
                uex = (Kmat.K11[kj]*uu[k][j][i].fx[7]+Kmat.K12[kj]*uu[k][j][i].fx[8]+Kmat.K13[zk]*uu[k][j][i].fx[9]);
                uey = (Kmat.K21[kj]*uu[k][j][i].fx[7]+Kmat.K22[kj]*uu[k][j][i].fx[8]+Kmat.K23[zk]*uu[k][j][i].fx[9]);
                uez = (Kmat.K31[yj]*uu[k][j][i].fx[7]+Kmat.K32[yj]*uu[k][j][i].fx[8]);

                //Cartesian components of the electric field: -ue X B part
                Epartx=By*uez - Bz*uey, Eparty=Bz*uex - Bx*uez, Epartz=Bx*uey - By*uex;

                //----------- E^r equation
                ff[k][j][i].fx[34]=Jinv.Jiv11[kj] *Epartx+Jinv.Jiv21[kj] *Eparty+Jinv.Jiv31[yj]*Epartz;

                //----------- E^theta equation
                ff[k][j][i].fx[35]=Jinv.Jiv12[kji]*Epartx+Jinv.Jiv22[kji]*Eparty+Jinv.Jiv32[ji]*Epartz;

                //----------- E^phi equation
                ff[k][j][i].fx[36]=Jinv.Jiv13[kji]*Epartx+Jinv.Jiv23[kji]*Eparty;

                for (s=0; s<a4; s++) {
                    if (isnan(ff[k][j][i].fx[s]) || isinf(ff[k][j][i].fx[s])) {
                        cout<<"Stiff function is Nan or inf at ("<<i<<", "<<j<<", "<<k<<", "<<s
                            <<") in stiffunction "<<ff[k][j][i].fx[s]<<endl;
                        exit(-1);
                    }
                }
            }
        }
    }

    DMDAVecRestoreArray(da, localX, &xx);
    DMRestoreLocalVector(da, &localX);

    DMDAVecRestoreArray(da, localU, &uu);
    DMRestoreLocalVector(da, &localU);

    DMDAVecRestoreArray(da, Xdt, &dxdt);
    DMDAVecRestoreArray(da, F, &ff);

    return 0;
}