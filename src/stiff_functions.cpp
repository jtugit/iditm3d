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
    int    s, t, s3, s14;
    double Br, Bt, Bp, Ec_VxB[3];
    double usr[3], ust[3], usp[3], ui[3], unr, unt, unp, ns[7], nn[7], rhos[7];
    double ne, Te, Nn, rhon, Tn, TiO, TiH, TiHe; //ne = ni
    double sum_nusq[7], mtnust_msmt, nusq_msmq[7], sum_nueq_div_mq, nusttemp, nuis;
    double mqnusq_msmq[3], nust_msmt, rhossum_nusq_rhon[3], sum_nues, sum_nueq;
    double uiminusun_sq[3], uOi_uHi_sq, uHi_uHei_sq, uOi_uHei_sq;
    double neme, nuis_ms, deni; //Ce, 
    const double two3rd=2.0/3.0;

    Vec    localX, localU;
    AppCtx *params = (AppCtx*)ctx;
    Field  ***xx, ***dxdt, ***uu, ***ff;
    DM     da;
    int    xs, xm, ys, ym, zs, zm, zk, yj, xi, i, j, k, im, ip, jm, jp, km, kp;
    uint64_t kji, kj, ji;

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

    for (k = zs; k < zs+zm; k++) {
        zk=k-zs;

        km = k-1; kp = k+1;

        for (j = ys; j < ys+ym; j++) {
            if (j == 0 || j == Nth) {
                for (i = xs; i < xs+xm; i++) {
                    for (s = 0; s < a4; s++) ff[k][j][i].fx[s] = xx[k][j][i].fx[s];
                }
                continue;
            }

            yj=j-ys; jm = j-1; jp = j+1; kj=(uint64_t)(zk*ym+yj);

            for (i = xs; i < xs+xm; i++) {
                if (i == 0 || i == Nr) {
                    for (s = 0; s< a4; s++) ff[k][j][i].fx[s] = xx[k][j][i].fx[s];
                    continue;
                }

                xi=i-xs; im = i-1; ip = i+1; ji=(uint64_t)(yj*xm+xi); kji=(uint64_t)(zk*ym*xm+yj*xm+xi);

                ne=0.0; Nn=0.0; rhon=0.0; ui[0]=0.0; ui[1]=0.0; ui[2]=0.0; sum_nues=0.0; sum_nueq=0.0;
                for (s = 0; s < 7; s++) {
                    
                }
                sum_nueq_div_mq=0.0; 

                for (s = 0; s < sl; s++) {
                    ff[k][j][i].fx[s] = dxdt[k][j][i].fx[s];

                    ff[k][j][i].fx[20+s]=dxdt[k][j][i].fx[20+s];

                    sum_nues += nust[zk][yj][xi][s];
                    sum_nueq += nust[zk][yj][xi][s+7];
                    sum_nueq_div_mq += nust[zk][yj][xi][7+s]/ams[s];

                    ns[s] = exp(xx[k][j][i].fx[s]); ne += ns[s]; rhos[s]= ns[s]*ams[s];
                    nn[s] = exp(xx[k][j][i].fx[20+s]); Nn += nn[s]; rhon += nn[s]*ams[s];

                    s14=14*s; sum_nusq[s]=0.0; nusq_msmq[s]=0.0;
                    for (t = 0; t< sm; t++) {
                        nusttemp=nust[zk][yj][xi][21+s14+t]; sum_nusq[s] += nusttemp;
                        if (s < 3) mqnusq_msmq[s] += ams[t]*nusttemp/(ams[s]+ams[t]);
                        nusq_msmq[s] += nusttemp/(ams[s]+ams[t]);
                    }

                    if (s < 3) {
                        s3=s*3;
                        usr[s]=xx[k][j][i].fx[7+s3]; ust[s]=xx[k][j][i].fx[8+s3]; usp[s]=xx[k][j][i].fx[9+s3];

                        ui[0] += ns[s]*usr[s]; ui[1] += ns[s]*ust[s]; ui[2] += ns[s]*usp[s];
                    }
                }

                //ion bulk velocity
                deni=ns[3]+ns[4]+ns[5]+ns[6];
                ui[0]=(ui[0]+deni*xx[k][j][i].fx[7])/ne;
                ui[1]=(ui[1]+deni*xx[k][j][i].fx[8])/ne;
                ui[2]=(ui[2]+deni*xx[k][j][i].fx[9])/ne;

                //background magnetic field in normal spherical components
                Br=uu[k][j][i].fx[0]; Bt=uu[k][j][i].fx[1]; Bp=uu[k][j][i].fx[2];
                unr=xx[k][j][i].fx[27]; unt=xx[k][j][i].fx[28]; unp=xx[k][j][i].fx[29];

                //ion momentum equation of O+
                ff[k][j][i].fx[7] = dxdt[k][j][i].fx[7] + qms[0]*((ui[1]-ust[0])*Bp-(ui[2]-usp[0])*Bt)
                                   +sum_nusq[0]*(usr[0]-unr) + nust[zk][yj][xi][15]*(usr[0]-usr[1])
                                   +nust[zk][yj][xi][16]*(usr[0]-usr[2]);
                ff[k][j][i].fx[8] = dxdt[k][j][i].fx[8] + qms[0]*((ui[2]-usp[0])*Br-(ui[0]-usr[0])*Bp)
                                   +sum_nusq[0]*(ust[0]-unt) + nust[zk][yj][xi][15]*(ust[0]-ust[1])
                                   +nust[zk][yj][xi][16]*(ust[0]-ust[2]);
                ff[k][j][i].fx[9] = dxdt[k][j][i].fx[9] + qms[0]*((ui[0]-usr[0])*Bt-(ui[1]-ust[0])*Br)
                                   +sum_nusq[0]*(usp[0]-unp) + nust[zk][yj][xi][15]*(usp[0]-usp[1])
                                   +nust[zk][yj][xi][16]*(usp[0]-usp[2]);

                //ion momentum equation of H+
                nuis= nust[zk][yj][xi][29]+nust[zk][yj][xi][31]+nust[zk][yj][xi][32]+nust[zk][yj][xi][33]
                     +nust[zk][yj][xi][34];

                ff[k][j][i].fx[10]= dxdt[k][j][i].fx[10]+qms[1]*((ui[1]-ust[1])*Bp-(ui[2]-usp[1])*Bt)
                                   +sum_nusq[1]*(usr[1]-unr)+nuis*(usr[1]-usr[0])+nust[zk][yj][xi][30]*(usr[1]-usr[2]);
                ff[k][j][i].fx[11]= dxdt[k][j][i].fx[11]+qms[1]*((ui[2]-usp[1])*Br-(ui[0]-usr[1])*Bp)
                                   +sum_nusq[1]*(ust[1]-unt)+nuis*(ust[1]-ust[0])+nust[zk][yj][xi][30]*(ust[1]-ust[2]);
                ff[k][j][i].fx[12]= dxdt[k][j][i].fx[12]+qms[1]*((ui[0]-usr[1])*Bt-(ui[1]-ust[1])*Br)
                                   +sum_nusq[1]*(usp[1]-unp)+nuis*(usp[1]-usp[0])+nust[zk][yj][xi][30]*(usp[1]-usp[2]);

                //ion momentum equation of He+
                nuis= nust[zk][yj][xi][43]+nust[zk][yj][xi][45]+nust[zk][yj][xi][46]+nust[zk][yj][xi][47]
                     +nust[zk][yj][xi][48];

                ff[k][j][i].fx[13]= dxdt[k][j][i].fx[13]+qms[2]*((ui[1]-ust[2])*Bp-(ui[2]-usp[2])*Bt)
                                   +sum_nusq[2]*(usr[2]-unr)+nuis*(usr[2]-usr[0])+nust[zk][yj][xi][44]*(usr[2]-usr[1]);
                ff[k][j][i].fx[14]= dxdt[k][j][i].fx[14]+qms[2]*((ui[2]-usp[2])*Br-(ui[0]-usr[2])*Bp)
                                   +sum_nusq[2]*(ust[2]-unt)+nuis*(ust[2]-ust[0])+nust[zk][yj][xi][44]*(ust[2]-ust[1]);
                ff[k][j][i].fx[15]= dxdt[k][j][i].fx[15]+qms[2]*((ui[0]-usr[2])*Bt-(ui[1]-ust[2])*Br)
                                   +sum_nusq[2]*(usp[2]-unp)+nuis*(usp[2]-usp[0])+nust[zk][yj][xi][44]*(usp[2]-usp[1]);

                //----- O+ temperature equation
                Tn=xx[k][j][i].fx[30]; Te=xx[k][j][i].fx[19];
                TiO=xx[k][j][i].fx[16]; TiH=xx[k][j][i].fx[17]; TiHe=xx[k][j][i].fx[18];
                neme=ne*ame;

                uiminusun_sq[0]=(usr[0]-unr)*(usr[0]-unr)+(ust[0]-unt)*(ust[0]-unt)+(usp[0]-unp)*(usp[0]-unp);
                uOi_uHi_sq= (usr[0]-usr[1])*(usr[0]-usr[1])+(ust[0]-ust[1])*(ust[0]-ust[1])
                           +(usp[0]-usp[1])*(usp[0]-usp[1]);
                uOi_uHei_sq= (usr[0]-usr[2])*(usr[0]-usr[2])+(ust[0]-ust[2])*(ust[0]-ust[2])
                            +(usp[0]-usp[2])*(usp[0]-usp[2]);

                ff[k][j][i].fx[16] = dxdt[k][j][i].fx[16]
                                    -two3rd*ams[0]*( mqnusq_msmq[0]*uiminusun_sq[0]
                                                    +ams[1]*nust[zk][yj][xi][15]/(ams[0]+ams[1])*uOi_uHi_sq
                                                    +ams[2]*nust[zk][yj][xi][16]/(ams[0]+ams[2])*uOi_uHei_sq)
                                    +2.0*ams[0]*( nusq_msmq[0]*(TiO-Tn)+nust[zk][yj][xi][15]/(ams[0]+ams[1])*(TiO-TiH)
                                                 +nust[zk][yj][xi][16]/(ams[0]+ams[2])*(TiO-TiHe))
                                    +2.0*neme/(ns[0]*ams[0])*nust[zk][yj][xi][0]*(TiO-Te);

                //----- H+ temperature equation
                uiminusun_sq[1]=(usr[1]-unr)*(usr[1]-unr)+(ust[1]-unt)*(ust[1]-unt)+(usp[1]-unp)*(usp[1]-unp);
                mtnust_msmt= ams[0]*nust[zk][yj][xi][29]/(ams[1]+ams[0])+ams[3]*nust[zk][yj][xi][31]/(ams[1]+ams[3])
                            +ams[4]*nust[zk][yj][xi][32]/(ams[1]+ams[4])+ams[5]*nust[zk][yj][xi][33]/(ams[1]+ams[5])
                            +ams[6]*nust[zk][yj][xi][34]/(ams[1]+ams[6]);
                nust_msmt= nust[zk][yj][xi][29]/(ams[1]+ams[0])+nust[zk][yj][xi][31]/(ams[1]+ams[3])
                          +nust[zk][yj][xi][32]/(ams[1]+ams[4])+nust[zk][yj][xi][33]/(ams[1]+ams[5])
                          +nust[zk][yj][xi][34]/(ams[1]+ams[6]);

                uHi_uHei_sq= (usr[1]-usr[2])*(usr[1]-usr[2])+(ust[1]-ust[2])*(ust[1]-ust[2])
                            +(usp[1]-usp[2])*(usp[1]-usp[2]);
                ff[k][j][i].fx[17] = dxdt[k][j][i].fx[17]
                                    -two3rd*ams[1]*( mqnusq_msmq[1]*uiminusun_sq[1]+mtnust_msmt*uOi_uHi_sq
                                                    +ams[2]*nust[zk][yj][xi][30]/(ams[1]+ams[2])*uHi_uHei_sq)
                                    +2.0*ams[1]*( nusq_msmq[1]*(TiH-Tn)+nust_msmt*(TiH-TiO)
                                                 +nust[zk][yj][xi][30]/(ams[1]+ams[2])*(TiH-TiHe))
                                    +2.0*neme/(ns[1]*ams[1])*nust[zk][yj][xi][1]*(TiH-Te);

                //----- He+ temperature equation
                uiminusun_sq[2]=(usr[2]-unr)*(usr[2]-unr)+(ust[2]-unt)*(ust[2]-unt)+(usp[2]-unp)*(usp[2]-unp);
                mtnust_msmt= ams[0]*nust[zk][yj][xi][43]/(ams[2]+ams[0])+ams[3]*nust[zk][yj][xi][45]/(ams[2]+ams[3])
                            +ams[4]*nust[zk][yj][xi][46]/(ams[2]+ams[4])+ams[5]*nust[zk][yj][xi][47]/(ams[2]+ams[5])
                            +ams[6]*nust[zk][yj][xi][48]/(ams[2]+ams[6]);
                nust_msmt= nust[zk][yj][xi][43]/(ams[2]+ams[0])+nust[zk][yj][xi][45]/(ams[2]+ams[3])
                          +nust[zk][yj][xi][46]/(ams[2]+ams[4])+nust[zk][yj][xi][47]/(ams[2]+ams[5])
                          +nust[zk][yj][xi][48]/(ams[2]+ams[6]);

                ff[k][j][i].fx[18] = dxdt[k][j][i].fx[18]
                                    -two3rd*ams[2]*( mqnusq_msmq[2]*uiminusun_sq[2]+mtnust_msmt*uOi_uHei_sq
                                                    +ams[1]*nust[zk][yj][xi][44]/(ams[2]+ams[1])*uHi_uHei_sq)
                                    +2.0*ams[2]*( nusq_msmq[2]*(TiHe-Tn)+nust_msmt*(TiHe-TiO)
                                                 +nust[zk][yj][xi][44]/(ams[2]+ams[1])*(TiHe-TiH))
                                    +2.0*neme/(ns[2]*ams[2])*nust[zk][yj][xi][2]*(TiHe-Te);

                //----- electron temperature equation
                nuis_ms= nust[zk][yj][xi][0]/ams[0]+nust[zk][yj][xi][3]/ams[3]+nust[zk][yj][xi][4]/ams[4]
                        +nust[zk][yj][xi][5]/ams[5]+nust[zk][yj][xi][6]/ams[6];
                //Ce=ele_cooling_rate(xx, Te, Tn, ne, i, j, k);

                ff[k][j][i].fx[19] = dxdt[k][j][i].fx[19]
                                    +2.0*ame*( nuis_ms*(Te-TiO)+nust[zk][yj][xi][1]/ams[1]*(Te-TiH)
                                              +nust[zk][yj][xi][2]/ams[2]*(Te-TiHe)+sum_nueq_div_mq*(Te-Tn));
                                    //-two3rd*sum_nueq*ame
                                    // *((ui[0]-unr)*(ui[0]-unr)+(ui[1]-unr)*(ui[1]-unt)+(ui[2]-unp)*(ui[2]-unp))
                                    //+two3rd/ne*Ce;

                //----- neutral momentum equation
                rhossum_nusq_rhon[0]=( rhos[0]*sum_nusq[0]+rhos[3]*sum_nusq[3]+rhos[4]*sum_nusq[4]
                                      +rhos[5]*sum_nusq[5]+rhos[6]*sum_nusq[6])/rhon;
                rhossum_nusq_rhon[1]=rhos[1]*sum_nusq[1]/rhon;
                rhossum_nusq_rhon[2]=rhos[2]*sum_nusq[2]/rhon;

                ff[k][j][i].fx[27] = dxdt[k][j][i].fx[27]+rhossum_nusq_rhon[0]*(unr-usr[0])
                                                         +rhossum_nusq_rhon[1]*(unr-usr[1])
                                                         +rhossum_nusq_rhon[2]*(unr-usr[2]);
                ff[k][j][i].fx[28] = dxdt[k][j][i].fx[28]+rhossum_nusq_rhon[0]*(unt-ust[0])
                                                         +rhossum_nusq_rhon[1]*(unt-ust[1])
                                                         +rhossum_nusq_rhon[2]*(unt-ust[2]);
                ff[k][j][i].fx[29] = dxdt[k][j][i].fx[29]+rhossum_nusq_rhon[0]*(unp-usp[0])
                                                         +rhossum_nusq_rhon[1]*(unp-usp[1])
                                                         +rhossum_nusq_rhon[2]*(unp-usp[2]);

                //neutral temperature equation
                nuis_ms = rhos[0]*ams[0]*nusq_msmq[0]+rhos[3]*ams[3]*nusq_msmq[3]+rhos[4]*ams[4]*nusq_msmq[4]
                         +rhos[5]*ams[5]*nusq_msmq[5]+rhos[6]*ams[6]*nusq_msmq[6];
                nuis = rhos[0]*nusq_msmq[0]+rhos[3]*nusq_msmq[3]+rhos[4]*nusq_msmq[4]+rhos[5]*nusq_msmq[5]
                      +rhos[6]*nusq_msmq[6];
                ff[k][j][i].fx[30] =( 2.0*(nuis*(Tn-TiO)+rhos[1]*nusq_msmq[1]*(Tn-TiH)+rhos[2]*nusq_msmq[2]*(Tn-TiHe))
                                     -two3rd*( nuis_ms*uiminusun_sq[0]+rhos[1]*ams[1]*nusq_msmq[1]*uiminusun_sq[1]
                                              +rhos[2]*ams[2]*nusq_msmq[2]*uiminusun_sq[2]))/Nn;

                //----------- B^r equation
                ff[k][j][i].fx[31]= dxdt[k][j][i].fx[31]+0.5*( (xx[k][jp][i].fx[36]-xx[k][jm][i].fx[36])/dth
                                                              -(xx[kp][j][i].fx[35]-xx[km][j][i].fx[35])/dph);

                //----------- B^theta equation
                ff[k][j][i].fx[32]= dxdt[k][j][i].fx[32]+0.5*( (xx[kp][j][i].fx[34]-xx[km][j][i].fx[34])/dph
                                                              -(xx[k][j][ip].fx[36]-xx[k][j][im].fx[36])/dr);

                //----------- B^phi equation
                ff[k][j][i].fx[33]= dxdt[k][j][i].fx[33]+0.5*( (xx[k][j][ip].fx[35]-xx[k][j][im].fx[35])/dr
                                                              -(xx[k][jp][i].fx[34]-xx[k][jm][i].fx[34])/dth);

                electric_field_vxB(xx, uu, ui, i, j, k, xi, yj, zk, xm, ym, Ec_VxB);

                //----------- E^r equation
                ff[k][j][i].fx[34]= xx[k][j][i].fx[23]
                                   -( Jinv.Jiv11[kj]*Ec_VxB[0]+Jinv.Jiv21[kj]*Ec_VxB[1]
                                     +Jinv.Jiv31[(uint64_t)yj]*Ec_VxB[2]);

                //----------- E^theta equation
                ff[k][j][i].fx[35]= xx[k][j][i].fx[35]
                                   -(Jinv.Jiv12[kji]*Ec_VxB[0]+Jinv.Jiv22[kji]*Ec_VxB[1]+Jinv.Jiv32[ji]*Ec_VxB[2]);

                //----------- E^phi equation
                ff[k][j][i].fx[36]= xx[k][j][i].fx[36]-(Jinv.Jiv13[kji]*Ec_VxB[0]+Jinv.Jiv23[kji]*Ec_VxB[1]);

                for (s=0; s<a4; s++) {
                    if (isnan(ff[k][j][i].fx[s]) || isinf(ff[k][j][i].fx[s])) {
                        cout<<"Stiff function is Nan or inf at ("<<i<<", "<<j<<", "<<k<<", "<<s
                            <<") in stiffunction "<<ff[k][j][i].fx[s]<<" "<<a4<<endl;
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