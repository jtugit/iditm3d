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
    int    s;
    double Br, Bt, Bp, Ec_VxB[3];
    double u7, u10, u13, u8, u11, u14, u9, u12, u15, ui[3], unr, unt, unp, ns[7], rhos[7];
    double ne, Te, Nn, rhon, Tn, TiO, TiH, TiHe, nuis, rhossum_nusq_rhon[3]; //, sum_nues
    double ui_un_sq[3], uOi_uHi_sq, uHi_uHei_sq, uOi_uHei_sq, neme, nuis_ms, nsmore; //Ce, 
    const double two3rd=2.0/3.0;

    Vec    localX, localU;
    AppCtx *params = (AppCtx*)ctx;
    Field  ***xx, ***dxdt, ***uu, ***ff;
    DM     da;
    int    xs, xm, ys, ym, zs, zm, zk, yj, xi, i, j, k, im, ip, jm, jp, km, kp;

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

            yj=j-ys; jm = j-1; jp = j+1;

            for (i = xs; i < xs+xm; i++) {
                if (i == 0 || i == Nr) {
                    for (s = 0; s< a4; s++) ff[k][j][i].fx[s] = xx[k][j][i].fx[s];
                    continue;
                }

                xi=i-xs; im = i-1; ip = i+1;

                ui[0]=0.0; ui[1]=0.0; ui[2]=0.0; //sum_nues=0.0;

                for (s = 0; s < sl; s++) {
                    ff[k][j][i].fx[s] = dxdt[k][j][i].fx[s];

                    ff[k][j][i].fx[20+s]=dxdt[k][j][i].fx[20+s];

                    //ion densities are from last time step
                    ns[s] = uu[k][j][i].fx[s+27]; rhos[s]= ns[s]*ams[s]; 
                }

                u7=xx[k][j][i].fx[7]; u10=xx[k][j][i].fx[10]; u13=xx[k][j][i].fx[13];
                u8=xx[k][j][i].fx[8]; u11=xx[k][j][i].fx[11]; u14=xx[k][j][i].fx[14];
                u9=xx[k][j][i].fx[9]; u12=xx[k][j][i].fx[12]; u15=xx[k][j][i].fx[15];

                //ion bulk velocity
                nsmore=ns[0]+ns[3]+ns[4]+ns[5]+ns[6]; ne=uu[k][j][i].fx[6];
                ui[0]=(nsmore*u7+ns[1]*u10+ns[2]*u13)/ne;
                ui[1]=(nsmore*u8+ns[1]*u11+ns[2]*u14)/ne;
                ui[2]=(nsmore*u9+ns[1]*u12+ns[2]*u15)/ne;

                //background magnetic field in normal spherical components
                Br=uu[k][j][i].fx[0]; Bt=uu[k][j][i].fx[1]; Bp=uu[k][j][i].fx[2];

                //neutral velocity
                unr=xx[k][j][i].fx[27]; unt=xx[k][j][i].fx[28]; unp=xx[k][j][i].fx[29];

                //ion momentum equation of O+
                ff[k][j][i].fx[7] = dxdt[k][j][i].fx[7] + qms[0]*((ui[1]-u8)*Bp-(ui[2]-u9)*Bt)
                                   +nust[zk][yj][xi][11]*(u7-unr)+nust[zk][yj][xi][7]*(u7-u10)
                                   +nust[zk][yj][xi][8]*(u7-u13);
                ff[k][j][i].fx[8] = dxdt[k][j][i].fx[8] + qms[0]*((ui[2]-u9)*Br-(ui[0]-u7)*Bp)
                                   +nust[zk][yj][xi][11]*(u8-unt)+nust[zk][yj][xi][7]*(u8-u11)
                                   +nust[zk][yj][xi][8]*(u8-u14);
                ff[k][j][i].fx[9] = dxdt[k][j][i].fx[9] + qms[0]*((ui[0]-u7)*Bt-(ui[1]-u8)*Br)
                                   +nust[zk][yj][xi][11]*(u9-unp)+nust[zk][yj][xi][7]*(u9-u12)
                                   +nust[zk][yj][xi][8]*(u9-u15);

                //ion momentum equation of H+
                ff[k][j][i].fx[10]= dxdt[k][j][i].fx[10]+qms[1]*((ui[1]-u11)*Bp-(ui[2]-u12)*Bt)
                                   +nust[zk][yj][xi][19]*(u10-unr)+nust[zk][yj][xi][14]*(u10-u7)
                                   +nust[zk][yj][xi][15]*(u10-u13);
                ff[k][j][i].fx[11]= dxdt[k][j][i].fx[11]+qms[1]*((ui[2]-u12)*Br-(ui[0]-u10)*Bp)
                                   +nust[zk][yj][xi][19]*(u11-unt)+nust[zk][yj][xi][14]*(u11-u8)
                                   +nust[zk][yj][xi][15]*(u11-u14);
                ff[k][j][i].fx[12]= dxdt[k][j][i].fx[12]+qms[1]*((ui[0]-u10)*Bt-(ui[1]-u11)*Br)
                                   +nust[zk][yj][xi][19]*(u12-unp)+nust[zk][yj][xi][14]*(u12-u9)
                                   +nust[zk][yj][xi][15]*(u12-u15);

                //ion momentum equation of He+
                ff[k][j][i].fx[13]= dxdt[k][j][i].fx[13]+qms[2]*((ui[1]-u14)*Bp-(ui[2]-u15)*Bt)
                                   +nust[zk][yj][xi][27]*(u13-unr)+nust[zk][yj][xi][22]*(u13-u7)
                                   +nust[zk][yj][xi][23]*(u13-u10);
                ff[k][j][i].fx[14]= dxdt[k][j][i].fx[14]+qms[2]*((ui[2]-u15)*Br-(ui[0]-u13)*Bp)
                                   +nust[zk][yj][xi][27]*(u14-unt)+nust[zk][yj][xi][22]*(u14-u8)
                                   +nust[zk][yj][xi][23]*(u14-u11);
                ff[k][j][i].fx[15]= dxdt[k][j][i].fx[15]+qms[2]*((ui[0]-u13)*Bt-(ui[1]-u14)*Br)
                                   +nust[zk][yj][xi][27]*(u15-unp)+nust[zk][yj][xi][22]*(u15-u9)
                                   +nust[zk][yj][xi][23]*(u15-u12);

                //----- O+ temperature equation
                Tn=xx[k][j][i].fx[30]; Te=xx[k][j][i].fx[19];
                TiO=xx[k][j][i].fx[16]; TiH=xx[k][j][i].fx[17]; TiHe=xx[k][j][i].fx[18];
                neme=ne*ame;

                ui_un_sq[0]=(u7-unr)*(u7-unr)+(u8-unt)*(u8-unt)+(u9-unp)*(u9-unp);
                uOi_uHi_sq =(u7-u10)*(u7-u10)+(u8-u11)*(u8-u11)+(u9-u12)*(u9-u12);
                uOi_uHei_sq=(u7-u13)*(u7-u13)+(u8-u14)*(u8-u14)+(u9-u15)*(u9-u15);

                ff[k][j][i].fx[16] = dxdt[k][j][i].fx[16]
                                    -two3rd*ams[0]*( nust[zk][yj][xi][13]*ui_un_sq[0]
                                                    +ams[1]*nust[zk][yj][xi][9] *uOi_uHi_sq
                                                    +ams[2]*nust[zk][yj][xi][10]*uOi_uHei_sq)
                                    +2.0*ams[0]*( nust[zk][yj][xi][12]*(TiO-Tn)
                                                 +nust[zk][yj][xi][9] *(TiO-TiH)
                                                 +nust[zk][yj][xi][10]*(TiO-TiHe))
                                    +2.0*neme/(ns[0]*ams[0])*nust[zk][yj][xi][0]*(TiO-Te);

                //----- H+ temperature equation
                ui_un_sq[1]=(u10-unr)*(u10-unr)+(u11-unt)*(u11-unt)+(u12-unp)*(u12-unp);

                uHi_uHei_sq= (u10-u13)*(u10-u13)+(u11-u14)*(u11-u14)+(u12-u15)*(u12-u15);
                ff[k][j][i].fx[17] = dxdt[k][j][i].fx[17]
                                    -two3rd*ams[1]*( nust[zk][yj][xi][21]*ui_un_sq[1]
                                                    +nust[zk][yj][xi][18]*uOi_uHi_sq
                                                    +ams[2]*nust[zk][yj][xi][17]*uHi_uHei_sq)
                                    +2.0*ams[1]*( nust[zk][yj][xi][20]*(TiH-Tn)
                                                 +nust[zk][yj][xi][16]*(TiH-TiO)
                                                 +nust[zk][yj][xi][17]*(TiH-TiHe))
                                    +2.0*neme/(ns[1]*ams[1])*nust[zk][yj][xi][1]*(TiH-Te);

                //----- He+ temperature equation
                ui_un_sq[2]=(u13-unr)*(u13-unr)+(u14-unt)*(u14-unt)+(u15-unp)*(u15-unp);

                ff[k][j][i].fx[18] = dxdt[k][j][i].fx[18]
                                    -two3rd*ams[2]*( nust[zk][yj][xi][29]*ui_un_sq[2]
                                                    +nust[zk][yj][xi][26]*uOi_uHei_sq
                                                    +ams[1]*nust[zk][yj][xi][25]*uHi_uHei_sq)
                                    +2.0*ams[2]*( nust[zk][yj][xi][28]*(TiHe-Tn)
                                                 +nust[zk][yj][xi][24]*(TiHe-TiO)
                                                 +nust[zk][yj][xi][25]*(TiHe-TiH))
                                    +2.0*neme/(ns[2]*ams[2])*nust[zk][yj][xi][2]*(TiHe-Te);

                //----- electron temperature equation
                //Ce=ele_cooling_rate(xx, Te, Tn, ne, i, j, k);

                ff[k][j][i].fx[19] = dxdt[k][j][i].fx[19]
                                    +2.0*ame*( nust[zk][yj][xi][4]*(Te-TiO)+nust[zk][yj][xi][1]/ams[1]*(Te-TiH)
                                              +nust[zk][yj][xi][2]/ams[2]*(Te-TiHe)+nust[zk][yj][xi][6]*(Te-Tn))
                                    -two3rd*ame*nust[zk][yj][xi][5]
                                     *((ui[0]-unr)*(ui[0]-unr)+(ui[1]-unr)*(ui[1]-unt)+(ui[2]-unp)*(ui[2]-unp));
                                    //+two3rd/ne*Ce;

                //----- neutral momentum equation
                Nn=uu[k][j][i].fx[10]; rhon=uu[k][j][i].fx[34];
                rhossum_nusq_rhon[0]=( rhos[0]*nust[zk][yj][xi][11]+rhos[3]*nust[zk][yj][xi][30]
                                      +rhos[4]*nust[zk][yj][xi][32]+rhos[5]*nust[zk][yj][xi][34]
                                      +rhos[6]*nust[zk][yj][xi][36])/rhon;
                rhossum_nusq_rhon[1]=rhos[1]*nust[zk][yj][xi][19]/rhon;
                rhossum_nusq_rhon[2]=rhos[2]*nust[zk][yj][xi][27]/rhon;

                ff[k][j][i].fx[27] = dxdt[k][j][i].fx[27]+rhossum_nusq_rhon[0]*(unr-u7)
                                                         +rhossum_nusq_rhon[1]*(unr-u10)
                                                         +rhossum_nusq_rhon[2]*(unr-u13);
                ff[k][j][i].fx[28] = dxdt[k][j][i].fx[28]+rhossum_nusq_rhon[0]*(unt-u8)
                                                         +rhossum_nusq_rhon[1]*(unt-u11)
                                                         +rhossum_nusq_rhon[2]*(unt-u14);
                ff[k][j][i].fx[29] = dxdt[k][j][i].fx[29]+rhossum_nusq_rhon[0]*(unp-u9)
                                                         +rhossum_nusq_rhon[1]*(unp-u12)
                                                         +rhossum_nusq_rhon[2]*(unp-u15);

                //neutral temperature equation
                nuis_ms = rhos[0]*ams[0]*nust[zk][yj][xi][12]+rhos[3]*ams[3]*nust[zk][yj][xi][31]
                         +rhos[4]*ams[4]*nust[zk][yj][xi][33]+rhos[5]*ams[5]*nust[zk][yj][xi][35]
                         +rhos[6]*ams[6]*nust[zk][yj][xi][37];
                nuis = rhos[0]*nust[zk][yj][xi][12]+rhos[3]*nust[zk][yj][xi][31]
                      +rhos[4]*nust[zk][yj][xi][33]+rhos[5]*nust[zk][yj][xi][35]
                      +rhos[6]*nust[zk][yj][xi][37];
                ff[k][j][i].fx[30] =( 2.0*( nuis*(Tn-TiO)+rhos[1]*nust[zk][yj][xi][20]*(Tn-TiH)
                                           +rhos[2]*nust[zk][yj][xi][28]*(Tn-TiHe)
                                           +neme*nust[zk][yj][xi][6]*(Tn-Te))
                                     -two3rd*( nuis_ms*ui_un_sq[0]
                                              +rhos[1]*ams[1]*nust[zk][yj][xi][20]*ui_un_sq[1]
                                              +rhos[2]*ams[2]*nust[zk][yj][xi][28]*ui_un_sq[2]))/Nn;

                //----------- B^r equation
                ff[k][j][i].fx[31]= dxdt[k][j][i].fx[31]+0.5*( (xx[k][jp][i].fx[36]-xx[k][jm][i].fx[36])/dth
                                                              -(xx[kp][j][i].fx[35]-xx[km][j][i].fx[35])/dph);

                //----------- B^theta equation
                ff[k][j][i].fx[32]= dxdt[k][j][i].fx[32]+0.5*( (xx[kp][j][i].fx[34]-xx[km][j][i].fx[34])/dph
                                                              -(xx[k][j][ip].fx[36]-xx[k][j][im].fx[36])/dr);

                //----------- B^phi equation
                ff[k][j][i].fx[33]= dxdt[k][j][i].fx[33]+0.5*( (xx[k][j][ip].fx[35]-xx[k][j][im].fx[35])/dr
                                                              -(xx[k][jp][i].fx[34]-xx[k][jm][i].fx[34])/dth);

                electric_field_vxB(xx, uu, ui, i, j, k, xi, yj, zk, Ec_VxB);

                //----------- E^r equation
                ff[k][j][i].fx[34]= xx[k][j][i].fx[34]
                                   -(Jiv11[zk][yj]*Ec_VxB[0]+Jiv21[zk][yj]*Ec_VxB[1]+Jiv31[yj]*Ec_VxB[2]);

                //----------- E^theta equation
                ff[k][j][i].fx[35]= xx[k][j][i].fx[35]
                                   -(Jiv12[zk][yj][xi]*Ec_VxB[0]+Jiv22[zk][yj][xi]*Ec_VxB[1]+Jiv32[yj][xi]*Ec_VxB[2]);

                //----------- E^phi equation
                ff[k][j][i].fx[36]= xx[k][j][i].fx[36]-(Jiv13[zk][yj][xi]*Ec_VxB[0]+Jiv23[zk][yj][xi]*Ec_VxB[1]);

                for (s=0; s<a4; s++) {
                    if (isnan(ff[k][j][i].fx[s]) || isinf(ff[k][j][i].fx[s])) {
                        cout<<"Stiff function is Nan or inf at ("<<i<<", "<<j<<", "<<k<<", "<<s
                            <<") in stiff_functions "<<ff[k][j][i].fx[s]<<endl;
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