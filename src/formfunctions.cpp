#include "funcdef.h"
#include "boundary.h"

#include <fstream>
#include <iostream>

using namespace std;

#include "electric_field.h"
#include "ele_cooling_rate.h"
#include "operators.h"
#include "neu_cooling_rate.h"

int formfunctions(SNES snes, Vec X, Vec F, void* ctx)
{
    Vec    localX, localXn, localU;
    Field  ***xx, ***xn, ***uu, ***ff;
    PetscInt xs, xm, ys, ym, zs, zm;
    AppCtx *params = (AppCtx*)ctx;
    DM     da;

    SNESGetDM(snes, &da);

    DMDAGetCorners(da, &xs ,&ys, &zs, &xm, &ym, &zm);

    DMGetLocalVector(da, &localX);
    DMGlobalToLocalBegin(da, X, INSERT_VALUES,localX);
    DMGlobalToLocalEnd(da, X, INSERT_VALUES, localX);
    DMDAVecGetArray(da, localX, &xx);

    DMGetLocalVector(da, &localU);
    DMGlobalToLocalBegin(da,params->U, INSERT_VALUES, localU);
    DMGlobalToLocalEnd(da,params->U, INSERT_VALUES, localU);
    DMDAVecGetArray(da, localU, &uu);

    DMGetLocalVector(da, &localXn);
    DMGlobalToLocalBegin(da,params->Xn, INSERT_VALUES, localXn);
    DMGlobalToLocalEnd(da,params->Xn, INSERT_VALUES, localXn);
    DMDAVecGetArray(da, localXn, &xn);

    DMDAVecGetArray(da, F, &ff);

    //functions(xx, xn, uu, xs, xm, ys, ym, zs, zm, params, ff);

    int      i, j, k, xi, yj, zk, s, ip, im, jp, jm, kp, km;
    double   Nn, rhon, Tn, ns[7], TiO, TiH, TiHe, TiO_n, TiH_n, TiHe_n, Te_n, Tn_n;    
    vector3D gradns[7], gradNq, gradTe, gradne, gradTs, gradNn, gradTn;
    vector3D grad_usr, grad_usth, grad_usphi, grad_unr, grad_untheta, grad_unphi;
    double   usr[3], usth[3], usph[3], usrr, ustht, usphp, unr, unth, unph, uer, ueth, ueph;
    double   div_ui[3], div_un, div_ue, divui, ne, Te, Ec_gradPe[3], Ec_VxB[3];
    double   jr, jtheta, jphi, Qeuv, Qphoto, Cn, Ce;
    double   ui[3], uir, uith, uiph, Br, Bt, Bp, neme, rhos[7], rhossum_nusq_rhon[3];
    double   usr_n[3], usth_n[3], usph_n[3], nsmore, uir_n, uith_n, uiph_n, unr_n, unth_n, unph_n;
    double   ui_un_sq[3], uOi_uHi_sq, uOi_uHei_sq, uHi_uHei_sq, nuis_ms, nuis;
    double   ui_un_sq_n[3], uOi_uHi_sq_n, uOi_uHei_sq_n, uHi_uHei_sq_n, Ps[14], Ls[14];
    double   dt_quart=0.25*dt;

    const double one3rd=1.0/3.0, two3rd=2.0/3.0;

    //if (xs+xm == a1) top_bc_vel(params, ys, ym, zs, zm);

    for (k = zs; k < zs+zm; k++) {
        zk=k-zs;

        for (j = ys; j< ys+ym; j++) {
            //set boundary conditions at j=0 and j=Nth
            if (j == 0) {
                for (i = xs; i < xs+xm; i++) {
                    for (s = 0; s < a4; s++) {
                        if (s == 8 || s == 11 || s == 14 || s == 28 || s ==32 || s == 35)
                            ff[k][0][i].fx[s]=xx[k][0][i].fx[s]+xn[(k+a3/2) % a3][1][i].fx[s];
                        else ff[k][0][i].fx[s]=xx[k][0][i].fx[s]-xn[(k+a3/2) % a3][1][i].fx[s];
                    }
                }

                continue;
            }
            else if (j == Nth) {
                for (i = xs; i < xs+xm; i++) {
                    for (s = 0; s < a4; s++) {
                        if (s == 8 || s == 11 || s == 14 || s == 28 || s ==32 || s == 35)
                            ff[k][Nth][i].fx[s]=xx[k][Nth][i].fx[s]+xn[(k+a3/2) % a3][Nthm][i].fx[s];
                        else ff[k][Nth][i].fx[s]=xx[k][Nth][i].fx[s]-xn[(k+a3/2) % a3][Nthm][i].fx[s];
                    }
                }

                continue;
            }

            yj=j-ys;

            for (i = xs; i < xs+xm; i++) {
                if (i == 0) {
                    lower_boundary_bc(xx, xn, ff, j, k); continue;
                }
                else if (i == Nr) {
                    upper_boundary_bc(xx, xn, ff, j, k, yj, zk); continue;
                }

                xi=i-xs;

                usr_n[0]=xn[k][j][i].fx[7]; usth_n[0]=xn[k][j][i].fx[8]; usph_n[0]=xn[k][j][i].fx[9];
                usr_n[1]=xn[k][j][i].fx[10]; usth_n[1]=xn[k][j][i].fx[11]; usph_n[1]=xn[k][j][i].fx[12];
                usr_n[2]=xn[k][j][i].fx[13]; usth_n[2]=xn[k][j][i].fx[14]; usph_n[2]=xn[k][j][i].fx[15];
                unr_n=xn[k][j][i].fx[27]; unth_n=xn[k][j][i].fx[28]; unph_n=xn[k][j][i].fx[29];

                div_ui[0]=divergence(xn, i, j, k, yj, xi, 7);
                div_ui[1]=divergence(xn, i, j, k, yj, xi, 10);
                div_ui[2]=divergence(xn, i, j, k, yj, xi, 13);
                div_un   =divergence(xn, i, j, k, yj, xi, 27);

                ne=uu[k][j][i].fx[6]; Te=xx[k][j][i].fx[19]; Te_n=xn[k][j][i].fx[19];

                //production and loss rates
                prod_loss_rates(xn, uu, i, j, k, zk, yj, xi, Qeuv, Qphoto, Ps, Ls);

                rhon=0.0;
                for (s = 0; s < sl; s++) {
                    ns[s] = xn[k][j][i].fx[s];
                    rhos[s] = ns[s]*ams[s];

                    rhon += exp(xn[k][j][i].fx[s+20])*ams[s];

                    gradns[s] = gradient(xn, i, j, k, yj, xi, s);
                    gradNq = gradient(xn, i, j, k, yj, xi, s+20);

                    if (s < 3) {
                        divui=div_ui[s]; usrr=usr_n[s]; ustht=usth_n[s]; usphp=usph_n[s];
                    }
                    else {
                        divui=div_ui[0]; usrr=usr_n[0]; ustht=usth_n[0]; usphp=usph_n[0];
                    }

                    ff[k][j][i].fx[s]= (1.0+dt*Ls[s])*xx[k][j][i].fx[s]-xn[k][j][i].fx[s]
                                      +dt*( divui*xn[k][j][i].fx[s]+usrr*gradns[s].r
                                           +ustht*gradns[s].t+usphp*gradns[s].p - Ps[s]);
                    ff[k][j][i].fx[s+20]= xx[k][j][i].fx[s+20] - xn[k][j][i].fx[s+20]
                                         +dt*( div_un + unr_n*gradNq.r + unth_n*gradNq.t + unph_n*gradNq.p
                                              -Ps[s+7]*exp(-xn[k][j][i].fx[s+20])+Ls[s+7]);
                }

                usr[0]=xx[k][j][i].fx[7]; usth[0]=xx[k][j][i].fx[10]; usph[0]=xx[k][j][i].fx[13];
                usr[1]=xx[k][j][i].fx[8]; usth[1]=xx[k][j][i].fx[11]; usph[1]=xx[k][j][i].fx[14];
                usr[2]=xx[k][j][i].fx[9]; usth[2]=xx[k][j][i].fx[12]; usph[2]=xx[k][j][i].fx[15];
                unr=xx[k][j][i].fx[27]; unth=xx[k][j][i].fx[28]; unph=xx[k][j][i].fx[29];

                //ion bulk velocity
                nsmore=ns[0]+ns[3]+ns[4]+ns[5]+ns[6];
                uir_n=(nsmore*usr_n[0]+ns[1]*usr_n[1]+ns[2]*usr_n[2])/ne;
                uith_n=(nsmore*usth_n[0]+ns[1]*usth_n[1]+ns[2]*usth_n[2])/ne;
                uiph_n=(nsmore*usph_n[0]+ns[1]*usph_n[1]+ns[2]*usph_n[2])/ne;

                uir=(nsmore*usr[0]+ns[1]*usr[1]+ns[2]*usr[2])/ne;
                uith=(nsmore*usth[0]+ns[1]*usth[1]+ns[2]*usth[2])/ne;
                uiph=(nsmore*usph[0]+ns[1]*usph[1]+ns[2]*usph[2])/ne;

                jr=uu[k][j][i].fx[16]; jtheta=uu[k][j][i].fx[17]; jphi=uu[k][j][i].fx[18];

                //background magnetic field in normal spherical components
                Br=uu[k][j][i].fx[0]; Bt=uu[k][j][i].fx[1]; Bp=uu[k][j][i].fx[2];

                gradTe=gradient(xn, i, j, k, yj, xi, 19);
                gradne=gradient(uu, i, j, k, yj, xi, 6);
                neme=ne*ame;

                TiO=xx[k][j][i].fx[16]; TiH=xx[k][j][i].fx[17]; TiHe=xx[k][j][i].fx[18];
                Tn =xx[k][j][i].fx[30];

                ui_un_sq[0]= (usr[0]-unr)*(usr[0]-unr)+(usth[0]-unth)*(usth[0]-unth)
                            +(usph[0]-unph)*(usph[0]-unph);
                ui_un_sq[1]= (usr[1]-unr)*(usr[1]-unr)+(usth[1]-unth)*(usth[1]-unth)
                            +(usph[1]-unph)*(usph[1]-unph);
                ui_un_sq[2]= (usr[2]-unr)*(usr[2]-unr)+(usth[2]-unth)*(usth[2]-unth)
                            +(usph[2]-unph)*(usph[2]-unph);
                uOi_uHi_sq = (usr[0]-usr[1])*(usr[0]-usr[1])+(usth[0]-usth[1])*(usth[0]-usth[1])
                            +(usph[0]-usph[1])*(usph[0]-usph[1]);
                uOi_uHei_sq= (usr[0]-usr[2])*(usr[0]-usr[2])+(usth[0]-usth[2])*(usth[0]-usth[2])
                            +(usph[0]-usph[2])*(usph[0]-usph[2]);
                uHi_uHei_sq= (usr[1]-usr[2])*(usr[1]-usr[2])+(usth[1]-usth[2])*(usth[1]-usth[2])
                            +(usph[1]-usph[2])*(usph[1]-usph[2]);

                TiO_n=xn[k][j][i].fx[16]; TiH_n=xn[k][j][i].fx[17]; TiHe_n=xn[k][j][i].fx[18];
                Tn_n =xn[k][j][i].fx[30];

                ui_un_sq_n[0]= (usr_n[0]-unr_n)*(usr_n[0]-unr_n)
                              +(usth_n[0]-unth_n)*(usth_n[0]-unth_n)
                              +(usph_n[0]-unph_n)*(usph_n[0]-unph_n);
                ui_un_sq_n[1]= (usr[1]-unr_n)*(usr[1]-unr_n)
                              +(usth_n[1]-unth_n)*(usth_n[1]-unth_n)
                              +(usph_n[1]-unph_n)*(usph_n[1]-unph_n);
                ui_un_sq_n[2]= (usr_n[2]-unr_n)*(usr_n[2]-unr_n)
                              +(usth_n[2]-unth_n)*(usth_n[2]-unth_n)
                              +(usph_n[2]-unph_n)*(usph_n[2]-unph_n);
                uOi_uHi_sq_n = (usr_n[0]-usr_n[1])*(usr_n[0]-usr_n[1])
                              +(usth_n[0]-usth_n[1])*(usth_n[0]-usth_n[1])
                              +(usph_n[0]-usph_n[1])*(usph_n[0]-usph_n[1]);
                uOi_uHei_sq_n= (usr_n[0]-usr_n[2])*(usr_n[0]-usr_n[2])
                              +(usth_n[0]-usth_n[2])*(usth_n[0]-usth_n[2])
                              +(usph_n[0]-usph_n[2])*(usph_n[0]-usph_n[2]);
                uHi_uHei_sq_n= (usr_n[1]-usr_n[2])*(usr_n[1]-usr_n[2])
                              +(usth_n[1]-usth_n[2])*(usth_n[1]-usth_n[2])
                              +(usph_n[1]-usph_n[2])*(usph_n[1]-usph_n[2]);

//------------ for O+ ---------------------------------------------------------------------------------------
                grad_usr=gradient(xn, i, j, k, yj, xi, 7);
                grad_usth=gradient(xn, i, j, k, yj, xi, 8);
                grad_usphi=gradient(xn, i, j, k, yj, xi, 9);
                gradTs = gradient(xn, i, j, k, yj, xi, 16);

                //----------- uO+_r equation
                ff[k][j][i].fx[7]= xx[k][j][i].fx[7] - xn[k][j][i].fx[7]
                                  +dt_half*( qms[0]*((uith-usth[0])*Bp-(uiph-usph[0])*Bt)
                                            +nust[zk][yj][xi][11]*(usr[0]-unr)
                                            +nust[zk][yj][xi][7]*(usr[0]-usr[1])
                                            +nust[zk][yj][xi][8]*(usr[0]-usr[2])   //implicit part C-N scheme
                                            +qms[0]*((uith_n-usth_n[0])*Bp-(uiph_n-usph_n[0])*Bt)
                                            +nust[zk][yj][xi][11]*(usr_n[0]-unr_n)
                                            +nust[zk][yj][xi][7]*(usr_n[0]-usr_n[1])
                                            +nust[zk][yj][xi][8]*(usr_n[0]-usr_n[2]))
                                       +dt*( (usr_n[0]*grad_usr.r+usth_n[0]*grad_usr.t+usph_n[0]*grad_usr.p)
                                            +(gradTs.r+TiO_n/ns[0]*gradns[0].r+gradTe.r+Te_n/ne*gradne.r)/ams[0]
                                            -(jtheta*Bp-jphi*Bt)/rhos[0] + gr[xi]
                                            -(usth_n[0]*usth_n[0]+usph_n[0]*usph_n[0])/rr[i]);

                //----------- uO+_{theta} equation
                ff[k][j][i].fx[8]= xx[k][j][i].fx[8] - xn[k][j][i].fx[8]
                                  +dt_half*( qms[0]*((uiph-usph[0])*Br-(uir-usr[0])*Bp)
                                            +nust[zk][yj][xi][11]*(usth[0]-unth)
                                            +nust[zk][yj][xi][7]*(usth[0]-usth[1])
                                            +nust[zk][yj][xi][8]*(usth[0]-usth[2])
                                            +qms[0]*((uiph_n-usph_n[0])*Br-(uir_n-usr_n[0])*Bp)
                                            +nust[zk][yj][xi][11]*(usth_n[0]-unth_n)
                                            +nust[zk][yj][xi][7]*(usth_n[0]-usth_n[1])
                                            +nust[zk][yj][xi][8]*(usth_n[0]-usth_n[2]))
                                   +dt*( (usr_n[0]*grad_usth.r+usth_n[0]*grad_usth.t+usph_n[0]*grad_usth.p)
                                        +(gradTs.t+TiO_n/ns[0]*gradns[0].t+gradTe.t+Te_n/ne*gradne.t)/ams[0]
                                        -(jphi*Br-jr*Bp)/rhos[0]
                                        +usth_n[0]*usr_n[0]/rr[i]-usph_n[0]*usph_n[0]*cot_div_r[yj][xi]);

                //----------- uO+_{phi} equation
                ff[k][j][i].fx[9]= xx[k][j][i].fx[9] - xn[k][j][i].fx[9]
                                  +dt_half*( qms[0]*((uir-usr[0])*Bt-(uith-usth[0])*Br)
                                            +nust[zk][yj][xi][11]*(usph[0]-unph)
                                            +nust[zk][yj][xi][7]*(usph[0]-usph[1])
                                            +nust[zk][yj][xi][8]*(usph[0]-usph[2])
                                            +qms[0]*((uir_n-usr_n[0])*Bt-(uith_n-usth_n[0])*Br)
                                            +nust[zk][yj][xi][11]*(usph_n[0]-unph_n)
                                            +nust[zk][yj][xi][7]*(usph_n[0]-usph_n[1])
                                            +nust[zk][yj][xi][8]*(usph_n[0]-usph_n[2]))
                                  +dt*( (usr_n[0]*grad_usphi.r+usth_n[0]*grad_usphi.t+usph_n[0]*grad_usphi.p)
                                       +(gradTs.p+TiO_n/ns[0]*gradns[0].p+gradTe.p+Te_n/ne*gradne.p)/ams[0]
                                       -(jr*Bt-jtheta*Br)/rhos[0]
                                       +usph_n[0]*usr_n[0]/rr[i]+usth_n[0]*usph_n[0]*cot_div_r[yj][xi]);

                //----------- TO+ equation
                ff[k][j][i].fx[16]= xx[k][j][i].fx[16] - xn[k][j][i].fx[16]
                                   +dt*( ams[0]*( nust[zk][yj][xi][12]*(TiO-Tn)+nust[zk][yj][xi][9]*(TiO-TiH)
                                                 +nust[zk][yj][xi][10]*(TiO-TiHe))
                                        +neme/rhos[0]*nust[zk][yj][xi][0]*(TiO-Te)
                                        -one3rd*ams[0]*( nust[zk][yj][xi][13]*ui_un_sq[0]
                                                        +ams[1]*nust[zk][yj][xi][9] *uOi_uHi_sq
                                                        +ams[2]*nust[zk][yj][xi][10]*uOi_uHei_sq) //implicit part
                                        +ams[0]*( nust[zk][yj][xi][12]*(TiO_n-Tn_n)
                                                 +nust[zk][yj][xi][9] *(TiO_n-TiH_n)
                                                 +nust[zk][yj][xi][10]*(TiO_n-TiHe_n))
                                        +neme/rhos[0]*nust[zk][yj][xi][0]*(TiO_n-Te_n)
                                        -one3rd*ams[0]*( nust[zk][yj][xi][13]*ui_un_sq_n[0]
                                                        +ams[1]*nust[zk][yj][xi][9] *uOi_uHi_sq_n
                                                        +ams[2]*nust[zk][yj][xi][10]*uOi_uHei_sq_n)
                                        +(usr_n[0]*gradTs.r + usth_n[0]*gradTs.t + usph_n[0]*gradTs.p)
                                        +two3rd*(TiO_n*div_ui[0] + uu[k][j][i].fx[11]/ns[0]));

//------------ for H+ ---------------------------------------------------------------------------------------
                grad_usr=gradient(xn, i, j, k, yj, xi, 10);
                grad_usth=gradient(xn, i, j, k, yj, xi, 11);
                grad_usphi=gradient(xn, i, j, k, yj, xi, 12);
                gradTs = gradient(xn, i, j, k, yj, xi, 17);

                //----------- uH+_r equation
                ff[k][j][i].fx[10]= xx[k][j][i].fx[10] - xn[k][j][i].fx[10]
                                   +dt_half*( qms[1]*((uith-usth[1])*Bp-(uiph-usph[1])*Bt)
                                             +nust[zk][yj][xi][19]*(usr[1]-unr)
                                             +nust[zk][yj][xi][14]*(usr[1]-usr[0])
                                             +nust[zk][yj][xi][15]*(usr[1]-usr[2])   //implicit part C-N scheme
                                             +qms[1]*((uith_n-usth_n[1])*Bp-(uiph_n-usph_n[1])*Bt)
                                             +nust[zk][yj][xi][19]*(usr_n[1]-unr_n)
                                             +nust[zk][yj][xi][14]*(usr_n[1]-usr_n[0])
                                             +nust[zk][yj][xi][15]*(usr_n[1]-usr_n[2]))
                                   +dt*( (usr_n[1]*grad_usr.r+usth_n[1]*grad_usr.t+usph_n[1]*grad_usr.p)
                                        +(gradTs.r+TiH_n/ns[1]*gradns[1].r+gradTe.r+Te_n/ne*gradne.r)/ams[1]
                                        -(jtheta*Bp-jphi*Bt)/rhos[1] + gr[xi]
                                        -(usth_n[1]*usth_n[1]+usph_n[1]*usph_n[1])/rr[i]);

                //----------- uH+_{theta} equation
                ff[k][j][i].fx[11]= xx[k][j][i].fx[11] - xn[k][j][i].fx[11]
                                  +dt_half*( qms[1]*((uiph-usph[1])*Br-(uir-usr[1])*Bp)
                                            +nust[zk][yj][xi][19]*(usth[1]-unth)
                                            +nust[zk][yj][xi][14]*(usth[1]-usth[0])
                                            +nust[zk][yj][xi][15]*(usth[1]-usth[2])
                                            +qms[1]*((uiph_n-usph_n[1])*Br-(uir_n-usr_n[1])*Bp)
                                            +nust[zk][yj][xi][19]*(usth_n[1]-unth_n)
                                            +nust[zk][yj][xi][14]*(usth_n[1]-usth_n[0])
                                            +nust[zk][yj][xi][15]*(usth_n[1]-usth_n[2]))
                                  +dt*( (usr_n[1]*grad_usth.r+usth_n[1]*grad_usth.t+usph_n[1]*grad_usth.p)
                                        +(gradTs.t+TiH_n/ns[1]*gradns[1].t+gradTe.t+Te_n/ne*gradne.t)/ams[1]
                                        -(jphi*Br-jr*Bp)/rhos[1]
                                        +usth_n[1]*usr_n[1]/rr[i]-usph_n[1]*usph_n[1]*cot_div_r[yj][xi]);

                //----------- uH+_{phi} equation
                ff[k][j][i].fx[12]= xx[k][j][i].fx[12] - xn[k][j][i].fx[12]
                                  +dt_half*( qms[1]*((uir-usr[1])*Bt-(uith-usth[1])*Br)
                                            +nust[zk][yj][xi][19]*(usph[1]-unph)
                                            +nust[zk][yj][xi][14]*(usph[1]-usph[0])
                                            +nust[zk][yj][xi][15]*(usph[1]-usph[2])
                                            +qms[1]*((uir_n-usr_n[1])*Bt-(uith_n-usth_n[1])*Br)
                                            +nust[zk][yj][xi][19]*(usph_n[1]-unph_n)
                                            +nust[zk][yj][xi][14]*(usph_n[1]-usph_n[0])
                                            +nust[zk][yj][xi][15]*(usph_n[1]-usph_n[2]))
                                  +dt*( (usr_n[1]*grad_usphi.r+usth_n[1]*grad_usphi.t+usph_n[1]*grad_usphi.p)
                                       +(gradTs.p+TiH_n/ns[1]*gradns[1].p+gradTe.p+Te_n/ne*gradne.p)/ams[1]
                                       -(jr*Bt-jtheta*Br)/rhos[1]
                                       +usph_n[1]*usr_n[1]/rr[i]+usth_n[1]*usph_n[1]*cot_div_r[yj][xi]);

                //----------- TH+ equation
                ff[k][j][i].fx[17]= xx[k][j][i].fx[17] - xn[k][j][i].fx[17]
                                   +dt*( ams[1]*( nust[zk][yj][xi][20]*(TiH-Tn)
                                                 +nust[zk][yj][xi][16]*(TiH-TiO)
                                                 +nust[zk][yj][xi][17]*(TiH-TiHe))
                                        +neme/rhos[1]*nust[zk][yj][xi][1]*(TiH-Te)
                                        -one3rd*ams[1]*( nust[zk][yj][xi][21]*ui_un_sq[1]
                                                        +nust[zk][yj][xi][18]*uOi_uHi_sq
                                                        +ams[2]*nust[zk][yj][xi][17]*uHi_uHei_sq) //implicit part
                                        +ams[1]*( nust[zk][yj][xi][20]*(TiH_n-Tn_n)
                                                 +nust[zk][yj][xi][16]*(TiH_n-TiO_n)
                                                 +nust[zk][yj][xi][17]*(TiH_n-TiHe_n))
                                        +neme/rhos[1]*nust[zk][yj][xi][1]*(TiH_n-Te_n)
                                        -one3rd*ams[1]*( nust[zk][yj][xi][21]*ui_un_sq_n[1]
                                                        +nust[zk][yj][xi][18]*uOi_uHi_sq_n
                                                        +ams[2]*nust[zk][yj][xi][17]*uHi_uHei_sq_n)
                                        +(usr_n[1]*gradTs.r + usth_n[1]*gradTs.t + usph_n[1]*gradTs.p)
                                        +two3rd*(TiH_n*div_ui[1] + uu[k][j][i].fx[12]/ns[1]));

//------------ for He+ ---------------------------------------------------------------------------------------
                grad_usr=gradient(xn, i, j, k, yj, xi, 13);
                grad_usth=gradient(xn, i, j, k, yj, xi, 14);
                grad_usphi=gradient(xn, i, j, k, yj, xi, 15);
                gradTs = gradient(xn, i, j, k, yj, xi, 18);

                //----------- uHe+_r equation
                ff[k][j][i].fx[13]= xx[k][j][i].fx[13] - xn[k][j][i].fx[13]
                                   +dt_half*( qms[2]*((uith-usth[2])*Bp-(uiph-usph[2])*Bt)
                                            +nust[zk][yj][xi][27]*(usr[2]-unr)
                                            +nust[zk][yj][xi][22]*(usr[2]-usr[0])
                                            +nust[zk][yj][xi][23]*(usr[2]-usr[1])   //implicit part C-N scheme
                                            +qms[2]*((uith_n-usth_n[2])*Bp-(uiph_n-usph_n[2])*Bt)
                                            +nust[zk][yj][xi][27]*(usr_n[2]-unr_n)
                                            +nust[zk][yj][xi][22]*(usr_n[2]-usr_n[0])
                                            +nust[zk][yj][xi][23]*(usr_n[2]-usr_n[1]))
                                       +dt*( (usr_n[2]*grad_usr.r+usth_n[2]*grad_usr.t+usph_n[2]*grad_usr.p)
                                            +(gradTs.r+TiHe_n/ns[2]*gradns[2].r+gradTe.r+Te_n/ne*gradne.r)/ams[2]
                                            -(jtheta*Bp-jphi*Bt)/rhos[2] + gr[xi]
                                            -(usth_n[2]*usth_n[2]+usph_n[2]*usph_n[2])/rr[i]);

                //----------- uHe+_{theta} equation
                ff[k][j][i].fx[14]= xx[k][j][i].fx[14] - xn[k][j][i].fx[14]
                                   +dt_half*( qms[2]*((uiph-usph[2])*Br-(uir-usr[2])*Bp)
                                             +nust[zk][yj][xi][27]*(usth[2]-unth)
                                             +nust[zk][yj][xi][22]*(usth[2]-usth[0])
                                             +nust[zk][yj][xi][23]*(usth[2]-usth[1])
                                             +qms[2]*((uiph_n-usph_n[2])*Br-(uir_n-usr_n[2])*Bp)
                                             +nust[zk][yj][xi][27]*(usth_n[2]-unth_n)
                                             +nust[zk][yj][xi][22]*(usth_n[2]-usth_n[0])
                                             +nust[zk][yj][xi][23]*(usth_n[2]-usth_n[1]))
                                   +dt*( (usr_n[2]*grad_usth.r+usth_n[2]*grad_usth.t+usph_n[2]*grad_usth.p)
                                        +(gradTs.t+TiHe_n/ns[2]*gradns[2].t+gradTe.t+Te_n/ne*gradne.t)/ams[2]
                                        -(jphi*Br-jr*Bp)/rhos[2]
                                        +usth_n[2]*usr_n[2]/rr[i]-usph_n[2]*usph_n[2]*cot_div_r[yj][xi]);

                //----------- uHe+_{phi} equation
                ff[k][j][i].fx[15]= xx[k][j][i].fx[15] - xn[k][j][i].fx[15]
                                  +dt_half*( qms[2]*((uir-usr[2])*Bt-(uith-usth[2])*Br)
                                            +nust[zk][yj][xi][27]*(usph[2]-unph)
                                            +nust[zk][yj][xi][22]*(usph[2]-usph[0])
                                            +nust[zk][yj][xi][23]*(usph[2]-usph[1])
                                            +qms[2]*((uir_n-usr_n[2])*Bt-(uith_n-usth_n[2])*Br)
                                            +nust[zk][yj][xi][27]*(usph_n[2]-unph_n)
                                            +nust[zk][yj][xi][22]*(usph_n[2]-usph_n[0])
                                            +nust[zk][yj][xi][23]*(usph_n[2]-usph_n[1]))
                                  +dt*( (usr_n[2]*grad_usphi.r+usth_n[2]*grad_usphi.t+usph_n[2]*grad_usphi.p)
                                       +(gradTs.p+TiHe_n/ns[2]*gradns[2].p+gradTe.p+Te_n/ne*gradne.p)/ams[2]
                                       -(jr*Bt-jtheta*Br)/rhos[2]
                                       +usph_n[2]*usr_n[2]/rr[i]+usth_n[2]*usph_n[2]*cot_div_r[yj][xi]);

                //----------- THe+ equation
                ff[k][j][i].fx[18]= xx[k][j][i].fx[18] - xn[k][j][i].fx[18]
                                   +dt*( ams[2]*( nust[zk][yj][xi][28]*(TiHe-Tn)
                                                 +nust[zk][yj][xi][24]*(TiHe-TiO)
                                                 +nust[zk][yj][xi][25]*(TiHe-TiH))
                                        +neme/rhos[2]*nust[zk][yj][xi][2]*(TiHe-Te)
                                        -one3rd*ams[2]*( nust[zk][yj][xi][29]*ui_un_sq[2]
                                                        +nust[zk][yj][xi][26]*uOi_uHei_sq
                                                        +ams[1]*nust[zk][yj][xi][25]*uHi_uHei_sq) //implicit part
                                        +ams[2]*( nust[zk][yj][xi][28]*(TiH_n-Tn_n)
                                        +nust[zk][yj][xi][24]*(TiH_n-TiO_n)
                                                 +nust[zk][yj][xi][25]*(TiH_n-TiHe_n))
                                        +neme/rhos[2]*nust[zk][yj][xi][2]*(TiH_n-Te_n)
                                        -one3rd*ams[2]*( nust[zk][yj][xi][29]*ui_un_sq_n[2]
                                                        +nust[zk][yj][xi][26]*uOi_uHei_sq_n
                                                        +ams[1]*nust[zk][yj][xi][25]*uHi_uHei_sq_n)
                                        +(usr_n[2]*gradTs.r + usth_n[2]*gradTs.t + usph_n[2]*gradTs.p)
                                        +two3rd*(TiHe_n*div_ui[2] + uu[k][j][i].fx[13]/ns[2]));

//----------- Te equation --------------------------------------------------------------------------------
                uer = uu[k][j][i].fx[7]; ueth = uu[k][j][i].fx[8]; ueph = uu[k][j][i].fx[9];
                div_ue  = divergence(uu, i, j, k, yj, xi, 7);

                Ce=ele_cooling_rate(xn, Te_n, Tn_n, ne, i, j, k);

                ff[k][j][i].fx[19]= xx[k][j][i].fx[19] - xn[k][j][i].fx[19]
                                   +dt*( ame*( nust[zk][yj][xi][4]*(Te-TiO)
                                              +nust[zk][yj][xi][1]/ams[1]*(Te-TiH)
                                              +nust[zk][yj][xi][2]/ams[2]*(Te-TiHe)
                                              +nust[zk][yj][xi][6]*(Te-Tn))
                                        -one3rd*ame*nust[zk][yj][xi][5]
                                         *( (uir-unr)*(uir-unr)+(uith-unth)*(uith-unth)
                                           +(uiph-unph)*(uiph-unph))
                                        +ame*( nust[zk][yj][xi][4]*(Te_n-TiO_n)
                                              +nust[zk][yj][xi][1]/ams[1]*(Te_n-TiH_n)
                                              +nust[zk][yj][xi][2]/ams[2]*(Te_n-TiHe_n)
                                              +nust[zk][yj][xi][6]*(Te_n-Tn_n))
                                        -one3rd*ame*nust[zk][yj][xi][5]
                                         *( (uir_n-unr_n)*(uir_n-unr_n)
                                           +(uith_n-unth_n)*(uith_n-unth_n)
                                           +(uiph_n-unph_n)*(uiph_n-unph_n))
                                        +(uer*gradTe.r + ueth*gradTe.t + ueph*gradTe.p)
                                        +two3rd*(Te_n*div_ue+(uu[k][j][i].fx[14]-Qeuv + Ce)/ne));

//------------------ for neutrals
                //-----------  un_r equation
                Nn=uu[k][j][i].fx[10];

                gradTn=gradient(xn, i, j, k, yj, xi, 30);
                gradNn=gradient(uu, i, j, k, yj, xi, 10);
                grad_unr=gradient(xn, i, j, k, yj, xi, 27);
                grad_untheta=gradient(xn, i, j, k, yj, xi, 28);
                grad_unphi=gradient(xn, i, j, k, yj, xi, 29);

                rhossum_nusq_rhon[0]=( rhos[0]*nust[zk][yj][xi][11]+rhos[3]*nust[zk][yj][xi][30]
                                             +rhos[4]*nust[zk][yj][xi][32]+rhos[5]*nust[zk][yj][xi][34]
                                             +rhos[6]*nust[zk][yj][xi][36])/rhon;
                rhossum_nusq_rhon[1]=rhos[1]*nust[zk][yj][xi][19]/rhon;
                rhossum_nusq_rhon[2]=rhos[2]*nust[zk][yj][xi][27]/rhon;

                ff[k][j][i].fx[27]= xx[k][j][i].fx[27] - xn[k][j][i].fx[27]
                                   +dt_half*( rhossum_nusq_rhon[0]*(unr-usr[0])
                                             +rhossum_nusq_rhon[1]*(unr-usr[1])
                                             +rhossum_nusq_rhon[2]*(unr-usr[2])  //implicit part
                                             +rhossum_nusq_rhon[0]*(unr_n-usr_n[0])
                                             +rhossum_nusq_rhon[1]*(unr_n-usr_n[1])
                                             +rhossum_nusq_rhon[2]*(unr_n-usr_n[2]))
                                   +dt*( (unr_n*grad_unr.r + unth_n*grad_unr.t + unph_n*grad_unr.p)
                                        -(unth_n*unth_n+unph_n*unph_n)/rr[i] + gr[xi]
                                        +(Nn*gradTn.r + Tn_n*gradNn.r)/rhon
                                       +2.0*(rotat_t[zk][yj]*unph_n-rotat_p[zk][yj]*unth_n)
                                       +cenf_r[zk][yj][xi]);

                //----------- un_{theta} equation
                ff[k][j][i].fx[28]= xx[k][j][i].fx[28] - xn[k][j][i].fx[28]
                                   +dt_half*( rhossum_nusq_rhon[0]*(unth-usth[0])
                                             +rhossum_nusq_rhon[1]*(unth-usth[1])
                                             +rhossum_nusq_rhon[2]*(unth-usth[2])  //implicit part
                                             +rhossum_nusq_rhon[0]*(unth_n-usth_n[0])
                                             +rhossum_nusq_rhon[1]*(unth_n-usth_n[1])
                                             +rhossum_nusq_rhon[2]*(unth_n-usth_n[2]))
                                   +dt*( (unr_n*grad_untheta.r+unth_n*grad_untheta.t+unph_n*grad_untheta.p)
                                        -unph_n*unph_n*cot_div_r[yj][xi]+unth_n*unr_n/rr[i]
                                        +(Nn*gradTn.t + Tn_n*gradNn.t)/rhon
                                        +2.0*(rotat_p[zk][yj]*unr_n-rotat_r[zk][yj]*unph_n)
                                        +cenf_t[zk][yj][xi]);

                //----------- un_{phi} equation
                ff[k][j][i].fx[29]= xx[k][j][i].fx[29] - xn[k][j][i].fx[29]
                                   +dt_half*( rhossum_nusq_rhon[0]*(unph-usph[0])
                                             +rhossum_nusq_rhon[1]*(unph-usph[1])
                                             +rhossum_nusq_rhon[2]*(unph-usph[2])  //implicit part
                                             +rhossum_nusq_rhon[0]*(unph_n-usph_n[0])
                                             +rhossum_nusq_rhon[1]*(unph_n-usph_n[1])
                                             +rhossum_nusq_rhon[2]*(unph_n-usph_n[2]))
                                   +(unr_n*grad_unphi.r + unth_n*grad_unphi.t + unph_n*grad_unphi.p)
                                   +unth_n*unph_n*cot_div_r[yj][xi] + unph_n*unr_n/rr[i]
                                   +(Nn*gradTn.p + Tn_n*gradNn.p)/rhon
                                   +2.0*(rotat_r[zk][yj]*unth_n-rotat_t[zk][yj]*unr_n)+cenf_p[zk][yj][xi];

                //----------- non-stiff terms of Tn equation
                Cn=neu_cooling_rate(xn, i, j, k);

                nuis_ms = rhos[0]*ams[0]*nust[zk][yj][xi][12]+rhos[3]*ams[3]*nust[zk][yj][xi][31]
                         +rhos[4]*ams[4]*nust[zk][yj][xi][33]+rhos[5]*ams[5]*nust[zk][yj][xi][35]
                         +rhos[6]*ams[6]*nust[zk][yj][xi][37];
                nuis = rhos[0]*nust[zk][yj][xi][12]+rhos[3]*nust[zk][yj][xi][31]
                      +rhos[4]*nust[zk][yj][xi][33]+rhos[5]*nust[zk][yj][xi][35]
                      +rhos[6]*nust[zk][yj][xi][37];

                ff[k][j][i].fx[30]= xx[k][j][i].fx[30] - xn[k][j][i].fx[30]
                                   +dt*( ( nuis*(Tn-TiO)+rhos[1]*nust[zk][yj][xi][20]*(Tn-TiH)
                                        +rhos[2]*nust[zk][yj][xi][28]*(Tn-TiHe)
                                        +neme*nust[zk][yj][xi][6]*(Tn-Te))
                                        -one3rd*( nuis_ms*ui_un_sq[0]
                                                 +rhos[1]*ams[1]*nust[zk][yj][xi][20]*ui_un_sq[1]
                                                 +rhos[2]*ams[2]*nust[zk][yj][xi][28]*ui_un_sq[2]))/Nn 
                                   +dt*( ( nuis*(Tn_n-TiO_n)+rhos[1]*nust[zk][yj][xi][20]*(Tn_n-TiH_n)
                                        +rhos[2]*nust[zk][yj][xi][28]*(Tn_n-TiHe_n)
                                        +neme*nust[zk][yj][xi][6]*(Tn_n-Te_n))
                                        -one3rd*( nuis_ms*ui_un_sq_n[0]
                                                 +rhos[1]*ams[1]*nust[zk][yj][xi][20]*ui_un_sq_n[1]
                                                 +rhos[2]*ams[2]*nust[zk][yj][xi][28]*ui_un_sq_n[2]))/Nn
                                   +dt*( (unr_n*gradTn.r + unth_n*gradTn.t + unph_n*gradTn.p)
                                        +two3rd*(Tn_n*div_un + (uu[k][j][i].fx[15] -Qeuv+Cn)/Nn));

//----------- magnetic induction equation
                ip=i+1; im=i-1; jp=j+1; jm=j-1; kp=k+1; km=k-1;

                ff[k][j][i].fx[31] = xx[k][j][i].fx[31] - xn[k][j][i].fx[31]
                                    +dt_quart*( ( (xx[k][jp][i].fx[36]-xx[k][jm][i].fx[36])/dth
                                                 -(xx[kp][j][i].fx[35]-xx[km][j][i].fx[35])/dph)
                                               +( (xn[k][jp][i].fx[36]-xn[k][jm][i].fx[36])/dth
                                                 -(xn[kp][j][i].fx[35]-xn[km][j][i].fx[35])/dph));

                ff[k][j][i].fx[32] = xx[k][j][i].fx[32] - xn[k][j][i].fx[32]
                                    +dt_quart*( ( (xx[kp][j][i].fx[34]-xx[km][j][i].fx[34])/dph
                                                 -(xx[k][j][ip].fx[36]-xx[k][j][im].fx[36])/dr)
                                               +( (xn[kp][j][i].fx[34]-xn[km][j][i].fx[34])/dph
                                                 -(xn[k][j][ip].fx[36]-xn[k][j][im].fx[36])/dr));

                ff[k][j][i].fx[33] = xx[k][j][i].fx[33] - xn[k][j][i].fx[33]
                                    +dt_quart*( ( (xx[k][j][ip].fx[35]-xx[k][j][im].fx[35])/dr
                                                 -(xx[k][jp][i].fx[34]-xx[k][jm][i].fx[34])/dth)
                                               +( (xn[k][j][ip].fx[35]-xn[k][j][im].fx[35])/dr
                                                 -(xn[k][jp][i].fx[34]-xn[k][jm][i].fx[34])/dth));

                //----------- electric field equation
                ui[0]=uir; ui[1]=uith; ui[2]=uiph;
                electric_field_vxB(xx, uu, ui, i, j, k, xi, yj, zk, Ec_VxB);
                E_gradPe(xn, uu, i, j, k, xi, yj, zk, Ec_gradPe);

                ff[k][j][i].fx[34]= xx[k][j][i].fx[34] -xn[k][j][i].fx[34]
                                   -( Jiv11[zk][yj]*(Ec_VxB[0]-Ec_gradPe[0])
                                     +Jiv21[zk][yj]*(Ec_VxB[1]-Ec_gradPe[1])
                                     +Jiv31[yj]*(Ec_VxB[2]-Ec_gradPe[2]));

                ff[k][j][i].fx[35]= xx[k][j][i].fx[35]-xn[k][j][i].fx[35]
                                   -( Jiv12[zk][yj][xi]*(Ec_VxB[0]-Ec_gradPe[0])
                                     +Jiv22[zk][yj][xi]*(Ec_VxB[1]-Ec_gradPe[1])
                                     +Jiv32[yj][xi]*(Ec_VxB[2]-Ec_gradPe[2]));

                ff[k][j][i].fx[36]= xx[k][j][i].fx[36]-xn[k][j][i].fx[36]
                                   -( Jiv13[zk][yj][xi]*(Ec_VxB[0]-Ec_gradPe[0])
                                     +Jiv23[zk][yj][xi]*(Ec_VxB[1]-Ec_gradPe[1]));

                for (s=0; s<a4; s++) {
                    if (isnan(ff[k][j][i].fx[s]) || isinf(ff[k][j][i].fx[s])) {
                        std::cout<<"function is Nan or inf at (i, j, k, s) = ("<<i<<", "<<j<<", "<<k<<", "<<s
                            <<") in formfunctions "<<ff[k][j][i].fx[s]<<endl;
                        exit(-1);
                    }
                }
            }
        }
    }

/*#include <iostream>
#include <fstream>

fstream fstr;
    fstr.open("run0.log", fstream::out);

    for (k = zs; k < zs+zm; k++) {
        for (j = ys; j< ys+ym; j++) {
            for (i = xs; i < xs+xm; i++) {
                for (s=0; s< a4; s++) fstr<<k<<" "<<j<<" "<<i<<" "<<s
                <<" "<<ff[k][j][i].fx[s]<<" "<<xx[k][j][i].fx[s]<<" "<<xn[k][j][i].fx[s]<<endl;
            }
        }
    }
    fstr.close();*/

    DMDAVecRestoreArray(da, localX, &xx);
    DMDAVecRestoreArray(da, localXn, &xn);
    DMDAVecRestoreArray(da, localU, &uu);

    DMRestoreLocalVector(da, &localX);
    DMRestoreLocalVector(da, &localXn);
    DMRestoreLocalVector(da, &localU);

    DMDAVecRestoreArray(da, F, &ff);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout<<"------- OK --------------"<<endl;

    return 0;
}
