#include "funcdef.h"
#include "boundary.h"

#include <fstream>
#include <iostream>

using namespace std;

#include "electric_field.h"
#include "ele_cooling_rate.h"
#include "operators.h"

int rhsfunctions(TS ts, double ftime, Vec X, Vec G, void* ctx)
{
    Vec    localX, localU;
    Field  ***xx, ***uu, ***gg;
    PetscInt xs, xm, ys, ym, zs, zm;
    AppCtx *params = (AppCtx*)ctx;
    DM     da;

    TSGetDM(ts, &da);

    params->sec = ftime;

    DMDAGetCorners(da, &xs ,&ys, &zs, &xm, &ym, &zm);

    DMGetLocalVector(da, &localX);
    DMGlobalToLocalBegin(da, X, INSERT_VALUES,localX);
    DMGlobalToLocalEnd(da, X, INSERT_VALUES, localX);
    DMDAVecGetArray(da, localX, &xx);

    DMGetLocalVector(da, &localU);
    DMGlobalToLocalBegin(da,params->U,INSERT_VALUES,localU);
    DMGlobalToLocalEnd(da,params->U,INSERT_VALUES,localU);
    DMDAVecGetArray(da, localU, &uu);

    DMDAVecGetArray(da, G, &gg);

    int    i, j, k, xi, yj, zk, s, s7, s8, s9, s16;
    int    im, ip, jm, jp, km, kp, kcm, kcp, kji, kj, ji;

    double   Nn, rhon, Tn;    
    vector3D gradNs, gradNq, gradTe, gradne, gradTs, graNs, gradNn, gradTn;
    vector3D grad_usr, grad_ustheta, grad_usphi, grad_unr, grad_untheta, grad_unphi;
    double   dTi_dr, dni_dr, dTi_dth, dni_dth, dTi_dph, dni_dph, dTe_dr, dTe_dth, dTe_dph;
    double   duir_dr, duir_dth, duir_dph, duith_dr, duith_dth, duith_dph, duiph_dr, duiph_dth, duiph_dph;
    double   uir[3], uith[3], uiph[3], uirr, uitht, uiphp, unr, unth, unph, uer, ueth, ueph;
    double   dunr_dr, dunr_dth, dunr_dph, dunth_dr, dunth_dth, dunth_dph, dunph_dr, dunph_dth, dunph_dph;
    double   dTn_dr, dnn_dr, dTn_dth, dnn_dth, dTn_dph, dnn_dph;
    double   div_ui[3], div_un, div_ue, divui, ne, Te, Ns, Ts;
    double   jr, jtheta, jphi, Br0, Btheta0, Bphi0;
    double   ms_ne, ene;

    double   nuinOiHi_r, nuinOiHi_th, nuinOiHi_ph, nuinOiHei_r, nuinOiHei_th, nuinOiHei_ph;
    double   nuinHiHei_r, nuinHiHei_th, nuinHiHei_ph;
    double   nuinHi_part, nuinHei_part;

    vector3D efd_gradPe;

    const double two3rd=2.0/3.0;

    for (k = zs; k < zs+zm; k++) {
        zk=k-zs; km=k-1; kp = k+1;

        for (j = ys; j< ys+ym; j++) {

            //boundary conditions at j=0 and j=Nth
            if (j == 0) {
                for (i == xs; i < xs+xm; i++) {
                    for (s = 0; s < a4; s++) {
                        if (s == 8 || s == 11 || s == 14 || s == 28 || s ==32 || s == 35)
                            gg[k][0][i].fx[s]=-xx[(k+a3/2) % a3][1][i].fx[s];
                        else gg[k][0][i].fx[s]=xx[(k+a3/2) % a3][1][i].fx[s];
                    }
                }
                continue;
            }
            else if (j == Nth) {
                for (i == xs; i < xs+xm; i++) {
                    for (s = 0; s < a4; s++) {
                        if (s == 8 || s == 11 || s == 14 || s == 28 || s ==32 || s == 35)
                            gg[k][Nth][i].fx[s]=-xx[(k+a3/2) % a3][Nthm][i].fx[s];
                        else gg[k][Nth][i].fx[s]=xx[(k+a3/2) % a3][Nthm][i].fx[s];
                    }
                }
                continue;
            }

            yj=j-ys; jm=j-1; jp=j+1; kj=zk*ym+yj;

            for (i = xs; i < xs+xm; i++) {
                xi=i-xs; ip=i+1; yj=yj*xm+xi; kji=zk*ym*xm+yj*xm+xi;
                im = i-1; ip = i+1;

                if (i == 0 || i == Nr) continue;

                uir[0]=xx[k][j][i].fx[7]; uith[0]=xx[k][j][i].fx[8]; uiph[0]=xx[k][j][i].fx[9];
                uir[1]=xx[k][j][i].fx[10]; uith[1]=xx[k][j][i].fx[11]; uiph[1]=xx[k][j][i].fx[12];
                uir[2]=xx[k][j][i].fx[13]; uith[2]=xx[k][j][i].fx[14]; uiph[2]=xx[k][j][i].fx[15];

                div_ui[0]=divergence(xx, i, j, k, yj, xi, 7);
                div_ui[1]=divergence(xx, i, j, k, yj, xi, 10);
                div_ui[2]=divergence(xx, i, j, k, yj, xi, 13);
                div_un = divergence(xx, i, j, k, yj, xi, 27);

                gradTe=gradient(xx, i, j, k, yj, xi, 19);
                gradne=gradient(uu, i, j, k, yj, xi, 6);

                ne=uu[k][j][i].fx[6]; Te=xx[k][j][i].fx[19];

                jr=uu[k][j][i].fx[20]; jtheta=uu[k][j][i].fx[21]; jphi=uu[k][j][i].fx[22];
                Br0=uu[k][j][i].fx[0]; Btheta0=uu[k][j][i].fx[1]; Bphi0=uu[k][j][i].fx[2];

                for (s = 0; s < sl; s++) {
                    gradNs = gradient(xx, i, j, k, yj, xi, s);
                    gradNq = gradient(xx, i, j, k, yj, xi, s+20);

                    if (s < 3) {
                        divui=div_ui[s]; uirr=uir[s]; uitht=uith[s]; uiphp=uiph[s];
                    }
                    else {
                        divui=div_ui[0]; uirr=uir[0]; uitht=uith[0]; uiphp=uiph[0];
                    }

                    gg[k][j][i].fx[s]= Ps[zk][yj][xi][s]*1.0e6*exp(-xx[k][j][i].fx[s])-Ls[zk][yj][xi][s]
                                      -divui - uirr*gradNs.r-uitht*gradNs.t-uiphp*gradNs.p;
                    gg[k][j][i].fx[s+20]= Ps[zk][yj][xi][s+7]*1.0e6*exp(-xx[k][j][i].fx[s+20])-Ls[zk][yj][xi][s+7]
                                         -div_un - unr*gradNq.r-unth*gradNq.t-unph*gradNq.p;

                    if (s >= 3) continue;

                    Ns=xx[k][j][i].fx[s]; Ts=xx[k][j][i].fx[s+16]; ms_ne=ams[s]*ne;

                    //----------- non-stiff terms of uir equation
                    s7=7+s*3; s8=8+s*3; s9=9+s*3;
                    grad_usr=gradient(xx, i, j, k, yj, xi, s7);
                    grad_ustheta=gradient(xx, i, j, k, yj, xi, s8);
                    grad_usphi=gradient(xx, i, j, k, yj, xi, s9);

                    gradTs = gradient(xx, i, j, k, yj, xi, s16);
                    gradNs = gradient(uu, i, j, k, yj, xi, 28+s);

                    gg[k][j][i].fx[s7]= (jtheta*Bphi0-jphi*Btheta0)/ms_ne-gr[xi]+(uith[s]*uith[s]+uiph[s]*uiph[s])/rr[i]
                                       -(uir[s]*grad_usr.r + uith[s]*grad_usr.t + uiph[s]*grad_usr.p)
                                       -(gradTs.r + Ts*gradNs.r + gradTe.r + Te/ne*gradne.r)/ams[s];

                    //----------- non-stiff terms of uitheta equation
                    gg[k][j][i].fx[s8]= (jphi*Br0-jr*Bphi0)/ms_ne-uith[s]*uir[s]/rr[i]+uiph[s]*uiph[s]*cot_div_r[yj][xi]
                                       -(uir[s]*grad_ustheta.r + uith[s]*grad_ustheta.t + uiph[s]*grad_ustheta.p)
                                       -(gradTs.t + Ts*gradNs.t + gradTe.t + Te/ne*gradne.t)/ams[s];

                    //----------- non-stiff terms of uiphi equation
                    gg[k][j][i].fx[s9]= (jr*Btheta0-jtheta*Br0)/ms_ne-uiph[s]*uir[s]/rr[i]-uith[s]*uiph[s]*cot_div_r[yj][xi]
                                       -(uir[s]*grad_usphi.r + uith[s]*grad_usphi.t + uiph[s]*grad_usphi.p)
                                       -(gradTs.p + Ts*gradNs.p + gradTe.p + Te/ne*gradne.p)/ams[s];

                    //----------- non-stiff terms of Ti equation
                    gg[k][j][i].fx[s16]=-(uir[s]*gradTs.r + uith[s]*gradTs.t + uiph[s]*gradTs.p)
                                        -two3rd*(div_ui[s] + uu[k][j][i].fx[s+15]*exp(-xx[k][j][i].fx[s]));
                }

                nuinOiHi_r = nust[zk][yj][xi][15]*(uir[1] - uir[0]);
                nuinOiHi_th = nust[zk][yj][xi][15]*(uith[1] - uith[0]);
                nuinOiHi_ph = nust[zk][yj][xi][15]*(uiph[1] - uiph[0]);

                nuinOiHei_r = nust[zk][yj][xi][15]*(uir[2] - uir[0]);
                nuinOiHei_th = nust[zk][yj][xi][15]*(uith[2] - uith[0]);
                nuinOiHei_ph = nust[zk][yj][xi][15]*(uiph[2] - uiph[0]);

                gg[k][j][i].fx[7] += nuinOiHi_r + nuinOiHei_r;
                gg[k][j][i].fx[8] += nuinOiHi_th + nuinOiHei_th;
                gg[k][j][i].fx[9] += nuinOiHi_ph + nuinOiHei_ph;

                nuinHiHei_r = nust[zk][yj][xi][30]*(uir[2] - uir[1]);
                nuinHiHei_th = nust[zk][yj][xi][30]*(uith[2] - uith[1]);
                nuinHiHei_ph = nust[zk][yj][xi][30]*(uiph[2] - uiph[1]);

                nuinHi_part = nust[zk][yj][xi][31]+nust[zk][yj][xi][32]+nust[zk][yj][xi][33]+nust[zk][yj][xi][34];

                gg[k][j][i].fx[10] += nuinHiHei_r - nuinOiHi_r + nuinHi_part*(uir[0] - uir[1]);
                gg[k][j][i].fx[11] += nuinHiHei_th - nuinOiHi_th + nuinHi_part*(uith[0] - uith[1]);
                gg[k][j][i].fx[12] += nuinHiHei_ph - nuinOiHi_ph + nuinHi_part*(uiph[0] - uiph[1]);

                nuinHei_part = nust[zk][yj][xi][45]+nust[zk][yj][xi][46]+nust[zk][yj][xi][47]+nust[zk][yj][xi][48];

                gg[k][j][i].fx[13] += nuinHei_part*(uir[0] - uir[2]) - nuinHiHei_r - nuinOiHei_r;
                gg[k][j][i].fx[14] += nuinHei_part*(uith[0] - uith[2]) - nuinHiHei_th - nuinOiHei_th;
                gg[k][j][i].fx[15] += nuinHei_part*(uiph[0] - uiph[2]) - nuinHiHei_ph - nuinOiHei_ph;

                //----------- non-stiff terms of Te equation
                uer = uu[k][j][i].fx[7]; ueth = uu[k][j][i].fx[8]; ueph = uu[k][j][i].fx[9];
                div_ue  = divergence(uu, i, j, k, yj, xi, 7);

                gg[k][j][i].fx[19]=-(uer*gradTe.r + ueth*gradTe.t + ueph*gradTe.p)
                                   -two3rd*(div_ue + uu[k][j][i].fx[18]/ne);

                //----------- non-stiff terms of unr equation
                Nn=uu[k][j][i].fx[10]; rhon=uu[k][j][i].fx[11]; Tn=xx[k][j][i].fx[30];
                unr=xx[k][j][i].fx[27]; unth=xx[k][j][i].fx[28]; unph=xx[k][j][i].fx[29];

                gradTn=gradient(xx, i, j, k, yj, xi, 30);
                gradNn=gradient(uu, i, j, k, yj, xi, 10);
                grad_unr=gradient(xx, i, j, k, yj, xi, 27);
                grad_untheta=gradient(xx, i, j, k, yj, xi, 28);
                grad_unphi=gradient(xx, i, j, k, yj, xi, 29);

                gg[k][j][i].fx[27]= (unth*unth+unph*unph)/rr[i] - gr[xi]
                                   -(unr*grad_unr.r + unth*grad_unr.t + unph*grad_unr.p)
                                   -(Nn*gradTn.r + Tn*gradNn.r)/rhon
                                   -2.0*(rotat_t[zk][yj]*unph-rotat_p[zk][yj]*unth)-cenf_r[zk][yj][xi];

                //----------- non-stiff terms of untheta equation
                gg[k][j][i].fx[28]= unph*unph*cot_div_r[yj][xi]-unth*unr/rr[i]
                                   -(unr*grad_untheta.r + unth*grad_untheta.t + unph*grad_untheta.p)
                                   -(Nn*gradTn.t + Tn*gradNn.t)/rhon
                                   -2.0*(rotat_p[zk][yj]*unr-rotat_r[zk][yj]*unph)-cenf_t[zk][yj][xi];

                //----------- non-stiff terms of unphi equation
                gg[k][j][i].fx[29]=-unth*unph*cot_div_r[yj][xi] - unph*unr/rr[i]
                                   -(unr*grad_unphi.r + unth*grad_unphi.t + unph*grad_unphi.p)
                                   -(Nn*gradTn.p + Tn*gradNn.p)/rhon
                                   -2.0*(rotat_r[zk][yj]*unth-rotat_t[zk][yj]*unr)-cenf_p[zk][yj][xi];

                //----------- non-stiff terms of Tn equation
                div_un=divergence(uu, i, j, k, yj, xi, 27);

                gg[k][j][i].fx[30]= two3rd*((uu[k][j][i].fx[12]-uu[k][j][i].fx[13]-uu[k][j][i].fx[19])/Nn-Tn*div_un)
                                   -(unr*gradTn.r + unth*gradTn.t + unph*gradTn.p);

                //----------- non-stiff magnetic field equation
                gg[k][j][i].fx[31] = 0.0; gg[k][j][i].fx[32] = 0.0; gg[k][j][i].fx[33] = 0.0;

                //----------- non-stiff electric field equation
                efd_gradPe = E_gradPe(xx, uu, i, j, k, xi, yj, zk, xm, ym, zm);
                ene=e*ne;

                gg[k][j][i].fx[34]=-efd_gradPe.r/ene;
                gg[k][j][i].fx[35]=efd_gradPe.t/ene;
                gg[k][j][i].fx[36]=efd_gradPe.p/ene;
            }

            if (xs == 0) lower_boundary_bc(xx, gg, j, k);
            if (xs+xm == a1) upper_boundary_bc(xx, gg, j, k, yj, zk);
        }
    }

    DMDAVecRestoreArray(da, localX, &xx);
    DMDAVecRestoreArray(da, localU, &uu);

    DMRestoreLocalVector(da, &localX);
    DMRestoreLocalVector(da, &localU);

    DMDAVecRestoreArray(da, G, &gg);

    return 0;
}
