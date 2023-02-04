#include "funcdef.h"
#include "boundary.h"

#include <fstream>
#include <iostream>

using namespace std;

#include "electric_field.h"
#include "ele_cooling_rate.h"
#include "operators.h"
#include "neu_cooling_rate.h"

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
    double   Nn, rhon, Tn;    
    vector3D gradns, gradNq, gradTe, gradne, gradTs, gradNn, gradTn;
    vector3D grad_usr, grad_ustheta, grad_usphi, grad_unr, grad_untheta, grad_unphi;
    double   uir[3], uith[3], uiph[3], uirr, uitht, uiphp, unr, unth, unph, uer, ueth, ueph;
    double   div_ui[3], div_un, div_ue, divui, ne, Te, Ts, Ec_gradPe[3];
    double   jr, jtheta, jphi, B0r, B0t, B0p, ms_ne, Qeuv, Qphoto, Cn, Ce;

    const double two3rd=2.0/3.0;

    if (xs+xm == a1) top_bc_vel(params, ys, ym, zs, zm);

    for (k = zs; k < zs+zm; k++) {
        zk=k-zs;

        for (j = ys; j< ys+ym; j++) {
            //boundary conditions at j=0 and j=Nth are set later
            if (j == 0 || j == Nth) continue;

            yj=j-ys;

            for (i = xs; i < xs+xm; i++) {
                if (i == 0 || i == Nr) continue;

                xi=i-xs;

                uir[0]=xx[k][j][i].fx[7];  uith[0]=xx[k][j][i].fx[8];  uiph[0]=xx[k][j][i].fx[9];
                uir[1]=xx[k][j][i].fx[10]; uith[1]=xx[k][j][i].fx[11]; uiph[1]=xx[k][j][i].fx[12];
                uir[2]=xx[k][j][i].fx[13]; uith[2]=xx[k][j][i].fx[14]; uiph[2]=xx[k][j][i].fx[15];

                div_ui[0]=divergence(xx, i, j, k, yj, xi, 7);
                div_ui[1]=divergence(xx, i, j, k, yj, xi, 10);
                div_ui[2]=divergence(xx, i, j, k, yj, xi, 13);
                div_un   =divergence(xx, i, j, k, yj, xi, 27);

                gradTe=gradient(xx, i, j, k, yj, xi, 19);
                gradne=gradient(uu, i, j, k, yj, xi, 6);

                ne=uu[k][j][i].fx[6]; Te=xx[k][j][i].fx[19];

                jr=uu[k][j][i].fx[16]; jtheta=uu[k][j][i].fx[17]; jphi=uu[k][j][i].fx[18];
                B0r=uu[k][j][i].fx[0]; B0t=uu[k][j][i].fx[1]; B0p=uu[k][j][i].fx[2];
                unr=xx[k][j][i].fx[27]; unth=xx[k][j][i].fx[28]; unph=xx[k][j][i].fx[29];

                //production and loss rates
                prod_loss_rates(xx, uu, i, j, k, zk, yj, xi, Qeuv, Qphoto);

                rhon=0.0;
                for (s = 0; s < sl; s++) {
                    rhon += exp(xx[k][j][i].fx[s+20])*ams[s];

                    gradns = gradient(xx, i, j, k, yj, xi, s);
                    gradNq = gradient(xx, i, j, k, yj, xi, s+20);

                    if (s < 3) {
                        divui=div_ui[s]; uirr=uir[s]; uitht=uith[s]; uiphp=uiph[s];
                    }
                    else {
                        divui=div_ui[0]; uirr=uir[0]; uitht=uith[0]; uiphp=uiph[0];
                    }

                    gg[k][j][i].fx[s]= Ps[zk][yj][xi][s]-(Ls[zk][yj][xi][s]+divui)*xx[k][j][i].fx[s]
                                      -uirr*gradns.r-uitht*gradns.t-uiphp*gradns.p;
                    gg[k][j][i].fx[s+20]= Ps[zk][yj][xi][s+7]*exp(-xx[k][j][i].fx[s+20])-Ls[zk][yj][xi][s+7]
                                         -div_un - unr*gradNq.r-unth*gradNq.t-unph*gradNq.p;

                    if (s >= 3) continue;

                    Ts=xx[k][j][i].fx[s+16]; ms_ne=ams[s]*ne;

                    //----------- non-stiff terms of uir equation
                    s7=7+s*3; s8=8+s*3; s9=9+s*3; s16=16+s;
                    grad_usr=gradient(xx, i, j, k, yj, xi, s7);
                    grad_ustheta=gradient(xx, i, j, k, yj, xi, s8);
                    grad_usphi=gradient(xx, i, j, k, yj, xi, s9);

                    gradTs = gradient(xx, i, j, k, yj, xi, s16);

                    gg[k][j][i].fx[s7]= (jtheta*B0p-jphi*B0t)/ms_ne-gr[xi]+(uith[s]*uith[s]+uiph[s]*uiph[s])/rr[i]
                                       -(uir[s]*grad_usr.r + uith[s]*grad_usr.t + uiph[s]*grad_usr.p)
                                       -(gradTs.r + Ts*gradns.r + gradTe.r + Te/ne*gradne.r)/ams[s];

                    //----------- non-stiff terms of uitheta equation
                    gg[k][j][i].fx[s8]= (jphi*B0r-jr*B0p)/ms_ne-uith[s]*uir[s]/rr[i]+uiph[s]*uiph[s]*cot_div_r[yj][xi]
                                       -(uir[s]*grad_ustheta.r + uith[s]*grad_ustheta.t + uiph[s]*grad_ustheta.p)
                                       -(gradTs.t + Ts*gradns.t + gradTe.t + Te/ne*gradne.t)/ams[s];

                    //----------- non-stiff terms of uiphi equation
                    gg[k][j][i].fx[s9]= (jr*B0t-jtheta*B0r)/ms_ne-uiph[s]*uir[s]/rr[i]-uith[s]*uiph[s]*cot_div_r[yj][xi]
                                       -(uir[s]*grad_usphi.r + uith[s]*grad_usphi.t + uiph[s]*grad_usphi.p)
                                       -(gradTs.p + Ts*gradns.p + gradTe.p + Te/ne*gradne.p)/ams[s];

                    //----------- non-stiff terms of Ti equation
                    gg[k][j][i].fx[s16]= two3rd*(uu[k][j][i].fx[s+11]/xx[k][j][i].fx[s] - Ts*div_ui[s])
                                        -(uir[s]*gradTs.r + uith[s]*gradTs.t + uiph[s]*gradTs.p);
                }

                //----------- non-stiff terms of Te equation
                uer = uu[k][j][i].fx[7]; ueth = uu[k][j][i].fx[8]; ueph = uu[k][j][i].fx[9];
                div_ue  = divergence(uu, i, j, k, yj, xi, 7);

                Tn=xx[k][j][i].fx[30];
                Ce=ele_cooling_rate(xx, Te, Tn, ne, i, j, k);

                gg[k][j][i].fx[19]= two3rd*((Qphoto + uu[k][j][i].fx[14]-Ce)/ne - Te*div_ue)
                                   -(uer*gradTe.r + ueth*gradTe.t + ueph*gradTe.p);

                //----------- non-stiff terms of unr equation
                Nn=uu[k][j][i].fx[10]; Tn=xx[k][j][i].fx[30];
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
                Cn=neu_cooling_rate(xx, i, j, k);

                gg[k][j][i].fx[30]= two3rd*((Qeuv-Cn + uu[k][j][i].fx[15])/Nn-Tn*div_un)
                                   -(unr*gradTn.r + unth*gradTn.t + unph*gradTn.p);

                //----------- non-stiff magnetic field equation
                gg[k][j][i].fx[31] = 0.0; gg[k][j][i].fx[32] = 0.0; gg[k][j][i].fx[33] = 0.0;

                //----------- non-stiff electric field equation
                E_gradPe(xx, uu, i, j, k, xi, yj, zk, Ec_gradPe);

                gg[k][j][i].fx[34]=-( Jiv11[zk][yj]*Ec_gradPe[0]+Jiv21[zk][yj]*Ec_gradPe[1]
                                     +Jiv31[yj]*Ec_gradPe[2]);

                gg[k][j][i].fx[35]=-( Jiv12[zk][yj][xi]*Ec_gradPe[0]+Jiv22[zk][yj][xi]*Ec_gradPe[1]
                                     +Jiv32[yj][xi]*Ec_gradPe[2]);

                gg[k][j][i].fx[36]=-(Jiv13[zk][yj][xi]*Ec_gradPe[0]+Jiv23[zk][yj][xi]*Ec_gradPe[1]);

                for (s=0; s<a4; s++) {
                    if (isnan(gg[k][j][i].fx[s]) || isinf(gg[k][j][i].fx[s])) {
                        cout<<"Rightside function is Nan or inf at ("<<i<<", "<<j<<", "<<k<<", "<<s
                            <<") in rhs_functions "<<gg[k][j][i].fx[s]<<endl;
                        exit(-1);
                    }
                }
            }

            if (xs == 0) lower_boundary_bc(xx, gg, j, k);
            if (xs+xm == a1) upper_boundary_bc(xx, gg, j, k, yj, zk);
        }
    }

        //set boundary conditions at j=0 and j=Nth
    if (ys == 0) {
        for (k = zs; k < zs+zm; k++) {
            for (i = xs; i < xs+xm; i++) {
                for (s = 0; s < a4; s++) {
                    if (s == 8 || s == 11 || s == 14 || s == 28 || s ==32 || s == 35)
                        gg[k][0][i].fx[s]=-xx[(k+a3/2) % a3][1][i].fx[s];
                    else gg[k][0][i].fx[s]=xx[(k+a3/2) % a3][1][i].fx[s];
                }
            }
        }
    }
    if (ys+ym == a2) {
        for (k = zs; k < zs+zm; k++) {
            for (i = xs; i < xs+xm; i++) {
                for (s = 0; s < a4; s++) {
                    if (s == 8 || s == 11 || s == 14 || s == 28 || s ==32 || s == 35)
                        gg[k][Nth][i].fx[s]=-xx[(k+a3/2) % a3][Nthm][i].fx[s];
                    else gg[k][Nth][i].fx[s]=xx[(k+a3/2) % a3][Nthm][i].fx[s];
                }
            }
        }
    }

    DMDAVecRestoreArray(da, localX, &xx);
    DMDAVecRestoreArray(da, localU, &uu);

    DMRestoreLocalVector(da, &localX);
    DMRestoreLocalVector(da, &localU);

    DMDAVecRestoreArray(da, G, &gg);

    return 0;
}
