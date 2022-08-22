#include "fluxes.h"
#include "funcdef.h"
#include "check_positivity.h"

#include <fstream>
#include <iostream>

using namespace std;

/** evaluate source terms excluding collision terms */
void source_terms(Field ***,Field ***,Field ***,Field ***,int,int,int,int,int,int,int,double [], double[]);

int forward_scheme(DM da, Vec Xn, Vec Xn1, AppCtx *params)
{
    double source[nvar-3], Ls[14];
    Vec    localX, localU, localV, localW, localZ;
    Field  ***xn, ***xn1, ***uu, ***vv, ***ww, ***zz;
    int    xs, xm, ys, ym, zs, zm;

    DMDAGetCorners(da, &xs ,&ys, &zs, &xm, &ym, &zm);

    DMGetLocalVector(da, &localX);
    DMGlobalToLocalBegin(da, Xn1, INSERT_VALUES,localX);
    DMGlobalToLocalEnd(da, Xn1, INSERT_VALUES, localX);
    DMDAVecGetArray(da, localX, &xn1);

    /*----------------------------------------------------------------------
     * Calculate various parameters and store in vectors U, V, W, Z.
     * xn1 used is local array.
     * This routine must be called before declaring local uu, vv, ww, & zz
     *---------------------------------------------------------------------*/
    parameters(da, xn1, params);
    /*---------------------------------------------------------------------
     *---------------------------------------------------------------------*/

    //now get local arrays uu, vv, ww, zz
    DMGetLocalVector(da, &localU);
    DMGlobalToLocalBegin(da,params->U,INSERT_VALUES,localU);
    DMGlobalToLocalEnd(da,params->U,INSERT_VALUES,localU);
    DMDAVecGetArray(da, localU, &uu);

    DMGetLocalVector(da, &localV);
    DMGlobalToLocalBegin(da,params->V,INSERT_VALUES,localV);
    DMGlobalToLocalEnd(da,params->V,INSERT_VALUES,localV);
    DMDAVecGetArray(da, localV, &vv);

    DMGetLocalVector(da, &localW);
    DMGlobalToLocalBegin(da,params->W,INSERT_VALUES,localW);
    DMGlobalToLocalEnd(da,params->W,INSERT_VALUES,localW);
    DMDAVecGetArray(da, localW, &ww);

    DMGetLocalVector(da, &localZ);
    DMGlobalToLocalBegin(da,params->Z,INSERT_VALUES,localZ);
    DMGlobalToLocalEnd(da,params->Z,INSERT_VALUES,localZ);
    DMDAVecGetArray(da, localZ, &zz);

    DMDAVecGetArray(da, Xn, &xn);

    int    i, j, k, xi, yj, zk, s, nvarm3=nvar-3, kc;
    int    ip, jm, jp, kp, kprime;
    double flux_rthph[nvar-3], Fr_atfaces[2][nvar-3], Ftheta_Lface[nvar-3], Fphi_Lface[nvar-3];
    double Erss;

    for (k = zs; k < zs+zm; k++) {
        zk=k-zs;
        if (k < Np) kp=k+1; else kp = 0;

        for (j = ys; j< ys+ym; j++) {
            yj=j-ys; jm=j-1; jp=j+1;

            if (j == Nth) {
                for (i = xs; i < xs+xm; i++) {
                    if (i == 0 || i == Nr) continue;
                    else {
                        kc=(k+a3/2) % a3;
                        if (kc + 1 > Np) kprime=0; else kprime=kc+1;

                        xn[k][j][i].fx[24]= xn1[k][j][i].fx[24]
                                           +dt_half*( ( vv[kprime][jm][i].fx[23]-vv[k][jm][kp].fx[23]
                                                     -vv[kc][jm][i].fx[23]+vv[k][jm][k].fx[23])/(2.0*rfavg_dth[i]*dph)
                                                     +(rh[i+1]*vv[k][j][i+1].fx[25]-rh[i]*vv[k][j][i].fx[25])/rh_d2[i]);

                        if (isnan(xn[k][j][i].fx[24]) || isinf(xn[k][j][i].fx[24])) {
                            cout<<"Solution is Nan or inf at ("<<i<<", "<<j<<", "<<k
                                <<", "<<24<<") in forwars_scheme"<<endl;
                            exit(-1);
                        }
                    }
                }

                continue;
            }

            for (i = xs; i < xs+xm; i++) {
                xi=i-xs; ip=i+1;

                if (i == 0 || i == Nr) continue;

                source_terms(xn1, uu, ww, zz, xs, i, j, k, zk, yj, xi, source, Ls);

                /*-------numerical fluxes at the bottom face i-1/2 of the cell i
                 *-------------------------------------------------------------------*/
                if (i==1 || (i==xs && i > 0)) {
                    fluxes_r(xn1, uu, vv, zz, i, j, k, flux_rthph);
                    for (s=0; s< nvarm3; s++) Fr_atfaces[0][s]=flux_rthph[s];
                }
                else {
                    for (s=0; s< nvarm3; s++) Fr_atfaces[0][s]=Fr_atfaces[1][s];
                }

                //numerical fluxes at the top face i+1/2 of the cell (i, j, k)
                fluxes_r(xn1, uu, vv, zz, ip, j, k, flux_rthph);
                for (s=0; s< nvarm3; s++) Fr_atfaces[1][s]=flux_rthph[s];

                /*--------numerical fluxes at the left face j-1/2 of the cell (i, j, k)
                 *-------------------------------------------------------------------*/
                if (j == ys) {
                    if (j == 0) for (s=0; s<nvarm3; s++) Ftheta_Lface[s] = 0.0;
                    else {
                        fluxes_theta(xn1, uu, ww, zz, i, j, k, flux_rthph);
                        for (s=0; s<nvarm3; s++) Ftheta_Lface[s]=flux_rthph[s];
                    }
                }
                else {
                    for (s=0; s<nvarm3; s++) Ftheta_Lface[s]=Ftheta_Rface[xi][s];
                }

                //numerical fluxes at the right face j+1/2 of the cell (i, j, k)
                if (j < Nthm) {
                    fluxes_theta(xn1, uu, ww, zz, i, jp, k, flux_rthph);
                    for (s=0; s<nvarm3; s++) Ftheta_Rface[xi][s]=flux_rthph[s];
                }
                else for (s=0; s<nvarm3; s++) Ftheta_Rface[xi][s]=0.0;

                /*---------numerical fluxes at the back face k-1/2 of the cell (i, j, k)
                 *---------------------------------------------------------------------*/
                if (k == zs) {
                    fluxes_phi(xn1, uu, zz, i, j, k, flux_rthph);
                    for (s=0; s<nvarm3; s++) Fphi_Lface[s]=flux_rthph[s];
                }
                else {
                    for (s=0; s<nvarm3; s++) Fphi_Lface[s]=Fphi_Rface[yj][xi][s];
                }

                //numerical fluxes at the front face k+1/2 of the cell (i, j, k)
                fluxes_phi(xn1, uu, zz, i, j, kp, flux_rthph);
                for (s=0; s<nvarm3; s++) Fphi_Rface[yj][xi][s]=flux_rthph[s];

                //solution of fluid equations
                for (s=0; s<nvarm3; s++) {
                    xn[k][j][i].fx[s]= xn1[k][j][i].fx[s]
                                      -dt_half*( (rh2[ip]*Fr_atfaces[1][s]-rh2[i]*Fr_atfaces[0][s])/rh_d3[i]
                                                +(sinth_h[jp]*Ftheta_Rface[xi][s]-sinth_h[j]*Ftheta_Lface[s])
                                                 /rfavg_costh[j][i]
                                                +(Fphi_Rface[yj][xi][s]-Fphi_Lface[s])/rfavg_costh_dth_dph[j][i]
                                                -source[s]);
                }

                //solution of magnetic equation
                xn[k][j][i].fx[23]= xn1[k][j][i].fx[23]
                                   +dt_half*( (vv[kp][j][i].fx[24]-vv[k][j][i].fx[24])/rh_costh_dth_dph[j][i]
                                             -(sinth_h[jp]*vv[k][jp][i].fx[25]-sinth_h[j]*vv[k][j][i].fx[25])
                                              /rh_costh[j][i]);

                if (j == 0) {
                    kc = (k+a3/2) % a3;
                    if (kc + 1 > Np) kprime=0; else kprime=kc+1;
                    xn[k][j][i].fx[24]= xn1[k][j][i].fx[24]
                                       +dt_half*( (rh[ip]*vv[k][j][ip].fx[25]-rh[i]*vv[k][j][i].fx[25])/rh_d2[i]
                                                 -( vv[kp][j][i].fx[23]-vv[kprime][j][i].fx[23]-vv[k][j][i].fx[23]
                                                   +vv[kc][j][i].fx[23])/(2.0*rfavg_dth[i]*dph));

                    Erss = 0.0;
                    for (s = 0; s < a3; s++) Erss += vv[s][j][i].fx[23];
                    Erss = Erss/(double)a3;

                    xn[k][j][i].fx[25]= xn1[k][j][i].fx[25]
                                       +dt_half*( (vv[k][jp][i].fx[23]-Erss)/rfavg_dth[i]
                                                 -(rh[ip]*vv[k][j][ip].fx[24]-rh[i]*vv[k][j][i].fx[24])/rh_d2[i]);
                }
                else {
                    xn[k][j][i].fx[24]= xn1[k][j][i].fx[24]
                                       +dt_half*( (rh[ip]*vv[k][j][ip].fx[25]-rh[i]*vv[k][j][i].fx[25])/rh_d2[i]
                                                 -(vv[kp][j][i].fx[23]-vv[k][j][i].fx[23])/rfavg_sinth_dph[j][i]);

                    if (j < Nthm)
                        xn[k][j][i].fx[25]= xn1[k][j][i].fx[25]
                                           +dt_half*( (vv[k][jp][i].fx[23]-vv[k][j][i].fx[23])/rfavg_dth[i]
                                                     -(rh[ip]*vv[k][j][ip].fx[24]-rh[i]*vv[k][j][i].fx[24])/rh_d2[i]);
                    else {
                        Erss = 0.0;
                        for (s = 0; s < a3; s++) Erss += vv[s][jm][i].fx[23];
                        Erss = Erss/(double)a3;

                        xn[k][j][i].fx[25]= xn1[k][j][i].fx[25]
                                           +dt_half*( (Erss - vv[k][jm][i].fx[23])/rfavg_dth[i]
                                                     -(rh[ip]*vv[k][j][ip].fx[24]-rh[i]*vv[k][j][i].fx[24])/rh_d2[i]);
                    }
                }

                for (s=0; s<nvar; s++) {
                    if (isnan(xn[k][j][i].fx[s]) || isinf(xn[k][j][i].fx[s])) {
                        cout<<"Solution is Nan or inf at ("<<i<<", "<<j<<", "<<k
                            <<", "<<s<<") in forward_scheme"<<endl;
                        exit(-1);
                    }
                }

                if(check_positivity(xn, i, j, k, 0) < 0) exit(-12);
            }
        }
    }

    boundary_bc(xn, xs, xm, ys, ym, zs, zm, params);
    DMDAVecRestoreArray(da, Xn, &xn);

    DMDAVecRestoreArray(da, localX, &xn1);
    DMDAVecRestoreArray(da, localU, &uu);
    DMDAVecRestoreArray(da, localV, &vv);
    DMDAVecRestoreArray(da, localW, &ww);
    DMDAVecRestoreArray(da, localZ, &zz);

    DMRestoreLocalVector(da, &localX);
    DMRestoreLocalVector(da, &localU);
    DMRestoreLocalVector(da, &localV);
    DMRestoreLocalVector(da, &localW);
    DMRestoreLocalVector(da, &localZ);

    params->sec += dt_half;

    return 0;
}