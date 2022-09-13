#include "fluxes.h"
#include "funcdef.h"
#include "check_positivity.h"
#include "boundary.h"

#include <fstream>
#include <iostream>

using namespace std;

#include "ele_cooling_rate.h"

/** evaluate source terms excluding collision terms */
void source_terms(Field ***,Field ***,Field ***,Field ***,int,int,int,int,int,int,int,double []);

int rhsfunctions(TS ts, double ftime, Vec X, Vec G, void* ctx)
{
    double source[nvar-3];
    Vec    localX, localU, localV, localW, localZ;
    Field  ***xx, ***uu, ***vv, ***ww, ***zz, ***gg;
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

    DMDAVecGetArray(da, G, &gg);

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
                    for (s=0; s<nvar; s++) {
                        if (s != 24) gg[k][j][i].fx[s] = xx[(k+a3/2) % a3][Nthm][i].fx[s];
                    }

                    kc=(k+a3/2) % a3;
                    if (kc + 1 > Np) kprime=0; else kprime=kc+1;

                    gg[k][j][i].fx[24]= ( vv[kprime][jm][i].fx[23]-vv[k][jm][kp].fx[23]
                                         -vv[kc][jm][i].fx[23]+vv[k][jm][k].fx[23])/(2.0*rfavg_dth[i]*dph)
                                       +(rh[i+1]*vv[k][j][i+1].fx[25]-rh[i]*vv[k][j][i].fx[25])/rh_d2[i];

                    if (isnan(xx[k][j][i].fx[24]) || isinf(xx[k][j][i].fx[24])) {
                        cout<<"Solution is Nan or inf at ("<<i<<", "<<j<<", "<<k
                            <<", "<<24<<") in leap_frog advace"<<endl;
                        exit(-1);
                    }
                }

                continue;
            }

            for (i = xs; i < xs+xm; i++) {
                xi=i-xs; ip=i+1;

                if (i == 0) {
                    lower_boundary_bc(xx, gg, j, k);
                    continue;
                }
                else if (i == Nr) {
                    upper_boundary_bc(xx, gg, j, k, yj, zk);
                    continue;
                }

                source_terms(xx, uu, ww, zz, xs, i, j, k, zk, yj, xi, source);

                /*-------numerical fluxes at the bottom face i-1/2 of the cell i
                 *-------------------------------------------------------------------*/
                if (i==1 || (i==xs && i > 0)) {
                    fluxes_r(xx, uu, vv, i, j, k, flux_rthph);
                    for (s=0; s< nvarm3; s++) Fr_atfaces[0][s]=flux_rthph[s];
                }
                else {
                    for (s=0; s< nvarm3; s++) Fr_atfaces[0][s]=Fr_atfaces[1][s];
                }

                //numerical fluxes at the top face i+1/2 of the cell (i, j, k)
                fluxes_r(xx, uu, vv, ip, j, k, flux_rthph);
                for (s=0; s< nvarm3; s++) Fr_atfaces[1][s]=flux_rthph[s];

                /*--------numerical fluxes at the left face j-1/2 of the cell (i, j, k)
                 *-------------------------------------------------------------------*/
                if (j == ys) {
                    if (j == 0) for (s=0; s<nvarm3; s++) Ftheta_Lface[s] = 0.0;
                    else {
                        fluxes_theta(xx, uu, ww, i, j, k, flux_rthph);
                        for (s=0; s<nvarm3; s++) Ftheta_Lface[s]=flux_rthph[s];
                    }
                }
                else {
                    for (s=0; s<nvarm3; s++) Ftheta_Lface[s]=Ftheta_Rface[xi][s];
                }

                //numerical fluxes at the right face j+1/2 of the cell (i, j, k)
                if (j < Nthm) {
                    fluxes_theta(xx, uu, ww, i, jp, k, flux_rthph);
                    for (s=0; s<nvarm3; s++) Ftheta_Rface[xi][s]=flux_rthph[s];
                }
                else for (s=0; s<nvarm3; s++) Ftheta_Rface[xi][s]=0.0;

                /*---------numerical fluxes at the back face k-1/2 of the cell (i, j, k)
                 *---------------------------------------------------------------------*/
                if (k == zs) {
                    fluxes_phi(xx, uu, zz, i, j, k, flux_rthph);
                    for (s=0; s<nvarm3; s++) Fphi_Lface[s]=flux_rthph[s];
                }
                else {
                    for (s=0; s<nvarm3; s++) Fphi_Lface[s]=Fphi_Rface[yj][xi][s];
                }

                //numerical fluxes at the front face k+1/2 of the cell (i, j, k)
                fluxes_phi(xx, uu, zz, i, j, kp, flux_rthph);
                for (s=0; s<nvarm3; s++) Fphi_Rface[yj][xi][s]=flux_rthph[s];

                //RHS of fluid equations
                for (s=0; s<nvarm3; s++) {
                    gg[k][j][i].fx[s]= source[s]-(rh2[ip]*Fr_atfaces[1][s]-rh2[i]*Fr_atfaces[0][s])/rh_d3[i]
                                      -(sinth_h[jp]*Ftheta_Rface[xi][s]-sinth_h[j]*Ftheta_Lface[s])/rfavg_costh[j][i]
                                      -(Fphi_Rface[yj][xi][s]-Fphi_Lface[s])/rfavg_costh_dth_dph[j][i];
                }

                //right hand side of magnetic equation
                gg[k][j][i].fx[23]= (vv[kp][j][i].fx[24]-vv[k][j][i].fx[24])/rh_costh_dth_dph[j][i]
                                   -(sinth_h[jp]*vv[k][jp][i].fx[25]-sinth_h[j]*vv[k][j][i].fx[25])/rh_costh[j][i];

                if (j == 0) {
                    kc = (k+a3/2) % a3;
                    if (kc + 1 > Np) kprime=0; else kprime=kc+1;
                    gg[k][j][i].fx[24]= (rh[ip]*vv[k][j][ip].fx[25]-rh[i]*vv[k][j][i].fx[25])/rh_d2[i]
                                       -( vv[kp][j][i].fx[23]-vv[kprime][j][i].fx[23]-vv[k][j][i].fx[23]
                                         +vv[kc][j][i].fx[23])/(2.0*rfavg_dth[i]*dph);

                    Erss = 0.0;
                    for (s = 0; s < a3; s++) Erss += vv[s][j][i].fx[23];
                    Erss = Erss/(double)a3;

                    gg[k][j][i].fx[25]= (vv[k][jp][i].fx[23]-Erss)/rfavg_dth[i]
                                       -(rh[ip]*vv[k][j][ip].fx[24]-rh[i]*vv[k][j][i].fx[24])/rh_d2[i];
                }
                else {
                    gg[k][j][i].fx[24]= (rh[ip]*vv[k][j][ip].fx[25]-rh[i]*vv[k][j][i].fx[25])/rh_d2[i]
                                       -(vv[kp][j][i].fx[23]-vv[k][j][i].fx[23])/rfavg_sinth_dph[j][i];

                    if (j < Nthm)
                        gg[k][j][i].fx[25]= (vv[k][jp][i].fx[23]-vv[k][j][i].fx[23])/rfavg_dth[i]
                                           -(rh[ip]*vv[k][j][ip].fx[24]-rh[i]*vv[k][j][i].fx[24])/rh_d2[i];
                    else {
                        Erss = 0.0;
                        for (s = 0; s < a3; s++) Erss += vv[s][jm][i].fx[23];
                        Erss = Erss/(double)a3;

                        gg[k][j][i].fx[25]= (Erss - vv[k][jm][i].fx[23])/rfavg_dth[i]
                                           -(rh[ip]*vv[k][j][ip].fx[24]-rh[i]*vv[k][j][i].fx[24])/rh_d2[i];
                    }
                }

                for (s=0; s<nvar; s++) {
                    if (isnan(gg[k][j][i].fx[s]) || isinf(gg[k][j][i].fx[s])) {
                        cout<<"function is Nan or inf at ("<<i<<", "<<j<<", "<<k
                            <<", "<<s<<") in rhsfunctions"<<endl;
                        exit(-1);
                    }
                }
            }
        }
    }

    DMDAVecRestoreArray(da, localX, &xx);
    DMDAVecRestoreArray(da, localU, &uu);
    DMDAVecRestoreArray(da, localV, &vv);
    DMDAVecRestoreArray(da, localW, &ww);
    DMDAVecRestoreArray(da, localZ, &zz);

    DMRestoreLocalVector(da, &localX);
    DMRestoreLocalVector(da, &localU);
    DMRestoreLocalVector(da, &localV);
    DMRestoreLocalVector(da, &localW);
    DMRestoreLocalVector(da, &localZ);

    DMDAVecRestoreArray(da, G, &gg);

    return 0;
}

int stifffunction(TS ts, double ftime, Vec X, Vec Xdt, Vec F, void* ctx)
{
    int  s, t;
    double rhos, rhoi, rhon, sum_nues, sum_nueq, sum_rhonusq, sum_nues_div_ms, sum_nueq_div_mq;
    double rhos_nusq_msmq, rhos_nusqmq_msmq, rhos_nusqms_msmq, temp;
    double ne, rhoe, Nn; //ne = ni
    double rhoe_sum_nues, rhoe_sum_nueq, uiminusun_sq;
    double uer, uet, uep, uir, uit, uip, unr, unt, unp;
    double Qefric, Qifric, Qnfric, Te, Tn;
    const double two3rd=2.0/3.0, one6th=1.0/6.0;

    Vec    localX, localU, localZ;
    AppCtx *params = (AppCtx*)ctx;
    Field  ***xx, ***dxdt, ***uu, ***zz, ***ff;
    DM     da;
    PetscInt xs, xm, ys, ym, zs, zm, zk, yj, xi, i, j, k, im, ip, jm, kcm, jp, kcp, km, kp;
    double fsf_imjk, fsf_ipjk;
    double a_imjk[2], a_imjk_multi_a_imjk_subtr, a_ipjk[2], a_ipjk_multi_a_ipjk_subtr;
    double fsf_ijmk, fsf_ijpk;
    double b_ijmk[2], b_ijmk_multi_b_ijmk_subtr, b_ijpk[2], b_ijpk_multi_b_ijpk_subtr;
    double fsf_ijkm, fsf_ijkp;
    double c_ijkm[2], c_ijkm_multi_c_ijkm_subtr, c_ijkp[2], c_ijkp_multi_c_ijkp_subtr;

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

    DMGetLocalVector(da, &localZ);
    DMGlobalToLocalBegin(da, params->Z, INSERT_VALUES,localZ);
    DMGlobalToLocalEnd(da, params->Z, INSERT_VALUES, localZ);
    DMDAVecGetArray(da, localZ, &zz);

    DMDAVecGetArray(da, Xdt, &dxdt);
    DMDAVecGetArray(da, F, &ff);

    DMDAGetCorners(da, &xs ,&ys, &zs, &xm, &ym, &zm);

    double dr2=dr*dr;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (k = zs; k < zs+zm; k++) {
        zk=k-zs;

        km = k-1; kp = k+1;

        for (j = ys; j < ys+ym; j++) {
            yj=j-ys;

            if (j == Nth) {
                for (i = xs; i < xs+xm; i++) {
                    for (s = 0; s< nvar; s++) {
                        if (s != 24) ff[k][j][i].fx[s] = xx[k][j][i].fx[s];
                    }

                    ff[k][j][i].fx[24]=dxdt[k][j][i].fx[24];
                }

                continue;
            }

            if (j == 0) {jm = j; kcm = (k+a3/2) % a3;}
            else {jm = j-1; kcm = k;}
            if (j < Nthm) {jp = j+1; kcp = k;}
            else {jp = Nthm; kcp = (k+a3/2) % a3;}

            for (i = xs; i < xs+xm; i++) {
                if (i == 0 or i == Nr) {
                    for (s = 0; s< nvar; s++) ff[k][j][i].fx[s] = xx[k][j][i].fx[s];
                    continue;
                }

                xi=i-xs;
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
                b_ijpk[0]=uu[kcp][jp][i].fx[5];
                b_ijpk[1]=uu[kcp][jp][i].fx[6];
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
                    fsf_imjk = a_imjk_multi_a_imjk_subtr*(xx[k][j][i].fx[s]-xx[k][j][im].fx[s]);
                    fsf_ipjk = a_ipjk_multi_a_ipjk_subtr*(xx[k][j][ip].fx[s]-xx[k][j][i].fx[s]);

                    fsf_ijmk = b_ijmk_multi_b_ijmk_subtr*(xx[k][j][i].fx[s]-xx[kcm][jm][i].fx[s]);
                    fsf_ijpk = b_ijpk_multi_b_ijpk_subtr*(xx[kcp][jp][i].fx[s]-xx[k][j][i].fx[s]);

                    fsf_ijkm = c_ijkm_multi_c_ijkm_subtr*(xx[k][j][i].fx[s]-xx[km][j][i].fx[s]);
                    fsf_ijkp = c_ijkp_multi_c_ijkp_subtr*(xx[kp][j][i].fx[s]-xx[k][j][i].fx[s]);

                    ff[k][j][i].fx[s] = dxdt[k][j][i].fx[s] + Ls[zk][yj][xi][s]*xx[k][j][i].fx[s]
                                       +(rh2[ip]*fsf_ipjk - rh2[i]*fsf_imjk)/rh_d3[i]
                                       +(sinth_h[jp]*fsf_ijpk - sinth_h[j]*fsf_ijmk)/rfavg_costh[j][i]
                                       +(fsf_ijkp - fsf_ijkm)/rfavg_costh_dth_dph[j][i];

                    ff[k][j][i].fx[12+s] = dxdt[k][j][i].fx[12+s] + Ls[zk][yj][xi][7+s]*xx[k][j][i].fx[12+s];

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

                uer=xx[k][j][i].fx[7]/rhoi; uet=xx[k][j][i].fx[8]/rhoi; uep=xx[k][j][i].fx[9]/rhoi;
                uir=uer; uit=uet; uip=uep;
                unr=xx[k][j][i].fx[19]/rhon; unt=xx[k][j][i].fx[20]/rhon; unp=xx[k][j][i].fx[21]/rhon;

                ff[k][j][i].fx[7] = dxdt[k][j][i].fx[7]-sum_nues*(uer-uir)-sum_rhonusq*(unr-uir);
                ff[k][j][i].fx[8] = dxdt[k][j][i].fx[8]-sum_nues*(uet-uit)-sum_rhonusq*(unt-uit);
                ff[k][j][i].fx[9] = dxdt[k][j][i].fx[9]-sum_nues*(uep-uip)-sum_rhonusq*(unp-uip);

                rhoe=ne*me;
                rhoe_sum_nues=rhoe*sum_nues;
                rhoe_sum_nueq=rhoe*sum_nueq;

                uiminusun_sq=(uir-unr)*(uir-unr)+(uit-unt)*(uit-unt)+(uip-unp)*(uip-unp);
                Qifric=two3rd*rhos_nusqmq_msmq*uiminusun_sq;
                Qefric=two3rd*( rhoe_sum_nues*((uer-uir)*(uer-uir)+(uet-uit)*(uet-uit)+(uep-uip)*(uep-uip))
                   +rhoe_sum_nueq*((uer-unr)*(uer-unr)+(uet-unt)*(uet-unt)+(uep-unp)*(uep-unp)));
                Qnfric=two3rd*rhos_nusqms_msmq*uiminusun_sq;

                ff[k][j][i].fx[10] = dxdt[k][j][i].fx[10]
                                    +2.0*( me*sum_nues_div_ms*(xx[k][j][i].fx[10]-xx[k][j][i].fx[11])
                                          +rhos_nusq_msmq*(xx[k][j][i].fx[10]/ne-xx[k][j][i].fx[22]/Nn))
                                    +two3rd*Qifric;

                Te=xx[k][j][i].fx[11]/(ne*kb);
                Tn=xx[k][j][i].fx[22]/(Nn*kb);
                ff[k][j][i].fx[11] = dxdt[k][j][i].fx[11]
                                    +2.0*me*( sum_nues_div_ms*(xx[k][j][i].fx[11]-xx[k][j][i].fx[10])
                                             +sum_nueq_div_mq*(xx[k][j][i].fx[11]-ne/Nn*xx[k][j][i].fx[22]))
                                    -two3rd*(Qefric+ele_cooling_rate(xx, Te, Tn, ne, i, j, k))
                                    +one6th/kb*( (zz[k][j][ip].fx[24]-zz[k][j][im].fx[24])/dr2
                                                 *( xx[k][j][ip].fx[11]/uu[k][j][ip].fx[17]
                                                   -xx[k][j][im].fx[11]/uu[k][j][im].fx[17])
                                                +(zz[kcp][jp][i].fx[24]-zz[kcm][jm][i].fx[24])/(rfavg_dth[i]*rfavg_dth[i])
                                                 *( xx[kcp][jp][i].fx[11]/uu[kcp][jp][i].fx[17]
                                                   -xx[kcm][jm][i].fx[11]/uu[kcm][jm][i].fx[17])
                                                +(zz[kp][j][i].fx[24]-zz[km][j][i].fx[24])
                                                 /(rfavg_sinth_dph[j][i]*rfavg_sinth_dph[j][i])
                                                 *( xx[kp][j][i].fx[11]/uu[kp][j][i].fx[17]
                                                   -xx[km][j][i].fx[11]/uu[km][j][i].fx[17]))
                                    +two3rd/kb*zz[k][j][i].fx[24]
                                              *( ( xx[k][j][ip].fx[11]/uu[k][j][ip].fx[17]
                                                  -xx[k][j][im].fx[11]/uu[k][j][im].fx[17])/(rfavg[i]*dr)
                                                +( xx[k][j][ip].fx[11]/uu[k][j][ip].fx[17]
                                                  -2.0*xx[k][j][i].fx[11]/uu[k][j][i].fx[11]
                                                  +xx[k][j][im].fx[11]/uu[k][j][im].fx[17])/dr2
                                                +cotth[j]*( xx[kcp][jp][i].fx[11]/uu[kcp][jp][i].fx[11]
                                                           -xx[kcm][jm][i].fx[11]/uu[kcm][jm][i].fx[11])
                                                          /(rfavg_dth[i]*rfavg[i])
                                                +( xx[kcp][jp][i].fx[11]/uu[kcp][jp][i].fx[17]
                                                  -2.0*xx[k][j][i].fx[11]/uu[k][j][i].fx[11]
                                                  +xx[kcm][jm][i].fx[11]/uu[kcm][jm][i].fx[17])
                                                 /(rfavg_dth[i]*rfavg_dth[i])
                                                +( xx[kp][j][i].fx[11]/uu[kp][j][i].fx[17]
                                                  -2.0*xx[k][j][i].fx[11]/uu[k][j][i].fx[11]
                                                  +xx[km][j][i].fx[11]/uu[km][j][i].fx[17])
                                                 /(rfavg_sinth_dph[j][i]*rfavg_sinth_dph[j][i]));

                for (s = 7; s < 12; s++) {
                    fsf_imjk = a_imjk_multi_a_imjk_subtr*(xx[k][j][i].fx[s]-xx[k][j][im].fx[s]);
                    fsf_ipjk = a_ipjk_multi_a_ipjk_subtr*(xx[k][j][ip].fx[s]-xx[k][j][i].fx[s]);

                    fsf_ijmk = b_ijmk_multi_b_ijmk_subtr*(xx[k][j][i].fx[s]-xx[kcm][jm][i].fx[s]);
                    fsf_ijpk = b_ijpk_multi_b_ijpk_subtr*(xx[kcp][jp][i].fx[s]-xx[k][j][i].fx[s]);

                    fsf_ijkm = c_ijkm_multi_c_ijkm_subtr*(xx[k][j][i].fx[s]-xx[km][j][i].fx[s]);
                    fsf_ijkp = c_ijkp_multi_c_ijkp_subtr*(xx[kp][j][i].fx[s]-xx[k][j][i].fx[s]);

                    ff[k][j][i].fx[s] += (rh2[ip]*fsf_ipjk - rh2[i]*fsf_imjk)/rh_d3[i]
                                        +(sinth_h[jp]*fsf_ijpk - sinth_h[j]*fsf_ijmk)/rfavg_costh[j][i]
                                        +(fsf_ijkp - fsf_ijkm)/rfavg_costh_dth_dph[j][i];
                }

                ff[k][j][i].fx[19] = dxdt[k][j][i].fx[19]-rhoe*sum_nueq*(uer-unr)-sum_rhonusq*(uir-unr);
                ff[k][j][i].fx[20] = dxdt[k][j][i].fx[20]-rhoe*sum_nueq*(uet-unt)-sum_rhonusq*(uit-unt);
                ff[k][j][i].fx[21] = dxdt[k][j][i].fx[21]-rhoe*sum_nueq*(uep-unp)-sum_rhonusq*(uip-unp);
                ff[k][j][i].fx[22] = dxdt[k][j][i].fx[22] + two3rd*Qnfric
                                    +2.0*rhos_nusq_msmq*(xx[k][j][i].fx[22]/Nn-xx[k][j][i].fx[10]/ne);

                for (s = 23; s< 26; s++) ff[k][j][i].fx[s] = dxdt[k][j][i].fx[s];
            }
        }
    }

    DMDAVecRestoreArray(da, localX, &xx);
    DMRestoreLocalVector(da, &localX);

    DMDAVecRestoreArray(da, localU, &uu);
    DMRestoreLocalVector(da, &localU);

    DMDAVecRestoreArray(da, localZ, &zz);
    DMRestoreLocalVector(da, &localZ);

    DMDAVecRestoreArray(da, Xdt, &dxdt);
    DMDAVecRestoreArray(da, F, &ff);

    return 0;
}