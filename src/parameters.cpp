#include <cmath>

#include "electric_field.h"
#include "funcdef.h"
#include "collision_freq.h"
#include "operators.h"

void solar_zenith(AppCtx *params, int ys, int ym, int zs, int zm);

/* ----------------------------------------------------------------------------
 *-------- Compute negative of \nabla \mathbf{q} for O+, H+, He+, and electron
 /* ----------------------------------------------------------------------------*/
inline double heat_flow_divergence(Field ***xx, Field ***localuu, int i, int j, int k, int xi, int yj, int s)
{
    double d2Tdth2;

    double dlamdadr = difference_r(localuu, i, j, k, s+7);
    double dlamdadth = difference_theta(localuu, i, j, k, s+7);
    double dlamdadph = difference_phi(localuu, i, j, k, s+7);

    double dTdr = difference_r(xx, i, j, k, s);
    double dTdth = difference_theta(xx, i, j, k, s);
    double dTdph = difference_phi(xx, i, j, k, s);

    double d2Tdr2 = (xx[k][j][i+1].fx[s]-2.0*xx[k][j][i].fx[s]+xx[k][j][i-1].fx[s])/(dr*dr);

    if (j==1) d2Tdth2=(xx[k][j+1][i].fx[s]-2.0*xx[k][j][i].fx[s]+xx[(k+a3/2) % a3][j][i].fx[s])/(dth*dth);
    if (j > 1 && j < Nthm) d2Tdth2=(xx[k][j+1][i].fx[s]-2.0*xx[k][j][i].fx[s]+xx[k][j-1][i].fx[s])/(dth*dth);
    else d2Tdth2=(xx[(k+a3/2) % a3][j][i].fx[s]-2.0*xx[k][j][i].fx[s]+xx[k][j-1][i].fx[s])/(dth*dth);

    double d2Tdph2=(xx[k+1][j][i].fx[s]-2.0*xx[k][j][i].fx[s]+xx[k][j][i].fx[s])/(dph*dph);

    double r2sin2=rsin[yj][xi]*rsin[yj][xi], rr2=rr[i]*rr[i];

    double div_q = dlamdadr*dTdr+dlamdadth*dTdth/rr2+dlamdadph*dTdph/r2sin2
                  +localuu[k][j][i].fx[s+7]*(2.0*dTdr/rr[i]+cot_div_r[yj][xi]/rr[i]*dTdth
                  +d2Tdr2 + d2Tdth2/rr2 + d2Tdph2/r2sin2);

    return div_q;
}

/* ----------------------------------------------------------------------------
 *-------- Compute negative of \nabla \mathbf{q} for neutrals
 /* ----------------------------------------------------------------------------*/
inline double neu_heat_flow_divergence(Field ***xx, Field ***localuu, int i, int j, int k, int xi, int yj)
{
    double d2Tdth2;

    double dlamdadr = difference_r(localuu, i, j, k, 27);
    double dlamdadth = difference_theta(localuu, i, j, k, 27);
    double dlamdadph = difference_phi(localuu, i, j, k, 27);

    double dTdr = difference_r(xx, i, j, k, 30);
    double dTdth = difference_theta(xx, i, j, k, 30);
    double dTdph = difference_phi(xx, i, j, k, 30);

    double d2Tdr2 = (xx[k][j][i+1].fx[30]-2.0*xx[k][j][i].fx[30]+xx[k][j][i-1].fx[30])/(dr*dr);

    if (j==1) d2Tdth2=(xx[k][j+1][i].fx[30]-2.0*xx[k][j][i].fx[30]+xx[(k+a3/2) % a3][j][i].fx[30])/(dth*dth);
    if (j > 1 && j < Nthm) d2Tdth2=(xx[k][j+1][i].fx[30]-2.0*xx[k][j][i].fx[30]+xx[k][j-1][i].fx[30])/(dth*dth);
    else d2Tdth2=(xx[(k+a3/2) % a3][j][i].fx[30]-2.0*xx[k][j][i].fx[30]+xx[k][j-1][i].fx[30])/(dth*dth);

    double d2Tdph2=(xx[k+1][j][i].fx[30]-2.0*xx[k][j][i].fx[30]+xx[k][j][i].fx[30])/(dph*dph);

    double r2sin2=rsin[yj][xi]*rsin[yj][xi], rr2=rr[i]*rr[i];

    double div_q = dlamdadr*dTdr+dlamdadth*dTdth/rr2+dlamdadph*dTdph/r2sin2
                  +localuu[k][j][i].fx[27]*(2.0*dTdr/rr[i]+cot_div_r[yj][xi]/rr[i]*dTdth
                  +d2Tdr2 + d2Tdth2/rr2 + d2Tdph2/r2sin2);

    return div_q;
}

/*----------------------------------------------------------------------
 * Calculate various parameters and store in vectors U, V, W, Z.
 * xx used is local array, corresponding to solution X.
 * This routine must be called before declaring local uu, vv, ww, & zz
 *---------------------------------------------------------------------*/
int parameters(DM da, Vec X, AppCtx *params)
{
    Vec      localX, localU;
    int      i, j, k, s, m, tt, s3;
    PetscInt xs, ys, zs, xm, ym, zm;
    Field    ***xx, ***uu, ***localuu;

    double Te, Te12, ne, ni[7], nn[7], Ti, Tn, fq, nimole;
    int    zk, yj, xi, s0;
    double Td, rhoi, rhon;
    double qn[5], nqd, Dst, amt, Te2, nuss[7], nuin, lambdai[7], lambdan, lambdae;
    double Br, Btheta, Bphi, uir, uitheta, uiphi, PiPe, unr, untheta, unphi;

    const double two3rdmu=2.0e-6/3.0, n00=n0*1.0e-6;

    //------------- xx local array ---------------------
    DMGetLocalVector(da, &localX);
    DMGlobalToLocalBegin(da, X, INSERT_VALUES, localX);
    DMGlobalToLocalEnd(da, X, INSERT_VALUES, localX);
    DMDAVecGetArray(da, localX, &xx);

    /*------- uu - global array -------------*/
    DMDAVecGetArray(da, params->U, &uu);

    DMDAGetCorners(da, &xs ,&ys, &zs, &xm, &ym, &zm);
    solar_zenith(params, (int)ys, (int)ym, (int)zs, (int)zm);

    for (k=zs; k<zs+zm; k++) {
        zk=k-zs;

        for (j=ys; j<ys+ym; j++) {
            if (j == 0 || j == Nth) continue;

            yj=j-ys;

            for (i=xs; i<xs+xm; i++) {
                xi=i-xs;

                for (s=6; s< 12; s++) uu[k][j][i].fx[s]=0.0;
                for (s = 0; s < sl; s++) {
                    uu[k][j][i].fx[28+s] = xx[k][j][i].fx[s];
                    uu[k][j][i].fx[6] += exp(xx[k][j][i].fx[s]);  //normalized electron density

                    nn[s] = exp(xx[k][j][i].fx[20+s]);
                    uu[k][j][i].fx[10] += nn[s];  //normalized neutral density
                    uu[k][j][i].fx[11] += nn[s]*ams[s]; //normalized neutral mass density

                    if (s < 3) {
                        s3=s*3;
                        uu[k][j][i].fx[7] += ni[s]*xx[k][j][i].fx[7+s3];
                        uu[k][j][i].fx[8] += ni[s]*xx[k][j][i].fx[8+s3];
                        uu[k][j][i].fx[9] += ni[s]*xx[k][j][i].fx[9+s3];
                    }
                }

                //normalized ion bulk velocity. Note that velocity of molecule ions is that of O+
                nimole=ni[3]+ni[4]+ni[5]+ni[6];
                uu[k][j][i].fx[7] = (uu[k][j][i].fx[7]+nimole*xx[k][j][i].fx[7])/uu[k][j][i].fx[6];
                uu[k][j][i].fx[8] = (uu[k][j][i].fx[8]+nimole*xx[k][j][i].fx[8])/uu[k][j][i].fx[6];
                uu[k][j][i].fx[9] = (uu[k][j][i].fx[9]+nimole*xx[k][j][i].fx[9])/uu[k][j][i].fx[6];

                /*----------- normazlied collision frequencies -----------*/
                collision_freq(xx, i, j, k, xi, yj, zk);

                //production and loss rates
                prod_loss_rates(xx, uu, i, j, k, zk, yj, xi);

                neu_cooling_rate(xx, uu, i, j, k);

                //current density by Ampere's law
                currents(xx, uu, i, j, k, yj, xi);

/*-----------------------------------------------------------------------------*/
//--------------- thermal conductivities
/*-----------------------------------------------------------------------------*/
                ne=uu[k][j][i].fx[6]*n00;   //electron density in cm^{-3}
                Te=uu[k][j][i].fx[19]*T0;
                Te12=sqrt(Te);
                Te2=Te*Te;

                nn[0]=exp(xx[k][j][i].fx[20])*n00;  //neutral density in cm^{-3}
                nn[1]=exp(xx[k][j][i].fx[23])*n00;
                nn[2]=exp(xx[k][j][i].fx[24])*n00;
                nn[3]=exp(xx[k][j][i].fx[21])*n00;
                nn[4]=exp(xx[k][j][i].fx[22])*n00;

                qn[0]=1.1e-16*(1.0+5.7e-4*Te);
                qn[1]=2.2e-16*(1.0+0.036*Te12);
                qn[2]=2.82e-17*Te12*(1.0-1.21e-4*Te);
                qn[3]=5.47e-15*(1.0-1.35e-4*Te);
                qn[4]=5.6e-16;

                /* normalized electron thermal conductivity */
                nqd=nn[0]*qn[0]+nn[3]*qn[1]+nn[4]*qn[2]+nn[1]*qn[3]+nn[2]*qn[4];
                uu[k][j][i].fx[26]=1.233694e-11*Te2*Te12/(1.0+3.32e4*Te2/ne*nqd)/lamda0;

                amt=uu[k][j][i].fx[11]/uu[k][j][i].fx[10];  //averaged neutral mass in amu

                /* calculate normalized O+, H+, He+ thermal conductivity Schunk & Nagy eq(5.167) */
                for (m = 0; m < 3; m++) {
                    s0=14*m;

                    Td=xx[k][j][i].fx[16+m]*T0; Td=Td*sqrt(Td);

                    nqd=0.0;
                    for (s = 0; s < sl+sm; s++) {
                        // ion - electron Coulomb interaction
                        if (s==0) nqd += nust[zk][yj][xi][s0]*(3.0+1.5*ame/ams[m]);
                        else if (s > 0 && s < sl) {
                            // ion - ion Coulomb interactions Schunk & Nagy eq (4.141b)
                            if (s <= m) tt=s-1; else tt=s;
                            Dst=(3.0*ams[m]*ams[m]-0.2*ams[tt]*ams[tt]+0.1*ams[m]*ams[tt])
                                /((ams[m]+ams[tt])*(ams[m]+ams[tt])); 
                            nqd += nust[zk][yj][xi][s0+s]*(Dst+1.5*ams[tt]/(ams[m]+ams[tt]));
                        }
                        else {
                            // ion - neutral interaction Schunk & Nagy eq (4.147b)
                            Dst=(3.0*ams[m]*ams[m]+amt*amt+1.6*ams[m]*amt)/((ams[m]+amt)*(ams[m]+amt));
                            nqd += nust[zk][yj][xi][s0+s+7]*(Dst+1.5*amt/(ams[m]+amt));
                        }
                    }

                    nuss[m]=1.27*ni[m]*n0/(sqrt(ams[m])*Td); //Schunk & Nagy eq(4.142)
                    nqd=1.0+1.25*nqd/nuss[m];

                    uu[k][j][i].fx[23+m]=4.96682e-13*Td*xx[k][j][i].fx[16+m]/(sqrt(ams[m])*nqd)/lamda0;
                }

                /* normalized neutral thermal conductivity */
                Tn=uu[k][j][i].fx[30]*T0;
                Td=pow(Tn, 0.69);
                uu[k][j][i].fx[27]=( 7.59e-4*Td + (3.93e-4*Td+0.255e-4*Tn-9.27e-4)    //O, O2
                                    +(3.82e-4*Td+0.190e-4*Tn+5.14e-4)                 //N2
                                    +3.79e-3*Td +2.99e-3*Td)/lamda0;                  //H, He
            }
        }
    }

    DMDAVecRestoreArray(da,params->U,&uu);

    // get local array localuu with ghost values updated for each process
    DMGetLocalVector(da, &localU);
    DMGlobalToLocalBegin(da,params->U,INSERT_VALUES,localU);
    DMGlobalToLocalEnd(da,params->U,INSERT_VALUES,localU);
    DMDAVecGetArray(da, localU, &localuu);

    //also get global array again
    DMDAVecGetArray(da, params->U, &uu);

    //calcuation of heat fluxes needs local uu because non-local grids required
    //therefre, they are evaluated in separate loops
    for (k=zs; k<zs+zm; k++) {
        zk=k-zs;

        for (j=ys; j<ys+ym; j++) {
            if (j == 0 || j == Nth) continue;

            yj=j-ys;

            //j = Nth is the grids on the other side (different k) of south pole so skip.
            for (i=xs; i<xs+xm; i++) {
                if (i == 0 || i == Nr) continue;

                xi=i-xs;

                //calculate electron pressure gradient for special spheric components of e-field
                //xx and uu are local arrays
                //efd_gradPe = E_gradPe(xx, localuu, i, j, k, xi, yj, zk, xm, ym, zm);
                //uu[k][j][i].fx[12]=efd_gradPe.r; uu[k][j][i].fx[13]=efd_gradPe.t; uu[k][j][i].fx[14]=efd_gradPe.p;

                electric_field_vxB(xx, localuu, i, j, k, xi, yj, zk, xm, ym, zm);

                //negative of divergence of heat flows
                uu[k][j][i].fx[15]=heat_flow_divergence(xx, localuu, i, j, k, xi, yj, 16);
                uu[k][j][i].fx[16]=heat_flow_divergence(xx, localuu, i, j, k, xi, yj, 17);
                uu[k][j][i].fx[17]=heat_flow_divergence(xx, localuu, i, j, k, xi, yj, 18);
                uu[k][j][i].fx[18]=heat_flow_divergence(xx, localuu, i, j, k, xi, yj, 19);
                uu[k][j][i].fx[19]=neu_heat_flow_divergence(xx, localuu, i, j, k, xi, yj);
            }

            /*--------- simple numerical boundary condisitons -----------------*/
            if (xs == 0) {
                for (s = 12; s < 28; s++) {
                    uu[k][j][0].fx[s]=2.0*uu[k][j][1].fx[s]-uu[k][j][2].fx[s];
                }
            }
            if (xs+xm == a1) {
                for (s = 12; s < 28; s++) {
                    uu[k][j][Nr].fx[s]=2.0*uu[k][j][Nrm].fx[s]-uu[k][j][Nrm-1].fx[s];
                }
            }
        }
    }

    // restore local arrays and local vectors
    DMDAVecRestoreArray(da,localU,&localuu);
    DMRestoreLocalVector(da,&localU);

    DMDAVecRestoreArray(da,params->U,&uu);

    // restore local arrays and local vectors
    DMDAVecRestoreArray(da, localX, &xx);
    DMRestoreLocalVector(da, &localX);

    /*-----------------------------------------------------------------------------------------
     *----- smooth gradPe Catesian components (see note) for gradient Pe par tof efd,
     * ---- divergence of q, and current density
     *----- multidomensional Shapiro filter is used to conduct smoothing (Falissard, JCP 2013)
     -----------------------------------------------------------------------------------------*/
    /*VecScatter   Vsc;
    VecScatterCreateToAll(params->U, &Vsc, &localU);
    VecScatterBegin(Vsc, params->U, localU, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(Vsc, params->U, localU, INSERT_VALUES, SCATTER_FORWARD);

    DMDAVecGetArray(da, localU, &localuu);*/

    // get local array localuu with ghost values updated for each process
    DMGetLocalVector(da, &localU);
    DMGlobalToLocalBegin(da,params->U,INSERT_VALUES,localU);
    DMGlobalToLocalEnd(da,params->U,INSERT_VALUES,localU);
    DMDAVecGetArray(da, localU, &localuu);

    //also get global array again
    DMDAVecGetArray(da, params->U, &uu);

    int imin, imax, jm, jp, jm2, jp2, kcm, kcp, kcm2, kcp2;
    double sgn;

    if (xs == 0) imin=xs+2; else imin = xs;
    if (xs+xm==a1) imax=xs+xm-2; else imax=xs+xm;

    double delta2_r, delta2_th, delta2_ph, delta4_r, delta4_th, delta4_ph;

    for (k = zs ; k < zs+zm; k++) {
        for (j = ys; j < ys+ym; j++) {
            if (j==1) {kcm=(k+a3/2) % a3; jm=j; kcm2=kcm; jm2=1;}
            else {kcm=k; jm=j-1; kcm2=k; jm2=j-2;}

            if (j < Nth) {kcp=k; jp=j+1; kcp2=k; jp2=j+2;}
            else {kcp=(k+a3/2) % a3; jm=Nth; kcp2=kcp; jp2=Nthm;}

            for (i = imin; i < imax; i++) {
                for (s = 12; s < 23; s++) {
                    delta2_r = localuu[k][j][i+1].fx[s]-2.0*localuu[k][j][i].fx[s]+localuu[k][j][i-1].fx[s];
                    delta2_th =localuu[kcp][jp][i].fx[s]-2.0*localuu[k][j][i].fx[s]+localuu[kcm][jm][i].fx[s];
                    delta2_ph =localuu[k+1][j][i].fx[s]-2.0*localuu[k][j][i].fx[s]+localuu[k-1][j][i].fx[s];

                    delta4_r = localuu[k][j][i+2].fx[s]-4.0*(localuu[k][j][i+1].fx[s]+localuu[k][j][i-1].fx[s])
                              +6.0*localuu[k][j][i].fx[s]+localuu[k][j][i-2].fx[s];
                    delta4_th =localuu[kcp2][jp2][i].fx[s]-4.0*(localuu[kcp][jp][i].fx[s]+localuu[kcm][jm][i].fx[s])
                              +6.0*localuu[k][j][i].fx[s]+localuu[kcm2][jm2][i].fx[s];
                    delta4_ph =localuu[k+2][j][i].fx[s]-4.0*(localuu[k+1][j][i].fx[s]+localuu[k-1][j][i].fx[s])
                              +6.0*localuu[k][j][i].fx[s]+localuu[k-2][j][i].fx[s];

                    uu[k][j][i].fx[s]= localuu[k][j][i].fx[s]
                                      +( delta2_r*delta2_th+delta2_th*delta2_ph+delta2_r*delta2_ph
                                        -delta4_r-delta4_th-delta4_ph)/16.0+delta2_r*delta2_th*delta2_ph/64.0;
                }
            }
        }
    }

    //VecScatterDestroy(&Vsc);
    //DMRestoreLocalVector(da,&localU);
    //VecDestroy(&localU);

    // restore local arrays and local vectors
    DMDAVecRestoreArray(da,localU,&localuu);
    DMRestoreLocalVector(da,&localU);

    DMDAVecRestoreArray(da,params->U,&uu);

    return 0;
}
