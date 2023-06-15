#include <cmath>

#include "electric_field.h"
#include "funcdef.h"
#include "collision_freq.h"
#include "operators.h"

void solar_zenith(AppCtx *params, int ys, int ym, int zs, int zm);

/* ----------------------------------------------------------------------------
 *-------- Compute negative of \nabla \mathbf{q} for O+, H+, He+, and electron
  ----------------------------------------------------------------------------*/
inline double heat_flow_divergence(Field ***xx, Field ***localuu, int i, int j, int k, int xi, int yj, int s)
{
    double dlamdadr, dlamdadth, dlamdadph, dTdr, dTdth, dTdph;
    double d2Tdr2, d2Tdth2, d2Tdph2;
    int    im=i-1, ip=i+1, jm=j-1, jp=j+1, km=k-1, kp=k+1, kc, s3=s+3;

    //dlamdadr = difference_r(localuu, i, j, k, s3);
    //dlamdadth = difference_theta(localuu, i, j, k, s3);
    //dlamdadph = difference_phi(localuu, i, j, k, s3);

    //dTdr = difference_r(xx, i, j, k, s);
    //dTdth = difference_theta(xx, i, j, k, s);
    //dTdph = difference_phi(xx, i, j, k, s);

    dlamdadr=(localuu[k][j][ip].fx[s3]-localuu[k][j][im].fx[s3])/dr;
    dlamdadph=(localuu[kp][j][i].fx[s3]-localuu[km][j][i].fx[s3])/dph;

    dTdr=(xx[k][j][ip].fx[s]-xx[k][j][im].fx[s])/dr;
    dTdph=(xx[kp][j][i].fx[s]-xx[km][j][i].fx[s])/dph;

    d2Tdr2 = (xx[k][j][ip].fx[s]-2.0*xx[k][j][i].fx[s]+xx[k][j][im].fx[s])/(dr*dr);

    if (j==1) {
        kc=(k+a3/2) % a3;
        dlamdadth=(localuu[k][jp][i].fx[s3]-localuu[kc][j][i].fx[s3])/dth;
        dTdth=(xx[k][jp][i].fx[s]-xx[kc][j][i].fx[s])/dth;

        d2Tdth2=(xx[k][jp][i].fx[s]-2.0*xx[k][j][i].fx[s]+xx[kc][j][i].fx[s])/(dth*dth);
    }
    else if (j > 1 && j < Nthm) {
        dlamdadth=(localuu[k][jp][i].fx[s3]-localuu[k][jm][i].fx[s3])/dth;
        dTdth=(xx[k][jp][i].fx[s]-xx[k][jm][i].fx[s])/dth;

        d2Tdth2=(xx[k][jp][i].fx[s]-2.0*xx[k][j][i].fx[s]+xx[k][jm][i].fx[s])/(dth*dth);
    }
    else {
        kc=(k+a3/2) % a3;
        dlamdadth=(localuu[kc][j][i].fx[s3]-localuu[k][jm][i].fx[s3])/dth;
        dTdth=(xx[kc][j][i].fx[s]-xx[k][jm][i].fx[s])/dth;

        d2Tdth2=(xx[kc][j][i].fx[s]-2.0*xx[k][j][i].fx[s]+xx[k][jm][i].fx[s])/(dth*dth);
    }

    d2Tdph2=(xx[kp][j][i].fx[s]-2.0*xx[k][j][i].fx[s]+xx[km][j][i].fx[s])/(dph*dph);

    double r2sin2=rsin[yj][xi]*rsin[yj][xi]*dph*dph, rr2=rr[i]*rr[i];

    double div_q = dlamdadr*dTdr+dlamdadth*dTdth/rr2+dlamdadph*dTdph/r2sin2
                  +localuu[k][j][i].fx[s3]*(2.0*dTdr/rr[i]+cot_div_r[yj][xi]/rr[i]*dTdth
                  +d2Tdr2 + d2Tdth2/rr2 + d2Tdph2/r2sin2);

    return div_q;
}

/* ----------------------------------------------------------------------------
 *-------- Compute negative of \nabla \mathbf{q} for neutrals
  ----------------------------------------------------------------------------*/
inline double neu_heat_flow_divergence(Field ***xx, Field ***localuu, int i, int j, int k, int xi, int yj)
{
    double dlamdadr, dlamdadth, dlamdadph, dTdr, dTdth, dTdph;
    double d2Tdr2, d2Tdth2, d2Tdph2;
    int    im=i-1, ip=i+1, jm=j-1, jp=j+1, km=k-1, kp=k+1, kc;

    //dlamdadr = difference_r(localuu, i, j, k, 23);       //\Delta{lambda_s}/dr
    //dlamdadth = difference_theta(localuu, i, j, k, 23);  //\Delta{lambda_s}/dtheta
    //dlamdadph = difference_phi(localuu, i, j, k, 23);    //\Delta{lambda_s}/dphi

    //dTdr = difference_r(xx, i, j, k, 30);               //Delta{T_s}/dr
    //dTdth = difference_theta(xx, i, j, k, 30);          //Delta{T_s}/dtheta
    //dTdph = difference_phi(xx, i, j, k, 30);            //Delta{T_s}/dphi

    dlamdadr=(localuu[k][j][ip].fx[23]-localuu[k][j][im].fx[23])/dr;
    dlamdadph=(localuu[kp][j][i].fx[23]-localuu[km][j][i].fx[23])/dph;

    dTdr=(xx[k][j][ip].fx[30]-xx[k][j][im].fx[30])/dr;
    dTdph=(xx[kp][j][i].fx[30]-xx[km][j][i].fx[30])/dph;

    d2Tdr2 = (xx[k][j][ip].fx[30]-2.0*xx[k][j][i].fx[30]+xx[k][j][im].fx[30])/(dr*dr);

    if (j==1) {
        kc=(k+a3/2) % a3;
        dlamdadth=(localuu[k][jp][i].fx[23]-localuu[kc][j][i].fx[23])/dth;
        dTdth=(xx[k][jp][i].fx[30]-xx[kc][j][i].fx[30])/dth;

        d2Tdth2=(xx[k][jp][i].fx[30]-2.0*xx[k][j][i].fx[30]+xx[kc][j][i].fx[30])/(dth*dth);
    }
    else if (j > 1 && j < Nthm) {
        dlamdadth=(localuu[k][jp][i].fx[23]-localuu[k][jm][i].fx[23])/dth;
        dTdth=(xx[k][jp][i].fx[30]-xx[k][jm][i].fx[30])/dth;

        d2Tdth2=(xx[k][jp][i].fx[30]-2.0*xx[k][j][i].fx[30]+xx[k][jm][i].fx[30])/(dth*dth);
    }
    else {
        kc=(k+a3/2) % a3;
        dlamdadth=(localuu[kc][j][i].fx[23]-localuu[k][jm][i].fx[23])/dth;
        dTdth=(xx[kc][j][i].fx[30]-xx[k][jm][i].fx[30])/dth;

        d2Tdth2=(xx[kc][j][i].fx[30]-2.0*xx[k][j][i].fx[30]+xx[k][jm][i].fx[30])/(dth*dth);
    }

    d2Tdph2=(xx[kp][j][i].fx[30]-2.0*xx[k][j][i].fx[30]+xx[km][j][i].fx[30])/(dph*dph);

    double r2sin2=rsin[yj][xi]*rsin[yj][xi]*dph*dph, rr2=rr[i]*rr[i];

    double div_q = dlamdadr*dTdr+dlamdadth*dTdth/rr2+dlamdadph*dTdph/r2sin2
                  +localuu[k][j][i].fx[23]*( 2.0*dTdr/rr[i]+cot_div_r[yj][xi]/rr[i]*dTdth
                                            +d2Tdr2 + d2Tdth2/rr2 + d2Tdph2/r2sin2);

    return div_q;
}

inline void Bspecial_sphereto_Bpshere(Field ***xx, Field ***uu, int i, int j, int k, 
    int xi, int yj, int zk)
{
    double    Bx, By, Bz;

    Bx=( Jiv11[zk][yj]*xx[k][j][i].fx[31]+Jiv12[zk][yj][xi]*xx[k][j][i].fx[32]
        +Jiv13[zk][yj][xi]*xx[k][j][i].fx[33]);
    By=( Jiv21[zk][yj]*xx[k][j][i].fx[31]+Jiv22[zk][yj][xi]*xx[k][j][i].fx[32]
        +Jiv23[zk][yj][xi]*xx[k][j][i].fx[33]);
    Bz=(Jiv31[yj]*xx[k][j][i].fx[31]+Jiv32[yj][xi]*xx[k][j][i].fx[32]);

    uu[k][j][i].fx[24]=(K11[zk][yj]*Bx+K21[zk][yj]*By+K31[yj]*Bz)/r2sintheta[yj][xi];
    uu[k][j][i].fx[25]=(K12[zk][yj]*Bx+K22[zk][yj]*By+K32[yj]*Bz)/r2sintheta[yj][xi];
    uu[k][j][i].fx[26]=(K13[zk]*Bx+K23[zk]*By)/r2sintheta[yj][xi];
}

/*----------------------------------------------------------------------
 * Calculate various parameters and store in vectors U, V, W, Z.
 * xx used is local array, corresponding to solution X.
 * This routine must be called before declaring local uu, vv, ww, & zz
 *---------------------------------------------------------------------*/
int parameters(DM da, Vec X, AppCtx *params)
{
    Vec      localX, localU;
    int      i, j, k, s, m, s3;
    PetscInt xs, ys, zs, xm, ym, zm;
    Field    ***xx, ***uu, ***localuu;

    double Te, Te12, ne, ni[7], nn[7], Tn, nimole, nO, nO2, nN2, nH, nHe;
    int    zk, yj, xi;
    double Td;
    double qn[5], nqd[4], Te2, nuss[7];

    const double n00=n0*1.0e-6;

    //------------- xx local array ---------------------
    DMGetLocalVector(da, &localX);
    DMGlobalToLocalBegin(da, X, INSERT_VALUES, localX);
    DMGlobalToLocalEnd(da, X, INSERT_VALUES, localX);
    DMDAVecGetArray(da, localX, &xx);

    /*------- uu - global array -------------*/
    DMDAVecGetArray(da, params->U, &uu);

    DMDAGetCorners(da, &xs ,&ys, &zs, &xm, &ym, &zm);

    solar_zenith(params, (int)ys, (int)ym, (int)zs, (int)zm);
    top_bc_vel(params, ys, ym, zs, zm);

    for (k=zs; k<zs+zm; k++) {
        zk=k-zs;

        for (j=ys; j<ys+ym; j++) {
            if (j == 0 || j == Nth) continue;

            yj=j-ys;

            for (i=xs; i<xs+xm; i++) {
                xi=i-xs;

                for (s=6; s< 10; s++) uu[k][j][i].fx[s]=0.0;
                uu[k][j][i].fx[34]=0.0;  //rhon - normalized neutral mass density
                for (s = 0; s < sl; s++) {
                    ni[s] = xx[k][j][i].fx[s];
                    uu[k][j][i].fx[6] += ni[s];   //normalized electron density

                    uu[k][j][i].fx[s+27]=ni[s];

                    nn[s] = exp(xx[k][j][i].fx[20+s]);
                    uu[k][j][i].fx[10] += nn[s];  //normalized neutral density
                    uu[k][j][i].fx[34] += nn[s]*ams[s];

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
                collision_freq(xx, ni, nn, uu[k][j][i].fx[6], i, j, k, xi, yj, zk, nqd);

                //convert B in special spherical components to spherical components
                Bspecial_sphereto_Bpshere(xx, uu, i, j, k, xi, yj, zk);

/*-----------------------------------------------------------------------------*/
//--------------- thermal conductivities
/*-----------------------------------------------------------------------------*/
                ne=uu[k][j][i].fx[6]*n00;   //electron density in cm^{-3}
                Te=xx[k][j][i].fx[19]*T0;
                Te12=sqrt(Te);
                Te2=Te*Te;

                //neutral O, O2, N2, H, and He density in cm^{-3}
                nO =nn[0]*n00; nO2=nn[3]*n00; nN2=nn[4]*n00; nH =nn[1]*n00; nHe=nn[2]*n00;

                qn[0]=1.1e-16*(1.0+5.7e-4*Te);           //O contribution
                qn[1]=2.2e-16*(1.0+0.036*Te12);          //O2
                qn[2]=2.82e-17*Te12*(1.0-1.21e-4*Te);    //N2
                qn[3]=5.47e-15*(1.0-1.35e-4*Te);         //H
                qn[4]=5.6e-16;                           //He

                /* normalized electron thermal conductivity */
                nqd[0]=nO*qn[0]+nO2*qn[1]+nN2*qn[2]+nH*qn[3]+nHe*qn[4];
                uu[k][j][i].fx[22]=1.233694e-11*Te2*Te12/(1.0+3.32e4*Te2/ne*nqd[0])/lamda0;

                /* calculate normalized O+, H+, He+ thermal conductivity Schunk & Nagy eq(5.167) */
                for (m = 0; m < 3; m++) {
                    Td=xx[k][j][i].fx[16+m]*T0; Td=Td*sqrt(Td);

                    nuss[m]=1.27*ni[m]*n0/(sqrt(ams[m])*Td); //Schunk & Nagy eq(4.142)
                    nqd[m]=1.0+1.25*nqd[m]/nuss[m];

                    uu[k][j][i].fx[19+m]=4.96682e-13*Td*xx[k][j][i].fx[16+m]/(sqrt(ams[m])*nqd[m])/lamda0;
                }

                /* normalized neutral thermal conductivity */
                Tn=xx[k][j][i].fx[30]*T0; Td=pow(Tn, 0.69);
                uu[k][j][i].fx[23]=( 7.59e-4*Td + (3.93e-4*Td+0.255e-4*Tn-9.27e-4)    //O, O2
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

    for (k=zs; k<zs+zm; k++) {
        zk=k-zs;

        for (j=ys; j<ys+ym; j++) {
            if (j == 0 || j == Nth) continue;

            yj=j-ys;

            for (i=xs; i<xs+xm; i++) {
                if (i == 0 || i == Nr) continue;

                xi=i-xs;

                //current density by Ampere's law
                currents(localuu, uu, i, j, k, yj, xi);

                //divergence of heat flows
                uu[k][j][i].fx[11]=heat_flow_divergence(xx, localuu, i, j, k, xi, yj, 16);
                uu[k][j][i].fx[12]=heat_flow_divergence(xx, localuu, i, j, k, xi, yj, 17);
                uu[k][j][i].fx[13]=heat_flow_divergence(xx, localuu, i, j, k, xi, yj, 18);
                uu[k][j][i].fx[14]=heat_flow_divergence(xx, localuu, i, j, k, xi, yj, 19);
                uu[k][j][i].fx[15]=neu_heat_flow_divergence(xx, localuu, i, j, k, xi, yj);
            }
        }
    }

    // restore local array/local vectors, and global array uu
    DMDAVecRestoreArray(da,localU,&localuu);
    DMRestoreLocalVector(da,&localU);
    DMDAVecRestoreArray(da,params->U,&uu);

    // restore local array and local vector
    DMDAVecRestoreArray(da, localX, &xx);
    DMRestoreLocalVector(da, &localX);

/*-----------------------------------------------------------------------------------------
 *----- smooth divergence of q and and current density
 *----- multidomensional Shapiro filter is used to conduct smoothing (Falissard, JCP 2013)
 -----------------------------------------------------------------------------------------*/
    smooth_multi_dim_U(da, params->U, 11, 18);

    return 0;
}
