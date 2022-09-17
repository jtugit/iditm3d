#include <cmath>

#include "electric_field.h"
#include "funcdef.h"
#include "collision_freq.h"

void solar_zenith(AppCtx *params, int ys, int ym, int zs, int zm);

/*----------------------------------------------------------------------
 * Calculate various parameters and store in vectors U, V, W, Z.
 * xx used is local array, corresponding to solution X.
 * This routine must be called before declaring local uu, vv, ww, & zz
 *---------------------------------------------------------------------*/

int parameters(DM da, Vec X, AppCtx *params)
{
    Vec      localX, localU;
    int      i, j, k, s, m, tt;
    PetscInt xs, ys, zs, xm, ym, zm;
    Field    ***xx, ***uu, ***localuu, ***vv, ***ww, ***zz;

    double Te, Te12, ne, ni[7], nn[7], Ti, Tn, fq;
    int    zk, yj, xi, s0;
    double Td, rhoi, rhon;
    double qn[5], nqd, Dst, amt, Te2, nuss[7], nuin, lambdai, lambdan;
    double Br, Btheta, Bphi, uir, uitheta, uiphi, PiPe, unr, untheta, unphi;
    double a_imjk[2], b_ijmk[2], c_ijkm[2];

    const double two3rdmu=2.0e-6/3.0;

    int    Nrm2=Nr-2;
    double rrb=(rC[0]-rC[1])/(rC[2]-rC[1]);
    double rrt=(rC[Nr]-rC[Nrm])/(rC[Nrm]-rC[Nrm2]);
    double y1, y2;

    //------------- xx local array ---------------------
    DMGetLocalVector(da, &localX);
    DMGlobalToLocalBegin(da, X, INSERT_VALUES, localX);
    DMGlobalToLocalEnd(da, X, INSERT_VALUES, localX);
    DMDAVecGetArray(da, localX, &xx);

    /*------- uu, vv, ww, zz are global arrays -------------*/
    DMDAVecGetArray(da, params->U, &uu);
    DMDAVecGetArray(da, params->V, &vv);
    DMDAVecGetArray(da, params->W, &ww);
    DMDAVecGetArray(da, params->Z, &zz);

    DMDAGetCorners(da, &xs ,&ys, &zs, &xm, &ym, &zm);
    solar_zenith(params, (int)ys, (int)ym, (int)zs, (int)zm);

    /* compute point quantities at rfavg, thetaC, phi_k */
    for (k=zs; k<zs+zm; k++) {
        zk=k-zs;

        for (j=ys; j<ys+ym; j++) {
            if (j == 0 || j == Nth) continue;

            yj=j-ys;

            for (i=xs; i<xs+xm; i++) {
                xi=i-xs;

                //cell averaged ne and Nn (~ those at rC, thetaC, phi_k)
                rhoi=0.0; rhon=0.0; uu[k][j][i].fx[17]=0.0; uu[k][j][i].fx[18]=0.0;
                for (s = 0; s < sl; s++) {
                    ni[s] = xx[k][j][i].fx[s];
                    rhoi += ms[s]*ni[s];          //ion mass density
                    uu[k][j][i].fx[17] += ni[s];  //electron density

                    nn[s] = xx[k][j][i].fx[12+s];
                    rhon += ms[s]*nn[s];          //neutral mass density
                    uu[k][j][i].fx[18] += nn[s];  //neutral density
                }

                //cell averaged uir=(rhoi*uir)/rhoi, uitheta=(rhoi*uitheta)/rhoi, uiphi=(rhoi*uiphi)/rhoi
                for (s = 7; s < 10; s++) {
                    uu[k][j][i].fx[s] = xx[k][j][i].fx[s]/rhoi;
                    uu[k][j][i].fx[s-4]=uu[k][j][i].fx[s]; //electron velocity without current term

                    //unr=(rhon*unr)/rhon, untheta=(rhon*untheta)/rhon, unphi=(rhon*unphi)/rhon
                    uu[k][j][i].fx[12+s] = xx[k][j][i].fx[12+s]/rhon;
                }

                //Ti = Pi/ne, Te = pe/ne (ne = ni), Tn = Pn/Nn
                uu[k][j][i].fx[10] = xx[k][j][i].fx[10]/(uu[k][j][i].fx[17]*kb); 
                uu[k][j][i].fx[11] = xx[k][j][i].fx[11]/(uu[k][j][i].fx[17]*kb);
                uu[k][j][i].fx[22] = xx[k][j][i].fx[22]/(uu[k][j][i].fx[18]*kb);

/*--------------------------------------------------*/
/*---------------- collision frequencies -----------*/
/*--------------------------------------------------*/
                collision_freq(xx, i, j, k, xi, yj, zk);

                //production and loss rates
                prod_loss_rates(xx, uu, i, j, k, zk, yj, xi);

                //calculate ratio of ion-neutral coolision freqency over electron gyro-frequency
                nuin = 0.0;
                for (s = 0; s < sl; s++) {
                    for (m = 0; m < sm; m++) nuin += nust[zk][yj][xi][s*14+21+m];
                }
                zz[k][j][i].fx[23]=nuin/Omegae[zk][yj][xi];

                ne=uu[k][j][i].fx[17];

                uir=uu[k][j][i].fx[7];
                uitheta=uu[k][j][i].fx[8];
                uiphi=uu[k][j][i].fx[9];

                unr=uu[k][j][i].fx[19];
                untheta=uu[k][j][i].fx[20];
                unphi=uu[k][j][i].fx[21];

                //cell averaged flux = that at (rC, thetaC, phi_k) of continuity equations
                for (s=0; s<sl; s++) {
                    vv[k][j][i].fx[s]=ni[s]*uir;
                    ww[k][j][i].fx[s]=ni[s]*uitheta;
                    zz[k][j][i].fx[s]=ni[s]*uiphi;

                    vv[k][j][i].fx[12+s]=nn[s]*unr;
                    ww[k][j][i].fx[12+s]=nn[s]*untheta;
                    zz[k][j][i].fx[12+s]=nn[s]*unphi;
                }

                //assume Br(rC[Nr],theta, phi) = Br[rC[Nr-1], theta, phi)]
                if (i==Nr) Br=reconstructed_Br(xx, i-1, j, k, rC[i-1], thetaC[j], phi[k])+ww[k][j][i].fx[23];
                else Br=reconstructed_Br(xx, i, j, k, rC[i], thetaC[j], phi[k])+ww[k][j][i].fx[23];
                Btheta=reconstructed_Btheta(xx, i, j, k, rC[i], thetaC[j], phi[k])+ww[k][j][i].fx[24];
                Bphi=reconstructed_Bphi(xx, i, j, k, rC[i], thetaC[j], phi[k]);

                PiPe=xx[k][j][i].fx[10]+xx[k][j][i].fx[11];

                //flux F_r(rC, thetaC, phi_k)
                vv[k][j][i].fx[7]= xx[k][j][i].fx[7]*uir+PiPe+0.5e-6*(Btheta*Btheta-Br*Br+Bphi*Bphi);
                vv[k][j][i].fx[8]= xx[k][j][i].fx[8]*uir-1.0e-6*Br*Btheta;
                vv[k][j][i].fx[9]= xx[k][j][i].fx[9]*uir-1.0e-6*Br*Bphi;
                vv[k][j][i].fx[10]=xx[k][j][i].fx[10]*uir; //heat flux has to be calculated and included later
                vv[k][j][i].fx[11]=xx[k][j][i].fx[11]*uu[k][j][i].fx[7];

                vv[k][j][i].fx[19]=xx[k][j][i].fx[19]*unr+xx[k][j][i].fx[22];
                vv[k][j][i].fx[20]=xx[k][j][i].fx[20]*unr;
                vv[k][j][i].fx[21]=xx[k][j][i].fx[21]*unr;
                vv[k][j][i].fx[22]=xx[k][j][i].fx[22]*unr;   //heat flux has to be calculated and included later

                //flux F_theta(rC, thetaC, phi_k)
                ww[k][j][i].fx[7]= xx[k][j][i].fx[7]*uitheta-1.0e-6*Br*Btheta;
                ww[k][j][i].fx[8]= xx[k][j][i].fx[8]*uitheta+PiPe+0.5e-6*(Br*Br-Btheta*Btheta+Bphi*Bphi);
                ww[k][j][i].fx[9]= xx[k][j][i].fx[9]*uitheta-1.0e-6*Btheta*Bphi;
                ww[k][j][i].fx[10]=xx[k][j][i].fx[10]*uitheta; //heat flux has to be calculated and included later
                ww[k][j][i].fx[11]=xx[k][j][i].fx[11]*uu[k][j][i].fx[8];

                ww[k][j][i].fx[19]=xx[k][j][i].fx[19]*untheta;
                ww[k][j][i].fx[20]=xx[k][j][i].fx[20]*untheta+xx[k][j][i].fx[22];
                ww[k][j][i].fx[21]=xx[k][j][i].fx[21]*untheta;
                ww[k][j][i].fx[22]=xx[k][j][i].fx[22]*untheta;   //heat flux has to be calculated and included later

                //flux F_phi(rC, thetaC, phi_k)
                zz[k][j][i].fx[7]= xx[k][j][i].fx[7]*uiphi-1.0e-6*Br*Bphi;
                zz[k][j][i].fx[8]= xx[k][j][i].fx[8]*uiphi-1.0e-6*Btheta*Bphi;
                zz[k][j][i].fx[9]= xx[k][j][i].fx[9]*uiphi+PiPe+0.5e-6*(Br*Br+Btheta*Btheta-Bphi*Bphi);
                zz[k][j][i].fx[10]=xx[k][j][i].fx[10]*uiphi;     //heat flux has to be calculated and included later
                zz[k][j][i].fx[11]=xx[k][j][i].fx[11]*uu[k][j][i].fx[9];

                zz[k][j][i].fx[19]=xx[k][j][i].fx[19]*unphi;
                zz[k][j][i].fx[20]=xx[k][j][i].fx[20]*unphi;
                zz[k][j][i].fx[21]=xx[k][j][i].fx[21]*unphi+xx[k][j][i].fx[22];
                zz[k][j][i].fx[22]=xx[k][j][i].fx[22]*unphi;  //heat flux has to be calculated and included later

/**-----------------------------------------------------------------------------------
 * -------------- cell volume averaged Electric field - reconstructed
 * -----------------------------------------------------------------------------------*/
                //xx local arrays but ww global array. uu global to store cell averaged efd at (rC, thetaC, phi_k)
                electric_field(xx, ww, uu, i, j, k);
            }
        }
    }

    // also get local array localuu with ghost values updated for each process
    DMGetLocalVector(da, &localU);
    DMGlobalToLocalBegin(da,params->U,INSERT_VALUES,localU);
    DMGlobalToLocalEnd(da,params->U,INSERT_VALUES,localU);
    DMDAVecGetArray(da, localU, &localuu);

    Vec localZ;
    Field ***localzz;
    DMGetLocalVector(da, &localZ);
    DMGlobalToLocalBegin(da,params->Z,INSERT_VALUES,localZ);
    DMGlobalToLocalEnd(da,params->Z,INSERT_VALUES,localZ);
    DMDAVecGetArray(da, localZ, &localzz);

    //calcuation of heat fluxes needs local uu because non-local grids required
    //therefre, they are evaluated in separate loops
    for (k=zs; k<zs+zm; k++) {
        zk=k-zs;

        for (j=ys; j<ys+ym; j++) {
            if (j == 0) continue;

            yj=j-ys;

            //j = Nth is the grids on the other side (different k) of south pole so skip.
            for (i=xs; i<xs+xm; i++) {
                if (j > 0 && j < Nth) {
                    xi=i-xs;

                    if (i > 0 && i < Nr) {
                        //maximum speed at cell interfaces
                        max_speed_r_face(xx, localzz, i, j, k, rh[i], thetaC[j], phi[k], a_imjk);
                        max_speed_theta_face(xx, localzz, i, j, k, rfavg[i], thetah[j], phi[k], b_ijmk);
                        max_speed_phi_face(xx, localzz, i, j, k, rfavg[i], theta[j], phih[k], c_ijkm);

                        uu[k][j][i].fx[3] = a_imjk[0];
                        uu[k][j][i].fx[4] = a_imjk[1];
                        uu[k][j][i].fx[5] = b_ijmk[0];
                        uu[k][j][i].fx[6] = b_ijmk[1];
                        uu[k][j][i].fx[15]= c_ijkm[0];
                        uu[k][j][i].fx[16]= c_ijkm[1];
                    }

/*-----------------------------------------------------------------------------*/
//--------------- thermal conductivities
/*-----------------------------------------------------------------------------*/
                    ne=uu[k][j][i].fx[17];
                    Te=uu[k][j][i].fx[11];
                    Te12=sqrt(Te);
                    Te2=Te*Te;

                    nn[0]=xx[k][j][i].fx[12];
                    nn[1]=xx[k][j][i].fx[15];
                    nn[2]=xx[k][j][i].fx[16];
                    nn[3]=xx[k][j][i].fx[13];
                    nn[4]=xx[k][j][i].fx[14];

                    qn[0]=1.1e-16*(1.0+5.7e-4*Te);
                    qn[1]=2.2e-16*(1.0+0.036*Te12);
                    qn[2]=2.82e-17*Te12*(1.0-1.21e-4*Te);
                    qn[3]=5.47e-15*(1.0-1.35e-4*Te);
                    qn[4]=5.6e-16;

                    fq=1.0/(1.0+sinh(zh[xi]*1.0e3/(rr[0]-Re)-1.0));

                    /* electron thermal conductivity */
                    nqd=nn[0]*qn[0]+nn[3]*qn[1]+nn[4]*qn[2]+nn[1]*qn[3]+nn[2]*qn[4];
                    zz[k][j][i].fx[24]=1.233694e-11*Te2*Te12/(1.0+3.32e4*Te2/ne*nqd);

                    rhon = 0.0;
                    for (s = 0; s < sm; s++) rhon += xx[k][j][i].fx[12+s]*ms[s];
                    amt=rhon/(uu[k][j][i].fx[18]*mp);  //averaged neutral mass in amu

                    Ti=uu[k][j][i].fx[10];
                    Td=Ti*sqrt(Ti);

                    /* calculate ion thermal conductivity & heat fluxes */
                    for (s = 12; s < 15; s++) uu[k][j][i].fx[s]=0.0; 
                    for (m = 0; m < sl; m++) {
                        s0=14*m;

                        nqd=0.0;
                        for (s = 0; s < sl+sm; s++) {
                            // ion - electron Coulomb interaction
                            if (s==0) nqd += nust[zk][yj][xi][s0]*(3.0+1.5*ame/ams[m]);
                            else if (s > 0 && s < sl) {
                                // ion - ion Coulomb interactions
                                if (s <= m) tt=s-1; else tt=s;
                                Dst=(3.0*ams[m]*ams[m]-0.2*ams[tt]*ams[tt]+0.1*ams[m]*ams[tt])
                                     /((ams[m]+ams[tt])*(ams[m]+ams[tt]));
                                nqd += nust[zk][yj][xi][s0+s]*(Dst+1.5*ams[tt]/(ams[m]+ams[tt]));
                            }
                            else {
                                // ion - neutral interaction
                                Dst=(3.0*ams[m]*ams[m]+amt*amt+1.6*ams[m]*amt)/((ams[m]+amt)*(ams[m]+amt));
                                nqd += nust[zk][yj][xi][s0+s+7]*(Dst+1.5*amt/(ams[m]+amt));
                            }
                        }

                        nuss[m]=1.27*ni[m]/(sqrt(ams[m])*Td);
                        nqd=1.0+1.25*nqd/nuss[m];

                        lambdai=4.96682e-13*Td/(sqrt(ams[m])*nqd);

                        // summ of all normalized ion heat conductivities in J m^-1 s^-1 K^-1
                        uu[k][j][i].fx[14] -=lambdai;
                    }

                    //ion heat flux vector
                    uu[k][j][i].fx[12] =uu[k][j][i].fx[14]*limited_slope_r(localuu, i, j, k, 10);
                    uu[k][j][i].fx[13] =uu[k][j][i].fx[14]*limited_slope_theta(localuu, i, j, k, 10)/rC[i];
                    uu[k][j][i].fx[14] =uu[k][j][i].fx[14]*limited_slope_phi(localuu, i, j, k, 10)/rCsinC[j][i];
            
                    /* neutral thermal conductivity */
                    //fq=1.0/(1.0+0.05*sinh(zh[xi]*1.0e3/(rr[0]-Re)-1.0));

                    Tn=uu[k][j][i].fx[22];
                    Td=pow(Tn, 0.69);
                    lambdan=( 7.59e-4*Td +1.0e-4*(3.93*Td+0.255*Tn-9.27)    //O, O2
                                         +1.0e-4*(3.82*Td+0.190*Tn+5.14)    //N2
                                         +3.79e-3*Td                        //H
                                         +2.99e-3*Td)*fq;                   //He

                    /* neutral heat flux. */
                    uu[k][j][i].fx[23] =-lambdan*limited_slope_r(localuu, i, j, k, 22);
                    uu[k][j][i].fx[24] =-lambdan*limited_slope_theta(localuu, i, j, k, 22)/rC[i];
                    uu[k][j][i].fx[25] =-lambdan*limited_slope_phi(localuu, i, j, k, 22)/rCsinC[j][i];

                    // now add heat fluxes to the ion and neutral pressure fluxes
                    vv[k][j][i].fx[10] += two3rdmu*uu[k][j][i].fx[12];
                    vv[k][j][i].fx[22] += two3rdmu*uu[k][j][i].fx[23];

                    ww[k][j][i].fx[10] += two3rdmu*uu[k][j][i].fx[13];
                    ww[k][j][i].fx[22] += two3rdmu*uu[k][j][i].fx[24];

                    zz[k][j][i].fx[10] += two3rdmu*uu[k][j][i].fx[14];
                    zz[k][j][i].fx[22] += two3rdmu*uu[k][j][i].fx[25];
                }
                else {
                    if (i > 0 && i < Nr) {
                        //maximum speed at cell interfaces j=Nth-1/2, thetah[Nth]=180 deg
                        max_speed_theta_face(xx, localzz, i, j, k, rfavg[i], thetah[j], phi[k], b_ijmk);

                        uu[k][j][i].fx[15]= c_ijkm[0];
                        uu[k][j][i].fx[16]= c_ijkm[1];
                    }
                }
            }


            //bottom BC for maximum/minimum speeds at r- and phi-intefaces
            if (j > 0 && j < Nth) {
                if (xs == 0) {
                    for (s = 3; s < 7; s++) {
                        y1= uu[k][j][1].fx[s];
                        y2= uu[k][j][2].fx[s];
                        uu[k][j][0].fx[s]= y1 + rrb*(y2-y1);
                    }
                }

                //top BC for maximum/minimum speeds at r- and phi-intefaces
                if (xs+xm == a1) {
                    for (s = 3; s < 7; s++) {
                        y1=uu[k][j][Nrm].fx[s];
                        y2=uu[k][j][Nrm2].fx[s];
                        uu[k][j][Nr].fx[s] = y1+rrt*(y1-y2);
                    }
                }
            }

            //bottom BC for maximum/minimum speeds at theta interfaces
            if (xs == 0) {
                for (s = 15; s < 17; s++) {
                    y1= uu[k][j][1].fx[s];
                    y2= uu[k][j][2].fx[s];
                    uu[k][j][0].fx[s]= y1 + rrb*(y2-y1);
                }
            }

            //top BC for maximum/minimum speeds at theta-interfaces
            if (xs+xm == a1) {
                for (s = 15; s < 17; s++) {
                    y1=uu[k][j][Nrm].fx[s];
                    y2=uu[k][j][Nrm2].fx[s];
                    uu[k][j][Nr].fx[s] = y1+rrt*(y1-y2);
                }
            }
        }
    }

    DMDAVecRestoreArray(da,params->U,&uu);
    DMDAVecRestoreArray(da,params->W,&ww);
    DMDAVecRestoreArray(da,params->Z,&zz);

    // restore local arrays and local vectors
    DMDAVecRestoreArray(da,localU,&localuu);
    DMRestoreLocalVector(da,&localU);

    DMDAVecRestoreArray(da,localZ,&localzz);
    DMRestoreLocalVector(da,&localZ);

    DMGetLocalVector(da, &localU);
    DMGlobalToLocalBegin(da,params->U,INSERT_VALUES,localU);
    DMGlobalToLocalEnd(da,params->U,INSERT_VALUES,localU);
    DMDAVecGetArray(da, localU, &localuu);

/**-----------------------------------------------------------------------------------
 * ---------- Electric field at edges - reconstructed and stored in Field vv
 * -----------------------------------------------------------------------------------*/
    for (k=zs; k<zs+zm; k++) {
        for (j=ys; j<ys+ym; j++) {
            if (j == 0) continue;

            for (i=xs; i<xs+xm; i++) {
                if (i == 0 || i == Nr) continue;

                Estar(xx, localuu, i, j, k, vv);
            }

            //bottom BC for edge-averaged electric field
            if (xs == 0) {
                for (s = 0; s < 3; s++) {
                    y1= vv[k][j][1].fx[s];
                    y2= vv[k][j][2].fx[s];
                    vv[k][j][0].fx[s]= y1 + rrb*(y2-y1);
                }
            }

            //top BC for edge-averaged electric field
            if (xs+xm == a1) {
                for (s = 0; s < 3; s++) {
                    y1=vv[k][j][Nrm].fx[s];
                    y2=vv[k][j][Nrm2].fx[s];
                    vv[k][j][Nr].fx[s] = y1+rrt*(y1-y2);
                }
            }
        }
    }

    DMDAVecRestoreArray(da,params->V,&vv);

    // restore local arrays and local vectors
    DMDAVecRestoreArray(da, localX, &xx);
    DMRestoreLocalVector(da, &localX);

    DMDAVecRestoreArray(da,localU,&localuu);
    DMRestoreLocalVector(da,&localU);

    return 0;
}
