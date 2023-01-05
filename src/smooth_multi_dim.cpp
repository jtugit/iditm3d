/*-----------------------------------------------------------------------------------------
 *----- multidomensional Shapiro filter is used to conduct smoothing (Falissard, JCP 2013)
 -----------------------------------------------------------------------------------------*/
#include "param.h"

void smooth_multi_dim_U(DM da, Vec U, int startIndex, int endIndex)
{
    Vec   localU;
    Field ***uu, ***localuu;
    int   xs, xm, ys, ym, zs, zm;
    int i, j, k, s, imin, imax, jm, jp, jm2, jp2, kcm, kcp, kcm2, kcp2;

    DMGetLocalVector(da, &localU);
    DMGlobalToLocalBegin(da, U, INSERT_VALUES, localU);
    DMGlobalToLocalEnd(da, U, INSERT_VALUES, localU);
    DMDAVecGetArray(da, localU, &localuu);

    //also get global array
    DMDAVecGetArray(da, U, &uu);
    
    DMDAGetCorners(da, &xs ,&ys, &zs, &xm, &ym, &zm);

    if (xs == 0) imin=xs+2; else imin = xs;
    if (xs+xm==a1) imax=xs+xm-2; else imax=xs+xm;

    double delta2_r, delta2_th, delta2_ph, delta4_r, delta4_th, delta4_ph;
    double sgn0, sgn1, sgnm1, sgnm2;

    for (k = zs ; k < zs+zm; k++) {
        for (j = ys; j < ys+ym; j++) {
            if (j==0 || j == Nth) continue;

            jm=j-1; jm2=j-2; jp=j+1; jp2=j+2;
            kcm=k; kcm2=k; kcp=k; kcp2=k;
            if (j==1) {kcm=(k+a3/2) % a3; jm=j; kcm2=kcm; jm2=2;}
            else if (j==2) {kcm=k; jm=j-1; kcm2=(k+a3/2) % a3; jm2=1;}

            else if (j==Nthm-1) {kcp=k; jp=Nthm; kcp2=(k+a3/2) % a3; jp2=j;}
            else if (j==Nthm) {kcp=(k+a3/2) % a3; jm=j; kcp2=kcp; jp2=Nthm-1;}

            for (i = imin; i < imax; i++) {
                for (s = startIndex; s < endIndex; s++) {
                    sgn0=1.0; sgn1=1.0; sgnm1=1.0; sgnm2=1.0;
                    if (s == 17) {
                        if (j == 1) {sgn0=-1.0; sgn1=-1.0;}
                        else if (j == 2) sgn1=-1.0;
                        else if (j == Nthm-1) sgnm2=-1.0;
                        else if (j == Nthm) {sgnm1=-1.0; sgnm2=-1.0;}
                    }

                    delta2_r = localuu[k][j][i+1].fx[s]-2.0*localuu[k][j][i].fx[s]+localuu[k][j][i-1].fx[s];
                    delta2_th=sgnm1*localuu[kcp][jp][i].fx[s]-2.0*localuu[k][j][i].fx[s]+sgn0*localuu[kcm][jm][i].fx[s];
                    delta2_ph=localuu[k+1][j][i].fx[s]-2.0*localuu[k][j][i].fx[s]+localuu[k-1][j][i].fx[s];

                    delta4_r = localuu[k][j][i+2].fx[s]-4.0*(localuu[k][j][i+1].fx[s]+localuu[k][j][i-1].fx[s])
                              +6.0*localuu[k][j][i].fx[s]+localuu[k][j][i-2].fx[s];
                    delta4_th= sgnm2*localuu[kcp2][jp2][i].fx[s]-4.0*sgnm1*(localuu[kcp][jp][i].fx[s]
                              +sgn0*localuu[kcm][jm][i].fx[s])+6.0*localuu[k][j][i].fx[s]
                              +sgn1*localuu[kcm2][jm2][i].fx[s];
                    delta4_ph= localuu[k+2][j][i].fx[s]-4.0*(localuu[k+1][j][i].fx[s]+localuu[k-1][j][i].fx[s])
                              +6.0*localuu[k][j][i].fx[s]+localuu[k-2][j][i].fx[s];

                    uu[k][j][i].fx[s]= localuu[k][j][i].fx[s]
                                      +( delta2_r*delta2_th+delta2_th*delta2_ph+delta2_r*delta2_ph
                                        -delta4_r-delta4_th-delta4_ph)/16.0+delta2_r*delta2_th*delta2_ph/64.0;
                }
            }
        }
    }

    // restore local arrays and local vectors
    DMDAVecRestoreArray(da, localU, &localuu);
    DMRestoreLocalVector(da, &localU);

    DMDAVecRestoreArray(da, U, &uu);
}

/*-------- for smoothing soultions xx ------------------------------------*/
void smooth_multi_dim_X(DM da, Vec X)
{
    Vec   localX;
    Field ***xx, ***localxx;
    int   xs, xm, ys, ym, zs, zm;
    int i, j, k, s, imin, imax, jm, jp, jm2, jp2, kcm, kcp, kcm2, kcp2;

    DMGetLocalVector(da, &localX);
    DMGlobalToLocalBegin(da, X, INSERT_VALUES, localX);
    DMGlobalToLocalEnd(da, X, INSERT_VALUES, localX);
    DMDAVecGetArray(da, localX, &localxx);

    //also get global array again
    DMDAVecGetArray(da, X, &xx);
    
    DMDAGetCorners(da, &xs ,&ys, &zs, &xm, &ym, &zm);

    if (xs == 0) imin=xs+2; else imin = xs;
    if (xs+xm==a1) imax=xs+xm-2; else imax=xs+xm;

    double delta2_r, delta2_th, delta2_ph, delta4_r, delta4_th, delta4_ph;
    double sgn0, sgn1, sgnm1, sgnm2;

    for (k = zs ; k < zs+zm; k++) {
        for (j = ys; j < ys+ym; j++) {
            if (j==0 || j == Nth) continue;

            jm=j-1; jm2=j-2; jp=j+1; jp2=j+2;
            kcm=k; kcm2=k; kcp=k; kcp2=k;
            if (j==1) {kcm=(k+a3/2) % a3; jm=j; kcm2=kcm; jm2=2;}
            else if (j==2) {kcm=k; jm=j-1; kcm2=(k+a3/2) % a3; jm2=1;}

            else if (j==Nthm-1) {kcp=k; jp=Nthm; kcp2=(k+a3/2) % a3; jp2=j;}
            else if (j==Nthm) {kcp=(k+a3/2) % a3; jm=j; kcp2=kcp; jp2=Nthm-1;}

            for (i = imin; i < imax; i++) {
                for (s = 0; s < a4; s++) {
                    sgn0=1.0; sgn1=1.0; sgnm1=1.0; sgnm2=1.0;
                    if (s==8 || s==11 || s==14 || s==28 || s==32 || s==35) {
                        if (j == 1) {sgn0=-1.0; sgn1=-1.0;}
                        else if (j == 2) sgn1=-1.0;
                        else if (j == Nthm-1) sgnm2=-1.0;
                        else if (j == Nthm) {sgnm1=-1.0; sgnm2=-1.0;}
                    }

                    delta2_r = localxx[k][j][i+1].fx[s]-2.0*localxx[k][j][i].fx[s]+localxx[k][j][i-1].fx[s];
                    delta2_th=sgnm1*localxx[kcp][jp][i].fx[s]-2.0*localxx[k][j][i].fx[s]+sgn0*localxx[kcm][jm][i].fx[s];
                    delta2_ph=localxx[k+1][j][i].fx[s]-2.0*localxx[k][j][i].fx[s]+localxx[k-1][j][i].fx[s];

                    delta4_r = localxx[k][j][i+2].fx[s]-4.0*(localxx[k][j][i+1].fx[s]+localxx[k][j][i-1].fx[s])
                              +6.0*localxx[k][j][i].fx[s]+localxx[k][j][i-2].fx[s];
                    delta4_th= sgnm2*localxx[kcp2][jp2][i].fx[s]-4.0*sgnm1*(localxx[kcp][jp][i].fx[s]
                              +sgn0*localxx[kcm][jm][i].fx[s])+6.0*localxx[k][j][i].fx[s]
                              +sgn1*localxx[kcm2][jm2][i].fx[s];
                    delta4_ph= localxx[k+2][j][i].fx[s]-4.0*(localxx[k+1][j][i].fx[s]+localxx[k-1][j][i].fx[s])
                              +6.0*localxx[k][j][i].fx[s]+localxx[k-2][j][i].fx[s];

                    xx[k][j][i].fx[s]= localxx[k][j][i].fx[s]
                                      +( delta2_r*delta2_th+delta2_th*delta2_ph+delta2_r*delta2_ph
                                        -delta4_r-delta4_th-delta4_ph)/16.0+delta2_r*delta2_th*delta2_ph/64.0;
                }
            }
        }
    }

    // restore local arrays and local vectors
    DMDAVecRestoreArray(da, localX, &localxx);
    DMRestoreLocalVector(da, &localX);

    DMDAVecRestoreArray(da, X, &xx);
}
