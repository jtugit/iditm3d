/*-----------------------------------------------------------------------------------------
 *----- multidomensional Shapiro filter is used to conduct smoothing (Falissard, JCP 2013)
 -----------------------------------------------------------------------------------------*/
#include "param.h"

void smooth_multi_dim(DM da, Vec U, int startIndex, int endIndex)
{
    Vec   localU;
    Field ***uu, ***localuu;
    int   xs, xm, ys, ym, zs, zm;
    int i, j, k, s, imin, imax, jm, jp, jm2, jp2, kcm, kcp, kcm2, kcp2;

    DMGetLocalVector(da, &localU);
    DMGlobalToLocalBegin(da, U, INSERT_VALUES, localU);
    DMGlobalToLocalEnd(da, U, INSERT_VALUES, localU);
    DMDAVecGetArray(da, localU, &localuu);

    //also get global array again
    DMDAVecGetArray(da, U, &uu);
    
    DMDAGetCorners(da, &xs ,&ys, &zs, &xm, &ym, &zm);

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
                for (s = startIndex; s <= endIndex; s++) {
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

    // restore local arrays and local vectors
    DMDAVecRestoreArray(da, localU, &localuu);
    DMRestoreLocalVector(da, &localU);

    DMDAVecRestoreArray(da, U, &uu);
}
