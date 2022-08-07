#include <algorithm>
#include <cmath>
using namespace std;

#include "param.h"

/*--- Following code was adapted from source code of Sami2 simulation model ---*/

/*******************************************
*******************************************

!            f1026, f584, f304, f1216

*******************************************
*******************************************/

void nightflux(double ***f, int line, PetscInt xs, PetscInt xm)
{
    int    xi, i, j, k;
    int    k90, ji, ki, jip1, kip1;
    double delk, flog, atanhf, zz, z90, z93, z94, z170, z138;

//-- determine f for the 4 known values of theta
//-- emission flux in photons m^-2 s^-1

    for (i = xs; i < xs+xm; i++) {
        xi=i-xs;

        if ( zh[xi] < zaltnt[line][0] ) {
            for (k = 1; k <= 4; k++) f[xi][(int)(thetant[line][k])+1-90][line]=1.0;
        }
        else {
            if ( zaltnt[line][0] <= zh[xi] && zh[xi] <= zaltnt[line][1]) {
                zz=zh[xi];
            }
            else zz=zaltnt[line][1];

            if (line==0) {
                if (zz-90.0 < 0.0) z90=0.0;
                else z90=zz-90.0;
                if (zz-93.0 < 0.0) z93=0.0;
                else z93=zz-90.3;
                if (zz-94.0 < 0.0) z94=0.0;
                else z94=zz-94.0;
                f[xi][(int)(thetant[line][1])+1-90][line]=1.4e12*tanh(z90/50.);
                f[xi][(int)(thetant[line][2])+1-90][line]=3.8e11*tanh(z90/50.);
                f[xi][(int)(thetant[line][3])+1-90][line]=1.4e11*tanh(z93/55.);
                f[xi][(int)(thetant[line][4])+1-90][line]=9.2e10*tanh(z94/55.);
            }
            else if (line==1) {
                if (zz-170.0 < 0.0) z170=0.0;
                else z170=zz-170.0;
                f[xi][(int)(thetant[line][1])+1-90][line]=1.85e9*pow(z170, 1.20);
                f[xi][(int)(thetant[line][2])+1-90][line]=2.60e8*pow(z170, 1.25);
                f[xi][(int)(thetant[line][3])+1-90][line]=2.60e7*pow(z170, 1.20);
                f[xi][(int)(thetant[line][4])+1-90][line]=2.60e6*pow(z170, 1.20);
            }
            else if (line==2) {
                if (zz-138.0 < 0.0) z138=0.0;
                else z138=zz-138.0;
                atanhf=tanh(z138 / 80.);
                f[xi][(int)(thetant[line][1])+1-90][line] = 3.8e10*atanhf;
                f[xi][(int)(thetant[line][2])+1-90][line] = 3.0e10*atanhf;
                f[xi][(int)(thetant[line][3])+1-90][line] = 2.5e10*atanhf;
                f[xi][(int)(thetant[line][4])+1-90][line] = 2.5e10*atanhf;
            }
            else {
                f[xi][(int)(thetant[line][1])+1-90][line]
                        = 1.2e14*tanh((zz-80.)/50.) + 3.0e13;
                f[xi][(int)(thetant[line][2])+1-90][line]
                        = 4.0e13*tanh((zz-80.)/50.) +1.0e13;
                f[xi][(int)(thetant[line][3])+1-90][line]
                        = 2.0e13*tanh((zz-65.)/50.) +1.0e12;
                f[xi][(int)(thetant[line][4])+1-90][line]
                        = 1.5e13*tanh((zz-75.)/50.) +1.0e12;
            }
        }
    }

    for (k = 1; k <= 4; k++) {
        for (i = xs; i < xs+xm; i++) {
            xi=i-xs;
            f[xi][(int)(thetant[line][k])+1-90][line]
                    = max(1., f[xi][(int)(thetant[line][k])+1-90][line]);
        }
    }

/*--now interpolate to all values of theta (90 - 180) */
 
    if (line==0 || line==3) {
        for (k=1; k<=91; k++) {
            k90 = 90 + k - 1;
            ji  = 1;
            ki  = (int)thetant[line][1];
            for (j = 1; j <= 4; j++) {
                if (k90 > (int)thetant[line][j]) {
                    ji = j;
                    ki = (int)thetant[line][ji];
                }
            }
            jip1 = ji + 1;
            kip1 = (int)thetant[line][jip1];
            delk = (thetant[line][jip1] - thetant[line][ji]);
            for (i = xs; i < xs+xm; i++) {
                xi=i-xs;
                flog = log10(f[xi][ki+1-90][line])
                      + (k90 - ki) / delk 
                                   *( log10(f[xi][kip1+1-90][line])
                                     -log10(f[xi][ki+1-90][line])); 
                f[xi][k][line] = pow(10.0, flog);
            }
        }
    }
    else{
        for (k=1; k<=91; k++) {
            k90 = 90 + k - 1;
            ji  = 1;
            ki  = (int)thetant[line][1];
            for (j = 1; j <= 4; j++) {
                if (k90 > (int)thetant[line][j]) {
                    ji = j;
                    ki = (int)thetant[line][ji];
                }
            }
            if (ji != 4) {
                jip1 = ji + 1;
                kip1 = (int)thetant[line][jip1];
                delk = (thetant[line][jip1] - thetant[line][ji]);
                for (i = xs; i < xs+xm; i++) {
                    xi=i-xs;
                    flog = log10(f[xi][ki+1-90][line])
                          + (k90 - ki) / delk 
                                   *( log10(f[xi][kip1+1-90][line])
                                     -log10(f[xi][ki+1-90][line])); 
                    f[xi][k][line] = pow(10.0, flog);
                }
            }
            else {
                delk = 180.0 - thetant[line][ji];
                for (i = xs; i < xs+xm; i++) {
                    xi=i-xs;
                    flog = log10(f[xi][ki+1-90][line]) 
                          + (k90 - ki) / delk 
                                       * (log10(1.)-log10(f[xi][ki+1-90][line])); 
                    f[xi][k][line] = pow(10.0, flog);
                }
            }
        }
    }

    return;
}
