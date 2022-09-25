#include <iostream>

using namespace std;

#include "param.h"

int check_positivity(DM da, Field ***xx)
{
    int i, j, k, s, ngnp=0, xs, ys, zs, xm, ym, zm;
    string  disp[2]={"forward_scheme", "imex_leap_frog"};

    DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);

//check if any negative densit or temperature
    for (k = zs; k < zs+zm; k++) {
        for (j = ys; j < ys+ym; j++) {
            for (i = xs; i < xs+xm; i++) {
                if (j == 0 || j == Nth) continue;

                for (s=0; s<nvar; s++) {
                    if (isnan(xx[k][j][i].fx[s]) || isinf(xx[k][j][i].fx[s])) {
                        cout<<"Solution is Nan or inf at ("<<i<<", "<<j<<", "<<k<<", "<<s
                            <<") in rhsfunctions"<<endl;
                        exit(-1);
                    }
                }

                for (s = 0; s < sl; s++) {
                    if (xx[k][j][i].fx[s] <= 0.0) {
                        cout<<"Ion density <= 0 of species "<<s<<" = "<<xx[k][j][i].fx[s]
                            <<" at (i, j, k) = ("<<i<<", "<<j<<", "<<k<<")"<<endl;
                            ngnp=-1;
                    }
                    if (xx[k][j][i].fx[12+s] <= 0.0) {
                        cout<<"Neutral density <= 0 of species "<<s<<" = "<<xx[k][j][i].fx[12+s]
                            <<" at (i, j, k) = ("<<i<<", "<<j<<", "<<k<<")"<<endl;
                        ngnp=-2;
                    }
                }

                if (xx[k][j][i].fx[10] <= 0.0) {
                    cout<<"Ion pressure <= 0 "<<xx[k][j][i].fx[10]
                        <<" at (i, j, k) = (" << i << ", " << j << ", "<< k << ")" <<endl;
                    ngnp=-3;
                }

                if (xx[k][j][i].fx[11] <= 0.0) {
                    cout<<"Electron pressure <= 0 "<<xx[k][j][i].fx[11]
                        <<" at (i, j, k) = (" << i << ", " << j << ", "<< k << ")" <<endl;
                    ngnp=-4;
                }

                if (xx[k][j][i].fx[22] <= 0.0) {
                    cout<<"Neutral pressure <= 0 "<<xx[k][j][i].fx[22]
                        <<" at (i, j, k) = (" << i << ", " << j << ", "<< k << ")" <<endl;
                    ngnp=-5;
                }
                if (ngnp < 0) exit(ngnp);
            }
        }
    }

    return ngnp;
}