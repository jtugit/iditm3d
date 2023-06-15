#include <iostream>

using namespace std;

#include "param.h"

int check_positivity(DM da, Field ***xx)
{
    int i, j, k, s, ngnp=0, xs, ys, zs, xm, ym, zm;
    string  spec[4]={"O+", "H+", "He+", "electron"};

    DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);

//check if any negative densit or temperature
    for (k = zs; k < zs+zm; k++) {
        for (j = ys; j < ys+ym; j++) {
            if (j == 0 || j == Nth) continue;

            for (i = xs; i < xs+xm; i++) {
                for (s=0; s<a4; s++) {
                    if (isnan(xx[k][j][i].fx[s]) || isinf(xx[k][j][i].fx[s])) {
                        cout<<"Solution is Nan or inf at (i, j, k, s) = ("<<i<<", "<<j<<", "<<k<<", "<<s
                            <<")"<<endl;
                        exit(-1);
                    }
                }

                for (s = 0; s < 7; s++) {
                    if (xx[k][j][i].fx[s] <= 0.0) {
                        cout<<"Density of ion species "<<s<<" negative "<<xx[k][j][i].fx[s]
                            <<" at (i, j, k) = (" << i << ", " << j << ", "<< k << ")" <<endl;
                        ngnp=-3;
                    }
                }

                for (s = 16; s <= 19; s++) {
                    if (xx[k][j][i].fx[s] <= 0.0) {
                        cout<<spec[s-16]<<" temperature negative "<<xx[k][j][i].fx[s]
                            <<" at (i, j, k) = (" << i << ", " << j << ", "<< k << ")" <<endl;
                        ngnp=-3;
                    }
                }

                if (xx[k][j][i].fx[30] <= 0.0) {
                    cout<<"Neutral temperature negative "<<xx[k][j][i].fx[30]
                        <<" at (i, j, k) = (" << i << ", " << j << ", "<< k << ")" <<endl;
                    ngnp=-5;
                }
                if (ngnp < 0) exit(ngnp);
            }
        }
    }

    return ngnp;
}