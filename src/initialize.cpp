/* -------------------------------------------------------------------------------
 * function initialize 
 *   Call house keeping routines to input parameters, set up initial conditions,
 *   specifically, input state variables from files and so on.
 *
 * Jiannan Tu 5/21/2022
 * -------------------------------------------------------------------------------*/
#include <string.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>

using namespace std;

#include <cmath>
#include "param.h"
#include "funcdef.h"

void nightflux(double ***, int, PetscInt, PetscInt);
void print_top_bc_vel(AppCtx *, int, int, int, int, int);
int input_iri_msis(DM da, Vec X, Field ***, AppCtx *params);

PetscErrorCode initialize(DM da, Vec X, AppCtx *params)
{
    Field     ***xx, ***xn;
    PetscInt  xs, ys, zs, xm, ym, zm;

    DMDAGetCorners(da, &xs ,&ys, &zs, &xm, &ym, &zm);

    /* allocate memory for arrays */
    array_allocate(xm, ym, zm);

    if(euvflux_seg(da, params)<0) return -1;

    const_normalize(params);
    
    /* calculate mesh system in geomagnetic spherical coordinates */
    if(grids(da, params) < 0) return -1;

    dipole_magnetic(da, params);

    //if (xs+xm == a1) print_top_bc_vel(params, xs+xm, ys, ym, zs, zm);

    params->ndt=params->npre;

    DMDAVecGetArray(da, X, &xx);
    DMDAVecGetArray(da, params->Xn, &xn);

    if (params->smod == 0) {
        /* initialize solution at time step 0 from IRI and MSIS model data */
        if (input_iri_msis(da, X, xx, params) < 0) exit(-1);

        if (output_solution(da, xx, params) < 0) exit(-1);

        cout<<"Finished reading IRI-MSIS data, which are stored in output/uvbenp0.h5"<<endl;

        return -1;
    }
    else {
    /* input values of state variables at time step 0 */
        if (input_psolutions(da, xx, params) < 0) exit(-1);

        if (output_solution(da, xx, params) < 0) exit(-1);
    }

    int aa=zm*ym*xm;
    copy(&xx[zs][ys][xs], &xx[zs][ys][xs]+aa, &xn[zs][ys][xs]);

    DMDAVecRestoreArray(da, X, &xx);
    DMDAVecRestoreArray(da, params->Xn, &xn);

    /*----- night time EUV flux interpolated to grid altitude and -------------*/
    /*----- integer zenith angle, from Stroble et al. [1974] ------------------*/
    nightflux(fluxn, 0, xs, xm);   //Lyman-alpha  1026
    nightflux(fluxn, 1, xs, xm);   //He I         584
    nightflux(fluxn, 2, xs, xm);   //He II        304
    nightflux(fluxn, 3, xs, xm);   //Lyman-beta   1216
    
    return 0;
}
