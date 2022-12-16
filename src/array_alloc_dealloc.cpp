/*******************************************************************************
 *  array_alloc
 *  allocate memory for multi-dimension arrays
 * 
 *  Jiannan Tu
 *  5/21/2022
 ******************************************************************************/
#include <cmath>
using namespace std;

#include "param.h"
#include "funcdef.h"

void array_allocate(PetscInt xm, PetscInt ym, PetscInt zm)
{
    PetscInt  i, j, k;

//mean spherical coordinates, every process has a copy
    rr = new double[a1];
    rh = new double[a1+1];

    theta =new double[a2];
    thetah =new double[a2];

    phi = new double[a3];
    phih =new double[a3+1];

    zh    =new double[a1];

// gravity
    gr  = new double[xm];
    
// ---- solar zenith angles
    zenith = new double*[zm];

    vt = new double*[zm];
    vp = new double*[zm];

// ---- Earth's rotation vector in geomagnetic coordinates
    rotat_r =new double*[zm];
    rotat_t=new double*[zm];
    rotat_p=new double*[zm];
    cenf_r =new double**[zm];
    cenf_t=new double**[zm];
    cenf_p=new double**[zm];

    for (k = 0; k < zm; k++) {
        zenith[k]=new double[ym];

        vt[k] = new double[ym];
        vp[k] = new double[ym];

        rotat_r[k]=new double[ym];
        rotat_t[k]=new double[ym];
        rotat_p[k]=new double[ym];

        cenf_r[k]=new double*[ym];
        cenf_t[k]=new double*[ym];
        cenf_p[k]=new double*[ym];

        for (j = 0; j < ym; j++) {
            cenf_r[k][j]=new double[xm];
            cenf_t[k][j]=new double[xm];
            cenf_p[k][j]=new double[xm];
        }
    }

    r2sintheta = new double*[ym];
    cot_div_r = new double*[ym];
    rsin = new double*[ym];
    for (j = 0; j < ym; j++) {
        r2sintheta[j] =  new double[xm];
        cot_div_r[j] =  new double[xm];
        rsin[j] =  new double[xm];
    }

    grad_pe = new vector3D**[a3];
    for (k = 0; k < a3; k++) {
        grad_pe[k] = new vector3D*[a2];

        for (j = 0; j < a2; j++) grad_pe[k][j] = new vector3D[a2];
    }

    nust=new double***[zm];
    Ps = new double***[zm];
    Ls = new double***[zm];
    //Qee= new double**[zm];
    //Qeuv=new double**[zm];

    for (k = 0; k< zm; k++) {
        nust[k]=new double**[ym];

        Ps[k] = new double**[ym];
        Ls[k] = new double**[ym];

        //Qee[k]= new double*[ym];
        //Qeuv[k]=new double*[ym];

        for (j = 0; j < ym; j++) {
            nust[k][j]=new double*[xm];

            Ps[k][j] = new double*[xm];
            Ls[k][j] = new double*[xm];

            //Qee[k][j]= new double[xm];
            //Qeuv[k][j]=new double[xm];

            for (i = 0 ; i < xm; i++) {
                nust[k][j][i]=new double[112];

                Ps[k][j][i] = new double[14];
                Ls[k][j][i] = new double[14];
            }
        }
    }
 
    for (k=0; k<zm; k++) {
        for (j=0; j<ym; j++) {
            vt[k][j]=-1.0e6;
            vp[k][j]=-1.0e6;
        }
    }

    fluxn=new double**[xm];
    for (i =0 ; i < xm; i++) {
        fluxn[i]=new double*[92];
        for (j=0; j<92; j++) fluxn[i][j]=new double[4];
    }

    return;
}

void array_deallocate(PetscInt xm, PetscInt ym, PetscInt zm)
{
    PetscInt i, j, k;

    delete[] rr;
    delete[] rh;

    delete[] theta;
    delete[] thetah;

    delete[] phi;
    delete[] phih;

    delete[] zh;

    delete[] gr;

    for (k = 0; k < zm; k++) {
        delete[] zenith[k];

        delete[] vt[k];
        delete[] vp[k];

        for (j = 0; j < ym; j++) {
            delete[] cenf_r[k][j];
            delete[] cenf_t[k][j];
            delete[] cenf_p[k][j];
        }

        delete[] rotat_r[k];
        delete[] rotat_t[k];
        delete[] rotat_p[k];

        delete[] cenf_r[k];
        delete[] cenf_t[k];
        delete[] cenf_p[k];
    }
    delete[] zenith;
    delete[] vt;
    delete[] vp;

    delete[] rotat_r;
    delete[] rotat_t;
    delete[] rotat_p;
    delete[] cenf_r;
    delete[] cenf_t;
    delete[] cenf_p;

    for (j = 0; j < ym; j++) {
        delete[] r2sintheta[j]; delete[] cot_div_r[j]; delete[] rsin[j];
    }
    delete[] r2sintheta; delete[] cot_div_r; delete[] rsin;

    for (k = 0; k < a3; k++) {
        for (j = 0; j < a2; j++) delete[] grad_pe[k][j];
        delete[] grad_pe[k];
    }
    delete[] grad_pe;

    for (k = 0; k < zm; k++) {
        for (j = 0; j < ym; j++) {
            for (i = 0; i < xm; i++) {
                delete[] nust[k][j][i];

                delete[] Ps[k][j][i];
                delete[] Ls[k][j][i];
            }
            delete[] nust[k][j];

            delete[] Ps[k][j];
            delete[] Ls[k][j];
        }
        delete[] nust[k];

        delete[] Ps[k];
        delete[] Ls[k];
    }
    delete[] nust;
    delete[] Ps;
    delete[] Ls;

    for (i=0; i< xm; i++) {
        for (j=0; j<92; j++) delete[] fluxn[i][j];
        delete[] fluxn[i];
    }
    delete[] fluxn;
}
