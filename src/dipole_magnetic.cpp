/***************************************************************************** 
 *   dipole_magnetic
 *   Calculate background magnetic field components with magnetic meridian in
 *   curvilinear coordinates
 *
 *   Jiannan Tu
 *   5/21/2022
 *****************************************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
using namespace std;

#include "hdf5.h"
#include "param.h"
#include "funcdef.h"

#define tag1  12
void RECALC_08(int,int,int,int,int,double,double,double,double[],double[],double[],double[]);
void GEOMAG_08(double&, double&, double&, double&, double&, double&, int, double[]);
void SPHCAR_08(double&, double&, double&, double&, double&, double&, int);
void BCARSP_08(double, double, double, double, double,double, double&,double&, double&);

void dipole_magnetic(DM da, AppCtx *params)
{
    Field        ***uu;
    int          i, j, k;
    PetscInt     xs, xm, ys, ym, zs, zm;
    double       xmag, ymag, zmag;
    double       r3, Br=0.0, Bt=0.0, Bp=0.0;
    double       AAP[106],G[106],H[106],REC[106];
    const double Mz=8.0e15;   //dipole moment at the center of the Earth (Tesla m^3)
    char         dsetn[6]="B_BGD";

    double       wx, wy, wz, wxm, wym, wzm, wr, wt, wp;

    hsize_t  dimsf[data_dim];     /* dataset dimensions */
    hsize_t  offset[data_dim];    /* start location in each dimension */
    hsize_t  count[data_dim];     /* number of hyperslab dimension */
    hsize_t  block[data_dim];     /* number of blocks in each dimension */
    double   *xdata;                 /* pointer to data buffer for writing */
    int vsize, zk, yj, xi, s;
    double   B0x, B0y, B0z, rr_unorm;

    /*
     Get pointers to vector data
    */
    DMDAVecGetArray(da, params->U, &uu);
    DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);

    /* Initialize a one dimensional data buffer for storing 3-component B-field*/
    vsize=3*xm*ym*zm;
    xdata=new double[vsize];

    int iyr=params->iyr, mn=params->mon, dd=params->idate;
    int iday=dayno(params->iyr, mn, dd);
    int ihour=(int)((double)(params->sec)/3600.0);
    int imin=(int)((double)(params->sec-ihour*3600)/60.0);
    int isec=(int)(params->sec-ihour*3600-imin*60);

    RECALC_08(iyr,iday,ihour,imin,isec,params->vgsex,params->vgsey,params->vgsez,AAP,G,H,REC);

    /* convert Earth's rotation rate (normalized) to geomagnetic Cartesian coordinates*/
    wx=0.0; wy=0.0, wz=w0n;
    GEOMAG_08(wx, wy, wz, wxm, wym, wzm, 1, AAP);

    int kj, kji, ji;
    /* calculate background magnetic field, Earth's rotation rate, centrifugal force */
    for (k = zs; k < zs+zm; k++) {
        zk=k-zs;

        for (j=ys; j<ys+ym; j++) {
            if (j < 1) j = 1;
            if (j > Nthm) j = Nthm;

            yj=j-ys; kj=zk*ym+yj;

            for (i=xs; i<xs+xm; i++) {
                xi=i-xs; kji=zk*(ym*xm)+yj*xm+xi; ji=yj*xm+xi;

                /* dipole magnetic field in geomagnetic coordinates at (rC_i, thetaC_j, phi_k) */
                // Note thetaC[Nth] = thetaC[Nth-1] on kc=(k+a3/2) % a3
                rr_unorm=rr[i]*r0;
                r3=rr_unorm*rr_unorm*rr_unorm;
                Br=-2.0*Mz*cos(theta[j])/r3;
                Bt=-Mz*sin(theta[j])/r3;
                Bp=0.0;

                uu[k][j][i].fx[0]=Br/B0; uu[k][j][i].fx[1]=Bt/B0; uu[k][j][i].fx[2]=Bp/B0;

                /* dipole magnetic field in special spherical coordinates */
                B0x=(Kmat.K11[kj]*Br+Kmat.K12[kj]*Bt+Kmat.K13[zk]*Bp);
                B0y=(Kmat.K21[kj]*Br+Kmat.K22[kj]*Bt+Kmat.K23[zk]*Bp);
                B0z=(Kmat.K31[yj]*Br+Kmat.K32[yj]*Bt);

                uu[k][j][i].fx[28]=B0x/B0; uu[k][j][i].fx[29]=B0y/B0; uu[k][j][i].fx[30]=B0z/B0;

                s=3*(zk*ym*xm+yj*xm+xi);
                xdata[s]  =Br;
                xdata[s+1]=Bt;
                xdata[s+2]=Bp;

                /* convert spherical coordinates to Cartesian coordinates */
                SPHCAR_08(rr[i], theta[j], phi[k], xmag, ymag, zmag, 1);
                BCARSP_08(xmag, ymag, zmag, wxm, wym, wzm, wr, wt, wp);

                /*--- Earth's rotation rate in magnetic spherical coordinates */
                rotat_r[zk][yj]=wr;
                rotat_t[zk][yj]=wt;
                rotat_p[zk][yj]=wp;

                /*---- centrifugal force = Omega x (Omega x r) ----*/
                /*--- in rotating coordinates (for fluid needed to multiply rho) */
                cenf_r[zk][yj][xi]=-rr[i]*(wt*wt+wp*wp);
                cenf_t[zk][yj][xi]= rr[i]*wr*wt;
                cenf_p[zk][yj][xi]= rr[i]*wr*wp;
            }
        }
    }

    /* release resources not needed anymore */
    DMDAVecRestoreArray(da, params->U, &uu);

    /* file name of hdf5 file and write buffer xdata to the file in parallel*/
    string strfname = params->outpdir + "/B0.h5";
    char *fname = &strfname[0];

    hyperslab_set(xs, xm, ys, ym, zs, zm, 3, dimsf, offset, count, block);

    hdf5parallelwrite(MPI_COMM_WORLD, dimsf, offset, count, block, fname, dsetn, xdata);

    delete[] xdata;

    return;
}
