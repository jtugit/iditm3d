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
    int          i, j, k, kc;
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
    vsize=6*xm*ym*zm;
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

    /* calculate background magnetic field, Earth's rotation rate, centrifugal force */
    for (k = zs; k < zs+zm; k++) {
        zk=k-zs; 

        for (j=ys; j<ys+ym; j++) {
            if (j < 1 || j > Nthm) continue;

            yj=j-ys;

            for (i=xs; i<xs+xm; i++) {
                xi=i-xs;
                /* dipole magnetic field in geomagnetic coordinates at (r_i, theta_j, phi_k) */
                rr_unorm=rr[i]*r0;
                r3=rr_unorm*rr_unorm*rr_unorm;
                Br=-2.0*Mz*cos(theta[j])/r3;
                Bt=-Mz*sin(theta[j])/r3;
                Bp=0.0;

                uu[k][j][i].fx[0]=Br/B0; uu[k][j][i].fx[1]=Bt/B0; uu[k][j][i].fx[2]=Bp/B0;

                /* dipole magnetic field in Cartesian coordinates */
                B0x=(K11[zk][yj]*Br+K12[zk][yj]*Bt+K13[zk]*Bp);
                B0y=(K21[zk][yj]*Br+K22[zk][yj]*Bt+K23[zk]*Bp);
                B0z=(K31[yj]*Br+K32[yj]*Bt);

                uu[k][j][i].fx[3]=B0x/B0; uu[k][j][i].fx[4]=B0y/B0; uu[k][j][i].fx[5]=B0z/B0;

                s=6*(zk*ym*xm+yj*xm+xi);
                xdata[s]=rr[i]*r0; xdata[s+1]=theta[j]; xdata[s+2]=phi[k];
                xdata[s+3]=Br; xdata[s+4]=Bt; xdata[s+5]=Bp; //background mfd in Tesla

                /* convert spherical coordinates to Cartesian coordinates */
                SPHCAR_08(rr[i], theta[j], phi[k], xmag, ymag, zmag, 1);
                BCARSP_08(xmag, ymag, zmag, wxm, wym, wzm, wr, wt, wp);

                /*--- Earth's rotation rate in magnetic spherical coordinates */
                rotat_r[zk][yj]=wr; rotat_t[zk][yj]=wt; rotat_p[zk][yj]=wp;

                /*---- centrifugal force = Omega x (Omega x r) ----*/
                /*--- in rotating coordinates (for fluid needed to multiply rho) */
                cenf_r[zk][yj][xi]=-rr[i]*(wt*wt+wp*wp);
                cenf_t[zk][yj][xi]= rr[i]*wr*wt;
                cenf_p[zk][yj][xi]= rr[i]*wr*wp;
            }
        }
    }

    for (k = zs; k < zs+zm; k++) {
        zk=k-zs; kc = (k+a3/2) % a3;

        //boundary condition at j = 0
        if (ys == 0) {
            j=0; yj=0;

            for (i=xs; i<xs+xm; i++) {
                xi=i-xs;

                s=6*(zk*ym*xm+yj*xm+xi);
                xdata[s]=rr[i]*r0; xdata[s+1]=theta[j]; xdata[s+2]=phi[k];
                xdata[s+3]=uu[kc][1][i].fx[0]*B0;
                xdata[s+4]=uu[kc][1][i].fx[1]*B0;
                xdata[s+5]=uu[kc][1][i].fx[2]*B0; //background mfd in Tesla
            }
        }

        //boundary condition at j = Nth
        if (ys+ym == a2) {
            j=Nth; yj=j-ys;

            for (i=xs; i<xs+xm; i++) {
                xi=i-xs;

                s=6*(zk*ym*xm+yj*xm+xi);
                xdata[s]=rr[i]*r0; xdata[s+1]=theta[j]; xdata[s+2]=phi[k];
                xdata[s+3]=uu[kc][Nthm][i].fx[0]*B0;
                xdata[s+4]=uu[kc][Nthm][i].fx[1]*B0;
                xdata[s+5]=uu[kc][Nthm][i].fx[2]*B0;
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
