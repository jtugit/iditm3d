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
    Field        ***ww, ***zz;
    int          i, j, k;
    PetscInt     xs, xm, ys, ym, zs, zm;
    double       xmag, ymag, zmag;
    double       r3, Br=0.0, Bt=0.0, Bp=0.0;
    double       AAP[106],G[106],H[106],REC[106];
    const double Mz=8.0e15;   //dipole moment at the center of the Earth (Tesla m^3)
    char         dsetn[6]="B_BGD";

    double       wx, wy, wz, wxm, wym, wzm, wr, wt, wp;

    //int          file_free;
    //fstream      divbfs;
    //MPI_Status   status;
    //int     nproc, rank;

    hsize_t  dimsf[data_dim];     /* dataset dimensions */
    hsize_t  offset[data_dim];    /* start location in each dimension */
    hsize_t  count[data_dim];     /* number of hyperslab dimension */
    hsize_t  block[data_dim];     /* number of blocks in each dimension */
    double   *xdata;                 /* pointer to data buffer for writing */
    int vsize, zk, yj, xi, s;

    //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    /*
     Get pointers to vector data
    */
    DMDAVecGetArray(da, params->W, &ww);
    DMDAVecGetArray(da, params->Z, &zz);
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
    wx=0.0; wy=0.0, wz=w0;
    GEOMAG_08(wx, wy, wz, wxm, wym, wzm, 1, AAP);

    /* calculate background magnetic field, Earth's rotation rate, centrifugal force */
    for (k = zs; k < zs+zm; k++) {
        zk=k-zs;

        for (j=ys; j<ys+ym; j++) {
            yj=j-ys;

            for (i=xs; i<xs+xm; i++) {
                xi=i-xs;

                /* dipole magnetic field in geomagnetic coordinates at (rC_i, thetaC_j, phi_k) */
                // Note thetaC[Nth] = thetaC[Nth-1] on kc=(k+a3/2) % a3
                r3=rC[i]*rC[i]*rC[i];
                Br=-2.0*Mz*cos(thetaC[j])/r3;
                Bt=-Mz*sin(thetaC[j])/r3;  

                /* normalized dipole magnetic field in spherical coordinates */
                ww[k][j][i].fx[23]=Br;  //normalized
                ww[k][j][i].fx[24]=Bt;  //normalized

                Omegae[zk][yj][xi] = e*sqrt(Br*Br+Bt*Bt)/me;

                s=3*(zk*ym*xm+yj*xm+xi);
                xdata[s]  =Br;
                xdata[s+1]=Bt;
                xdata[s+2]=Bp;

                /* dipole magnetic field in geomagnetic coordinates at (<r>^{phi}_i, theta_j, phi_{km}) */
                // Note theta[Nth] = theta[Nth-1] on kc=(k+a3/2) % a3
                r3=rfavg[i]*rfavg[i]*rfavg[i];
                Br=-2.0*Mz*cos(theta[j])/r3;
                Bt=-Mz*sin(theta[j])/r3;
                ww[k][j][i].fx[25]=Br;
                zz[k][j][i].fx[25]=Bt;

                /* dipole magnetic field in geomagnetic coordinates at (<r>^{theta}_i, theta_{jm}, phi_k) */
                // Note thetah[0] = 0 & theta[Nth]=pi
                Br=-2.0*Mz*cos(thetah[j])/r3;
                Bt=-Mz*sin(thetah[j])/r3;
                zz[k][j][i].fx[23]=Br;
                zz[k][j][i].fx[24]=Bt;

                /* convert spherical coordinates to Cartesian coordinates */
                SPHCAR_08(rC[i], thetaC[j], phi[k], xmag, ymag, zmag, 1);
                BCARSP_08(xmag, ymag, zmag, wxm, wym, wzm, wr, wt, wp);

                /*--- normalized Earth's rotation rate in magnetic spherical coordinates */
                rotat_r[zk][yj]=wr;
                rotat_t[zk][yj]=wt;
                rotat_p[zk][yj]=wp;

                /*---- normalized centrifugal force = Omega x (Omega x r) ----*/
                /*--- in rotating coordinates (for fluid needed to multiply rho) */
                cenf_r[zk][yj][xi]=-rr[i]*(wt*wt+wp*wp);
                cenf_t[zk][yj][xi]= rr[i]*wr*wt;
                cenf_p[zk][yj][xi]= rr[i]*wr*wp;
            }
        }
    }

    /* release resources not needed anymore */
    DMDAVecRestoreArray(da, params->W, &ww);
    DMDAVecRestoreArray(da, params->Z, &zz);

    /* file name of hdf5 file and write buffer xdata to the file in parallel*/
    string strfname = params->outpdir + "/B0.h5";
    char *fname = &strfname[0];

    hyperslab_set(xs, xm, ys, ym, zs, zm, 3, dimsf, offset, count, block);

    hdf5parallelwrite(MPI_COMM_WORLD, dimsf, offset, count, block, fname, dsetn, xdata);

/* write text file --- diagnostic use only */
    /*strncpy(fname, params->outpdir, 150);
    strcat(fname, "/output/B0divB0.dat");

    if (rank==0) {
      file_free=1;

      divbfs.open(fname, fstream::out);
      if(!divbfs) {
        cout << "Can't open file " << fname << endl;
        MPI_Abort(comm,-1);
      }
    }
    else {
      MPI_Recv(&file_free,1,MPI_INT,rank-1,tag1,comm,&status);

      //open the file to append contents
      divbfs.open(fname);
      if(!divbfs) {        cout << "Can't open file " << fname << endl;
        MPI_Abort(comm,-1);
      }
      divbfs.seekp(0,fstream::end);
    }

    if(file_free==1) {
        for (k=zs; k<zs+zm; k++) {
            zk=k-zs;
            for (j=ys; j<ys+ym; j++) {
                yj=j-ys;
                for (i=xs; i<xs+xm; i++) {
                    xi=i-xs;
                    s=3*(zk*ym*xm+yj*xm+xi);
                    divbfs<<scientific<<setw(12)<<setprecision(4)
                          <<(rr[i]*r0-Re)/1.0e3
                          <<scientific<<setw(14)<<setprecision(5)<<xdata[s]
                          <<scientific<<setw(14)<<setprecision(5)<<xdata[s+1]
                          <<scientific<<setw(14)<<setprecision(5)<<xdata[s+2]<<endl;
                }
            }
        }
        divbfs.close();
    }

    if(rank != nprocs-1) MPI_Send(&file_free,1,MPI_INT,rank+1,tag1,comm);*/

    delete[] xdata;

    return;
}
