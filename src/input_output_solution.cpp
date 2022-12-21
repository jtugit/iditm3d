/*****************************************************************************
 *  output_final
 *   Output simulation results at given time step.
 *   Parallel output 
 *
 *   Jiannan Tu
 *   2/24/2017
 ****************************************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <cmath>
#include <stdlib.h>
using namespace std;

#include "hdf5.h"

#include "param.h"
#include "funcdef.h"

int input_psolutions(DM da, Field ***xx, AppCtx *params)
{
    int i, j, k, m, ierr, xs, xm, ys, ym, zs, zm, vsize, s, zk, yj, xi;
    uint64_t kj, kji, ji;
    double   *xdata;                 /* pointer to data buffer for writing */
    double   Bx, By, Bz, Bsph[3], Ex, Ey, Ez, Esph[3];
    string   spec[4]={"O+", "H+", "He+", "electron"};

    hsize_t  dimsf[data_dim];     /* dataset dimensions */
    hsize_t  offset[data_dim];    /* start location in each dimension */
    hsize_t  count[data_dim];     /* number of hyperslab dimension */
    hsize_t  block[data_dim];     /* number of blocks in each dimension */

    ierr=DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

    vsize=a4*xm*ym*zm;
    xdata=new double[vsize];

    hyperslab_set(xs, xm, ys, ym, zs, zm, a4, dimsf, offset, count, block);

    string strfname = params->workdir + "/inp/" + params->prefln;
    char *fname = &strfname[0];
    hdf5parallelread(MPI_COMM_WORLD, offset, count, block, fname, params->dset_name, xdata);

    for (k=zs; k<zs+zm; k++) {
        zk=k-zs;

        for (j=ys; j<ys+ym; j++) {
            yj=j-ys; kj = (uint64_t)(zk*ym+yj);

            for (i=xs; i<xs+xm; i++) {
                xi=i-xs; kji=(uint64_t)(zk*ym*xm+yj*xm+xi); ji=(uint64_t)(yj*xm+xi);
                s=a4*(int)kji;

                //ion density (O+, H+, He+, O2+, N2+, NO+, N+), ion mass density, and electron density
                for (m = 0; m < sl; m++) xx[k][j][i].fx[m] = log(xdata[s+m]*1.0e6/n0); //density in m^{-3} converted to log scale 

                //ion velocity in m/s, ion and electron temperatuer in K 
                for (m = sl; m < 16; m++) xx[k][j][i].fx[m] = xdata[s+m]/v0;
                for (m = 16; m < 20; m++) xx[k][j][i].fx[m] = xdata[s+m]/T0;

                /* O, H, He, O2, N2, NO, N m^{-3} in log scale */
                for (m = 20; m < 27; m++) xx[k][j][i].fx[m]=log(xdata[s+m]*1.0e6/n0); //density m^{-3} in log scale

                //neutral velocity unr, untheta, unphi (m/s) and temperature in K
                for (m = 27; m < 30; m++) xx[k][j][i].fx[m] = xdata[s+m]/v0;
                xx[k][j][i].fx[30] = xdata[s+30]/T0;
            
                /* normalized delta_B in spherical coordinates */
                for (m = 31; m < 34; m++) Bsph[m-31]=xdata[s+m]*1.0e-9/B0;

                //converte to Cartesian coordinates
                Bx=Kmat.K11[kj]*Bsph[0]+Kmat.K12[kj]*Bsph[1]+Kmat.K13[(uint64_t)zk]*Bsph[2];
                By=Kmat.K21[kj]*Bsph[0]+Kmat.K22[kj]*Bsph[1]+Kmat.K23[(uint64_t)zk]*Bsph[2];
                Bz=Kmat.K31[(uint64_t)yj]*Bsph[0]+Kmat.K32[(uint64_t)yj]*Bsph[1];

                /* normalized delta_B in special spherical coordinates */
                xx[k][j][i].fx[31]=r2sintheta[yj][xi]*(Jmat.J11[kj]*Bx+Jmat.J12[kj]*By+Jmat.J13[(uint64_t)yj]*Bz);
                xx[k][j][i].fx[32]=r2sintheta[yj][xi]*(Jmat.J21[kji]*Bx+Jmat.J22[kji]*By+Jmat.J23[ji]*Bz);
                xx[k][j][i].fx[33]=r2sintheta[yj][xi]*(Jmat.J31[kji]*Bx+Jmat.J32[kji]*By);

                /* normalized delta_B in spherical coordinates */
                for (m = 34; m < 37; m++) Esph[m-34]=xdata[s+m]*1.0e-3/E0;

                //converte to Cartesian coordinates
                Ex=Kmat.K11[kj]*Esph[0]+Kmat.K12[kj]*Esph[1]+Kmat.K13[(uint64_t)zk]*Esph[2];
                Ey=Kmat.K21[kj]*Esph[0]+Kmat.K22[kj]*Esph[1]+Kmat.K23[(uint64_t)zk]*Esph[2];
                Ez=Kmat.K31[(uint64_t)yj]*Esph[0]+Kmat.K32[(uint64_t)yj]*Esph[1];

                /* normalized E in special spherical coordinates */
                xx[k][j][i].fx[34]=(Jinv.Jiv11[kj]*Ex+Jinv.Jiv21[kj]*Ey+Jinv.Jiv31[(uint64_t)yj]*Ez);
                xx[k][j][i].fx[35]=(Jinv.Jiv12[kji]*Ex+Jinv.Jiv22[kji]*Ey+Jinv.Jiv32[ji]*Ez);
                xx[k][j][i].fx[36]=(Jinv.Jiv13[kji]*Ex+Jinv.Jiv23[kji]*Ey);

/*----- check if any variables ar not a number (Nan) or infinity (inf) or ------
 *----- negative density or temperature -----------*/
                for (m=0; m<a4; m++) {
                    if (isnan(xx[k][j][i].fx[m]) || isinf(xx[k][j][i].fx[m])) {
                        cout<<"variable is Nan or inf at ("<<i<<", "<<j<<", "<<k<<", "<<m<<") from input file"<<endl;
                        MPI_Abort(MPI_COMM_WORLD,ierr);
                    }
                }
                for (s = 16; s < 20; s++) {
                    if (xx[k][j][i].fx[s] <= 0.0) {
                        cout<<"Negative or zero temperature of "<<spec[s]<<" = "<<xx[k][j][i].fx[s]
                            <<" at ("<<i<<", "<<j<<", "<<k<<")" <<" from input file"<<endl;
                        MPI_Abort(MPI_COMM_WORLD,ierr);
                    }
                }
                if (xx[k][j][i].fx[30] <= 0.0) {
                    cout<<"Negative or zero neutral temperature = "<<xx[k][j][i].fx[30] << " at ("<<i<< ", "
                        << j << ", " << k << ") from input file" << endl;
                    MPI_Abort(MPI_COMM_WORLD,ierr);
                }
            }
        }
    }

    delete[] xdata;

    return 0;
}

/*****************************************************************************
 *  output_solution
 *   Output simulation results 
 *
 *   Jiannan Tu
 *   12/11/2013, 1/11/2014, 3/26/2020
 ****************************************************************************/
int output_solution(DM da, Field ***xx, AppCtx *params)
{
    PetscInt xs, xm, ys, ym, zs, zm, vsize, s;
    int i, j, k, m, zk, yj, xi, ierr;
    uint64_t kj, kji, ji;
    double   Bx, By, Bz, Ex, Ey, Ez;
    double   *xdata;                 /* pointer to data buffer for writing */
    const double B00=B0*1.0e9, n00=n0*1.0e-6, E00=E0*1.0e-3;

    hsize_t  dimsf[data_dim];     /* dataset dimensions */
    hsize_t  offset[data_dim];    /* start location in each dimension */
    hsize_t  count[data_dim];     /* number of hyperslab dimension */
    hsize_t  block[data_dim];     /* number of blocks in each dimension */
    PetscMPIInt rank;

    ierr=DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

    vsize=a4*xm*ym*zm;
    xdata=new double[vsize];

    /* all quantities in SI unit */
    for (k=zs; k<zs+zm; k++) {
        zk=k-zs;

        for (j=ys; j<ys+ym; j++) {
            yj=j-ys; kj = (uint64_t)(zk*ym+yj);

            for (i=xs; i<xs+xm; i++) {
                xi=i-xs; kji=(uint64_t)(zk*ym*xm+yj*xm+xi); ji=(uint64_t)(yj*xm+xi);
                s=a4*(int)kji;

                //output ion density in cm^{-3}
                for (m = 0; m < sl; m++) xdata[s+m]=exp(xx[k][j][i].fx[m])*n00;

                //ion velocity components (m/s) and ion and ele temperatures (K)
                for (m = 7; m < 16; m++) xdata[s+m]=xx[k][j][i].fx[m]*v0;
                for (m = 16; m < 20; m++) xdata[s+m]=xx[k][j][i].fx[m]*T0;

                /* neutral number density (cm^{-3}), velocity (m/s), and temperature (K)*/
                for (m = 20; m < 27; m++) xdata[s+m]=exp(xx[k][j][i].fx[m])*n00;
                for (m = 27; m < 30; m++) xdata[s+m]=xx[k][j][i].fx[m]*v0;
                xdata[s+30]=xx[k][j][i].fx[30]*T0;

                //convert perturbation magentic field to Cartesian coordinates
                Bx=( Jinv.Jiv11[kj]*xx[k][j][i].fx[31]+Jinv.Jiv12[kji]*xx[k][j][i].fx[32]
                    +Jinv.Jiv13[kji]*xx[k][j][i].fx[33])/r2sintheta[yj][xi];
                By=( Jinv.Jiv21[kj]*xx[k][j][i].fx[31]+Jinv.Jiv22[kji]*xx[k][j][i].fx[32]
                    +Jinv.Jiv23[kji]*xx[k][j][i].fx[33])/r2sintheta[yj][xi];
                Bz=( Jinv.Jiv31[(uint64_t)yj]*xx[k][j][i].fx[31]+Jinv.Jiv32[ji]*xx[k][j][i].fx[32])/r2sintheta[yj][xi];

                /* perturbation magnetic field (nT) in spherical coordinates */
                xdata[s+31]=(Kmat.K11[kj]*Bx+Kmat.K21[kj]*By+Kmat.K31[(uint64_t)yj]*Bz)*B00;
                xdata[s+32]=(Kmat.K12[kj]*Bx+Kmat.K22[kj]*By+Kmat.K32[(uint64_t)yj]*Bz)*B00;
                xdata[2+33]=(Kmat.K13[(uint64_t)zk]*Bx+Kmat.K23[(uint64_t)zk]*By)*B00;

                //convert electric field to Cartesian coordinates
                Ex=( Jmat.J11[kj]*xx[k][j][i].fx[34]+Jmat.J21[kji]*xx[k][j][i].fx[35]
                    +Jmat.J31[kji]*xx[k][j][i].fx[36]);
                Ey=( Jmat.J12[kj]*xx[k][j][i].fx[34]+Jmat.J22[kji]*xx[k][j][i].fx[35]
                    +Jmat.J32[kji]*xx[k][j][i].fx[36]);
                Ez=( Jmat.J13[(uint64_t)yj]*xx[k][j][i].fx[34]+Jmat.J23[ji]*xx[k][j][i].fx[35]);

                /* E (mV/m) in spherical coordinates */
                xdata[s+34]=(Kmat.K11[kj]*Ex+Kmat.K21[kj]*Ey+Kmat.K31[(uint64_t)yj]*Ez)*E00;
                xdata[s+35]=(Kmat.K12[kj]*Ex+Kmat.K22[kj]*Ey+Kmat.K32[(uint64_t)yj]*Ez)*E00;
                xdata[2+36]=(Kmat.K13[(uint64_t)zk]*Ex+Kmat.K23[(uint64_t)zk]*Ey)*E00;
            }
        }
    }

    string strfname = params->outpdir + "/uvbenp" + to_string(params->ndt) + ".h5";
    char *fname = &strfname[0];

    hyperslab_set(xs, xm, ys, ym, zs, zm, a4, dimsf, offset, count, block);
    hdf5parallelwrite(MPI_COMM_WORLD, dimsf, offset, count, block, fname, params->dset_name, xdata);

    delete[] xdata;

    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);CHKERRQ(ierr);

    if (!rank) data_attribute(fname, params, 0);

    return ierr;
}

int data_attribute(char *fname, AppCtx *params, int dset)
{
    char     attr_str[55]="a1 a2 a3 Year Month Date sl sm UTsec F107 F107A Ap dt";
    int attr_grp1[8], ierr;
    double   attr_grp2[5];
    hid_t    file_id, dset_id, attr_dataspace, attr_id, attr_type;
    hsize_t  dims;
    char     dsetname[11];

    if (dset==0) strncpy(dsetname, params->dset_name, 11);
    else strncpy(dsetname, params->dset_diag, 11);

    file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
    dset_id = H5Dopen2(file_id, dsetname, H5P_DEFAULT);

    attr_grp1[0]=a1;
    attr_grp1[1]=a2;
    attr_grp1[2]=a3;
    attr_grp1[3]=params->iyr;
    attr_grp1[4]=params->mon;
    attr_grp1[5]=params->idate;
    attr_grp1[6]=sl;
    attr_grp1[7]=sm;

    attr_grp2[0]=params->sec;
    attr_grp2[1]=params->f107;
    attr_grp2[2]=params->f107a;
    attr_grp2[3]=params->Ap;
    attr_grp2[4]=dt;

    /* create string attribute and write it to the hdf5 file*/
    dims=1;
    attr_dataspace=H5Screate_simple(1, &dims, NULL);
    attr_type=H5Tcopy(H5T_C_S1);
    ierr=H5Tset_size(attr_type, 49);
    ierr=H5Tset_strpad(attr_type,H5T_STR_NULLTERM);
    attr_id=H5Acreate(dset_id, "Parameter_Names", attr_type, attr_dataspace,
            H5P_DEFAULT, H5P_DEFAULT);
    ierr=H5Awrite(attr_id, attr_type, attr_str);
    H5Sclose(attr_dataspace);
    H5Aclose(attr_id);
        
    dims=8;
    attr_dataspace=H5Screate_simple(1, &dims, NULL);
    if (sizeof(int)==4) {
        attr_id=H5Acreate(dset_id, "Parameters_Group_1", H5T_NATIVE_INT, 
            attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
        ierr=H5Awrite(attr_id, H5T_NATIVE_INT, attr_grp1);
    }
    else {
        attr_id=H5Acreate(dset_id, "Parameters_Group_1", H5T_NATIVE_LONG, 
            attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
        ierr=H5Awrite(attr_id, H5T_NATIVE_LONG, attr_grp1);
    }
    H5Sclose(attr_dataspace);
    H5Aclose(attr_id);

    dims=5;
    attr_dataspace=H5Screate_simple(1, &dims, NULL);
    attr_id=H5Acreate(dset_id, "Parameters_Group_2", H5T_NATIVE_DOUBLE, 
            attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    ierr=H5Awrite(attr_id, H5T_NATIVE_DOUBLE, attr_grp2);
    H5Sclose(attr_dataspace);
    H5Aclose(attr_id);

    H5Dclose(dset_id);
    H5Fclose(file_id);

    return ierr;
}
