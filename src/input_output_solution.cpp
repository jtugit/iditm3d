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
    int i, j, k, m, ierr;
    PetscInt xs, xm, ys, ym, zs, zm, vsize, s;
    int zk, yj, xi;
    double   *xdata;                 /* pointer to data buffer for writing */

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

    double rhoi, ne, rhon, Nn, ni00=ni_0/1.0e6, nn00=nn_0/1.0e6;

    for (k=zs; k<zs+zm; k++) {
        zk=k-zs;

        for (j=ys; j<ys+ym; j++) {
            yj=j-ys;

            for (i=xs; i<xs+xm; i++) {
                xi=i-xs;
                s=a4*(zk*ym*xm+yj*xm+xi);

                //ion density (O+, H+, He+, O2+, N2+, NO+, N+), ion mass density, and electron density
                rhoi = 0.0; ne = 0.0;
                for (m = 0; m < sl; m++) {
                    xx[k][j][i].fx[m] = xdata[s+m]/ni00; //density in m^{-3}/ni_0 
                    rhoi += ms[m]*xx[k][j][i].fx[m];     //ion mass density in kg m^{-3}/ni_0
                    ne += xx[k][j][i].fx[m];             //electron density in m^{-3}/ni_0
                }

                //momentum density (rhoi*uir), (rhoi*uitheta), (rhoi*uiphi) in kg/m^{3} (m/s) / ni_0 
                xx[k][j][i].fx[7] = rhoi*xdata[s+7];
                xx[k][j][i].fx[8] = rhoi*xdata[s+8];
                xx[k][j][i].fx[9] = rhoi*xdata[s+9];

                //pi=ni*kb*Ti and pe=ne*kb*Te (ni=ne) in J m^{-3} / ni_0
                xx[k][j][i].fx[10] = ne*kb*xdata[s+10];
                xx[k][j][i].fx[11] = ne*kb*xdata[s+11];

                /* O, H, He, O2, N2, NO, N normalized density in cm^{-3} */
                rhon = 0.0; Nn = 0.0;
                for (m = 0; m < sm; m++) {
                    xx[k][j][i].fx[12+m]=xdata[s+12+m]/nn00; //density in m^{-3}/nn_0
                    rhon += ms[m]*xx[k][j][i].fx[12+m];  //neutral mass density in kg m^{-3} / nn_0
                    Nn += xx[k][j][i].fx[12+m];          //neutral number density in m^{-3}/nn_0
                }

                //neutral momentum (rhon*unr), (rhon*untheta), (rhon*unphi) in kg/m^{3} (m/s) / nn_0
                xx[k][j][i].fx[19] = rhon*xdata[s+19];
                xx[k][j][i].fx[20] = rhon*xdata[s+20];
                xx[k][j][i].fx[21] = rhon*xdata[s+21];

                //pn=Nn*kb*Tn in J/m^{3} / nn_0
                xx[k][j][i].fx[22] = Nn*kb*xdata[s+22];
            
                /* delta_B at face center in Tesla */
                xx[k][j][i].fx[23]=xdata[s+23];
                xx[k][j][i].fx[24]=xdata[s+24];
                xx[k][j][i].fx[25]=xdata[s+25];

/*----- check if any variables ar not a number (Nan) or infinity (inf) or ------
 *----- negative density or temperature -----------*/
                for (m=0; m<nvar; m++) {
                    if (isnan(xx[k][j][i].fx[m]) || isinf(xx[k][j][i].fx[m])) {
                        cout<<"variable is Nan or inf at ("<<i<<", "<<j<<", "<<k
                            <<", "<<m<<") from input file"<<endl;
                        MPI_Abort(MPI_COMM_WORLD,ierr);
                    }
                }
                for (m = 0; m < sl; m++) {
                    if (xx[k][j][i].fx[m] <= 0.0) {
                        cout<<"Negative or zero ion density of species "<< m
                            <<" = "<<xx[k][j][i].fx[m] <<" at (" << i << ", " << j << ", "
                            << k << ") from input file" <<endl;
                        MPI_Abort(MPI_COMM_WORLD,ierr);
                    }
                }
                if (xx[k][j][i].fx[10] <= 0.0) {
                    cout<<"Negative or zero ion pressure = "<<xx[k][j][i].fx[10]
                        <<" at ("<<i<<", "<<j<<", "<<k<<")" <<" from input file"<<endl;
                    MPI_Abort(MPI_COMM_WORLD,ierr);
                }
                if (xx[k][j][i].fx[11] <= 0.0) {
                    cout<<"Negative or zero electron pressure = "<<xx[k][j][i].fx[11]
                        <<" at (i, j, k) = ("<<i<<", "<<j<<", "<<k<<")" <<" from input file"<<endl;
                    MPI_Abort(MPI_COMM_WORLD,ierr);
                }
                for (m=0; m<sm; m++) {
                    if (xx[k][j][i].fx[m+12] <= 0.0) {
                        cout<<"Negative or zero neutral density of species "<< m
                            <<" = "<<xx[k][j][i].fx[m+12] <<" at (" << i << ", " << j << ", "
                            << k << ") from input file" <<endl;
                        MPI_Abort(MPI_COMM_WORLD,ierr);
                    }
                }
                if (xx[k][j][i].fx[22] <= 0.0) {
                    cout<<"Negative or zero neutral pressure = "<<xx[k][j][i].fx[22] << " at ("<<i<< ", "
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
    double   *xdata;                 /* pointer to data buffer for writing */

    hsize_t  dimsf[data_dim];     /* dataset dimensions */
    hsize_t  offset[data_dim];    /* start location in each dimension */
    hsize_t  count[data_dim];     /* number of hyperslab dimension */
    hsize_t  block[data_dim];     /* number of blocks in each dimension */
    PetscMPIInt rank;

    ierr=DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

    vsize=a4*xm*ym*zm;
    xdata=new double[vsize];

    double rhoi, ne, rhon, Nn, ni00=ni_0/1.0e6, nn00=nn_0/1.0e6;;

    /* all quantities in SI unit */
    for (k=zs; k<zs+zm; k++) {
        zk=k-zs;

        for (j=ys; j<ys+ym; j++) {
            yj=j-ys;

            for (i=xs; i<xs+xm; i++) {
                xi=i-xs;
                s=a4*(zk*ym*xm+yj*xm+xi);

                /* ion number density (cm^{-3}), velocity (m/s), and temperature (K) */
                rhoi = 0.0; ne = 0.0;
                for (m = 0; m < sl; m++) {
                    xdata[s+m]=xx[k][j][i].fx[m]*ni00;  //ion density in cm^{-3}
                    rhoi += xx[k][j][i].fx[m]*ms[m];          //ion mass density in kg m^{-3}/ni_0
                    ne += xx[k][j][i].fx[m];                  //electron mass density in m^{-3}/ni_0
                }
                xdata[s+7]=xx[k][j][i].fx[7]/rhoi;  //uir     //velocity in m/s
                xdata[s+8]=xx[k][j][i].fx[8]/rhoi;  //uitheta
                xdata[s+9]=xx[k][j][i].fx[9]/rhoi;  //uiphi

                xdata[s+10] = xx[k][j][i].fx[10]/(ne*kb);     //temperature in K
                xdata[s+11] = xx[k][j][i].fx[11]/(ne*kb);

                /* neutral number density (cm^{-3}), velocity (m/s), and temperature (K)*/
                rhon = 0.0; Nn = 0.0;
                for (m = 0; m < sm; m++) {
                    xdata[s+12+m]=xx[k][j][i].fx[12+m]*nn00;
                    rhon += xx[k][j][i].fx[12+m]*ms[m];
                    Nn += xx[k][j][i].fx[12+m];
                }
                xdata[s+19]=xx[k][j][i].fx[19]/rhon;  //unr
                xdata[s+20]=xx[k][j][i].fx[20]/rhon;  //untheta
                xdata[s+21]=xx[k][j][i].fx[21]/rhon;  //unphi

                xdata[s+22] = xx[k][j][i].fx[22]/(Nn*kb);  //Tn

                /* perturbation magnetic field (Tesla) at face interfaces*/
                xdata[s+23]=xx[k][j][i].fx[23];
                xdata[s+24]=xx[k][j][i].fx[24];
                xdata[s+25]=xx[k][j][i].fx[25];
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
