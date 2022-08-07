/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <cmath>
using namespace std;

#include "hdf5.h"
#include "param.h"
#include "funcdef.h"

void hyperslab_set(PetscInt xs, PetscInt xm, PetscInt ys, PetscInt ym, PetscInt zs, 
        PetscInt zm, int d4, hsize_t dimsf[], hsize_t offset[], hsize_t count[], hsize_t block[])
{
    dimsf[0]=(hsize_t)a3;
    dimsf[1]=(hsize_t)a2;
    dimsf[2]=(hsize_t)a1;
    dimsf[3]=(hsize_t)d4;

    offset[0]=(hsize_t)zs;  /*start location in z dimension*/
    offset[1]=(hsize_t)ys;  /*start location in y dimension*/
    offset[2]=(hsize_t)xs;  /*start location in x dimension*/
    offset[3]=0;   /*start location 4th components dimension*/

    count[0]=1;    /*number of block in z dimension */
    count[1]=1;
    count[2]=1;
    count[3]=1;

    block[0]=(hsize_t)zm;   /*size of the hyperslab in z dimension */
    block[1]=(hsize_t)ym;
    block[2]=(hsize_t)xm;
    block[3]=(hsize_t)d4;
}

int hdf5parallelread(MPI_Comm comm, hsize_t *offset, hsize_t *count, hsize_t *block, char *fname, 
     char *dset_name, double *xdata)
{
    hid_t    acc_tpl;       /* file access template */
    hid_t    fid;           /* file and data space identifiers */
    hid_t    filespace;     /* file data space identifiers */
    hid_t    memspace;      /* memory data space ID */
    hid_t    dataset;       /* data set ID */
    hid_t    xfer_plt;      /* data set transfer property list */
    herr_t   hstatus;       /* generic return value */

    MPI_Info info=MPI_INFO_NULL;

    /* -------------------------
     * OPEN AN HDF5 FILE
     * ------------------------- */
    /* set up file access template with parallel I/O access */
    acc_tpl = H5Pcreate(H5P_FILE_ACCESS);

    /* set parallel access with communicator */
    hstatus = H5Pset_fapl_mpio(acc_tpl, comm, info);

    /* open the file collectively and release property list ID */
    fid = H5Fopen(fname, H5F_ACC_RDONLY, acc_tpl);
    H5Pclose(acc_tpl);

    /* ------------------------------------
     * OPEN THE DATA SET IN THE HDF5 FILE
     * ------------------------------------ */
    /* open the data set collectively */
    dataset = H5Dopen2(fid, dset_name, H5P_DEFAULT);

    /* get file data space independently. Each process takes a hyper slab */
    filespace = H5Dget_space(dataset);
    hstatus = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, block);

    /* create a memory data space independently */
    memspace = H5Screate_simple(data_dim, block, NULL);

    /* set up collective transfer properties lis */
    xfer_plt = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(xfer_plt, H5FD_MPIO_INDEPENDENT);

    /* read data collectively */
    hstatus = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, filespace, xfer_plt, xdata);
    if (hstatus<0) return -1;

    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Dclose(dataset);
    H5Pclose(xfer_plt);
    H5Fclose(fid);

    return 0;
}

int hdf5parallelwrite(MPI_Comm comm, hsize_t *dimsf, hsize_t *offset, 
    hsize_t *count, hsize_t *block, char *fname, char *dset_name, double *xdata)
{
    /* HDF5 APIs definitions */
    hid_t    file_id, dset_id;    /* file and dataset identifiers */
    hid_t    dataspace, memspace; /* file and memory dataspace identifiers */
    hid_t    plist_id;            /* property list identifier */
    herr_t   hstatus;

    MPI_Info info=MPI_INFO_NULL;

    /*
    * set up file access property list with parallel I/O access
    */
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, comm, info);

    /*
     * create a new file collectively and release property list identifier
     */
    file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    /*
     * Create the dataspace for the dataset
     */
    dataspace = H5Screate_simple(data_dim, dimsf, NULL);

    /*
     * Each process select a hyperslab in the file, and defines the memory space
     * for writing it to the hyperslab in the file
     */
    memspace = H5Screate_simple(data_dim, block, NULL);

    /* Create the dataset */
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id, data_dim, block);

    /* define property list as compressed with level 6 */
    H5Pset_deflate(plist_id, 6);

    dset_id = H5Dcreate(file_id, dset_name, H5T_NATIVE_DOUBLE, dataspace, 
            H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    H5Sclose(dataspace);

    dataspace=H5Dget_space(dset_id);
    hstatus = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, 
            count, block);

    /* Create property list for collective dataset write and write buffer
     * to the file */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    hstatus=H5Dwrite(dset_id,H5T_NATIVE_DOUBLE,memspace,dataspace,plist_id,xdata);
    if (hstatus<0) return -1;

    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Pclose(plist_id);
    H5Dclose(dset_id);
    H5Fclose(file_id);

    return 0;
}