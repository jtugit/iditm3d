#include "ele_cooling_rate.h"
#include "collision_freq.h"

using namespace std;

// Gaussian elimination without pivot
/*inline void GaussianElimination(double mat[][10], int n, double res[])
{
    int    i, j, k;
    double f;

    for(i=0; i<n; i++) {                   
        for(j=i+1; j<n; j++) {
            if(abs(mat[i][i]) < abs(mat[j][i])) {
                for(k=0; k<n+1; k++) {
                    // swapping mat[i][k] and mat[j][k]
                    mat[i][k]=mat[i][k]+mat[j][k];
                    mat[j][k]=mat[i][k]-mat[j][k];
                    mat[i][k]=mat[i][k]-mat[j][k];
                }
            }
        }
    }
   
    // performing Gaussian elimination
    for(i=0; i<n-1; i++) {
        for(j=i+1; j<n; j++) {
            f=mat[j][i]/mat[i][i];
            for(k=0; k<n+1; k++) {
                mat[j][k]=mat[j][k]-f*mat[i][k];
            }
        }
    }

    // Backward substitution for discovering values of unknowns
    for(i=n-1; i>=0; i--)          {                     
        res[i]=mat[i][n];
                    
        for(j=i+1; j<n; j++) {
            if(i != j) res[i]=res[i]-mat[i][j]*res[j];
        }
  
        res[i]=res[i]/mat[i][i];  
    }
}*/

inline void GaussianElimination(double a[][4], int n, double res[])
{
    int    i, j, k;
    double multiplier;

    //------------ to get the diagonal matrix -------------------|
    for (k = 0;k < n;k++){
        for (i=0; i < n; i++){
            if(i!= k){
                multiplier = a[i][k]/a[k][k];
                for (j=0; j < n+1; j++){
                    a[i][j] = a[i][j] - a[k][j]*multiplier;
                }
            }
        }
    }

    // -------------------Solutions---------------------|
    for (i = 0;i < n;i++){
        res[i] = a[i][n]/a[i][i];
    }
}

inline void implicit_solver(Field ***xx, Field ***xn, Field ***xn1, int i, int j, int k, int xi, int yj, int zk)
{
    int  s, t;
    double vmat[6][7]={{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                       {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                       {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    double tmat[3][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
    double bb[6], cc[3];
    double rhos, rhoi, rhon, sum_nues, sum_nueq, sum_rhonusq, sum_nues_div_ms;
    double sum_nueq_div_mq, rhos_nusq_msmq, rhos_nusqmq_msmq, rhos_nusqms_msmq, temp;
    double ne, rhoe, uir_minus_unr, uit_minus_unt, uip_minus_unp, Nn; //ne = ni
    double Ce, dCedTe, dCedTn, Te, Tn, rhoe_sum_nues, rhoe_sum_nueq;
    double collr, collt, collp, coll_Tn_Ti;
    double uer, uet, uep, uir, uit, uip, unr, unt, unp;
    const double two3rd=2.0/3.0, two3rd_dt = 2.0/3.0*dt;

    rhoi=0.0; rhon=0.0;
    sum_nues=0.0; sum_nueq=0.0; sum_rhonusq=0.0; sum_nues_div_ms=0.0;
    rhos_nusq_msmq=0.0; rhos_nusqmq_msmq=0.0; sum_nueq_div_mq=0.0;
    rhos_nusqms_msmq=0.0; ne=0.0; Nn=0.0;

    //compute quantities at tn time step
    for (s = 0; s < sl; s++) {
        sum_nues += nust[zk][yj][xi][s];
        sum_nueq += nust[zk][yj][xi][s+7];

        ne += xx[k][j][i].fx[s];
        rhos= xx[k][j][i].fx[s]*ms[s];
        rhoi += rhos;

        Nn += xx[k][j][i].fx[12+s];
        rhon += xx[k][j][i].fx[12+s]*ms[s];

        for (t = 0; t< sm; t++) {
            temp=rhos*nust[zk][yj][xi][21+14*s+t];
            sum_rhonusq += temp;

            rhos_nusq_msmq += temp/(ms[s]+ms[t]);
            rhos_nusqmq_msmq += temp*ms[t]/(ms[s]+ms[t]);
            rhos_nusqms_msmq += temp*ms[s]/(ms[s]+ms[t]);
        }

        sum_nues_div_ms =+ nust[zk][yj][xi][s]/ms[s];
        sum_nueq_div_mq += nust[zk][yj][xi][s+7]/ms[s];
    }

    rhoe=ne*me;
    rhoe_sum_nues=rhoe*sum_nues;
    rhoe_sum_nueq=rhoe*sum_nueq;

    //electron heating terms evaluated at tn step
    Te=xn[k][j][i].fx[11]/(ne*kb);
    Tn=xn[k][j][i].fx[22]/(Nn*kb);
    ele_cooling_rate(xn, ne, Te, Tn, i, j, k, Ce, dCedTe, dCedTn);

/*-------------------------------------------------------------------------
 *------------------- calculate matrix elements
 *-----------------------------------------------------------------------*/
    //for velocity components
    vmat[0][0]=1.0+dt*(rhoe_sum_nues+sum_rhonusq)/rhoi;
    vmat[0][3]=-dt*sum_rhonusq/rhon;

    vmat[1][1]=vmat[0][0];
    vmat[1][4]=vmat[0][3];

    vmat[2][2]=vmat[0][0];
    vmat[2][5]=vmat[0][3];

    vmat[3][0]=-dt*sum_rhonusq/rhoi;
    vmat[3][3]=1.0 + dt*(rhoe_sum_nueq + sum_rhonusq)/rhon;

    vmat[4][1]=vmat[3][0];
    vmat[4][4]=vmat[3][3];

    vmat[5][2]=vmat[3][0];
    vmat[5][5]=vmat[3][3];

    //for pressures
    tmat[0][0]=1.0 + dt2*(me*sum_nues_div_ms + rhos_nusq_msmq/ne);
    tmat[0][1]=-dt2*me*sum_nues_div_ms;
    tmat[0][2]=-dt2*rhos_nusq_msmq/Nn;

    tmat[1][1]=1.0 - two3rd_dt*dCedTe/(kb*ne);
    tmat[1][2]=-two3rd_dt*dCedTn/(kb*Nn);

    tmat[2][0]=-dt2*rhos_nusq_msmq/ne;
    tmat[2][1]=-dt2*me*sum_nueq_div_mq;
    tmat[2][2]=1.0 + dt2*(rhoe*sum_nueq_div_mq + rhos_nusq_msmq)/Nn;

/*-------------------------------------------------------------------------
 *----- right hand side vector; using quantities at tn-dt nut nust at tn
 *-----------------------------------------------------------------------*/
    rhoi=0.0; rhon=0.0; ne=0.0; Nn=0.0;
    for (s = 0; s < sl; s++) {
        ne += xn1[k][j][i].fx[s];
        rhos= xn1[k][j][i].fx[s]*ms[s];
        rhoi += rhos;

        Nn += xn1[k][j][i].fx[12+s];
        rhon += xn1[k][j][i].fx[12+s]*ms[s];
    }

    rhoe=ne*me;
    rhoe_sum_nues=rhoe*sum_nues;
    rhoe_sum_nueq=rhoe*sum_nueq;

    //electron heating terms evaluated at tn-dt step
    double Ce1;
    double Te1=xn1[k][j][i].fx[11]/(ne*kb);
    double Tn1=xn1[k][j][i].fx[22]/(Nn*kb);
    ele_cooling_rate(xn1, ne, Te1, Tn1, i, j, k, Ce1, dCedTe, dCedTn);

    uir=xn1[k][j][i].fx[7]/rhoi; uit=xn1[k][j][i].fx[8]/rhoi; uip=xn1[k][j][i].fx[9]/rhoi;
    uer=uir; uet=uit; uep=uip;
    unr=xn1[k][j][i].fx[19]/rhon; unt=xn1[k][j][i].fx[20]/rhon; unp=xn1[k][j][i].fx[21]/rhon;

    uir_minus_unr=uir - unr;
    uit_minus_unt=uit - unt;
    uip_minus_unp=uip - unp;

    collr=sum_rhonusq*uir_minus_unr;
    collt=sum_rhonusq*uit_minus_unt;
    collp=sum_rhonusq*uip_minus_unp;

    vmat[0][6]= xx[k][j][i].fx[7]+dt*(rhoe_sum_nues*(2.0*uer-uir) - collr);
    vmat[1][6]= xx[k][j][i].fx[8]+dt*(rhoe_sum_nues*(2.0*uet-uit) - collt);
    vmat[2][6]= xx[k][j][i].fx[9]+dt*(rhoe_sum_nues*(2.0*uep-uip) - collp);
    vmat[3][6]= xx[k][j][i].fx[19]+dt*(rhoe_sum_nueq*(2.0*uer-unr) + collr);
    vmat[4][6]= xx[k][j][i].fx[20]+dt*(rhoe_sum_nueq*(2.0*uet-unt) + collt);
    vmat[5][6]= xx[k][j][i].fx[21]+dt*(rhoe_sum_nues*(2.0*uep-unp) + collp);

    coll_Tn_Ti=rhos_nusq_msmq*(xn1[k][j][i].fx[22]/Nn - xn1[k][j][i].fx[10]/ne);
    double veldiff2=two3rd*( uir_minus_unr*uir_minus_unr+uit_minus_unt*uit_minus_unt
                            +uip_minus_unp*uip_minus_unp);

    tmat[0][3]=xx[k][j][i].fx[10]+dt2*( me*sum_nues_div_ms*(xn1[k][j][i].fx[11]-xn1[k][j][i].fx[10])
                                       +coll_Tn_Ti + rhos_nusqmq_msmq*veldiff2);
    tmat[1][3]=xx[k][j][i].fx[11]+two3rd_dt*(Ce+Ce1-dCedTe*Te - dCedTn*Tn);
    tmat[2][3]=xx[k][j][i].fx[22]+dt2*( me*sum_nueq_div_mq*(xn1[k][j][i].fx[11]-rhoe/Nn*xn1[k][j][i].fx[22])
                                       -coll_Tn_Ti + rhos_nusqms_msmq*veldiff2);

/*-------------------------------------------------------------------------
 *------------------- solve linear equations
 *-----------------------------------------------------------------------*/
    //for velocity components
    double delta=(vmat[0][0]*vmat[3][3]-vmat[3][0]*vmat[0][3]);
    bb[0] = (vmat[0][6]*vmat[3][3]-vmat[3][6]*vmat[0][3])/delta;
    bb[3] = (vmat[0][0]*vmat[3][6]-vmat[3][0]*vmat[0][6])/delta;

    delta=(vmat[1][1]*vmat[4][4]-vmat[4][1]*vmat[1][4]);
    bb[1] = (vmat[1][6]*vmat[4][4]-vmat[4][6]*vmat[1][4])/delta;
    bb[4] = (vmat[1][1]*vmat[4][6]-vmat[4][1]*vmat[1][6])/delta;

    delta=(vmat[2][2]*vmat[5][5]-vmat[5][2]*vmat[2][5]);
    bb[2] = (vmat[2][6]*vmat[5][5]-vmat[5][6]*vmat[2][5])/delta;
    bb[5] = (vmat[2][2]*vmat[5][6]-vmat[5][2]*vmat[2][6])/delta;

    //for pressures
    GaussianElimination(tmat, 3, cc);

    for (s = 0; s < 6; s++) {
        if (isnan(bb[s]) || isinf(bb[s])) {
            cout<<"Velocity is Nan or inf at ("<<i<<", "<<j<<", "<<k<<"), s = "<<s<<" in implicit_solver"<<endl;
            exit(-2);
        }
    }

    for (s = 0; s < 3; s++) {
        if (isnan(cc[s]) || isinf(cc[s])) {
            cout<<"Pressure is Nan or inf at ("<<i<<", "<<j<<", "<<k<<"), s = "<<s<<" in implicit_solver"<<endl;
            exit(-2);
        }
        if (cc[s] <= 0.0) {
            cout << "Negative pressure in implicit_solver at ("<<i<<", "<<j<<", "<<k<<"), s = "<<s<<endl;
            exit(-2);
        }
    }
    for (s = 0; s < 3; s++) {
        xx[k][j][i].fx[7+s]=bb[s];
    }
    for (s = 3; s < 6; s++) xx[k][j][i].fx[16+s]=bb[s];
    xx[k][j][i].fx[10]=cc[0];
    xx[k][j][i].fx[11]=cc[1];
    xx[k][j][i].fx[22]=cc[2];
}