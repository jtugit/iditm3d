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
    double tmat[3][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
    double bb[6], cc[3];
    double rhos, rhoi, rhon, sum_nues, sum_nueq, sum_rhonusq, sum_nues_div_ms;
    double rhos_nusq_msmq, rhos_nusqmq_msmq, rhos_nusqms_msmq, temp;
    double ne, rhoe, uir_minus_unr, uit_minus_unt, uip_minus_unp, Nn; //ne = ni
    double Len_part, Te, Tn, Ti, rhoe_sum_nues, rhoe_sum_nueq;
    double collr, collt, collp, uiminusun_sq;
    double uer, uet, uep, uir, uit, uip, unr, unt, unp;
    double Qefric, Qifric, Qnfric, Te12;
    const double two3rd=2.0/3.0, two3rd_dt = 2.0/3.0*dt;

    rhoi=0.0; rhon=0.0;
    sum_nues=0.0; sum_nueq=0.0; sum_rhonusq=0.0; sum_nues_div_ms=0.0;
    rhos_nusq_msmq=0.0; rhos_nusqmq_msmq=0.0;
    rhos_nusqms_msmq=0.0; ne=0.0; Nn=0.0;

    const double Ei[3]={228.0, 326.0, 98.0};
    const double Ai[3]={7.833e-6, 9.466e-6, 1.037e-8};
    const double Bi[3]={1.021, 0.8458, 1.633};
    const double Ci[3]={1.009, 0.9444, 1.466};

    //compute quantities at tn+dt ne, rhoi, rhon, nust at tn time step
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
    }

    uer=xx[k][j][i].fx[7]/rhoi; uet=xx[k][j][i].fx[8]/rhoi; uep=xx[k][j][i].fx[9]/rhoi;
    uir=uer; uit=uet; uip=uep;
    unr=xx[k][j][i].fx[19]/rhon; unt=xx[k][j][i].fx[20]/rhon; unp=xx[k][j][i].fx[21]/rhon;

    rhoe=ne*me;
    rhoe_sum_nues=rhoe*sum_nues;
    rhoe_sum_nueq=rhoe*sum_nueq;

    double rhoe_sum_nues_uer = rhoe_sum_nues*uer;
    double rhoe_sum_nues_uet = rhoe_sum_nues*uet;
    double rhoe_sum_nues_uep = rhoe_sum_nues*uep;
    double rhoe_sum_nueq_uer = rhoe_sum_nueq*uer;
    double rhoe_sum_nueq_uet = rhoe_sum_nueq*uet;
    double rhoe_sum_nueq_uep = rhoe_sum_nueq*uep;

    Qefric=two3rd*( rhoe_sum_nues*((uer-uir)*(uer-uir)+(uet-uit)*(uet-uit)+(uep-uip)*(uep-uip))
                   +rhoe_sum_nueq*((uer-unr)*(uer-unr)+(uet-unt)*(uet-unt)+(uep-unp)*(uep-unp)));

    uiminusun_sq=(uir-unr)*(uir-unr)+(uit-unt)*(uit-unt)+(uip-unp)*(uip-unp);
    Qifric=two3rd*rhos_nusqmq_msmq*uiminusun_sq;
    Qnfric=two3rd*rhos_nusqms_msmq*uiminusun_sq;

    //electron heating terms evaluated at tn step
    Te=xn[k][j][i].fx[11]/(ne*kb);
    Tn=xn[k][j][i].fx[22]/(Nn*kb);

    //part of thermal energy change between electrons and neutrals are treated explicitly
    Len_part = 0.0; // ele_cooling_rate(xn, ne, Te, Tn, i, j, k);

    //part of fine structure of O treated implicitly: at tn time step
    double Z=5.0+3.0*exp(-228.0/Tn)+exp(-326.0/Tn), Dxi[3], Exi[3];
    Dxi[0]=exp(-228.0/Tn); Dxi[1]=exp(-326.0/Tn); Dxi[2]=exp(-326.0/Tn);
    Exi[0]=exp(-228.0/Te); Exi[1]=exp(-326.0/Te); Exi[2]=exp(-(98.0/Te+228.0/Tn));

    double Le_coeff=0.0;
    for (int m = 0; m < 3; m++) {
        Le_coeff += Ai[m]*Ci[m]*pow(Te, Bi[m]-0.5)*((1.0+Bi[m])*Dxi[m]+(Ei[m]/Te+1.0+Bi[m])*Exi[m]);
    }
    Le_coeff = 5.099739e-14*xn[k][j][i].fx[12]/Z*Le_coeff;

    Te12=sqrt(Te);

    //add contribution from other cooling processes
    Le_coeff = ( Le_coeff+(2.9e-14/Te12 + 1.77e-19*(1.0-1.21e-4*Te)*Te)*xx[k][j][i].fx[16]
                +(6.9e-14/Te12+1.21e-18*(1.0+3.6e-2*Te12)*Te12)*xx[k][j][i].fx[15]
                +7.9e-19*xx[k][j][i].fx[12]*(1.0+5.7e-4*Te)*Te12
                +9.63e-16*xx[k][j][i].fx[13]*(1.0-1.35e-4*Te)*Te12
                +2.46e-17*xx[k][j][i].fx[14]*Te12)*ne*e;

/*-------------------------------------------------------------------------
 *------------------- calculate matrix elements
 *-----------------------------------------------------------------------*/
    //for velocity components
    double vmatr[2][3], vmatt[2][3], vmatp[2][3];
    vmatr[0][0]=1.0+dt*(rhoe_sum_nues+sum_rhonusq)/rhoi;      //a1
    vmatr[0][1]=-dt*sum_rhonusq/rhon;                         //b1
    vmatr[1][0]=-dt*sum_rhonusq/rhoi;                         //a2
    vmatr[1][1]=1.0 + dt*(rhoe_sum_nueq + sum_rhonusq)/rhon;  //b2

    vmatt[0][0]=vmatr[0][0];
    vmatt[0][1]=vmatr[0][1];
    vmatt[1][0]=vmatr[1][0];
    vmatt[1][1]=vmatr[1][1];

    vmatp[0][0]=vmatr[0][0];
    vmatp[0][1]=vmatr[0][1];
    vmatp[1][0]=vmatr[1][0];
    vmatp[1][1]=vmatr[1][1];

    //for pressures
    tmat[0][0]=1.0 + dt2*(me*sum_nues_div_ms + rhos_nusq_msmq/ne);
    tmat[0][1]=-dt2*me*sum_nues_div_ms;
    tmat[0][2]=-dt2*rhos_nusq_msmq/Nn;

    tmat[1][0]=-dt2*me*sum_nues_div_ms;
    tmat[1][1]=1.0+dt2*(me*sum_nues_div_ms + Le_coeff/(3.0*ne*kb));
    tmat[1][2]=-two3rd_dt*Le_coeff/(Nn*kb);

    tmat[2][0]=-dt2*rhos_nusq_msmq/ne;
    tmat[2][1]=-two3rd_dt*Le_coeff/(ne*kb);
    tmat[2][2]=1.0 + dt2*(rhos_nusq_msmq/Nn + Le_coeff/(3.0*Nn*kb));

/*-------------------------------------------------------------------------
 *----- right hand side vector; using quantities at tn-dt but nust at tn
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

    uir=xn1[k][j][i].fx[7]/rhoi; uit=xn1[k][j][i].fx[8]/rhoi; uip=xn1[k][j][i].fx[9]/rhoi;
    uer=uir; uet=uit; uep=uip;
    unr=xn1[k][j][i].fx[19]/rhon; unt=xn1[k][j][i].fx[20]/rhon; unp=xn1[k][j][i].fx[21]/rhon;

    uir_minus_unr=uir - unr;
    uit_minus_unt=uit - unt;
    uip_minus_unp=uip - unp;

    collr=sum_rhonusq*uir_minus_unr;
    collt=sum_rhonusq*uit_minus_unt;
    collp=sum_rhonusq*uip_minus_unp;

    vmatr[0][2]= xx[k][j][i].fx[7]+dt*(rhoe_sum_nues*(uer-uir)+rhoe_sum_nues_uer - collr);  //c1
    vmatr[1][2]=xx[k][j][i].fx[19]+dt*(rhoe_sum_nueq*(uer-unr)+rhoe_sum_nueq_uer + collr);  //c2
    vmatt[0][2]= xx[k][j][i].fx[8]+dt*(rhoe_sum_nues*(uet-uit)+rhoe_sum_nues_uet - collt);
    vmatt[1][2]=xx[k][j][i].fx[20]+dt*(rhoe_sum_nueq*(uet-unt)+rhoe_sum_nueq_uet + collt);
    vmatp[0][2]= xx[k][j][i].fx[9]+dt*(rhoe_sum_nues*(uep-uip)+rhoe_sum_nues_uep - collp);
    vmatp[1][2]=xx[k][j][i].fx[21]+dt*(rhoe_sum_nueq*(uep-unp)+rhoe_sum_nueq_uep + collp);

    //electron heating terms evaluated at tn-dt step
    Te=xn1[k][j][i].fx[11]/(ne*kb);
    Ti=xn1[k][j][i].fx[10]/(ne*kb);
    Tn=xn1[k][j][i].fx[22]/(Nn*kb);
    Te12=sqrt(Te);

    //part of fine structure of O treated implicitly: at tn-dt time step
    Z=5.0+3.0*exp(-228.0/Tn)+exp(-326.0/Tn);
    Dxi[0]=exp(-228.0/Tn); Dxi[1]=exp(-326.0/Tn); Dxi[2]=exp(-326.0/Tn);
    Exi[0]=exp(-228.0/Te); Exi[1]=exp(-326.0/Te); Exi[2]=exp(-(98.0/Te+228.0/Tn));

    Le_coeff=0.0;
    for (int m = 0; m < 3; m++) {
        Le_coeff += Ai[m]*Ci[m]*pow(Te, Bi[m]-0.5)*((1.0+Bi[m])*Dxi[m]+(Ei[m]/Te+1.0+Bi[m])*Exi[m]);
    }
    Le_coeff = 5.099739e-14*xn1[k][j][i].fx[12]/Z*Le_coeff;

    Le_coeff=( Le_coeff+ (2.9e-14/Te12 + 1.77e-19*(1.0-1.21e-4*Te)*Te)*xx[k][j][i].fx[16]
              +(6.9e-14/Te12+1.21e-18*(1.0+3.6e-2*Te12)*Te12)*xx[k][j][i].fx[15]
              +7.9e-19*xx[k][j][i].fx[12]*(1.0+5.7e-4*Te)*Te12
              +9.63e-16*xx[k][j][i].fx[13]*(1.0-1.35e-4*Te)*Te12
              +2.46e-17*xx[k][j][i].fx[14]*Te12)*ne*e;

    double coll_Te_Ti=me*sum_nues_div_ms*kb*(Te-Ti);
    double coll_Tn_Ti=rhos_nusq_msmq*kb*(Tn-Ti);
    double coll_Tn_Te=two3rd_dt*Le_coeff*(Tn-Te);

    tmat[0][3]=xx[k][j][i].fx[10]+dt2*(Qifric+coll_Te_Ti+coll_Tn_Ti);
    tmat[1][3]=xx[k][j][i].fx[11]+dt2*(Qefric-coll_Te_Ti+two3rd*Len_part)+coll_Tn_Te;
    tmat[2][3]=xx[k][j][i].fx[22]+dt2*(Qnfric-coll_Tn_Ti-two3rd*Len_part)-coll_Tn_Te;

/*-------------------------------------------------------------------------
 *------------------- solve linear equations
 *-----------------------------------------------------------------------*/
    //for velocity components
    //vmat[0][0] - a1, vmat[0][1] - b1, vmat[1][0] - a2, vmat[1][1] - b2
    double delta=(vmatr[0][0]*vmatr[1][1]-vmatr[1][0]*vmatr[0][1]);
    bb[0] = (vmatr[0][2]*vmatr[1][1]-vmatr[1][2]*vmatr[0][1])/delta;
    bb[1] = (vmatr[0][0]*vmatr[1][2]-vmatr[1][0]*vmatr[0][2])/delta;

    delta=(vmatt[0][0]*vmatt[1][1]-vmatt[1][0]*vmatt[0][1]);
    bb[2] = (vmatt[0][2]*vmatt[1][1]-vmatt[1][2]*vmatt[0][1])/delta;
    bb[3] = (vmatt[0][0]*vmatt[1][2]-vmatt[1][0]*vmatt[0][2])/delta;

    delta=(vmatp[0][0]*vmatp[1][1]-vmatp[1][0]*vmatp[0][1]);
    bb[4] = (vmatp[0][2]*vmatp[1][1]-vmatp[1][2]*vmatp[0][1])/delta;
    bb[5] = (vmatp[0][0]*vmatp[1][2]-vmatp[1][0]*vmatp[0][2])/delta;

    //for pressures
    GaussianElimination(tmat, 3, cc);

    //if (i==20 && j==13 && k==68) cout <<sqrt(bb[0]*bb[0]+bb[2]*bb[2]+bb[4]*bb[4]) <<" "<<cc[1]<<endl;
    for (s = 0; s < 6; s++) {
        if (isnan(bb[s]) || isinf(bb[s])) {
            cout<<"Velocity is Nan or inf at ("<<i<<", "<<j<<", "<<k<<"), s = "<<s
                <<" in implicit_solver"<<endl;
            exit(-2);
        }
    }

    for (s = 0; s < 3; s++) {
        if (isnan(cc[s]) || isinf(cc[s])) {
            cout<<"Pressure is Nan or inf at ("<<i<<", "<<j<<", "<<k<<"), s = "<<s
                <<" in implicit_solver"<<endl;
            exit(-2);
        }
        if (cc[s] <= 0.0) {
            cout << "Negative pressure in implicit_solver at ("<<i<<", "<<j<<", "<<k<<"), s = "<<s<<endl;
            exit(-2);
        }
    }
    xx[k][j][i].fx[7]=bb[0];
    xx[k][j][i].fx[8]=bb[2];
    xx[k][j][i].fx[9]=bb[4];
    xx[k][j][i].fx[19]=bb[1];
    xx[k][j][i].fx[20]=bb[3];
    xx[k][j][i].fx[21]=bb[5];
    xx[k][j][i].fx[10]=cc[0];
    xx[k][j][i].fx[11]=cc[1];
    xx[k][j][i].fx[22]=cc[2];
}