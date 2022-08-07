def divB(xx, i, j, k, alt, cr, rdth, cos_rsinth, rsinth_dph, a1, a2, a3):
    Nr=a1-1
    Nth=a2-1
    Np=a3-1

    rr=alt[i]*1e3+6371.2e3
    r0=1.0e6

    if i==0:
        div_r=div_r= 1.5*xx[k][j][i+1][0]-0.25*xx[k][j][i-1][0] \
              -0.83333333*xx[k][j][i][0]-0.5*xx[k][j][i+2][0] \
              +0.08333333*xx[k][j][i+3][0]
    elif i > 0 and i < Nr-2:
        div_r= 0.08333333*xx[k][j][i-2][0]-0.66666667*xx[k][j][i-1][0] \
              +0.66666667*xx[k][j][i+1][0]-0.08333333*xx[k][j][i+2][0]
    elif i==Nr-1:
        div_r= 0.5*xx[k][j][i-2][0]-1.0/12.0*xx[k][j][i-3][0] \
              -1.5*xx[k][j][i-1][0]+10.0/12.0*xx[k][j][i][0] \
              +0.25*xx[k][j][i+1][0]
    else:
        div_r= 0.25*xx[k][j][i-4][0]-4.0/3.0*xx[k][j][i-3][0] \
              +3.0*xx[k][j][i-2][0]-4.0*xx[k][j][i-1][0] \
              +25.0/12.0*xx[k][j][i][0]

    if j == 0:
        kc=(k+int(a3/2)) % a3
        div_th=( 0.66666667*xx[kc][0][i][1]-0.08333333*xx[kc][1][i][1]
                +0.66666667*xx[k][j+1][i][1]-0.08333333*xx[k][j+2][i][1])
    elif j == 1:
        kc=(k+int(a3/2)) % a3
        div_th=( 0.66666667*xx[k][j+1][i][1]-0.08333333*xx[kc][0][i][1]
                -0.66666667*xx[k][j-1][i][1]-0.08333333*xx[k][j+2][i][1])
    elif j > 1 and j < Nth-1:
        div_th=( 0.08333333*xx[k][j-2][i][1]-0.66666667*xx[k][j-1][i][1]
                +0.66666667*xx[k][j+1][i][1]-0.08333333*xx[k][j+2][i][1])
    elif j == Nth-1:
        kc=(k+int(a3/2)) % a3
        div_th=( 0.08333333*xx[k][j-2][i][1]-0.66666667*xx[k][j-1][i][1]
                +0.66666667*xx[k][j+1][i][1]+0.08333333*xx[kc][Nth][i][1])
    else:
        kc=(k+int(a3/2)) % a3
        div_th=( 0.08333333*xx[k][j-2][i][1]-0.66666667*xx[k][j-1][i][1]
                -0.66666667*xx[kc][Nth][i][1]+0.08333333*xx[kc][Nth-1][i][1])

    if k==0:
        km2=Np-1
        km1=Np
        kp1=1
        kp2=2
    elif k==1:
        km2=Np
        km1=0
        kp1=2
        kp2=3
    elif k>1 and k<Np-1:
        km2=k-2
        km1=k-1
        kp1=k+1
        kp2=k+2
    elif k==Np-1:
        km2=k-2
        km1=k-1
        kp1=Np
        kp2=0
    else:
        km2=k-2
        km1=k-1
        kp1=0
        kp2=1

    div_ph=( 0.08333333*xx[km2][j][i][2]-0.66666667*xx[km1][j][i][2]
            +0.66666667*xx[kp1][j][i][2]-0.08333333*xx[kp2][j][i][2])

    div= 2.0*xx[k][j][i][0]/rr+(cr[i]*div_r \
        +xx[k][j][i][1]*cos_rsinth[j][i]+div_th/rdth[i] \
        +div_ph/rsinth_dph[j][i])/r0

    return div