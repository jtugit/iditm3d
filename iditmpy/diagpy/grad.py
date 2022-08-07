#------------ fourth order finite difference approximation of -----------------
#---------------------------- pressure gradient -------------------------------
import numpy as np

def grad_f(xx, i, j, k, s, a1, a2, a3, cr, r_dth, rsinth_dph):
    grad=np.zeros(3)

    Nr=a1-1
    Nth=a2-1

    if (i == 0):
        grad[0]=cr[7][i]*( 4.0*xx[k][j][i+1].fx[s]-2.08333333*xx[k][j][i].fx[s] \
                          -3.0*xx[k][j][i+2].fx[s]+1.33333333*xx[k][j][i+3].fx[s] \
                          -0.25*xx[k][j][i+4].fx[s])
    elif (i == 1):
        grad[0]=cr[7][i]*( 1.5*xx[k][j][i+1].fx[s]-0.25*xx[k][j][i-1].fx[s] \
                          -0.83333333*xx[k][j][i].fx[s]-0.5*xx[k][j][i+2].fx[s] \
                          +0.08333333*xx[k][j][i+3].fx[s])
    elif (i > 1 and i < Nr):
        grad[0]=cr[7][i]*( 0.08333333*xx[k][j][i-2].fx[s] \
                          -0.66666667*xx[k][j][i-1].fx[s] \
                          +0.66666667*xx[k][j][i+1].fx[s] \
                          -0.08333333*xx[k][j][i+2].fx[s])
    else:
        #Te depends on BC dTe/dr=Fe/lambda_e and involves value at i = Nr+1
        grad[0]=cr[7][i]*( 0.5*xx[k][j][i-2].fx[s]-0.08333333*xx[k][j][i-3].fx[s] \
                          -1.5*xx[k][j][i-1].fx[s]+0.83333333*xx[k][j][i].fx[s] \
                          +0.25*xx[k][j][i+1].fx[s])

    if (j == 0):
        kc=(k+a3/2) % a3
        grad[1]=( 0.08333333*xx[kc][1][i].fx[s]-0.66666667*xx[kc][0][i].fx[s] \
                 +0.66666667*xx[k][j+1][i].fx[s]-0.08333333*xx[k][j+2][i].fx[s])

    elif (j == 1):
        kc=(k+a3/2) % a3
        grad[1]=( 0.08333333*xx[kc][0][i].fx[s]-0.66666667*xx[k][j-1][i].fx[s] \
                 +0.66666667*xx[k][j+1][i].fx[s]-0.08333333*xx[k][j+2][i].fx[s])       

    elif (j > 1 and j < Nth-1):
        grad[1]=( 0.08333333*xx[k][j-2][i].fx[s]-0.66666667*xx[k][j-1][i].fx[s] \
                 +0.66666667*xx[k][j+1][i].fx[s]-0.08333333*xx[k][j+2][i].fx[s])        

    elif (j == Nth-1):
        kc=(k+a3/2) % a3
        grad[1]=( 0.08333333*xx[k][j-2][i].fx[s]-0.66666667*xx[k][j-1][i].fx[s] \
                 +0.66666667*xx[k][j+1][i].fx[s]-0.08333333*xx[kc][Nth][i].fx[s])        

    else:
        kc=(k+a3/2) % a3
        grad[1]=( 0.08333333*xx[k][j-2][i].fx[s]-0.66666667*xx[k][j-1][i].fx[s] \
                 +0.66666667*xx[kc][Nth][i].fx[s]-0.08333333*xx[kc][Nth-1][i].fx[s])

    grad[1]=grad[1]/r_dth[i]

    grad[2]=( 0.08333333*xx[k-2][j][i].fx[s]-0.66666667*xx[k-1][j][i].fx[s] \
             +0.66666667*xx[k+1][j][i].fx[s]-0.08333333*xx[k+2][j][i].fx[s]) \
            /rsinth_dph[j][i]       

    return grad
