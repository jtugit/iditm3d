#------------ fourth order finite difference approximation of -----------------
#---------------------------- pressure gradient -------------------------------
import numpy as np

def evelocity(xx, i, j, k):
    evel=np.zeros(3)

    ne=0.0

    for l in range(6):
        s0=5*l

        ns=xx[k][j][i][4+s0]

        ne=ne+ns

        evel[0]=evel[0]+ns*xx[k][j][i][5+s0]
        evel[1]=evel[1]+ns*xx[k][j][i][6+s0]
        evel[2]=evel[2]+ns*xx[k][j][i][7+s0]

    evel[0]=evel[0]/ne
    evel[1]=evel[1]/ne
    evel[2]=evel[2]/ne

    return evel
