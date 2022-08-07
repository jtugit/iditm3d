#------------ fourth order finite difference approximation of -----------------
#---------------------------- pressure gradient -------------------------------
import numpy as np

def evelocity(xx, i, j, k):
    evel=np.zeros(3)

    ne=0.0

    for s in range(7):
        ns=xx[k][j][i][s]

        ne=ne+ns

        evel[0]=evel[0]+ns*xx[k][j][i][7]
        evel[1]=evel[1]+ns*xx[k][j][i][8]
        evel[2]=evel[2]+ns*xx[k][j][i][9]

    evel[0]=evel[0]/ne
    evel[1]=evel[1]/ne
    evel[2]=evel[2]/ne

    return evel
