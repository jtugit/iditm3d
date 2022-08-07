#------------ fourth order finite difference approximation of -----------------
#---------------------------- pressure gradient -------------------------------
import numpy as np

def efield(xx, B0, i, j, k):
    efield=np.zeros(3)

    ver=0.0
    vet=0.0
    vep=0.0
    ne=0.0

    for s in range(7):
        s0=5*l

        ns=xx[k][j][i][s]

        ne=ne+ns

        ver=ver+ns*xx[k][j][i][7]
        vet=vet+ns*xx[k][j][i][8]
        vep=vep+ns*xx[k][j][i][9]

    ver=ver/ne
    vet=vet/ne
    vep=vep/ne

    Br=xx[k][j][i][0]+B0[k][j][i][0]
    Bt=xx[k][j][i][1]+B0[k][j][i][1]
    Bp=xx[k][j][i][2]+B0[k][j][i][2]

    #print(i,ver,vet,vep)
    efield[0]=-(vet*Bp - vep*Bt)
    efield[1]=-(vep*Br - ver*Bp)
    efield[2]=-(ver*Bt - vet*Br)

    return efield
