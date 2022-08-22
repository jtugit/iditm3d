import math

def fast_speed(xx, B0, i, j, k, sl):
    ms=[16.0,1.0,4.0,32.0,28,30.0,14.0]
    me=9.1e-31
    mp=1.67216219e-27
    kb=1.380649e-23
    mu0=4.0e-7*math.pi

    ne=0.0
    rho=0.0
    p=0.0
    for s in range(sl):
        ns = xx[k][j][i][s]

        ne = ne + ns
        rho = rho + ns*ms[s]

    #p in J cm^{-3}
    p = ne*kb*(xx[k][j][i][10] + xx[k][j][i][11])

    #rho in kg cm^{-3}
    rho=rho*mp

    #sound speed square in (m/s)^2
    gama=5.0/3.0
    Cs2=gama*p/rho

    #Alfven speed square in (m/s)^2
    BB2=B0[k][j][i][0]**2+B0[k][j][i][1]**2+B0[k][j][i][2]**2
    BB=math.sqrt(BB2)
    Va2=BB2/(mu0*rho*1.0e6)

    #cosin of angle between magnetic field and coordinate r, theta and phi
    costh_r=B0[k][j][i][0]/BB
    costh_t=B0[k][j][i][1]/BB
    costh_f=B0[k][j][i][2]/BB
    Va2Cs2=Va2+Cs2

    #fast wave speed along each of three individual coordinate in m/s
    vf_r=math.sqrt(0.5*(Va2Cs2+math.sqrt(Va2Cs2**2-4.0*Va2*Cs2*costh_r**2)))
    vf_t=math.sqrt(0.5*(Va2Cs2+math.sqrt(Va2Cs2**2-4.0*Va2*Cs2*costh_t**2)))
    vf_f=math.sqrt(0.5*(Va2Cs2+math.sqrt(Va2Cs2**2-4.0*Va2*Cs2*costh_f**2)))
    Cs=math.sqrt(Cs2)

    return vf_r, vf_t, vf_f, Cs