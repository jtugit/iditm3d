import math
import numpy as np

def paramaters(iditm_arr, i, j, k):
    ms=[16.0,32.0,28.0,1.0,4.0,30.0,14.0]
    me=9.10938e-31/1.67262e-27

# coefficients of all collision frequencies
    coe=[0.0]*15
    coe[0]=54.5e-6               # e - ions
    coe[1]=8.90e-17              # e - O
    coe[2]=1.82e-16              # e - O2
    coe[3]=2.33e-17              # e - N2
    coe[4]=4.50e-15              # e - H
    coe[5]=4.60e-16              # e - He
    coe[6]=1.0e-6                # e - NO
    coe[7]=1.0e-6                # e - N

     #for these used in Te equation, multiplied by by me/(me+ms)
    coe[8]= me/ms[0]*8.90e-17    # e - O
    coe[9]= me/ms[1]*1.82e-16    # e - O2
    coe[10]=me/ms[2]*2.33e-17    # e - N2
    coe[11]=me/ms[3]*4.50e-15    # e - H
    coe[12]=me/ms[4]*4.60e-16    # e - He
    coe[13]=me/ms[5]*1.0e-6      # e - NO
    coe[14]=me/ms[6]*1.0e-6      # e - N

    coiO=[0.0]*12
    coiO[0]=2.60e-7        # O+ - O2+
    coiO[1]=2.50e-7        # O+ - N2+
    coiO[2]=7.70e-8        # O+ - H+
    coiO[3]=1.40e-7        # O+ - H+e
    coiO[4]=2.60e-7        # O+ - NO+
    coiO[5]=1.5*3.67e-17         # O+ - O
    coiO[6]=6.64e-16             # O+ - O2
    coiO[7]=6.82e-16             # O+ - N2
    coiO[8]=6.61e-17             # O+ - H
    coiO[9]=1.32e-16             # O+ - He
    coiO[10]=2.69e-15            # O+ - NO
    coiO[11]=4.62e-16            # O+ - N

    coiO2=[0.0]*12
    coiO2[0]=1.30e-7           # O2+ - O+
    coiO2[1]=1.50e-7           # O2+ - N2+
    coiO2[2]=3.90e-8           # O2+ - H+
    coiO2[3]=7.50e-8           # O2+ - He+
    coiO2[4]=1.60e-7           # O2+ - NO+
    coiO2[5]=2.31e-16          # O2+ - O
    coiO2[6]=2.59e-17          # O2+ - O2
    coiO2[7]=4.13e-16          # O2+ - N2
    coiO2[8]=6.50e-17          # O2+ - H
    coiO2[9]=7.00e-17          # O2+ - He
    coiO2[10]=2.69e-15         # O2+ - NO
    coiO2[11]=2.64e-16         # O2+ - N

    coiN2=[0.0]*12
    coiN2[0]=1.50e-7         # N2+ - O+
    coiN2[1]=1.80e-7         # N2+ - O2+
    coiN2[2]=4.50e-8         # N2+ - H+
    coiN2[3]=8.50e-8         # N2+ - He+
    coiN2[4]=1.70e-7         # N2+ - NO+
    coiN2[5]=2.58e-16           # N2+ - O
    coiN2[6]=4.49e-16           # N2+ - O2
    coiN2[7]=5.14e-17           # N2+ - N2
    coiN2[8]=7.40e-17           # N2+ - H
    coiN2[9]=7.90e-17           # N2+ - He
    coiN2[10]=2.69e-15          # N2+ - NO
    coiN2[11]=2.95e-16          # N2+ - N
 
    coiH=[0.0]*12
    coiH[0]=1.23e-6       # H+ - O+
    coiH[1]=1.25e-6       # H+ - O2+
    coiH[2]=1.25e-6       # H+ - N2+
    coiH[3]=1.14e-6       # H+ - He+
    coiH[4]=1.25e-6       # H+ - NO+
    coiH[5]=6.61e-17            # H+ - O
    coiH[6]=3.20e-15            # H+ - O2
    coiH[7]=3.36e-15            # H+ - N2
    coiH[8]=2.65e-16            # H+ - H
    coiH[9]=10.6e-16            # H+ - He
    coiH[10]=2.69e-15           # H+ - NO
    coiH[11]=26.1e-16           # H+ - N

    coiHe=[0.0]*12
    coiHe[0]=0.57e-6       # He+ - O+
    coiHe[1]=0.60e-6       # He+ - O2+
    coiHe[2]=0.59e-6       # He+ - N2+
    coiHe[3]=0.28e-6       # He+ - H+
    coiHe[4]=0.60e-6       # He+ - NO+
    coiHe[5]=10.1e-16         # He+ - O
    coiHe[6]=15.3e-16         # He+ - O2
    coiHe[7]=16.0e-16         # He+ - N2
    coiHe[8]=4.71e-16         # He+ - H
    coiHe[9]=8.73e-17         # He+ - He
    coiHe[10]=2.69e-15        # He+ - NO
    coiHe[11]=11.9e-16        # He+ - N

    coiNO=[0.0]*12 
    coiNO[0]=0.14e-6         # NO+ - O+
    coiNO[1]=1.70e-7         # NO+ - O2+
    coiNO[2]=1.60e-7         # NO+ - N2+
    coiNO[3]=4.20e-8         # NO+ - H+
    coiNO[4]=8.00e-8         # NO+ - He+
    coiNO[5]=2.44e-16          # NO+ - O
    coiNO[6]=4.27e-16          # NO+ - O2
    coiNO[7]=4.34e-16          # NO+ - N2
    coiNO[8]=6.90e-17          # NO+ - H
    coiNO[9]=7.40e-17          # NO+ - He
    coiNO[10]=2.69e-25         # NO+ - NO
    coiNO[11]=2.79e-16         # NO+ - N
 
    con=[0.0]*10
    con[0]=2.26216e-17      # O - O2
    con[1]=7.60616e-18      # O - N2
    con[2]=5.15410e-18      # O2 - N2
    con[3]=9.05133e-18      # H - O
    con[4]=5.43880e-18      # H - O2
    con[5]=6.05369e-18      # H - N2
    con[6]=2.33451e-17      # H - He
    con[7]=1.49978e-17      # He - O
    con[8]=8.63023e-18      # He - O2
    con[9]=1.00277e-17      # He - N2

#--------------------------------------------------*/
#---------------- collision frequencies -----------*/
#--------------------------------------------------*/
    Te=iditm_arr[k][j][i][3]
    Te12=math.sqrt(Te)
    Te32=Te*Te12
    cee=coe[0]/Te32

    nuet=np.zeros(7)
    nust=np.zeros(81)

    temp=np.zeros(4)
    Ti=np.zeros(6)
    ni=np.zeros(6)
    nn=np.zeros(7)

    for l in range(7):
        nn[l]=iditm_arr[k][j][i][34+l]

    Tn = iditm_arr[k][j][i][48]

    ne=0.0
    for l in range(6):
        s0=5*l

        #electron density
        ni[l]=iditm_arr[k][j][i][4+s0]
        ne = ne+ni[l]

        #ion temperature
        Ti[l]=iditm_arr[k][j][i][8+s0]

        #electron Coulomb collision frequencies with O+, O2+, N2+, H+, He+, NO+
        nuet[l]=cee*ni[l]

    #electron - neutral collision frequencies
    tem=[0.0]*4
    tem[0]=nn[0]*(1.0+5.7e-4*Te)
    tem[1]=nn[1]*(1.0+3.6e-2*Te12)
    tem[2]=nn[2]*(1.0-1.21e-4*Te)*Te
    tem[3]=nn[3]*(1.0-1.35e-4*Te)

    nust[76]=coe[1]*tem[0]*Te12
    nust[77]=coe[2]*tem[1]*Te12
    nust[78]=coe[3]*tem[2]
    nust[79]=coe[4]*tem[3]*Te12
    nust[80]=coe[5]*nn[4]*Te12

    #total electron - neutral collision frequency nu_{e,n}
    nuet[6]= nust[76] +nust[77]+nust[78]+nust[79]+nust[80]

    #---- O+ collision frequencies --------------------------------------------------*/
    # O+ Coulomb collision frequencies
    nust[0]=ne*me/(ni[0]*ms[0])*nuet[0]           # O+ - e
    nust[1]=coiO[0]*ni[1]/(Ti[1]*math.sqrt(Ti[1]))     # O+ - O2+
    nust[2]=coiO[1]*ni[2]/(Ti[2]*math.sqrt(Ti[2]))     # O+ - N2+
    nust[3]=coiO[2]*ni[3]/(Ti[3]*math.sqrt(Ti[3]))     # O+ - H+
    nust[4]=coiO[3]*ni[4]/(Ti[4]*math.sqrt(Ti[4]))     # O+ - He+
    nust[5]=coiO[4]*ni[5]/(Ti[5]*math.sqrt(Ti[5]))     # O+ - NO+

    # O+ - neutral collision frequencies
    Tr=0.5*(Ti[0]+Tn)
    tem[0]=nn[0]*math.sqrt(Tr)*(1.0-0.064*math.log10(Tr))**2
    tem[1]=nn[3]*math.sqrt(Ti[0])*(1.0-0.047*math.log10(Ti[0]))**2
    nust[6]= coiO[5]*tem[0]                        # O+ - O
    nust[7]= coiO[6]*nn[1]                         # O+ - O2
    nust[8]= coiO[7]*nn[2]                         # O+ - N2
    nust[9]= coiO[8]*tem[1]                        # O+ - H
    nust[10]= coiO[9]*nn[4]                        # O+ - He

#---- O2+ collision frequencies --------------------------------------------------------
    # O2+ Coulomb collision frequencies
    nust[11]=ne*me/(ni[1]*ms[1])*nuet[1]           # O2+ - e
    nust[12]=coiO2[0]*ni[0]/(Ti[0]*math.sqrt(Ti[0]))    # O2+ - O+
    nust[13]=coiO2[1]*ni[2]/(Ti[2]*math.sqrt(Ti[2]))    # O2+ - N2+
    nust[14]=coiO2[2]*ni[3]/(Ti[3]*math.sqrt(Ti[3]))    # O2+ - H+
    nust[15]=coiO2[3]*ni[4]/(Ti[4]*math.sqrt(Ti[4]))    # O2+ - He+
    nust[16]=coiO2[4]*ni[5]/(Ti[5]*math.sqrt(Ti[5]))    # O2+ - NO+

    # O2+ - neutral collision frequencies
    Tr=0.5*(Ti[1]+Tn)
    tem[0]=nn[1]*math.sqrt(Tr)*(1.0-0.073*math.log10(Tr))**2
    nust[17]= coiO2[5]*nn[0]               # O2+ - O
    nust[18]= coiO2[6]*tem[0]              # O2+ - O2
    nust[19]= coiO2[7]*nn[2]               # O2+ - N2
    nust[20]= coiO2[8]*nn[3]               # O2+ - H
    nust[21]= coiO2[9]*nn[4]               # O2+ - He

# ---- N2+ collision frequencies --------------------------------------------------------
    # N2+ Coulomb collision frequencies
    nust[22]=ne*me/(ni[2]*ms[2])*nuet[2]             # N2+ - e
    nust[23]=coiN2[0]*ni[0]/(Ti[0]*math.sqrt(Ti[0]))      # N2+ - O+
    nust[24]=coiN2[1]*ni[1]/(Ti[1]*math.sqrt(Ti[1]))      # N2+ - O2+
    nust[25]=coiN2[2]*ni[3]/(Ti[3]*math.sqrt(Ti[3]))      # N2+ - H+
    nust[26]=coiN2[3]*ni[4]/(Ti[4]*math.sqrt(Ti[4]))      # N2+ - He+
    nust[27]=coiN2[4]*ni[5]/(Ti[5]*math.sqrt(Ti[5]))      # N2+ - NO+

    # N2+ - neutral collision frequencies
    Tr=0.5*(Ti[2]+Tn)
    tem[0]=nn[2]*math.sqrt(Tr)*(1.0-0.069*math.log10(Tr))**2
    nust[28]= coiN2[5]*nn[0]                 # N2+ - O
    nust[29]= coiN2[6]*nn[1]                 # N2+ - O2
    nust[30]= coiN2[7]*tem[0]                # N2+ - N2
    nust[31]= coiN2[8]*nn[3]                 # N2+ - H
    nust[32]= coiN2[9]*nn[4]                 # N2+ - He

# ---- H+ collision frequencies ---------------------------------------------------------
    # H+ Coulomb collision frequencies
    nust[33]=ne*me/(ni[3]*ms[3])*nuet[3]           # H+ - e
    nust[34]=coiH[0]*ni[0]/(Ti[0]*math.sqrt(Ti[0]))     # H+ - O+
    nust[35]=coiH[1]*ni[1]/(Ti[1]*math.sqrt(Ti[1]))     # H+ - O2+
    nust[36]=coiH[2]*ni[2]/(Ti[2]*math.sqrt(Ti[2]))     # H+ - N2+
    nust[37]=coiH[3]*ni[4]/(Ti[4]*math.sqrt(Ti[4]))     # H+ - He+
    nust[38]=coiH[4]*ni[5]/(Ti[5]*math.sqrt(Ti[5]))     # H+ - NO+

    # H+ - neutral collision frequencies
    Tr=0.5*(Ti[3]+Tn)
    tem[0]=nn[0]*math.sqrt(Ti[3])*(1.0-0.047*math.log10(Ti[3]))**2
    tem[1]=nn[3]*math.sqrt(Tr)*(1.0-0.083*math.log10(Tr))**2
    nust[39]= coiH[5]*tem[0]                  # H+ - O
    nust[40]= coiH[6]*nn[1]                   # H+ - O2
    nust[41]= coiH[7]*nn[2]                   # H+ - N2
    nust[42]= coiH[8]*tem[1]                  # H+ - H
    nust[43]= coiH[9]*nn[4]                   # H+ - He

# ---- He+ collision frequencies ---------------------------------------------------------
    # He+ Coulomb collision frequencies
    nust[44]=ne*me/(ni[4]*ms[4])*nuet[4]             # He+ - e
    nust[45]=coiHe[0]*ni[0]/(Ti[0]*math.sqrt(Ti[0]))      # He+ - O+
    nust[46]=coiHe[1]*ni[1]/(Ti[1]*math.sqrt(Ti[1]))      # He+ - O2+
    nust[47]=coiHe[2]*ni[2]/(Ti[2]*math.sqrt(Ti[2]))      # He+ - N2+
    nust[48]=coiHe[3]*ni[3]/(Ti[3]*math.sqrt(Ti[3]))      # He+ - H+
    nust[49]=coiHe[4]*ni[5]/(Ti[5]*math.sqrt(Ti[4]))      # He+ - NO+

    # He+ - neutral collision frequencies
    Tr=0.5*(Ti[4]+Tn)
    tem[0]=nn[4]*math.sqrt(Tr)*(1.0-0.083*math.log10(Tr))**2
    nust[50]= coiHe[5]*nn[0]                   # He+ - O
    nust[51]= coiHe[6]*nn[1]                   # He+ - O2
    nust[52]= coiHe[7]*nn[2]                   # He+ - N2
    nust[53]= coiHe[8]*nn[3]                   # He+ - H
    nust[54]= coiHe[9]*tem[0]                  # He+ - He

# ---- NO+ collision frequencies -------------------------------------------------------
    # NO+ Coulomb collision frequencies
    nust[55]=ne*me/(ni[5]*ms[5])*nuet[5]           # NO+ - e
    nust[56]=coiNO[0]*ni[0]/(Ti[0]*math.sqrt(Ti[0]))    # NO+ - O+
    nust[57]=coiNO[1]*ni[1]/(Ti[1]*math.sqrt(Ti[1]))    # NO+ - H+
    nust[58]=coiNO[2]*ni[2]/(Ti[2]*math.sqrt(Ti[2]))    # NO+ - H+e
    nust[59]=coiNO[3]*ni[3]/(Ti[3]*math.sqrt(Ti[3]))    # NO+ - O2+
    nust[60]=coiNO[4]*ni[4]/(Ti[4]*math.sqrt(Ti[4]))    # NO+ - N2+

    # NO+ - neutral collision frequencies nu_{NO+,n}
    nust[61]= coiNO[5]*nn[0]                   # NO+ - O
    nust[62]= coiNO[6]*nn[1]                   # NO+ - O2
    nust[63]= coiNO[7]*nn[2]                   # NO+ - N2
    nust[64]= coiNO[8]*nn[3]                   # NO+ - H
    nust[65]= coiNO[9]*nn[4]                   # NO+ - He

    # collision frequency between neutral species
    Tn226=pow(Tn,0.226)
    nust[66]=con[0]*(nn[0]+2.0*nn[1])*Tn226             # O - O2
    nust[67]=con[1]*(4.0*nn[0]+7.0*nn[2])*Tn226         # O - N2 
    nust[68]=con[2]*(8.0*nn[1]+7.0*nn[2])*Tn**0.25      # O2 - N2
    nust[69]=con[3]*(nn[3]+16.0*nn[0])*Tn**0.292        # H - O
    nust[70]=con[4]*(nn[3]+32.0*nn[1])*Tn**0.289        # H - O2
    nust[71]=con[5]*(nn[3]+28.0*nn[2])*Tn**0.302        # H - N2
    nust[72]=con[6]*(nn[3]+4.0*nn[4])*Tn**0.294         # H - He
    nust[73]=con[7]*(nn[4]+4.0*nn[0])*Tn**0.251         # He - O
    nust[74]=con[8]*(nn[4]+8.0*nn[1])*Tn**0.290         # He - O2
    nust[75]=con[9]*(nn[4]+7.0*nn[2])*Tn**0.282         # He - N2

#-----------------------------------------------------------------------------
#--------------- thermal conductivities --------------------------------------
#-----------------------------------------------------------------------------
    lamda=np.zeros(8)

    Te2=Te*Te

    qn=np.zeros(5)
    qn[0]=1.1e-16*(1.0+5.7e-4*Te)
    qn[1]=2.2e-16*(1.0+0.036*Te12)
    qn[2]=2.82e-17*Te12*(1.0-1.21e-4*Te)
    qn[3]=5.47e-15*(1.0-1.35e-4*Te)
    qn[4]=5.6e-16

    # electron thermal conductivity
    nqd=nn[0]*qn[0]+nn[1]*qn[1]+nn[2]*qn[2]+nn[3]*qn[3]+nn[4]*qn[4]
    lamda[0]=1.233694e-11*Te2*Te12/(1.0+3.32e4*Te2/ne*nqd)

    # ion thermal conductivity

    #first calculate average neutral mass
    Nn=0.0
    rhon=0.0
    for l in range(7):
        Nn=Nn+nn[l]
        rhon=rhon+nn[l]*ms[l]

    mt=rhon/Nn

    nuss=np.zeros(6)

    for m in range(6):
        s0=11*m

        nqd=0.0
        for l in range(6):
            # ion - electron Coulomb interaction
            if l==0:
                nqd = nqd+nust[s0]*(3.0+1.5*me/ms[m])
            elif l > 0 and l < 6:
                # ion - ion Coulomb interactions
                if l<=m:
                    ll=l-1
                else:
                    ll=l
                            
                Dst=(3.0*ms[m]*ms[m]-0.2*ms[ll]*ms[ll]+0.1*ms[m]*ms[ll]) \
                    /((ms[m]+ms[ll])*(ms[m]+ms[ll]))

                nqd = nqd+nust[s0+l]*(Dst+1.5*ms[ll]/(ms[m]+ms[ll]))
            else:
                # ion - neutral interaction
                Dst=(3.0*ms[m]*ms[m]+mt*mt+1.6*ms[m]*mt)/((ms[m]+mt)*(ms[m]+mt))

                nqd = nqd+nust[s0+l]*(Dst+1.5*mt/(ms[m]+mt))

            Td=Ti[m]*math.sqrt(Ti[m])
            nuss[m]=1.27*ni[m]/(math.sqrt(ms[m])*Td)
            nqd=1.0+1.25*nqd/nuss[m]

            lamda[1+m]=4.96682e-13*Td/(math.sqrt(ms[m])*nqd)
            
        # normalized neutral thermal conductivity 
        Td=Tn**0.69
        lamda[7]= 7.59e-4*Td \
                  +1.0e-4*(3.93*Td+0.255*Tn-9.27)  \
                  +1.0e-4*(3.82*Td+0.190*Tn+5.14) \
                  +3.79e-3*Td \
                  +2.99e-3*Td

    return nuet, nust, lamda
