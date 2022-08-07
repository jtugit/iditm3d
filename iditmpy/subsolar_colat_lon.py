#######################################################################
# Evaluate subsolar colatitude and longitude with given year (4 digits)
# day of the year, and UT (seconds)
#----------------------------------------------------------------------
from math import floor , sin , cos , pi , atan2 , asin

def subsolar_colat_lon(year, doy, ut):
    yr = year - 2000
    nleap = floor ((year-1601/4.0))
    nleap = nleap - 99
    if year <= 1900: 
        ncent = floor ( ( year - 1601)/100.)
        ncent = 3 - ncent
        nleap = nleap + ncent

    l0 = -79.549- ( 0.238699*( yr - 4*nleap )- 3.08514e-2*nleap )
    g0 = -2.472 - ( 0.2558905*( yr  -4*nleap)+ 3.79617e-2*nleap )
    df = ( ut /86400. - 1.5)+ doy
    lf = .9856474 *df
    gf = .9856003*df
    l = l0 + lf 
    g = g0 + gf 
    grad = g*pi / 180.
    lmbda = l + 1.915*sin ( grad ) + .020 *sin ( 2. *grad )
    lmrad = lmbda*pi / 180.
    sinlm = sin ( lmrad )
    n = df + 365. *yr + nleap
    epsilon = 23.439  -4.0e-7*n
    epsrad = epsilon*pi / 180.
    alpha = atan2 ( cos ( epsrad )*sinlm , cos ( lmrad ) ) *180. / pi 
    delta = asin ( sin ( epsrad ) *sinlm ) *180. / pi 
    subsolar_colat = 90.0  -delta 

    etdeg = l - alpha
    nrot = round ( etdeg / 360. ) 
    etdeg = etdeg -360.*nrot 
    aptime = ut /240. + etdeg 
    subsolar_lon = 180. - aptime
    nrot = round ( subsolar_lon / 360. ) 
    subsolar_lon = subsolar_lon - 360.*nrot

    return subsolar_colat, subsolar_lon