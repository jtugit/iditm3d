def shapiro_r(i0, i1, cshap, z):
    y=[0.0]*(i1-i0)

    for i in range(i0, i1):
        y[i-i0]=z[i-i0]

    for i in range(i0,i1):
        if i>0 and i < i1-1:
            z[i-i0]=(1.0-cshap)*y[i-i0]+0.5*cshap*(y[i-i0-1]+y[i-i0+1])

    return z
