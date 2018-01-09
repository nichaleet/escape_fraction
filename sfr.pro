function sfr, beta,betaerr,f1600,f1600err,z
mpc    = 3.086d24
clight = 3.e18
A16    = 2.31*beta+4.85 ;Calzetti,2000
A16err = 2.31*betaerr
f1600dered    = f1600*10.^(0.4*A16)  ;flambda
f1600derederr = sqrt((f1600dered*0.4*alog(10.)*A16err)^2+(10.^(0.4*A16)*f1600err)^2)

f1600dered    = 1600.^2*f1600dered/clight    ;fnu in erg/s/hz/cm^2
f1600derederr = 1600.^2*f1600derederr/clight ;fnu in erg/s/hz/cm^2

dist          = lumdist(z)*mpc ;cm
L1600dered    = f1600dered*4.*!dpi*dist^2
L1600derederr = f1600derederr*4.*!dpi*dist^2

SFR    = 1.4d-28*L1600dered
SFRerr = 1.4d-28*L1600derederr
print, 'SFR =', SFR,SFRerr
return,[SFR,SFRerr]
end
