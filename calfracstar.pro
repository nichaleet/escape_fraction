function Fmt,m,t
common param, lmparams,lfparams,mionmax
  logm = alog10(m)
  lfuv_m = 10.^(lmparams[0]+lmparams[1]*logm+lmparams[2]*logm^2+lmparams[3]*logm^3)
  IMF = m^(-2.35) ;salpeter
  SFR = 10      
  return, SFR*IMF*lfuv_m
end
function ionizing_Limits, t
  common param, lmparams, lfparams,mionmax
  mmax = alog((t-lfparams[2])/lfparams[0])/alog(lfparams[1]) ;msun
  if mmax gt mionmax then mrange= [mionmax,mmax] else mrange= [mionmax,mionmax]
  return, mrange
end

function nonionizing_Limits, t
  common param, lmparams, lfparams,mionmax
  mmax = alog((t-lfparams[2])/lfparams[0])/alog(lfparams[1]) ;msun
  if mmax gt mionmax then mrange= [1.,mionmax] else mrange= [1.,mmax]
  return,mrange
end

pro calfracstar
common param, lmparams,lfparams,mionmax
mionmax =2.5
;function of L(M) from Parravano2003
logm = [0.075,0.12,0.17,0.24,0.29,0.38,0.46,0.58,0.68,0.76,0.83,0.93,1.07,1.16,1.46,1.6,1.76,2.,2.08]
logLFUV = [-2.8,-2.3,-1.5,-0.5,0,0.8,1.5,2.2,2.5,2.8,3.,3.5,3.8,4.2,4.7,5,5.3,5.6,5.8]
LMparams = poly_fit(logm,loglfuv,3,yfit=yfitlogfuv)
!p.multi=[0,1,2]
plot,10.^logm,10.^loglfuv,xtitle='log M',ytitle='log LFUV'
oplot, 10.^logm, 10.^yfitlogfuv,color=fsc_color('green')
;function of lifetime from Encyclopedia of astronomy http://www.astro.caltech.edu/~george/ay20/eaa-stellarmasses.pdf
mass     = [120.,60.,25,12,5.] ;solarmass 
lifetime = [2.56,3.45,6.51,16.0,94.5] ;Myr
a=[500.,0.9,10]
lfparams = comfit(mass,lifetime,a,/exponential)
plot,mass,lifetime,xtitle='mass',ytitle='life time in Myr'
massarr = findgen(120)
model = lfparams[0]*lfparams[1]^massarr+lfparams[2]
oplot,massarr,model,color=fsc_color('green')

;integration
AB_Limits=[5,100] ;Myr
lfuv_ionizing = int_2d('Fmt',AB_Limits,'ionizing_Limits',48,/double)
lfuv_nonionizing = int_2d('Fmt',AB_Limits,'nonionizing_Limits',48,/double)
frac=lfuv_nonionizing/(lfuv_ionizing+lfuv_ionizing)
print, frac
stop
end
