function caluvslope,fluxarc,z,namegal,namearc,namestar,fixbeta=fixbeta,noconfit=noconfit,wlmin=wlmin,wlmax=wlmax
  clight = 3.d18                ;angstrom/s
  mpc    = 3.086d24             ;cm/mpc

  ;convert to rest frame values
  lambda = fluxarc.lambda/(1.+z)
  flux   = fluxarc.flux*(1.+z)
  fluxerr = fluxarc.fluxerr*(1.+z)
  ;convert to magnitude
  fnu    = flux*lambda^2/clight
  fnuerr = fluxerr*lambda^2/clight
  mAB    = -2.5*alog10(fnu)-48.6
  mABerr = abs(2.5*fnuerr/fnu/alog(10.))
  
  ;make a mask 
  mask   = bytarr(n_elements(lambda))+1
  ;take only wavelength longer than lya+2000km/s(1225) and less than 1650.
  lyc = where(lambda lt 1225.,clyc)
  if clyc eq 0 then stop, 'Something wrong with wavelength'
  mask(lyc) = 0
  if ~keyword_set(wlmax) then wlmax=1600.
  toored = where(lambda gt wlmax,cred)
  if cred gt 0 then mask(toored) = 0
  if keyword_set(wlmin) then begin
     tooblue = where(lambda lt wlmin,cblu)
     if cblu gt 0 then mask(tooblue) = 0.
  endif
  ;block off the absorption features
  readcol, '/scr2/nichal/workspace3/lines.txt',linewave,fosc,linename,format ='D,D,A,X',comment='#'
  linestart= (linewave*(-500.)/3.e5+linewave)
  lineend  = (linewave*(300)/3.e5+linewave)
  for i=0,n_elements(linestart)-1 do begin
     w = where(lambda ge linestart[i] and lambda le lineend[i], c)
     if c gt 0 then mask[w] = 0
  endfor
  ;block off where error is infinity
  baddata = where(~finite(mABerr) or ~finite(mAB),cbaddata)
  if cbaddata gt 0 then mask[baddata] = 0
  ;block off telluric features
 ; readcol,'/scr2/nichal/workspace2/telluric/telluric.mask', tellstart, tellend, format='D,D', /silent, comment='#'
  readcol,'/scr2/nichal/workspace3/telluric.mask', tellstart, tellend, format='D,D', /silent, comment='#'
  for i=0,n_elements(tellstart)-1 do begin
     woff = where(lambda*(1.+z) gt tellstart(i) and lambda*(1.+z) lt tellend(i),cwoff)
     ;if cwoff gt 0 then mask(woff) = 0
  endfor

  won = where(mask ne 0,cwon)
  if cwon eq 0 then stop,'BOO'
  
  ;fit the function to the continuum fit (linear fit)
  ;fit continuum
  upperlim = 1.5
  lowerlim = 1.5
  if namegal eq 'macs0940' and namearc eq 'arcB' then upperlim = 0.5 
  if namegal eq 'macs0940' and namearc eq 'arcB' then lowerlim = 3.
  if namegal eq 'a2219' then upperlim = 0.5 
  if namegal eq 'a2219' then lowerlim = 7.
  if namegal eq 'h3a' then upperlim = 5. 
  if namegal eq 'h3a' then lowerlim = .5
  bkpt = slatec_splinefit(lambda[won], mAB[won], coeff, invvar=1./(mABerr[won])^2, bkspace=150, upper=upperlim, lower=lowerlim, /silent)
  if bkpt[0] eq -1 then stop,'CANNOT FIT CONTINUUM'
  cont = slatec_bvalu(lambda,bkpt,coeff)

  xmp = alog10(lambda[won])
  if keyword_set(noconfit) then ymp=mAB[won] else ymp = cont[won]
  dymp = mABerr[won]
  
  if n_Elements(fixbeta) eq 0 then begin 
     params = linfit(xmp,ymp,chisqr=chisq,covar=covar,/double,measure_errors=dymp,prob=prob,sigma=paramserr)
     if prob lt 0.1 then print, 'WARNING:LINEAR FIT IS POOR' else print,'GOOD LINEAR FIT'
     beta = (params[1]/(-2.5))-2.0
     betaerr = paramserr[1]/2.5
     yfit = -2.5*(beta+2.)*alog10(lambda)+params[0]
  endif
  if n_Elements(fixbeta) ne 0 then begin
     beta = fixbeta[0]
     betaerr = fixbeta[1]
     params  = poly_fit(xmp,ymp+2.5*(beta+2.)*xmp,0,measure_errors=dymp)
     yfit =  -2.5*(beta+2.)*alog10(lambda)+params[0]
  endif

  ;find the magnitude and flux at 1600A
  m1600     = -2.5*(Beta+2.)*alog10(1600.)+params[0]
  m1600err  = abs(2.5*alog10(1600.)*betaerr)
  f1600     = 10.^((m1600+48.6)/(-2.5))*clight/1600.^2 ;flambda
  f1600err  = f1600*(alog(10.)*m1600err/2.5)
  okrange   = where(lambda gt 1598 and lambda lt 1602.,cokrange)
  meanerr,flux(okrange),fluxerr(okrange),meanflux,sigmam,sigmad,sigmas
  f1600err  = sqrt(f1600err^2+sigmam^2)  

  ;SFR!!! YAY
  A16    = 2.31*beta+4.85       ;Calzetti,2000
  A16err = 2.31*betaerr
  f1600dered    = f1600*10.^(0.4*A16) ;flambda
  f1600derederr = sqrt((f1600dered*0.4*alog(10.)*A16err)^2+(10.^(0.4*A16)*f1600err)^2)
  
  f1600dered    = 1600.^2*f1600dered/clight  ;fnu in erg/s/hz/cm^2
  f1600derederr = 1600.^2*f1600derederr/clight ;fnu in erg/s/hz/cm^2
  distmpc       = lumdist(z)
  dist          = distmpc*mpc ;cm
  L1600dered    = f1600dered*4.*!dpi*dist^2
  L1600derederr = f1600derederr*4.*!dpi*dist^2

  MUV  = m1600-5.*alog10(distmpc*1.e6/10.)
  MUVerr = m1600err
  
  SFR    = 1.4d-28*L1600dered
  SFRerr = 1.4d-28*L1600derederr

  ;PLOT
  loglambda = alog10(lambda)
  set_plot,'ps'
  !p.font = 0 
  !P.multi = [0,1,1]
  !p.charsize = 1
  psname='/scr2/nichal/workspace3/flux_calib/'+namegal+'_'+namearc+'_'+namestar+'_beta.eps'
  device, filename = psname,xsize = 15,ysize = 10, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  plot,loglambda,mAB,xtitle='log(lambda)',ytitle='AB Magnitude',yrange=[20,40],xrange=minmax(loglambda),xstyle=1
  plotmask,loglambda,mask
  oplot,loglambda,mAB+MABerr,color=fsc_color('rose')
  oplot,loglambda,mAB-MABerr,color=fsc_color('rose')
  oplot,loglambda,mAB
  oplot,loglambda,cont,color=fsc_color('orange')
  oplot,loglambda,yfit,color=fsc_color('dark green')
  betaletter = "142B
  pmletter   = "261B
  muletter   = "155B
  xyouts,0.9,0.25,'!9'+string(betaletter)+'!X'+'='+strtrim(string(beta,format='(F5.2)'))+'!9'+string(pmletter)+'!X'+strtrim(string(Betaerr,format='(F4.2)')),alignment=1,/normal
  xyouts,0.9,0.18,'SFRx'+'!9'+string(muletter)+'!X'+'='+strtrim(string(fix(SFR),format='(I3)'))+'!9'+string(pmletter)+'!X'+strtrim(string(fix(SFRerr),format='(I3)')),alignment=1.,/normal

  device,/close
  print, namearc
  print, 'BETA =',beta,betaerr
  print, 'A1600 =',A16,A16err
  print, 'SFR =', SFR,SFRerr
  print, 'MUV =', MUV,MUVerr
  return,[beta,betaerr,f1600,f1600err,SFR,SFRerr]
end
