pro runsfr_jones
  clight = 3.d18                ;angstrom/s
  mpc    = 3.086d24             ;cm/mpc
  betaarr = fltarr(5)           ;H3a, H3b, H5a, H5b, J1261
  betaerrarr = fltarr(5)
  magarr = fltarr(5)
  magerrarr = fltarr(5)
;Abell2390
  zH3 = 4.043
  zH5 = 4.0448
  wlobs = [641.,669.,799.,832,1270.]
;H3 Pello
  wl    = alog10(wlobs/(zH3+1))
  mH3A  = [22.8,22.9,22.1,22.1,21.0]
  mH3B  = [23,23.4,22.5,22.3,21.6]
  mH3 = [[mH3A],[mH3B]]
  for i=0,1 do begin
     mnow = mh3[*,i]
     params = linfit(wl,mnow,sigma=paramserr)
     beta = (params[1]/(-2.5))-2.0
     betaerr = paramserr[1]/2.5
     plot,wl,mH3A,psym=1
     oplot,!x.crange,params[0]+params[1]*!x.crange
     print, 'H3 Beta=',beta,betaerr
  endfor
;H3a Bunker2000
  wl = alog10([1585,3300])
  mag = [23.05,22.85]
  magerr = [0.05,0.06]
  params = linfit(wl,mag,sigma=paramserr,measure_error=magerr)
  beta = (params[1]/(-2.5))-2.0
  betaerr = paramserr[1]/2.5
  plot,wl,mag,psym=1
  oplot,!x.crange,params[0]+params[1]*!x.crange
  print, 'H3 Beta=',beta,betaerr
  ;use this beta
  betaarr[0:1] = beta
  betaerrarr[0:1] = betaerr
  magarr[0] = 22.4 ;J13
;  magarr[0] = mag[0] 
  magerrarr[0] = 0.05
  magarr[1] = magarr[0]+0.4 ;compare to Pello1990
  magerrarr[1] = magerr[0]

;H5 Pello
  wl    = alog10(wlobs/(zH5+1))
  mH5A  = [23.5,23.9,23.7,23.3,23.0]
  mH5B  = [23.1,23.1,23.0,22.5,22.3]
  mH5 = [[mH5A],[mH5B]]
  for i=0,1 do begin
     mnow = mh5[*,i]
     params = linfit(wl,mnow,sigma=paramserr)
     beta = (params[1]/(-2.5))-2.0
     betaerr = paramserr[1]/2.5
     plot,wl,mH5A,psym=1
     oplot,!x.crange,params[0]+params[1]*!x.crange
     print, 'H5 Beta=',beta,betaerr
  endfor
;H5a Sextracted Cat
  wl = alog10([12584.8,8186.4]/(4.045+1.))
  mag = [24.871,24.279]
  magerr = [0.05,0.05]
  params = linfit(wl,mag,sigma=paramserr,measure_error=magerr)
  beta = (params[1]/(-2.5))-2.0
  betaerr = paramserr[1]/2.5
  plot,wl,mag,psym=1
  oplot,!x.crange,params[0]+params[1]*!x.crange
  print, 'H5a Beta=',beta,betaerr
  magarr[2] = 22.8+0.7 ;j13 for H5b and Pello for mag diff
  magerrarr[2] = magerr[1]

;H5b Sextracted Cat
  wl = alog10([12584.8,8186.4]/(4.045+1.))
  mag = [23.656,23.279]
  magerr = [0.05,0.05]
  params = linfit(wl,mag,sigma=paramserr,measure_error=magerr)
  beta = (params[1]/(-2.5))-2.0
  betaerr = paramserr[1]/2.5
  plot,wl,mag,psym=1
  oplot,!x.crange,params[0]+params[1]*!x.crange
  print, 'H5b Beta=',beta,betaerr
  ;use this beta
  betaarr[2:3] = beta
  betaerrarr[2:3] = betaerr
  magarr[3] = 22.8 ;J13
  magarr[3] = mag[1]
  magerrarr[3] = 0.05

;J1261 Sextracted Cat
  wl = alog10([11029.6,15235.9]/(4.13+1.))
  mag = [23.05,23.017]
  magerr = [0.05,0.05]
  params = linfit(wl,mag,sigma=paramserr,measure_error=magerr)
  beta = (params[1]/(-2.5))-2.0
  betaerr = paramserr[1]/2.5
  plot,wl,mag,psym=1
  oplot,!x.crange,params[0]+params[1]*!x.crange
  print, 'J1261 Beta=',beta,betaerr
  ;use this beta
  betaarr[4] = beta
  betaerrarr[4] = betaerr
  magarr[4] = 22.6
  magerrarr[4] = 0.05

;J1621 Tucker and Bayliss
  wl = alog10([799.,630.]/(4.13+1.))
  mag = [22.6,22.28]
  magerr = [0.05,0.05]
  params = linfit(wl,mag,sigma=paramserr,measure_error=magerr)
  beta = (params[1]/(-2.5))-2.0
  betaerr = paramserr[1]/2.5
  plot,wl,mag,psym=1
  oplot,!x.crange,params[0]+params[1]*!x.crange
  print, 'J1261 Beta=',beta,betaerr


  zarr=[zh3,zh3,zh5,zh5,4.13]
  name = ['H3a','H3b','H5a','H5b','J1261']
  mag = [19.5,9.4,4.8,6.5,38.]
  magerr = 0.2*mag
  magerr(where(mag gt 10.)) = 0.5*mag(where(mag gt 10.))
  print, 'beta',betaarr
  for i=0,4 do begin
     A16    = 2.31*betaarr[i]+4.85    ;Calzetti,2000
     A16err = 2.31*betaerrarr[i]
     m1600  = magarr[i]
     m1600err = magerrarr[i]
     f1600  = double(10.^((m1600+48.6)/(-2.5))) ;fnu in cgs  
     f1600err = abs(f1600*0.4*alog(10)*0.05)
     f1600dered    = f1600*10.^(0.4*A16) 
     f1600derederr = sqrt((f1600dered*0.4*alog(10.)*A16err)^2+(10.^(0.4*A16)*f1600err)^2)
     distmpc       = lumdist(zarr[i],/silent)
     dist          = distmpc*mpc ;cm
     L1600dered    = f1600dered*4.*!dpi*dist^2
     L1600derederr = f1600derederr*4.*!dpi*dist^2

     MUV  = m1600-5.*alog10(distmpc*1.e6/10.)
     MUVerr = m1600err
  
     SFR    = 1.4d-28*L1600dered/mag[i]
     SFRerr = SFR*sqrt((L1600derederr/L1600dered)^2+(magerr[i]/mag[i])^2)
     print,name[i]+' lensed SFR=',1.4d-28*L1600dered,1.4d-28*L1600derederr

     print,name[i]+' unlensed SFR=',SFR,SFRerr
  endfor

stop

end
