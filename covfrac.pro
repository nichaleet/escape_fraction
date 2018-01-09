pro covfrac,objname,returnave,returnave_si,lineuse=lineuse,zspec=zspec,rangevel=rangevel,redshiftdomain=redshiftdomain

;Note: This calculate the weighted mean flux in each velocity bin. See covfracnew which interpolate velocity to velocities of each pixel of Si1260 line

;INPUT:
;lineuse = according to line list below, what line to use for cov fraction and systemic redshift from low-ionization absorption line


;;;;;;;;;;;;;;;;;;;;;
;Constants
linelist =['SiII1260','OI1302','SiII1304','CII1334','SiII1526']
linewl   = [1260.42,1302.168,1304.37,1334.532,1526.72]
;osc = [1.22,0.052,0.0928,0.129,0.133] ;nist
;tau is greater than 1 when N>2.45039e+11  5.56468e+12  3.11288e+12  2.18873e+12  1.85567e+12
if n_Elements(lineuse) eq 0 then lineuse = [1,1,1,1,1]
if n_elements(rangevel) eq 0 then rangevel = [-1000.,500.]
use_zlya  = 0 
use_zis   = 1
clight = 3.e5
shiftlya = -330. ;km Jones12
shiftISM = 190.  ;km Jones12
delv_is_sys = 125. ;km/s Jones13
restlya  = 1215.67
zoommin = 1220*(zspec+1)  ;SiII
zoommax = 1283*(zspec+1)
zoommin2 = 1500*(zspec+1) ;SiII
zoommax2 = 1600*(zspec+1)
fullmin = 1200*(zspec+1)
fullmax = 1560*(zspec+1)
fscience = '/scr2/nichal/workspace3/SCIENCE/'+objname+'_science.fits'
fscienceall = '/scr2/nichal/workspace3/SCIENCE/'+objname+'_scienceall.fits'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
setplot,14
science=mrdfits(fscience,1,/silent)
allscience = mrdfits(fscienceall,1,/silent)
science.lineuse = lineuse

;-------------------------------------------------------------------------
;PLOTS
!p.multi=[0,1,3]
plot,science.lambda,science.contdiv,yrange=[0,8],xrange=[fullmin,fullmax],xstyle=1
plotuvlines,science.zspec
plot,science.lambda,science.contdiv*sqrt(science.contdivivar),xrange=[fullmin,fullmax],xstyle=1,ytitle='SN per pixel'
plotuvlines,science.zspec
plot,science.lambda,science.sn,xrange=[fullmin,fullmax],xstyle=1,ytitle='SN per dlam'
plotuvlines,science.zspec
print, 'Median S/N =', median(science.sn)
wait,1.

linewl   = [1260.42,1302.168,1334.532,1526.72]
window,1,xsize=1200,ysize=900
!p.multi=[0,1,4]
for jj=0,3 do begin
   zoomlinemin = (1.-1500/3.e5)*linewl(jj)*(science.zspec+1)
   zoomlinemax = (1.+1500/3.e5)*linewl(jj)*(science.zspec+1)
   toplot = where(science.lambda lt zoomlinemax and science.lambda gt zoomlinemin)
   plot,science.lambda,science.contdiv,xrange=[zoomlinemin,zoomlinemax],/nodata,title=science.objname[0]
   skyspec = science.skyspec[*,0]
   skytoplot = skyspec(toplot)
   skytoplot = skytoplot/max(skytoplot)*(max(!y.crange)-min(!y.crange))
   oplot,science.lambda(toplot),skytoplot,color=fsc_color('pink')
   oplot,science.lambda,science.contdiv
   oplot,!x.crange,[1.,1.],color=fsc_color('green')
   plotuvlines,science.zspec
endfor
linewl   = [1260.42,1302.168,1304.37,1334.532,1526.72]
;-------------------------------------------------------------------------
;-------------------------------------------------------------------------
;REDSHIFT 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Find REDSHIFTS

if n_Elements(redshiftdomain) eq 0 then findredshift,science, objname, restlya, linelist,linewl, lineuse, clight, shiftlya, shiftISM, delv_is_sys, zoommin, zoommax,zoommin2,zoommax2,fullmin,fullmax, fscience else copyredshift,science,redshiftdomain

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
set_plot,'x'
mwrfits,science,fscience,/create,/silent
stop
;-------------------------------------------------------------------------
;-------------------------------------------------------------------------
;-------------------------------------------------------------------------
;COVERING FRACTION
;Average Low-Ionization Absortion profile
if use_zlya eq 1 then science.zspec = science.zsys_lya
if use_zis  eq 1 then science.zspec = science.zsys_ism_ave 
if use_zlya eq 1 and use_zis eq 1 then stop,'You have to choose type of systemic redshift'
if use_zlya eq 0 and use_zis eq 0 then stop,'You have to choose type of systemic redshift'

z = science.zspec
set_plot,'ps'
psname='output/'+objname+'_covfrac'+strjoin(strtrim(string(lineuse),2))+'.eps'
device, filename = psname,xsize = 22,ysize = 22, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
!p.multi = [0,2,3]
!p.font = 0
!p.charsize = 2
for i=0,n_elements(linelist)-1 do begin
   wl   = linewl(i)
   minwl = (1.-1000/3.e5)*wl*(z+1)
   maxwl = (1.+1000/3.e5)*wl*(z+1)
   goodwl = where(science.lambda lt maxwl and science.lambda gt minwl)
   lambdanow  = science.lambda(goodwl)
   contdivnow = science.contdiv(goodwl)
   ivarnow    = science.contdivivar(goodwl)
   snnow      = science.sn(goodwl)
   plot,lambdanow,contdivnow,title=linelist(i),ytitle='Normalized Spectrum'
   oplot,lambdanow,snnow/5.,color=fsc_color('cyan')
   plotuvlines,z
   if i eq 0 then  al_Legend,['Spec','SN/5'],psym=[0,0],color=[fsc_color('black'),fsc_color('cyan')],box=0,thick=2,charsize=0,symsize=1.5,/right,/bottom

   if lineuse(i) eq 0 then begin
      ;draw a cross to mark out the line
      oplot,!x.crange,!y.crange,color=fsc_color('dark red'),thick=2
      oplot,!x.crange,reverse(!y.crange),color=fsc_color('dark red'),thick=2
   endif
endfor
 mwrfits,science,fscience,/create,/silent

;make velocity step in mean (fwhm of dlam)
 goodlines = linewl(where(lineuse eq 1))*(z+1)
 goodlines_name = linelist(where(lineuse eq 1))
 nlines = n_Elements(goodlines)
 stepv  = mean(science.dlam)/mean(science.lambda)*clight*2.355
 nstepv = round(abs(rangevel[1]-rangevel[0])/stepv)
 velarr = dblarr(nstepv)
 absarr = dblarr(nstepv)
 absarr_wivar = dblarr(nstepv)
 Si1260       = dblarr(nstepv)
 Si1260_wivar  = dblarr(nstepv)
 wlSiII = 1260.42*(z+1.)
 velsi = (science.lambda-wlSiII)/wlSiII*clight
 vminSi1304 = -200. ;km
 ls = 1304.37
 lo = 1302.168
 vmaxOi1302 = vminsi1304*ls/lo+clight*(ls-lo)/lo
 for nv=0,nstepv-1 do begin
    contdivnow = []
    ivarnow    = []
    for nl=0,nlines-1 do begin
       vnowmin = rangevel[0]+nv*stepv
       vnowmax = rangevel[0]+(nv+1.)*stepv
       velarr(nv) = 0.5*(vnowmin+vnowmax)
       if goodlines_name(nl) eq 'SiII1304' then vnowmin = -1000.+nv*stepv > vminsi1304
       if goodlines_name(nl) eq 'OI1302' then vnowmax = -1000.+(nv+1.)*stepv < vmaxOi1302 
       wlnow  = (science.lambda-goodlines(nl))/goodlines(nl)*clight;translate wl to velocity
       goodwl = where(wlnow gt vnowmin and wlnow lt vnowmax,cwl)
       if cwl gt 0 then begin
          contdivnow = [contdivnow,science.contdiv(goodwl)]
          ivarnow    = [ivarnow,science.contdivivar(goodwl)]
       endif else print, 'no wl found for',goodlines_name(nl),'at velocity ',strtrim(string(velarr(nv)),2)
    endfor
    if n_Elements(contdivnow) gt 0 then begin
       absarr(nv) = total(contdivnow*ivarnow)/total(ivarnow)
       ;absarr_wivar(nv) = total(ivarnow)/total(ivarnow*(contdivnow-absarr(nv))^2) ;weighted standard deviation
       absarr_wivar(nv) = total(ivarnow) ;variance of weighted mean $
       ;(This has stronger assumption on that the flux is constant over the velocity bin of dlam)           
    endif
    vnowmin = rangevel[0]+nv*stepv
    vnowmax = rangevel[0]+(nv+1.)*stepv
    goodsi = where(velsi gt vnowmin and velsi lt vnowmax,csi)
    si1260(nv) = total(science.contdiv(goodsi)*science.contdivivar(goodsi))/total(science.contdivivar(goodsi))
    ;si1260_wivar(nv) = total(science.contdivivar(goodsi))/total(science.contdivivar(goodsi)*(science.contdiv(goodsi)-si1260(nv))^2)
    si1260_wivar(nv) = total(science.contdivivar(goodsi))
 endfor
 ;make finer grids to shade uncertainties in the plots
 nfine = n_elements(velarr)*2-1
 velarr_fine = dblarr(nfine)
 upperabs    = dblarr(nfine)
 lowerabs    = dblarr(nfine)
 Siupperabs    = dblarr(nfine)
 Silowerabs    = dblarr(nfine)
 velarr_fine[0] = velarr[0]
 for i=0,nfine-1 do begin
    if i mod 2 eq 0 then begin
       velarr_fine(i) = velarr_fine(i-1>0)
       upperabs[i]    = absarr(i/2)+0.5*sqrt(1./absarr_wivar(i/2))
       lowerabs[i]    = absarr(i/2)-0.5*sqrt(1./absarr_wivar(i/2))
       Siupperabs[i]    = Si1260(i/2)+0.5*sqrt(1./Si1260_wivar(i/2))
       Silowerabs[i]    = Si1260(i/2)-0.5*sqrt(1./Si1260_wivar(i/2))
    endif
    if i mod 2 eq 1 then begin
       velarr_fine(i) = 0.5*(velarr((i-1)/2)+velarr((i+1)/2))
       upperabs[i]    = upperabs[i-1]
       lowerabs[i]    = lowerabs[i-1]
       Siupperabs[i]  = Siupperabs[i-1]
       Silowerabs[i]  = Silowerabs[i-1]
    endif
 endfor
 plot,velarr,absarr,psym=10,title=objname,ytitle='Average Absorbtion profile',xtitle='velocity(km/s)',/nodata,yrange=[0,1.5]
 cgcolorfill,[velarr_fine,reverse(velarr_fine),velarr_fine[0]],[Siupperabs,reverse(Silowerabs),Siupperabs[0]],color=fsc_color('cyan')
 cgcolorfill,[velarr_fine,reverse(velarr_fine),velarr_fine[0]],[upperabs,reverse(lowerabs),upperabs[0]],color=fsc_color('rose')
 oplot,velarr,si1260,psym=10,color=fsc_color('blue')
 oplot,velarr,absarr,psym=10,color=fsc_color('red')
 al_Legend,['SIII1260','Average line profile'],psym=[15,15],color=[fsc_color('blue'),fsc_color('red')],box=0,thick=2,charsize=0,symsize=1.5,/right,/bottom
 device,/close
 set_plot,'x'

 returnave=[[velarr],[absarr],[absarr_wivar]]
 returnave_si = [[velarr],[si1260],[si1260_wivar]]
 if objname eq 'ms1358arcac' then begin
    badmark = where(returnave[*,1] lt 0.)
    returnave[badmark,1] = 0.190
 endif
 if tag_exist(science,'ave_covfrac') then begin
    science.ave_covfrac=returnave
    science.ave_covfrac_si=returnave_si
 endif else science=create_struct(science,'ave_covfrac',returnave,'ave_covfrac_si',returnave_si)

 mwrfits,science,fscience,/create,/silent
stop
end
