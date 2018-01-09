pro covfracnew,objname,returnave,returnave_si,lineuse=lineuse,zspec=zspec

;Note: This program interpolates velocity to velocities of each pixel of Si1260 line

;INPUT:
;lineuse = according to line list below, what line to use for cov fraction and systemic redshift from low-ionization absorption line


;;;;;;;;;;;;;;;;;;;;;
;Constants
linelist =['SiII1260','OI1302','SiII1304','CII1334','SiII1526']
linewl   = [1260.42,1302.168,1304.37,1334.532,1526.72]
if n_Elements(lineuse) eq 0 then lineuse = [1,1,1,1,1]
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
set_plot,'x'
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
findredshift,science, objname, restlya, linelist,linewl, lineuse, clight, shiftlya, shiftISM, delv_is_sys, zoommin, zoommax,zoommin2,zoommax2,fullmin,fullmax, fscience 

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
psname='output/'+objname+'_covfracnew'+strjoin(strtrim(string(lineuse),2))+'.eps'
device, filename = psname,xsize = 22,ysize = 22, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
!p.multi = [0,2,3]
!p.font = 0
!p.charsize = 2
;Draw spectra at each line
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
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;finish drawing 5 lines. Now going to combine the lines. 


;Calculate the velocities of each pixel for a referent line (SiII1260 unless specified)
 if n_Elements(lineref) eq 0 then lineref = 'SiII1260'
 ref = where(linelist eq lineref,cref)
 if cref eq 0 then stop, 'CANNOT FIND REFERENCE LINE'
 linerefwl = (linewl(ref)*(z+1))[0]
 velrefarr = (science.lambda-linerefwl)/linerefwl*clight ;translate wl to velocity
 goodvel   = where(velrefarr gt -1000. and velrefarr lt 500., cgoodvel)
 if cgoodvel ne 0 then velrefarr = velrefarr(goodvel) else stop, 'THERE IS A PROBLEM'

 ;Do SiII1260
 wlSiII = 1260.42*(z+1.)
 velsi = (science.lambda-wlSiII)/wlSiII*clight
 Si1260       = interpol(science.contdiv,velsi,velrefarr,/spline)
 Si1260_ivar  = interpol(Science.contdiv,velsi,velrefarr,/spline)
 dSi1260 = 1./sqrt(Si1260_ivar)

 ;going to start the loop
 goodlines = linewl(where(lineuse eq 1))*(z+1)
 goodlines_name = linelist(where(lineuse eq 1))
 nlines = n_Elements(goodlines)
 nvel   = n_elements(velrefarr)
 absarr_all = dblarr(nvel,nlines)
 absarrivar_all = dblarr(nvel,nlines)

 vminSi1304 = -200. ;km
 vbadsi = where(velrefarr lt vminsi1304,cbadsi)
 ls = 1304.37
 lo = 1302.168
 vmaxOi1302 = vminsi1304*ls/lo+clight*(ls-lo)/lo
 vbadoi = where(velrefarr gt vmaxoi1302,cbadoi)

 ;begin loop over lines
 for nl = 0, nlines-1 do begin
    ;translate wl to velocity
    curvel = (science.lambda-goodlines(nl))/goodlines(nl)*clight
    goodvel = where(curvel gt -1100. and curvel lt 600.,cgoodvel)
    if cgoodvel eq 0 then stop, 'NO wl FOUND FOR THE WHOLE VELOCITY'

    ;interpolate to velrefarr
    contdivmatch = interpol(science.contdiv,curvel,velrefarr,/spline)
    contdivivarmatch = interpol(science.contdivivar,curvel,velrefarr,/spline)

    ;mark the ivar of overlap region of siII1304 with OI1302 to be 0.
    if goodlines_name(nl) eq 'SiII1304' and cbadsi gt 0 then contdivivarmatch(vbadsi) = 0
    if goodlines_name(nl) eq 'OI1302' and cbadoi gt 0 then contdivivarmatch(vbadoi) = 0

    ;add the values to the big array
    absarr_all[*,nl] = contdivmatch
    absarrivar_all[*,nl] = contdivivarmatch  
 endfor

 ;combine the lines
 absarr_ivar = total(absarrivar_all,2)
 absarr = total(absarr_all*absarrivar_all,2)/absarr_ivar
 dabsarr = 1./sqrt(absarr_ivar)

 ;PLOTS (the last panel in the covfracnew11111.fits)
 plot,velrefarr,absarr,psym=10,title=objname,ytitle='Average Absorbtion profile',xtitle='velocity(km/s)',/nodata,yrange=[0,1.5]
 ;shade uncertainties 
 ;cgcolorfill,[velrefarr,reverse(velrefarr)],[Si1260+dSi1260,reverse(Si1260-dSi1260)],color=fsc_color('cyan')
 cgcolorfill,[velrefarr,reverse(velrefarr)],[absarr+dabsarr,reverse(absarr-dabsarr)],color=fsc_color('rose') 

 ;overplot with lines
 oplot,velrefarr,si1260,color=fsc_color('blue')
 oplot,velrefarr,absarr,psym=10,color=fsc_color('red')
 al_Legend,['SIII1260','Average line profile'],psym=[15,15],color=[fsc_color('blue'),fsc_color('red')],box=0,thick=2,charsize=0,symsize=1.5,/right,/bottom
 device,/close
 set_plot,'x'

 returnave=[[velrefarr],[absarr],[absarr_ivar]]
 returnave_si = [[velrefarr],[si1260],[si1260_ivar]]
 if tag_exist(science,'ave_covfracnew') then begin
    science.ave_covfracnew=returnave
    science.ave_covfracnew_si=returnave_si
 endif else science=create_struct(science,'ave_covfracnew',returnave,'ave_covfracnew_si',returnave_si)
 mwrfits,science,fscience,/create,/silent
stop
end
