pro covfrac,file1ds,objname,returnave,returnave_si,zspec=zspec,smoothfactor=smoothfactor,redo=redo
linelist =['SiII1260','OI1302','SiII1304','CII1334','SiII1526']
linewl   = [1260.42,1302.168,1304.37,1334.532,1526.72]
lineuse = [1,1,1,1,0]
use_zlya  = 0 
use_zis   = 1
clight = 3.e5
shiftlya = -330. ;km Jones12
shiftISM = 190.  ;km Jones12
restlya  = 1215.67
zoommin = 1220*(zspec+1)  ;SiII
zoommax = 1283*(zspec+1)
zoommin2 = 1500*(zspec+1) ;SiII
zoommax2 = 1600*(zspec+1)
fullmin = 1200*(zspec+1)
fullmax = 1560*(zspec+1)
if ~keyword_set(smoothfactor) then smoothfactor = 2.
nmasks = n_Elements(file1ds)
fscience = '/scr2/nichal/workspace/SCIENCE/'+objname+'_science.fits'
fscienceall = '/scr2/nichal/workspace/SCIENCE/'+objname+'_scienceall.fits'

scienceall = []
if ~file_test(fscience) or keyword_set(redo) then begin
   for n=0,nmasks-1 do begin
      file1d = file1ds(n)
                                ;read spec1d file
      blu = mrdfits(file1d,1,hdrblue)
      red = mrdfits(file1d,2,hdrred)
      ;find overlapping
      badred = where(red.lambda lt max(blu.lambda),cbadred)
      mingood = max(badred)+1
      if cbadred gt 0 then begin
         tagnames = tag_names(red)
         nold  = n_elements(red.lambda)
         for i=0,n_tags(red)-1 do begin
            if n_elements(red.(i)) gt 1 then newvals=(red.(i))[mingood:nold-1] else newvals = red.(i)
            if i eq 0 then newstruct = create_struct(tagnames(i),newvals) else newstruct = create_struct(newstruct,tagnames(i),newvals)
         endfor
         red=newstruct
      endif
      
                                ;make science structure
      npix = n_elements(blu.lambda)+n_elements(red.lambda)
      nsky = 200
      science = {objname:objname, $
                 spec1dfile:file1d, $
                 objpos:[-999d,-999d], $
                 fwhm:[-999d,-999d], $
                 exptime:[-999d,-999d],$
                 lambda:dblarr(npix), $
                 spec:dblarr(npix), $
                 origspec:dblarr(npix), $
                 smoothed:dblarr(npix), $
                 ivar:dblarr(npix), $
                 skyspec:dblarr(npix), $
                 telldiv:dblarr(npix), $
                 telldivivar:dblarr(npix), $
                 contmask:bytarr(npix), $
                 cont:dblarr(npix), $
                 contdiv:dblarr(npix), $
                 contdivivar:dblarr(npix), $
                 fitmask:bytarr(npix), $
                 dlam:dblarr(npix), $
                 skyfit:[[dblarr(nsky)], [dblarr(nsky)], [dblarr(nsky)]], $
                 skylinemask:lonarr(nsky),$
                 goodsky:0, $
                 airmass:-999d, $
                 zspec:-999d, $
                 sn:-999d, $
                 smoothfactor:smoothfactor}
      science.objpos = [blu.objpos,red.objpos]
      science.fwhm   = [blu.fwhm,red.fwhm]
      science.exptime= [sxpar(hdrblue,'exptime'),sxpar(hdrred,'exptime')]
      science.lambda = [blu.lambda,red.lambda]
      science.spec   = [blu.spec,red.spec]
      science.ivar   = [blu.ivar,red.ivar]
      science.skyspec= [blu.skyspec,red.skyspec]
      science.zspec  = zspec
      science.airmass = sxpar(hdrblue,'AIRMASS')

      ;spectral resolution
      science.skylinemask = -1
      specres,science,qf=specres_poly
      save,specres_poly,filename='/scr2/nichal/workspace/SCIENCE/'+objname+'_specres_poly.sav'

      ;find signal to noise for each lambda in a spectral resolution
      sn = calsn(science,smoothfactor)
      loadct,2
      set_plot,'ps'
      psname=objname+strtrim(string(n),2)+'_signal_noise.eps'
      device, filename = psname,xsize = 12,ysize = 12, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      multiplot, /init
      multiplot, [1,4]
      !p.font = 0
                                ;plot, science.lambda,1./sqrt(science.ivar),xrange=[6500,8000],yrange=[10,20],ytitle='S/N'
      plot, science.lambda,sn,xrange=[6500,8000],yrange=[0,10],ytitle='S/N'
      multiplot
      plot, science.lambda,science.spec,yrange=[-500,500],xrange=[fullmin,fullmax],ytitle='spec'
      plotlines,science.zspec
      multiplot
      plot, science.lambda,science.dlam*smoothfactor,xrange=[fullmin,fullmax],ytitle='Angstrom'
      oplot,science.lambda,science.dlam,linestyle=2
      xyouts,[0.92*mean(!x.crange)],[mean(science.dlam)*smoothfactor],'smoothing range used for SN'
      xyouts,[mean(!x.crange)],[mean(science.dlam)],'spectral resolution'
          
      ;smooth with twice spectral resolution
      ;science.smoothed = median_smooth(science.lambda,science.telldiv,science.lambda,2.*science.dlam)
      science.origspec = science.spec     
      science.spec = smooth_gauss_wrapper(science.lambda,science.spec,science.lambda,smoothfactor*science.dlam)
      multiplot
      plot, science.lambda, science.spec,ytitle='smoothed spec',xrange=[fullmin,fullmax],xtitle='wavelength(A)'
      plotlines,science.zspec
      multiplot,/reset,/default
      device,/close
      set_plot,'x'
      setplot,14
      plot,science.lambda,science.origspec,xrange=[zoommin,zoommax],yrange=[-1000,1000],title=file1d,/nodata,xtitle='CHECK SMOOTH and SKY'
      oplot,science.lambda,science.skyspec,color=fsc_color('pink')
      oplot,science.lambda,science.origspec
      oplot,science.lambda,science.spec,color=fsc_color('green')      
      plotlines,science.zspec
      stop

      ;telluric
      ;tellfile = file_Search('/scr2/nichal/workspace/SPEC1D/'+objname+'_telluric.fits',count=ctell)
      ;if ctell ne 0 then tellfile=tellfile(0) else stop, 'Stop: No telluric file found.'      
      ;telluric,science,tellfile
      ;telluric_enk,science
      science.telldiv = science.spec
      science.telldivivar = science.ivar

      ;fit continuum
      fitcontinuum, science
      !p.multi=[0,2,3]
      plot,science.lambda,science.spec,xrange=[fullmin,fullmax],yrange=[-300,300],title='CHECK CONTINUUM'
      oplot,science.lambda,science.cont,color=50
      plotlines,science.zspec
      linewl   = [1260.42,1302.168,1304.37,1334.532,1526.72]
      for jj=0,4 do begin
         zoomlinemin = (1.-1500/3.e5)*linewl(jj)*(science.zspec+1)
         zoomlinemax = (1.+1500/3.e5)*linewl(jj)*(science.zspec+1)
         plot,science.lambda,science.spec,xrange=[zoomlinemin,zoomlinemax],/nodata
         oplot,science.lambda,science.skyspec,color=fsc_color('pink')
         oplot,science.lambda,science.spec
         oplot,science.lambda,science.cont,color=fsc_color('green')
         plotlines,science.zspec
      endfor

      scienceall = [scienceall,science]
      stop
      !p.multi=[0,1,1]

   endfor

   ;combine each frames
   if nmasks gt 1 then begin
      nlamb = n_Elements(scienceall[0].lambda)
      contdiv = dblarr(nlamb,nmasks)
      ivar    = dblarr(nlamb,nmasks)
      dlam    = dblarr(nlamb,nmasks)
      cont    = dblarr(nlamb,nmasks)
      for i=0,nmasks-1 do begin
         science = scienceall[i]
         contdiv[*,i] = science.contdiv
         ivar[*,i]    = science.ivar
         dlam[*,i]    = science.dlam
         cont[*,i]    = science.cont
      endfor
      s_ivar = total(ivar,2)
      fin_dlam = mean(dlam,dimension=2)
      fin_contdiv = total(contdiv*ivar,2)/s_ivar
      fin_cont    = total(cont*ivar,2)/s_ivar
      result={contdiv:fin_contdiv,lambda:science.lambda,ivar:s_ivar,dlam:fin_dlam,cont:fin_cont,zspec:science.zspec}
      science = result
      !p.multi=[0,1,2]
      plot,science.lambda,science.contdiv,xrange=[zoommin,zoommax]
      oplot,science.lambda,science.contdiv
      plotlines,science.zspec
      plot,science.lambda,fin_dlam,ytitle='mean spectral resolution'
      stop
   endif 
   mwrfits,science,fscience,/create,/silent
   mwrfits,scienceall,fscienceall,/create,/silent
endif  
;-------------------------------------------------------------------------
;FINISH PREPARING SCIENCE (THE REDO PART)
;-------------------------------------------------------------------------
setplot,14
science=mrdfits(fscience,1,/silent)
allscience = mrdfits(fscienceall,1,/silent)
;Plot Continuum 
set_plot,'ps'
psname=objname+'_continuum.eps'
device, filename = psname,xsize = 25,ysize = 8, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
!p.font = 0
!p.charsize = 1
!p.multi=[0,1,1]
if nmasks gt 1 then spec = total(allscience.spec*allscience.ivar,2)/total(allscience.ivar,2) else spec=allscience.spec
plot,science.lambda,spec,xrange=[fullmin,fullmax],xstyle=1
plotlines,science.zspec
oplot,science.lambda,science.cont,color=fsc_Color('red')
device,/close
set_plot,'x'
;-------------------------------------------------------------------------
;-------------------------------------------------------------------------
;-------------------------------------------------------------------------
;REDSHIFT 

!p.multi=[0,1,2]
plot,science.lambda,science.contdiv,yrange=[0,8],xrange=[fullmin,fullmax],xstyle=1
plotlines,science.zspec
plot,science.lambda,science.cont*sqrt(science.ivar),xrange=[fullmin,fullmax],xstyle=1,ytitle='SN'
plotlines,science.zspec
wait,1.
;stop

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Find Redshifts
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
peak     = 0
centroid = 1
gausspeak = 0

;Systematic Redshift from Lyman alpha redshift 
z  = science.zspec
minlya = (1.-1000/3.e5)*restlya*(z+1)
maxlya = (1.+1000/3.e5)*restlya*(z+1)
set_plot,'ps'
psname=objname+'_redshift.eps'
device, filename = psname,xsize = 15,ysize = 15, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
;set_plot,'x'
!p.font = 0
!p.charsize = 1
goodlya = where(science.lambda gt minlya and science.lambda lt maxlya,cgoodlya)
yrange = minmax(science.contdiv(goodlya))
plot,science.lambda,science.contdiv,yrange=yrange,xrange=[minlya,maxlya],xstyle=1
;plotlines,science.zspec
if cgoodlya gt 0 then begin
   case 1 of
      peak: begin
                                ;peak location
         peakflux = max(science.contdiv[goodlya],peakloc)
         wllya  = (science.lambda[goodlya])(peakloc) ;peak location 
         zsys_lya   = (wllya/restlya+shiftlya/clight)-1. 
      end
      centroid: begin
         wllya    = total(science.contdiv[goodlya]*science.lambda[goodlya])/total(science.contdiv[goodlya])
         zsys_lya = wllya/restlya-1.+shiftlya/clight 
      end
      gausspeak: begin 
                                ;fit gaussian to the red side
         if objname eq 'ms1358' then minlya = wllya
         goodlya =  where(science.lambda gt minlya and science.lambda lt maxlya,cgoodlya)
         lambnow = science.lambda(goodlya)
         specnow = science.contdiv(goodlya)
         param   = [6,wllya,20.,0]
         gauss   = gaussfit(lambnow,specnow,param,sigma=sigma,measure_errors=sigmanow,nterms=4,estimates=param)
         wllya   = param(1)
         oplot,lambnow,gauss,color=fsc_color('pink')
         zsys_lya   = (wllya/restlya+shiftlya/clight)-1.
      end
   endcase

   plotlines,zsys_lya,color='blue',delv=35
   if tag_exist(science,'zsys_lya') then science.zsys_lya = zsys_lya else science=create_struct(science,'zsys_lya',zsys_lya)
endif else stop,'Stop at line 91: Cannot find Lymanalpha redshift'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Systematic Redshift from absorption lines
;SiII 1260
minsi = (1.-1000/3.e5)*1260.42*(z+1)
maxsi = (1.+1000/3.e5)*1264.73*(z+1)
;CII1334
minsi = (1.-1000/3.e5)*1332*(z+1)
maxsi = (1.-1000/3.e5)*1336*(z+1)

plot,science.lambda,science.contdiv,yrange=[0,2],xrange=[minsi,maxsi],xstyle=1
plotlines,zsys_lya,color='blue',delv=35
if tag_exist(science,'minSiabs') then minline = science.minSiabs else begin
   read,minline,prompt='enter minimum region'
   science=create_struct(science,'minSiabs',minline)
endelse
if tag_exist(science,'maxSiabs') then maxline = science.maxSiabs else begin
   read,maxline,prompt='enter maximum region'
   science=create_struct(science,'maxSiabs',maxline)
endelse

goodlamb = where(science.lambda gt minline and science.lambda lt maxline,cgood)
if cgood gt 0 then begin
   lambnow = science.lambda(goodlamb)
   specnow = 1.-science.contdiv(goodlamb)
   ;wlISM  = total(specnow*lambnow)/total(specnow)
   sigmanow= 1./sqrt((science.cont(goodlamb))^2*science.ivar(goodlamb))
   ;param   = [0.5,7460.,20.,-0.4] ms1358
   param   = [0.5,0.5*(minline+maxline),20.,-0.4]
   gauss   = gaussfit(lambnow,specnow,param,sigma=sigma,measure_errors=sigmanow,nterms=4,estimates=param)
   wlISM   = param(1)
   zsys_ISM = (wlISM/1260.42+shiftISM/clight)-1.
   delv_ISM = sqrt(70^2+(sigma(1)/param(1)*3e5)^2)
   oplot,lambnow,1.-gauss,color=fsc_color('dark green')
   plotlines,zsys_ISM,color='dark green',delv=delv_ISM
   al_Legend,['Lya Redshift','SiII1260 Absorption Redshift'],psym=[0,0],color=[fsc_color('Blue'),fsc_color('dark green')],box=0,thick=5,linestyle=2,charsize=0,/right,/bottom
   if tag_exist(science,'zsys_ISM') then science.zsys_ISM = zsys_ISM else science=create_struct(science,'zsys_ISM',zsys_ISM)
endif

device,/close
set_plot,'x'
mwrfits,science,fscience,/create,/silent
;-------------------------------------------------------------------------
;-------------------------------------------------------------------------
;-------------------------------------------------------------------------
;COVERING FRACTION
;Average Low-Ionization Absortion profile
if use_zlya eq 1 then z = science.zsys_lya
if use_zis  eq 1 then z = science.zsys_ism 
if use_zlya eq 1 and use_zis eq 1 then stop,'You have to choose type of systemic redshift'
if use_zlya eq 0 and use_zis eq 0 then stop,'You have to choose type of systemic redshift'

set_plot,'ps'
psname=objname+'_covfrac'+strjoin(strtrim(string(lineuse),2))+'.eps'
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
   contnow    = science.cont(goodwl)
   ivarnow    = contnow^2*science.ivar(goodwl)
   snnow      = sqrt(ivarnow)
   plot,lambdanow,contdivnow,title=linelist(i),ytitle='Normalized Spectrum'
   oplot,lambdanow,snnow/5.,color=fsc_color('cyan')
   plotlines,science.zspec
   ;draw a cross
   if i eq 0 then  al_Legend,['Spec','~SN/5.'],psym=[0,0],color=[fsc_color('black'),fsc_color('cyan')],box=0,thick=2,charsize=0,symsize=1.5,/right,/bottom

   if lineuse(i) eq 0 then begin
      oplot,!x.crange,!y.crange,color=fsc_color('dark red'),thick=2
      oplot,!x.crange,reverse(!y.crange),color=fsc_color('dark red'),thick=2
   endif
endfor
 mwrfits,science,fscience,/create,/silent
; read,lineuse,prompt='input good lines 1=y [SiII1260,OI1302,SiII1304,CII1334,SiII1526] in [1,0,0,1,0]:'
;make velocity step in mean fwhm of dlam
 goodlines = linewl(where(lineuse eq 1))*(science.zsys_lya+1)
 goodlines_name = linelist(where(lineuse eq 1))
 nlines = n_Elements(goodlines)
 stepv  = mean(science.dlam)/mean(science.lambda)*clight*2.355
 nstepv = round(1500/stepv)
 velarr = dblarr(nstepv)
 absarr = dblarr(nstepv)
 absarr_ivar = dblarr(nstepv)
 absarr_wivar = dblarr(nstepv)
 Si1260       = dblarr(nstepv)
 Si1260_ivar  = dblarr(nstepv)
 Si1260_wivar  = dblarr(nstepv)
 wlSiII = 1260.42*(science.zsys_lya+1.)
 velsi = (science.lambda-wlSiII)/wlSiII*clight
 vminsi1304 = -200. ;km
 ls = 1304.37
 lo = 1302.168
 vmaxOi1302 = vminsi1304*ls/lo+clight*(ls-lo)/lo
 for nv=0,nstepv-1 do begin
    contdivnow = []
    ivarnow    = []
    contnow    = []
    for nl=0,nlines-1 do begin
       vnowmin = -1000.+nv*stepv
       vnowmax = -1000.+(nv+1.)*stepv
       velarr(nv) = 0.5*(vnowmin+vnowmax)
       if goodlines_name(nl) eq 'SiII1304' then vnowmin = -1000.+nv*stepv > vminsi1304
       if goodlines_name(nl) eq 'OI1302' then vnowmax = -1000.+(nv+1.)*stepv < vmaxOi1302 
       wlnow  = (science.lambda-goodlines(nl))/goodlines(nl)*clight
       goodwl = where(wlnow gt vnowmin and wlnow lt vnowmax,cwl)
       if cwl gt 0 then begin
          contdivnow = [contdivnow,science.contdiv(goodwl)]
          ivarnow    = [ivarnow,science.ivar(goodwl)]
          contnow    = [contnow,science.cont(goodwl)]
       endif else print, 'no wl found for',goodlines_name(nl),'at velocity ',strtrim(string(velarr(nv)),2)
    endfor
    if n_Elements(contdivnow) ne 0 then begin
       absarr(nv) = total(contdivnow*ivarnow)/total(ivarnow)
       absarr_ivar(nv) = total(contnow^2*ivarnow) ;error of the mean
       absarr_wivar(nv) = total(ivarnow)/total(ivarnow*(contdivnow-absarr(nv))^2);weighted standard deviation       
    endif
    vnowmin = -1000.+nv*stepv
    vnowmax = -1000.+(nv+1.)*stepv
    goodsi = where(velsi gt vnowmin and velsi lt vnowmax,csi)
    si1260(nv) = total(science.contdiv(goodsi)*science.ivar(goodsi))/total(science.ivar(goodsi))
    si1260_ivar(nv) = total((science.cont(goodsi))^2*science.ivar(goodsi))
    si1260_wivar(nv) = total(science.ivar(goodsi))/total(science.ivar(goodsi)*(science.contdiv(goodsi)-si1260(nv))^2)
 endfor
 ;make finer grids for shaded uncertainties
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
 if tag_exist(science,'ave_covfrac') then begin
    science.ave_covfrac=returnave
    science.ave_covfrac_si=returnave_si
 endif else science=create_struct(science,'ave_covfrac',returnave,'ave_covfrac_si',returnave_si)
 mwrfits,science,fscience,/create,/silent

end
