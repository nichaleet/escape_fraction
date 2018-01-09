pro prepspec,file1ds,objname,zspec=zspec,smoothfactor=smoothfactor,stardiv=stardiv,enkdiv=enkdiv,deep2div=deep2div,twicecontdiv=twicecontdiv,plot=plot
  
;PURPOSE: telluric correction, continuum normalization, and combine specra

;OPTIONAL INPUTS
;Telluric type
;stardiv: use observed star spectrum as telluric correction. The star must be in the /SPEC1D folder with the name of objname_telluric.fits
;enkdiv: This is DEFAULT. Use the code from Evan. It uses the workspace2/telluric/deimos_telluric_1.0.fits Feige110 star to do telluric correction
;deep2div: This is from Tucker. It uses deep2's pipeline

;OUTPUT:
;fscience    ='/scr2/nichal/workspace3/SCIENCE/'+objname+'_science.fits'
;fscienceall ='/scr2/nichal/workspace3/SCIENCE/'+objname+'_scienceall.fits'
;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;
;Constants
linelist =['SiII1260','OI1302','SiII1304','CII1334','SiII1526']
linewl   = [1260.42,1302.168,1304.37,1334.532,1526.72]
clight = 3.e5
restlya  = 1215.67
zoommin = 1220*(zspec+1)  ;SiII
zoommax = 1283*(zspec+1)
zoommin2 = 1500*(zspec+1) ;SiII
zoommax2 = 1600*(zspec+1)
fullmin = 1200*(zspec+1)
fullmax = 1560*(zspec+1)
if ~keyword_set(smoothfactor) then smoothfactor = 1.
nmasks = n_Elements(file1ds)
fscience = '/scr2/nichal/workspace3/SCIENCE/'+objname+'_science.fits'
fscienceall = '/scr2/nichal/workspace3/SCIENCE/'+objname+'_scienceall.fits'

;check keywords
if n_elements(enkdiv) then enkdiv = 1 else enkdiv = 0
if n_Elements(deep2div) then deep2div = 1 else deep2div = 0
if n_Elements(stardiv) then stardiv = 1 else stardiv = 0 
if (deep2div or stardiv) eq 0 then enkdiv = 1
if n_elements(plot) then plotstop = 1 else plotstop = 0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
scienceall = []
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
      nline = n_Elements(linelist)
      science = {objname:objname, $
                 spec1dfile:file1d, $
                 objpos:[-999d], $
                 fwhm:[-999d], $
                 exptime:[-999d],$
                 lambda:dblarr(npix), $
                 spec:dblarr(npix), $
                 smoothed:dblarr(npix), $
                 ivar:dblarr(npix), $
                 sn:dblarr(npix), $
                 skyspec:dblarr(npix), $
                 telldiv:dblarr(npix), $
                 telldivivar:dblarr(npix), $
                 telltype:'telltype',$
                 contmask:bytarr(npix), $
                 cont:dblarr(npix), $
                 contdiv:dblarr(npix), $
                 contdivivar:dblarr(npix), $
                 fitmask:bytarr(npix)+1, $
                 dlam:dblarr(npix), $
                 skyfit:[[dblarr(nsky)], [dblarr(nsky)], [dblarr(nsky)]], $
                 skylinemask:lonarr(nsky),$
                 goodsky:0, $
                 sn_all:-999d, $
                 airmass:-999d, $
                 zspec:-999d, $
                 zsys_lya:-999d,$
                 zsys_ism_arr:dblarr(nline),$
                 zsys_ism_ave:-999d,$
                 zsys_ism_std:-999d,$
                 method:'method',$
                 linelist:linelist,$
                 linewl:linewl,$
                 smoothfactor:smoothfactor}
      science.objpos = blu.objpos
      science.fwhm   = blu.fwhm
      science.exptime= sxpar(hdrblue,'exptime')
      science.lambda = [blu.lambda,red.lambda]
      science.spec   = [blu.spec,red.spec]
      if objname eq 'a2219' then science.spec = science.spec+20.;make the lyc=0
      science.ivar   = [blu.ivar,red.ivar]
      science.skyspec= [blu.skyspec,red.skyspec]
      science.zspec  = zspec
      science.airmass = sxpar(hdrblue,'AIRMASS')

      ;SPECTRAL RESOLUTION
      science.skylinemask = -1
      specres,science,qf=specres_poly
      save,specres_poly,filename='/scr2/nichal/workspace3/SCIENCE/'+objname+'_specres_poly.sav'
     
      ;TELLURIC CORRECTION
      case 1 of 
         stardiv: begin
            ;telluric with observed star
            tellfile = file_Search('/scr2/nichal/workspace/Telluric/'+objname+'_telluric.fits',count=ctell)
            if ctell ne 0 then tellfile=tellfile(0) else stop, 'Stop: No telluric file found.'      
            telluric,science,tellfile
         end
         enkdiv: begin
            ;telluric with Evan's telluric file
            telluric_enk,science,plotstop
            
         end
         deep2div: begin
            ;telluric with deep2
            specmerge = fill_gap(file1d, head=hdr)
            tell_corrected = qe_correction2(specmerge, hdr, /telluric, ncorrect=7)
            science.telldiv = tell_corrected
            science.telldivivar = science.ivar
            plot,science.lambda,science.telldiv
            science.telltype = 'DEEP2'
         end
      endcase

      ;SMOOTH
      ;smooth with smoothfactor*spectral resolution
      ;science.smoothed = median_smooth(science.lambda,science.telldiv,science.lambda,2.*science.dlam)
      science.smoothed = smooth_gauss_wrapper(science.lambda,science.telldiv,science.lambda,smoothfactor*science.dlam)

      ;PLOT TELLDIV SMOOTH AND SKY
      !p.font = 0
      set_plot,'x'
      setplot,14
      plot,science.lambda,science.spec,xrange=[zoommin,zoommax],yrange=[-300,1000],title=file1d,/nodata,xtitle='CHECK Telldiv+SMOOTH and SKY'
      oplot,science.lambda,science.skyspec,color=fsc_color('pink')
      oplot,science.lambda,science.spec
      oplot,science.lambda,science.smoothed,color=fsc_color('green')
      plotuvlines,science.zspec
      al_Legend,['Original Spec',science.telltype+' Telldiv&Smoothed spectra','sky spectra'],psym=[0,0,0],color=fsc_color(['black','green','pink']),box=0,linestyle=2,/right,/bottom
      if plotstop then stop

      set_plot,'ps'
      !p.font = 1
      psname='output/'+objname+strtrim(string(n),2)+'_telldiv_smooth_sky.eps'
      device, filename = psname,xsize = 25,ysize = 12,/encapsulated,/color
      plot, science.lambda, science.spec,ytitle='SPEC',xrange=[fullmin,fullmax],xtitle='wavelength(A)',yrange=[-500,1000]
      oplot,science.lambda,science.skyspec,color=fsc_color('pink')
      oplot,science.lambda,science.spec
      oplot,science.lambda,science.smoothed,color=fsc_color('green')
      al_Legend,['Original Spec',science.telltype+' Telldiv&Smoothed spectra','sky spectra'],psym=[0,0,0],color=fsc_color(['black','green','pink']),box=0,linestyle=2,/right,/bottom
      plotuvlines,science.zspec
      device,/close
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      science.telldiv = science.smoothed ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ;FIT CONTINUUM AND DO CONTINUUM NORMALIZATION
      fitcontinuum, science

      set_plot,'x'
      window,0,xsize=1200,ysize=400
      !p.multi=[0,1,2]
      plot,science.lambda,science.telldiv,xrange=[fullmin,fullmax],yrange=[-300,300],title='CHECK CONTINUUM'
      oplot,science.lambda,science.contmask*50-300.,color=fsc_color('cyan'),thick=2
      oplot,science.lambda,science.cont,color=50
      plotuvlines,science.zspec
      plot,science.lambda,science.contdiv,xrange=[fullmin,fullmax],yrange = [-10,10]
      oplot,!x.crange,[1,1],color=fsc_color('cyan'),linestyle=2
  
      linewl   = [1260.42,1302.168,1304.37,1334.532,1526.72]
      window,1,xsize=1200,ysize=1000
      !p.multi=[0,1,5]
      for jj=0,4 do begin
         zoomlinemin = (1.-1500/3.e5)*linewl(jj)*(science.zspec+1)
         zoomlinemax = (1.+1500/3.e5)*linewl(jj)*(science.zspec+1)
         toplot = where(science.lambda lt zoomlinemax and science.lambda gt zoomlinemin)
         plot,science.lambda,science.telldiv,xrange=[zoomlinemin,zoomlinemax],/nodata,title=science.spec1dfile
         skytoplot = science.skyspec(toplot)
         skytoplot = skytoplot/max(skytoplot)*(max(!y.crange)-min(!y.crange))
         oplot,science.lambda(toplot),skytoplot,color=fsc_color('pink')
         oplot,science.lambda,science.telldiv
         oplot,science.lambda,science.cont,color=fsc_color('green')
         oplot,science.lambda,science.contmask*50.+!y.crange[0],color=fsc_color('cyan')
         plotuvlines,science.zspec
      endfor
      stop
      !p.multi = [0,1,1]
      window,0
      scienceall = [scienceall,science]
 
      ;find signal to noise for each lambda in a spectral resolution
      sn = calsn(science,smoothfactor)
      science.sn = sn
      
      ;PLOT 
      set_plot,'ps'
      psname='output/'+objname+strtrim(string(n),2)+'_signal2noise.eps'
      device, filename = psname,xsize = 12,ysize = 20, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      !p.multi = [0,1,4]
      plot, science.lambda,science.telldiv*sqrt(science.telldivivar),ytitle='telldiv*sqrt(telldivivar)'
      plot, science.lambda,sn,ytitle='wmean/wstdev'
      plot, science.lambda,science.telldiv,yrange=[-500,500],xrange=[fullmin,fullmax],ytitle='telldiv'
      oplot,science.lambda,science.spec,color=fsc_color('pink')
      plotuvlines,science.zspec
      plot, science.lambda,science.dlam*smoothfactor,xrange=[fullmin,fullmax],ytitle='Angstrom'
      oplot,science.lambda,science.dlam,linestyle=2
      xyouts,[0.92*mean(!x.crange)],[mean(science.dlam)*smoothfactor],'smoothing range used for SN'
      xyouts,[mean(!x.crange)],[mean(science.dlam)],'spectral resolution'
      device,/close
      !p.multi = [0,1,1]

   endfor
  
          
   ;combine each frames
   if nmasks gt 1 then begin
      nlamb = n_Elements(scienceall[0].lambda)
      contdiv = dblarr(nlamb,nmasks)
      ivar    = dblarr(nlamb,nmasks)
      dlam    = dblarr(nlamb,nmasks)
      for i=0,nmasks-1 do begin
         sciencenow = scienceall[i]
         contdiv[*,i] = sciencenow.contdiv
         ivar[*,i]    = sciencenow.contdivivar
         dlam[*,i]    = sciencenow.dlam
      endfor
      s_ivar = total(ivar,2)
      fin_dlam = mean(dlam,dimension=2)
      fin_contdiv = total(contdiv*ivar,2)/s_ivar
      sn_cont = calsn_cont(sciencenow.lambda,fin_contdiv,s_ivar,fin_dlam*smoothfactor)
      result={contdiv:fin_contdiv,lambda:science.lambda,ivar:s_ivar,dlam:fin_dlam,zspec:science.zspec}
      science = {objname:strjoin(scienceall.objname,' ; '), $
                 spec1dfile:strjoin(scienceall.spec1dfile,' ; '), $
                 objpos:scienceall.objpos, $
                 fwhm:scienceall.fwhm, $
                 exptime:scienceall.exptime,$
                 lambda:sciencenow.lambda, $
                 spec:scienceall.spec, $
                 smoothed:scienceall.smoothed, $
                 ivar:scienceall.ivar, $
                 skyspec:scienceall.skyspec, $
                 telldiv:scienceall.telldiv, $
                 telldivivar:scienceall.telldivivar, $
                 contmask:scienceall.contmask, $
                 fitmask:total(scienceall.fitmask,2), $
                 skylinemask:scienceall.skylinemask, $
                 contdiv:fin_contdiv, $
                 contdivivar:s_ivar, $
                 dlam:fin_dlam, $
                 skyfit:scienceall.skyfit, $
                 airmass:scienceall.airmass, $
                 zspec:science.zspec ,$
                 zsys_lya:-999d,$
                 zsys_ism_arr:dblarr(nline),$
                 zsys_ism_ave:-999d,$
                 zsys_ism_std:-999d,$
                 method:'method',$
                 linelist:linelist,$
                 lineuse:[1,1,1,1,1], $
                 linewl:linewl,$
                 sn:sn_cont, $
                 smoothfactor:smoothfactor}
      if n_elements(twicecontdiv) eq 1 then fitcontinuum2,science
      set_plot,'x'
      !p.multi=[0,1,3]
      setplot,14
      plot,science.lambda,science.contdiv,xrange=[fullmin,fullmax],title='combined contdiv',yrange=[-5,5]
      oplot,scienceall[0].lambda,scienceall[0].contdiv,color=fsc_color('red')
      oplot,science.lambda,science.contdiv,color=fsc_color('black')
      plotuvlines,science.zspec
      plot,science.lambda/(1.+science.zspec),science.contdiv,title='combined contdiv',yrange=[-5,50]
      colors = fsc_color(['purple','cyan','green','yellow','orange','red','pink'])
      for jj=1,nmasks do begin
         cur = scienceall[nmasks-jj]
         oplot, cur.lambda/(1.+science.zspec),cur.contdiv,color=colors(jj)
      endfor
      plotuvlines, 0.,color='black'
      plot,science.lambda,fin_dlam,ytitle='mean spectral resolution'
   endif else scienceall=science
   mwrfits,science,fscience,/create,/silent
   mwrfits,scienceall,fscienceall,/create,/silent
   stop  
;-------------------------------------------------------------------------
;FINISH PREPARING SCIENCE 
;-------------------------------------------------------------------------
end
