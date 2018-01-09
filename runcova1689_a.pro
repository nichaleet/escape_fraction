pro runcova1689_a,redo=redo

;Combine the long slit data with the mask data. The main is long slit. The output is fout which is the same as normal spec1d from deep2 pipeline but with an extra tag in the structure called skydiv because we combine the sky_div data from each mask. the skyspec in the output is from the long slit only.

fa = '/scr2/nichal/deimos/rawdata/2015may/a1689/spec1d.a1689.007.arcz4868a.fits'
fb = '/scr2/nichal/deimos/rawdata/2015may/a1689/spec1d.a1689.012.arcz4868b.fits'
fout = '/scr2/nichal/deimos/reduced_data/spec1d.a1689.007.arcz4868.fits'
faout = '/scr2/nichal/workspace/SPEC1D/spec1d.a1689.007.arcz4868a.matched.fits'
fbout = '/scr2/nichal/workspace/SPEC1D/spec1d.a1689.012.arcz4868b.matched.fits'
;if ~file_test(faout) and ~file_test(fbout) or keyword_set(redo) then begin
if ~file_test(faout) and ~file_test(fbout) then begin
;loop over blue and red side
   for i=1,2 do begin
      a = mrdfits(fa,i,hdra)
      b = mrdfits(fb,i,hdrb)
      minwl = max([min(a.lambda),min(b.lambda)])
      maxwl = min([max(a.lambda),max(b.lambda)])
 ;trim a slit
      posmin = min(where(a.lambda gt minwl))>0
      posmax = max(where(a.lambda lt maxwl))
 ;find tags to trim and trim
 ;trim only if the number of n_Elements(tag) =  n_elements(lambda)
      tnames = tag_names(a)
      trim_arr =  bytarr(n_Elements(tnames))
      for t=0,n_elements(tnames)-1 do begin
         if size(a.(t),/dimensions) eq n_elements(a.lambda) then trim_arr(t) = 1
         if t eq 0 then begin
            if trim_arr(t) eq 0 then tag_now = a.(t) else tag_now=(a.(t))[posmin:posmax]
            tmp=create_struct(tnames[t],tag_now)
         endif else begin
            if trim_arr(t) eq 0 then tag_now = a.(t) else tag_now=(a.(t))[posmin:posmax]
            tmp=create_struct(tmp,tnames[t],tag_now)
         endelse
      endfor
      a = tmp
      
 ;interpolate b to match a 
      spec_int=spline(b.lambda,b.spec,a.lambda,/double)
      ivar_int=spline(b.lambda,b.ivar,a.lambda,/double)
      skyspec_int = spline(b.lambda,b.skyspec,a.lambda,/double)
      bnew = create_struct('SPEC',spec_int,'LAMBDA',a.lambda,'IVAR',ivar_int,'OBJPOS',b.objpos,'FWHM',b.fwhm,'NSIGMA',b.nsigma,'R1',b.R1,'R2',b.R2,'SKYSPEC',skyspec_int,'IVARFUDGE',b.ivarfudge)
      
 ;write output for b and a
      if i eq 1 then begin
         mwrfits,a,faout,hdra,/create
         mwrfits,bnew,fbout,hdrb,/create 
      endif else begin
         mwrfits,a,faout,hdra
         mwrfits,bnew,fbout,hdrb
      endelse
      b = bnew
 ;combine spectra with weighted arithmatic mean
      spec = [[a.spec],[b.spec]]
      ivar = [[a.ivar],[b.ivar]]
      s_ivar = total(ivar,2)
   
      fin_spec = total(spec*ivar,2)/s_ivar
      result={spec:fin_spec,lambda:a.lambda,ivar:s_ivar}
      for j = 3,n_elements(tnames)-1 do result=create_struct(result,tnames[j],a.(j))
   ;fix header
      expa = sxpar(hdra,'exptime',comment=comments)
      expb = sxpar(hdrb,'exptime')
      sxaddpar,hdra,'exptime',expa+expb,comments
      mskname1 = sxpar(hdra,'slmsknam')
      mskname2 = sxpar(hdrb,'slmsknam')
      sxaddpar,hdra,'slmsknam',mskname1+'&'+mskname2
      slitno1 = sxpar(hdra,'slitno',comment=comments)
      slitno2 = sxpar(hdrb,'slitno')
      sxaddpar,hdra,'slitno',strtrim(string(slitno1),2)+'&'+strtrim(string(slitno2),2),comments
;write output
      if i eq 1 then mwrfits,result,fout,hdra,/create else mwrfits,result,fout,hdra

   endfor

;plot
   window,0,xsize=700,ysize=900
   !p.multi = [0,1,3]
   setplot,14
   restore,'line.sav'
   z=4.868
   file =[fout,fa,fb]
   title=['All data','A slit data','B data']
   for f=0,2 do begin
      blu = mrdfits(file(f),1,hdrblu)
      red = mrdfits(file(f),2,hdrred)
      linewave = line_wl*(z+1.)
      linename = line_name
      nlines= n_elements(line_wl)
      lambda = [blu.lambda,red.lambda]
   ;spec=[blu.skyspec,red.skyspec]
      spec=[blu.spec,red.spec]
      skyspec=[blu.skyspec,red.skyspec]
      plot,lambda,spec,xrange=[min(lambda),7500],xstyle=1,title=title[f],font=0,ystyle=1,/nodata
      oplot,lambda,skyspec,color=fsc_color('pink')
      oplot,lambda,spec
      
      for j=0,nlines-1 do begin ;plot markers for emission lines
         oplot,[linewave[j],linewave[j]],[(!y.crange)[0],((!y.crange)[0])*.25],color=55,linestyle=2
         xyouts,linewave[j],((!y.crange)[0])/2.,linename[j]
      endfor
   endfor
stop
endif

;window,0,xsize=1100,ysize=700

covfrac,[faout,fbout],'a1689a',returnave,returnave_si,zspec=4.868,smoothfactor=1.,redo=redo
; returnave=[[velarr],[absarr],[absarr_wivar]]
; returnave_si = [[velarr],[si1260],[si1260_wivar]]

covfrac_si,'a1689a',velarr,fc_mcmc,clmd_mcmc,fc_grid,clmd_grid
; fc_grid = size[4,nstepvel] in peak,first quartile,median, third quartiles

set_plot,'ps'
psname='a1689a_covfrac_all.eps'
device, filename = psname,xsize = 14,ysize = 10, $
             xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
plot,returnave[*,0],returnave[*,1],psym=10,title=objname,ytitle='Average Absorbtion profile',xtitle='velocity(km/s)',/nodata,yrange=[0,1.5],color=fsc_color('black'),font=0,charsize=1
shaded_uncertainties,returnave[*,0],returnave[*,1],returnave[*,2],color='rose'
shaded_uncertainties_ul,velarr,1.-fc_grid[0,*],1.-fc_grid[1,*],1.-fc_grid[3,*],color='grey'
oplot,returnave[*,0],returnave[*,1],psym=10,color=fsc_color('red')
oplot,velarr,1.-fc_grid[0,*],psym=10,linestyle=1,color=fsc_color('black')
al_Legend,['Average line profiles','Si II covering fraction'],psym=[15,15],color=[fsc_color('red'),fsc_color('black')],linestyle=[1,2],box=0,thick=2,charsize=1,symsize=1.5,/right,/bottom,font=0
device,/close
stop
end
