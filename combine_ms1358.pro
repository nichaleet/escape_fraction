pro combine_ms1358
;Combine the long slit data with the mask data. The main is long slit. The output is fout which is the same as normal spec1d from deep2 pipeline but with an extra tag in the structure called skydiv because we combine the sky_div data from each mask. the skyspec in the output is from the long slit only.

flong = '/scr2/nichal/deimos/reduced_data/ms1358_long/boxspec1d.Long1.0B.007.N3b.fits'
fmask = '/scr2/nichal/deimos/reduced_data/ms1358_try3/boxspec1d.ms1358.009.arcBC.fits'
fout = '/scr2/nichal/deimos/reduced_data/boxspec1d.ms1358.combine.fits'
fmaskout = '/scr2/nichal/workspace/SPEC1D/boxspec1d.ms1358.009.arcBC.matched.fits'
flongout = '/scr2/nichal/workspace/SPEC1D/boxspec1d.Long1.0B.007.N3b.matched.fits'

;loop over blue and red side
for i=1,2 do begin
   long = mrdfits(flong,i,hdrlong)
   mask = mrdfits(fmask,i,hdrmask)
   minwl = max([min(long.lambda),min(mask.lambda)])
   maxwl = min([max(long.lambda),max(mask.lambda)])
 ;trim long slit
   posmin = min(where(long.lambda gt minwl))>0
   posmax = max(where(long.lambda lt maxwl))
 ;find tags to trim and trim
 ;trim only if the number of n_Elements(tag) =  n_elements(lambda)
   tnames = tag_names(long)
   trim_arr = bytarr(n_Elements(tnames))
   for t=0,n_elements(tnames)-1 do begin
      if size(long.(t),/dimensions) eq n_elements(long.lambda) then trim_arr(t) = 1
      if t eq 0 then begin
         if trim_arr(t) eq 0 then tag_now = long.(t) else tag_now=(long.(t))[posmin:posmax]
         tmp=create_struct(tnames[t],tag_now)
      endif else begin
         if trim_arr(t) eq 0 then tag_now = long.(t) else tag_now=(long.(t))[posmin:posmax]
         tmp=create_struct(tmp,tnames[t],tag_now)
      endelse
   endfor
   long = tmp
 ;finished trimming now long has a new size for lambda
 ;make skydiv spectra
   mask_skydiv = mask.spec/mask.skyspec
   long_skydiv = long.spec/long.skyspec

 ;interpolate mask to match long 
   spec_int=spline(mask.lambda,mask.spec,long.lambda,/double)
   ivar_int=spline(mask.lambda,mask.ivar,long.lambda,/double)
   skydiv_int = spline(mask.lambda,mask_skydiv,long.lambda,/double)
   skyspec_int = spline(mask.lambda,mask.skyspec,long.lambda,/double)
   masknew = create_struct('SPEC',spec_int,'LAMBDA',long.lambda,'IVAR',ivar_int,'skydiv',skydiv_int,'OBJPOS',mask.objpos,'FWHM',mask.fwhm,'NSIGMA',mask.nsigma,'R1',mask.R1,'R2',mask.R2,'SKYSPEC',skyspec_int,'IVARFUDGE',mask.ivarfudge)
   
 ;write output for mask and long
   if i eq 1 then begin
      mwrfits,long,flongout,hdrlong,/create
      mwrfits,masknew,fmaskout,hdrmask,/create 
   endif else begin
      mwrfits,long,flongout,hdrlong
      mwrfits,masknew,fmaskout,hdrmask
   endelse
   mask = masknew
 ;combine spectra with weighted arithmatic mean
   spec = [[long.spec],[mask.spec]]
   skydiv = [[long_skydiv],[mask.skydiv]]
   ivar = [[long.ivar],[mask.ivar]]
   s_ivar = total(ivar,2)
   
   fin_spec = total(spec*ivar,2)/s_ivar
   fin_skydiv = total(skydiv*ivar,2)/s_ivar
   result={spec:fin_spec,lambda:long.lambda,ivar:s_ivar,skydiv:fin_skydiv}
   for j = 3,n_elements(tnames)-1 do result=create_struct(result,tnames[j],long.(j))
   ;fix header
   explong = sxpar(hdrlong,'exptime',comment=comments)
   expmask = sxpar(hdrmask,'exptime')
   sxaddpar,hdrlong,'exptime',explong+expmask,comments
   mskname1 = sxpar(hdrlong,'slmsknam')
   mskname2 = sxpar(hdrmask,'slmsknam')
   sxaddpar,hdrlong,'slmsknam',mskname1+'&'+mskname2
   slitno1 = sxpar(hdrlong,'slitno',comment=comments)
   slitno2 = sxpar(hdrmask,'slitno')
   sxaddpar,hdrlong,'slitno',strtrim(string(slitno1),2)+'&'+strtrim(string(slitno2),2),comments
;write output
   if i eq 1 then mwrfits,result,fout,hdrlong,/create else mwrfits,result,fout,hdrlong

endfor

;plot
window,0,xsize=700,ysize=900
!p.multi = [0,1,3]
setplot,14
restore,'line.sav'
z=4.93
file =[fout,flong,fmask]
title=['All data','Long slit data','Mask data']
for f=0,2 do begin
   blu = mrdfits(file(f),1,hdrblu)
   red = mrdfits(file(f),2,hdrred)
   linewave = line_wl*(z+1.)
   linename = line_name
   nlines= n_elements(line_wl)
   lambda = [blu.lambda,red.lambda]
   ;spec=[blu.skyspec,red.skyspec]
   if f eq 0 then spec   =[blu.skydiv,red.skydiv] else spec=[blu.spec/blu.skyspec,red.spec/red.skyspec]
   plot,lambda,spec,yrange=[-10,30],xrange=minmax(lambda),xstyle=1,title=title[f],font=0,ystyle=1
   for j=0,nlines-1 do begin    ;plot markers for emission lines
      oplot,[linewave[j],linewave[j]],[(!y.crange)[0],((!y.crange)[0])*.25],color=55,linestyle=2
      xyouts,linewave[j],((!y.crange)[0])/2.,linename[j]
   endfor
endfor
stop
end
