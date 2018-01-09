pro match_lambda,fmain,fsub,z=z
;Purpose: Interpolate the second mask to match the wavelengths of the main mask 
;Inputs : fmain and fsub are the spec1dfiles (strings). The files must be in /scr2/nichal/workspace/SPEC1D
;Outputs: Strings of the names of the new spec1d files in the same folder.

dir = '/scr2/nichal/workspace3/SPEC1D';/rcs0224/'
fsubout  = dir+strmid(fsub,0,strpos(fsub,'.fits'))+'.matched.fits'
fmain = dir+fmain
fsub  = dir+fsub


;loop over blue and red side

for i=1,4 do begin
   ;Read and merge sub data
   sub   = mrdfits(fsub,i,hdrsub)
   if i eq 1 or i eq 2 then sub = fill_gap(fsub,header=hdr,/boxsprof)
   if i eq 3 or i eq 4 then sub = fill_gap(fsub,header=hdr,/horne)
   ;readmain
   main  = mrdfits(fmain,i,hdrmain)
   minwl = min(main.lambda)
   maxwl = max(main.lambda)
   nlamb = n_elements(main.lambda)
   ;interpolate sub to match main 
   spec_int=interpol(sub.spec,sub.lambda,main.lambda,/spline)
   ivar_int=interpol(sub.ivar,sub.lambda,main.lambda,/quadratic)
   skyspec_int = interpol(sub.skyspec,sub.lambda,main.lambda,/quadratic)
   crmask_int = interpol(sub.crmask,sub.lambda,main.lambda,/quadratic)
   bitmask_int = interpol(sub.bitmask,sub.lambda,main.lambda,/quadratic)
   ormask_int = interpol(sub.ormask,sub.lambda,main.lambda,/quadratic)
   nbadpix_int = interpol(sub.nbadpix,sub.lambda,main.lambda,/quadratic)
   infomask_int = interpol(sub.infomask,sub.lambda,main.lambda,/quadratic)
   ;mask the region where it's extrapolate for the sub mask
   masksub  = bytarr(nlamb)+1
   badsubwl = where(main.lambda lt min(sub.lambda) or main.lambda gt max(sub.lambda), cbadsub)
   if cbadsub gt 0 then begin
      masksub(badsubwl) = 0 
      ivar_int(badsubwl) = 0.
      spec_int(badsubwl) = -99.
   endif
 ;create new structures
   newsub = create_struct('SPEC',spec_int,'LAMBDA',main.lambda,'IVAR',ivar_int,'skyspec',skyspec_int,'MASK_INTPOL',masksub,'CRMASK',crmask_int,'BITMASK',bitmask_int,'ORMASK',ormask_int,'NBADPIX',nbadpix_int,'INFOMASK',infomask_int,'OBJPOS',sub.objpos,'FWHM',sub.fwhm,'NSIGMA',sub.nsigma,'R1',sub.r1,'R2',sub.r2,'IVARFUDGE',sub.ivarfudge)
   if i eq 1 then mwrfits,newsub,fsubout,hdrsub,/create else mwrfits,newsub,fsubout,hdrsub

endfor

;plot
window,0,xsize=700,ysize=900
!p.multi = [0,1,3]
setplot,14
restore,'line.sav'
file =[fmain,fsubout,fsub]
title=['Main slit data','New Sub data','Sub data']
for i=0,1 do begin
   for f=0,2 do begin
      blu = mrdfits(file(f),(i*2)+1,hdrblu)
      red = mrdfits(file(f),(i*2)+2,hdrred)
      linewave = line_wl*(z+1.)
      linename = line_name
      nlines= n_elements(line_wl)
      lambda = [blu.lambda,red.lambda]
      spec=[blu.skyspec,red.skyspec]
      spec=[blu.spec,red.spec]
      plot,lambda,spec,yrange=[-100,2000],title=title[f],font=0,ystyle=1
      vline,linewave,color=fsc_color('red')
      for j=0,nlines-1 do begin ;plot markers for emission lines
         xyouts,linewave[j],((!y.crange)[1])/2.,linename[j]
      endfor
   endfor
stop
endfor
end

