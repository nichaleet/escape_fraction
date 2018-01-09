pro match_lambda,fmain,fsub,fmainout,fsubout
;Purpose: Interpolate the second mask to match the wavelengths of the main mask 
;Inputs : fmain and fsub are the spec1dfiles (strings). The files must be in /scr2/nichal/workspace/SPEC1D
;Outputs: Strings of the names of the new spec1d files in the same folder.

dir = '/scr2/nichal/workspace3/SPEC1D/'
fmainout = dir+strmid(fmain,0,strpos(fmain,'.fits'))+'.matched.fits'
fsubout  = dir+strmid(fsub,0,strpos(fsub,'.fits'))+'.matched.fits'

;loop over blue and red side
for i=1,2 do begin
   main  = mrdfits(fmain,i,hdrmain)
   sub   = mrdfits(fsub,i,hdrsub)
   minwl = max([min(main.lambda),min(sub.lambda)])
   maxwl = min([max(main.lambda),max(sub.lambda)])
 ;trim main slit
   posmin = min(where(main.lambda gt minwl))>0
   posmax = max(where(main.lambda lt maxwl))
 ;find tags to trim and trim
 ;trim only if the number of n_Elements(tag) =  n_elements(lambda)
   tnames = tag_names(main)
   trim_arr = bytarr(n_Elements(tnames))
   for t=0,n_elements(tnames)-1 do begin
      if size(main.(t),/dimensions) eq n_elements(main.lambda) then trim_arr(t) = 1
      if t eq 0 then begin
         if trim_arr(t) eq 0 then tag_now = main.(t) else tag_now=(main.(t))[posmin:posmax]
         tmp=create_struct(tnames[t],tag_now)
      endif else begin
         if trim_arr(t) eq 0 then tag_now = main.(t) else tag_now=(main.(t))[posmin:posmax]
         tmp=create_struct(tmp,tnames[t],tag_now)
      endelse
   endfor
   main = tmp
 ;finished trimming now main has a new size for lambda
 
 ;interpolate sub to match main 
   spec_int=interpol(sub.lambda,sub.spec,main.lambda,/spline)
   ivar_int=interpol(sub.lambda,sub.ivar,main.lambda,/spline)
   skyspec_int = interpol(sub.lambda,sub.skyspec,main.lambda,/spline)

   sub = create_struct('SPEC',spec_int,'LAMBDA',main.lambda,'IVAR',ivar_int,'skyspec',skyspec_int)

   for j = 3,n_elements(tnames)-1 do result=create_struct(result,tnames[j],main.(j))
   ;fix header
   expmain = sxpar(hdrmain,'exptime',comment=comments)
   expsub = sxpar(hdrsub,'exptime')
   sxaddpar,hdrmain,'exptime',expmain+expsub,comments
   mskname1 = sxpar(hdrmain,'slmsknam')
   mskname2 = sxpar(hdrsub,'slmsknam')
   sxaddpar,hdrmain,'slmsknam',mskname1+'&'+mskname2
   slitno1 = sxpar(hdrmain,'slitno',comment=comments)
   slitno2 = sxpar(hdrsub,'slitno')
   sxaddpar,hdrmain,'slitno',strtrim(string(slitno1),2)+'&'+strtrim(string(slitno2),2),comments
;write output
   if i eq 1 then mwrfits,result,fout,hdrmain,/create else mwrfits,result,fout,hdrmain

endfor

;plot
window,0,xsize=700,ysize=900
!p.multi = [0,1,3]
setplot,14
restore,'line.sav'
z=4.93
file =[fout,fmain,fsub]
title=['All data','Main slit data','Sub data']
for f=0,2 do begin
   blu = mrdfits(file(f),1,hdrblu)
   red = mrdfits(file(f),2,hdrred)
   linewave = line_wl*(z+1.)
   linename = line_name
   nlines= n_elements(line_wl)
   lambda = [blu.lambda,red.lambda]
   spec=[blu.skyspec,red.skyspec]
   ;if f eq 0 then spec   =[blu.skydiv,red.skydiv] else spec=[blu.spec/blu.skyspec,red.spec/red.skyspec]
   plot,lambda,spec,yrange=[-100,1000],xrange=minmax(lambda),xstyle=1,title=title[f],font=0,ystyle=1
   for j=0,nlines-1 do begin    ;plot markers for emission lines
      oplot,[linewave[j],linewave[j]],[(!y.crange)[0],((!y.crange)[0])*.25],color=55,linestyle=2
      xyouts,linewave[j],((!y.crange)[0])/2.,linename[j]
   endfor
endfor
stop
end
return,newsub
end
