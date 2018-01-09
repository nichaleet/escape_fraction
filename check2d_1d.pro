pro check2d_1d,directory,slit,objname,z
;e.g. IDL> check2d_1d,'ms1358_long/oneobj','007','N3b',4.93
  setplot,14
  restore,'line.sav'
  linewave = line_wl*(z+1.)
  linename = line_name
  nlines= n_elements(line_wl)
;  parentdir = '/scr2/nichal/deimos/reduced_data/'
;  parentdir = '/scr2/nichal/deimos/rawdata/2015may/'
  parentdir = '/scr2/nichal/deimos/rawdata/2016jan/'
  directory = parentdir+directory
  file_2d_blue = file_search(directory+'/slit.*.'+slit+'B.{fits,fits.gz}', count=c)
  if c eq 0 then stop, 'No files found'
  twodblue = mrdfits(file_2d_blue[0],1)

  file_2d_red = file_search(directory+'/slit.*.'+slit+'R.{fits,fits.gz}', count=c)
  if c eq 0 then stop, 'No files found'
  twodred = mrdfits(file_2d_red[0],1)

  file1d = file_search(directory+'/*spec1d.*.'+slit+'.'+objname+'.{fits,fits.gz}', count=c)
  if c eq 0 then stop, 'No files found'

  for i=0,c-1 do begin
     print, file1d[i]
     onedblue =mrdfits(file1d[i],1,hdr)
     onedred  =mrdfits(file1d[i],2)
     ;ymin = max([-1000,min(onedblue.spec)])
     ymin = min(onedblue.spec)
     ymin =-1000
     ;ymax = min([1000,max(onedblue.spec)])
     ymax = max(onedblue.spec)
     ymax = 1000
     ;window,i
     plot,[onedblue.lambda,onedred.lambda],[onedblue.spec,onedred.spec],yrange=[ymin,ymax],xrange=[7400,8500]
     for j=0,nlines-1 do begin ;plot markers for emission lines
        oplot,[linewave[j],linewave[j]],[(!y.crange)[0],((!y.crange)[0])/2.],color=fsc_color('red'),linestyle=2
        xyouts,linewave[j],((!y.crange)[0])/2.,linename[j]
     endfor
     print, 'obj position', onedblue.objpos
     print, 'obj fwhm',onedblue.fwhm

     bflux = twodblue.flux
     rflux = twodred.flux
     width = min([(size(bflux,/dimensions))[1],(size(rflux, /dimensions))[1]])
     flux = [bflux[*,0:width-1],rflux[*,0:width-1]]
     atv,flux
     ;atv,[twodred.flux]
  endfor
end
