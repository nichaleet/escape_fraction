pro tv_slit2d, directory, slit
  ;parentdir = '/scr2/nichal/deimos/rawdata/2015mar/'
  parentdir = '/scr2/nichal/deimos/reduced_data/'
  directory = parentdir+directory
  file_2d_blue = file_search(directory+'/slit.*.'+slit+'B.{fits,fits.gz}', count=c)
  if c eq 0 then stop, 'No files found'
  bslit = mrdfits(file_2d_blue[0],1)
  file_2d_red = file_search(directory+'/slit.*.'+slit+'R.{fits,fits.gz}', count=c)
  if c eq 0 then stop, 'No files found'
  
  loadct, 0, /silent
  window, 0, xsize=1600,ysize=300
  tv,dblarr(1600,300),0,0
  
  rslit = mrdfits(file_2d_red[0],1)  
  bflux = bslit.flux
  rflux = rslit.flux
  width = min([(size(bflux,/dimensions))[1],(size(rflux, /dimensions))[1]])
  flux = [bflux[*,0:width-1],rflux[*,0:width-1]]

  med = median(flux)
  sig = stddev(flux)
  max_value = (med + (10 * sig)) < max(flux)
  min_value = (med - (2 * sig))  > min(flux)
  if (finite(min_value) EQ 0) then min_value = image_min
  if (finite(max_value) EQ 0) then max_value = image_max
  if (min_value GE max_value) then begin
     min_value = min_value - 1
     max_value = max_value + 1
  endif

  bflux = bytscl(bslit.flux, /nan, min=min_value, max=max_value, top=255) + 8
  bflux = -1.0*bflux + 255.0
  bheight = 1600./2047.*(size(bslit.dlambda, /dimensions))[1]
  if bheight gt 74 then begin
     bwidth = round(74./bheight*1600.)
     bheight = 74
  endif else begin
     bwidth = 1600
     bheight = round(bheight)
  endelse
  bflux1 = congrid(bflux[0:2047,*], bwidth, bheight, /interp)
  bflux2 = congrid(bflux[2048:4095,*], bwidth, bheight, /interp)
  tv, bflux1, 0, 225
  tv, bflux2, 0, 150

  rflux = bytscl(rslit.flux, /nan, min=min_value, max=max_value, top=247) + 8
  rflux = -1.0*rflux + 255.0
  rheight = 1600./2047.*(size(rslit.dlambda, /dimensions))[1]
  if rheight gt 74 then begin
     rwidth = round(74./rheight*1600.)
     rheight = 74
  endif else begin
     rwidth = 1600
     rheight = round(rheight)
  endelse
  rflux1 = congrid(rflux[0:2047,*], rwidth, rheight, /interp)
  rflux2 = congrid(rflux[2048:4095,*], rwidth, rheight, /interp)
  tv, rflux1, 0, 75
  tv, rflux2, 0, 0
  rgood = 1
  atv,flux
end
