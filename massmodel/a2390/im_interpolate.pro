function im_interpolate,img
  imgnew = img
  fluxpix = where(img gt 0,cflux)
  xyfluxind = array_indices(img,fluxpix)
  xmin = min(xyfluxind(0,*))
  xmax = max(xyfluxind(0,*))
  ymin = min(xyfluxind(1,*))
  ymax = max(xyfluxind(1,*))

  nofluxpix = where(img eq 0,cnoflux)
  xynofluxind = array_indices(img,nofluxpix)
  xind = xynofluxind(0,*)
  yind = xynofluxind(1,*)
  xind = reform(xind)
  yind = reform(yind)

  tobeint = where(xind ge xmin and xind le xmax and yind ge ymin and yind le ymax,ctobeint)
  xind = xind(tobeint)
  yind = yind(tobeint)
  countchange = 0.
  for jj=0, n_elements(xind)-1 do begin
     x = xind(jj)
     y = yind(jj)
     surroundpix = [img(x-1,y-1),img(x-1,y),img(x-1,y+1),img(x,y-1),img(x,y+1),img(x+1,y-1),img(x+1,y),img(x+1,y+1)]
     goodpix = where(surroundpix gt 0.,goodpixcount)
     if goodpixcount ge 5. then begin
        imgnew(x,y) = mean(surroundpix(goodpix))
        countchange = countchange+1.
     endif
  endfor
  print, 'npix changed=',countchange
return,imgnew
end
