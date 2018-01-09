function cal_flux_conversion,phot,segment,imflux

  nobjs = n_elements(phot.id)
  flux  = fltarr(nobjs)
  markgood = bytarr(nobjs)
  for i=0,nobjs-1 do begin
     id = phot.id(i)
     pix = where(segment eq id,cpix)
     if cpix gt 0 then begin
        markgood(i) = 1
        flux(i)     = total(imflux(pix))
     endif
  endfor
  ok = where(markgood eq 1, cok)
  print,'From', nobjs, 'objects, find', cok,'goodmatch'

  magarr = phot.magauto(ok)
  magerr = phot.magerr(ok)
  fluxcgs = 10.^((magarr+48.6)/(-2.5))/1.e-29;fnu/1.e-29
  flux   = flux(ok)
  
  param = linfit(flux,fluxcgs,sigma=sigma,yfit=yfit)
  conversion = [param[1]*1.e-29,sigma[1]*1.e-29]
  print, 'params A, B are', param
  print, 'conversion is', conversion ,'fnu/count'
  set_plot,'x'
  !p.multi=[0,1,1]
  plot,flux,fluxcgs,psym=1,xtitle='HST Counts',ytitle='fnu(cgs)'
  oplot,flux,yfit,color=fsc_color('orange')
  return, conversion

end
