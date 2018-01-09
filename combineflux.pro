function combineflux, fluxarcs,plot=plot
clight = 3.d18 ;angstrom per second

  lambda = fluxarcs.lambda
  if total(lambda[*,0]-lambda[*,1]) ne 0 then stop, 'lambda are diff'
  lambda = lambda[*,0]
  newflux = fltarr(n_Elements(lambda))
  newfluxerr = newflux
  flux = fluxarcs.flux
  fluxerr = fluxarcs.fluxerr
  for i=0,n_elements(lambda)-1 do begin
     newflux[i]=wmean(flux[i,*],fluxerr[i,*],error=ferr)
     newfluxerr[i] = ferr
  endfor
  newmag     = -2.5*alog10(lambda^2*newflux/clight)-48.6
  newmagerr  = abs(2.5*newfluxerr/newflux/alog(10.))
  str =  {lambda:lambda,flux:newflux,fluxerr:newfluxerr,mag:newmag,magerr:newmagerr}
  if keyword_Set(plot) then begin
     set_plot,'x'
     !p.multi=[0,1,1]
     colorlist = ['purple','cyan','green','orange','red']
     plot,lambda,flux[*,0],title='combineflux.pro'
     for i=0,n_elements(flux[0,*])-1 do oplot,lambda,flux[*,i],color=fsc_color(colorlist[i])
     oplot,lambda,newflux
     wait,2
  endif
  return,str
end
