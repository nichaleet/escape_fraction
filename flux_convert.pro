function flux_convert,scienceall,param,namegal,namestar
  clight = 3.d18 ;angstrom/s
  set_plot,'ps'
  !p.font = 0 
  !P.multi = [0,2,n_elements(scienceall)]
  !p.charsize = 1
  psname='/scr2/nichal/workspace3/flux_calib/'+namegal+'_'+namestar+'_magnitude.eps'
  device, filename = psname,xsize = 25,ysize = 15, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  strout = []
  for i=0,n_Elements(scienceall)-1 do begin
     sci = scienceall[i]
     spec = sci.telldiv
     lambda = sci.lambda
     ;interpolate the param lambda to galaxy lambda
     conv = interpol(param.conv,param.lambda,lambda)
     ;calculate flux
     flux  = spec*conv
     fluxerr = 1./sqrt(sci.telldivivar)*conv
     plot,lambda,flux,xtitle='lambda',ytitle='flux (erg/s/cm2/A)'
     oplot,lambda,flux+fluxerr,color=fsc_color('rose')
     oplot,lambda,flux-fluxerr,color=fsc_color('rose')
     oplot,lambda,flux
     ;convert back to magnitude
     fnu    = flux*lambda^2/clight
     fnuerr = fluxerr*lambda^2/clight
     mAB    = -2.5*alog10(fnu)-48.6
     mABerr = 2.5*fnuerr/fnu/alog(10.)
     plot,lambda,mAB,xtitle='lambda',ytitle='AB Magnitude'
     oplot,lambda,mAB+MABerr,color=fsc_color('rose')
     oplot,lambda,mAB-MABerr,color=fsc_color('rose')
     oplot,lambda,mAB
     ;save to structure
     str = {lambda:lambda,flux:flux,fluxerr:fluxerr,mag:mAB,magerr:mABerr}
     strout = [strout,str] 
  endfor
  xyouts,0.5,0.97,namegal,alignment=[0.5],/normal
  device,/close
  return,strout
end
