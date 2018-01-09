pro telluric_enk, science,plotstop
  ;TELLURIC CORRECTION WITH TELLURIC CURVE FROM EVAN

  ;Decide which file to use based on zspec
  if science.zspec lt 4.5 then tellfile = '/scr2/nichal/workspace2/telluric/deimos_telluric_1.0.fits'$
             else tellfile = '/scr2/nichal/workspace2/telluric/deimos_telluric_1.0_1200G_8300.fits'
  tell = mrdfits(tellfile, 1, /silent) 
  wtell = n_elements(tell)-1
  if wtell gt 0 then tell = tell[wtell]
  
  specresfile = '/scr2/nichal/workspace/SCIENCE/'+science.objname+'_specres_poly.sav'
  if file_test(specresfile) then begin
     restore, specresfile
     dlam = poly(science.lambda/1000 - 7.8, specres_poly) / 2.35
  endif
  
  aratio = science.airmass/(tell.airmass)
  telllambda = tell.lambda
  tellspec = (tell.spec)^aratio
  tellivar = (tell.ivar)*((tell.spec)/(aratio*tellspec))^2.
  ivarmissing = 10d10
  w = where(tellivar ge 1d8, c)
  if c gt 0 then tellivar[w] = ivarmissing
  f = where(finite(tellspec) and tellspec gt 0 and finite(tellivar) and tellivar gt 0 and tellivar lt 1d8)
  telllambda = telllambda[f]
  tellspec = tellspec[f]
  tellivarfull = tellivar
  tellivar = tellivar[f]

  tellspecnew = interpolate(tellspec, findex(telllambda, science.lambda), missing=1.)
  tellivarnew = interpolate(tellivarfull, findex(telllambda, science.lambda), missing=ivarmissing)
  science.telldiv = science.spec / tellspecnew
  science.telldivivar = (science.ivar^(-1.) + (science.telldiv)^2.*tellivarnew^(-1.))^(-1.) * tellspecnew^2.
  set_plot,'x'
  plot,science.lambda,science.telldiv,title=science.spec1dfile
  oplot,science.lambda,science.spec,color=fsc_color('pink')
  science.telltype = 'ENK:'+tellfile
  if plotstop then stop
end

