pro telluric, science,tellfile
  tell = mrdfits(tellfile,1,hdr)
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
  plot,science.lambda,science.telldiv
  science.telltype = 'OBS STAR:'+tellfile
  stop
end

