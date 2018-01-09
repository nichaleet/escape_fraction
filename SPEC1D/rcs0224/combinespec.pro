pro combinespec
for i=1,4 do begin
  a11=mrdfits('spec1d.rcs0224.007.arc1.013016_try2.matched.fits',i,hdr)
  a12=mrdfits('spec1d.rcs0224.007.arc1.013116_try2.fits',i)
  a21=mrdfits('spec1d.rcs0224.007.arc2.013016_try2.matched.fits',i)
  a22=mrdfits('spec1d.rcs0224.007.arc2.013116_try2.matched.fits',i)
  a31=mrdfits('spec1d.rcs0224.008.arc3.013016_try2.matched.fits',i)
  a32=mrdfits('spec1d.rcs0224.008.arc3.013116_try2.matched.fits',i)
  maxflux = max([median(a11.spec),median(a12.spec),median(a21.spec),median(a22.spec),median(a31.spec),median(a32.spec)])
  a11.spec = a11.spec*maxflux/median(a11.spec)
  a12.spec = a12.spec*maxflux/median(a12.spec)
  a21.spec = a21.spec*maxflux/median(a21.spec)
  a22.spec = a22.spec*maxflux/median(a22.spec)
  a31.spec = a31.spec*maxflux/median(a31.spec)
  a32.spec = a32.spec*maxflux/median(a32.spec)
  spec = fltarr(n_elements(a11.spec))
  ivar = spec
  for j=0,n_Elements(spec)-1 do begin
     farr = [a11.spec[j],a12.spec[j],a21.spec[j],a22.spec[j],a31.spec[j],a32.spec[j]]
     fivar =[a11.ivar[j],a12.ivar[j],a21.ivar[j],a22.ivar[j],a31.ivar[j],a32.ivar[j]]
     good = where(fivar ne 0.,cgood)
     if cgood eq 0 then begin
        spec[j] = 0.
        ivar[j] = 0.
     endif else begin
        meanerr,farr(good),1./sqrt(fivar(good)),fmean,sigmam,sigmad,sigmas
        spec[j] = fmean
        if cgood gt 1 then ivar[j] = 1./sigmas^2 else ivar[j] = fivar(good)
     endelse
  endfor
  comb = a11
  comb.spec = spec
  comb.ivar = ivar
  if i eq 1 then mwrfits,comb,'spec1d.rcs0224.combine.fits',hdr,/create
  if i gt 1 then mwrfits,comb,'spec1d.rcs0224.combine.fits',hdr
endfor
end
