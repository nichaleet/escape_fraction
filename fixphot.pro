function fixphot, fluxarcs, photmag, photmagerr, photmagwl, name, namestar
clight = 3.d18 ;angstrom per second
strout = []
for i = 0,n_Elements(fluxarcs)-1 do begin
   sci = fluxarcs(i)
   lambda = sci.lambda
   flux   = sci.flux
   fluxerr= sci.fluxerr

   centralwl = mean(photmagwl)
   realflux = 10.^((photmag(i)+48.6)/(-2.5)) ;fnu
   realflux = clight*realflux/centralwl^2 ;f lambda
   realfluxerr = realflux*abs(0.4*alog(10.)*photmagerr(i))

   ;calculate average flux in the range
   infilter = where(lambda gt photmagwl[0] and lambda lt photmagwl[1] and finite(fluxerr), cinfilter)
   if cinfilter eq 0 then stop,'PHOTOMETRY WL DOES NOT MATCH THE OBSERVED WL'
   fluxmeasure = wmean(flux(infilter),fluxerr(infilter),error=fluxmeasure_err)
   measuredmag = -2.5*alog10(centralwl^2*fluxmeasure/clight)-48.6
   measuredmagerr = abs(2.5*fluxmeasure_err/fluxmeasure/alog(10.))
   print,'real magnitude is ',photmag(i),photmagerr(i),' measured to be ', measuredmag,measuredmagerr

   ;calculate the shift
   multfactor = realflux/fluxmeasure
   multfactorerr = abs(multfactor*sqrt((realfluxerr/realflux)^2+(fluxmeasure_err/fluxmeasure)^2))
   
   newflux    = flux*multfactor
   newfluxerr = sqrt(multfactorerr^2*flux^2+multfactor*fluxerr^2)
   newmag     = -2.5*alog10(lambda^2*newflux/clight)-48.6
   newmagerr  = abs(2.5*newfluxerr/newflux/alog(10.))
   str =  {lambda:lambda,flux:newflux,fluxerr:newfluxerr,mag:newmag,magerr:newmagerr}
   strout = [strout,str]
   set_plot,'x'
   !p.multi=[0,1,2]
   plot,lambda,newflux,ytitle='flux'
   plot,lambda,newmag,ytitle='mag'
   stop
endfor
return,strout
end
