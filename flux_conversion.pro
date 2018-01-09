pro flux_conversion,starblu,starred,reflam,refmag,namegal,namestar
  flux_conversionall=[]
  flambdaall=[]
  contall=[]
  for j=0,1 do begin
     if j eq 0 then star = starblu
     if j eq 1 then star = starred
     lambda = star.lambda
     spec   = star.spec
  ;interpolate the reference to obs_lambda
     refmag_match = interpol(refmag,reflam,lambda)
     refmag_match = double(refmag_match)
  ;convert AB mag to cgs flux
     clight = 3.d18                     ;angstrom/s
     fnu = 10.^(-(refmag_match+48.6)/2.5d) ;erc/s/cm^2/hz
     flambda = fnu*clight/lambda^2         ;erg/s/cm^2/Angstrom
  
  ;;;fit continuum to the observed flux with the telluric masked out
  ;masking
     readcol,'/scr2/nichal/workspace3/telluric.mask', tellstart, tellend, format='D,D', /silent, comment='#'
     tellmask = bytarr(n_elements(lambda))+1
     for i=0,n_elements(tellstart)-1 do begin
        woff = where(lambda gt tellstart(i) and lambda lt tellend(i),cwoff)
        if cwoff gt 0 then tellmask(woff) = 0
     endfor
     won = where(tellmask)
  ;fit continuum
     bkpt = slatec_splinefit(lambda[won], spec[won], coeff, invvar=star.ivar[won], bkspace=150, upper=10, lower=0.05, /silent)
     if bkpt[0] eq -1 then stop,'CANNOT FIT CONTINUUM'
     cont = slatec_bvalu(lambda,bkpt,coeff)

  ;find the conversion
     flux_conversionall = [flux_conversionall,flambda/cont]
     flambdaall = [flambdaall,flambda]
     contall    = [contall,cont]
  endfor
     
  ;save the conversion
  lambda = [starblu.lambda,starred.lambda]
  spec   = [starblu.spec,starred.spec]
  conversion = {lambda:lambda,conv:flux_conversionall,unit:'flux(cgs)/Deimos flux',star:namestar,gal:namegal}
  mwrfits,conversion,'/scr2/nichal/workspace3/flux_calib/'+namegal+'_'+namestar+'.fits',/create
  ;plot
  set_plot,'ps'
  !p.font = 0 
  !P.multi = [0,1,3]
  !p.charsize = 1.5
  psname='/scr2/nichal/workspace3/flux_calib/'+namegal+'_'+namestar+'_fluxconversion.eps'
  xrange = minmax(lambda)
  device, filename = psname,xsize = 15,ysize = 15, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  plot,lambda,spec,xtitle='lambda',ytitle='observed star spectra',title=namegal+', '+namestar,/nodata,xrange=xrange,xstyle=1
  masktell
  oplot,lambda,spec
  oplot,lambda,contall,color=fsc_color('red')
  plot,lambda,flambdaall,xtitle='lambda',ytitle='Real flambda(cgs)',xrange=xrange,xstyle=1
  plot,lambda,flux_conversionall,xtitle='lambda',ytitle='flux conversion (flux(cgs)/deimos flux)',xrange=xrange,xstyle=1
  device,/close
end
