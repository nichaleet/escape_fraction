pro fitcontinuum_tell,tfile,name=name,spline=spline,poly=poly,horne=horne,boxsprof=boxsprof

;FIT RED AND BLUE SEPARATELY CAUSE THE SPECTRUM MIGHT HAVE A DISCONTIINUOUS JUMP AND THAT MESS THE FITTING
;INPUT: 
;       tfile = spec1d from the spec2d reduction pipeline
;       name  = name of the star. see the output
;OPTIONAL INPUT
;       /Spline: Spline fit is the default
;       /Poly  : Polynomial fitting
;       /horne : Use the star spectrum extracted from horne method --> DEFAULT
;       /boxprof: USE the star spectrum extracted with box car

;OUTPUT: is a telluric file with the name of 'name_telluric.fits'. It's a continuum normoallized for the telluric features of the star.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Check Parameters
  if n_Elements(poly) ge 1 then poly =1 else poly = 0
  if n_Elements(spline) ge 1 then spline = 1 else spline = 0
  if poly eq 0 then spline = 1

  if n_Elements(boxsprof) ne 0 then initibr = 1 else initibr = 3
  if n_elements(horne) ne 0 then initibr=3

  fout = '/scr2/nichal/workspace3/Telluric/'+name+'_telluric.fits'
  tfile = '/scr2/nichal/workspace3/SPEC1D/'+tfile

  niter = 4

  ;Loop Over Blue and Red side
  for ibr = initibr,initibr+1 do begin  
     ;Loop for ITERATIONS (Keep blocking out the absorption lines (those with cont level less than 0.9 out)
     For iit = 0,niter do begin
        ;Read Data
        tfile_exist = file_search(tfile,count=cf)
        if cf eq 0 then STOP, 'NO FILE '+tfile+' FOUND.'
        tell = mrdfits(tfile,ibr,hdr)
     
        ;masking the expected telluric regions 
        n        = n_elements(tell.lambda)
        tellmask = bytarr(n)+1

        velout = -700.          ;km/s
        velin  = 300.           ;km/s

        tellstart = [6864., 7591., 8938.]
        tellend = [6935., 7694., 9000.] 
        for i=0,n_elements(tellstart)-1 do begin
           w = where(tell.lambda ge tellstart[i] and tell.lambda le tellend[i], c)
           if c gt 0 then tellmask[w] = 0
        endfor

        ;mask the first/last few pixels
        tellmask[0:2] = 0
        tellmask[n-3:n-1] = 0

        ;mask the extra tell absorption lines
        tellmask2 = tellmask
        if iit gt 0 then begin
           w = where(tellold lt 0.9,c)
           if c gt 0 then tellmask2[w]=0
        endif

        spec = tell.spec
        lambda = tell.lambda
        won = where(tellmask eq 1 and tellmask2 eq 1, complement=woff, con) 
 
        ;fitting
        case 1 of
           spline: begin
              bkpt = slatec_splinefit(lambda[won], spec[won], coeff, invvar=tell.ivar[won], bkspace=150, upper=5, lower=1.5, /silent)
              if bkpt[0] eq -1 then return
              cont = slatec_bvalu(lambda, bkpt, coeff)
           end
           poly: begin
              degree = 12
              norm = median(spec[won])
              a = [norm, replicate(0.0, degree-1)]
              p = lmfit(lambda[won], spec[won], a, measure_errors=(tell.ivar[won])^(-0.5), /double, function_name='legendre_poly')
              cont = legendre_poly(lambda, a, /noderiv)
           end
        endcase
        if ibr eq 1 or ibr eq 3 then begin
           specb = spec
           contb = cont
           lambb = lambda
           tellb = spec/cont
           tellbivar = tell.ivar*cont^2
           skyb      = tell.skyspec
           tellmaskb = tellmask
        endif
        if ibr eq 2 or ibr eq 4 then begin
           specr = spec
           contr = cont  
           lambr = lambda
           tellr = spec/cont
           tellrivar = tell.ivar*cont^2
           skyr      = tell.skyspec
           tellmaskr = tellmask
        endif
        tellold=spec/cont
        !p.multi = [0,1,1]
        plot,lambda,spec
        oplot,lambda,cont,color=fsc_color('red')
        wait,0.5
     endfor
  endfor

  ;MAKE OUTPUT STRUCTURE 
  spec   = [specb,specr]
  cont   = [contb, contr] 
  lambda = [lambb,lambr]
  tell   = [tellb,tellr]
  tellivar= [tellbivar,tellrivar]
  sky     = [skyb,skyr]
  tellmask = [tellmaskb,tellmaskr]
  object  = sxpar(hdr,'OBJECT')
  airmass = sxpar(hdr,'AIRMASS')
  exptime = sxpar(hdr,'EXPTIME')
  date    = sxpar(hdr,'DATE-OBS')
  UTC     = sxpar(hdr,'UTC')
  RA      = sxpar(hdr,'RA')
  DEC     = sxpar(hdr,'DEC')
  
  dev = abs((spec-cont) / cont)
  avgdev = mean(dev)
  w = where(dev lt 3.0*avgdev, c)
  if c gt 0 then sn = 1.0/mean(dev[w])
  print, 'SIGNAL TO NOISE: ', sn

  ;MAKE THE TELLURIC FREE REGIONS TO BE 1
  ;Read in the template from Evan's telluric curve
  telltemp1 = mrdfits('/scr2/nichal/workspace2/telluric/deimos_telluric_1.0_1200G_8300.fits',1,hdrtemp)
  good = where(telltemp1.spec ne -1.)
  tempspec = interpol(telltemp1.spec(good),telltemp1.lambda(good),lambda,/spline)
  won = where(tempspec eq 1. and lambda gt min(telltemp1.lambda(good)) and lambda lt max(telltemp1.lambda(good)) and tellmask ne 0,cwon)
  if cwon gt 0 then tell(won) = 1.

  ;FIX SOME PASCHEN LINES. USE THE TEMPLATE
  wpaschen_start = [8980,9200]
  wpaschen_end   = [9040,9260]
  for i=0,n_Elements(wpaschen_start)-1 do begin
     wpaschen = where(lambda gt wpaschen_start[i] and lambda lt wpaschen_end[i],cpas)
     if cpas gt 0 then tell(wpaschen)= tempspec(wpaschen)
  endfor

  ;SAVE FILE
  tellstr = {OBJECT:object,NAME:name,AIRMASS:airmass,EXPTIME:exptime,DATE:date,UTC:UTC,RA:RA,DEC:DEC,LAMBDA:lambda,SPEC:tell,IVAR:tellivar}
  mwrfits,tellstr,fout,hdr
  
  ;PLOT
  !p.font = 0
  !p.multi=[0,1,3]
  spec1d = fill_gap(tfile)
  plot,spec1d.lambda,spec1d.spec
  oplot,lambda,cont,color=fsc_color('red')
  plot,tellstr.lambda,tellstr.spec
  plot,tellstr.lambda,tellstr.spec
  oplot,telltemp1.lambda,telltemp1.spec,color=fsc_color('yellow')
  stop

end
	
