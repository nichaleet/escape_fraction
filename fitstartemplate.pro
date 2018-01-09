pro fitstartemplate,tfile,startemplate,name=name,horne=horne,boxsprof=boxsprof

;FIT RED AND BLUE SEPARATELY CAUSE THE SPECTRUM MIGHT HAVE A DISCONTIINUOUS JUMP AND THAT MESS THE FITTING
;INPUT: 
;       tfile = spec1d from the spec2d reduction pipeline
;       name  = name of the star. see the output
;       startemplate = name of the template to be fitted in the folder 
;                      /scr2/nichal/UVKLIBspectral_lib/
;                      Should be corresponding to the star type
;OPTIONAL INPUT
;       /horne : Use the star spectrum extracted from horne method --> DEFAULT
;       /boxprof: USE the star spectrum extracted with box car

;OUTPUT: is a telluric file with the name of 'name_telluric.fits'. It's a continuum normoallized for the telluric features of the star.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Check Parameters

  if n_Elements(boxsprof) ne 0 then initibr = 1 else initibr = 3
  if n_elements(horne) ne 0 then initibr=3

  fout = '/scr2/nichal/workspace3/Telluric/'+name+'_telluric.fits'
  tfile = '/scr2/nichal/workspace3/SPEC1D/'+tfile


  ;Loop Over Blue and Red side
  for ibr = initibr+1,initibr+1 do begin  
     ;Read Data
     tfile_exist = file_search(tfile,count=cf)
     if cf eq 0 then STOP, 'NO FILE '+tfile+' FOUND.'
     tell = mrdfits(tfile,ibr,hdr)
     
     ;masking the expected telluric regions 
     n        = n_elements(tell.lambda)
     tellmask = bytarr(n)+1
     
     velout = -700.             ;km/s
     velin  = 300.              ;km/s
     
     tellstart = [6864., 7591., 8938.]
     tellend = [6935., 7694., 9000.] 
     for i=0,n_elements(tellstart)-1 do begin
        w = where(tell.lambda ge tellstart[i] and tell.lambda le tellend[i], c)
        if c gt 0 then tellmask[w] = 0
     endfor

     ;mask the first/last few pixels
     tellmask[0:2] = 0
     tellmask[n-3:n-1] = 0

     spec = tell.spec
     lambda = tell.lambda
     won = where(tellmask eq 1, complement=woff, con) 
 
     ;Fitting star template
        ;guess the continuum level
        ;Read in star template & interpolate the lambda to match with the DEIMOS
     starfile = '/scr2/nichal/UVKLIBspectral_lib/'+startemplate
     readcol,starfile,slamb,sspec,v3,v4,v5,v6,v7
     starspec = interpol(sspec,slamb,lambda,/spline)
     contlevel = median(spec(won))/median(starspec)
     stop
     if ibr eq 1 or ibr eq 3 then begin
        specb = spec
        contb = cont
        lambb = lambda
        skyb      = tell.skyspec
     endif
     if ibr eq 2 or ibr eq 4 then begin
        specr = spec
        contr = cont  
        lambr = lambda
        skyr      = tell.skyspec
     endif
     tellold=spec/cont
     !p.multi = [0,1,2]
     plot,lambda,spec
     oplot,lambda,cont,color=fsc_color('red')
     
  endfor

  ;MAKE OUTPUT STRUCTURE 
  spec   = [specb,specr]
  cont   = [contb, contr] 
  lambda = [lambb,lambr]
  tell   = [tellb,tellr]
  tellivar= [tellbivar,tellrivar]
  sky     = [skyb,skyr]
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

  ;SAVE FILE
  tellstr = {OBJECT:object,NAME:name,AIRMASS:airmass,EXPTIME:exptime,DATE:date,UTC:UTC,RA:RA,DEC:DEC,LAMBDA:lambda,SPEC:tell,IVAR:tellivar}
  mwrfits,tellstr,fout,hdr
  
  ;PLOT
  !p.multi=[0,1,2]
  spec1d = fill_gap(tfile)
  plot,spec1d.lambda,spec1d.spec
  oplot,lambda,cont,color=fsc_color('red')
  plot,tellstr.lambda,tellstr.spec
  stop

end
	
