pro fitcontinuum,science
  fft = 0
  spline = 1
  poly = 0
  usesmooth = 0

  contmask = science.contmask
  w = where(contmask ne 0, c)
  lambda = science.lambda / (1d + science.zspec)
  n = n_elements(lambda)
  nhalf = round(double(n) / 2.)
  ;MASKING
  if c eq 0 then begin
;     velout = -700. ;km/s ;TO BE TWEAK
;     velin  = 300.  ;km/s ;TO BE TWEAK
     velout = -250. ;km/s ;TO BE TWEAK
     velin  = 200.  ;km/s ;TO BE TWEAK
     readcol, '/scr2/nichal/workspace3/lines.txt',linewave,fosc,linename,format ='D,D,A,X',comment='#'
     linestart= (linewave*velout/3.e5+linewave)
     lineend  = (linewave*velin/3.e5+linewave)
     contmask = bytarr(n_elements(lambda))+1
     for i=0,n_elements(linestart)-1 do begin
        w = where(lambda ge linestart[i] and lambda le lineend[i], c)
        if c gt 0 then contmask[w] = 0
     endfor
     ;make extra wide mask for Lya
     w = where(lambda gt -1215.*500./3.e5+1215. and lambda lt 1215.*1500./3.e5+1215)
     contmask[w]=0
     ;telluric mask
     tellmask = bytarr(n_elements(lambda))
     tellstart = [6864., 7591., 8938.]
     tellend = [6935., 7694., 9000.] 
     for i=0,n_elements(tellstart)-1 do begin
        w = where(science.lambda ge tellstart[i] and science.lambda le tellend[i], c)
        if c gt 0 then begin
    ;       contmask[w] = 0   ;do not mask
           tellmask[w] = 1
        endif
     endfor
     contmask[0:2] = 0
     contmask[nhalf-3:nhalf+2] = 0
     contmask[n-3:n-1] = 0
     ;if science.objname eq 'ms1358' then begin
     ;   badskyline = where(science.skylinemask eq 1 and science.skyfit[*,0] gt 8300 and science.skyfit[*,0] lt 8600,cbadskyline)
     ;   meanwidth = mean(science.dlam)*2.35
     ;   for csky=0,cbadskyline-1 do begin
     ;      badlam = science.skyfit[badskyline[csky],0]
     ;      w = where(science.lambda ge badlam-meanwidth and science.lambda le badlam+meanwidth, c)
     ;      if c gt 0 then begin
     ;         contmask[w] = 0
     ;         tellmask[w] = 1
     ;      endif
     ;   endfor
     ;endif
  endif

  satbands = [6870, 7650]
  niter = 5

  if n lt 8193 then begin
     wwhole = lindgen(nhalf)
     ccd2 = 2
  endif else begin
     wwhole = lindgen(n)
     ccd2 = 1
  endelse
  spec = science.telldiv
  for ccd=1,ccd2 do begin
     contmask[wwhole[0:3]] = 0
     contmask[wwhole[nhalf-4:nhalf-1]] = 0
     won = where(contmask[wwhole] eq 1, complement=woff, con) + wwhole[0]
     woff += wwhole[0]
     
     case 1 of
        spline: begin
           lower = 1.5          ;1.5
           upper = 10.
           bkpt = slatec_splinefit(lambda[won], spec[won], coeff, invvar=science.telldivivar[won], bkspace=5, upper=upper, lower=lower, /silent)
           if bkpt[0] eq -1 then return
           cont = slatec_bvalu(lambda[wwhole], bkpt, coeff)
        end
        poly: begin
           degree = 12
           norm = median(spec[won])
           a = [norm, replicate(0.0, degree-1)]
           p = lmfit(lambda[won], spec[won], a, measure_errors=(science.telldivivar[won])^(-0.5), /double, function_name='legendre_poly')
           cont = legendre_poly(lambda[wwhole], a, /noderiv)
        end
        fft: begin
           wfft = wwhole[10:nhalf-11]
           nfft = n_elements(wfft)
           nkeep = 10
           ft = fft(spec[wfft], -1, /double)
           
           ft[nkeep:nfft-1] = 0
           cont = real_part(fft(ft, 1, /double))
           plot, cont
           stop
        end
        usesmooth: begin
           ww = wwhole[3:nhalf-4]
           wcont = where(contmask[ww] eq 1, ccont) + ww[0]
           for i=1,5 do begin
              wcontold = wcont
              if ccont lt 50 then begin
                 science.cont = replicate(-999, n_elements(science.cont))
                 return
              endif
              cont = smooth_gauss_wrapper(lambda[wcont], spec[wcont], lambda[wwhole], 30.0, ivar1=science.telldivivar[wcont])
              wcont = where(abs(spec[ww]-cont[ww-wwhole[0]]) lt (science.telldivivar[ww])^(-0.5) and contmask[ww] eq 1, ccont) + ww[0]
              
              if array_equal(wcont, wcontold) then break
           endfor
        end
     endcase

     if ccd eq 1 then contb = cont
     if ccd eq 2 then contr = cont        
     wwhole += nhalf        
  endfor
  
  if n lt 8193 then cont = [contb, contr] else cont = contb
  if n_elements(cont) gt n then cont = cont[0:n-1]

  science.contmask = contmask
  science.cont = cont
  science.contdiv = science.telldiv / cont
  science.contdivivar = science.telldivivar * cont^2.

 ;CHECK WHERE CONTINUUM IS NEGATIVE ==> make it zero and make condivivar to 0
  noflux = where(cont lt 0. or science.spec eq -99,cnoflux)
  if cnoflux gt 0 then begin
     science.contmask(noflux) = 0
     science.cont(noflux) = 0
     science.contdiv(noflux) = 0.
     science.contdivivar(noflux)= 0.
  endif
  
  wcont = where(contmask[3:n-4] eq 1)+3
  wcont = wcont[where(finite(science.telldiv[wcont]) and finite(science.cont[wcont]) and science.cont[wcont] ne 0)]
  dev = abs((science.telldiv[wcont] - science.cont[wcont]) / science.cont[wcont])
  avgdev = mean(dev)
  w = where(dev lt 3.0*avgdev, c)
  if c gt 0 then science.sn_all = 1.0/mean(dev[w])
stop
end
	
