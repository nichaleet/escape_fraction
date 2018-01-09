pro fitcontinuum2,science
  fft = 0
  spline = 1
  poly = 0
  usesmooth = 0

  contmask = science.contmask
  w = where(contmask ne 0, c)
  lambda = science.lambda / (1d + science.zspec)
  n = n_elements(lambda)
  nhalf = round(double(n) / 2.)

  satbands = [6870, 7650]
  niter = 5

  if n lt 8193 then begin
     wwhole = lindgen(nhalf)
     ccd2 = 2
  endif else begin
     wwhole = lindgen(n)
     ccd2 = 1
  endelse
  spec = science.contdiv
  for ccd=1,ccd2 do begin
     won = where(contmask[wwhole] eq 1, complement=woff, con) + wwhole[0]
     woff += wwhole[0]
     
     case 1 of
        spline: begin
           lower = 1.5
           upper = 5
           bkpt = slatec_splinefit(lambda[won], spec[won], coeff, invvar=science.telldivivar[won], bkspace=150, upper=upper, lower=lower, /silent)
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
                 cont = replicate(-999, n_elements(science.lambda))
                 return
              endif
              cont = smooth_gauss_wrapper(lambda[wcont], spec[wcont], lambda[wwhole], 5., ivar1=science.contdivivar[wcont])
              wcont = where(abs(spec[ww]-cont[ww-wwhole[0]]) lt (science.contdivivar[ww])^(-0.5) and contmask[ww] eq 1, ccont) + ww[0]
              
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

  set_plot,'x'
  !p.multi = [0,1,2]
  plot,science.lambda,science.contdiv,yrange=[-5,5]
  oplot,science.lambda,cont,color=fsc_color('green')
  plot,science.lambda,science.contdiv,yrange=[-5,5]
  oplot,science.lambda,science.contdiv/cont,color=fsc_color('green')
  oplot,science.lambda,science.contmask-5,color=fsc_color('cyan')

  science.contmask = contmask
  science.contdiv = science.contdiv / cont
  science.contdivivar = science.contdivivar * cont^2.

 ;CHECK WHERE CONTINUUM IS NEGATIVE ==> make it zero and make condivivar to 0
  noflux = where(cont lt 0. or science.spec eq -99,cnoflux)
  if cnoflux gt 0 then begin
     science.contmask(noflux) = 0
     science.contdiv(noflux) = 0.
     science.contdivivar(noflux)= 0.
  endif

end
	
