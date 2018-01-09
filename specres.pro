function fitskylines,science
  lambda = science.lambda
  skyspec = science.skyspec
  w = where(~finite(skyspec), c)
  if c gt 0 then skyspec[w] = 0
  skyspec = skyspec/max(skyspec)
  medskyspec = median(skyspec)
  ;plot, lambda, skyspec
  ;oplot, minmax(lambda), 2.0*replicate(medskyspec, 2), linestyle=1
  ;stop
  
  deriv1skyspec = deriv(lambda, skyspec)
  deriv2skyspec = deriv(lambda, deriv1skyspec)
  ;plot, lambda, deriv1skyspec
  ;oplot, minmax(lambda), replicate(0.2, 2), linestyle=1
  ;stop
  ;plot, lambda, deriv2skyspec, yrange=[-0.1, 0.1]
  ;oplot, minmax(lambda), replicate(-0.01, 2), linestyle=1
  ;stop
  thresh = 1.5
  nlines = 1000
  while nlines gt 200 do begin
     w = where(abs(deriv1skyspec) lt 0.2 and deriv2skyspec lt -0.01 and skyspec gt thresh*medskyspec)
     w = [w, n_elements(lambda)+1]
     wstep = round(-1*ts_diff(w, 1))
     linestart = where(wstep gt 1, nlines)
     if nlines lt 5 then begin
        message, 'Fewer than 5 sky lines in this spectrum.', /info
        return, [[-1], [-1], [-1]]
     endif
     linepix = round(-1*ts_diff(linestart, 1))
     nlines -= 1
     thresh *= 2
  endwhile
  linepix = linepix[0:nlines-1]
  wlocmax = lindgen(nlines)
  sigma = dblarr(nlines)
  sigmaerr = dblarr(nlines)
  linelambda = dblarr(nlines)
  nskyspec = n_elements(skyspec)
  wgoodline = bytarr(nlines)+1
  for i=0,nlines-1 do begin        
     if linepix[i] eq 1 then wlocmax[i] = w[linestart[i]] else begin
        junk = min(abs(deriv1skyspec[w[linestart[i]]:w[linestart[i]]+linepix[i]-1]), wmin)
        wlocmax[i] = w[linestart[i]]+wmin ;location of peak
     endelse
     if wlocmax[i]-10 lt 0 or wlocmax[i]+10 ge nskyspec then begin
        wgoodline[i] = 0
        continue
     endif
     skyspecmax = max(skyspec[wlocmax[i]-10:wlocmax[i]+10], wlocmaxtemp)
     wlocmax[i] += wlocmaxtemp-10
     wfit = lindgen(20)+wlocmax[i]-10
     wfit = wfit[where(wfit ge 0 and wfit lt n_elements(lambda), nfit)]
     lambdafit = lambda[wfit]
     skyspecfit = skyspec[wfit]
 
     skyspecmax = max(skyspecfit, wmax)
     
     guess = [skyspecmax, lambdafit[wmax], 0, medskyspec]
     yfit = gaussfit(lambdafit, skyspecfit, a, estimates=guess, sigma=aerr, nterms=4, chisq=chisq)
     sigma[i] = abs(a[2])
     sigmaerr[i] = aerr[2]
     linelambda[i] = a[1]

     ;print, sigma[i], sigmaerr[i]/sigma[i], chisq
     ;plot, lambdafit, skyspecfit
     ;oplot, lambdafit, yfit, color=fsc_color('red')
     ;wait, 0.1

     if chisq gt 1d-3 or sigmaerr[i] gt 0.8 then wgoodline[i] = 0
  endfor
  wgood = where(sigma gt 0 and wgoodline eq 1, c)
  if c gt 0 then begin
     linelambda = linelambda[wgood]
     sigma = sigma[wgood]
     sigmaerr = sigmaerr[wgood]
     ploterr, linelambda, sigma, sigmaerr, psym=1
        
     plot, lambda, skyspec
     oplot, lambda[wlocmax], skyspec[wlocmax], psym=1, color=fsc_color('red')
    
     return, [[linelambda], [sigma], [sigmaerr]]
  endif else return, [[-1], [-1], [-1]]
end

pro specres, science, qf=qf, goodoverride=goodoverride
    lambda = science.lambda
    wsky = where(science.skylinemask ne -1, cw)
    if cw eq 0 then begin
        fit = fitskylines(science)
        n = (size(fit))[1]
        if n gt 200 then message, 'Too many sky lines!'
        if n gt 0 then science.skyfit[0:n-1,*] = fit
        if n le 3 then begin
            science.dlam = replicate(1.37/2.35, n_elements(lambda))
            science.skylinemask = lonarr(n_elements(science.skylinemask))-1
            science.goodsky = 0
            message, 'Unstable sky line fit.  Using FWHM = 1.37 A', /info
            return
        endif

        w = where(2.35*fit[*,1] gt 0.8 and 2.35*fit[*,1] lt 7.0, cwprev)
        if cwprev lt 3 then begin
            science.dlam = replicate(1.37/2.35, n_elements(lambda))
            science.skylinemask = lonarr(n_elements(science.skylinemask))-1
            science.goodsky = 0
            message, 'Unusuable arc lines.  Using FWHM = 1.37 A', /info
            return
        endif

        ;quadratic fit
        qf = poly_fit(fit[w,0]/1000.0 - 7.8, 2.35*fit[w,1], 2, measure_errors=2.35*fit[w,2], chisq=chisq, /double, yfit=yfit)

        for j=0,4 do begin
            wnow = where(abs(2.35*fit[w,1] - yfit) lt 2.*2.35*fit[w,2], cw)
            if cw eq cwprev then break
            cwprev = cw
            if cw lt 3 then begin
                science.goodsky = 0
                message, 'The spectral resolution fit is very poor.', /info
                break
            endif
            w = w[wnow]
            qf = poly_fit(fit[w,0]/1000.0 - 7.8, 2.35*fit[w,1], 2, measure_errors=2.35*fit[w,2], chisq=chisq, /double, yfit=yfit)
        endfor
        n = (size(fit))[1]
        science.skylinemask = 0
        science.skylinemask[w] = 1
        if n lt 200 then science.skylinemask[n:n_elements(science.skylinemask)-1] = -1
    endif else begin
        wsky = where(science.skylinemask eq 1, csky)
        if csky lt 3 then begin
            science.dlam = replicate(1.37/2.35, n_elements(lambda))
            science.goodsky = 0
            message, 'Too few arc lines for new fit.', /info
            return
        endif
        fit = science.skyfit[wsky,*]
        qf = poly_fit(fit[wsky,0]/1000.0 - 7.8, 2.35*fit[wsky,1], 2, measure_errors=2.35*fit[wsky,2], chisq=chisq, /double, yfit=yfit)        
    endelse

    l = lambda / 1000. - 7.8
    dlam = poly(l, qf)
    dlam /= 2.35
    science.dlam = dlam
    if ~keyword_set(goodoverride) then science.goodsky = 1
end
