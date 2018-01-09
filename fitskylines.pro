pro fitskylines,science
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
            wlocmax[i] = w[linestart[i]]+wmin
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

        ;lambdarange = max(lambda[wfit])-min(lambda[wfit])
        ;lambdastep = 0.03
        ;lambdafit = dindgen(round(lambdarange/lambdastep))*lambdastep+min(lambda[wfit])
        ;skyspecfit = interpol(skyspec[wfit], lambda[wfit], lambdafit, /quadratic)

        skyspecmax = max(skyspecfit, wmax)

        ;a = [lambdafit[wmax], 1.3/2.35, skyspecmax]
        ;yfit = curvefit(lambdafit, skyspecfit, lindgen(n_elements(lambdafit))+1.0, a, aerr, function_name='fitskyspec_gauss', chisq=chisq, status=status, /double, tol=1d-12)
        ;a = [skyspecmax, lambdafit[wmax], 1.3/2.35]

        guess = [skyspecmax, lambdafit[wmax], 0, medskyspec]
        yfit = gaussfit(lambdafit, skyspecfit, a, estimates=guess, sigma=aerr, nterms=4, chisq=chisq)
        sigma[i] = abs(a[2])
        sigmaerr[i] = aerr[2]
        linelambda[i] = a[1]

        ;print, sigma[i], sigmaerr[i]/sigma[i], chisq
        ;plot, lambdafit, skyspecfit
        ;oplot, lambdafit, yfit, color=fsc_color('red')
        ;wait, 0.5

        if chisq gt 1d-3 or sigmaerr[i] gt 0.8 then wgoodline[i] = 0
    endfor
    wgood = where(sigma gt 0 and wgoodline eq 1, c)
    if c gt 0 then begin
        linelambda = linelambda[wgood]
        sigma = sigma[wgood]
        sigmaerr = sigmaerr[wgood]
        ;ploterr, linelambda, sigma, sigmaerr, psym=1
        
        ;plot, lambda, skyspec
        ;oplot, lambda[wlocmax], skyspec[wlocmax], psym=1, color=fsc_color('red')

        return, [[linelambda], [sigma], [sigmaerr]]
    endif else return, [[-1], [-1], [-1]]
end



