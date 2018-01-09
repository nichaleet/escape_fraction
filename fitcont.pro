pro fitcont,filesprof
;filesprof is produced in nl_do_extract.pro e.g. profile.slit.rcs0224.007R.fits
;Fit with Polynomial to remove the background sky from 2dslit.
  sprof = mrdfits(filesprof,0)
  objpos = mrdfits(filesprof,1) ;desipos
  fwhm = mrdfits(filesprof,2)*2.1   ;fwtm  ;THIS CAN BE TWEAK
  nspatial = n_elements(sprof)
  x = findgen(nspatial)
  mask = bytarr(nspatial)+1
  good = 'n'
  while good ne 'y' do begin
     for i = 0,n_elements(objpos)-1 do begin
        xmin = round(objpos[i]-fwhm[i])>0
        xmax = round(objpos[i]+fwhm[i])<(nspatial-1)
        mask[xmin:xmax] = 0
     endfor
                                ;fix end
     mask[0:10] = 0
     mask[nspatial-10:nspatial-1]=0
     bck = where(mask,nbck)
     ndeg = 1                   ;THIS CAN BE CHANGED (degree of polynomial)
     if nbck gt 0 then begin
        params = poly_fit(x(bck),sprof(bck),ndeg)
        case ndeg of
           3: ymodel = params(0)+params(1)*x+params(2)*x^2+params(3)*x^3
           2: ymodel = params(0)+params(1)*x+params(2)*x^2
           1: ymodel = params(0)+params(1)*x
           else: stop,'CHECK YOUR DEGREE'
        endcase 
        
        plot,x,sprof
        oplot,x(bck),sprof(bck),psym=1,color=fsc_color('red')
        oplot,x,ymodel,color=fsc_color('green')
        oplot,!x.crange,[0,0],color=fsc_color('pink')
     endif
     
     read,good,prompt='good? (Y/N)'
     if good eq 'y' then save,params,filename='polynomial.'+filesprof+'.sav'
     stop
  endwhile
end
