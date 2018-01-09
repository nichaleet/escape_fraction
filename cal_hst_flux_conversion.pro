function cal_hst_flux_conversion,photfile
photall = read_ascii(photfile,comment_symbol='#')
photall = photall.field001
phot = {xcenter:reform(photall[0,*]),ycenter:reform(photall[1,*]),ra:reform(photall[2,*]),dec:reform(photall[3,*]),magauto:reform(photall[30,*]),magerr:reform(photall[31,*]),fluxauto:reform(photall[32,*]),fluxerr:reform(photall[33,*]),bck:reform(photall[17,*]),flag:reform(photall[18,*]),pix:reform(photall[19,*])}


fluxcgs =  10.^((phot.magauto+48.6)/(-2.5))/1.e-29;fnu
fluxcgserr = fluxcgs*0.4*alog(10.)*phot.magerr
gs = where(phot.flag le 1 and fluxcgserr gt 0, cgoodsource)
set_plot,'x'
!p.multi=[0,1,1]
ploterror,phot.fluxauto(gs),fluxcgs(gs),phot.fluxerr(gs),fluxcgserr(gs)
param=linfit(phot.fluxauto(gs),fluxcgs(gs),measure_errors=fluxcgserr(gs),sigma=sigma)
median_background = median(phot.bck/phot.pix)
conversion = [param[1]*1.e-29,sigma[1]*1.e-29,median_background]
print, 'conversion is', conversion ,'fnu/count'
wait, 1
return, conversion
end
