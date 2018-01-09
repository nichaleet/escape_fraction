pro runcov0224,redo=redo

;match_lambda,'spec1d.rcs0224.007.arc1.013116_try2_adj.fits','spec1d.rcs0224.007.arc1.013016_try2_adj.fits',z=4.88
;match_lambda,'spec1d.rcs0224.007.arc1.013116_try2_adj.fits','spec1d.rcs0224.007.arc2.013016_try2_adj.fits',z=4.88
;match_lambda,'spec1d.rcs0224.007.arc1.013116_try2_adj.fits','spec1d.rcs0224.008.arc3.013016_try2_adj.fits',z=4.88
;match_lambda,'spec1d.rcs0224.007.arc1.013116_try2_adj.fits','spec1d.rcs0224.007.arc2.013116_try2_adj.fits',z=4.88
;match_lambda,'spec1d.rcs0224.007.arc1.013116_try2_adj.fits','spec1d.rcs0224.008.arc3.013116_try2_adj.fits',z=4.88

;;make telluric file
;tfile = 'spec1d.ms1358_Feige98.007.serendip1.fits'
;fitcontinuum_tell, tfile,name='rcs0224',/spline

dir = '/scr2/nichal/workspace3/SPEC1D/rcs0224/'

;if n_Elements(redo) eq 1 then prepspec,dir+['spec1d.rcs0224.007.arc1.013016_try2_adj.matched.fits','spec1d.rcs0224.007.arc1.013116_try2_adj.fits','spec1d.rcs0224.007.arc2.013016_try2_adj.matched.fits','spec1d.rcs0224.007.arc2.013116_try2_adj.fits','spec1d.rcs0224.008.arc3.013016_try2_adj.matched.fits','spec1d.rcs0224.008.arc3.013116_try2_adj.fits'],'rcs0224',zspec=4.88
if n_Elements(redo) eq 1 then prepspec,dir+['spec1d.rcs0224.008.arc3.013016_try2.matched.fits','spec1d.rcs0224.008.arc3.013116_try2.matched.fits'],'rcs0224',zspec=4.88

covfrac,'rcs0224',returnave,returnave_si,lineuse=[1,1,1,1,1],zspec=4.88
covfrac_si,'rcs0224',velarr,fc_mcmc,clmd_mcmc,fc_grid,clmd_grid,lineuse_in=[1,1,1];,rangevel=[-1000.,600.]
set_plot,'ps'
psname='output/rcs0224_covfrac_all.eps'
device, filename = psname,xsize = 14,ysize = 10, $
             xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
plot,returnave[*,0],returnave[*,1],psym=10,title=objname,ytitle='Average Absorbtion profile',xtitle='velocity(km/s)',/nodata,yrange=[0,1.5],color=fsc_color('black'),font=0,charsize=1,xstyle=1,xrange=minmax(returnave[*,0])
shaded_uncertainties,returnave[*,0],returnave[*,1],returnave[*,2],color='rose'
shaded_uncertainties_ul,velarr,1.-fc_grid[0,*],1.-fc_grid[1,*],1.-fc_grid[3,*],color='grey'
oplot,returnave[*,0],returnave[*,1],psym=10,color=fsc_color('red')
oplot,velarr,1.-fc_grid[0,*],psym=10,linestyle=1,color=fsc_color('black')
al_Legend,['Average line profiles','Si II covering fraction'],psym=[15,15],color=[fsc_color('red'),fsc_color('black')],linestyle=[1,2],box=0,thick=2,charsize=1,symsize=1.5,/right,/bottom,font=0
device,/close
stop

end
