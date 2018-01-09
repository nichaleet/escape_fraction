pro runcov0940,redo=redo
;match_lambda,'spec1d.Long1.5.007.arcA.fits','spec1d.Long1.5.007.arcB.fits',z=4.04

dir='/scr2/nichal/workspace3/SPEC1D/'

if n_Elements(redo) eq 1 then prepspec,dir+['spec1d.Long1.5.007.arcA.fits','spec1d.Long1.5.007.arcB.matched.fits'],'macs0940',zspec=4.03

;model1
covfrac,'macs0940',returnave,returnave_si,lineuse=[1,1,1,1,1],zspec=4.03,rangevel=[-1000.,600.]

;model2 (full model1)
covfrac_si,'macs0940',velarr,fc_mcmc,clmd_mcmc,fc_grid,clmd_grid,lineuse_in=[1,1,1],rangevel=[-1000.,600.]
; fc_grid = size[4,nstepvel] in peak,first quartile,median, third quartiles

set_plot,'ps'
psname='output/macs0940_covfrac_all.eps'
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
