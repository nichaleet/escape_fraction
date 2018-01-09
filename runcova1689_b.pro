pro runcova1689_b,redo=redo

;Combine the long slit data with the mask data. The main is long slit. The output is fout which is the same as normal spec1d from deep2 pipeline but with an extra tag in the structure called skydiv because we combine the sky_div data from each mask. the skyspec in the output is from the long slit only.

fa = '/scr2/nichal/deimos/rawdata/2015may/a1689/spec1d.a1689.009.arcz5120.fits'

covfrac,[fa],'a1689b',returnave,returnave_si,zspec=5.120,smoothfactor=1.,redo=redo
; returnave=[[velarr],[absarr],[absarr_wivar]]
; returnave_si = [[velarr],[si1260],[si1260_wivar]]

covfrac_si,'a1689b',velarr,fc_mcmc,clmd_mcmc,fc_grid,clmd_grid
; fc_grid = size[4,nstepvel] in peak,first quartile,median, third quartiles

set_plot,'ps'
psname='a1689b_covfrac_all.eps'
device, filename = psname,xsize = 14,ysize = 10, $
             xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
plot,returnave[*,0],returnave[*,1],psym=10,title=objname,ytitle='Average Absorbtion profile',xtitle='velocity(km/s)',/nodata,yrange=[0,1.5],color=fsc_color('black'),font=0,charsize=1
shaded_uncertainties,returnave[*,0],returnave[*,1],returnave[*,2],color='rose'
shaded_uncertainties_ul,velarr,1.-fc_grid[0,*],1.-fc_grid[1,*],1.-fc_grid[3,*],color='grey'
oplot,returnave[*,0],returnave[*,1],psym=10,color=fsc_color('red')
oplot,velarr,1.-fc_grid[0,*],psym=10,linestyle=1,color=fsc_color('black')
al_Legend,['Average line profiles','Si II covering fraction'],psym=[15,15],color=[fsc_color('red'),fsc_color('black')],linestyle=[1,2],box=0,thick=2,charsize=1,symsize=1.5,/right,/bottom,font=0
device,/close
stop
end
