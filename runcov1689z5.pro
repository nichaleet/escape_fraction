pro runcov1689z5,redo=redo
;match_lambda,'spec1d.a1689.007.obja1.fits','spec1d.a1689.012.obja2.fits',z=4.868
;!!!!!!!!!!!!!!!!
; IF RUN REDO, GO CHANGE THE MASK REGION IN FITCONTINUUM.PRO TO -200 TO 200
; ALSO CHANGE THE BKSPACE TO LOW 10 -15 PIX
;!!!!!!!!!!!!!!!!!
dir='/scr2/nichal/workspace3/SPEC1D/'

if n_Elements(redo) eq 1 then prepspec,dir+['spec1d.a1689.009.arcz5120.fits','spec1d.a1689.009.arcz5120.fits'],'a1689z5',zspec=5.120

covfrac,'a1689z5',returnave,returnave_si,lineuse=[1,1,1,1,0],zspec=5.120,rangevel=[-1000.,600.]

covfrac_si,'a1689z5',velarr,fc_mcmc,clmd_mcmc,fc_grid,clmd_grid,lineuse_in=[1,1,1],rangevel=[-1000.,600.]
; fc_grid = size[4,nstepvel] in peak,first quartile,median, third quartiles

set_plot,'ps'
psname='output/a1689z5_covfrac_all.eps'
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
end
