pro runcov1358_clump,redo=redo
;match_lambda,'spec1d.1358J.009.clump0.fits','spec1d.Long1.0B.007.clump0.fits',z=4.925
;match_lambda,'spec1d.1358J.009.clump1.fits','spec1d.Long1.0B.007.clump1.fits',z=4.925
;match_lambda,'spec1d.1358J.009.clump4.fits','spec1d.Long1.0B.007.clump4.fits',z=4.925
;match_lambda,'spec1d.1358J.009.clump5.fits','spec1d.Long1.0B.007.clump5.fits',z=4.925

dir='/scr2/nichal/workspace3/SPEC1D/'

name=['clump0','clump1','clump4','clump5']

for i=0,2 do begin

   if n_Elements(redo) eq 1 then prepspec,dir+['spec1d.1358J.009.'+name(i)+'.fits','spec1d.Long1.0B.007.'+name(i)+'.matched.fits'],'ms1358_'+name(i),zspec=4.925

;model1
   covfrac,'ms1358_'+name(i),returnave,returnave_si,lineuse=[1,1,1,1,0],zspec=4.925,rangevel=[-1000.,500.],redshiftdomain='ms1358arcac_science.fits'

;model2 (full model1)
   covfrac_si,'ms1358_'+name(i),velarr,fc_mcmc,clmd_mcmc,fc_grid,clmd_grid,lineuse_in=[1,1,1],rangevel=[-1000.,500.]

   set_plot,'ps' 
   psname='output/clump/ms1358_'+name(i)+'covfrac_all.eps'
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
endfor
end
