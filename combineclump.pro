pro combineclump
set_plot,'ps'
psname='output/combine_covfrac_clumps.eps'
device, filename = psname,xsize = 12,ysize = 8, $
             xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
names =    ['clump0','clump1','clump4']
namereal = ['C0','C1+C2+C3','C4']

!p.charsize=3
!p.font=0
;multiplot,[1,1],mytitle='Average Absorption Profile'
xrange= [-1000.,300]
yrange= [0,1.2] 
position = [[0.15,0.69,0.85,0.96],[0.15,0.42,0.85,0.69],[0.15,0.15,0.85,0.42]]
ytickv=[0.0,0.5,1.0]

for i=0,2 do begin
   file = '/scr2/nichal/workspace3/SCIENCE/ms1358_'+names(i)+'_science.fits'
   sci  = mrdfits(file,1,hdr)  
   returnave = sci.ave_covfrac
   returnave_si = sci.ave_covfrac_si
   if i eq 0 then plot,returnave[*,0],returnave[*,1],psym=10,/nodata,yrange=yrange,color=fsc_color('black'),font=0,charsize=1,xrange=xrange,xstyle=1,position=position[*,i],ystyle=1,ytickinterval=0.5,xtickformat="(A1)"
   if i gt 0 and i lt 4 then plot,returnave[*,0],returnave[*,1],psym=10,/nodata,yrange=yrange,color=fsc_color('black'),font=0,charsize=1,xrange=xrange,xstyle=1,/noerase,position=position[*,i],ytickv=ytickv,ystyle=1,ytickinterval=0.5,xtickformat="(A1)"
   if i eq 2 then begin
      plot,returnave[*,0],returnave[*,1],psym=10,/nodata,yrange=yrange,color=fsc_color('black'),font=0,charsize=1,xrange=xrange,xstyle=1,/noerase,position=position[*,i],xtitle='velocity(km/s)',ytickv=ytickv,ystyle=1,ytickinterval=0.5
      al_Legend,[' Average line profiles',' Si II covering fraction'],psym=[15,15],color=[fsc_color('rose'),fsc_color('grey')],linestyle=[1,2],box=0,charsize=0.8,symsize=1,font=0,position=[-1000,0.45],/data
      oplot,[-980,-900],[0.3,0.3],color=fsc_color('red')
      oplot,[-970,-890],[0.1,0.1],color=fsc_color('black'),linestyle=1
   endif
   shaded_uncertainties,returnave[*,0],returnave[*,1],returnave[*,2],color='rose'
   velarr = returnave[*,0]
   fc_grid = sci.fc_grid
   
   lim = where(fc_grid[1,*] gt fc_grid[0,*],clim)
   if clim gt 0 then for j=0,clim-1 do fc_grid[1,lim(j)]=fc_grid[0,lim(j)]
   shaded_uncertainties_ul,velarr,1.-fc_grid[0,*],1.-fc_grid[1,*],1.-fc_grid[3,*],color='grey'
   oplot,returnave[*,0],returnave[*,1],psym=10,color=fsc_color('red')
   oplot,velarr,1.-fc_grid[0,*],psym=10,linestyle=1,color=fsc_color('black') 
   cgplot,[-2000,2000],[1,1],linestyle=2,color=fsc_color('black'),/overplot   
   cgplot,[0,0],[0.,1.2],linestyle=2,color=fsc_color('black'),/overplot   
   xyouts,280,0.1,namereal(i),/data,alignment=1,font=0,charsize=0
   axis,yaxis=1,yrange=[100,-20],ystyle=1,yticks=4,ytickv=[0,50,100],font=0,charsize=0

endfor
xyouts,0.05,0.5,'Average Absorption Profile',/normal,orientation=90,alignment=0.5,font=0,charsize=0
xyouts,0.95,0.5,'Neutral covering fraction(%)',/normal,orientation=90,alignment=0.5,font=0,charsize=0
device,/close
stop
end
