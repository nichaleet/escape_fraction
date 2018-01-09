pro combinefig
set_plot,'ps'
psname='output/combine_covfrac.eps'
device, filename = psname,xsize = 12,ysize = 10, $
             xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
names =    ['macs0940','a2219','a1689','ms1358arcac']
namereal = ['M0940, z=4.03','A2219, z=4.45','A1689, z=4.88','MS1358, z=4.93','RCS0224']

!p.charsize=3
!p.font=0
;multiplot,[1,1],mytitle='Average Absorption Profile'
xrange= [-1000.,600]
yrange= [0,1.2] 
position = [[0.15,0.77,0.85,0.97],[0.15,0.57,0.85,0.77],[0.15,0.37,0.85,0.57],[0.15,0.12,0.85,0.37]]
ytickv=[0.0,0.5,1.0]

for i=0,3 do begin
   file = '/scr2/nichal/workspace3/SCIENCE/'+names(i)+'_science.fits'
   sci  = mrdfits(file,1,hdr)  
   returnave = sci.ave_covfrac
   returnave_si = sci.ave_covfrac_si
   if i eq 0 then plot,returnave[*,0],returnave[*,1],psym=10,/nodata,yrange=yrange,color=fsc_color('black'),font=0,charsize=1,xrange=xrange,xstyle=1,position=position[*,i],ystyle=9,ytickinterval=0.5,xtickformat="(A1)"
   if i gt 0 and i lt 4 then plot,returnave[*,0],returnave[*,1],psym=10,/nodata,yrange=yrange,color=fsc_color('black'),font=0,charsize=1,xrange=xrange,xstyle=9,/noerase,position=position[*,i],ytickv=ytickv,ystyle=9,ytickinterval=0.5,xtickformat="(A1)"

   shaded_uncertainties,returnave[*,0],returnave[*,1],returnave[*,2],color='rose'
   velarr = returnave[*,0]
   fc_grid = sci.fc_grid
   
   lim = where(fc_grid[1,*] gt fc_grid[0,*],clim)
   if clim gt 0 then for j=0,clim-1 do fc_grid[1,lim(j)]=fc_grid[0,lim(j)]
   shaded_uncertainties_ul,velarr,1.-fc_grid[0,*],1.-fc_grid[1,*],1.-fc_grid[3,*],color='grey'
   oplot,returnave[*,0],returnave[*,1],psym=10,color=fsc_color('red')
   oplot,velarr,1.-fc_grid[0,*],psym=10,linestyle=1,color=fsc_color('black') 

   cgplot,[-2000,2000],[1,1],linestyle=2,color=fsc_color('darkgray'),/overplot   
   cgplot,[0,0],[0.,1.2],linestyle=2,color=fsc_color('darkgray'),/overplot   
   xyouts,580,0.1,namereal(i),/data,alignment=1,font=0,charsize=0
   axis,yaxis=1,yrange=[100,-20],ystyle=1,yticks=4,ytickv=[0,50,100],font=0,charsize=0,yminor=5
   if i eq 0 then begin
      al_Legend,[' Average line profiles',' Si II covering fraction'],psym=[15,15],color=[fsc_color('rose'),fsc_color('grey')],linestyle=[1,2],box=0,charsize=0.8,symsize=1,font=0,position=[-1000,0.45],/data
      oplot,[-980,-900],[0.3,0.3],color=fsc_color('red')
      oplot,[-970,-890],[0.1,0.1],color=fsc_color('black'),linestyle=1
   endif
endfor
axis,xrange=xrange,xtitle='velocity(km/s)',xstyle=1,font=0,charsize=0


xyouts,0.05,0.5,'Average low-ionization profile',/normal,orientation=90,alignment=0.5,font=0,charsize=0
xyouts,0.95,0.5,'Neutral covering fraction(%)',/normal,orientation=90,alignment=0.5,font=0,charsize=0

device,/close
stop
end
