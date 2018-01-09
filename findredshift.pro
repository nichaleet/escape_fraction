pro findredshift, science,objname, restlya, linelist ,linewl, lineuse, clight, shiftlya, shiftISM, delv_is_sys, zoommin, zoommax,zoommin2,zoommax2,fullmin,fullmax, fscience 

peak     = 0
centroid = 1
gausspeak = 0

;Systematic Redshift from Lyman alpha redshift 
z  = science.zspec
minlya = (1.-1000/3.e5)*restlya*(z+1) ;guess lya region
maxlya = (1.+1000/3.e5)*restlya*(z+1)

set_plot,'ps'
psname='output/'+objname+'_redshift.eps'
device, filename = psname,xsize = 30,ysize = 25, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
!p.font = 0
!p.charsize = 1
!p.multi=[0,2,3]
goodlya = where(science.lambda gt minlya and science.lambda lt maxlya,cgoodlya)
yrange = minmax(science.contdiv(goodlya))
plot,science.lambda,science.contdiv,yrange=yrange,xrange=[minlya,maxlya],xstyle=1
;plotlines,science.zspec
if cgoodlya gt 0 then begin
   case 1 of
      peak: begin
         ;peak location
         peakflux = max(science.contdiv[goodlya],peakloc)
         wllya  = (science.lambda[goodlya])(peakloc) ;peak location 
         zsys_lya   = (wllya/restlya+shiftlya/clight)-1. 
         science.method = 'peak'
      end
      centroid: begin
         wllya    = total(science.contdiv[goodlya]*science.lambda[goodlya])/total(science.contdiv[goodlya])
         zsys_lya = wllya/restlya-1.+shiftlya/clight
         science.method = 'centroid'
      end
      gausspeak: begin 
                                ;fit gaussian to the red side
         if objname eq 'ms1358' then minlya = wllya
         goodlya =  where(science.lambda gt minlya and science.lambda lt maxlya,cgoodlya)
         lambnow = science.lambda(goodlya)
         specnow = science.contdiv(goodlya)
         param   = [6,wllya,20.,0]
         gauss   = gaussfit(lambnow,specnow,param,sigma=sigma,measure_errors=sigmanow,nterms=4,estimates=param)
         wllya   = param(1)
         oplot,lambnow,gauss,color=fsc_color('pink')
         zsys_lya   = (wllya/restlya+shiftlya/clight)-1.
         science.method = 'gausspeak'
      end
   endcase
   plotuvlines,science.zspec,color='black'
   plotuvlines,zsys_lya,color='blue',delv=delv_is_sys
   al_Legend,['Lya Redshift','average ISM Redshift','individual ISM redshift'],psym=0,color=fsc_color(['Blue','dark green','yellow']),box=0,thick=5,linestyle=2,charsize=0,/right,/bottom
   science.zsys_lya = zsys_lya 
endif else stop,'Cannot find Lymanalpha redshift'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Systematic Redshift from absorption lines
;linelist =['SiII1260','OI1302','SiII1304','CII1334','SiII1526']
;linewl   = [1260.42,1302.168,1304.37,1334.532,1526.72]
minrange = (1.-1000/3.e5)*linewl*(z+1)
maxrange = (1.+1000/3.e5)*linewl*(z+1)

for i=0,n_Elements(linelist)-1 do begin
   if lineuse(i) eq 1 then begin
      good = where(science.lambda gt minrange(i) and science.lambda lt maxrange(i),cgood) 
      lambnow = science.lambda(good)
      weightnow = science.contdivivar(good)
      fluxnow = 2.-science.contdiv(good)
      case 1 of
         peak: begin
            peakflux = max(fluxnow,peakloc)
            wlis   = lambnow(peakloc) ;peak location
         end
         centroid: begin
            wlis    = total(fluxnow*weightnow*lambnow)/total(fluxnow*weightnow)
         end
         gausspeak:begin
            param   = [0.5,median(lambnow),20.,0]
            gauss   = gaussfit(lambnow,fluxnow,param,sigma=sigma,measure_errors=1./sqrt(weightnow),nterms=4,estimates=param)
            wlis   = param(1)
         end
      endcase
   science.zsys_ism_arr(i)   = (wlis/linewl(i)+shiftISM/clight)-1.
   print, 'Doing Absorption Lines Redshift with method: ',science.method
   print, 'At ',linewl(i),'A, z_sys_ism =',science.zsys_ism_arr(i)
   endif
endfor
science.zsys_ism_ave = mean(science.zsys_ism_arr(where(lineuse eq 1)))
science.zsys_ism_std = stdev(science.zsys_ism_arr(where(lineuse eq 1)))
delv = science.zsys_ism_std*3.e5
print, 'mean =',science.zsys_ism_ave, '  SD=',science.zsys_ism_std, ' delta V=',delv
print, 'delta Z total = ', sqrt(delv^2+delv_is_sys^2)/clight
print, 'z_sys_lya=',science.zsys_lya

for i=0,n_Elements(linelist)-1 do begin
   if lineuse(i) eq 1 then begin
      good = where(science.lambda gt minrange(i) and science.lambda lt maxrange(i),cgood) 
      lambnow = science.lambda(good)
      fluxnow = science.contdiv(good)
      skynow = 1./sqrt(science.contdivivar(good))
      plot,lambnow,fluxnow,title=linelist(i),/nodata
      plotuvlines,zsys_lya,color='blue',delv=delv_is_sys
      delwl = linewl(i)*sqrt(delv^2+delv_is_sys^2)/clight
      ymax = max(!y.crange)
      ymin = min(!y.crange)
      linewlshift = linewl[i]*(1.+science.zsys_ism_ave)
      cgcolorfill,[linewlshift-delwl,linewlshift+delwl,linewlshift+delwl,linewlshift-delwl,linewlshift-delwl],[ymax,ymax,ymin,ymin,ymax],color=fsc_color('cyan')
      vline,linewl(i)*(1.+science.zsys_ism_arr(i)),color=fsc_color('yellow')
      vline,linewl(i)*(1.+science.zsys_ISM_ave),color=fsc_color('dark green')
      oplot,lambnow,skynow/max(skynow)*max(!y.crange),color=fsc_color('pink')
      xyouts,!x.crange[0]+1,!y.crange[1]*0.9,'pink: contdiv_err/'+strtrim(string(max(skynow)*max(!y.crange)),2)
      oplot,lambnow,fluxnow
   endif
endfor
device,/close
end
