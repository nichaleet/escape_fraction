pro lyaequiv,redo=redo
;CALCULATE EQUIVALENT WIDTH OF FILES IN NAMES
;MAKE PLOT OF EW and Max-lowionization absorption depth
;REDO IF WANT To calculate again if not, it'll just plot
  namegalcap= ['MACS0940','A2219','A1689','MS1358']

names =    ['macs0940','a2219','a1689','ms1358arcac']
namereal = ['M0940, z=4.03','A2219, z=4.45','A1689, z=4.88','MS1358, z=4.93','RCS0224']
nameclump =    'ms1358_'+['clump0','clump1','clump4']
nameclumpreal = ['C0','C1+C2+C3','C4']

clight = 300000.
Lya    = 1215.67
set_plot,'x'
if n_Elements(redo) gt 0 then begin
   minvel = -800.
   maxvel = +1300.
   for j=0,1 do begin
      if j eq 0 then namenow = names
      if j eq 1 then namenow = nameclump
      for i=0,n_Elements(namenow)-1 do begin
;read in spec
         file = '/scr2/nichal/workspace3/SCIENCE/'+namenow(i)+'_science.fits'
         sci  = mrdfits(file,1,hdr,/silent)
         lambda = sci.lambda
         zlya   = sci.zsys_lya
         lambda = sci.lambda/(1.+zlya)
         minlam = (1.+minvel/clight)*lya 
         maxlam = (1.+maxvel/clight)*lya
         peakflux = max(sci.contdiv(where(lambda gt minlam and lambda lt maxlam)),cpeakflux)
         lambdac = (lambda(where(lambda gt minlam and lambda lt maxlam)))(cpeakflux)
         ok = 'n'
         nloop = 0
         while ok ne 'y' do begin
            if nloop gt 0 then read,minlam,prompt='min lambda:'
            if nloop gt 0 then read,maxlam,prompt='max lambda:'
            goodreg = where(lambda gt minlam and lambda lt maxlam,cgoodreg)
            if cgoodreg eq 0 then stop,'PROBLEM WITH REDSHIFT'
            lambdanow = lambda(goodreg)
            spec      = sci.contdiv(goodreg)
            specerr   = 1./sqrt(sci.contdivivar(goodreg))
            EW        = int_tabulated(lambdanow,1.-spec)
            EWerr     = sqrt(total(specerr^2))*median(abs(ts_diff(lambda,1)))
            plot,lambda,sci.contdiv,xtitle='lambda',ytitle='spec',xrange=[minlam-30,maxlam+30],yrange=[-2,peakflux+2],title=namenow(i)
            polyfill,[lambdac+EW/2.,lambdac+EW/2.,lambdac-EW/2.,lambdac-EW/2.],[0,1,1,0],color=fsc_color('cyan')
            oplot,lambdanow,spec+specerr,color=fsc_color('orange')
            oplot,lambdanow,spec-specerr,color=fsc_color('orange')
            oplot,lambda,sci.contdiv
            print,namenow(i),'EW=',EW,'+-',EWerr
            read,ok,prompt='ok? y or n:'
            nloop = nloop+1
         endwhile
         if tag_exist(sci,'LyaEW') then begin
            sci.LyaEW    = EW
            sci.LyaEWerr = EWerr
         endif else sci=create_struct(sci,'LyaEW',EW,'LyaEWerr',EWerr)
         mwrfits,sci,file,/create,/silent
      endfor
   endfor
endif

;OLD DATA
J13EW = [7.3,19.4,10.4]
J13fc = [0.68,0.33,0.89]
J13fcerr = [0.1,0.08,0.18]

;READ DATA
NL_EW    = fltarr(4)
NL_EWerr = fltarr(4)
NL_fc    = fltarr(4)
NL_fcerr = fltarr(4)
for i=0,3 do begin
   file = '/scr2/nichal/workspace3/SCIENCE/'+names(i)+'_science.fits'
   sci  = mrdfits(file,1,hdr,/silent)
   NL_EW(i)    = abs(sci.lyaew)
   NL_EWerr(i) = sci.lyaewerr
   NL_fc(i)    = 1.-min(sci.ave_covfrac[*,1],minpos)
   NL_fcerr(i) = 1./sqrt(sci.ave_covfrac[minpos,2]) 
endfor
;READ DATA CLUMP
NL_EWcl    = fltarr(3)
NL_EWclerr = fltarr(3)
NL_fccl    = fltarr(3)
NL_fcclerr = fltarr(3)
for i=0,2 do begin
   file = '/scr2/nichal/workspace3/SCIENCE/'+nameclump(i)+'_science.fits'
   sci  = mrdfits(file,1,hdr,/silent)
   NL_EWcl(i)    = abs(sci.lyaew)
   NL_EWclerr(i) = sci.lyaewerr
   NL_fccl(i)    = 1.-min(sci.ave_covfrac[*,1],minpos)
   NL_fcclerr(i) = 1./sqrt(sci.ave_covfrac[minpos,2])
endfor 

;PLOT 
colors=['orange','dark green','purple','dark red']
set_plot,'ps'
psname='output/lyaew.eps'
!p.font = 0
!x.margin=[10,6]
device, filename = psname,xsize = 11,ysize = 8, $
             xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
xrange = minmax([[J13EW,NL_EW,NL_EWcl]-5,[J13EW,NL_EW,NL_EWcl]+5])
plot, NL_EW,NL_fc,psym=1,xtitle='W!ILya!N('+cgsymbol('angstrom')+')',ytitle='Maximum absorption depth',xrange=[0,25],yrange=[0,1],/nodata,ystyle=8,position=[0.15,0.15,0.85,0.95]
for i=0,3 do begin
   oploterror,NL_EW[i],NL_fc[i],NL_Ewerr[i],NL_fcerr[i],color=fsc_color(colors[i]),psym=1
   cgplot, NL_EW[i],NL_fc[i],psym=14,color=fsc_color(colors[i]),/overplot,symsize=1.2
endfor
oploterror,NL_EWcl,NL_fccl,NL_Ewclerr,NL_fcclerr,color=fsc_color('darkred'),psym=1
cgplot, NL_EWcl,NL_fccl,psym=15,color=fsc_color('darkred'),/overplot
oploterror,j13ew,j13fc,j13fcerr,psym=1
cgplot, J13EW,J13fc,psym=14,color=fsc_color('navy'),/overplot,symsize=1.2
oplot,[0,25],[1,0.16],linestyle=2
axis,yaxis=1,yrange=[1,0]
xyouts,0.96,0.5,'f!Iesc',charsize=1.2,alignment=0.5,orientation=90,/normal
al_legend,namegalcap[0:2],psym=[14,14,14],color=fsc_color(colors[0:2]),position=[2,0.2],box=0,symsize=1,charsize=0.9
al_legend,[namegalcap[3],'Jones et al. 2013','MS1358 clumps'],psym=[14,14,15],color=fsc_color([colors[3],'blue',colors[3]]),position=[13,0.2],box=0,symsize=[1,1,0],charsize=0.9

device,/close
stop
end
