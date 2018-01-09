pro sfrplot

;;;;;;;;;;;;;;;;;;;;;;;SFR;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  clustername = ['MS1358','MS1358','MS1358','MS1358','MS1358','A1689','A1689','a2219','m0940','m0940']
  arcname     = ['arcC','arcA','c0','c123','c4','a2','a1','arc2','arcA','arcB']
  beta        = [-1.67,-1.60,-1.74,-1.89,-1.99,-1.92,-1.76,-1.76,-1.49,-1.35]
  betaerr     = [0.04,0.11,0.07,0.08,0.09,0.21,0.07,0.09,0.11,0.26]
  SFRim       = [804,215,333,182,79,98,248,143,275,363]
  SFRimerr    = [254,184,168,110,54,160,142,99,231,713]
  MUV         = [-24.3,-22.72,-23.51,-23.2,-22.53,-22.6,-23.23,-22.65,-22.74,-22.71]
  MUVerr      = [-.34,0.89,0.52,0.62,0.71,1.7,0.59,0.72,0.88,2.05]
  mu          = [12.1,2.9,10.0,11.5,13.1,2.2,5.3,1.6,182,38]
  muerr       = [2.6,0.05,0.8,1.1,1.5,0.002,0.01,0.01,411,2]
  muerrnew    = muerr


;adopt 20% and 50% uncertainty in magnification according to fig 6 in Zitrin2015
  toolow      = where(muerr/mu lt 0.2 and mu le 10,ctoolow)
  if ctoolow gt 0 then muerrnew(toolow) = 0.2*mu(toolow)
  toolow2     = where(muerr/mu lt 0.5 and mu gt 10 and mu lt 40, ctoolow2)
  if ctoolow2 gt 0 then muerrnew(toolow2) = 0.5*mu(toolow2)
  SFR     = SFRim/mu
  SFRerr  = SFR*sqrt((SFRimerr/SFRim)^2+(muerrnew/mu)^2)
  MUV = MUV+2.5*alog10(mu)
  MUVerr = sqrt(MUVerr^2+(2.5*muerrnew/mu/alog(10.))^2)
  print, 'Cluster      Arc       SFR          SFRerr       MUV      MUVerr'
  for i=0,n_Elements(arcname)-1 do print,clustername(i),'',arcname(i),sfr(i),sfrerr(i),MUV(i),MUVerr(i)


;If there are multiple images, take the weighted average of SFR
  namegal   = ['ms1358arcac','a1689','a2219','macs0940']
  namegalcap= ['MS1358','A1689','A2219','MACS0940']
  redshiftgal  = [4.927,4.875,4.450,4.031]
  SFRgal    = fltarr(4)
  SFRgalerr = fltarr(4)
  SFRgal[0] = wmean(SFR[0:1],SFRerr[0:1],error=errnow)
  SFRgalerr[0] = 0.5*sqrt(total((SFRerr[0:1])^2))
  SFRgal[1] = wmean(SFR[5:6],SFRerr[5:6],error=errnow)
  SFRgalerr[1] = 0.5*sqrt(total((SFRerr[5:6])^2))
  SFRgal[2] = SFR[7]
  SFRgalerr[2] = SFRerr[7]
  SFRgal[3] = wmean(SFR[8:9],SFRerr[8:9],error=errnow)
  SFRgalerr[3] = 0.5*sqrt(total((SFRerr[8:9])^2)) ;errnow

  MUVgal    = fltarr(4)
  MUVgalerr = fltarr(4)
  MUVgal[0] = wmean(MUV[0:1],MUVerr[0:1],error=errnow)
  MUVgalerr[0] = errnow
  MUVgal[1] = wmean(MUV[5:6],MUVerr[5:6],error=errnow)
  MUVgalerr[1] = errnow
  MUVgal[2] = MUV[7]
  MUVgalerr[2] = MUVerr[7]
  MUVgal[3] = wmean(MUV[8:9],MUVerr[8:9],error=errnow)
  MUVgalerr[3] = errnow


  betagal    = fltarr(4)
  betagalerr = fltarr(4)
  betagal[0] = wmean(beta[0:1],betaerr[0:1],error=errnow)
  betagalerr[0] = 0.5*sqrt(total((betaerr[0:1])^2))
  betagal[1] = wmean(beta[5:6],betaerr[5:6],error=errnow)
  betagalerr[1] = 0.5*sqrt(total((betaerr[5:6])^2))
  betagal[2] = beta[7]
  betagalerr[2] = betaerr[7]
  betagal[3] = wmean(beta[8:9],betaerr[8:9],error=errnow)
  betagalerr[3] = 0.5*sqrt(total((betaerr[8:9])^2)) ;errnow

  for i=0,n_Elements(namegal)-1 do print,namegal(i),sfrgal(i),sfrgalerr(i),muvgal(i),muvgalerr(i)

  nameclump = ['clump0','clump1','clump4']
  nameclumpcap = ['C0','C1+C2+C3','C4']
  SFRclump = SFR[2:4]
  SFRclumperr = SFRerr[2:4]
  redshiftclump = [4.927,4.927,4.927]
  betaclump = beta[2:4]
  betaclumperr = betaerr[2:4]
;lllllllllllllllllllllllll  Covering Fraction llllllllllllllllllllllllllllll
  fcovgal    = fltarr(4)
  fcovgalerr = fltarr(4)
  fcovclump  = fltarr(3)
  fcovclumperr = fltarr(3)
  for i=0,3 do begin
     file = '/scr2/nichal/workspace3/SCIENCE/'+namegal(i)+'_science.fits'
     sci  = mrdfits(file,1,hdr,/silent)
     fcovgal(i)    = 1.-min(sci.ave_covfrac[*,1],minpos)
     fcovgalerr(i) = 1./sqrt(sci.ave_covfrac[minpos,2]) 
  endfor
  for i=0,2 do begin
     file = '/scr2/nichal/workspace3/SCIENCE/ms1358_'+nameclump(i)+'_science.fits'
     sci  = mrdfits(file,1,hdr,/silent)
     fcovclump(i)    = 1.-min(sci.ave_covfrac[*,1],minpos)
     fcovclumperr(i) = 1./sqrt(sci.ave_covfrac[minpos,2])
  endfor 

;Jones13 galaxies
  Jclustername = ['Abell2390','Abell2390','Abell2390','Abell2390','J1261']
  Jarcname = ['H3a','H3b','H5a','H5b','thearc']
  Jbeta = [-1.75,-1.75,-2.8,-2.8,-1.9]
  Jbetaerr = [0.1,0.1,0.15,0.15,0.2]
  JSFRim = [1854,1282,71,87,1163]
  JSFRimerr = [396,274,23,28,501]
  Jmu    = [19.5,9.4,4.8,6.5,40]
  Jmuerr = 0.2*Jmu
  Jmuerr[4] = 0.5*Jmu[4]
  JSFRarcs    = JSFRim/Jmu
  JSFRerrarcs = JSFRarcs*sqrt((Jmuerr/Jmu)^2+(JSFRimerr/JSFRim)^2)
  JSFR = fltarr(3)
  JSFRerr = fltarr(3)
  jfcov    = [0.68,0.33,0.89]
  jfcoverr = [0.1,0.08,0.18]
  JSFR[0] = wmean(JSFRarcs[0:1],JSFRerrarcs[0:1],error=err)
  JSFRerr[0]=0.5*sqrt(total((JSFRerrarcs[0:1])^2))
  JSFR[1] = wmean(JSFRarcs[2:3],JSFRerrarcs[2:3],error=err)
  JSFRerr[1] = 0.5*sqrt(total((JSFRerrarcs[2:3])^2))
  JSFR[2] = JSFRarcs[4]
  JSFRerr[2] = JSFRerrarcs[4]
;;;;;;;;;;;PRINT OUT STATISTICS;;;;;;;;;;;;;;;;;;;;
fcovarr = [jfcov,fcovgal]
fcoverrarr = [jfcoverr,fcovgalerr]
betaarr = [betagal,jbeta]
betaerrarr = [betagalerr,jbetaerr]
meanerr,fcovarr,fcoverrarr,meanfcov,sigmam,sigmad,meanfcoverr
print,'average observed covering fraction = ',meanfcov,meanfcoverr
print, 'median',median(fcovarr),sigmad
;calculate fcovrel (relative cov frac)
;1) dust screen model
;find tau of lyc due to dust
A16 = (2.31*betaarr+4.85)  ;Calzetti2000
A16err = 2.31*betaerrarr
;e^-tau = 10^(-0.4*A16)
meanerr,a16,a16err,meana16,sigmam,sigmad,meana16err
print,'average extinction in magnitude = ',meana16,meana16err
fescabs = 1.-fcovarr/(fcovarr+(1.-fcovarr)*10.^(-0.4*A16))
sigma_2ndlowerterm = (1.-fcovarr)*10.^(-0.4*A16)*sqrt((fcoverrarr/fcovarr)^2+(0.4*alog(10.)*A16err)^2)
sigma_lowerterm = sqrt((fcoverrarr/fcovarr)^2+sigma_2ndlowerterm^2)
fescabserr = fescabs*sqrt((fcoverrarr/fcovarr)^2+(sigma_lowerterm/(fcovarr+(1.-fcovarr)*10.^(-0.4*A16)))^2)
meanerr,fescabs,fescabserr,meanfescabs,sigmam,sigmad,meanfescabserr
print, 'average absolute escape fraction = ', meanfescabs, meanfescabserr
print, 'median = ',median(fescabs),sigmad

jfescabs = fescabs[0:2]
fescabsgal = fescabs[3:6]
jfescabserr = fescabserr[0:2]
fescabsgalerr = fescabserr[3:6]

A16clump = (2.31*betaclump+4.85)  ;Calzetti2000
A16clumperr = 2.31*betaclumperr
fescclumpabs = 1.-fcovclump/(fcovclump+(1.-fcovclump)*10.^(-0.4*A16clump))
sigma_2ndlowerterm = (1.-fcovclump)*10.^(-0.4*A16clump)*sqrt((fcovclumperr/fcovclump)^2+(0.4*alog(10.)*A16clumperr)^2)
sigma_lowerterm = sqrt((fcovclumperr/fcovclump)^2+sigma_2ndlowerterm^2)
fescclumpabserr = fescclumpabs*sqrt((fcovclumperr/fcovclump)^2+(sigma_lowerterm/(fcovclump+(1.-fcovclump)*10.^(-0.4*A16clump)))^2)


;;CORRELATION
allcor = correlate(alog10([sfrgal,sfrclump,jsfr]),1.-[fcovgal,fcovclump,jfcov])
galcor = correlate(alog10([sfrgal,jsfr]),1.-[fcovgal,jfcov])
clumpcor =  correlate(alog10([sfrclump]),1.-[fcovclump])
;bootstrap
;sfrall = [sfrgal,sfrclump,jsfr]
;sfrallerr = [sfrgalerr,sfrclumperr,jsfrerr]
;fcovall = [fcovgal,fcovclump,jfcov]
;fcovallerr = [fcovgalerr,fcovclumperr,jfcoverr]

;sfrall = [sfrgal,jsfr]
;sfrallerr = [sfrgalerr,jsfrerr]
;fcovall = [fcovgal,jfcov]
;fcovallerr = [fcovgalerr,jfcoverr]

sfrall = [sfrclump]
sfrallerr = [sfrclumperr]
fcovall = [fcovclump]
fcovallerr = [fcovclumperr]

corarr = fltarr(1000)
for i=0,999 do begin
   sfrnow = sfrall+randomn(seed,n_Elements(sfrall))*sfrallerr
   badsfr = where(sfrnow lt 0.,cbadsfr)
   if cbadsfr gt 0 then sfrnow(badsfr)=1.
   fcovnow = fcovall+randomn(seed,n_elements(sfrall))*fcovallerr
   corarr[i] = correlate(alog10(sfrnow),1.-fcovnow)
endfor
print,'correlation coefficient for (all, galaxies, clumps): ',allcor,galcor,clumpcor
print,'bootstrap result (mean,median,stdev):',mean(corarr),median(corarr),stdev(corarr)
set_plot,'x'
plothist,corarr,xhist,yhist,bin=0.05
aa=gaussfit(xhist,yhist,a,nterms=4)
print, a[1:2]
stop
;;;;;;;;;;; PLOT SFR VS COV FRAC ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  sfrrange = [1.,130]
  set_plot,'ps'
  psname='/scr2/nichal/workspace3/output/sfr_vs_covfrac.eps'
  device, filename = psname,xsize = 11,ysize = 8, $
             xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  !p.font = 0
  plot,SFRgal,fcovgal,xtitle='SFR(Msun/yr)',ytitle='Maximum absorption depth',xrange=sfrrange,psym=1,yrange=[0,1],/xlog,/nodata,xstyle=1,ystyle=8,position=[0.15,0.15,0.87,0.95]
  colors=['dark red','purple','dark green','orange']
  for i=0,3 do begin
     oploterror,sfrgal[i],fcovgal[i],sfrgalerr[i],fcovgalerr[i],color=fsc_color(colors[i]),errcolor=fsc_color(colors[i])
     cgplot,sfrgal[i],fcovgal[i],color=fsc_color(colors[i]),psym=14,/overplot,symsize=1
  end  
  axis,yaxis=1,yrange=[1,0],ytitle='f!Desc!N'
  oploterror,sfrclump,fcovclump,sfrclumperr,fcovclumperr,color=fsc_color(colors[0]),errcolor=fsc_color(colors[0]),psym=1
  cgplot,sfrclump,fcovclump,color=fsc_color(colors[0]),psym=15,/overplot
  oploterror,jsfr,jfcov,jsfrerr,jfcoverr,errcolor=fsc_color('blue'),psym=1
  cgplot,jsfr,jfcov,color=fsc_color('blue'),psym=14,/overplot,symsize=1
  al_legend,namegalcap[0:2],psym=[14,14,14],color=fsc_color(colors[0:2]),position=[2,0.2],box=0,symsize=1,charsize=0.9
  al_legend,[namegalcap[3],'Jones et al. 2013','MS1358 clumps'],psym=[14,14,15],color=fsc_color([colors[3],'blue',colors[0]]),position=[10,0.2],box=0,symsize=[1,1,0],charsize=0.9
  device,/close

;;;;;;;;;;; PLOT A16 VS COV FRAC ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  set_plot,'ps'
  psname='/scr2/nichal/workspace3/output/EBV_vs_covfrac.eps'
  device, filename = psname,xsize = 11,ysize = 8, $
             xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  !p.font = 0
  EBVgal = (betagal*2.31+4.85)/4.39
  EBVgalerr = betagalerr*2.31/4.39
  EBVclump = (betaclump*2.31+4.85)/4.39
  EBVclumperr = betaclumperr*2.31/4.39
  jEBV = (jbeta*2.31+4.85)/4.39
  jEBVerr = jbetaerr*2.31/4.39
  EBVrange = [0,0.35]
  plot,EBVgal,fcovgal,xtitle='E(B-V)',ytitle='Maximum absorption depth',xrange=EBVrange,psym=1,yrange=[0,1],/nodata,xstyle=1,ystyle=8,position=[0.15,0.15,0.87,0.95]
  colors=['dark red','purple','dark green','orange']
  for i=0,3 do begin
     oploterror,EBVgal[i],fcovgal[i],EBVgalerr[i],fcovgalerr[i],color=fsc_color(colors[i]),errcolor=fsc_color(colors[i])
     cgplot,EBVgal[i],fcovgal[i],color=fsc_color(colors[i]),psym=14,/overplot,symsize=1
  end  
  axis,yaxis=1,yrange=[1,0],ytitle='f!Desc,obs!N'
  oploterror,EBVclump,fcovclump,EBVclumperr,fcovclumperr,color=fsc_color(colors[0]),errcolor=fsc_color(colors[0]),psym=1
  cgplot,EBVclump,fcovclump,color=fsc_color(colors[0]),psym=15,/overplot
  oploterror,jEBV,jfcov,jEBVerr,jfcoverr,errcolor=fsc_color('blue'),psym=1
  cgplot,jEBV,jfcov,color=fsc_color('blue'),psym=14,/overplot,symsize=1
  al_legend,namegalcap[0:2],psym=[14,14,14],color=fsc_color(colors[0:2]),box=0,symsize=1,charsize=0.9,/bottom,/left
  al_legend,[namegalcap[3],'Jones et al. 2013','MS1358 clumps'],psym=[14,14,15],color=fsc_color([colors[3],'blue',colors[0]]),box=0,symsize=[1,1,0],charsize=0.9,/bottom,/right
  device,/close

;;;;;;;;;;; PLOT A16 VS ABS COV FRAC ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  set_plot,'ps'
  psname='/scr2/nichal/workspace3/output/absEBV_vs_covfrac.eps'
  device, filename = psname,xsize = 11,ysize = 8, $
             xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  !p.font = 0
  fcovgalabs = 1-fescabsgal
  fcovgalabserr = fescabsgalerr
  jfcovabs = 1.-jfescabs
  jfcovabserr = jfescabserr
  fcovclumpabs = 1.-fescclumpabs
  fcovclumpabserr = fescclumpabserr

  EBVrange = [0,0.35]
  plot,EBVgal,fcovgalabs,xtitle='E(B-V)',ytitle='Absolute Covering Fraction',xrange=EBVrange,psym=1,yrange=[0,1],/nodata,xstyle=1,ystyle=8,position=[0.15,0.15,0.87,0.95]
  colors=['dark red','purple','dark green','orange']
  for i=0,3 do begin
     oploterror,EBVgal[i],fcovgalabs[i],EBVgalerr[i],fcovgalabserr[i],color=fsc_color(colors[i]),errcolor=fsc_color(colors[i])
     cgplot,EBVgal[i],fcovgalabs[i],color=fsc_color(colors[i]),psym=14,/overplot,symsize=1
  end  
  axis,yaxis=1,yrange=[1,0],ytitle='f!Desc,abs!N'
  oploterror,EBVclump,fcovclumpabs,EBVclumperr,fcovclumpabserr,color=fsc_color(colors[0]),errcolor=fsc_color(colors[0]),psym=1
  cgplot,EBVclump,fcovclumpabs,color=fsc_color(colors[0]),psym=15,/overplot
  oploterror,jEBV,jfcovabs,jEBVerr,jfcovabserr,errcolor=fsc_color('blue'),psym=1
  cgplot,jEBV,jfcovabs,color=fsc_color('blue'),psym=14,/overplot,symsize=1
  al_legend,namegalcap[0:2],psym=[14,14,14],color=fsc_color(colors[0:2]),box=0,symsize=1,charsize=0.9,/bottom,/left
  al_legend,[namegalcap[3],'Jones et al. 2013','MS1358 clumps'],psym=[14,14,15],color=fsc_color([colors[3],'blue',colors[0]]),box=0,symsize=[1,1,0],charsize=0.9,/bottom,/right
  device,/close

;;;;;;;;;WRITE OUT TABLE;;;;;;;;;;;;;;;;;;;;;;
print, '\tablehead{\colhead{Cluster} & \colhead{Arc} & \colhead{$\mu$} & \colhead{$\beta$} & \colhead{SFR}   & \colhead{SFR_{\textrm{err}}} &  \colhead{t$_{\text{exp}}$(hr)}}'   
print,'\startdata'

print, format='(A0,"&     &    &    &$",F10.2,"\pm",F6.2,"$\\")', clustername(0),SFRgal(0),SFRgalerr(0)
for i=0,1 do print, format='(A0,"&",A0,"&$",F10.1,"\pm",F5.1,"$&$",F10.2,"\pm",F6.2,"$&$",I,"\pm",I,"$\\")', clustername(i),arcname(i),mu(i),muerrnew(i),beta(i),betaerr(i),sfr(i),sfrerr(i)

print, format='(A0,"&     &    &    &$",F10.2,"\pm",F6.2,"$\\")', clustername(5),SFRgal(1),SFRgalerr(1)
for i=5,6 do print, format='(A0,"&",A0,"&$",F10.1,"\pm",F5.1,"$&$",F10.2,"\pm",F6.2,"$&$",I,"\pm",I,"$\\")', clustername(i),arcname(i),mu(i),muerrnew(i),beta(i),betaerr(i),sfr(i),sfrerr(i)

for i=7,7 do print, format='(A0,"&",A0,"&$",F10.1,"\pm",F5.1,"$&$",F10.2,"\pm",F6.2,"$&$",I,"\pm",I,"$\\")', clustername(i),arcname(i),mu(i),muerrnew(i),beta(i),betaerr(i),sfr(i),sfrerr(i)

print, format='(A0,"&     &    &    &$",F10.2,"\pm",F6.2,"$\\")', clustername(8),SFRgal(3),SFRgalerr(3)
for i=8,9 do print, format='(A0,"&",A0,"&$",F10.1,"\pm",F5.1,"$&$",F10.2,"\pm",F6.2,"$&$",I,"\pm",I,"$\\")', clustername(i),arcname(i),mu(i),muerrnew(i),beta(i),betaerr(i),sfr(i),sfrerr(i)

for i=2,4 do print, format='(A0,"&",A0,"&$",F10.1,"\pm",F5.1,"$&$",F10.2,"\pm",F6.2,"$&$",I,"\pm",I,"$\\")', clustername(i),nameclumpcap(i-2),mu(i),muerrnew(i),beta(i),betaerr(i),sfr(i),sfrerr(i)

;;;;;;;;;;Escape Fraction Redshift Evolution;;;;;;;;;;;;;;;;;;;;;;
;OLD DATA
  J13EW = [7.3,19.4,10.4]
  J13fc = [0.68,0.33,0.89,1.,0.91,0.69,1,0.4]
  jfcoverr = [0.1,0.08,0.18]
  J13z  = [4.04,4.015,4.13,2.73,3.07,2.38,2.73,2.96]
  J13age = galage(J13z,1000)/1.e9
  agegal = galage(redshiftgal,1000)/1.e9
  ageclump = galage(redshiftclump,1000)/1.e9

  psname='/scr2/nichal/workspace3/output/redshift_fc.eps'
  device, filename = psname,xsize = 11,ysize = 8, $
             xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  !p.font = 0
  plot,agegal,fcovgal,xtitle='Time since big bang (Gyr)',ytitle='Maximum absorption depth',xrange=[3.3,0.5],psym=1,yrange=[0,1.1],ystyle=9,xstyle=9,/nodata,position=[0.15,0.15,0.85,0.85]
  oploterror,j13age,j13fc[0:2],jfcoverr,color=fsc_color('blue'),errcolor=fsc_color('blue'),psym=1

  cgplot,J13age,J13fc,color=fsc_color('blue'),psym=14,/overplot,symsize=1.2
  oploterror,ageclump,fcovclump,fcovclumperr,color=fsc_color(colors[0]),errcolor=fsc_color('darkred'),psym=1
  cgplot,ageclump,fcovclump,color=fsc_color('darkred'),psym=15,/overplot
  for i=0,3 do begin
     oploterror,agegal[i],fcovgal[i],fcovgalerr[i],color=fsc_color('darkred'),errcolor=fsc_color(colors[i]),psym=1
     cgplot,agegal[i],fcovgal[i],color=fsc_color(colors[i]),psym=14,/overplot,symsize=1.2  
  endfor
  
  cgerrplot,[galage(3.,1000)/1.e9],[0.3],[0.65],color='black'
  cgplot,[galage(3.5,1000.)/1.e9],[0.4],psym=14,color='black',symsize=1.2,/overplot
  axis,yaxis=1,ytitle='f!Desc!N',yrange=[100,-10],yticks=5,ystyle=1,ytickv=[100,80,60,40,20,0],yminor=5
  z=[2,3,4,5,6,7]
  t=galage(z,1000.)/1.e9
  axis,xaxis=1,xrange=[3.3,0.5],xstyle=1,xticks=6,xtickv=t,xtickname=string(z,format='(I2)'),xtitle='Redshift'
  device,/close
stop

end

