pro sfrplot

;;;;;;;;;;;;;;;;;;;;;;;SFR;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  clustername = ['MS1358','MS1358','MS1358','MS1358','MS1358','A1689','A1689','a2219','m0940','m0940']
  arcname     = ['arcBC','arcA','c0','c123','c4','a2','a1','arc2','arcA','arcB']
  beta        = [-1.59,-1.60,-1.74,-1.89,-1.99,-1.92,-1.76,-1.76,-1.49,-1.35]
  betaerr     = [0.03,0.11,0.07,0.08,0.09,0.21,0.07,0.09,0.11,0.26]
  SFRim       = [1189,215,333,182,79,98,248,143,275,363]
  SFRimerr    = [312,184,168,110,54,160,142,99,231,713]
  MUV         = [-24.56,-22.72,-23.51,-23.2,-22.53,-22.6,-23.23,-22.65,-22.74,-22.71]
  MUVerr      = [-.27,0.89,0.52,0.62,0.71,1.7,0.59,0.72,0.88,2.05]
  mu          = [14.2,2.9,10.0,11.5,13.1,2.2,5.3,1.6,182,38]
  muerr       = [1.6,0.05,0.8,1.1,1.5,0.002,0.01,0.01,411,2]
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
  namegal   = ['ms1358','a1689','a2219','macs0940']
  namegalcap= ['MS1358','A1689','A2219','MACS0940']
  redshiftgal  = [4.927,4.875,4.450,4.031]
  SFRgal    = fltarr(4)
  SFRgalerr = fltarr(4)
  SFRgal[0] = wmean(SFR[0:1],SFRerr[0:1],error=errnow)
  SFRgalerr[0] = errnow
  SFRgal[1] = wmean(SFR[5:6],SFRerr[5:6],error=errnow)
  SFRgalerr[1] = errnow
  SFRgal[2] = SFR[7]
  SFRgalerr[2] = SFRerr[7]
  SFRgal[3] = wmean(SFR[8:9],SFRerr[8:9],error=errnow)
  SFRgalerr[3] = errnow

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

  for i=0,n_Elements(namegal)-1 do print,namegal(i),sfrgal(i),sfrgalerr(i),muvgal(i),muvgalerr(i)

  nameclump = ['clump0','clump1','clump4']
  nameclumpcap = ['C0','C1+C2+C3','C4']
  SFRclump = SFR[2:4]
  SFRclumperr = SFRerr[2:4]
  redshiftclump = [4.927,4.927,4.927]
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

;;;;;;;;;;; PLOT SFR VS COV FRAC ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  sfrrange = [1.,max([SFRgal,SFRclump])]
  set_plot,'ps'
  psname='/scr2/nichal/workspace3/output/sfr_vs_covfrac.eps'
  device, filename = psname,xsize = 11,ysize = 8, $
             xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  !p.font = 0
  plot,SFRgal,fcovgal,xtitle='SFR(Msun/yr)',ytitle='Maximum absorption depth',xrange=sfrrange,psym=1,yrange=[0,1],/xlog,/nodata
  colors=['dark red','purple','dark green','orange']
  for i=1,3 do begin
     oploterror,sfrgal[i],fcovgal[i],sfrgalerr[i],fcovgalerr[i],color=fsc_color(colors[i]),errcolor=fsc_color(colors[i])
     cgplot,sfrgal[i],fcovgal[i],color=fsc_color(colors[i]),psym=14,/overplot
  end  
  oploterror,sfrclump,fcovclump,sfrclumperr,fcovclumperr,color=fsc_color(colors[0]),errcolor=fsc_color(colors[0]),psym=1
  cgplot,sfrclump,fcovclump,color=fsc_color(colors[0]),psym=15,/overplot
  al_legend,[namegalcap[1:3],'MS1358 clumps'],psym=[14,14,14,15],color=fsc_color([colors[1:3],colors[0]]),position=[10,0.35],box=0
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
  J13z  = [4.04,1.015,4.13,2.73,3.07,2.38,2.73,2.96]
  J13age = galage(J13z,1000)/1.e9
  agegal = galage(redshiftgal,1000)/1.e9
  ageclump = galage(redshiftclump,1000)/1.e9

  psname='/scr2/nichal/workspace3/output/redshift_fc.eps'
  device, filename = psname,xsize = 11,ysize = 8, $
             xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  !p.font = 0
  plot,agegal,fcovgal,xtitle='Time since big bang (Gyr)',ytitle='Maximum absorption depth',xrange=[3.,0],psym=1,yrange=[0,1.1],/nodata
  oploterror,agegal[1:3],fcovgal[1:3],fcovgalerr[1:3],color=fsc_color('darkred'),errcolor=fsc_color('dark red'),psym=1
  cgplot,agegal[1:3],fcovgal[1:3],color=fsc_color('dark red'),psym=14,/overplot  
  oploterror,ageclump,fcovclump,fcovclumperr,color=fsc_color(colors[0]),errcolor=fsc_color('darkred'),psym=1
  cgplot,ageclump,fcovclump,color=fsc_color('darkred'),psym=15,/overplot
  cgplot,J13age,J13fc,color=fsc_color('orange'),psym=14,/overplot
  al_legend,['This paper','clumps','J13'],psym=[14,15,14],color=fsc_color(['darkred','darkred','orange']),/right,/bottom,box=0
  device,/close
stop

end

