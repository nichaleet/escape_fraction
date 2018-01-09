pro gridsearch_si,x,y,yerr,nstep,range,vel,probarr,narr,fcarr
;INPUT
;range is [2*nparam] in size for min max
;nsteps has nparam dimension
;OUTPUT:
; probarr is normallized probability according to narr and fcarr grids
;make paramsarr
chisqarr= dblarr(nstep[0],nstep[1])
n_step  = (range[1,0]-range[0,0])/nstep[0]
fc_step = (range[1,1]-range[0,1])/nstep[1]
narr    = range[0,0]+n_step*findgen(nstep[0])
fcarr   = range[0,1]+fc_step*findgen(nstep[1])
narr    = rebin(narr,nstep[0],nstep[1])
fcarr   = rebin(transpose(fcarr),nstep[0],nstep[1])
print,'-----------------------------'
print,'N    fc    tau              Imodel              Iobs'
for in=0,nstep[0]-1 do begin
   for ic=0,nstep[1]-1 do begin
      nnow = narr[in,ic]
      fcnow= fcarr[in,ic]
      tau  = x*nnow/3.768
      Imodel = 1.-fcnow*(1.-exp(-1.*tau))
      chisqarr[in,ic] = total((Imodel-y)^2/yerr^2)
      if in mod 200 eq 0 and ic mod 20 eq 0 then begin
         print, nnow,fcnow,tau(0),imodel(0),y(0)
         ;print, nnow,fcnow,tau(1),imodel(1),y(1)
      endif
   endfor   
endfor

Likelihoodarr = -0.5*chisqarr
likelihoodarr = likelihoodarr-max(likelihoodarr)
probarr  = exp(likelihoodarr)
;normallize
print, 'start inttab'
volume   = int_tabulated_2d(narr,fcarr,probarr)
print, 'end inttab'
probarr  = probarr/volume

;cgDisplay, 600, 650, Title='Velocity'+strtrim(string(fix(vel)),2)+'km/s'

cgLoadCT, 33, CLIP=[30,255]
position =   [0.125, 0.125, 0.9, 0.800]
cgImage, probarr, Stretch=1, MinValue=min(probarr), MaxValue=max(probarr), $
       /Axes, XTitle='N', YTitle='fc',title='Velocity'+strtrim(string(fix(vel)),2)+'km/s', Position=position,XRange=range[*,0], YRange=range[*,1];, /Keep_Aspect
xyouts,0.1,0.1,'Number of lines used: '+strtrim(string(n_elements(x)),2),/normal
end
