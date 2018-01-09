pro mcmc_si,x,y,yerr,nparam,nstep,p0,stepsizes,paraname,returnvalues,vel,objname
set_plot,'ps'
;before running mcmc pls make a new calymodel.pro
psname=objname+'_mcmc_vel'+strtrim(string(fix(vel)),2)+'.eps'
device, filename = psname,xsize = 22,ysize = 25, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
setplot,14
!x.omargin = [2,6]
!y.omargin = [2,6]
!p.multi = [0,2,3]
!p.font  = 0
;stop
ymodel0 = calymodel(x,p0)
oldchisq = calchisq(ymodel0,y,yerr)

print,'Chisq of the initial parameters:', oldchisq

randomparam = fix(randomu(seed,nstep)*nparam) 
randomstepsize = randomn(seed,nstep)
randomaccept = randomu(seed,nstep)

parr = fltarr(nparam,nstep+1)
chisqarr = fltarr(nstep+1)
pold = p0
parr(0) = p0
pcurrent = p0
chisqarr(0) = oldchisq
naccepts = 0UL
nrejects = 0UL

for ii=0UL,nstep-1 do begin
   if (ii mod 10000.) eq 0. then print, 'Doing step',ii
   chosenparam = randomparam(ii)
   pcurrent(chosenparam) = pcurrent(chosenparam)+stepsizes(chosenparam)*randomstepsize(ii)
   ymodel = calymodel(x,pcurrent)
   newchisq = calchisq(ymodel,y,yerr)
   if pcurrent[1] gt 1. or pcurrent[1] lt -0.5 then newchisq = 1.e8 ;covering fraction is less than 1
   if newchisq lt oldchisq then begin 
      naccepts = naccepts+1
   endif else begin
      prob = exp(-0.5*(newchisq-oldchisq))
      if randomaccept(ii) lt prob then begin ;accept
         naccepts = naccepts+1
      endif else begin    ;reject
         pcurrent = pold    ; do not accept new parameters
         newchisq = oldchisq 
         nrejects = nrejects+1
      endelse
   endelse
   ;stop
   parr(*,ii+1) = pcurrent
   chisqarr(ii+1) = newchisq
   pold = pcurrent
   oldchisq = newchisq
endfor
print, 'n accepts =', naccepts
print, 'n rejects =',nrejects
print, 'acceptance ratio', float(naccepts)/float(nrejects+naccepts)

plot,smooth(chisqarr,n_elements(chisqarr)/1000.),xtitle='step number', ytitle='chisq',psym=10,yrange=cgpercentiles(chisqarr,percentiles=[0,0.95])

medchi = median(chisqarr)
minstep = where(chisqarr le medchi)
minstep = minstep(0)

print, 'Trimming position is', minstep,' at chisq ', chisqarr(minstep)

chisqarr_trimmed = chisqarr
remainedind = indgen(nstep,/long)
if minstep gt 0 then begin
   removeind = indgen(minstep,/long)
   remove, removeind,chisqarr_trimmed,remainedind
endif
parr_trimmed = parr[*,[remainedind]]
plot,parr_trimmed[0,*],parr_trimmed[1,*],psym=1,xtitle=paraname(0),ytitle=paraname(1)
parr_sort = 0.*parr_trimmed

nelement = n_elements(remainedind)
sigmapos = long([0.16,0.5,0.84]*nelement)
returnvalues = fltarr(3,nparam)

print,'After trimming, there are ',nelement, ' elements left.'

;Calculating outputs and plotting
for jj=0,nparam-1 do begin
   sortind = bsort(parr_trimmed(jj,*),asort)
   parr_Sort(jj,*) = asort
   plot,parr(jj,*),xtitle='step number', ytitle=paraname(jj),yrange=[min(parr(jj,*))-.5,max(parr(jj,*))+.5]
   ;plothist,parr_sort(jj,*),xtitle = paraname(jj),/autobin
   bin = (max(parr_sort(jj,*))-min(parr_sort(jj,*)))/100.
   plothist,parr_sort(jj,*),xhist,yhist,bin=bin,/noplot,peak=1.
   yhist= yhist/tsum(xhist,yhist)
   plot,xhist,yhist,ytitle='posterior PDF',xtitle=paraname(jj)
   xloc = (!x.crange[1]-!x.crange[0])*0.1+!x.crange[0]
   xyouts,xloc,!y.crange[1]*0.9,strtrim(string(parr_sort[jj,sigmapos[0]]),2),charsize=1,color=fsc_color('blue')
   xyouts,xloc,!y.crange[1]*0.8,strtrim(string(parr_sort[jj,sigmapos[1]]),2),charsize=1,color=fsc_color('blue')
   xyouts,xloc,!y.crange[1]*0.7,strtrim(string(parr_sort[jj,sigmapos[2]]),2),charsize=1,color=fsc_color('blue')
   returnvalues(*,jj) =  parr_sort[jj,[sigmapos]]
endfor
xyouts, 0.5, 0.95,objname+' velocity='+strtrim(fix(vel),2)+'km/s DOF='+strtrim(string(n_elements(x)-nparam),2), ALIGNMENT=0.5, CHARSIZE=1.25, /NORMAL
;stop
device,/close
set_plot,'x'
end
