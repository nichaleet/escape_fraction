pro covfrac_si,objname,velarr,fc_mcmc,clmd_mcmc,fc_grid,clmd_grid,lineuse_in=lineuse_in,rangevel=rangevel
;This program should be preceded by covfrac.pro 
  fscience = '/scr2/nichal/workspace3/SCIENCE/'+objname+'_science.fits'
  science=mrdfits(fscience,1,/silent)
;science have tags: CONTDIV, LAMBDA, IVAR, DLAM, CONT, ZSPEC, ZSYS_LYA, MINSIABS,MAXSIABS, ZSYS_ISM
  clight   = 3.e5
  linelist = ['SiII1260','SiII1304','SiII1526']
  osc      = [1.22,0.0928,0.133] ;nist
  z  = science.zspec
  lambda   = science.lambda/(z+1.)
  linewl   = [1260.42,1304.37,1526.72]
  if n_elements(rangevel) eq 0 then rangevel = [-1000.,500.]
  stepv  = mean(science.dlam)/mean(science.lambda)*clight*2.355
  nstepv = round(abs(rangevel[1]-rangevel[0])/stepv)
  velarr = dblarr(nstepv)
  clmd_mcmc= dblarr(3,nstepv)    ;column density
  fc_mcmc  = dblarr(3,nstepv)    ;covering fraction
  rangen   = dblarr(2,nstepv)
  clmd_grid = dblarr(4,nstepv)
  fc_grid   = dblarr(4,nstepv)
  vminsi1304 = -200             ;km/s limit to overlap with OII1302
  checkrange = 'n'
  read,checkrange,prompt='Want to check range? (y/n):'
  for nv=0,nstepv-1 do begin
     vnowmin = rangevel[0]+nv*stepv  ;from -1000 km/s to 500 km/s 
     vnowmax = rangevel[0]+(nv+1.)*stepv 
     velarr(nv) = 0.5*(vnowmin+vnowmax)
     if n_elements(lineuse_in) eq 0 then lineuse = [1,1,1] else lineuse=lineuse_in ;default
     if velarr(nv) lt vminsi1304 then lineuse[1]=0
     absarr      = dblarr(3)
     absarr_ivar = dblarr(3)
     for il=0,2 do begin
        if lineuse(il) ne 0 then begin
           wlnow  = (lambda-linewl(il))/linewl(il)*clight
           goodwl = where(wlnow gt vnowmin and wlnow lt vnowmax,cwl)
           if cwl gt 0 then begin
              contdivnow = science.contdiv(goodwl)
              ivarnow    = science.contdivivar(goodwl)
              errnow     = 1./sqrt(ivarnow)
              absarr(il) = total(contdivnow*ivarnow)/total(ivarnow)
              ;absarr_ivar(il) = total(ivarnow)/total(ivarnow*(contdivnow-absarr(il))^2) ;weighted standard deviation                  
              absarr_ivar(il) = total(ivarnow) ;variance of weighted mean                    
           endif else print, 'no wl found for',linelist(il),'at velocity ',strtrim(string(velarr(nv)),2)
        endif
     endfor           
    ;Find the best column density and covering fraction at each velocity  
     goodline = where(lineuse ne 0) 
     flambda  = osc(goodline)*linewl(goodline)
     y        = absarr(goodline)
     y_err    = 1/sqrt(absarr_ivar(goodline))
     ;mcmc
     ;mcmc,x,y,yerr,nparam,nstep,p0,stepsizes,paraname,returnvalues
;     paraname = ['N','fc']
;     p0       = [3.768,0.5]  ;*10^14 per cm2 and covering fraction 
;     stepsizes = [0.1,0.1]
;     if total(~finite([y,y_err])) ne 0 then stop
;     mcmc_si,flambda,y,y_err,2,50000,p0,stepsizes,paraname,pfmcmc,velarr(nv),objname
;     ;pfmcmc is a 3x2 array of [first quartile, median, third quartile] of N and fc
;     fc_mcmc[*,nv]  = pfmcmc[*,1]
;     clmd_mcmc[*,nv]= pfmcmc[*,0]

;     science=create_struct(science,'fc_mcmc',fc_mcmc,'clmd_mcmc',clmd_mcmc,'velarr_mcmc',velarr)

     ;GRID SEARCH
     nsteps = [100,100] ; for N and fc
     range  = [[1.e-5,0.1],[0,1.0]] ;for N (units of *10^14 cm^-2) and fc
     ;if tag_exist(science,'range_N') then range[*,0] = science.range_n[*,nv]
     ;Calculate Probability
     gridsearch_si,flambda,y,y_err,nsteps,range,velarr(nv),probarr,narr,fcarr
     if checkrange eq 'y' then begin
        qfiner = 'n'
        read,qfiner,prompt='Want Finer Grid for N? (y/n):'
        while qfiner eq 'y' do begin
           read,minn,prompt='Input minimum N:'
           read,maxn,prompt='Input maximum N:'
           range[*,0] = [minn,maxn]
           gridsearch_si,flambda,y,y_err,nsteps,range,velarr(nv),probarr,narr,fcarr
           read,qfiner,prompt='Want Finer Grid? (y/n):'
        endwhile
     endif

     ;if good then marginalize,  keep data, and write plots
     ;marginalize
     p_n = dblarr(nsteps[0])
     p_fc= dblarr(nsteps[1])
     for i=0,nsteps[0]-1 do p_n(i) = int_tabulated(fcarr[i,*],probarr[i,*])
     for i=0,nsteps[1]-1 do p_fc(i)= int_tabulated(narr[*,i],probarr[*,i])
     clmd_grid[*,nv] = confidence_interval(narr[*,0],p_n)
     fc_grid[*,nv]   = confidence_interval(fcarr[0,*],p_fc)
     rangen[*,nv]    = range[*,0]
     set_plot,'ps'
     psname='/scr2/nichal/workspace3/output/'+objname+'_NFc_prob_vel'+strtrim(string(fix(velarr(nv))),2)+'.eps'
     device, filename = psname,xsize = 15,ysize = 15, $
             xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
     cgLoadCT, 33, CLIP=[30,255]
     !p.multi=[0,1,1]
     !p.charsize=1.3
     !p.font =0
     position =   [0.15, 0.1, 0.75, 0.75]
     cgImage, probarr, Stretch=1, MinValue=min(probarr), MaxValue=max(probarr), Position=position,XRange=range[*,0]*10, YRange=range[*,1]
   ;          /Axes, XTitle='N ($\tex10^{13}$ cm $\tex^{-2} (km s^{-1})^{-1}$)',YTitle='$\tex f_c$'
     axis,xaxis=0,XTitle='N(10!E13!Ncm!E-2!N(km s!E-1!N)!E-1!N)',color=fsc_color('black')
     axis,yaxis=0,YTitle='f!Dc!N',color=fsc_color('black')
     cgplot,narr[*,0]*10,p_n/10.,position=[0.15,0.75,0.75,0.95],xtickformat='(A1)',/noerase,ytickformat='(A1)',color=fsc_color('black')
     axis,yaxis=0,ytitle='P(N)',color=fsc_color('black'),yrange=range[*,1],ystyle=1,yticks=3,ytickv=[0.5,1.,1.5]
     vline,[clmd_grid[0,nv],clmd_grid[1,nv],clmd_grid[3,nv]]*10,color=fsc_color('indianred')
     cgplot,p_fc,fcarr[0,*],font=0,position=[0.75,0.1,0.95,0.75],xtickformat='(A1)',ytickformat='(A1)',/noerase,color=fsc_color('black')
     axis,xaxis=1,xtitle='P(f!Dc!N)',color=fsc_color('black'),xticks=4,xtickv=[1,2,3,4,5],xminor=1
     hline,[fc_grid[0,nv],fc_grid[1,nv],fc_grid[3,nv]],color=fsc_color('indianred'),linestyle=2
     device,/close
     set_plot,'x'
     !p.charsize=0
  endfor
  if tag_exist(science,'range_N') then science.range_N = rangen else science=create_struct(science,'range_N',rangen)
  if tag_exist(science,'clmd_grid') then science.clmd_grid = clmd_grid else science=create_struct(science,'clmd_grid',clmd_grid)
  if tag_exist(science,'fc_grid') then science.fc_grid = fc_grid else science=create_struct(science,'fc_grid',fc_grid)

  mwrfits,science,fscience,/create,/silent
  ;stop
end
