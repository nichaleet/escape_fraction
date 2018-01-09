pro covfrac3,objname,velarr,fc_grid,clmd_grid,lineuse_in=lineuse_in,rangevel=rangevel

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
  rangen   = dblarr(2,nstepv)
  clmd_grid = dblarr(4,nstepv)
  fc_grid   = dblarr(4,nstepv)
  vminsi1304 = -200             ;km/s limit to overlap with OII1302
  checkrange = 'y'
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

     ;GRID SEARCH
     nsteps = [100,100] ; for N and fc
     range  = [[1.e-7,1.e-3],[0,1.0]] ;for N (units of *10^14 cm^-2) and fc
     if range[1,0] gt 0.2 then range[1,0]=0.2 ;must be optically thin
     if tag_exist(science,'range_N3') then range[*,0] = science.range_n3[*,nv]
     ;Calculate Probability
     gridsearch_si3,flambda,y,y_err,nsteps,range,velarr(nv),probarr,narr,fcarr
     if checkrange eq 'y' then begin
        qfiner = 'n'
        read,qfiner,prompt='Want Finer Grid for N? (y/n):'
        while qfiner eq 'y' do begin
           read,minn,prompt='Input minimum N:'
           read,maxn,prompt='Input maximum N:'
           range[*,0] = [minn,maxn]
           gridsearch_si3,flambda,y,y_err,nsteps,range,velarr(nv),probarr,narr,fcarr
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
     psname='output/model3/'+objname+'_NFc_prob_vel'+strtrim(string(fix(velarr(nv))),2)+'.eps'
     device, filename = psname,xsize = 25,ysize = 6, $
             xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
     cgLoadCT, 33, CLIP=[30,255]
     !p.multi=[0,3,1]
     position =   [0.2, 0.2, 0.9, 0.85]
     cgImage, probarr, Stretch=1, MinValue=min(probarr), MaxValue=max(probarr), $
              /Axes, XTitle='N ($\tex10^{13}$ cm $\tex^{-2} (km s^{-1})^{-1}$)',YTitle='$\tex f_c$', Position=position,XRange=range[*,0]*10, YRange=range[*,1],font=0;,title=objname+' Velocity'+strtrim(string(fix(velarr(nv))),2)+'km/s'
     cgplot,narr[*,0]*10,p_n/10.,xtitle='N ($\tex 10^{13}cm^{-2} (km s^{-1})^{-1}$)',ytitle='P(N)',font=0;,title=objname
     vline,[clmd_grid[0,nv],clmd_grid[1,nv],clmd_grid[3,nv]]*10,color=fsc_color('pink')
     cgplot,fcarr[0,*],p_fc,xtitle='$\tex f_c$',ytitle='$\tex P(f_c)$',font=0,title=''
     vline,[fc_grid[0,nv],fc_grid[1,nv],fc_grid[3,nv]],color=fsc_color('pink')
     device,/close
     set_plot,'x'
     !p.multi=[0,1,1]
  endfor
  if tag_exist(science,'range_N3') then science.range_N = rangen else science=create_struct(science,'range_N3',rangen)
  if tag_exist(science,'clmd_grid3') then science.clmd_grid = clmd_grid else science=create_struct(science,'clmd_grid3',clmd_grid)
  if tag_exist(science,'fc_grid3') then science.fc_grid = fc_grid else science=create_struct(science,'fc_grid3',fc_grid)

  mwrfits,science,fscience,/create,/silent
  ;stop
end
