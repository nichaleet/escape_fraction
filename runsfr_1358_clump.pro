pro runsfr_1358_clump

;initial input
filestar = '/scr2/nichal/deimos/rawdata/2015mar/ms1358_standard/spec1d.Long1.0B.007.serendip1.fits'
filestarref  = 'standard_stars/feige98.dat'
namestar = 'Feige98'
namearcs = 'ms1358_'+['clump0','clump1','clump4']
arcmags  = [23.09,23.285,23.904]
arcmagserr = [0.2,0.2,0.2]
;read data
starblu     = mrdfits(filestar,1,hdr)
starred     = mrdfits(filestar,2,hdr)
readcol,filestarref,reflam,refmag,what

for i=0,2 do begin
   name = namearcs[i]
   file = '/scr2/nichal/workspace3/SCIENCE/'+name+'_science.fits'
   fileall = '/scr2/nichal/workspace3/SCIENCE/'+name+'_scienceall.fits'
   science  = mrdfits(file,1,hdr)
   scienceall  = mrdfits(fileall,1,hdr)
;get the conversion to convert the DEIMOS obs flux to magnitude
   flux_conversion,starblu,starred,reflam,refmag,name,namestar 

;Do the conversion to science
   param = mrdfits('/scr2/nichal/workspace3/flux_calib/'+name+'_'+namestar+'.fits',1)
   fluxarcs = flux_convert(scienceall,param,name,namestar)
   
;calculate photometric magnitude by hand
;from arc_mag.pro
   arcmag = arcmags[i]+fltarr(n_Elements(fluxarcs))
   arcmagerr = arcmagserr[i]+fltarr(n_Elements(fluxarcs))
   photmagwl  = [9033.-1031.,9033.+1031.] ;95% cummulative throughput width for F850lp
;fix the magnitude
   fluxarcsreal = fixphot(fluxarcs,arcmag,arcmagerr,photmagwl,name,namestar)
   fluxarctot = combineflux(fluxarcsreal,/plot)
;fluxarctot = fluxarcs
;cal UV slope(beta) and flux at 1600A and SFR
   if i eq 0 then results = caluvslope(fluxarctot,science.zspec,name,namearcs[i],namestar,/noconfit) else results = caluvslope(fluxarctot,science.zspec,name,namearcs[i],namestar)
endfor

end
