pro runsfr_1358arcac

;initial input
name = 'ms1358arcac'
filestar = '/scr2/nichal/deimos/rawdata/2015mar/ms1358_standard/spec1d.Long1.0B.007.serendip1.fits'
filestarref  = 'standard_stars/feige98.dat'
namestar = 'Feige98'
file = '/scr2/nichal/workspace3/SCIENCE/'+name+'_science.fits'
fileall = '/scr2/nichal/workspace3/SCIENCE/'+name+'_scienceall.fits'
namearcs = ['arcC','arcA']
arcmags   = [22.25,23.8]
arcmagserr = [0.3,0.3]/3.
photmagwl  = [9033.-1031.,9033.+1031.] ;95% cummulative throughput width for F850lp

;read data
science  = mrdfits(file,1,hdr)
scienceall  = mrdfits(fileall,1,hdr)
starblu     = mrdfits(filestar,1,hdr)
starred     = mrdfits(filestar,2,hdr)
readcol,filestarref,reflam,refmag,what
stop
;get the conversion to convert the DEIMOS obs flux to magnitude
flux_conversion,starblu,starred,reflam,refmag,name,namestar 

;Do the conversion to science
param = mrdfits('/scr2/nichal/workspace3/flux_calib/'+name+'_'+namestar+'.fits',1)
fluxarcs = flux_convert(scienceall,param,name,namestar)
for i=0,n_elements(namearcs)-1 do begin
   if namearcs[i] eq 'arcC' then fluxarcs = fluxarcs([0,1,3])
   if namearcs[i] eq 'arcA' then fluxarcs = fluxarcs[2]
   arcmag = arcmags[i]+fltarr(n_elements(fluxarcs))
   arcmagerr = arcmagserr[i]+fltarr(n_elements(fluxarcs))

   fluxarcsreal = fixphot(fluxarcs,arcmag,arcmagerr,photmagwl,name,namestar)
   if n_elements(fluxarcsreal) gt 1 then fluxarctot = combineflux(fluxarcsreal,/plot) else fluxarctot=fluxarcsreal
save,fluxarcs,namearcs,fluxarctot,filename='/scr2/nichal/workspace3/flux_calib/ms1358_flux.sav'
;fluxarctot = fluxarcs
;cal UV slope(beta) and flux at 1600A and SFR
   if i eq 0 then results = caluvslope(fluxarctot,science.zspec,name,namearcs[i],namestar)
   if i eq 1 then results = caluvslope(fluxarctot,science.zspec,name,namearcs[i],namestar)
endfor

end
