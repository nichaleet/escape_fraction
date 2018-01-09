pro runsfr_1689

;initial input
name = 'a1689'
filestar = '/scr2/nichal/deimos/rawdata/2015may/a1689_star/spec1d.Long1.0B.007.star1689.fits'
filestarref  = 'standard_stars/bd28d4211.dat'
namestar = 'BD28D4211'
file = '/scr2/nichal/workspace3/SCIENCE/'+name+'_science.fits'
fileall = '/scr2/nichal/workspace3/SCIENCE/'+name+'_scienceall.fits'
namearcs = ['arca1','arca2']
;read data
science  = mrdfits(file,1,hdr)
scienceall  = mrdfits(fileall,1,hdr)
starblu     = mrdfits(filestar,1,hdr)
starred     = mrdfits(filestar,2,hdr)
readcol,filestarref,reflam,refmag,what

;get the conversion to convert the DEIMOS obs flux to magnitude
flux_conversion,starblu,starred,reflam,refmag,name,namestar 

;Do the conversion to science
param = mrdfits('/scr2/nichal/workspace3/flux_calib/a1689_BD28D4211.fits',1)
fluxarcs = flux_convert(scienceall,param,name,namestar)

;get magnitude from photometry
photfile = '/scr2/nichal/deimos/May15_masks/abell1689/hst_9289_56_acs_wfc_f850lp_sexphot_trm.cat'
photall = read_ascii(photfile,comment_symbol='#')
photall = photall.field001
phot = {xcenter:reform(photall[0,*]),ycenter:reform(photall[1,*]),ra:reform(photall[2,*]),dec:reform(photall[3,*]),magauto:reform(photall[30,*]),magerr:reform(photall[31,*]),fluxauto:reform(photall[32,*]),fluxerr:reform(photall[33,*]),magbest:reform(photall[34,*]),magbesterr:reform(photall[35,*]),fluxbest:reform(photall[36,*]),fluxbesterr:reform(photall[37,*])}
bestguessx = [4292,2723]
bestguessy = [2230,2987]
photmag    = []
photmagerr = []
for i=0,n_elements(bestguessx)-1 do begin
   distance = sqrt((bestguessx[i]-phot.xcenter)^2+(bestguessy[i]-phot.ycenter)^2)
   min = min(distance,loc)      ;I check, it's correct 
   print, 'RA DEC', phot.ra[loc],phot.dec[loc]
   photmag = [photmag,phot.magauto[loc]]     ;ABmag
   photmagerr = [photmagerr,phot.magerr[loc]];ABmag
endfor

photmagwl  = [9033.-1031.,9033.+1031.] ;95% cummulative throughput width for F850lp

;fix the magnitude with photometry average
fluxarcsreal = fixphot(fluxarcs,photmag,photmagerr,photmagwl,name,namestar)
fluxarctot = combineflux(fluxarcsreal,/plot)
save,fluxarcs,namearcs,fluxarctot,filename='/scr2/nichal/workspace3/flux_calib/a1689flux.sav'
for i =0,1 do begin

   fluxarctot = fluxarcsreal[i]
;cal UV slope(beta) and flux at 1600A and SFR
   if i eq 0 then results = caluvslope(fluxarctot,science.zspec,name,namearcs[i],namestar)
   if i eq 1 then results = caluvslope(fluxarctot,science.zspec,name,namearcs[i],namestar,wlmax=1566.,wlmin=1290.)

endfor
stop
end
