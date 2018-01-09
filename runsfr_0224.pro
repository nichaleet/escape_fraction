pro runsfr_0224

;initial input
name = 'rcs0224'
filestar = '/scr2/nichal/deimos/rawdata/2016jan/013016/rcs0224_standard/spec1d.Long1.0B.007.feige25.fits'
filestarref  = 'standard_stars/feige25.dat'
namestar = 'Feige25'
file = '/scr2/nichal/workspace3/SCIENCE/'+name+'_science.fits'
fileall = '/scr2/nichal/workspace3/SCIENCE/'+name+'_scienceall.fits'
namearcs = ['arc1','arc2','arc3']
;read data
science  = mrdfits(file,1,hdr)
scienceall  = mrdfits(fileall,1,hdr)
starblu     = mrdfits(filestar,1,hdr)
starred     = mrdfits(filestar,2,hdr)
readcol,filestarref,reflam,refmag,what

;get the conversion to convert the DEIMOS obs flux to magnitude
flux_conversion,starblu,starred,reflam,refmag,name,namestar 

;Do the conversion to science
param = mrdfits('/scr2/nichal/workspace3/flux_calib/rcs0224_Feige25.fits',1)
fluxarcs = flux_convert(scienceall,param,name,namestar)

;get magnitude from photometry

photfile = '/scr2/nichal/deimos/August15_masks/hst_09135_02_wfpc2_f814w_wf_sexphot.cat'
photall = read_ascii(photfile,comment_symbol='#')
photall = photall.field001
phot = {id:reform(photall[4,*]),xcenter:reform(photall[0,*]),ycenter:reform(photall[1,*]),ra:reform(photall[2,*]),dec:reform(photall[3,*]),magauto:reform(photall[30,*]),magerr:reform(photall[31,*])}
;'spec1d.rcs0224.007.arc1.013016_try2.matched.fits','spec1d.rcs0224.007.arc1.013116_try2.fits','spec1d.rcs0224.007.arc2.013016_try2.matched.fits','spec1d.rcs0224.007.arc2.013116_try2.fits','spec1d.rcs0224.008.arc3.013016_try2.matched.fits','spec1d.rcs0224.008.arc3.013116_try2.fits'
photmag = fltarr(6)
photmagerr = fltarr(6)
arc1_id = [2920,2857]
arc2_id = [2730,2608]
arc3_id = [2523,2540]
arcids = [[arc1_id],[arc2_id],[arc3_id]]
for i=0,n_Elements(namearcs)-1 do begin
   loc1 = where(phot.id eq arcids[0,i],cloc1)
   if cloc1 ne 1 then stop,'Find 0 or more than 1 component for this ID'   
   loc2 = where(phot.id eq arcids[1,i],cloc2)
   if cloc2 ne 1 then stop,'Find 0 or more than 1 component for this ID'
   mtot = -2.5*alog10(10^(phot.magauto(loc1)/(-2.5))+10^(phot.magauto(loc2)/(-2.5)))
   mtoterr = sqrt((phot.magerr(loc1))^2+(phot.magerr(loc2))^2)
   photmag([i*2,i*2+1])=mtot
   photmagerr([i*2,i*2+1])=mtoterr
endfor
photmagwl = [8186-1250,8186+1250] 
;fix the magnitude with photometry average
fluxarcsreal = fixphot(fluxarcs,photmag,photmagerr,photmagwl,name,namestar)
for i =0,2 do begin
   fluxarctot = combineflux(fluxarcsreal([i*2,i*2+1]),/plot)
   results = caluvslope(fluxarctot,science.zspec,name,namearcs[i],namestar)
endfor
stop
end
