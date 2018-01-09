pro runsfr_0940

;initial input
name = 'macs0940'
filestar = '/scr2/nichal/deimos/rawdata/2015mar/m0940star2/spec1d.Long1.5.007.m0940star2.fits'
filestarref  = 'standard_stars/feige98.dat'
namestar = 'Feige98'
file = '/scr2/nichal/workspace3/SCIENCE/'+name+'_science.fits'
fileall = '/scr2/nichal/workspace3/SCIENCE/'+name+'_scienceall.fits'
namearcs = ['arcA','arcB']
;read data
science  = mrdfits(file,1,hdr)
scienceall  = mrdfits(fileall,1,hdr)
starblu     = mrdfits(filestar,1,hdr)
starred     = mrdfits(filestar,2,hdr)
readcol,filestarref,reflam,refmag,what

;get the conversion to convert the DEIMOS obs flux to magnitude
flux_conversion,starblu,starred,reflam,refmag,name,namestar 

;Do the conversion to science
param = mrdfits('/scr2/nichal/workspace3/flux_calib/macs0940_Feige98.fits',1)
fluxarcs = flux_convert(scienceall,param,name,namestar)

;get magnitude from photometry

photfile = '/scr2/nichal/deimos/Mar15_masks/macs0940/hst_11103_35_wfpc2_multiwave_wf_sexphot_trm.cat'
photall = read_ascii(photfile,comment_symbol='#')
photall = photall.field1
phot = {id:reform(photall[0,*]),xcenter:reform(photall[1,*]),ycenter:reform(photall[2,*]),ra:reform(photall[3,*]),dec:reform(photall[4,*]),magauto:reform(photall[5,*]),magerr:reform(photall[5,*])*0.}
fluxfile = '/scr2/nichal/deimos/Mar15_masks/macs0940/se_check_files/ck_objects_11103_35_f606w.fits'
im = mrdfits(fluxfile,0) 
segmentfile = '/scr2/nichal/deimos/Mar15_masks/macs0940/se_check_files/ck_segmentation_11103_35_f606w.fits'
segment = mrdfits(segmentfile,0)
conversion = cal_flux_conversion(phot,segment,im)
imsize = size(im)
xarr = rebin(indgen(imsize[1]),imsize[1],imsize[2])
yarr = rebin(transpose(indgen(imsize[2])),imsize[1],imsize[2])
xcorners = [556,565,703,712] 
ycorners = [640,625,725,708]
ifinslit     = inside(xarr,yarr,xcorners,ycorners)
photmag = [0.,0.]
for i=0,n_Elements(namearcs)-1 do begin
   if i eq 0 then ifsegment = where(segment eq 2165 or segment eq 2177 or  segment eq 2058 or  segment eq 2133 or segment eq 2136)
   if i eq 1 then ifsegment = where(segment eq 2176 or  segment eq 2221 or segment eq 2240)
   ifinsegment  = ifinslit*0.
   ifinsegment(ifsegment) = 1
   arcregion = where(ifinslit eq 1 and ifinsegment eq 1,carcpix)
   if carcpix eq 0 then stop,'Cannot find the arc area'
   fluxhstarc  = total(im(arcregion))
   fnuarc      = conversion(0)*fluxhstarc
   photmag(i)  = -2.5*alog10(fnuarc)-48.6
   print, carcpix, 'pixels in ',namearcs(i)
endfor
photmagerr = [0.1,0.1]
photmagwl  = [4670.,7349.] ;minmax for WFPC2 F606w

;fix the magnitude with photometry average
fluxarcsreal = fixphot(fluxarcs,photmag,photmagerr,photmagwl,name,namestar)
;fluxarcsreal = fluxarcs
fluxarctot = combineflux(fluxarcsreal,/plot)
save,fluxarcs,namearcs,fluxarctot,filename='/scr2/nichal/workspace3/flux_calib/m0940flux.sav'

for i =0,1 do begin

   fluxarctot = fluxarcsreal[i] 

;cal UV slope(beta) and flux at 1600A and SFR
   if i eq 0 then results = caluvslope(fluxarctot,science.zspec,name,namearcs[i],namestar)
   if i eq 1 then results = caluvslope(fluxarctot,science.zspec,name,namearcs[i],namestar,wlmin=1400.);,fixbeta=[-1.4521513,0.12335125])
;;;;;;;;;FIX SLOPE TO ARC B! MAYBE CONTAMINATION???
endfor
stop
end
