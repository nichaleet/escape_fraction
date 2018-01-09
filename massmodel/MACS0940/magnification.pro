pro magnification, name
;lightweighted magnification
;source plane
imsource = 'source'+name+'.fits'
ss = readfits(imsource)
ssnew = im_interpolate(ss)
writefits,'source'+name+'.new.fits',ssnew
sourceflux = total(ssnew)
goodsspix  = where(ssnew ne 0,csourcepix)

;image plane
reg    = readfits('Bw'+name+'.fits')
imim   = readfits('s606.fits')
if (size(reg))[1] lt (size(imim))[1] or (size(reg))[2] lt (size(imim))[2] then begin
   regnew = bytarr((size(imim))[1],(size(imim))[2])
   regnew[0:(size(reg))[1]-1,0:(size(reg))[2]-1] = reg
   reg = regnew
endif
good = where(reg eq 1,cimgood)
imflux = total(imim(good))
magflux = imflux/sourceflux*16.
magarea = cimgood/float(csourcepix)*16
;magnification map
mumap    = readfits('Magnification_m0940_4p04_20160412_LTM_Gc.fits')
muerrmap = readfits('Magnification_1sigmaErr_m0940_4p04_20160412_LTM_Gc.fits')
reg    = readfits('Bw'+name+'.fits')
good = where(reg eq 1,cgood)
mu = mumap(good)
muerr = muerrmap(good)
meanerr,mu,muerr,meanmag,sigmam,sigmad,sigmas
print,'magnification by area=',magarea
print,'magnification by flux=',magflux
print,'magnification by mu map(mean, median, sigmam, sigmad, sigmas)='
print,meanmag,median(mu),sigmam,sigmad,sigmas
stop
end
