pro tifftofits, imfile

im=read_tiff(imfile)
sizep = size(im)
imnew = bytarr(sizec[1],sizec[2])
for i=0,sizec[1]-1 do imnew[*,2500-i]=im[*,i]
writefits,'stiff.fits',imnew
end
