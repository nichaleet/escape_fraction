pro shiftmedian
;arc1
a1 = mrdfits('spec1d.rcs0224.007.arc1.013116_try2.fits',1,hdr)
toshift = median(a1.spec(where(a1.lambda lt 7100.)))
if toshift gt 0 then toshift=0
a1.spec = a1.spec-toshift
print, 'arc1', toshift
mwrfits,a1,'spec1d.rcs0224.007.arc1.013116_try2_adj.fits',hdr,/create
a1 = mrdfits('spec1d.rcs0224.007.arc1.013116_try2.fits',2,hdr)
a1.spec = a1.spec-toshift
mwrfits,a1,'spec1d.rcs0224.007.arc1.013116_try2_adj.fits',hdr

a1 = mrdfits('spec1d.rcs0224.007.arc1.013116_try2.fits',3,hdr)
toshift =  median(a1.spec(where(a1.lambda lt 7100.)))
if toshift gt 0 then toshift=0
print, 'arc1', toshift
a1.spec = a1.spec-toshift
mwrfits,a1,'spec1d.rcs0224.007.arc1.013116_try2_adj.fits',hdr
a1 = mrdfits('spec1d.rcs0224.007.arc1.013116_try2.fits',4,hdr)
a1.spec = a1.spec-toshift
mwrfits,a1,'spec1d.rcs0224.007.arc1.013116_try2_adj.fits',hdr

;arc2
a1 = mrdfits('spec1d.rcs0224.007.arc2.013116_try2.fits',1,hdr)
toshift = median(a1.spec(where(a1.lambda lt 7100.)))
if toshift gt 0 then toshift=0
print, 'arc2', toshift
a1.spec = a1.spec-toshift
mwrfits,a1,'spec1d.rcs0224.007.arc2.013116_try2_adj.fits',hdr,/create
a1 = mrdfits('spec1d.rcs0224.007.arc2.013116_try2.fits',2,hdr)
a1.spec = a1.spec-toshift
mwrfits,a1,'spec1d.rcs0224.007.arc2.013116_try2_adj.fits',hdr

a1 = mrdfits('spec1d.rcs0224.007.arc2.013116_try2.fits',3,hdr)
toshift =  median(a1.spec(where(a1.lambda lt 7100.)))
if toshift gt 0 then toshift=0
print, 'arc2', toshift
a1.spec = a1.spec-toshift
mwrfits,a1,'spec1d.rcs0224.007.arc2.013116_try2_adj.fits',hdr
a1 = mrdfits('spec1d.rcs0224.007.arc2.013116_try2.fits',4,hdr)
a1.spec = a1.spec-toshift
mwrfits,a1,'spec1d.rcs0224.007.arc2.013116_try2_adj.fits',hdr

;arc3
a1 = mrdfits('spec1d.rcs0224.008.arc3.013116_try2.fits',1,hdr)
toshift = median(a1.spec(where(a1.lambda lt 7100.)))
if toshift gt 0 then toshift=0
print, 'arc3', toshift
a1.spec = a1.spec-toshift
mwrfits,a1,'spec1d.rcs0224.008.arc3.013116_try2_adj.fits',hdr,/create
a1 = mrdfits('spec1d.rcs0224.008.arc3.013116_try2.fits',2,hdr)
a1.spec = a1.spec-toshift
mwrfits,a1,'spec1d.rcs0224.008.arc3.013116_try2_adj.fits',hdr

a1 = mrdfits('spec1d.rcs0224.008.arc3.013116_try2.fits',3,hdr)
toshift =  median(a1.spec(where(a1.lambda lt 7100.)))
if toshift gt 0 then toshift=0
print, 'arc3', toshift
a1.spec = a1.spec-toshift
mwrfits,a1,'spec1d.rcs0224.008.arc3.013116_try2_adj.fits',hdr
a1 = mrdfits('spec1d.rcs0224.008.arc3.013116_try2.fits',4,hdr)
a1.spec = a1.spec-toshift
mwrfits,a1,'spec1d.rcs0224.008.arc3.013116_try2_adj.fits',hdr

;firstnight shift median to match
;arc1
a1 = mrdfits('spec1d.rcs0224.007.arc1.013116_try2_adj.fits',1)
b1 = mrdfits('spec1d.rcs0224.007.arc1.013016_try2.fits',1,hdr)
toshift = median(a1.spec)-median(b1.spec)
if toshift lt 0 then toshift=0
print, 'arc1', toshift
b1.spec = b1.spec+toshift
mwrfits,b1,'spec1d.rcs0224.007.arc1.013016_try2_adj.fits',hdr,/create
b1 = mrdfits('spec1d.rcs0224.007.arc1.013016_try2.fits',2,hdr)
b1.spec = b1.spec+toshift
mwrfits,b1,'spec1d.rcs0224.007.arc1.013016_try2_adj.fits',hdr
a1 = mrdfits('spec1d.rcs0224.007.arc1.013116_try2_adj.fits',3)
b1 = mrdfits('spec1d.rcs0224.007.arc1.013016_try2.fits',3,hdr)
toshift = median(a1.spec)-median(b1.spec)
print, 'arc1', toshift
if toshift lt 0 then toshift=0
b1.spec = b1.spec+toshift
mwrfits,b1,'spec1d.rcs0224.007.arc1.013016_try2_adj.fits',hdr
b1 = mrdfits('spec1d.rcs0224.007.arc1.013016_try2.fits',4,hdr)
b1.spec = b1.spec+toshift
mwrfits,b1,'spec1d.rcs0224.007.arc1.013016_try2_adj.fits',hdr

;arc2
a1 = mrdfits('spec1d.rcs0224.007.arc2.013116_try2_adj.fits',1)
b1 = mrdfits('spec1d.rcs0224.007.arc2.013016_try2.fits',1,hdr)
toshift = median(a1.spec)-median(b1.spec)
if toshift lt 0 then toshift=0
print, 'arc2', toshift
b1.spec = b1.spec+toshift
mwrfits,b1,'spec1d.rcs0224.007.arc2.013016_try2_adj.fits',hdr,/create
b1 = mrdfits('spec1d.rcs0224.007.arc2.013016_try2.fits',2,hdr)
b1.spec = b1.spec+toshift
mwrfits,b1,'spec1d.rcs0224.007.arc2.013016_try2_adj.fits',hdr
a1 = mrdfits('spec1d.rcs0224.007.arc2.013116_try2_adj.fits',3)
b1 = mrdfits('spec1d.rcs0224.007.arc2.013016_try2.fits',3,hdr)
toshift = median(a1.spec)-median(b1.spec)
print, 'arc2', toshift
if toshift lt 0 then toshift=0
b1.spec = b1.spec+toshift
mwrfits,b1,'spec1d.rcs0224.007.arc2.013016_try2_adj.fits',hdr
b1 = mrdfits('spec1d.rcs0224.007.arc2.013016_try2.fits',4,hdr)
b1.spec = b1.spec+toshift
mwrfits,b1,'spec1d.rcs0224.007.arc2.013016_try2_adj.fits',hdr

;arc3
a1 = mrdfits('spec1d.rcs0224.008.arc3.013116_try2_adj.fits',1)
b1 = mrdfits('spec1d.rcs0224.008.arc3.013016_try2.fits',1,hdr)
toshift = median(a1.spec)-median(b1.spec)
if toshift lt 0 then toshift=0
print, 'arc3', toshift
b1.spec = b1.spec+toshift
mwrfits,b1,'spec1d.rcs0224.008.arc3.013016_try2_adj.fits',hdr,/create
b1 = mrdfits('spec1d.rcs0224.008.arc3.013016_try2.fits',2,hdr)
b1.spec = b1.spec+toshift
mwrfits,b1,'spec1d.rcs0224.008.arc3.013016_try2_adj.fits',hdr
a1 = mrdfits('spec1d.rcs0224.008.arc3.013116_try2_adj.fits',3)
b1 = mrdfits('spec1d.rcs0224.008.arc3.013016_try2.fits',3,hdr)
toshift = median(a1.spec)-median(b1.spec)
if toshift lt 0 then toshift=0
print, 'arc3', toshift
b1.spec = b1.spec+toshift
mwrfits,b1,'spec1d.rcs0224.008.arc3.013016_try2_adj.fits',hdr
b1 = mrdfits('spec1d.rcs0224.008.arc3.013016_try2.fits',4,hdr)
b1.spec = b1.spec+toshift
mwrfits,b1,'spec1d.rcs0224.008.arc3.013016_try2_adj.fits',hdr

end
