pro runsfr_j1621

;initial input
name = 'j1621'
file = '/scr2/nichal/workspace3/SPEC1D/spec1d.tuck.011.j1621_fcal.fits'
namearcs= 'j1621'
namestar = ['']
redshift = 4.1278
;multfac = 1.d-16
multfac = 1.d-8

;read data
science  = mrdfits(file,1)

;Do the conversion to science
;fluxarcs = {lambda:science.lambda,flux:science.flux*multfac,fluxerr:multfac/sqrt(science.ivar)}
fluxarcs = {lambda:science.lambda,flux:science.flux*multfac/science.lambda^2,fluxerr:multfac/science.lambda^2/sqrt(science.ivar)}
photmag = 22.6           ;ABmag
photmagerr = 0.05
photmagwl  = [8186.-1250.,8186.+1250.] ;95% cummulative throughput width for F814W

;fix the magnitude with photometry average
fluxarcsreal = fixphot(fluxarcs,photmag,photmagerr,photmagwl,name,namestar)
fluxarctot = fluxarcsreal

results = caluvslope(fluxarctot,redshift,name,namearcs,namestar)

stop
end
