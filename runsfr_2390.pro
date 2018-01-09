pro runsfr_2390

;initial input
name = 'a2390'
file = '/scr2/nichal/workspace3/SPEC1D/'+['spec1d.tuck.010.H3_fcal.fits','spec1d.tuck.009.H5_fcal.fits']
namearcs=['H3','H5']
namestar = ['']
redshift = [4.043,4.0448]
multfac = [1.d-10,1.d-12]
for i=0,1 do begin
;read data
   science  = mrdfits(file[i],1)

;Do the conversion to science
   fluxarcs = {lambda:science.lambda,flux:science.flux*multfac[i]/science.lambda^2,fluxerr:multfac[i]/science.lambda^2/sqrt(science.ivar)}
   photmag = [22.4,22.8]        ;ABmag
   photmagerr = [0.01,0.01]
   photmagwl  = [8186.-1250.,8186.+1250.] ;95% cummulative throughput width for F814W

;fix the magnitude with photometry average
   fluxarcsreal = fixphot(fluxarcs,photmag,photmagerr,photmagwl,name,namestar)
   fluxarctot = fluxarcsreal

   results = caluvslope(fluxarctot,redshift[i],name,namearcs[i],namestar)
stop
endfor
stop
end
