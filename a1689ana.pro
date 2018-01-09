pro a1689ana
!p.multi=[0,1,2]
path='/scr2/nichal/deimos/reduced_data/a1689/bothnights/'
files=['spec1d.a1689.007.arcz4868a.fits','spec1d.a1689.012.arcz4868b.fits','spec1d.a1689.009.arcz5120.fits']
name=['arcz487a','arcz487b','arcz512']
for nfile=0,n_elements(files)-1 do begin
window,nfile
blue=mrdfits(path+files(nfile),1)
red=mrdfits(path+files(nfile),2)
wl = [blue.lambda,red.lambda]
spec=[blue.spec,red.spec]
ivar = [blue.ivar,red.ivar]
plot,wl,spec,psym=10,title='Spectrum',yrange=[0,1000]
lya_wl=wl(where(spec eq max(spec(where(wl lt 7300)))))
if nfile eq 2 then lya_wl=wl(where(spec eq max(spec(where(wl lt 7600))))) 
z=lya_wl/1215.67-1.
restore,'line.sav'
line_z = line_wl*(z(0)+1.)
for i=0,n_elements(line_wl)-1 do oplot,[line_z(i),line_z(i)],[min(!y.crange),-500],color=255,linestyle=2
;calculate Signal to Noise
;Calculate resolution in angstrom
r=0.57 ;anamorphic factor for 1200g/mm grating at central wl 8360 is 0.57
w = 0.728 ;mm (scale at slit is 0.728 mm/arcsec) 
fcam=15. ;inches
fcoll = 86.5 ;inches
resolution = r*(fcam/fcoll)*w ;in mm at detector
resolution_pix = resolution/0.015 ;in pixel
dispersion = 0.33 ;angstrom per pix for 1200 line grating
resolution = resolution_pix*0.33 ;in Angstrom
;------------------------------------------------
;get the values per resolution element
numpix_in_res = round(resolution_pix)
n_res = floor(float(n_elements(wl))/numpix_in_res)
wl_res = fltarr(n_res)
signal = wl_res
noise   = wl_res
for ii=0,n_res-1 do begin
ind = findgen(numpix_in_res)+ii*numpix_in_res
print,ind
wl_res(ii) = mean(wl(ind))
;signal(ii)=total(spec(ind))
signal(ii) = total(spec(ind)*ivar(ind))/total(ivar(ind))
;noise(ii) = sqrt(total(1./ivar(ind)))
noise(ii)=sqrt(total(1./ivar(ind)))/float(numpix_in_res)
endfor
snr = signal/noise
;plot,wl_res,snr,psym=10,title='SNR'
plot,wl_res,snr,psym=10,title='SNR',yrange=[0,20]
wl=wl_res
save,wl,signal,noise,snr,filename=name(nfile)+'.sav'
stop
endfor
;add the two blobs of z487
restore,name(0)+'.sav'
wl_a=wl
signal_a=signal
noise_a=noise
restore,name(1)+'.sav'
wl_b=wl
signal_b=signal
noise_b=noise
;match the wavelength
start_a = where(wl_a-

end
