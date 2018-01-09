function median_smooth,xin,yin,xout,sigma
width = round(2.355*sigma/(-1.*ts_diff(xin,1)))
nan = where(finite(width) eq 0,c)
if c ne 0 then width(nan) = median(nan)
evenpos = where(width mod 2 eq 1, c)
if c ne 0 then width(evenpos) = width(evenpos)+1

npix = n_elements(xin)
hwidth = (width-1)/2
yout = dblarr(npix)
for i=0, npix-1 do begin
   posmin = (i-hwidth(i))>0
   posmax = (i+hwidth(i))<(npix-1)
   yout(i) = median(yin[posmin:posmax])
endfor

if total(abs(xin-xout)) ne 0 then yout = spline(xin,yout,xout,/double)  
return,yout
end
