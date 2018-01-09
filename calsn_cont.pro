function calsn_cont, x,y,ivar,dx

yerr = 1./sqrt(ivar)
nx = n_elements(x)
minx = min(x)
maxx = max(x)
sn   = dblarr(nx)

for i = 0,nx-1 do begin
 minbox = (x(i)-dx(i)/2.) > minx
 maxbox = (x(i)+dx(i)/2.) < maxx
 good   = where(x ge minbox and x le maxbox, cgood)
 if cgood ge 2 then begin
    yin = y(good)
    yerrin = yerr(good)
    meanerr,yin,yerrin,sigmean,sigmam,sigmad,sigmas
    signal  = sigmean
    noise   = 1./sqrt(total(ivar(good)))
    sn(i)   = signal/noise
 endif
endfor

return,sn
end
