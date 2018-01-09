pro shaded_uncertainties,x,y,yivar,color=color
;remove data outside range
  xmax =!x.crange[1]
  toolarge = where(x gt xmax,ctoolarge)
  if ctoolarge gt 0 then x=x[0:toolarge[0]-1]
  if ctoolarge gt 0 then y=y[0:toolarge[0]-1]
  if ctoolarge gt 0 then yivar=yivar[0:toolarge[0]-1]
;make finer grids for shaded uncertainties
  nfine = n_elements(x)*2-1
  x_fine = dblarr(nfine)
  uppery    = dblarr(nfine)
  lowery    = dblarr(nfine)
  x_fine[0] = x[0]
  for i=0,nfine-1 do begin
     if i mod 2 eq 0 then begin
        x_fine(i) = x_fine(i-1>0)
        uppery[i]    = y(i/2)+0.5*sqrt(1./yivar(i/2))
        lowery[i]    = y(i/2)-0.5*sqrt(1./yivar(i/2))
     endif
     if i mod 2 eq 1 then begin
        x_fine(i) = 0.5*(x((i-1)/2)+x((i+1)/2))
        uppery[i]    = uppery[i-1]
        lowery[i]    = lowery[i-1]
     endif
  endfor
  if keyword_Set(color) then cnow=color else cnow='rose'
  cgcolorfill,[x_fine,reverse(x_fine),x_fine[0]],[uppery,reverse(lowery),uppery[0]],color=fsc_color(cnow)

end
