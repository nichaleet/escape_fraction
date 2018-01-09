pro plotmask,lambda,mask

tdiff = ts_diff(mask,1)
maskbegin = where(tdiff eq -1.,nbegin)
maskend   = where(tdiff eq 1.,nend)
if abs(nbegin-nend) gt 1 then stop
if nbegin-nend eq 1 then maskend=[maskend,n_elements(lambda)-1]
yrange = !y.crange
yrangenew = [mean(yrange)-(yrange(1)-yrange(0))*0.4,mean(yrange)+(yrange(1)-yrange(0))*0.4]
yrange = yrangenew
for i =0,nbegin-1 do begin
   polyfill,[lambda(maskbegin(i)),lambda(maskbegin(i)),lambda(maskend(i)),lambda(maskend(i))],[yrange(0),yrange(1),yrange(1),yrange(0)],color=fsc_color('light cyan')
endfor

end
