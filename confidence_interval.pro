function confidence_interval,x,prob_x,noplot=noplot
;return peak,first q,median, third q
valx = dblarr(4)
numx = n_elements(x)

;peak
peakval = max(prob_x,pkloc)
valx[0] = x[pkloc]

;make CDF (cumulative distribution function)
CDF = dblarr(numx)
for i=1,numx-1 do CDF[i]=int_tabulated(x[0:i],prob_x[0:i]) 
;if ~keyword_set(noplot) then plot,x,cdf,ytitle='CDF'
valx[1:3] = interpol(x,cdf,[0.16,0.5,0.84])

return,valx
end
