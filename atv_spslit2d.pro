pro atv_spslit2d,slitnum,color

;have to cd into that directory
;color is 'B' or 'R'
  file = file_search('spSlit.*'+slitnum+color+'.{fits,fits.gz}',count=c)
  if c eq 0 then stop, 'no file found'
  twod = mrdfits(file[0],1)
  flux = twod.flux
  atv,flux

end
