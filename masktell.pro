pro masktell
  readcol,'/scr2/nichal/workspace3/telluric.mask', tellstart, tellend, format='D,D', /silent, comment='#'
  yrange = !y.crange
  ymin   = yrange[0]
  ymax   = yrange[1]
  for i=0,n_elements(tellstart)-1 do polyfill,[tellstart(i),tellend(i),tellend(i),tellstart(i)],[ymax,ymax,ymin,ymin],color=fsc_color('cyan'),/data
end
