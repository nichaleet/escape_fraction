pro plotlines,z,color=color,delv=delv
  readcol, '/scr2/nichal/workspace/lines.txt',linewave,fosc,linename,format ='D,D,A,X',comment='#'   
  nlines = n_elements(linewave)
  linewave = linewave*(z+1.)
 
  for j=0,nlines-1 do if linewave(j) ge !x.crange[0] and linewave(j) lt !x.crange[1] then begin
     ymin = (!y.crange)[1]*0.65
     ymax = (!y.crange)[1]
     if keyword_Set(delv) then begin
        delwl = linewave[j]*delv/3.e5
        cgcolorfill,[linewave[j]-delwl,linewave[j]+delwl,linewave[j]+delwl,linewave[j]-delwl,linewave[j]-delwl],[ymax,ymax,ymin,ymin,ymax],color=fsc_color('grey')
     endif
     if ~keyword_set(color) then cnow = 'red' else cnow = color
     oplot,[linewave[j],linewave[j]],[ymin,ymax],color=fsc_color(cnow),linestyle=2 
     xyouts,linewave[j],((!y.crange)[1])*0.7,linename[j],color=fsc_color(cnow),font=0,charsize=0
  endif
end
