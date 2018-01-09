pro check_spec1d_pos,directory
;input: dir - directory containing spec1d files
; It'll check for RA and Dec of all the identified objects and
; make a region mask
  parentdir = '/scr2/nichal/deimos/reduced_data/'
  directory = parentdir+directory
  files = file_search(directory+'*/spec1d*.{fits,fits.gz}', count=c)
  if c eq 0 then stop, 'No files found'
  objnames = strarr(c)
  slits    = strarr(c)

  for i=0, c-1 do begin
     basefile    = file_basename(files[i])
     extensions  = reverse(strsplit(basefile,'.',/extract))
     slits[i]    = extensions[2]
     objnames[i] = extensions[1]
     data1 = mrdfits(files(i),1,hdr,/silent)
     ras   = sxpar(hdr,'RA_OBJ')
     decs  = sxpar(hdr,'DEC_OBJ')
     print, ras,decs
  endfor
stop
end
