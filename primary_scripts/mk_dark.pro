pro mk_dark,  obsen,  cyclen

; Makes a dark for an IMaX observing run
; specificed by obsen and cyclen. These
; numbers are normally obtained by
; inspection of the corresponding set_hdr.save
; f2 mechanism tells when darks are obtained.
; But other options (close curtain) can be used as well.

; Output is an IDL saveset whose name must be changed
; according to input

; obsen and cyclen are obtained from inspection of int
; arrays from the set_hdr.save or from set_hdr.txt files

; This version is prepared to remove frame and 
; 936 images (or any mod 4) instead of 1024 !!!!!

; IDL>mk_dark, 168, [12, 87] 

; all these must be changed manually !!!!

; frame extraction

fx   = 49
fy   = 984
ndim = fy - fx + 1 ; all IMaX images will be 936 pixels (mod 4 is needed)
print, fx, fy

; camera 2 is flipped below,  before extracting frame 

;

;path  = '/data/sunrise/2009/IMaX/level0.0/2009_06_12/2009_06_12_18/'  ; data path
;path  = '/data/waller/IMaX.imax_pros_011110/imax8910hr'
;pathr = '/data/waller/IMaX/imax_pros_011110/sets_index/'
path = 'tempy/'

;restore, pathr + 'ind_12thJune91hr.sav'                ; index file from image_viewer_str
;restore, pathr + 'set_12thJune91hr.sav'                ; output file from image_viewer_str
restore, path + 'ind_12thJune_wdtest.sav'
restore,'12th_day_test.save'
; substitute path files to path variable above

;for ind = 0L, nima1-1L do begin
;	nlen = strlen(imas1[ind])
;	npos = strpos(imas1[ind], '2009')
;	imas1[ind] = path+strmid(imas1[ind], npos, nlen-npos)	
;endfor

;for ind = 0L, nima2-1L do begin
;	nlen = strlen(imas2[ind])
;	npos = strpos(imas2[ind], '2009')
;	imas2[ind] = path+strmid(imas2[ind], npos, nlen-npos)	
;endfor

; now we continue as usual

;outfilename = '/data/waller/IMaX/imax_pros_011110/darks/dar_12thJune91hr.sav'      ; idl save set with darks. Output file
outfilename = 'dar_12thJune_tr.sav'
; initialize darks

dc1  = fltarr(ss1, ss2)
dc2  = fltarr(ss1, ss2)
ndc1 = 0.
ndc2 = 0.

;  houskeeping data

dd1   = where(numobse1 eq obsen, nobsen1)
dd2   = where(numobse2 eq obsen, nobsen2)
nacd1 = mean(nacc1[dd1])
nacd2 = mean(nacc2[dd2])

print, 'Number of observations of Cam1 for this observing sequenece = ', nobsen1
print, 'Number of observations of Cam2 for this observing sequenece = ', nobsen2

print, '# of Accumulations is = ', mean(nacc1[dd1])  , mean(nacc2[dd2])
print, '# of Wavelengths   is = ', mean(nwavel1[dd1]), mean(nwavel2[dd2])
print, '# of Polarizations is = ', mean(npol1[dd1])  , mean(npol2[dd2])
wait, 1

cycl_0 = cyclen[0]
cycl_1 = cyclen[1]

for ind = cycl_0, cycl_1 do begin
	print, '  Image = ', ind
        dd1 = where(numobse1 eq obsen and numcycl1 eq ind, ndd1)
        print, '..................................Cycle '+ $
        strcompress(string(ind), /remove_all) +' has '+ $
        strcompress(string(ndd1), /remove_all) +' images for Cam1'
        for jnd = 0, ndd1 - 1 do begin
                ima  = float(readfits(imas1[dd1[jnd]]))                            ; Read in what we hope is a dark
;		tvwin, rebin(ima, 512, 512)                                        ; Plot it
		int1 = mean(ima[ss1/2-ss1/4:ss1/2+ss1/4, ss2/2-ss2/4:ss2/2+ss2/4]) ; Take centre of image for mean	
		print, 'Mean is  = ', int1
		if (int1 gt 50.*(nacd1-1) and int1 lt 50.*(nacd1+1)) then begin
			print, 'Image added to dark'
			dc1  = dc1 + ima
			ndc1 = ndc1 + 1.
		endif
	endfor
endfor

dc1 = dc1/ndc1

; Camera 2

for ind = cycl_0, cycl_1 do begin
	print, '  Image = ', ind
        dd2 = where(numobse2 eq obsen and numcycl2 eq ind, ndd2)
        print, '..................................Cycle '+ $
        strcompress(string(ind), /remove_all)+' has '+ $
        strcompress(string(ndd2), /remove_all)+' images for Cam2'
        for jnd = 0, ndd2-1 do begin
                ima = float(readfits(imas2[dd2[jnd]]))
;		tvwin, rebin(ima, 512, 512)
		int2 = mean(ima[ss1/2-ss1/4:ss1/2+ss1/4, ss2/2-ss2/4:ss2/2+ss2/4])	
		print, 'Mean is  = ', int2
		if (int2 gt 50.*(nacd2-1) and int2 lt 50.*(nacd2+3)) then begin ; used to be +1
			print, 'Image added to dark'
			dc2 = dc2+float(ima)
			ndc2 = ndc2+1.
		endif
	endfor
endfor

dc2 = dc2/ndc2

print, 'Images added for Cam1 & Cam2 = ', fix(ndc1), fix(ndc2)

; extracting frame

dc1 = dc1[fx:fy, fx:fy]
dc2 = rotate(dc2, 5)
dc2 = dc2[fx:fy, fx:fy]

save, dc1, dc2, ndc1, ndc2, nacd1, nacd2, fx, fy, ndim, filename = outfilename 

end
