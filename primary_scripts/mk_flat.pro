pro mk_flat, obsen, cyclen, outfilename

; Is line mk_flat_str adapted for the longitudinal cases
; that uses n_pol=4 when it should have been n_pol=2

;IDL>mk_flat_str,222,[1,52],'flats_DATA_10_set_l12-2.save'

; obsen and cyclen are guessed from image_viewer_str.pro and
; the ouput saveset. Output file is outfilename and a .txt log file


; open output log file

openw, 1, strmid(outfilename, 0, strpos(outfilename,'.s'))+'.txt'

; path for darks. Needs to be changed manually !!!

;pathd = 'darks/'
;restore , pathd + 'dar_12thJune91hr.sav', /v
restore, 'dar_12thJune_tr.sav',/v
print   ,    'Images are of size=', ndim, ndim
printf  , 1, 'Images are of size=', ndim, ndim

;data path. Needs to be changed manually !!!

;path  = '/data/sunrise/2009/IMaX/level0.0/2009_06_12/2009_06_12_18/' ; data path
;pathr = 'sets_index/' ;path to set and index

;restore, pathr + 'ind_12thJune91hr.sav', /v      ; index file from image_viewer_str
;restore, pathr + 'set_12thJune91hr.sav', /v                ; output file from image_viewer_str
restore,'tempy/ind_12thJune_wdtest.sav',/v
restore, '12th_day_test.save', /v
; substitute path files to path variable above

;for ind=0L,nima1-1L do begin
;        nlen=strlen(imas1[ind])
;        npos=strpos(imas1[ind],'2009')
;        imas1[ind]=path+strmid(imas1[ind],npos,nlen-npos)
;endfor

;for ind=0L,nima2-1L do begin
;        nlen=strlen(imas2[ind])
;        npos=strpos(imas2[ind],'2009')
;        imas2[ind]=path+strmid(imas2[ind],npos,nlen-npos)
;endfor

; Identify indexes for corresponding observing sequence

dd1 = where(numobse1 eq obsen, nobsen1)
;;print,'Number of observations of Cam1 for this observing sequenece=',nobsen1
;;printf,1,'Number of observations of Cam1 for this observing sequenece=',nobsen1
dd2 = where(numobse2 eq obsen, nobsen2)
;;print,'Number of observations of Cam2 for this observing sequenece=',nobsen2
;;printf,1,'Number of observations of Cam2 for this observing sequenece=',nobsen2
;;print,'# of Accumulations is=',mean(nacc1[dd1]),mean(nacc2[dd2])
;;printf,1,'# of Accumulations is=',mean(nacc1[dd1]),mean(nacc2[dd2])
;;print,'# of Wavelengths is=',mean(nwavel1[dd1]),mean(nwavel2[dd2])
;;printf,1,'# of Wavelengths is=',mean(nwavel1[dd1]),mean(nwavel2[dd2])
;;print,'# of Polarizations is=',mean(npol1[dd1]),mean(npol2[dd2])
;;printf,1,'# of Polarizations is=',mean(npol1[dd1]),mean(npol2[dd2])

;;n_accum  = reform(nacc1[dd1[0]])
;;n_accum  = fix(n_accum[0])
n_accum  = nacc1[dd1[0]];why reform a float? Des this even need to be an int?
;;n_lambda = reform(nwavel1[dd1[0]])
;;n_lambda = fix(n_lambda[0])
n_lambda = nwavel1[dd1[0]]
;;n_pol    = reform(npol1[dd1[0]])
;;n_pol    = fix(n_pol[0])
n_pol    = npol1[dd1[0]]

cycl_0   = cyclen[0]
cycl_1   = cyclen[1]
; inititalize flats

ff1 = fltarr(ndim, ndim, n_lambda, n_pol)
ff2 = fltarr(ndim, ndim, n_lambda, n_pol)
nf1 = fltarr(n_lambda, n_pol)
nf2 = fltarr(n_lambda, n_pol)
	
for ind = cycl_0, cycl_1 do begin
	dd1 = where(numobse1 eq obsen and numcycl1 eq ind, ndd1);can ndd1 ever not be 1?
;;	print,'..................................Cycle '+ $
;;	strcompress(string(ind),/remove_all)+' has '+ $
;;	strcompress(string(ndd1),/remove_all)+' images for Cam1'
;;	printf,1,'..................................Cycle '+ $
;;	strcompress(string(ind),/remove_all)+' has '+ $
;;	strcompress(string(ndd1),/remove_all)+' images for Cam1'
	for jnd = 0, ndd1 - 1 do begin
		ima    = float(readfits(imas1[dd1[jnd]]))
		ima    = ima[fx:fy, fx:fy]
		ima    = ima - dc1 * float(n_accum) / float(nacd1)
		;tvwin, rebin(ima, ndim/2, ndim/2)
		polind = fix(polst1[dd1[jnd]] - 1)
		lamind = fix(lamst1[dd1[jnd]] - 1)
		contrast=sss1[dd1[jnd]]/int1[dd1[jnd]]
		print,'Flat with contrast =',contrast
		printf,1,'Flat with contrast =',contrast
		nf1[lamind,polind]=nf1[lamind,polind]+1;why is this only done in one place everywhere else is xero, later we divide by it and get nans. Okay if ndd1 is greater than 1 we wont have this problem so bad
		ff1[*,*,lamind,polind]=ff1[*,*,lamind,polind]+ima
	endfor
	
	dd2=where(numobse2 eq obsen and numcycl2 eq ind,ndd2)
	print,'.................................Cycle '+ $
	strcompress(string(ind),/remove_all)+' has '+ $
	strcompress(string(ndd2),/remove_all)+' images for Cam2'
	printf,1,'..................................Cycle '+ $
	strcompress(string(ind),/remove_all)+' has '+ $
	strcompress(string(ndd2),/remove_all)+' images for Cam2'
	for jnd=0,ndd2-1 do begin
		ima=float(readfits(imas2[dd2[jnd]]))
                ima=rotate(ima,5) ; flip Camera 2
		ima=ima[fx:fy,fx:fy]
		ima=ima-dc2*float(n_accum)/float(nacd2)
		;tvwin,rebin(ima,ndim/2,ndim/2)
		polind=fix(polst2[dd2[jnd]]-1)
		lamind=fix(lamst2[dd2[jnd]]-1)
		contrast=sss2[dd2[jnd]]/int2[dd2[jnd]]
		print,'Flat with contrast =',contrast
		printf,1,'Flat with contrast =',contrast
		nf2[lamind,polind]=nf2[lamind,polind]+1
		ff2[*,*,lamind,polind]=ff2[*,*,lamind,polind]+ima
	endfor
endfor
print,'Files added for Camera 1=',nf1
printf,1,'Files added for Camera 1=',nf1
print,'Files added for Camera 2=',nf2
printf,1,'Files added for Camera 2=',nf2
for ind=0,n_lambda-1 do begin;NaNs appear here
	for jnd=0,n_pol-1 do begin
		ff1[*,*,ind,jnd]=ff1[*,*,ind,jnd]/float(nf1[ind,jnd])
		ff2[*,*,ind,jnd]=ff2[*,*,ind,jnd]/float(nf2[ind,jnd])
	endfor
endfor

save,ff1,ff2,nf1,nf2,obsen,cyclen,n_lambda,n_pol,n_accum,filename=outfilename 

close,1

end
