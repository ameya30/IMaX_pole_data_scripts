;Program: "prepare_deconvol.pro"

;**********************************************************************

@lib_wiener
@mean_power2
@masking_power
@filtering_fringes
@rubbish

;*********************************************************************************

pro prepare_deconvol,obsen,cyclen,save_in=save_in,save_out=save_out

; reading images, dark substraction, prefilter correction and
; flat fielding as in mk_magneto.pro

; obsen and cyclen have same meaninig as in mk_magneto. 
; Use obsen and cyclen = 0 if you are using the save_in 
; option to change any previously generated output of prepare_deconvol
; prepare_deconvol output is save_out

;IDL>prepare_deconvol,225,100,save_out='masks_225_100.save'

;Inputs through editor !!!!!!!!!!!

;Lower and upper (x,y) coordinates defining the image-box to be estracted.
;Si se ponen todos a 0 entonces toma la imagen completa
;lox=54 & hix=54+951 & loy=40 & hiy=40+951   
;lox=0 & hix=0 & loy=0 & hiy=0   ;para trabajar con la imagen completa

;_________________________________________________________________________

  ;Option: we may want to improve pre-existing masks included in "save_in"
  if keyword_set(save_in) ne 0 then begin
    print,'Recovering masks from the file '+save_in
    restore,save_in
  endif else begin

; images readout (from Level 0 fits !!)

if (obsen ne 0 and cyclen ne 0) then begin ; image readout

;path  = '/data/sunrise/2009/IMaX/level0.0/2009_06_12/2009_06_12_16/' ; path to data to restore
pathr = 'sets_index/'   ;path to sets and inds 
pathd = 'darks/'        ;path to darks
pathf = 'flats_stage4/' ;path to rem_lin flats

restore, pathd + 'dar_12thJune91hr.sav', /v
restore, pathf + 'stage4_flat_out.sav', /v           ;new flats from rem_lin_flat
restore, pathr + 'set_12thJune16hr.sav', /v           ;restore housekeepings of data not dark or flats?
restore, pathf + 'pref_factors_stage4_flat_out.sav', /v ;restore prefilter factors

;; overwrite ff1 and ff2 with flats that do not have low freq fringes
;ff1 = readfits('cam1_flat_full_no_med.fits')
;ff2 = readfits('cam2_flat_full_no_med.fits')

n_lambda = fix(n_lambda)
n_pol    = fix(n_pol)

;some index

dd1=where(cyclen eq numcycl1 and obsen eq numobse1,ndd1)
dd2=where(cyclen eq numcycl2 and obsen eq numobse2,ndd2)

; prefilter correction of housekeeping data

for ind=0,ndd1-1 do int1[dd1[ind]]=int1[dd1[ind]]/pref_fact[lamst1[dd1[ind]]-1]
for ind=0,ndd2-1 do int2[dd2[ind]]=int2[dd2[ind]]/pref_fact[lamst2[dd2[ind]]-1]

; plot intensities, should look like the spectral line

plot,time1[dd1],int1[dd1],psym=1,xtitle='Time obs.',title='Mean Intensities of data set'    
oplot,time2[dd2],int2[dd2]-1000,psym=2

; index data directory

restore, pathr + 'ind_12thJune16hr.sav'

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

; find for cycle cyclen all frames read them, flip Cam2, extract images,
; substract dc, flatfield data and correct prefilter. 

arr1=fltarr(ndim,ndim,n_lambda,n_pol)
arr2=fltarr(ndim,ndim,n_lambda,n_pol)

indm1=fltarr(n_lambda,n_pol)
indm2=fltarr(n_lambda,n_pol)

if ndd1 ne n_pol*n_lambda or ndd2 ne n_pol*n_lambda then begin 
	print,'Cycle not complete'
	close,1
	return
endif

; read files

for ind=0,n_lambda-1 do begin
	for jnd=0,n_pol-1 do begin
		indm1[ind,jnd]=where(numobse1 eq obsen and numcycl1 eq cyclen and lamst1 eq  $
		ind+1 and polst1 eq jnd+1,ncheck1)
		indm2[ind,jnd]=where(numobse2 eq obsen and numcycl2 eq cyclen and lamst2 eq  $
		ind+1 and polst2 eq jnd+1,ncheck2)
		if (ncheck1 ne 1) then begin
			print,'Image not found' 
			return
		endif else begin 
			ima=readfits(imas1[indm1[ind,jnd]])
			ima=ima[fx:fy,fx:fy]
			ima=ima-dc1*float(n_accum)/float(nacd1) ; dark substraction
			ima=ima/(ff1[*,*,ind,jnd]*pref_fact[ind])	     ; flat and prefilter
			arr1[*,*,ind,jnd]=ima	
			tvwin,rebin(ima,ndim/2,ndim/2)
			meanima=mean(ima)
			print,'Cam1 image contrast is=', $
			stdev(ima)/meanima
		endelse
		if (ncheck2 ne 1) then begin
			print,'Image not found' 
			return
		endif else begin 
			ima=readfits(imas2[indm2[ind,jnd]])
			ima=rotate(ima,5) ; flip Camera 2
			ima=ima[fx:fy,fx:fy]
			ima=ima-dc2*float(n_accum)/float(nacd2)
			ima=ima/(ff2[*,*,ind,jnd]*pref_fact[ind])	
			arr2[*,*,ind,jnd]=ima
			tvwin,rebin(ima,ndim/2,ndim/2)
			meanima=mean(ima)
			print,'Cam2 image contrast is=', $
			stdev(ima)/meanima 
		endelse
	endfor
endfor


endif

;Construct the average of normalized images and of their respective power spectra
;(for each camera).
    perc=12.5  ;10.
    print,''  &  print,'RUNNING PROGRAM mean_power2.pro'
    print,'CAMERA 1'
    mean_power2,arr1,perc,im1,pow1,poww
    print,'CAMERA 2'
    mean_power2,arr2,perc,im2,pow2,poww
  endelse

writefits, 'cam2_line160_meanimage.fits', im2

  ;Based on the mean power spectrum, construct the mask defining the areas including
  ;spurious signal (fringes) (for each camera). Also construct the weighting function
  ;to filter out this spurious signal.
  print,''  &  print,'RUNNING PROGRAM masking_power.pro'
  print,'INTERACTIVE PROCESS WITH THE MOUSE: "rubber tool"'
  print,'CAMERA 1'  
  print,'Hit Enter....' & PAUSE
  masking_power,pow1,powcorr1,mask1,weight1
  ; saving temporary arrays in case the program crashes later
  save,im1,im2,pow1,pow2,powcorr1,mask1,weight1,filename='prepare_deconvol_tmp.save'
  print,''  &  print,'CAMERA 2'  
  print,'Hit Enter....' & PAUSE
  masking_power,pow2,powcorr2,mask2,weight2
  save,im1,im2,pow1,pow2,powcorr1,powcorr2,mask1,mask2,weight1,weight2,filename='prepare_deconvol_tmp.save'

  ;Remove fringes from the mean image (for each camera)
  ;Lower and upper thresholds are required through display (reference: -0.01, 0.01)
  print,''  &  print,'RUNNING PROGRAM filtering_fringes.pro'
  print,'CAMERA 1'  
  print,'Hit Enter....' & PAUSE
  filtering_fringes,0,0,5,weight1,im1,imf1
  print,''  &  print,'CAMERA 2'  
  print,'Hit Enter....' & PAUSE
  filtering_fringes,0,0,5,weight2,im2,imf2
writefits, 'cam2_line185_filtered_fringes.fits', im2
  ;Construct the mask defining the position of dust-particules and scratches remaining
  ;after a defective flatfield (for each camera).
  print,''  &  print,'RUNNING PROGRAM rubbish.pro'
  print,'CAMERA 1'
  rubbish,imf1,5,grains1,maskg1,maskgd1
  print,'CAMERA 2'
  rubbish,imf2,5,grains2,maskg2,maskgd2

  ;Optional: Save all these calibration data (for both cameras together)
  print,''
  print,'Shall I save the new masks (the save file will be overwritten if it already exist)? (y/n)'
  testch=string(3)  &  testch=get_kbrd(1)
  if testch eq 'y' then begin
    print,'Writing masks into the file '+save_out
    save,im1,im2,pow1,pow2,powcorr1,powcorr2,mask1,mask2,weight1,weight2,imf1,imf2,grains1,grains2,maskg1,maskg2,maskgd1,maskgd2,filename=save_out
  print,'Deleting temporary files...'
  spawn,'rm prepare_deconvol_tmp.save'
  endif

end
