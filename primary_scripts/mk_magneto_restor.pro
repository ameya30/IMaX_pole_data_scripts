pro mk_magneto_restor,obsen,cyclen,maskfile,longitudinal=longitudinal,file_out

; obsen is observation number of IMaX logs
; cyclen is cycle to be analyzed

; define path to data and restore headers, flats and darks

;IDL>mk_magneto_restor,225,100,'masks_225_100.save',longitudinal=0 

; longitudinal=1 should be used for L12-2 or L3-2

path_out='saves_Oct11/' + file_out
if n_elements(longitudinal) eq 0 then longitudinal=0

path  = '/data/sunrise/2009/IMaX/level0.0/2009_06_12/2009_06_12_16/' ; path to data to restore
;pathr = 'sets_index/'   ;path to sets and inds
pathr = '/scratch/prabhu/HollyWaller/IMaX_pole_data_scripts/sets_inds/' 
;pathd = 'darks/'        ;path to darks
;pathf = 'flats_stage4/' ;path to rem_lin flats
pathf = 'at3/'
;restore, pathd + 'dar_12thJune91hr.sav', /v
restore, 'dar_12thJune_tr.sav', /v
;restore, pathf + 'srage4_flat_out.sav', /v
restore, pathf + 'stage4_flat_out_at3.sav', /v               ;new flats from rem_lin_flat
restore, pathr + 'set_12thJune16hr.sav', /v               ;restore housekeepings of data not dark or flats?
;restore, pathf + 'pref_factors_stage4_flat_out.sav', /v ;restore prefilter factors
restore, pathf + 'pref_factors_stage4_flat_out_at3.sav', /v ;restore prefilter factors

;; overwrite ff1 and ff2 with flats that do not have low freq fringes
;ff1 = readfits('cam1_flat_full_no_med.fits')
;ff2 = readfits('cam2_flat_full_no_med.fits')

n_lambda=fix(n_lambda)
n_pol=fix(n_pol)


; prefilter correction (3% for n_lambda 12)

precor=0.97
if (n_lambda eq 12) then pref_fact=pref_fact*poly(findgen(12),[precor,(1.-precor)/11.]) 

;some index

ss1h=ndim/2
ss1q=ndim/4
ss2h=ndim/2
ss2q=ndim/4

dd1=where(cyclen eq numcycl1 and obsen eq numobse1,ndd1)
dd2=where(cyclen eq numcycl2 and obsen eq numobse2,ndd2)

dd1 = dd1[0:19]
dd2 = dd2[0:19]

ndd1 = n_elements(dd1)
ndd2 = n_elements(dd2)

; prefilter correction of housekeeping data


for ind=0,ndd1-1 do int1[dd1[ind]]=int1[dd1[ind]]/pref_fact[lamst1[dd1[ind]]-1]
for ind=0,ndd2-1 do int2[dd2[ind]]=int2[dd2[ind]]/pref_fact[lamst2[dd2[ind]]-1]

; plot intensities, should look like the spectral line

plot,time1[dd1],int1[dd1],psym=1,xtitle='Time obs.',title='Mean Intensities of data set'    
oplot,time2[dd2],int2[dd2]-1000,psym=2

; index data directory

;imas1=file_search(path,'*Camera_1*.fits',count=nima1)
;imas2=file_search(path,'*Camera_2*.fits',count=nima2)
;restore, pathr + 'ind_12thJune16hr.sav'
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

iim1=fltarr(ndim,ndim,n_lambda,n_pol)  ; restored
iim2=fltarr(ndim,ndim,n_lambda,n_pol)  ; restored
iim1n=fltarr(ndim,ndim,n_lambda,n_pol) ; not restored
iim2n=fltarr(ndim,ndim,n_lambda,n_pol) ; not restored

indm1=fltarr(n_lambda,n_pol)
indm2=fltarr(n_lambda,n_pol)

; text file with log data from reduction procedure

textfilename='reduc_rnr_'+strcompress(string(fix(obsen)),/remove_all)+ $
'_'+strcompress(string(fix(cyclen)),/remove_all)+'.txt'

openw,1,textfilename

if ndd1 ne n_pol*n_lambda or ndd2 ne n_pol*n_lambda then begin 
	print,'Cycle not complete'
	printf,1,'Cycle not complete'
	close,1
	return
endif

; read files

date_arraycam1 = make_array(20, /string) ; Added this 31st May 2017
date_arraycam2 = make_array(20, /string) ; Added this 31st May 2017
count1 = 0
count2 = 0

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
			ima=readfits(imas1[indm1[ind,jnd]], cam1head)

			date_arraycam1[count1] = cam1head[11] ; Added 31st May 2017
			count1 = count1 + 1

			ima=ima[fx:fy,fx:fy]
			ima=ima-dc1*float(n_accum)/float(nacd1) ; dark substraction
			ima=ima/(ff1[*,*,ind,jnd]*pref_fact[ind])	     ; flat and prefilter
			iim1[*,*,ind,jnd]=ima	
			tvwin,rebin(ima,ndim/2,ndim/2)
			meanima=mean(ima)
			print,'Cam1 image contrast is=', $
			stdev(ima)/meanima
			printf,1,'Cam1 image contrast is=', $
			stdev(ima)/meanima
		endelse
		if (ncheck2 ne 1) then begin
			print,'Image not found' 
			return
		endif else begin 
			ima=readfits(imas2[indm2[ind,jnd]], cam2head)

			date_arraycam2[count2] = cam2head[11] ; Added 31st May 2017
			count2 = count2 + 1

			ima=rotate(ima,5) ; flip Camera 2
			ima=ima[fx:fy,fx:fy]
			ima=ima-dc2*float(n_accum)/float(nacd2)
			ima=ima/(ff2[*,*,ind,jnd]*pref_fact[ind])	
			iim2[*,*,ind,jnd]=ima
			tvwin,rebin(ima,ndim/2,ndim/2)
			meanima=mean(ima)
			print,'Cam2 image contrast is=', $
			stdev(ima)/meanima 
			printf,1,'Cam2 image contrast is=', $
			stdev(ima)/meanima 
		endelse
	endfor
endfor

; this is needed simply because IMaX tricks longitudinal
; mode into vector mode

if (longitudinal eq 1) then begin
	for ind=0,n_lambda-1 do begin
		iim1[*,*,ind,0]=0.5*(iim1[*,*,ind,0]+iim1[*,*,ind,2])
		iim1[*,*,ind,1]=0.5*(iim1[*,*,ind,1]+iim1[*,*,ind,3])
		iim2[*,*,ind,0]=0.5*(iim2[*,*,ind,0]+iim2[*,*,ind,2])
		iim2[*,*,ind,1]=0.5*(iim2[*,*,ind,1]+iim2[*,*,ind,3])
	endfor
	n_pol=2
	iim1=iim1[*,*,*,0:n_pol-1]
	iim2=iim2[*,*,*,0:n_pol-1]
endif

; Starting reconstruction (need mask generated previously) 

restore,maskfile ; created with prepare_deconvolve.pro
; psf file from Jose Antonio
;restore,'./ctel.save'
restore,'Zcoef_pd_promedios_varios_desde_c2_rad_188.save' ; valid for the night of the 10th
ctel=c_medt
;ctel(*)=0.

inter=0  ;1 interactive mode
option_fringes=2  ;2  ;1   ;0

; restoring loop

for ind=0,n_lambda-1 do begin
        for jnd=0,n_pol-1 do begin
          print,'Doing restoration for lambda and pol...............',ind,jnd
          im01=iim1(*,*,ind,jnd)

		      imr1=deconvolution_IMaX7(im01,ctel,inter,option_fringes,imnr1,maskf=mask1,maskg=maskgd1) ; global psf
     		
		      iim1n(*,*,ind,jnd)=imnr1 ; not restored
          iim1(*,*,ind,jnd)=imr1   ; restored

     		  im02=iim2(*,*,ind,jnd)

		      imr2=deconvolution_IMaX7(im02,ctel,inter,option_fringes,imnr2,maskf=mask2,maskg=maskgd2) ; global psf
          iim2n(*,*,ind,jnd)=imnr2 ; not restored
          iim2(*,*,ind,jnd)=imr2   ; restored
	endfor
endfor

; demodulation Lindau

;dmod1=[[0.250271 ,    0.249951  ,   0.249992  ,   0.250204], $
;[0.437954   ,  0.416056 ,   -0.367972  ,  -0.413209], $
;[0.712257   , -0.805102   ,  0.108710 ,  -0.0517977], $
;[0.0919034  ,  -0.146047  ,  -0.611572 ,    0.681498]]

;dmod2=[[0.249960 ,   0.249999   , 0.249612  ,  0.250079], $
;[-0.416352  , -0.392155  ,  0.409879  ,  0.431250], $
;[-0.730936  ,  0.820399  , -0.118402  , 0.0580606], $
;[-0.0848565 ,   0.136714 ,   0.617866 ,  -0.691365]]

; demodulation ESRANGE 

;dmod1=[[0.250154 ,   0.250123  ,  0.250084   , 0.250032], $
;[0.453535  ,  0.369986  , -0.373724 ,  -0.358070], $
;[0.678616  , -0.801690  ,  0.120243 , -0.0378298], $
;[0.0690960 ,  -0.165518 ,  -0.582772,    0.685490]]

;dmod2=[[0.250116  ,   0.250124   ,  0.250131  ,   0.250077], $
;[-0.395565  ,  -0.317963  ,   0.383667  ,   0.352761], $
;[-0.661183  ,   0.769807  ,  -0.125787  ,  0.0339422], $
;[-0.0586242 ,    0.151446 ,    0.568101 ,   -0.667675]]

; demodulation ESRANGE corrected for 30 degrees

dmod1=[[0.250109,    0.238974,    0.260214,    0.250423], $
[0.455083   , 0.356425   ,-0.377968   ,-0.342087], $
[0.678561   ,-0.827494   ,0.0977504   ,0.0104520], $
[0.0235024   ,-0.140358   ,-0.567420   , 0.690720]]


dmod2=[[0.250067 ,   0.239074  ,  0.260436  ,  0.250237],$
[-0.396238  , -0.307138  ,  0.388434 ,   0.338140],$
[-0.660712  ,  0.794409  , -0.104223 , -0.0122207],$
[-0.0139619  ,  0.126351  ,  0.552859 ,  -0.672138 ]]
	
iid1=fltarr(ndim,ndim,n_lambda,n_pol)
iid2=fltarr(ndim,ndim,n_lambda,n_pol)
iid1n=fltarr(ndim,ndim,n_lambda,n_pol)
iid2n=fltarr(ndim,ndim,n_lambda,n_pol)

if (longitudinal eq 1) then begin

	for ind=0,n_lambda-1 do begin
		iid1[*,*,ind,0]=0.5*(iim1[*,*,ind,0]+iim1[*,*,ind,1]) ; Stokes I
		iid1[*,*,ind,1]=0.5*(iim1[*,*,ind,0]-iim1[*,*,ind,1]) ; Stokes V

		iid2[*,*,ind,0]=0.5*(iim2[*,*,ind,0]+iim2[*,*,ind,1]) ; Stokes I
		iid2[*,*,ind,1]=-0.5*(iim2[*,*,ind,0]-iim2[*,*,ind,1]) ; Stokes V

		iid1n[*,*,ind,0]=0.5*(iim1n[*,*,ind,0]+iim1n[*,*,ind,1]) ; Stokes I
		iid1n[*,*,ind,1]=0.5*(iim1n[*,*,ind,0]-iim1n[*,*,ind,1]) ; Stokes V

		iid2n[*,*,ind,0]=0.5*(iim2n[*,*,ind,0]+iim2n[*,*,ind,1]) ; Stokes I
		iid2n[*,*,ind,1]=-0.5*(iim2n[*,*,ind,0]-iim2n[*,*,ind,1]) ; Stokes V
	endfor

endif else begin

	for ind=0,n_lambda-1 do begin
	
		ima11=reform(iim1[*,*,ind,0])	
		ima12=reform(iim1[*,*,ind,1])	
		ima13=reform(iim1[*,*,ind,2])	
		ima14=reform(iim1[*,*,ind,3])	
	
		ima21=reform(iim2[*,*,ind,0])	
		ima22=reform(iim2[*,*,ind,1])	
		ima23=reform(iim2[*,*,ind,2])	
		ima24=reform(iim2[*,*,ind,3])	
	

		iid1[*,*,ind,0]=ima11*dmod1[0,0]+ima12*dmod1[1,0]+ $
		ima13*dmod1[2,0]+ima14*dmod1[3,0]
		iid1[*,*,ind,1]=ima11*dmod1[0,1]+ima12*dmod1[1,1]+ $
		ima13*dmod1[2,1]+ima14*dmod1[3,1]
		iid1[*,*,ind,2]=ima11*dmod1[0,2]+ima12*dmod1[1,2]+ $
		ima13*dmod1[2,2]+ima14*dmod1[3,2]
		iid1[*,*,ind,3]=ima11*dmod1[0,3]+ima12*dmod1[1,3]+ $
		ima13*dmod1[2,3]+ima14*dmod1[3,3]
	
		iid2[*,*,ind,0]=ima21*dmod2[0,0]+ima22*dmod2[1,0]+ $
		ima23*dmod2[2,0]+ima24*dmod2[3,0]
		iid2[*,*,ind,1]=ima21*dmod2[0,1]+ima22*dmod2[1,1]+ $
		ima23*dmod2[2,1]+ima24*dmod2[3,1]
		iid2[*,*,ind,2]=ima21*dmod2[0,2]+ima22*dmod2[1,2]+ $
		ima23*dmod2[2,2]+ima24*dmod2[3,2]
		iid2[*,*,ind,3]=ima21*dmod2[0,3]+ima22*dmod2[1,3]+ $
		ima23*dmod2[2,3]+ima24*dmod2[3,3]
		
		ima11=reform(iim1n[*,*,ind,0])	
		ima12=reform(iim1n[*,*,ind,1])	
		ima13=reform(iim1n[*,*,ind,2])	
		ima14=reform(iim1n[*,*,ind,3])	
	
		ima21=reform(iim2n[*,*,ind,0])	
		ima22=reform(iim2n[*,*,ind,1])	
		ima23=reform(iim2n[*,*,ind,2])	
		ima24=reform(iim2n[*,*,ind,3])	
	
		iid1n[*,*,ind,0]=ima11*dmod1[0,0]+ima12*dmod1[1,0]+ $
		ima13*dmod1[2,0]+ima14*dmod1[3,0]
		iid1n[*,*,ind,1]=ima11*dmod1[0,1]+ima12*dmod1[1,1]+ $
		ima13*dmod1[2,1]+ima14*dmod1[3,1]
		iid1n[*,*,ind,2]=ima11*dmod1[0,2]+ima12*dmod1[1,2]+ $
		ima13*dmod1[2,2]+ima14*dmod1[3,2]
		iid1n[*,*,ind,3]=ima11*dmod1[0,3]+ima12*dmod1[1,3]+ $
		ima13*dmod1[2,3]+ima14*dmod1[3,3]
	
		iid2n[*,*,ind,0]=ima21*dmod2[0,0]+ima22*dmod2[1,0]+ $
		ima23*dmod2[2,0]+ima24*dmod2[3,0]
		iid2n[*,*,ind,1]=ima21*dmod2[0,1]+ima22*dmod2[1,1]+ $
		ima23*dmod2[2,1]+ima24*dmod2[3,1]
		iid2n[*,*,ind,2]=ima21*dmod2[0,2]+ima22*dmod2[1,2]+ $
		ima23*dmod2[2,2]+ima24*dmod2[3,2]
		iid2n[*,*,ind,3]=ima21*dmod2[0,3]+ima22*dmod2[1,3]+ $
		ima23*dmod2[2,3]+ima24*dmod2[3,3]
	
	endfor
endelse

; Stokes I ad-hoc xtalk removal. Removing only xtalk derived in continuum (n_lambda-1)

itos1=fltarr(n_lambda,n_pol)
itos2=fltarr(n_lambda,n_pol)

for ind=0,n_lambda-1 do begin
	for jnd=1,n_pol-1 do begin
		itos1[ind,jnd]=mean(iid1[*,*,ind,jnd])/ $
		mean(iid1[*,*,ind,0])
		itos2[ind,jnd]=mean(iid2[*,*,ind,jnd])/ $
		mean(iid2[*,*,ind,0])
	endfor
endfor

for ind=0,n_lambda-1 do begin
	for jnd=1,n_pol-1 do begin
		iid1[*,*,ind,jnd]=iid1[*,*,ind,jnd]-itos1[n_lambda-1,jnd]*iid1[*,*,ind,0]	
		iid2[*,*,ind,jnd]=iid2[*,*,ind,jnd]-itos2[n_lambda-1,jnd]*iid2[*,*,ind,0]	
	endfor
endfor

print,'Only used continuum xtalk (last column)'
print,'Continuum xtalk parameters for Cam1=',itos1
print,'Continuum xtalk parameters for Cam2=',itos2
printf,1,'Only used continuum xtalk (last column)'
printf,1,'Continuum xtalk parameters for Cam1=',itos1
printf,1,'Continuum xtalk parameters for Cam2=',itos2

for ind=0,n_lambda-1 do begin
	for jnd=1,n_pol-1 do begin
		itos1[ind,jnd]=mean(iid1n[*,*,ind,jnd])/ $
		mean(iid1n[*,*,ind,0])
		itos2[ind,jnd]=mean(iid2n[*,*,ind,jnd])/ $
		mean(iid2n[*,*,ind,0])
	endfor
endfor

for ind=0,n_lambda-1 do begin
	for jnd=1,n_pol-1 do begin
		iid1n[*,*,ind,jnd]=iid1n[*,*,ind,jnd]-itos1[n_lambda-1,jnd]*iid1n[*,*,ind,0]	
		iid2n[*,*,ind,jnd]=iid2n[*,*,ind,jnd]-itos2[n_lambda-1,jnd]*iid2n[*,*,ind,0]	
	endfor
endfor

print,'Only used continuum xtalk (last column)'
print,'Continuum xtalk parameters for Cam1 (not restored)=',itos1
print,'Continuum xtalk parameters for Cam2 (not restored)=',itos2
printf,1,'Only used continuum xtalk (last column)'
printf,1,'Continuum xtalk parameters for Cam1 (not restored)=',itos1
printf,1,'Continuum xtalk parameters for Cam2 (not restored)=',itos2

; equalize intensities Cam1 and Cam2

ratcam=fltarr(n_lambda)

for ind=0,n_lambda-1 do begin
	ratcam[ind]=mean(iid1[*,*,ind,0])/mean(iid2[*,*,ind,0])
	iid2[*,*,ind,*]=iid2[*,*,ind,*]*ratcam[ind]
endfor

print,'Camera 1/Camera 2 ratios=',ratcam
printf,1,'Camera 1/Camera 2 ratios=',ratcam

for ind=0,n_lambda-1 do begin
	ratcam[ind]=mean(iid1n[*,*,ind,0])/mean(iid2n[*,*,ind,0])
	iid2n[*,*,ind,*]=iid2n[*,*,ind,*]*ratcam[ind]
endfor

print,'Camera 1/Camera 2 ratios (not restored)=',ratcam
printf,1,'Camera 1/Camera 2 ratios (not restored)=',ratcam


; shift Cam2 to Cam1

; default values
shiftpl_def=[5.15,-8.05]

shiftpl=fltarr(n_lambda,2)

for ind=0,n_lambda-1 do begin
	ii1x=iid1[ss1h-ss1q:ss1h+ss1q,ss2h-ss2q:ss2h+ss2q,ind,0]
	ii2x=iid2[ss1h-ss1q:ss1h+ss1q,ss2h-ss2q:ss2h+ss2q,ind,0]
	tmpal = align(ii1x,ii2x)
	checkshift=where(finite(tmpal) eq 0)
	if (checkshift ne -1) then tmpal(checkshift)=shiftpl_def(checkshift)
	shiftpl[ind,*]=tmpal 
	for jnd=0,n_pol-1 do begin
		iid2[*,*,ind,jnd]=fft_shift(iid2[*,*,ind,jnd],shiftpl[ind,*])
	endfor
endfor


print,'Shift between Camera 2 and 1 is=',shiftpl
printf,1,'Shift between Camera 2 and 1 is=',shiftpl

shiftpln=fltarr(n_lambda,2)

for ind=0,n_lambda-1 do begin
	ii1x=iid1n[ss1h-ss1q:ss1h+ss1q,ss2h-ss2q:ss2h+ss2q,ind,0]
	ii2x=iid2n[ss1h-ss1q:ss1h+ss1q,ss2h-ss2q:ss2h+ss2q,ind,0]
	tmpal=align(ii1x,ii2x)
	checkshift=where(finite(tmpal) eq 0)
	if (checkshift ne -1) then tmpal(checkshift)=shiftpl_def(checkshift)
	shiftpln[ind,*]=tmpal 
	for jnd=0,n_pol-1 do begin
		iid2n[*,*,ind,jnd]=fft_shift(iid2n[*,*,ind,jnd],shiftpln[ind,*])
	endfor
endfor
print,'Shift between Camera 2 and 1 is (not restored)=',shiftpln
printf,1,'Shift between Camera 2 and 1 is (not restored)=',shiftpln

; add Cam1 and Cam2 to minimize stdevs.

iid=fltarr(ndim,ndim,n_lambda,n_pol) ; data array
ssadd=fltarr(n_lambda,n_pol,50)

mm=.35+0.3*findgen(50)/49.

print,'    Lambda ---  Polarization --- Cam1/Cam2  weight'
printf,1,'    Lambda ---  Polarization --- Cam1/Cam2  weight'
for ind=0,n_lambda-1 do begin
	for jnd=1,n_pol-1 do begin
		for knd=0,49 do begin
			iid[*,*,ind,jnd]=mm[knd]*iid1[*,*,ind,jnd]+(1.-mm[knd])*iid2[*,*,ind,jnd]
			ssadd[ind,jnd,knd]=stdev(iid[*,*,ind,jnd])
		endfor
;		plot,mm[*],ssadd[ind,jnd,*],/yno
		co=poly_fit(mm[*],ssadd[ind,jnd,*],2)
		aa=-co[1]/(2.*co[2])
;		ver,aa
		aa=mm(where(ssadd[ind,jnd,*] eq min(ssadd[ind,jnd,*])))
		print,ind,'       ',jnd,'       ',aa
		printf,1,ind,'       ',jnd,'       ',aa
		iid[*,*,ind,jnd]=aa[0]*iid1[*,*,ind,jnd]+(1.-aa[0])*iid2[*,*,ind,jnd]
		iim1[*,*,ind,jnd]=aa[0]*iid1[*,*,ind,jnd]-(1.-aa[0])*iid2[*,*,ind,jnd]
		;this array is motion induced xtalk. We reuse an old array to save memory
	endfor
	iid[*,*,ind,0]=0.5*iid1[*,*,ind,0]+0.5*iid2[*,*,ind,0]
	iim1[*,*,ind,0]=0.5*iid1[*,*,ind,0]-0.5*iid2[*,*,ind,0]

endfor

iidn=fltarr(ndim,ndim,n_lambda,n_pol) ; data array

print,'    Lambda ---  Polarization --- Cam1/Cam2  weight (not restored)'
printf,1,'    Lambda ---  Polarization --- Cam1/Cam2  weight (not restored)'
for ind=0,n_lambda-1 do begin
	for jnd=1,n_pol-1 do begin
		for knd=0,49 do begin
			iidn[*,*,ind,jnd]=mm[knd]*iid1n[*,*,ind,jnd]+(1.-mm[knd])*iid2n[*,*,ind,jnd]
			ssadd[ind,jnd,knd]=stdev(iidn[*,*,ind,jnd])
		endfor
;		plot,mm[*],ssadd[ind,jnd,*],/yno
		co=poly_fit(mm[*],ssadd[ind,jnd,*],2)
		aa=-co[1]/(2.*co[2])
;		ver,aa
		aa=mm(where(ssadd[ind,jnd,*] eq min(ssadd[ind,jnd,*])))
		print,ind,'       ',jnd,'       ',aa
		printf,1,ind,'       ',jnd,'       ',aa
		iidn[*,*,ind,jnd]=aa[0]*iid1n[*,*,ind,jnd]+(1.-aa[0])*iid2n[*,*,ind,jnd]
		iim1n[*,*,ind,jnd]=aa[0]*iid1n[*,*,ind,jnd]-(1.-aa[0])*iid2n[*,*,ind,jnd]
		;this array is motion induced xtalk. We reuse an old array to save memory
	endfor
	iidn[*,*,ind,0]=0.5*iid1n[*,*,ind,0]+0.5*iid2n[*,*,ind,0]
	iim1n[*,*,ind,0]=0.5*iid1n[*,*,ind,0]-0.5*iid2n[*,*,ind,0]

endfor

; remove any residual nonzero trend

smwin=5
ranfit=1

for ind=0,n_lambda-1 do begin
	meansti=mean(iid[*,*,ind,0])
	for jnd=1,n_pol-1 do begin
		tmp=smooth(iid[*,*,ind,jnd],smwin)
		sss=sfit(tmp,ranfit)
		print,'Residual non-zero trends=',ind,jnd,mean(abs(sss))/meansti
		printf,1,'Residual non-zero trends=',ind,jnd,mean(abs(sss))/meansti
		iid[*,*,ind,jnd]=iid[*,*,ind,jnd]-sss*iid[*,*,ind,0]/meansti 
		
		tmp=smooth(iim1[*,*,ind,jnd],smwin)
		sss=sfit(tmp,ranfit)
		iim1[*,*,ind,jnd]=iim1[*,*,ind,jnd]-sss 
	endfor
endfor

for ind=0,n_lambda-1 do begin
	meansti=mean(iidn[*,*,ind,0])
	for jnd=1,n_pol-1 do begin
		tmp=smooth(iidn[*,*,ind,jnd],smwin)
		sss=sfit(tmp,ranfit)
		print,'Residual non-zero trends (not restored)=',ind,jnd,mean(abs(sss))/meansti
		printf,1,'Residual non-zero trends (not restored)=',ind,jnd,mean(abs(sss))/meansti
		iidn[*,*,ind,jnd]=iidn[*,*,ind,jnd]-sss*iidn[*,*,ind,0]/meansti 
		
		tmp=smooth(iim1n[*,*,ind,jnd],smwin)
		sss=sfit(tmp,ranfit)
		iim1n[*,*,ind,jnd]=iim1n[*,*,ind,jnd]-sss 
	endfor
endfor

; mask uncommon area

fixshift=[fix(total(shiftpl[*,0])/float(n_lambda))+1*sign(fix(total(shiftpl[*,0])/float(n_lambda))), $
fix(total(shiftpl[*,1])/float(n_lambda))+1*sign(fix(total(shiftpl[*,1])/float(n_lambda)))]

if fixshift[0] gt 0 then begin
	mx1=0 & mx2=fixshift[0]
endif else begin
	mx1=ndim+fixshift[0] & mx2=ndim-1
endelse

if fixshift[1] gt 0 then begin
	my1=0 & my2=fixshift[1]
endif else begin
	my1=ndim+fixshift[1] & my2=ndim-1
endelse

for ind=0,n_lambda-1 do begin
	for jnd=0,n_pol-1 do begin
		if (jnd eq 0) then begin
			iid[mx1:mx2,*,ind,jnd]=mean(iid[*,*,ind,jnd])
			iid[*,my1:my2,ind,jnd]=mean(iid[*,*,ind,jnd])
		endif else begin
			iid[mx1:mx2,*,ind,jnd]=0.
			iid[*,my1:my2,ind,jnd]=0.
		endelse
	endfor 
endfor

for ind=0,n_lambda-1 do begin
	for jnd=0,n_pol-1 do begin
		if (jnd eq 0) then begin
			iidn[mx1:mx2,*,ind,jnd]=mean(iidn[*,*,ind,jnd])
			iidn[*,my1:my2,ind,jnd]=mean(iidn[*,*,ind,jnd])
		endif else begin
			iidn[mx1:mx2,*,ind,jnd]=0.
			iidn[*,my1:my2,ind,jnd]=0.
		endelse
	endfor 
endfor

; remove continuum from Q, U and V signals for lambdas 2 and 3 (specific of V5-6)

if (n_pol eq 4 and n_lambda eq 5) then begin

	for jnd=1,n_pol-1 do begin
		conpol=reform(iid[*,*,n_lambda-1,jnd])
		stpol=stdev(conpol)
		dd=where(abs(conpol) gt 3.*stpol,nns)
		if (nns gt 0) then conpol[dd]=0.
		conpol=smooth(conpol,6)
		for ind=2,3 do begin
			co=poly_fit(conpol[*,*],iid[*,*,ind,jnd],1)  
			iid[*,*,ind,jnd]=iid[*,*,ind,jnd]-conpol*co(1)
			print,'Substracted continuum polarization=',ind,jnd,co(1)
			printf,1,'Substracted continuum polarization=',ind,jnd,co(1)
		endfor
	endfor

	for jnd=1,n_pol-1 do begin
		conpol=reform(iidn[*,*,n_lambda-1,jnd])
		stpol=stdev(conpol)
		dd=where(abs(conpol) gt 3.*stpol,nns)
		if (nns gt 0) then conpol[dd]=0.
		conpol=smooth(conpol,6)
		for ind=2,3 do begin
			co=poly_fit(conpol[*,*],iidn[*,*,ind,jnd],1)  
			iidn[*,*,ind,jnd]=iidn[*,*,ind,jnd]-conpol*co(1)
			print,'Substracted continuum polarization (not restored)=',ind,jnd,co(1)
			printf,1,'Substracted continuum polarization (not restored)=',ind,jnd,co(1)
		endfor
	endfor
endif
; The above does not work for lambdas 0 and 1. Stokes I is
; too diferent from continuum

; remove ad-hoc xtalk between Q,U and V
; xtalk is estimated in non restored data (specific of V modes)

if (n_pol eq 4) then begin

	vtoqu=fltarr(2,n_lambda) 

	for ind=0,n_lambda-1 do begin
		co=poly_fit(iidn[30:ndim-30,30:ndim-30,ind,3], $
		iidn[30:ndim-30,30:ndim-30,ind,1],1)
		vtoqu[0,ind]=co[1]
		co=poly_fit(iidn[30:ndim-30,30:ndim-30,ind,3], $
		iidn[30:ndim-30,30:ndim-30,ind,2],1)
		vtoqu[1,ind]=co[1]
	endfor

	vtoq=mean(vtoqu[0,0:n_lambda-2]) ; exclude continuum
	vtou=mean(vtoqu[1,0:n_lambda-2]) ; exclude continuum

; we remove at each wavelength different residual xtalks, not the
; means computed in vtoq and vtou. This has been decided because
; these values (residual xtalk at different lambdas) are consistenly
; similar. Not in continuum ! (added in 21-10-10)

; continuum has not much info and values there are not as good. Use mean of
; last two wavelength points. Added in 21-10-10

	if (n_lambda le 5) then for ind=0,1 do vtoqu[ind,n_lambda-1]=mean(vtoqu[ind,n_lambda-3:n_lambda-2])
	if (n_lambda eq 12) then begin
		for ind=0,1 do begin
			vtoqu[ind,0:2]=mean(vtoqu[ind,3:5])
			vtoqu[ind,n_lambda-3:n_lambda-1]=mean(vtoqu[ind,6:8])
		endfor
	endif

	for ind=0,n_lambda-1 do begin
;		iid[*,*,ind,1]=iid[*,*,ind,1]-vtoq*iid[*,*,ind,3]
;		iid[*,*,ind,2]=iid[*,*,ind,2]-vtou*iid[*,*,ind,3]
;		iidn[*,*,ind,1]=iidn[*,*,ind,1]-vtoq*iidn[*,*,ind,3]
;		iidn[*,*,ind,2]=iidn[*,*,ind,2]-vtou*iidn[*,*,ind,3]
		iid[*,*,ind,1]=iid[*,*,ind,1]-vtoqu[0,ind]*iid[*,*,ind,3]
		iid[*,*,ind,2]=iid[*,*,ind,2]-vtoqu[1,ind]*iid[*,*,ind,3]
		iidn[*,*,ind,1]=iidn[*,*,ind,1]-vtoqu[0,ind]*iidn[*,*,ind,3]
		iidn[*,*,ind,2]=iidn[*,*,ind,2]-vtoqu[1,ind]*iidn[*,*,ind,3]
	endfor

	print,'ad-hoc xtalk from V to Q is (not restored)=',transpose(vtoqu[0,*])
	printf,1,'ad-hoc xtalk from V to Q is (not restored)=',transpose(vtoqu[0,*])
	print,'ad-hoc xtalk from V to U is (not restored)=',transpose(vtoqu[1,*])
	printf,1,'ad-hoc xtalk from V to U is (not restored)=',transpose(vtoqu[1,*])

	print,'ad-hoc xtalk from V to Q used is (to both)=',vtoq
	printf,1,'ad-hoc xtalk from V to Q used is (to both)=',vtoq
	print,'ad-hoc xtalk from V to U used is (to both)=',vtou
	printf,1,'ad-hoc xtalk from V to U used is (to both)=',vtou

endif

close,1

; testing further I and motion xtalk removal

cond_z=ndim mod (findgen(20)+1) 
dd=reverse(where(cond_z eq 0))+1
subfr=ndim/dd(0)

for ind=0,n_lambda-1 do begin
	iii=iid[*,*,ind,0]
	for jnd=1,n_pol-1 do begin
		pol=iid[*,*,ind,jnd]
		xta=iim1[*,*,ind,jnd]
		iid[*,*,ind,jnd]=test_local(pol,xta,iii,subfr)
	endfor
endfor

for ind=0,n_lambda-1 do begin
	iii=iidn[*,*,ind,0]
	for jnd=1,n_pol-1 do begin
		pol=iidn[*,*,ind,jnd]
		xta=iim1n[*,*,ind,jnd]
		iidn[*,*,ind,jnd]=test_local(pol,xta,iii,subfr)
	endfor
endfor

;path_out='/media/FREECOM HDD/longitudinal/'
save,iid,iidn, date_arraycam1, date_arraycam2,$
filename=path_out+'reduc_rnr_'+strcompress(string(fix(obsen)),/remove_all)+ $
'_'+strcompress(string(fix(cyclen)),/remove_all)+'.sav'
;save,iim1,iim1n, $
;filename=path_out+'reduc_rnr_'+strcompress(string(fix(obsen)),/remove_all)+ $
;'_'+strcompress(string(fix(cyclen)),/remove_all)+'_xtalk.save'
print,'Saving=', $
path_out+'reduc_rnr_'+strcompress(string(fix(obsen)),/remove_all)+ $
'_'+strcompress(string(fix(cyclen)),/remove_all)+'.sav'

nla=n_lambda<11
vvv=(total(iid[*,*,0:(nla-1)/2,1],3)/float((nla-1)/2)- $
total(iid[*,*,(nla-1)/2+1:nla-1,1],3)/float(nla-((nla-1)/2+1)))/ $
mean(iid[*,*,n_lambda-1,0])

window,0,xsize=900,ysize=800                                                            
tvframe,vvv<.005>(-.005),/bar,/aspect,xrange=[0,ndim*.055],yrange=[0,ndim*.055], $
title='IMaX/Sunrise Stokes V/Ic',charsize=2.
write_png, path_out + 'magneto_scale.png', tvrd(/true)

iii=iid[*,*,n_lambda-1,0]
iii=iii/mean(iii)
window,0,xsize=900,ysize=800                                                            
tvframe,iii>0.8<1.35,/bar,/aspect,xrange=[0,ndim*.055],yrange=[0,ndim*.055], $
title='IMaX/Sunrise Continuum',charsize=1.5
write_png, path_out + 'continuum_scale.png', tvrd(/true)

end
