pro image_viewer_str, path, outfilename

; This program analyzes an IMaX/SUNRISE directory 
; and generates a save file that can be used to 
; study the contents of that day (data mining). 

; Ends printing info on where flats, darks, pinholes
; and PD sets are available

; First run it with no flat and dark, then it can be run 
; with flats and darks included. 

;IDL>image_viewer_str,'/mnt/Iomega HDD/2009y_06m_09d/','set09_hdr.save'
;IDL>image_viewer_str,'/media/FREECOM HDD/2009y_06m_09d/','set09_hdr.save'

; NOTE: needs changing the name of indexing file manually

; path is directory for data
; outfilename is output save set

; Two save files come out of this routine:
;   index: - name given in procedure
;          - must be runfor each day seperately
;          - imas1/2: path to file
;          - nima1/2: number of camera1/2 files
;
;   set: - name given as outfile input
;        - basically a header 
;        - all files with same numobse1/2 will be pointing at the same
;          place (except when it breaks)
;        - numcycl1/2 gives cycle number (cycles over stokes params?)

;---------------------------------------------------------------------------;
; do indexing only first time is run. Output is indexing file

allFiles1 = file_search(path, '*Camera_1*.fits.gz', count=oriCount1)
allFiles2 = file_search(path, '*Camera_2*.fits.gz', count=oriCount2)
print, path+'*Camera_1 *.fits.gz'
save, allFiles1, oriCount1, allFiles2, oriCount2, filename = 'tempy/ind_11thJune_wdtest.sav' 
;---------------------------------------------------------------------------;
;Indexing takes some time. Once this saved is generated
;restore it instead of do the file_search again

;restore, 'imax_index_12thJune.save'
;---------------------------------------------------------------------------;

;print, 'Camera 1 images for this day is =', nima1
;print, 'Camera 2 images for this day is =', nima2

ind_first1 = 0 
ind_last1  = oriCount1
ind_first2 = 0  
ind_last2  = oriCount2  

nima1 = ind_last1 - ind_first1
nima2 = ind_last2 - ind_first2
print, oriCOunt1, ind_last1, ind_first1, nima1
; if no darks are available 
dc1 = fltarr(1024, 1024)
dc2 = fltarr(1024, 1024)

; if darks are available use this
;pathd='/home/vmp/IDLWorkspace/kiruna_flight/data/firstset/'
;restore,pathd+'0906509_dc_set1.save'

; restore flats (if available)

;restore,'flats_DATA_09_set1.save'

; initialize arrays (a description is available in imax_hdr_str.pro)

int1       = fltarr(nima1)           ; mean intensity central part
sss1       = fltarr(nima1)           ; standar deviation central part
time1      = fltarr(nima1)           ; time of given image
timerun1   = fltarr(nima1)           ; time at which observing sequence started
imgacst1   = fltarr(nima1)           ; time at which observing cycle started
imgacen1   = fltarr(nima1)           ; time at which observing cycle ended 
xpoint1    = fltarr(nima1)           ; x coordinate solar pointing
ypoint1    = fltarr(nima1)           ; y coordinate solar pointing
f2m1       = fltarr(nima1)           ; F2 Mechanism (identifies darks and pinholes)
cwl1       = fltarr(nima1)           ; CWS status (lock is good)
pss1       = fltarr(nima1)           ; pointing system (identifies flats)
pdp1       = fltarr(nima1)           ; phase diversity  
nwavel1    = fltarr(nima1)           ; number of wavelengths
nacc1      = fltarr(nima1)           ; numer of accumulations 
npol1      = fltarr(nima1)           ; number of polarizations
dwavel1    = fltarr(nima1)           ; positions of wavelengths from line center
polst1     = fltarr(nima1)           ; polarization state of current frame
lamst1     = fltarr(nima1)           ; wavelength state of current frame
ettemp1    = fltarr(nima1)           ; etalon temperature
lcvrtemp1  = fltarr(nima1)           ; LCVRs temperature
doubltemp1 = fltarr(nima1)           ; Doublets temperature
acoef1     = fltarr(nima1)           ; wavelength calibration coefficient
obsmod1    = fltarr(nima1)           ; observing mode
sscpar1    = fltarr(nima1)           ; SSC sequential counter
numcycl1   = fltarr(nima1)           ; cycle number
numobse1   = fltarr(nima1)           ; observing sequence number
rms1       = fltarr(nima1) 
; next is the same but for Cam 2
 
int2       = fltarr(nima2)
sss2       = fltarr(nima2)
time2      = fltarr(nima2)
timerun2   = fltarr(nima2)
imgacst2   = fltarr(nima2)
imgacen2   = fltarr(nima2)
xpoint2    = fltarr(nima2)
ypoint2    = fltarr(nima2)
f2m2       = fltarr(nima2)
cwl2       = fltarr(nima2)
pss2       = fltarr(nima2)
pdp2       = fltarr(nima2)
nwavel2    = fltarr(nima2)
nacc2      = fltarr(nima2)
npol2      = fltarr(nima2)
dwavel2    = fltarr(nima2)
polst2     = fltarr(nima2)
lamst2     = fltarr(nima2)
ettemp2    = fltarr(nima2)
lcvrtemp2  = fltarr(nima2)
doubltemp2 = fltarr(nima2)
acoef2     = fltarr(nima2)
obsmod2    = fltarr(nima2)
sscpar2    = fltarr(nima2)
numcycl2   = fltarr(nima2)
numobse2   = fltarr(nima2)
rms2       = fltarr(nima2)

for ind = ind_first1, ind_last1 - 1 do begin
	iind = ind - ind_first1
        print, '--------------------------------------'
	print, '  Image  = ', iind
        print, allFiles1[iind]
	ima = float(readfits(allFiles1[iind], hdr))
        imax_hdr_str, hdr, hdrstr
	ss1   = hdrstr.sz1 ; better be 1024
	ss2   = hdrstr.sz2
	ima   = ima - dc1
	nlamb = hdrstr.n_wave   - 1
	npola = hdrstr.pol_stat - 1

;	this is used only if flats are available
;	ima = ima/ff1[*,*,nlamb,npola]
        if iind EQ 0 then begin
           imas1 = allFiles1[iind]
        endif else begin
           imas1 = [[imas1], [allFiles1[iind]]]
        endelse 
	int1[iind] = mean(ima[ss1/2-ss1/4 : ss1/2+ss1/4, ss2/2-ss2/4 : ss2/2+ss2/4])	
	sss1[iind] = stdev(ima[ss1/2-ss1/4 : ss1/2+ss1/4, ss2/2-ss2/4 : ss2/2+ss2/4])
	rms1[iind] = sss1[iind] / int1[iind]	
	print, 'Contrast Cam1 is = ', rms1[iind]
        print, '--------------------------------------'
	time1[iind]      = hdrstr.time
	timerun1[iind]   = hdrstr.time_run
	imgacst1[iind]   = hdrstr.img_acst
	imgacen1[iind]   = hdrstr.img_acen
	xpoint1[iind]    = hdrstr.xcen
	ypoint1[iind]    = hdrstr.ycen
	f2m1[iind]       = hdrstr.f2_mech
	cwl1[iind]       = hdrstr.cw_loop
	pss1[iind]       = hdrstr.ps_state
	pdp1[iind]       = hdrstr.pd_plate
	nwavel1[iind]    = hdrstr.num_wave
	nacc1[iind]      = hdrstr.num_accu
	npol1[iind]      = hdrstr.num_pola
	dwavel1[iind]    = hdrstr.d_wave
	polst1[iind]     = hdrstr.pol_stat
	lamst1[iind]     = hdrstr.lam_stat
	ettemp1[iind]    = hdrstr.etal_tem
	lcvrtemp1[iind]  = hdrstr.rocl_tem
	doubltemp1[iind] = hdrstr.etad_tem
	acoef1[iind]     = hdrstr.etal_aco
	obsmod1[iind]    = hdrstr.obs_mode
	sscpar1[iind]    = hdrstr.ssc
	numobse1[iind]   = hdrstr.num_obse
	numcycl1[iind]   = hdrstr.num_cycl
;	tvwin, rebin(ima, 512, 512)

  
endfor


for ind = ind_first2, ind_last2 - 4 do begin
	iind = ind - ind_first2
        print, '--------------------------------------'
	print, '  Image  = ', iind
        print, allFiles1[iind]
	ima = float(readfits(allFiles2[iind], hdr))
	imax_hdr_str, hdr, hdrstr
	ima   = ima-dc2
	nlamb = hdrstr.n_wave   - 1
	npola = hdrstr.pol_stat - 1

;	this is used only if flats are available
;	ima = ima/ff2[*, *, nlamb, npola]
        
        if iind EQ 0 then begin
           imas2 = allFiles2[iind]
        endif else begin
           imas2 = [[imas2], [allFiles2[iind]]]
        endelse 

	int2[iind] = mean(ima[ss1/2-ss1/4 : ss1/2+ss1/4, ss2/2-ss2/4 : ss2/2+ss2/4])	
	sss2[iind] = stdev(ima[ss1/2-ss1/4 : ss1/2+ss1/4, ss2/2-ss2/4 : ss2/2+ss2/4])	
        rms2[iind] = sss2[iind] / int2[iind]
	print, 'Contrast Cam2 is = ', rms2[iind]
        print, '--------------------------------------'
	time2[iind]      = hdrstr.time
	timerun2[iind]   = hdrstr.time_run
	imgacst2[iind]   = hdrstr.img_acst
	imgacen2[iind]   = hdrstr.img_acen
	xpoint2[iind]    = hdrstr.xcen
	ypoint2[iind]    = hdrstr.ycen
	f2m2[iind]       = hdrstr.f2_mech
	cwl2[iind]       = hdrstr.cw_loop
	pss2[iind]       = hdrstr.ps_state
	pdp2[iind]       = hdrstr.pd_plate
	nwavel2[iind]    = hdrstr.num_wave
	nacc2[iind]      = hdrstr.num_accu
	npol2[iind]      = hdrstr.num_pola
	dwavel2[iind]    = hdrstr.d_wave
	polst2[iind]     = hdrstr.pol_stat
	lamst2[iind]     = hdrstr.lam_stat
	ettemp2[iind]    = hdrstr.etal_tem
	lcvrtemp2[iind]  = hdrstr.rocl_tem
	doubltemp2[iind] = hdrstr.etad_tem
	acoef2[iind]     = hdrstr.etal_aco
	obsmod2[iind]    = hdrstr.obs_mode
	sscpar2[iind]    = hdrstr.ssc
	numobse2[iind]   = hdrstr.num_obse
	numcycl2[iind]   = hdrstr.num_cycl
	
;	tvwin, rebin(ima, 512, 512)
endfor


save, ss1, ss2, $
int1, sss1, time1, timerun1, imgacst1, imgacen1, xpoint1, ypoint1, f2m1, cwl1, pss1, pdp1, $
nwavel1, nacc1, npol1, dwavel1, polst1, lamst1, ettemp1, lcvrtemp1, doubltemp1, acoef1,  $
obsmod1, sscpar1, numobse1, numcycl1, rms1, imas1,  $
int2, sss2, time2, timerun2, imgacst2, imgacen2, xpoint2, ypoint2, f2m2, cwl2, pss2, pdp2, $
nwavel2, nacc2, npol2, dwavel2, polst2, lamst2, ettemp2, lcvrtemp2, doubltemp2, acoef2,  $
obsmod2, sscpar2, numobse2, numcycl2, rms2, imas2, filename = outfilename

; do some data mining
; find obsenum for flats darks pinholes and PD
; writes outputfilename

openw, 1, strmid(outfilename, 0, strpos(outfilename, '.s')) + '.txt'

printf, 1, path
printf, 1, ' '

; flats

ddf = where(pss1 eq 3, nnf)
if (nnf gt 0) then begin
	obsef = numobse1[ddf[0]]
	for ind = 1, nnf - 1 do begin
		if (numobse1[ddf[ind]] ne  obsef[n_elements(obsef[*]) - 1]) then begin
			obsef = [obsef, numobse1[ddf[ind]]]
		endif
	endfor
	print, 'Flats in series (observations and PD)  = ', obsef
	printf, 1, 'Flats in series (observations and PD)  = ', obsef
endif else begin
	print, 'No flats in series'
	printf, 1, 'No flats in series'
endelse

;darks

ddd = where(f2m1 eq 3, nnd)
if (nnd gt 0) then begin
	obsed = numobse1[ddd[0]]
	for ind = 1, nnd - 1 do begin
		if (numobse1[ddd[ind]] ne  obsed[n_elements(obsed[*]) - 1]) then begin
         		 obsed = [obsed, numobse1[ddd[ind]]]
                endif
	endfor
	print, 'Darks in series (observations and PD)  = ', obsed
	printf, 1, 'Darks in series (observations and PD)  = ', obsed
endif else begin
	print, 'No darks in series'
	printf, 1, 'No darks in series'
endelse

;pinholes
ddp = where(f2m1 eq 2, nnp)
if (nnp gt 0) then begin
	obsep = numobse1[ddp[0]]
	for ind = 1, nnp - 1 do begin
		if (numobse1[ddp[ind]] ne  obsep[n_elements(obsep[*]) - 1]) then begin 
            		obsep = [obsep, numobse1[ddp[ind]]]
                endif
	endfor	
        print, 'Pinholes in series (observations and PD)  = ', obsep
	printf, 1, 'Pinholes in series (observations and PD)  = ', obsep
endif else begin
	printf, 1, 'No pinholes in series'
endelse

; PD sets
ddpd = where(pdp1 eq 1, nnpd)
if (nnpd gt 0) then begin
	obsepd = numobse1[ddpd[0]]
	for ind = 1, nnpd - 1 do begin
		if (numobse1[ddpd[ind]] ne obsepd[n_elements(obsepd[*]) - 1]) then begin ; changed this from sizeof to n_el... /01/17
                        obsepd = [obsepd, numobse1[ddpd[ind]]]
                endif
	endfor
	print, 'PD in series  = ', obsepd
	printf, 1, 'PD in series  = ', obsepd
endif else begin
	print, 'No PD in series'
	printf, 1, 'No PD in series'
endelse

close, 1	
if n_elements(obsef)  EQ 0 then obsef  = 'No flats'
if n_elements(obsed)  EQ 0 then obsed  = 'No darks'
if n_elements(obsep)  EQ 0 then obsep  = 'No pinholes'
if n_elements(obsepd) EQ 0 then obsepd = 'No PDs'
save, ss1, ss2, $
int1, sss1, time1, timerun1, imgacst1, imgacen1, xpoint1, ypoint1, f2m1, cwl1, pss1, pdp1, $
nwavel1, nacc1, npol1, dwavel1, polst1, lamst1, ettemp1, lcvrtemp1, doubltemp1, acoef1,  $
obsmod1, sscpar1, numobse1, numcycl1, rms1, imas1,  $
int2, sss2, time2, timerun2, imgacst2, imgacen2, xpoint2, ypoint2, f2m2, cwl2, pss2, pdp2, $
nwavel2, nacc2, npol2, dwavel2, polst2, lamst2, ettemp2, lcvrtemp2, doubltemp2, acoef2,  $
obsmod2, sscpar2, numobse2, numcycl2, rms2, imas2, obsef, obsed, obsep, obsepd, filename = outfilename
end
