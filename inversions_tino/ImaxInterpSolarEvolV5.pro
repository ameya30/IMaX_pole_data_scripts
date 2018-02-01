;*****************************************************************************
; ~/MyProjects/src/idlstuff/idl_jmb/ImaxInterpSolarEvolV5
;*****************************************************************************
;* Die Stokes-Inversion geht davon aus, dass alle Wellenlängenpositionen und
;* alle Polarisationszustaende gleichzeitig aufgenommen wurden, was aber bei
;* IMaX nicht der Fall ist. Diese Routine nimmt 3 IMaX-Aufnahmen einer
;* Zeitserie im V8-4 Mode und führt eine zeitliche Interpolation für den
;* mittleren Zeitpunkt durch, für den dann die Inversion gerechnet wird.
;* Die Input-Daten von IMaX liegen als 4D-Fits-Dateien vor.
; 
; Example:
; idl
; IDL> .r ImaxInterpSolarEvolV5.pro
; IDL> ImaxInterpSolarEvolV5
;******************************************************************************

;******************************************************************************
;* Reads all data of an IMaX cycle either from a 4D cube or from single
;* Stokes images.
;******************************************************************************
FUNCTION ReadData,SrcDir,SrcFile,Header,w,h,NWL,NrStokes

   if STRMID(SrcFile,0,1,/REVERSE_OFFSET) EQ '/' then begin ; single Stokes images
      ObsCyc = STRMID(SrcFile,0,STRLEN(SrcFile)-1)
      im = fltarr(w,h,NWL,NrStokes)
      FileName = SrcDir+SrcFile+'reduc_rnr_'+ObsCyc+'_m080_I.fits.gz' & print,'Reading file ',FileName & im[*,*,0,0] = readfits(FileName,Header)
      FileName = SrcDir+SrcFile+'reduc_rnr_'+ObsCyc+'_m080_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,0,1] = readfits(FileName,Header)
      FileName = SrcDir+SrcFile+'reduc_rnr_'+ObsCyc+'_m080_U.fits.gz' & print,'Reading file ',FileName & im[*,*,0,2] = readfits(FileName,Header)
      FileName = SrcDir+SrcFile+'reduc_rnr_'+ObsCyc+'_m080_V.fits.gz' & print,'Reading file ',FileName & im[*,*,0,3] = readfits(FileName,Header)
      FileName = SrcDir+SrcFile+'reduc_rnr_'+ObsCyc+'_m040_I.fits.gz' & print,'Reading file ',FileName & im[*,*,1,0] = readfits(FileName,Header)
      FileName = SrcDir+SrcFile+'reduc_rnr_'+ObsCyc+'_m040_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,1,1] = readfits(FileName,Header)
      FileName = SrcDir+SrcFile+'reduc_rnr_'+ObsCyc+'_m040_U.fits.gz' & print,'Reading file ',FileName & im[*,*,1,2] = readfits(FileName,Header)
      FileName = SrcDir+SrcFile+'reduc_rnr_'+ObsCyc+'_m040_V.fits.gz' & print,'Reading file ',FileName & im[*,*,1,3] = readfits(FileName,Header)
      FileName = SrcDir+SrcFile+'reduc_rnr_'+ObsCyc+'_p040_I.fits.gz' & print,'Reading file ',FileName & im[*,*,2,0] = readfits(FileName,Header)
      FileName = SrcDir+SrcFile+'reduc_rnr_'+ObsCyc+'_p040_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,2,1] = readfits(FileName,Header)
      FileName = SrcDir+SrcFile+'reduc_rnr_'+ObsCyc+'_p040_U.fits.gz' & print,'Reading file ',FileName & im[*,*,2,2] = readfits(FileName,Header)
      FileName = SrcDir+SrcFile+'reduc_rnr_'+ObsCyc+'_p040_V.fits.gz' & print,'Reading file ',FileName & im[*,*,2,3] = readfits(FileName,Header)
      FileName = SrcDir+SrcFile+'reduc_rnr_'+ObsCyc+'_p080_I.fits.gz' & print,'Reading file ',FileName & im[*,*,3,0] = readfits(FileName,Header)
      FileName = SrcDir+SrcFile+'reduc_rnr_'+ObsCyc+'_p080_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,3,1] = readfits(FileName,Header)
      FileName = SrcDir+SrcFile+'reduc_rnr_'+ObsCyc+'_p080_U.fits.gz' & print,'Reading file ',FileName & im[*,*,3,2] = readfits(FileName,Header)
      FileName = SrcDir+SrcFile+'reduc_rnr_'+ObsCyc+'_p080_V.fits.gz' & print,'Reading file ',FileName & im[*,*,3,3] = readfits(FileName,Header)
      FileName = SrcDir+SrcFile+'reduc_rnr_'+ObsCyc+'_p227_I.fits.gz' & print,'Reading file ',FileName & im[*,*,4,0] = readfits(FileName,Header)
      FileName = SrcDir+SrcFile+'reduc_rnr_'+ObsCyc+'_p227_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,4,1] = readfits(FileName,Header)
      FileName = SrcDir+SrcFile+'reduc_rnr_'+ObsCyc+'_p227_U.fits.gz' & print,'Reading file ',FileName & im[*,*,4,2] = readfits(FileName,Header)
      FileName = SrcDir+SrcFile+'reduc_rnr_'+ObsCyc+'_p227_V.fits.gz' & print,'Reading file ',FileName & im[*,*,4,3] = readfits(FileName,Header)
   endif else begin ; 4D cube
      print,'Reading file ',SrcDir,SrcFile
      im = readfits(SrcDir+SrcFile,Header)
   endelse
   
   RETURN,im
END

;******************************************************************************
;* Saves all data of an IMaX cycle either as a 4D cube or as single
;* Stokes images.
;******************************************************************************
FUNCTION SaveData,im,DstDir,SrcDir,SrcFile,Header

   if STRMID(SrcFile,0,1,/REVERSE_OFFSET) EQ '/' then begin ; single Stokes images
      ObsCyc = STRMID(SrcFile,0,STRLEN(SrcFile)-1)
      SPAWN,'mkdir -p '+DstDir+SrcFile & FileName = DstDir+SrcFile+'reduc_rnr_'+ObsCyc+'_m080_I.fits' & print,'Saving file ',FileName,'.gz' & DstImg = REFORM(im[*,*,0,0]) & writefits,FileName,DstImg,Header & SPAWN,'gzip -f '+FileName
      SPAWN,'mkdir -p '+DstDir+SrcFile & FileName = DstDir+SrcFile+'reduc_rnr_'+ObsCyc+'_m080_Q.fits' & print,'Saving file ',FileName,'.gz' & DstImg = REFORM(im[*,*,0,1]) & writefits,FileName,DstImg,Header & SPAWN,'gzip -f '+FileName
      SPAWN,'mkdir -p '+DstDir+SrcFile & FileName = DstDir+SrcFile+'reduc_rnr_'+ObsCyc+'_m080_U.fits' & print,'Saving file ',FileName,'.gz' & DstImg = REFORM(im[*,*,0,2]) & writefits,FileName,DstImg,Header & SPAWN,'gzip -f '+FileName
      SPAWN,'mkdir -p '+DstDir+SrcFile & FileName = DstDir+SrcFile+'reduc_rnr_'+ObsCyc+'_m080_V.fits' & print,'Saving file ',FileName,'.gz' & DstImg = REFORM(im[*,*,0,3]) & writefits,FileName,DstImg,Header & SPAWN,'gzip -f '+FileName
      SPAWN,'mkdir -p '+DstDir+SrcFile & FileName = DstDir+SrcFile+'reduc_rnr_'+ObsCyc+'_m040_I.fits' & print,'Saving file ',FileName,'.gz' & DstImg = REFORM(im[*,*,1,0]) & writefits,FileName,DstImg,Header & SPAWN,'gzip -f '+FileName
      SPAWN,'mkdir -p '+DstDir+SrcFile & FileName = DstDir+SrcFile+'reduc_rnr_'+ObsCyc+'_m040_Q.fits' & print,'Saving file ',FileName,'.gz' & DstImg = REFORM(im[*,*,1,1]) & writefits,FileName,DstImg,Header & SPAWN,'gzip -f '+FileName
      SPAWN,'mkdir -p '+DstDir+SrcFile & FileName = DstDir+SrcFile+'reduc_rnr_'+ObsCyc+'_m040_U.fits' & print,'Saving file ',FileName,'.gz' & DstImg = REFORM(im[*,*,1,2]) & writefits,FileName,DstImg,Header & SPAWN,'gzip -f '+FileName
      SPAWN,'mkdir -p '+DstDir+SrcFile & FileName = DstDir+SrcFile+'reduc_rnr_'+ObsCyc+'_m040_V.fits' & print,'Saving file ',FileName,'.gz' & DstImg = REFORM(im[*,*,1,3]) & writefits,FileName,DstImg,Header & SPAWN,'gzip -f '+FileName
      SPAWN,'mkdir -p '+DstDir+SrcFile & FileName = DstDir+SrcFile+'reduc_rnr_'+ObsCyc+'_p040_I.fits' & print,'Saving file ',FileName,'.gz' & DstImg = REFORM(im[*,*,2,0]) & writefits,FileName,DstImg,Header & SPAWN,'gzip -f '+FileName
      SPAWN,'mkdir -p '+DstDir+SrcFile & FileName = DstDir+SrcFile+'reduc_rnr_'+ObsCyc+'_p040_Q.fits' & print,'Saving file ',FileName,'.gz' & DstImg = REFORM(im[*,*,2,1]) & writefits,FileName,DstImg,Header & SPAWN,'gzip -f '+FileName
      SPAWN,'mkdir -p '+DstDir+SrcFile & FileName = DstDir+SrcFile+'reduc_rnr_'+ObsCyc+'_p040_U.fits' & print,'Saving file ',FileName,'.gz' & DstImg = REFORM(im[*,*,2,2]) & writefits,FileName,DstImg,Header & SPAWN,'gzip -f '+FileName
      SPAWN,'mkdir -p '+DstDir+SrcFile & FileName = DstDir+SrcFile+'reduc_rnr_'+ObsCyc+'_p040_V.fits' & print,'Saving file ',FileName,'.gz' & DstImg = REFORM(im[*,*,2,3]) & writefits,FileName,DstImg,Header & SPAWN,'gzip -f '+FileName
      SPAWN,'mkdir -p '+DstDir+SrcFile & FileName = DstDir+SrcFile+'reduc_rnr_'+ObsCyc+'_p080_I.fits' & print,'Saving file ',FileName,'.gz' & DstImg = REFORM(im[*,*,3,0]) & writefits,FileName,DstImg,Header & SPAWN,'gzip -f '+FileName
      SPAWN,'mkdir -p '+DstDir+SrcFile & FileName = DstDir+SrcFile+'reduc_rnr_'+ObsCyc+'_p080_Q.fits' & print,'Saving file ',FileName,'.gz' & DstImg = REFORM(im[*,*,3,1]) & writefits,FileName,DstImg,Header & SPAWN,'gzip -f '+FileName
      SPAWN,'mkdir -p '+DstDir+SrcFile & FileName = DstDir+SrcFile+'reduc_rnr_'+ObsCyc+'_p080_U.fits' & print,'Saving file ',FileName,'.gz' & DstImg = REFORM(im[*,*,3,2]) & writefits,FileName,DstImg,Header & SPAWN,'gzip -f '+FileName
      SPAWN,'mkdir -p '+DstDir+SrcFile & FileName = DstDir+SrcFile+'reduc_rnr_'+ObsCyc+'_p080_V.fits' & print,'Saving file ',FileName,'.gz' & DstImg = REFORM(im[*,*,3,3]) & writefits,FileName,DstImg,Header & SPAWN,'gzip -f '+FileName
      SPAWN,'mkdir -p '+DstDir+SrcFile & FileName = DstDir+SrcFile+'reduc_rnr_'+ObsCyc+'_p227_I.fits' & print,'Saving file ',FileName,'.gz' & DstImg = REFORM(im[*,*,4,0]) & writefits,FileName,DstImg,Header & SPAWN,'gzip -f '+FileName
      SPAWN,'mkdir -p '+DstDir+SrcFile & FileName = DstDir+SrcFile+'reduc_rnr_'+ObsCyc+'_p227_Q.fits' & print,'Saving file ',FileName,'.gz' & DstImg = REFORM(im[*,*,4,1]) & writefits,FileName,DstImg,Header & SPAWN,'gzip -f '+FileName
      SPAWN,'mkdir -p '+DstDir+SrcFile & FileName = DstDir+SrcFile+'reduc_rnr_'+ObsCyc+'_p227_U.fits' & print,'Saving file ',FileName,'.gz' & DstImg = REFORM(im[*,*,4,2]) & writefits,FileName,DstImg,Header & SPAWN,'gzip -f '+FileName
      SPAWN,'mkdir -p '+DstDir+SrcFile & FileName = DstDir+SrcFile+'reduc_rnr_'+ObsCyc+'_p227_V.fits' & print,'Saving file ',FileName,'.gz' & DstImg = REFORM(im[*,*,4,3]) & writefits,FileName,DstImg,Header & SPAWN,'gzip -f '+FileName
   endif else begin ; 4D cube
      BaseFileName = STRMID(SrcFile,0,STRLEN(SrcFile)-STRLEN('.fits')-1)
      FileName = DstDir+BaseFileName+'.fits'
      print,'Saving file ',FileName,'.gz'
      writefits,FileName,Dst,Header
      SPAWN,'gzip -f '+FileName
   endelse
   
   RETURN,0
END

;******************************************************************************
;******************************************************************************
PRO ImaxInterpSolarEvolV5

CLOSE,/ALL

NWL=8 ; number of wavelength points

SrcDir   = '/home/tino/nas/mhd/MURaM/MySimulations/MURaM_RUN_1/20171121.30G.ngrey.576x100x576/Stokes5250_degraded_0.145_Bin_4_WithStrayLight4_WithNoise_0.003/'
DstDir   = '/home/tino/nas/mhd/MURaM/MySimulations/MURaM_RUN_1/20171121.30G.ngrey.576x100x576/Stokes5250_degraded_0.145_Bin_4_WithStrayLight4_WithNoise_0.003_Interp/'
w        = 144
h        = 144
NWL      = 5
NrStokes = 4
SrcFile  = MAKE_ARRAY(59+2,/STRING)
FileNo   = 0
; Note: the first and the last file has to be duplicated in order to get as many output files as input files
SrcFile[FileNo] = '254160/' & FileNo++
SrcFile[FileNo] = '254160/' & FileNo++
SrcFile[FileNo] = '254400/' & FileNo++
SrcFile[FileNo] = '254640/' & FileNo++
SrcFile[FileNo] = '254880/' & FileNo++
SrcFile[FileNo] = '255120/' & FileNo++
SrcFile[FileNo] = '255360/' & FileNo++
SrcFile[FileNo] = '255600/' & FileNo++
SrcFile[FileNo] = '255840/' & FileNo++
SrcFile[FileNo] = '256080/' & FileNo++
SrcFile[FileNo] = '256320/' & FileNo++
SrcFile[FileNo] = '256560/' & FileNo++
SrcFile[FileNo] = '256800/' & FileNo++
SrcFile[FileNo] = '257040/' & FileNo++
SrcFile[FileNo] = '257280/' & FileNo++
SrcFile[FileNo] = '257520/' & FileNo++
SrcFile[FileNo] = '257760/' & FileNo++
SrcFile[FileNo] = '258000/' & FileNo++
SrcFile[FileNo] = '258240/' & FileNo++
SrcFile[FileNo] = '258480/' & FileNo++
SrcFile[FileNo] = '258720/' & FileNo++
SrcFile[FileNo] = '258960/' & FileNo++
SrcFile[FileNo] = '259200/' & FileNo++
SrcFile[FileNo] = '259440/' & FileNo++
SrcFile[FileNo] = '259680/' & FileNo++
SrcFile[FileNo] = '259920/' & FileNo++
SrcFile[FileNo] = '260160/' & FileNo++
SrcFile[FileNo] = '260400/' & FileNo++
SrcFile[FileNo] = '260640/' & FileNo++
SrcFile[FileNo] = '260880/' & FileNo++
SrcFile[FileNo] = '261120/' & FileNo++
SrcFile[FileNo] = '261360/' & FileNo++
SrcFile[FileNo] = '261600/' & FileNo++
SrcFile[FileNo] = '261840/' & FileNo++
SrcFile[FileNo] = '262080/' & FileNo++
SrcFile[FileNo] = '262320/' & FileNo++
SrcFile[FileNo] = '262560/' & FileNo++
SrcFile[FileNo] = '262800/' & FileNo++
SrcFile[FileNo] = '263040/' & FileNo++
SrcFile[FileNo] = '263280/' & FileNo++
SrcFile[FileNo] = '263520/' & FileNo++
SrcFile[FileNo] = '263760/' & FileNo++
SrcFile[FileNo] = '264000/' & FileNo++
SrcFile[FileNo] = '264240/' & FileNo++
SrcFile[FileNo] = '264480/' & FileNo++
SrcFile[FileNo] = '264720/' & FileNo++
SrcFile[FileNo] = '264960/' & FileNo++
SrcFile[FileNo] = '265200/' & FileNo++
SrcFile[FileNo] = '265440/' & FileNo++
SrcFile[FileNo] = '265680/' & FileNo++
SrcFile[FileNo] = '265920/' & FileNo++
SrcFile[FileNo] = '266160/' & FileNo++
SrcFile[FileNo] = '266400/' & FileNo++
SrcFile[FileNo] = '266640/' & FileNo++
SrcFile[FileNo] = '266880/' & FileNo++
SrcFile[FileNo] = '267120/' & FileNo++
SrcFile[FileNo] = '267360/' & FileNo++
SrcFile[FileNo] = '267600/' & FileNo++
SrcFile[FileNo] = '267840/' & FileNo++
SrcFile[FileNo] = '268080/' & FileNo++
SrcFile[FileNo] = '268080/' & FileNo++
NrFiles = FileNo

SPAWN,'mkdir -p '+DstDir

f=(-0.5+dindgen(NWL)/FLOAT(NWL)+1.0/(2.0*FLOAT(NWL)))

; lies FITS-Dateien, um auch an die Header-Informationen zu kommen
Src1 = ReadData(SrcDir,SrcFile[0],Header1,w,h,NWL,NrStokes) ; time t-1
Src2 = ReadData(SrcDir,SrcFile[1],Header2,w,h,NWL,NrStokes) ; time t

for FileNo=1,NrFiles-2 do begin
   Src3 = ReadData(SrcDir,SrcFile[FileNo+1],Header3,w,h,NWL,NrStokes) ; time t+1
   
   Dst = Src2
   for j=0,NWL-1 do begin
     if(f[j] gt 0) then begin
       Dst[*,*,j,*] = (1.0-f[j])*Src2[*,*,j,*] + f[j]*Src1[*,*,j,*]
     end else begin
       Dst[*,*,j,*] = (1.0+f[j])*Src2[*,*,j,*] - f[j]*Src3[*,*,j,*]
     end
   end

   ret = SaveData(Dst,DstDir,SrcDir,SrcFile[FileNo],Header2)

   ; prepare next loop run
   Src1 = Src2 & Header1 = Header2
   Src2 = Src3 & Header2 = Header3
end

end
