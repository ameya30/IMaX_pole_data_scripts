@lib_wiener
pro mean_power2,arr4D,perc,imam,pow,powm

;Given a 4D array including a complete spectropolarimetric scan (set of
;images for different wavelengths and polarimetric states),
;this program computes the average of the normalized images and of their
;respective power spectra. Also computes the power spectrum of the average
;image.

;INPUTS:
;       arr4D = 4D array including the images (dimx x dimy) corresponding
;               to several lambdas and polarimetric states:
;                        arr4D(dimx,dimy,lambdas,pol)
;       perc = percentage for apodization in the calculation of the power
;              spectra. For perc=0 no apodization is performed.
;              It is recommended perc=10 or 12 % to avoid the cross
;              in the Fourier domain produced by the replication pattern
;              due to the discrete Fourier Transform nature.

;!!!!!!! There are inputs TO BE CHANGED BY EDITOR !!!!!!!!!!!!!

;OUTPUT:
;       imam = average image
;       pow = average of the power spectra (origin of coords. in the center
;             of the matrix).
;       powm = power spectrum of the average image (origin of coords. in the center
;              of the matrix).
;_____________________________________________________________
;Edited on 13 Oct 2009 (IAC) from mean_power.pro
;Jose A. Bonet
;_____________________________________________________________

;!!!!!! PARAMETERS TO BE CHANGED BY EDITOR !!!!!!!!

;Vertices de la caja a extraer. Si se ponen todos a 0 entonces toma la imagen completa
;lox=54 & hix=54+951 & loy=40 & hiy=40+951   
lox=0 & hix=0 & loy=0 & hiy=0   ;para trabajar con la imagen completa
;________________________________________________________________________

  siz=size(arr4D)
  lamb=siz(3)  &  pola=siz(4)
  numima=lamb*pola

  ind=0
  for p=0,pola-1 do begin
    for l=0,lamb-1 do begin
      ima=float(arr4D(*,*,l,p))
      flagg= lox ne 0 or hix ne 0 or loy ne 0 or hiy ne 0
      if flagg then ima=ima(lox:hix,loy:hiy)
      if ind eq 0 then begin
        ind=1
        dimx=(size(ima))(1)  &  dimy=(size(ima))(2)
        dimx2=dimx/2  &  dimy2=dimy/2
        imam=fltarr(dimx,dimy)
        pow=fltarr(dimx,dimy)
        ;Create apodization window
        www=win_cos(dimx,dimy,perc)
        twww=total(www)
      endif
      ima=ima/mean(ima) ;Normalization of single images
      ;Compute average of normalized images and power spectra
      sus=total(www*ima)/twww
      fima=shift(fft((ima-sus)*www,-1),dimx2,dimy2)
      pow=pow+(abs(fima))^2/float(numima)
      imam=imam+ima/float(numima)
    endfor
  endfor

  ;Compute the power spectrum of the mean image
  print,'Statistics of the average image'
  statist,imam
  sus=total(www*imam)/twww
  fimam=shift(fft((imam-sus)*www,-1),dimx2,dimy2)
  powm=(abs(fimam))^2

;stop

end
