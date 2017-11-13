;deconvolution_IMaX7.pro

@lib_wiener
@zernike_annular_IMaX
;______________________________________________________________________

function weight_f,pow,mask

;Construccion del filtro de franjas, "weight", especifico para cada imagen 
;individual de cualquier scan polarimetrico. 
;Para ello se utiliza la mascara "mask" calculada (programa masking_power.pro) 
;a partir del espectro de potencia medio "pow" de un scan polarimetrico 
;completo en el caso mas desfavorable (i.e. cuando dicho scan esta muy 
;alejado del flatfield).

;Ejemplo:       weight = weight_f(pow,mask)
;_____________________________________________________
;Editado el 11 de Enero de 2010
; Jose A. Bonet
;_____________________________________________________

;!!!!!! PARAMETERS TO BE CHANGED BY EDITOR !!!!!!!!

fwhm = 7  ;0 ;17;0 ;1   ;7;17  ;smoothing param to smooth the input power
;_________________________________________________________

  dimx=(size(pow))(1)  &  dimy=(size(pow))(2)
  dimx2=dimx/2  &  dimy2=dimy/2
  xl=dimx2-25 & xh=dimx2+25 & yl=dimy2-25 & yh=dimy2+25 ;central area of "pow"

  if fwhm ne 0 then begin
    pows=sconvol(pow,fwhm=fwhm)
    ;pows(xl:xh,yl:yh)=pow(xl:xh,yl:yh) ;Linea alternativa a las 4 siguientes.
    masc=mask*0+1                       ;Preserva sin smoothing la potencia de
    masc(xl:xh,yl:yh)=mask(xl:xh,yl:yh) ;las franjas de baja freq. localizadas
    mascc=masc eq 0                     ;en areas muy peque�as. Asi se evita
    pows=pows*masc+pow*mascc            ;que la se�al se debilite en el smoothing.
    pows=(pows+especul(pows))/2.
  endif else begin
    pows=pow
  endelse
  pows=pows>(1.e-10)

;Create a weighting function to filter the power in the areas defined by "mask".

  ;Cartesian coords of the points within the areas masked with zeroes.
  wg=where(mask eq 0,nn)
  if nn ne 0 then begin
    x=lonarr(nn)  &  y=x 
    for i=0l,nn-1 do begin
       x(i)=wg(i) mod dimx
       y(i)=fix(wg(i)/dimx)
    endfor
  endif

  ;Masking the power spectrum
  ;poww=pows*float(mask)                   ;for linear scale !!!!!!!!
  lpows=alog10(pows)                       ;for log10 scale !!!!!!!!
  poww=lpows*float(mask)                   ;for log10 scale !!!!!!!!

  ;Interpolation in the gaps containing zeroes.
  if nn ne 0 then begin
    for i=0l,nn-1 do begin  ;filling the gaps in poww by interpolation
       xx=x(i)  &  yy=y(i)
       fill_gaps,float(poww),xx,yy,valor
       poww(xx,yy)=valor
    endfor
  endif
  poww=(poww+especul(poww))/2.

  ;Weigting factor for filtering the signal
  ;weight=poww/pows                      ;for linear scale !!!!!!!!
  weight=10.^(poww-lpows)                ;for log10 scale !!!!!!!!
  ;weight=median(weight,3)

  ;Values larger than 1 are set to 1.
  w=where(weight gt 1,nnn)
  print,'Number of points where the weigting function is >1 ',nnn
  print,'these points will be set to 1'
  if nnn ne 0 then begin
    weight(w)=1
  endif
;stop
return,weight
end
;____________________________________________________________________

function deconvolution_IMaX7,im0,c_medt,inter,flag,imfd,maskf=maskf,maskg=maskg
        ;,filter,filt_rest,t0,mtf0,tsupport

; Deconvoluciona una imagen de su PSF mediante un filtro de Wiener-Helstrom Modificado
; Difiere de "deconvolution_IMaX4.pro" en que aqui la funcion de peso para filtrar las 
; franjas se construye especificamente para cada imagen (rutina "weight_f.pro"). 
; Como en "deconvolution_IMaX4.pro", tanto el filtrado de franjas y/o de los granos
; del polvo, asi como como la restauracion se hacen sobre la imagen global.

;INPUTS:
;       im0 = imagen a deconvolucionar. Ha de ser cuadrada y con dimensiuones 
;             multiplo de 4.No se precisa que este normalizada; la imagen 
;             restaurada que resulta respeta el valor medio de la imagen im0.
;       c_medt = coeffs. de Zernike que caracterizan las aberraciones promedio
;              (i.e. en todo el FOV), a partir de c2 (tip). [radianes]
;       inter = Opcion de interacion para refinar los parametros de regularizacion.
;              (No interacion --> inter=0)
;       flag = (0) Only removal of fringes and/or dust remains (in combin.with the Keywords)
;              (1) the same as flag=0 plus noise filtering
;              (2) the same as flag=1 plus image deconvolution
;       imfd = Intermediate result: image after fringes and/or dust removal.
;              In case of flag=0 the main result of the routine ("imr") coincides with imfd

;       NOTE.- The removal of fringes and/or dust requires the setting of, at least, one
;              of the following keywords. If any of these keyword are set, then only
;              noise-filtering (for flag=1) or noise-filtering & restoration (for flag=2) 
;              will be performed.

;       KEYWORDS:  maskf  &  maskg
;       maskf = Binary mask set to zero in the areas with spurious power (fringes) and
;              to one in the rest. This mask is constructed from the mean power of a
;              polarimetric scan in the most unfavourable case (i.e. taken far away
;              from the nearest faltfield) (origin of coords in center of the matrix).
;       maskg  = Binary mask in the real domain, defining with 1s the image points
;                affected by rubbish (dust and scratches). It can be the simple
;                mask --> "maskg" or the dilated one by an "opening operation" --> "maskgd"

;OUTPUTS:
;           imr = deconvolution_IMaX7(im0,c_medt,inter,flag,imfd,maskf=mask,maskg=maskg)
;                             ;,filter,filt_rest,t0,mtf0,tsupport)
;                 siendo imr la imagen restaurada
;       filter = "filtro optimo" de ruido construido a partir del espectro de
;                potencia medio "pow".
;       filt_rest = filtro de restauracion que combina el filtro de ruido con la
;                   deconvolucion.
;       t0 y mtf0 = OTF y MTF del sistema optico, respectivamente
;       tsupport = mascara que define el dominio de frecuencias hasta la cut-off.

;             !!!! HAY PARAMETROS A CAMBIAR POR EDITOR !!!!
;________________________________________________________________________________
; Editado el 19 de Enero de 2010 
; Construido a partir de "deconvolution_IMaX4.pro" en el que se ha cambiado la
; keyword input "weight" por "maskf" y la funcion de peso para filtrar las franjas aqui
; se construye especificamente para cada imagen (rutina "weight_f.pro")
; Tambien se ha hecho una reorganizacion cosmetica con respecto al programa original
; ("deconvolution_IMaX4.pro")
; Jose A. Bonet
;_______________________________________________________________________________

; !!! A MODIFICAR POR EDITOR: SETTINGS !!!

  ;Initial values for the "regularization parameters"
  low_f=0.1 ;(0.3 M.L.);(Ref.values: 0.1,0.2) Lower boundary in construc.of noise filter.
  filterfactor=1.25;(Ref.values: 1.,1.1,1.2) Modify level of noise: filter more conservative
  rfac=0.9  ;1.0;reduction factor of the cutoff freq.to avoid noise amplification in the restoration
           ;It is only used for the purpose of constructing the filter of noise.
  regdelta=0.005 ;(Ref.values:0.,0.001) Modulate restoration at high freq.(Restoration filter)
  fwhm1=5   ;(Ref.value:5) FWHM convolving Gaussian (1st smooth.in construction of noise filter).
  fwhm2=11  ;(Ref.value:11) FWHM convolving Gaussian (2nd smooth.in construction of noise filter).
  fwhm = 7 ;(Ref.value: 7) smoothing param to smooth the input power. 
  rbox1 = 7  ;(Ref.value: 7) size of running box for smoothing the merging masks in dust removal.
  rbox2 = 9  ;(Ref.value: 9) size of running box for image smoothing in the dust removal.

  ;Other initial parameters
  perc=10.  ;12.5     ; Porcentage de apodizacion
  d_pix_cut=1.  ;0.  ;Cte.a restar para disminuir el valor efectivo de la frec.cut_off
  noiseapo=1 ;0  ;Determinacion potencia de ruido en imagen con/sin apodizacion (1/0)
  disp = 1;0       ; Optional display de nivel de ruido y contrastes: Yes (1), No (0)
  ccd = 1;0                 ;Optional inclusion of the CCD effect in the inversion
  telescope_D = 0.98 ;0.968(patch 128); 0.976(patch 256); 0.98(full image ~936);1. ; Telescope aperture width (m)
  cobs=0.327 ;(media alemanes+Carmen); 0.332 (alemanes) ;0.321 (Carmen) ;0.352 (Obsol);Central Obscuration (fraccion de telescope_D) (=0 for clear aperture)
  telescope_f = 45.00	    ; Telescope focal length (m) (at ~ 525 nm)
  lambda      = 5250.2e-10      ; Working wavelength (m)
  pix=0.055      ; Pixel-size en arcsec
  pix_CCD=pix    ; Pixel-size in the CCD (arcsec)
  s_CCD=pix_CCD  ; Sampling interval in the CCD (arcsec)

;Lee las mascaras de Mikey-mouse para deformaciones de soportes de M1 en SUNRISE
;restore,'/scratch/TALLERK/MICKEY_MASK/mickeymask952.save' ;carga mik_r45 y phi_m1 !!!!!
mik=0 ;mik_r45        ; =0 si no quiero topos !!!!!
spid=11  ;1(patch 128) ; 3(patch 256) ; ~ 11 (patch de 952 � 936 pix)!!!!! 
;__________________________________________________________________________

; Parameters and functions derived from the previous information

  siz=(size(im0))(1)
  jmax = 1 + n_elements(c_medt)
  pixrad = pix*!pi/180./3600.                   ; Pixel-size en radian.
  nu_cutoff=telescope_D/lambda                  ; Frecuencia de corte en rad^-1
  deltanu=1./siz/pixrad                         ; Intervalo de muestreo en rad^-1
  rpupil=nu_cutoff/deltanu/2.                   ; Radio de la pupila en pixels
  ;rpupil=nint(nu_cutoff/deltanu/2.)    ; Radio de la pupila en pixels (redondeo)
  pix_cut=2*rpupil-d_pix_cut    ;pixel en la frecuencia cut-off.(ojo!! -d_pix_cut)
  lim_freq = 2*rpupil*rfac  ;cut-off reducida solo a efectos de construccion filtro de ruido

  if cobs eq 0 then begin   ;clear aperture
    ze=zernike2(siz/4,rpupil,jmax)             ; Define the first jmax Zernike polynom.
  endif else begin
;    ze=zernike_annular(siz/4,rpupil,cobs,jmax) ;Define the first jmax Zernike polynom.
    ze=zernike_annular_IMaX(siz/4,rpupil,cobs,jmax,mik=mik,spid=spid) ;Zernikes in annular+masked aperture. Si mik=0 y spid=0 esta rutina devuelve simplemente los Zernikes anulares (igual que "Zernike_annular.pro".!!!!!
  endelse
  pupil=ze(*,*,0)                        ; Pupil-mask
  ;rd=radius_1(siz/2,pix_cut,tsupport)
  rd=radius_aper2(siz/2,lim_freq,tsupp) 
  tsupp = tsupp ne 0 ;Template defining freq.domain upto the reduced cutoff; this domain 
                     ;differs from that defined by tsupport, and it is only used in the 
                     ;construction of the filter of noise.
  rd=radius_aper2(siz/2,pix_cut,tsupport)
  tsupport=tsupport ne 0               ;Template defining freq. domain upto cutoff
  pmask = tsupport eq 0                  ; Mascara para determinar el nivel 
  n2=siz/2 & n6=siz/6                    ; de potencia del ruido en las imagenes
  pmask(n2-n6:n2+n6,*) = 0               ; enfocada y desenfocada.
  pmask(*,n2-n6:n2+n6) = 0
  www=win_cos(siz,siz,perc)              ; Ventana de apodizacion
  siz_central=siz-nint(2.*float(siz)*perc/100.) ;size of central area non-apodized
  low=(siz-siz_central)/2 & hi=low+siz_central-1   ; Def. sub-area central.
  ;tairy=float(otfunc(sqrt(pupil),norm))*tsupport     ;MTF of a diffrac.limited telescope

  print,''
  print,nu_cutoff,!radeg*3600./nu_cutoff,  $
      format='("Diffrac.cut-off freq. =",E11.4," rad^-1 (",F5.3," arcsec)")'
  print,rpupil,format='("Pupil-radius  =",F7.3," pixel")'

;stop

  otfccd=1.
  if ccd ne 0 then begin  ;Inclusion opcional del efecto de pixelacion CCD
    ax=pix_CCD/pix  &  ay=ax
    ;sx=s_CCD/pix
    sx=0  ;solo se considera la "OTF(footprint)"
    sy=sx
    otfccd = otf_CCD(siz,siz,ax,ay,sx,sy) * tsupport
  endif

;.............................................................................

; Calculo del contraste original.

  ima0=im0
  m0=mean(ima0(low:hi,low:hi))
  contrast_i=stdev(ima0(low:hi,low:hi))/m0  ;Contrast central part original image
  if disp ne 0 then begin
    print,''
    print,m0,contrast_i,$
     format='("Original image: Mean (central part)= ",F10.3,"  Contrast_rms= ",F7.5)'
  endif
;  print,'Statistics of the original image'
;  statist,ima0

; Transformada de Fourier y espectro de potencia.Desplazamiento origen al centro del array.

  twww=total(www)
  sus0=total(www*ima0)/twww ; cte.a restar para tener media 0 en imag.apodiz.(C)
  d0=fft((ima0-sus0)*www,-1) & d0(0,0)=complex(0.,0.)
  d0=shift(d0,siz/2,siz/2)
  power=(abs(d0))^2
  d1 = d0
; Filtrado de franjas optativo

  if keyword_set(maskf) then begin
    weight=weight_f(power,maskf)   ;construye filtro de franjas especifido para cada frame
    real=float(d0)*weight  ;removal of the fringes.
    imaginario=imaginary(d0)*weight
    d0=complex(real,imaginario)
    power=power * weight
  endif

  ima = float(fft(shift(d0,-siz/2,-siz/2),1))  ;image filtered from fringes
  
; Computation of the noise power level. Smoothing of the power spectrum

  pows=sconvol(power,fwhm=fwhm)
  ;pows=power
  pows=(pows+especul(pows))/2.
  ;print,''
  print,'Nivel de ruido a partir del espec.de potencia'
  level_of_noise_IMaX,pows,siz,disp,pmask,npow0,noise_sigma0

; Genera imagen o plot del espectro de potencia suavizado

  goto,noplot   ;Produccion de plot e imagen del espectro de potencia suavizado
    xnu=findgen(siz/2)/siz/pix
    plot_io,xnu,pows(siz/2:*,siz/2)>1.e-7,xtitle='frequency  [ arcsec!U -1!N ]',$
        ytitle='Power',yrange=[1.e-5,1.e3],/yst,charsize=1.6,back=255,col=0
    oplot,xnu,pows(siz/2,siz/2:*)>1.e-7,linestyle=3,col=0
    oplot,xnu,fltarr(siz/2)+filterfactor*npow0,linestyle=0,col=0  ;,thick=2
    oplot,xnu,fltarr(siz/2)+npow0,linestyle=0,col=0  ;,thick=2
    xyouts,4,10^2,'X-section:  ________',col=0,charsize=1.6
    xyouts,4,10^1.5,'Y-section:  . _ . _ .',col=0,charsize=1.6
    xyouts,4,10,'noise level: horizontal solid line',col=0,charsize=1.6
    imaplot=tvrd()
    write_tiff,'power.tif',reverse(imaplot,2)
    window,2,xs=siz,ys=siz
    tvscl,alog10(pows)*tsupport
    imaplot=tvrd()
    write_tiff,'power2.tif',reverse(imaplot,2)
  noplot:

;  Limpia de motas de polvo y scratches (Optativo)

  if keyword_set(maskg) then begin
    ;Crea mascara complementaria
    maskk = maskg eq 0
    mask_exterior=smooth(float(maskk),rbox1)
    mask_interior=smooth(float(maskg),rbox1)

    ;Coordenadas cartesianas de los puntos en las areas de la mascara con 0s.
    dimx=(size(maskk))(1) & dimy=(size(maskk))(2)
    wg=where(maskk eq 0,nn)
    if nn ne 0 then begin
      x=lonarr(nn)  &  y=x 
      for i=0l,nn-1 do begin
        x(i)=wg(i) mod dimx
        y(i)=fix(wg(i)/dimx)
      endfor
    endif

    ima=ima*maskk   ;ambos arrays han de tener las mismas dimensiones

    ;Interpola en los gaps que contienen ceros
    if nn ne 0 then begin
      for i=0l,nn-1 do begin
        xx=x(i)  &  yy=y(i)
        fill_gaps,ima,xx,yy,valor
        ima(xx,yy)=valor
        ;print,xx,yy,valor
      endfor
    endif
    imas= smooth(ima,rbox2)   ;sconvol(ima,7)
    ima=ima*mask_exterior+imas*mask_interior
    ima=median(ima,3)
  endif

; Intermediate result (up to here only findges and/or dust removal) always available
; as an output.
  imfd = ima+sus0

  if flag eq 0 then begin  ;only removal of fringes and/or dust in final result
    imr=imfd
    goto,fin
  endif

  if keyword_set(maskg) then begin
    d0=fft(ima,-1) & d0(0,0)=complex(0.,0.) ;apodiz.is not needed because it is included in ima.
    d0=shift(d0,siz/2,siz/2)
  endif

; Crea la fase aberrada

  ctel=c_medt(2:*)   ;extrae los coeffs a partir del desenfoque
  c=fltarr(jmax+1)
  c(4:*)=ctel
  phi_rms = sqrt(total(c(4:*)*c(4:*)))   ; rms de Phi sin tip/tilt
  if disp ne 0 then begin
    print,''
    print,phi_rms,format='("Wave-front aberration:  phi_rms [rad]=",F5.3,$)'
    print,1./(phi_rms/2./!pi),format='("    phi_rms [waves]= 1/",F5.3)'
  endif
  phi_0=fltarr(siz/2,siz/2)            ; aberraciones de fase
  for ii=4,jmax do begin
    phi_0=phi_0+c(ii)*ze(*,*,ii)
  endfor

; Calcula la OTF y MTF 

  ;h0=pupil*complex(cos(phi_0),sin(phi_0))      ;generalized pupil function (H)
  h0=sqrt(pupil)*complex(cos(phi_0),sin(phi_0))      ;generalized pupil function (H)
  t0=otfunc(h0,norm) * tsupport	             ;Hermitian function  (G.6)
  ;El multiplicar por tsupport en t0 es para suprimir valores de ~ 10^-17
  ;mas alla de la cutoff (ruido numerico en la FFT).
  t0=t0*otfccd      ;Efecto de pixelacion en la CCD (optativo)
  mtf0=abs(t0)
  ;lmtf0=alog10(mtf0 > 1.e-7) * tsupport

; ***** DETERMINACION ITERATIVA DE LOS PARAMETROS DE REGULARIZACION ********

start:

; Construye filtro de ruido

  filter=filtrop_IMaX(pows,siz,npow0,filterfactor,tsupport,fwhm1=fwhm1,fwhm2=fwhm2,low_f=low_f)
  filter = filter*tsupp         ;safety to avoid noise restoration

; Construye filtro de restauracion

  if flag eq 2 then begin
    filt_rest=t0*0.
    wfil=where(tsupport ne 0,nfil)
    ;Regularizacion con constante aditiva en la mtf�
    ;filt_rest(wfil)=filter(wfil)*conj(t0(wfil))/(mtf0(wfil)^2+regdelta)
    ;Regularizacion con cu�a aditiva en la mtf�
    nu_rel=(shift(dist(siz),siz/2,siz/2))/pix_cut ;normaliza distancias a la de corte.
    suma = nu_rel*regdelta*tsupport  ;genera una cu�a creciente hasta la cutoff
    filt_rest(wfil)=filter(wfil)*conj(t0(wfil))/(mtf0(wfil)^2+suma(wfil))
    ;filt_rest(wfil)=filter(wfil)*conj(t0(wfil))/(mtf0(wfil)+suma(wfil))^2 ;(filtra demasiado)
  endif

  if inter ne 0 then begin
    PAUSE
    print,'FILTRO OPTIMO DE RUIDO RESULTANTE'
    tvwinp,(filter* tsupport) + (tsupport eq 0),xr=0.5,yr=0.5,$
           title='FILTRO OPTIMO DE RUIDO'
    ;lfilter = alog10(filter > 1.e-3) * tsupport
    ;tvwinp,lfilter,/free
    if flag eq 2 then begin
      PAUSE
      print,'FILTRO DE RESTAURACION RESULTANTE'
      tvwinp,float(filt_rest)* tsupport + (tsupport eq 0),xr=0.5,yr=0.5,$
             title='FILTRO DE RESTAURACION (PARTE REAL)'
      wait,1
      tvwinp,imaginary(filt_rest)* tsupport + (tsupport eq 0),xr=0.5,yr=0.5,$
           title='FILTRO DE RESTAURACION (PARTE IMAGINARIA)'
    endif
  endif

; Deconvolucion y calculo del contraste imagen restaurada

  case flag of
    1: FFF=filter
    2: FFF=filt_rest
  endcase
  imr = float(fft(shift(d0*FFF,-siz/2,-siz/2),1)) ; total restoration
  imr = imr+sus0
  mr=mean(imr(low:hi,low:hi)) ;& print,'mr=',mr 
  contrast_f=stdev(imr(low:hi,low:hi))/mr ;Contrast central part restor.image
  if disp ne 0 then begin
    print,''
    print,mr,contrast_f,$
     format='("Restored image: Mean (central part)= ",F10.5,"  Contrast_rms= ",F7.5)'
  endif

; Opcion de modificar los parametros de regularizacion

;  wait,1
;  tvwinp,imr,title='VISUALIZE THE RESULT OF THE PRESENT INTERACTION'
  print,''
  print,low_f,filterfactor,regdelta,$
          format='("Actual regul.param.: low_f,filterfactor,regdelta= ",4F8.4)'
  print,''
  if inter ne 0 then begin
    print,'New regularization parameters ? : press y/n'
    tent=get_kbrd(1)
    if (tent eq 'y') then begin
      read,'New values for: low_f,filterfactor,regdelta= ',$
                               low_f,filterfactor,regdelta
      goto,start
    endif
  endif

; ********* FIN ITERACION PARAMETROS DE REGULARIZACION  ******************

fin:

;stop

return,imr

end 
