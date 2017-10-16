pro filtering_fringes,initima,finima,perc,weight,ima_single,imafff_single
                      ;real_ori,imaginario_ori,real,imaginario,ima,imafff

;Filtra las franjas en una secuencia de imagenes cuyos nos. de 
;identificacion van desde "initima" hasta "finima".
;Tambien puede filtrar una imagen individual introducida como argumento en la
;llamada de la rutina: "ima_single"
;Basados en la funcion "weight" generada con "masking_power.pro", 
;se filtran en la Transformada de Fourier de cada imagen (partes Real e 
;Imaginaria por separado).

;INPUTS:
;       initima, finima = nos. de identificacion de la primera
;                         y ultima imagenes que se van a filtrar.
;       perc = porcentaje de apodizacion en el calculo de los
;              espectros de potencia. Si perc=0 no se apodiza.
;              Se recomienda perc = 10 o 12 % para evitar el patron
;              de replicacion de la imagen (dreadfull cross).
;       weight = Weighting function (generated with masking_power.pro) 
;                to filter the spurious signal in the images 
;                This function must be < or =1 (origen de coords
;                en el centro de la matriz).
;       ima_single = Solo utilizar en caso de procesar una sola imagen
;              que se introducira a traves de este argumento.
;              En este caso los valores asignados a (initima,finima) son
;              ignorados por el programa ya que no se leen imagenes desde
;              el disco.
;              La salida en este caso es a traves del argum."imafff_single"
;              No utilizar en la llamada estos dos argumentos en caso de
;              que se procese una serie leida desde disco.

;OUTPUT:
;       Imagenes que resultan del filtrado de las franjas.
;       Se depositan en disco mediante ficheros .fits
;       Si solo se procesa una imagen introducida como argumento ("ima_single")
;       la salida no se deposita en disco y se llama: "imafff_single"

;!!! ATENCION: Hay parametros a cambiar por editor!!!!!

;_____________________________________________________________
;__________________________________________________________________
;Editado el 9 de Junio de 2009 (Esrange)
;Jose A. Bonet
;________________________________________________________________

;!!!!!! Parametros a cambiar por editor !!!!!!
; Estos parametros son ignorados en caso de procesar una sola imagen
; a traves de la llamada "ima_single"

path='/scratch/TALLERK/WIENER/FF_IMAGES/09-06-09_set1/'
namin='ffcam1.'
namout='fffcam_1.'
;namin='ffcam2.'
;namout='fffcam_2.'
;Vertices de la caja a extraer. Si se ponen todos a 0 entonces toma la imagen completa
;Las dimensiones finales de la caja extraida tiene que ser iguales a las de "weight".
;lox=54 & hix=54+951 & loy=40 & hiy=40+951   
lox=0 & hix=0 & loy=0 & hiy=0   ;para trabajar con la imagen completa
;_________________________________________________________________

;
;Determine dimensions and create apodization window
;
dimx=(size(weight))(1) & dimy=(size(weight))(2)
dimx2=dimx/2  &  dimy2=dimy/2
www=win_cos(dimx,dimy,perc)
twww=total(www)

;Filtrado de franjas

if keyword_set(ima_single) then begin  ;Caso de una sola imagen a traves de la llamada

  ima=ima_single
  ;Calculo de FFT de las imagenes y aplicacion de la mascara 
  sus=total(www*ima)/twww
  fima=shift(fft((ima-sus)*www,-1),dimx2,dimy2)
  ;real_ori=float(fima)
  real=float(fima)*weight
  ;imaginario_ori=imaginary(fima)
  imaginario=imaginary(fima)*weight
;  statist,fima  & statist,real  &  statist,imaginario

  ;Transformada inversa
  fima=complex(real,imaginario)
  imafff=float(fft(shift(fima,-dimx2,-dimy2),1))+sus
;  imafff=fix(round(imafff)) ;converts to single-precision integers to save space
  imafff_single=imafff

endif else begin           ;Caso de una serie de imags.leidas desde disco duro

  for k=initima,finima do begin

    ;Lectura de las imagenes y estraccion de una caja si procede
    n_image=string(k,format='(I3.3)')
    name=path+namin+n_image+'.fits'
;    print,name
    im_name=findfile(name)
    if (size(im_name))(1) ne 1 then begin 
      print,''
      print,'No existe imagen '+n_image+' o hay mas de una con este No.'
      goto,eend
    endif
    print,'Reading image: ',im_name
    ima=readfits(im_name,header)
    ima=float(ima)
    flagg= lox ne 0 or hix ne 0 or loy ne 0 or hiy ne 0
    if flagg then ima=float(ima(lox:hix,loy:hiy))
;   statist,ima

    ;Calculo de FFT de las imagenes y aplicacion de la mascara 
    sus=total(www*ima)/twww
    fima=shift(fft((ima-sus)*www,-1),dimx2,dimy2)
    ;real_ori=float(fima)
    real=float(fima)*weight
    ;imaginario_ori=imaginary(fima)
    imaginario=imaginary(fima)*weight
;    statist,fima  & statist,real  &  statist,imaginario

    ;Transformada inversa
    fima=complex(real,imaginario)
    imafff=float(fft(shift(fima,-dimx2,-dimy2),1))+sus

    ;Almacena imagen filtrada como integer
    imafff=fix(round(imafff)) ;converts to single-precision integers to save space
    nameo=path+namout+n_image+'.fits'
    print,'Writing image fits: ',nameo
    writefits,nameo,imafff,header

  endfor

endelse

eend:

end
