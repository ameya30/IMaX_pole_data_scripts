;pro masking_power
;_____________________________________________________________________________
;_____________________________________________________________________________

function sbox,image,sizbx,sizby,x0,y0,xinf,yinf,xsup,ysup

; Dada una imagen "image" se selecciona una caja de dimensiones prefijadas
; (sizbx,sizby). Para ello se selecciona con el raton el centro de la caja deseada.
; Si por efecto de las dimensiones prefijadas, la caja se sale de la imagen, la 
; rutina retorna el trozo de caja que queda dentro de la imagen y por lo tanto
; se redefinen sus dimensiones.
; La rutina tambien calcula la estadistica de la caja seleccionada.

;INPUTS:
;   image       = imagen en la que se va ha seleccionar la caja.
;   sizbx,sizby = dimensiones de la caja a seleccionar. Son inputs pero pueden
;                 ser outputs si la caja selecionada se sale de la imagen y
;                 hay que redimensionarla a la porcion que queda dentro de ella.

;OUTPUTS:
;	sizbx, sizby = dimensiones de la caja seleccionada en el caso de
;                      redimensionado citado arriba.
;	x0, y0       = coordenadas del centro de la caja seleccionada.
;	xinf,yinf    = coordenadas del vertice inferior izquierdo de la caja.
;	xsup,ysup    = coordenadas del vertice superior derecho de la caja.

;EJEMPLOS:
;               sizbx=40  &  sizby=74
;		box=sbox(image,sizbx,sizby,x0,y0,xinf,yinf,xsup,ysup)

;------------------------------------------------------------------------------
; Editado: 28 de Mayo de 2009 a partir de selectbox.pro
; Jose A. Bonet
;______________________________________________________________________________

  dimx=(size(image))(1) & dimy=(size(image))(2)
  window,0,xs=dimx,ys=dimy
  tvscl,image
  sx=fix(sizbx/2) & sy=fix(sizby/2)  ;mitad tamagno caja a extraer
  mmm=max(image)

  print,''
  print,'****** SELECT CENTER OF THE BOX TO BE EXTRACTED IN IMAGE'
  print,'left button to select and right to quit'
  cursor,x,y,/device,/down
  while (!err eq 1) do begin
    x0 = x & y0 = y
    imtv=image
;   Construimos los bordes brillantes de la caja, para elegirla.
    xinf=x0-sx & xsup=x0+sx-1 & yinf=y0-sy & ysup=y0+sy-1
    if xinf lt 0 then xinf=0  &  if xsup gt dimx-1 then xsup=dimx-1
    if yinf lt 0 then yinf=0  &  if ysup gt dimy-1 then ysup=dimy-1
    imtv(xinf:xsup,yinf)=mmm
    imtv(xinf:xsup,ysup)=mmm
    imtv(xinf,yinf:ysup)=mmm
    imtv(xsup,yinf:ysup)=mmm    
    tvscl,imtv
    print,'new selection or quit'
    cursor,x,y,/device,/down
  endwhile

  sizbx=xsup-xinf+1  &  sizby=ysup-yinf+1
  window,1,xs=sizbx,ys=sizby,xpos=0,ypos=600
  box =  image(xinf:xsup,yinf:ysup)
  tvscl,box

print,'Tamano caja: ',fix(sizbx),fix(sizby),'    Centro caja: ',fix(x0),fix(y0)
print,'Vert.inf.izq.: ',fix(xinf),fix(yinf),'  Vert.sup.dcho.:',fix(xsup),fix(ysup)
print,'Estadistica de la caja extraida'
statist,box

imtv=0b

return,box

end

;***********************************************************************************

;------------------------------------------------------------------------------
; RUBBER.PRO
;------------------------------------------------------------------------------
;
; Herramienta para borrar (poner a cero) o reponer areas previamente
; borradas en una imagen original "imori".
; Se trabaja por cajas de imori para poder utilizar factores de zoom 
; grandes y asi facilitar la tarea de borrado.
; El raton hace de borrador: 
;     left button para borrar; right button para reponer valor original
;     center button to quit.
; El tamaño del rubber-spot es regulable.
; La imagen que resulta del proceso, "imcorr", puede realimentar la rutina
; a traves de la Keyword IMCORR=imcorr, de manera que en otra sesion se 
; pueda continuar el proceso en donde se dejo.
;
;INPUTS:
;      imori = imagen original a retocar. Nunca se destruye ni sustituye. 
;              Incluso si se aplica la rutina reiteradamente haciendo 
;              feedbacks, imori es siempre la misma

; !!!! El tamaño de la caja de trabajo se puede MODIFICAR POR EDITOR!!!!!!

;KEYWORDS:
;      [IMCORR=imcorr] = hace de entrada y salida. La primera vez que se 
;              aplica la rutina, la entrada no aporta ningun dato y la salida
;              da la imagen con areas borradas. En proceso de realimentacion 
;              la entrada es la salida anterior y la salida es la nueva 
;              imagen retocada.
;      [MASK=mask] = mascara binaria que define las areas borradas en la 
;                    imagen "imori" con ceros. En procesos con feedback
;                    esta Keyword funciona de forma analoga a IMCORR.
;      [RUBW=rubw] = tamaño del rubber-spot que se asume de entrada. Se
;               puede modificar opcionalmente durante el proceso. En la 
;               salida guarda el ultimo valor utilizado por si se quiere 
;               utilizar en el feeback

; CALLING SEQUENCE:
;	               ruber,imori,[IMCORR=imcorr,MASK=mask,RUBW=rubw]

;________________________________________________________________________
;Editado el 29 Mayo 2009
;Modificacion de cosmetics3_b.pro de Michael Steinegger
;Jose A. Bonet
;______________________________________________________________________

pro rubber,imori,IMCORR=imcorr,MASK=mask,RUBW=rubw

;!!!! A MODIFICAR POR EDITOR !!!!!!!

sizbx0=200  &  sizby0=200  ;tamaño de la de trabajo a extraer
;______________________________________________

rimori=imori
nn=size(imori)  &  nx=nn(1)  &  ny=nn(2)
if keyword_set(imcorr) eq 0 then imcorr=imori
if keyword_set(mask) eq 0 then mask=intarr(nx,ny)+1
if keyword_set(rubw) eq 0 then rubw=1
if rubw lt 1 then rubw=1
;
; display image to be processed and selected box
;
newbox:
zoomflag=0
fzoom=1
sizbx=sizbx0  &  sizby=sizby0
box=sbox(imcorr,sizbx,sizby,x0,y0,xinf,yinf,xsup,ysup)
boxori=imori(xinf:xsup,yinf:ysup)
boxmask=mask(xinf:xsup,yinf:ysup)
rbox=box
rboxori=boxori
rboxmask=boxmask
rsx=sizbx  &  rsy=sizby
;
; Zoom image (optionally)
;
print,'*****************************************'
newzoom:
print,''
print,'CURRENT ZOOM FACTOR to rebin the window: ',fzoom
print,'TYPE New zoom factor (typically 3 and 0 to quit)'
testch=string(3)  &  testch=get_kbrd(1)
flag=testch
if flag ne '0' then begin
   zoomflag=1
   fzoom=fix(testch)
   rsx=fzoom*sizbx  &  rsy=fzoom*sizby
   window,6,xsize=rsx,ysize=rsy,title='ZOOMED BOX'
   rbox=rebin(box,rsx,rsy,/sample) 
   tvscl,rbox
   rboxori=rebin(boxori,rsx,rsy,/sample)
   rboxmask=rebin(boxmask,rsx,rsy,/sample)
   goto,newzoom
endif
mi=min(rboxori)  &  ma=max(rboxori)
;
; Change rubber width (Optionally)
;
print,'*****************************************'
newidth:
print,''
print,'CURRENT RUBBER WIDTH : ',rubw
print,'TYPE New rubber width (Typically 9 and 0 to start the rubber)'
testch=string(3)  &  testch=get_kbrd(1)
flag=testch
if flag ne '0' then begin rubw=fix(testch) & goto,newidth & endif
rubw=rubw*fzoom
;
; "rubber" (in accumulative mode)
;
print,'*****************************************'
print,''
accum:
if zoomflag eq 0 then begin wset,1 & endif else begin & wset,6 & endelse
print,'--- Press left mouse button to erase.'
print,'--- Press right mouse button to fill.'
print,'--- Press middle mouse button to quit.'
loop:  
cursor,x,y,/device
;x1=x  &  x2=x+rubw-1  &  y1=y  &  y2=y+rubw-1
x1=x-fix(rubw/2)  &  x2=x1+rubw-1  &  y1=y-fix(rubw/2)  &  y2=y1+rubw-1
;
; Mark erased/filled pixels in display
;
if y2 gt rsy-1 then y2=rsy-1
if x1 lt 0 then x1=0
if x2 gt rsx-1 then x2=rsx-1
if y1 lt 0 then y1=0

case !err of
1: begin   ;pone areas de la imagen a zero
   rbox(x1:x2,y1:y2)=mi
   rboxmask(x1:x2,y1:y2)=0
;   tvscl,rbox
   tv,bytscl(rbox,min=mi,max=ma)
   end
4: begin   ;repone el valor original de la imagen
   rbox(x1:x2,y1:y2)=rboxori(x1:x2,y1:y2)
   rboxmask(x1:x2,y1:y2)=1
;   tvscl,rbox
   tv,bytscl(rbox,min=mi,max=ma)
   end
2: goto,ende
endcase
goto,loop
ende:
;
; Checking the results by flickering
;
window,3,xs=rsx,ys=rsy
print,'Hit Enter to stop flickering..................................'
flick,bytscl(rbox,min=mi,max=ma),bytscl(rboxori,min=mi,max=ma),1
wdelete,3
print,'*****************************************'
print,''
print,'New changes in the current box ? (y/n)'
testch=string(3)  &  testch=get_kbrd(1)
if testch eq 'y' then goto,accum
;
;Box recover the original size and es encapsulated in the global image
;
if zoomflag eq 1 then begin
   box=rebin(rbox,sizbx,sizby,/sample)
   boxmask=rebin(rboxmask,sizbx,sizby,/sample)
   boxmask=boxmask eq 1
   rubw=rubw/fzoom
   wdelete,6
endif else begin
   box=rbox
   boxmask=rboxmask
endelse

imcorr(xinf:xsup,yinf:ysup)=box
mask(xinf:xsup,yinf:ysup)=boxmask
;
;Opcion de trabajo sobre otra caja de la imagen original
;
print,'*****************************************'
print,''
print,'New box ? (y/n)'
testch=string(3)  &  testch=get_kbrd(1)
if testch eq 'y' then goto,newbox
;
;Fin del proceso
;
imcorr=imcorr*mask+(mask eq 0)*min(imori)
tvwin,imcorr

;stop

end

;***************************************************************************

pro recover_even_function,f_half,f_par

;Dados el 1º y 4º cuadrantes de una funcion par, reconstruye
;la informacion en el 2º y 3er cuadrantes.
;
;INPUT:
;    f_half = imagen con 1º y 4º cuadrantes, de dimension
;             (dimx/2,dimy) y con el origen de coordenadas
;             en (0,dimy/2).
;OUTPUT:
;    f_par  = imagen completa de la funcion par con dimension
;             (dimx,dimy) y origen en (dimx/2,dimy/2)
;             Dado que a informacion en los bordes de "f_half"
;             es incompleta, la primera fila y columna de la
;             matriz "f_par" las llena con zeroes.
;_________________________________________________________
;11 de Junio de 2009  (Esrange)
;Jose A. Bonet
;__________________________________________________________

siz=size(f_half)  &  nx=siz(1)  &  ny=siz(2)
dimx=2*nx  &  dimy=ny
x0=dimx/2  &  y0=dimy/2

f_par=fltarr(dimx,dimy)
f_par(x0,1)=f_half(*,1:ny-1)

;print,f_par

a=f_half(0:nx-1,y0:ny-1)
;print,''
;print,a
a=reverse(reverse(a,2),1)
;print,''
;print,a
f_par(1,1)=a

;print,''
;print,f_par

b=f_half(1:nx-1,1:y0-1)
;print,''
;print,b
b=reverse(reverse(b,2),1)
;print,''
;print,b
f_par(1,y0+1)=b

;print,''
;print,f_par

end

;******************************************************************************

pro masking_power,pow,powcorr,mask,weight

; The problem consists in removing spurious structures in an image (e.g.
; interferential fringes).
; To that aim, we will identify in its power spectrum the areas 
; containing spurious signal and we will create by hand (with the mouse:
; program "rubber.pro") a binary mask filling up such areas with 0's,
; and with 1's the rest.
; From this mask and from the power spectrum a weight function to filter out
; the spurious signal in the images is built.

;INPUTS:
;       pow = power spectrum of the image to be filtered. It is always the
;             the same in an iterative process. The origin of coords.has to
;             be placed at the center of the matrix.

;INPUTS & OUTPUTS:
;       powcorr = power spectrum with areas of spurious signal erased.
;                 It can be an input when one wants to refine the result
;                 stemming from a previous iteration. (Origin of coords.
;                 placed at the center of the matrix).
;       mask = binary mask set to zero in the areas with spurious power and
;              to one in the rest. The same as in the case of "powcorr" this
;              argument may also be an input for iterative improvement of
;              the mask resulting from a previous iteration.(Origin of coords.
;                 placed at the center of the matrix).
;       weight = (only output) Weighting function to filter the
;                spurious signal in the images. This function must be < or =1 
;                (origin of coords in the center of the matrix). 
;                It is derived as the quotient of the interpolated power in 
;                the areas defined by "mask", and the original power. 
;                This function must be < or =1
;                NOTE that the construction of "weight" match the concept of
;                "Optimum-filter of noise" by Brault and White

;!!!!!!!WARNING: THERE ARE PARAMETERS TO BE CHANGED WITH EDITOR!!!!!!!!!

;__________________________________________________________________
;Edited on 9 June 2009 (Esrange)
;Jose A. Bonet
;________________________________________________________________

;!!!!!!!! PARAMETERS TO BE CHANGED WITH EDITOR !!!!!!!!!!!!

option_log=1 ;0  ;Display is wanted in scale: logarithmic/lineal = 1/0
threshold=0.05; 0.1  ;(It is ignored if option_log=1) to display the power.
option_smooth=1  ;0  ;If smoothing of "pow" is desired (YES/NO = 1/0)
fwhm = 7 ;smoothing param to smmoth the input power. It is ignored if option_smooth=0
;____________________________________________________________

dimx=(size(pow))(1) & dimy=(size(pow))(2)
x0=dimx/2  &  y0=dimy/2
xl=x0-25 & xh=x0+25 & yl=y0-25 & yh=y0+25 ;central area of "pow" (new line 12-1-2010)

;Optional smoothing of the input power spectrum "pow"

if option_smooth eq 1 then begin
  pows=sconvol(pow,fwhm=fwhm)
  pows(xl:xh,yl:yh)=pow(xl:xh,yl:yh) ;Preserva sin smoothing la potencia central (new line 12-1-2010)
  pows=(pows+especul(pows))/2.
endif else begin
  pows=pow
endelse

if option_log eq 0 then pows=pows<threshold

;Extracts 1º y 4º quadrants of the power.
powhalf=pows(x0:*,*)
if option_log ne 0 then powhalf=alog10(powhalf)
if option_log eq 0 then begin
  if keyword_set(powcorr) then powcorrhalf=powcorr(x0:*,*)
endif else begin
  if keyword_set(powcorr) then powcorrhalf=alog10(powcorr(x0:*,*))
endelse
if keyword_set(mask) then maskhalf=mask(x0:*,*)

;Create a mask defining the areas of spurious signal
;in the 1st and 4th quadrants of the power spectrum.
rubber,powhalf,IMCORR=powcorrhalf,Mask=maskhalf,RUBW=rubw
if option_log ne 0 then begin
  powcorrhalf=10.^(powcorrhalf)
  powhalf=10.^(powhalf)
endif

;Create even-functions from the 1st and 4th quadrants.
recover_even_function,maskhalf,maskk
mask=maskk*0.+1.
mask(1,1)=maskk(1:*,1:*)
mask=fix(mask)
recover_even_function,powcorrhalf,powc

;Recover the 1st row and column of the original power.
powcorr=pows
powcorr(1,1)=powc(1:*,1:*)

wdelete,0,1

;Create a weighting function to filter the power in the areas defined by the
;"mask" created. Values larger than 1 are set to 1.

  ;Cartesian coords of the points within the areas masked with zeroes.
  wg=where(mask eq 0,nn)
  if nn ne 0 then begin
    x=lonarr(nn)  &  y=x 
    for i=0l,nn-1 do begin
       x(i)=wg(i) mod dimx
       y(i)=fix(wg(i)/dimx)
    endfor
  endif

  ;Masking the smoothed mean power spectrum
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

  w=where(weight gt 1,nnn)
  print,'Number of points where the weigting function are >1 ',nnn
  print,'these points will be set to 1'
  if nnn ne 0 then begin
    weight(w)=1
  endif

;profiles_XY,alog10([[pow],[pow*weight]])
;wait,1
;profiles_XY,alog10((pows-pows*weight)>1.e-14)
profiles_XY,weight

;stop

end



