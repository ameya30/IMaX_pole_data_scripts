function zernike_annular_IMaX,r,r_c,e,jmax,mik=mik,spid=spid

; Retorna un array 3D de dimensiones (2r x 2r x jamx+1) que contiene los
; polinomios de Zernike desde j=1 hasta j=jmax. 
; La matriz para j=0 (es decir (*,*,0)) contiene la mascara que define la
; apertura.

; INPUTS: 
;	r  = semi-dimension of the matrices. 
;	r_c = radius of the aperture.
;	e   = obscurecimiento central expresado en fraccion de "rc" (p.ej. 0.25)
; KEYWORDS:
;       mik = mascara constituida por tres topos de ceros definiendo la posicion de 
;             los soportes de M1 en SUNRISE.
;             Poner mik = 0 u omitir esta keyword si no se quieren topos en la pupila.
;       spid = grueso de los brazos del spider [pixel]. Ha de ser impar.
;              Poner spid = 0 u omitir esta keyword si no se quiere spider en la pupila.

;     NOTA: Si mik=0 y spid=0 esta rutina devuelve simplemente los Zernikes anulares (igual 
;           que "Zernike_annular.pro". !!!!!

; OUTPUTS:
;	     zer = zernike_annular_IMaX(r,r_c,e,jmax,mik=mik,spid=spid)

;____________________________________________________________________________
; 22 de Octubre de 2009 y retocado el 9 de Noviembre de 2009
; Este programa es identico a zernike_annular.pro excepto porque se le han añadido en
; el input dos keywords: "mik" (mascara de los soportes de M1 en SUNRISE) y spid (grosor
; de los brazos del spider) y varias lineas (marcadas con !!!!!!!) en las que 
; opcionalmente se aplica la mascara "mik" como factor y/o se define un spider.
;
; Jose A. Bonet
;____________________________________________________________________________
   
  IF jmax gt 45 THEN BEGIN 
    print,'jmax too large: redefine jmax=45'
    jmax=45
  ENDIF

; Calcula las coordenadas (rd,theta) en cada pixel de la apertura, con
; origen en (r,r). Tambien calcula la mascara que define la apertura.

;  rd=radius_1(r,r_c,support)  ;when fractional pixel at the edge is not consid.
  rd = radius_aper_annular(r,r_c,e,support)     ; (DARIO y JOSE)
  th = theta(r)

  if keyword_set(mik) then begin    ;topos !!!!!!!!!!!!!
    rd = rd*mik                     ;topos !!!!!!!!!!!!!
    support = support*mik           ;topos !!!!!!!!!!!!!
  endif                             ;topos !!!!!!!!!!!!!
  if keyword_set(spid) then begin                   ;spider !!!!!!!!!!!!!
    as=fix(spid/2.)                                 ;spider !!!!!!!!!!!!!
    rd(*,r-as:r+as) = 0.  &  rd(r-as:r+as,*) = 0.   ;spider !!!!!!!!!!!!!
    support(*,r-as:r+as)=0 & support(r-as:r+as,*)=0 ;spider !!!!!!!!!!!!!
  endif	                                            ;spider !!!!!!!!!!!!!

  support_plus = support NE 0

; Prepara el calculo de los polinomios
			   
  mn = Zernike_mn(jmax)            ;Calculo de [m,n] para cada j.
  fac=fact(10)       ;Calculando 10 factoriales es suficiente hasta j=50 o mas.
  zer=fltarr(2*r,2*r,jmax+1)       ;Array almacen de los polinomios calculados.
  zer(0,0,0)=support   ;en j=0 se almacena la mascara que define la apertura.

; Comienza el gran loop para generar los polinomios *******************

FOR j=1,jmax DO BEGIN

  IF j EQ 1 THEN BEGIN

    z = support_plus

  ENDIF ELSE BEGIN

    m = mn(0,j-1)
    n = mn(1,j-1)

; Calculo de la parte angular del polinomio

;    print,''
    IF m EQ 0 THEN BEGIN
      c = sqrt(double(n+1))
;      print,j,m,n,n+1, $
;		     format='("j,m,n=",I2,",",I2,",",I2,"  sqrt(",F4.1,")*",$)'
    ENDIF ELSE BEGIN
      c = sqrt(double(2*(n+1)))
      IF odd(j) THEN BEGIN 
        c = c * sin(m*th)
;        print,j,m,n,2*(n+1),m, $
;	format='("j,m,n=",I2,",",I2,",",I2,"  sqrt(",F4.1,")*sin(",I1,"Th)*",$)'
      ENDIF ELSE BEGIN
        c = c * cos(m*th)
;        print,j,m,n,2*(n+1),m, $
;	format='("j,m,n=",I2,",",I2,",",I2,"  sqrt(",F4.1,")*cos(",I1,"Th)*",$)'
      ENDELSE
    ENDELSE

; Calculo de la parte radial del polinomio

    z = dblarr(2*r,2*r)
    FOR s=0,(n-m)/2 DO BEGIN
      IF n-2*s EQ 0 THEN BEGIN 
	rx = double(support_plus)
      ENDIF ELSE BEGIN 
	rx = double(rd)^(n-2*s)
      ENDELSE 
      tmp = ((-1.d0)^s*fac(n-s)) / (fac(s)*fac((n+m)/2-s)*fac((n-m)/2-s))
;;     print,max([n-s,s,(n+m)/2-s,(n-m)/2-s])
;      print,tmp,n-2*s,format='(F6.1,"r^",I1,$)'
      z = z + rx * tmp
    ENDFOR
;    print,''

    z = z * c
    
  ENDELSE

  
; Ortogonalization de Gram-Schmidt con los Zernikes ya creados (DARIO)

  FOR I=1,J-1 do begin
	Projec= total(support*Z*Zer(*,*,I) )/total(support)
	Z=Z-Projec*Zer(*,*,I)
  ENDFOR
  

; Reajuste de la rms(Z) a 1. (Ver apendice H)

;;;  zer(0,0,j) = float( z )  
  zer(0,0,j) = float( z / Sqrt(Total(support * z * z)/Total(support)))
; print,sqrt(total(z^2 * support)/total(support))      ;rms de Z_j sin reajuste
; print,sqrt(total(double(zer(*,*,j))^2 * support)/total(support))  ;rms de Z_j
; print,total(double(zer(*,*,j))*support)/total(support)          ;media de Z_j

ENDFOR

; Fin del gran loop *****************************************************

  Return,zer

end
 
