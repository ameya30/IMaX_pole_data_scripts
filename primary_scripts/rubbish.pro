pro rubbish,imm,smoow,grains,maskg,maskgd

; Detecta el patron de motas dobles y otra basurilla en una imagen
; defectuosamente flatfielded. Utiliza la tecnica de "unsharp-masking"
; y "thresholding". Los "lower and upper thresholds" los solicita el 
; programa tras una inspeccion de la imagen "grains" con contraste resaltado.

;INPUT:
;   imm = imagen sobre la que se va a detectar el patron de estruct.espureas.
;         Se recomienda que sea una imagen promedio de varias con el mismo patron.
;   smoow = smoothing window width for procedure "smooth.pro". Odd number whose
;           value has to be approx.the width of the structures we wish to detect.
;           For the case of IMaX, 5 pix is a reasonable value.
;
;OUTPUTS:
;   grains = Imagen con motas y basurilla resaldas sobre un fondo muy atenuado.
;   maskg  = Mascara binaria que define con 1s los puntos afectados de basurilla.
;            Se obtiene mediante "thresholding" en la imagen "grains".
;   maskgd = Dilated mask by an opening operation

;   INPUTS TO BE SELECTED BY EDITOR !!!!!!!!!!!!!
;   Available shape-operators (kernels) in the "closing" and "opening" operations,
;   respectively.
;
;        (1) 1 1  (2) 0 1 0  (3) 1 1 1  (4) 0 1 1 0  (5) 0 0 1 0 0  (6) 1 1 1 1
;            1 1      1 1 1      1 1 1      1 1 1 1      0 1 1 1 0      1 1 1 1
;                     0 1 0      1 1 1      1 1 1 1      1 1 1 1 1      1 1 1 1
;                                           0 1 1 0      0 1 1 1 0      1 1 1 1
;                                                        0 0 1 0 0
;
;        (7) 0 1 1 1 0  (8) 1 1 1 1 1
;            1 1 1 1 1      1 1 1 1 1
;            1 1 1 1 1      1 1 1 1 1
;            1 1 1 1 1      1 1 1 1 1
;            0 1 1 1 0      1 1 1 1 1

;For programming: Closing and opening operators for eroding and dilating structures
  ;kernels for opening and closing operations (shape operator)
  ;1: ker_o/c=replicate(1,2,2)
  ;2: ker_o/c=[[0,1,0],[1,1,1],[0,1,0]]
  ;3: ker_o/c=replicate(1,3,3)
  ;4: ker_o/c=[[0,1,1,0],[1,1,1,1],[1,1,1,1],[0,1,1,0]]
  ;5: ker_o/c=[[0,0,1,0,0],[0,1,1,1,0],[1,1,1,1,1],[0,1,1,1,0],[0,0,1,0,0]]
  ;6: ker_o/c=replicate(1,4,4)
  ;7: ker_o/c=[[0,1,1,1,0],[1,1,1,1,1],[1,1,1,1,1],[1,1,1,1,1],[0,1,1,1,0]]
  ;8: ker_o/c=replicate(1,5,5)
;__________________________________________________________________
; Edited 23-Sep-2009
; Jose A. Bonet
;________________________________________________________________

;Unsharp masking

grains=imm-smooth(imm,smoow)
grains=grains-smooth(grains,smoow)
;grains=grains-smooth(grains,3)
;tvwinp,grains>(-0.1)<0.1,title='DETECCION DE LOWER AND UPPER THRSHOLDS'
tvwin,grains>(-0.1)<0.1,title='DETECCION DE LOWER AND UPPER THRSHOLDS'
print,'The higher the thresholds the smaller the number of dust particles detected'
read,'lower and upper thresholds (e.g. -0.03,0.03)= ',lo,hi
maskg=grains*0.
maskg=grains ge hi or grains le lo

;Closing and opening operations
ker_c=replicate(1,2,2)   ;kernel for closing
ker_o=[[0,1,0],[1,1,1],[0,1,0]]  ;kernel for dilating
;maskg = dilate(maskg,ker_o)
;maskg = dilate(erode(maskg,ker_c),ker_o)
maskg = erode(dilate(maskg,ker_o),ker_c)
;maskg=smooth(float(maskg),3) gt 0.49
ker_o=[[0,1,1,1,0],[1,1,1,1,1],[1,1,1,1,1],[1,1,1,1,1],[0,1,1,1,0]]
maskgd=dilate(maskg,ker_o)

end
