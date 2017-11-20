;+
; NAME:
;       NINT
; PURPOSE:
;       returns nearest integer of real-value x
;*CATEGORY:            @CAT-# 19@
;       Mathematical Functions (General)
; CALLING SEQUENCE:
;       n = NINT(x [,/LONG])
; INPUTS:
;       x = real value or vector; values must be -32767 <= x <= 32767
;                   (set /LONG or use NLONG (KIS_LIB) or ROUND (IDL-
;                    Routine) if |x| larger).
; OPTIONAL KEYWORD INPUT:
;       LONG - If this keyword is set and non-zero, then the result of NINT
;               is of type LONG.   Otherwise, the result is of type INTEGER.
;               (Nearest long integer is also returned from NLONG (KIS_LIB)
;               or from ROUND (IDL- Routine Vers. 3.1) 
; OUTPUTS:
;       next integer (single value or integer-array);
; PROCEDURE:
;       
; MODIFICATION HISTORY:
;       nlte, 1990-01-05
;       nlte, 1992-05-06 : faster code, no indices 
;       nlte, 1993-08-04 : keyword /LONG
;-
FUNCTION NINT,x,long=long
on_error,1
if n_elements(x) lt 1 then message,'argument undefined'
;
if keyword_set(long) then return,long(x+float(x ge 0)-0.5) else $
                          return,fix(x+float(x ge 0)-0.5)
end

