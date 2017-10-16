pro imax_hdr_str,hdr,hdrstr

; input (hdr) is FITS header string format
; output (hdrstr) if FITS header structure format, all in numeric values

; structure is : datasize1, datasize2, time (hours from 00:00 8th july),
; time at which observing run started (last reset), time accumulations started
; time accumulations ended, IMaX observing Number, IMaX cycle ,
; xcen,ycen,f2_mechanism, cw_loop, ps_state, pd_plate, n_wavel,
; n_accum, npola, d_wavel, polariz. state, etalon temp, lcvr temp.,
; doubleits temp., a coeficient (wavel cal during flight), SSC

aa=fxpar(hdr, 'DATE_OBS') ; time of the observations
day=double(strmid(aa,8,2))
hour=double(strmid(aa,11,2))
min=double(strmid(aa,14,2))
sec=double(strmid(aa,17,2))
msec=double(strmid(aa,20,3))
timeobs=(day-8.d0)*24.+hour+min/60.d0+(sec+msec*1.d-3)/3600.d0

aa=fxpar(hdr, 'OBS_TIME') ; time at which observing run started
if (aa ne '           0') then begin ; sometimes this parameter is not set
	day=double(strmid(aa,8,2))
	hour=double(strmid(aa,11,2))
	min=double(strmid(aa,14,2))
	sec=double(strmid(aa,17,2))
	msec=double(strmid(aa,20,3))
	timerun=(day-8.d0)*24.+hour+min/60.d0+(sec+msec*1.d-3)/3600.d0
	obstime=0
endif else begin
	obstime=1
endelse

aa=fxpar( hdr, 'IMG_STIM') ; time at which accumulations started 
day=double(strmid(aa,8,2))
hour=double(strmid(aa,11,2))
min=double(strmid(aa,14,2))
sec=double(strmid(aa,17,2))
msec=double(strmid(aa,20,3))
imgacst=(day-8.d0)*24.+hour+min/60.d0+(sec+msec*1.d-3)/3600.d0

if obstime eq 1 then timerun=imgacst

aa=fxpar( hdr, 'IMG_ETIM') ; time at which accumulations ended
day=double(strmid(aa,8,2))
hour=double(strmid(aa,11,2))
min=double(strmid(aa,14,2))
sec=double(strmid(aa,17,2))
msec=double(strmid(aa,20,3))
imgacen=(day-8.d0)*24.+hour+min/60.d0+(sec+msec*1.d-3)/3600.d0

obsmod=strmid(fxpar( hdr, 'OBSMODE' ),0,4)
case obsmod of
	'V': obsmode=1
	'L': obsmode=2
	'D': obsmode=3
	'M': obsmode=4
else: obsmode=0
endcase

f2m=strmid(fxpar( hdr, 'F2_MECH' ),0,4)
case f2m of
	'fiel': f2state=1
	'pinh': f2state=2
	'dark': f2state=3
else: f2state=0
endcase

dummy=strmid(fxpar( hdr, 'AP_DOOR' ),0,4)
if (dummy eq 'clos') then f2state=3
dummy=strmid(fxpar( hdr, 'FLAT_MOD' ),0,4)
if (dummy eq 'on') then f2state=0



cwl=fxpar( hdr, 'CW_LOOP' )
case cwl of
	'open': cwlstate=0
else: cwlstate=1
endcase

pss=strmid(fxpar( hdr, 'PS_STATE' ),0,9)
case pss of
	'15 arcsec': psstate=1
	'15 arcmin': psstate=2
	'flatfield': psstate=3
else: psstate=0
endcase
if (dummy eq 'on') then psstate=3

pdpl=fxpar( hdr, 'PD_PLATE' )
case pdpl of
	'on': pdstate=1
	'off': pdstate=2
else: pdstate=0
endcase

numwave=fxpar( hdr, 'NUM_WAVE' )
naccu=fxpar( hdr, 'NUM_ACCU' )
npola=fxpar( hdr, 'NUM_POLA' )
dwave=fxpar(hdr, 'D_WAVE' )
nwave=fxpar(hdr, 'N_WAVE' )

polstat=fxpar( hdr, 'POL_STAT' )
case polstat of
	'I1': polstat=1
	'I2': polstat=2
	'I3': polstat=3
	'I4': polstat=4
else: polstat=0
endcase

lamstat=fxpar( hdr, 'N_WAVE' )
etaltem=fxpar(hdr, 'ETAL_TEM' )
rocltem=fxpar(hdr, 'ROCL_TEM' )
etadtem=fxpar(hdr, 'ETAD_TEM' )
etalaco=fxpar(hdr, 'ETAL_ACO' )
sscpar=fxpar(hdr, 'SSC' )

; create structure

hdrstr={                                                                                        $
sz1:fxpar( hdr, 'NAXIS1' ),sz2:fxpar( hdr, 'NAXIS2' ),                                          $
time:timeobs,time_run:timerun,img_acst:imgacst,img_acen:imgacen,                                $
num_obse:fxpar( hdr , 'NUM_OBSE' ),num_cycl:fxpar( hdr, 'NUM_CYCL' ),xcen:fxpar( hdr, 'XCEN' ), $
ycen:fxpar( hdr, 'YCEN' ),cw_loop:cwlstate,obs_mode:obsmode,f2_mech:f2state,ps_state:psstate,   $
pd_plate:pdstate,num_wave:numwave,num_accu:naccu,num_pola:npola,d_wave:dwave,n_wave:nwave,      $
pol_stat:polstat,lam_stat:lamstat,etal_tem:etaltem,rocl_tem:rocltem,etad_tem:etadtem,           $
etal_aco:etalaco,ssc:sscpar                                                                     $       
}

end
