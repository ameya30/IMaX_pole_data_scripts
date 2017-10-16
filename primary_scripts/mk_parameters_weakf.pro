pro mk_parameters_weakf,obsen,cyclen

;IDL>mk_parameters_weakf,163,346 

; Provides absolute velocity, line FWHM, line depth,
; continuum, averaged linear polarization, averaged circular 
; polarization and calibrated longitudinal and transverse 
; magnetograms. Continuum circular polarization included

; Variables ending with 'nr' are level 1 (not resotred)
; Variables starting with 's' are errors !

; Blueshift is fully taken care off
; Velocities are in an absolute scale

; restore, housekeepings

restore,'set09_hdr.save',/v

; restore data file

pathsi='~/IDLWorkspace71/kiruna_flight/DECON_MPS/new_jab/reduc9r/'
pathso='./'

infile=pathsi+'reduc_rnr_'+strcompress(string(fix(obsen)),/remove_all)+ $
'_'+strcompress(string(fix(cyclen)),/remove_all)+'.save'

restore,infile,/v

ndim=sizeof(iid,1)
n_lambda=sizeof(iid,3)
n_pol=sizeof(iid,4)

; restore blueshift file

restore,'blueshiftcal_remlin_09_1.save',/v

; recover blueshift 

xr = (findgen(ndim)) # replicate(1,ndim) 
yr = transpose(xr) 
lm=sk(0)+sk(1)*yr+sk(2)*yr^2+sk(3)*yr^3+                        $
sk(4)*xr+sk(5)*xr*yr+sk(6)*xr*yr^2+sk(7)*xr*yr^3+                $
sk(8)*xr^2+sk(9)*xr^2*yr+sk(10)*xr^2*yr^2+sk(11)*xr^2*yr^3+      $
sk(12)*xr^3+sk(13)*xr^3*yr+sk(14)*xr^3*yr^2+sk(15)*xr^3*yr^3     

; lm is blueshift array

lll=[-.08d0,-.04d0,.04d0,.08d0,.227d0]

llo=5250.225d0

lla=fltarr(ndim,ndim,n_lambda) ; wavelengths for all FOV pixels

for ind=0,n_lambda-1 do lla[*,*,ind]=-lm[*,*]+lll[ind]

; reconstructed data

lmin=fltarr(ndim,ndim) ; wavelength of each pixel (i.e. IMaX wavelength calibration) 
fwhm=fltarr(ndim,ndim) ; fwhm for each pixel
linc=fltarr(ndim,ndim) ; line center intensity for each pixel
slmin=fltarr(ndim,ndim) ; sigma of wavelength of each pixel (i.e. IMaX wavelength calibration) 
sfwhm=fltarr(ndim,ndim) ; sigma of fwhm for each pixel
slinc=fltarr(ndim,ndim) ; sigma of line center intensity for each pixel

; nr is non reconstructed

lminnr=fltarr(ndim,ndim) ; wavelength of each pixel (i.e. IMaX wavelength calibration) 
fwhmnr=fltarr(ndim,ndim) ; fwhm for each pixel
lincnr=fltarr(ndim,ndim) ; line center intensity for each pixel
slminnr=fltarr(ndim,ndim) ; sigma of wavelength of each pixel (i.e. IMaX wavelength calibration) 
sfwhmnr=fltarr(ndim,ndim) ; sigma of fwhm for each pixel
slincnr=fltarr(ndim,ndim) ; sigma of line center intensity for each pixel

; arrays for problematic areas

probl=fltarr(ndim,ndim)
problnr=fltarr(ndim,ndim)

; Stokes I errors are from S/N 800 (continuum)
; Stokes I errors are from S/N 800/3 for restored data (continuum)

yerrnr=fltarr(n_lambda-1)
yerr=fltarr(n_lambda-1)
for ind=0,n_lambda-2 do yerrnr[ind]=mean(iidn[*,*,ind,0])/800.
for ind=0,n_lambda-2 do yerr[ind]=mean(iid[*,*,ind,0])/267.

saa=fltarr(3)

; lambda array for weak field

llwf=(findgen(200)-99)/300.
blon_wf=fltarr(ndim,ndim)
blonnr_wf=fltarr(ndim,ndim)

; fitting loop

for ind=0,ndim-1 do begin
	print,ind	
	for jnd=0,ndim-1 do begin

;;;;;;;;;non-restored
		fint=reform(iidn[ind,jnd,*,0])
		syerr=yerrnr/fint[n_lambda-1]
		fint=1.-fint/fint[n_lambda-1]
		yyy=fint[0:n_lambda-2]
		dd=where(yyy lt 0.)
		if dd(0) ne -1 then begin 
;			yyy[dd]=0.01
			problnr[ind,jnd]=mean(yyy[dd])	
		endif
		co=gaussfit_status(lla[ind,jnd,0:n_lambda-2],yyy,aa,measure_errors=syerr, $
		nterms=3,sigma=saa,STATUS=st_fit)
		if (st_fit ne 0) then problnr[ind,jnd]=1.	
		lminnr[ind,jnd]=aa[1]
		fwhmnr[ind,jnd]=aa[2]*2.*SQRT(2.*ALOG(2.))
		lincnr[ind,jnd]=1.-aa[0]
		slminnr[ind,jnd]=saa[1]
		sfwhmnr[ind,jnd]=saa[2]*2.*SQRT(2.*ALOG(2.))
		slincnr[ind,jnd]=saa[0]	
;weak field non restored
                yarr=iidn[ind,jnd,n_lambda-1,0]*(1.-(aa[0])*exp(-(llwf-aa[1])^2/(2.*aa[2]^2)))
                dyarr=deriv(llwf,yarr)
                idyarr=interpol(dyarr,llwf,lla[ind,jnd,*])
                co=poly_fit(reform(idyarr),-reform(iidn[ind,jnd,*,3]),1)
                blonnr_wf[ind,jnd]=co(1)/(4.67d-13*3.*5250.^2) ; field in G


;;;;;;;;;restored
		fint=reform(iid[ind,jnd,*,0])
		syerr=yerr/fint[n_lambda-1]
		fint=1.-fint/fint[n_lambda-1]
		yyy=fint[0:n_lambda-2]
		dd=where(yyy lt 0.)
		if dd(0) ne -1 then begin 
;			yyy[dd]=0.01
			probl[ind,jnd]=mean(yyy[dd])	
		endif
		co=gaussfit_status(lla[ind,jnd,0:n_lambda-2],yyy,aa,measure_errors=syerr, $
		nterms=3,sigma=saa,STATUS=st_fit)
		if (st_fit ne 0) then probl[ind,jnd]=1.	
		lmin[ind,jnd]=aa[1]
		fwhm[ind,jnd]=aa[2]*2.*SQRT(2.*ALOG(2.))
		linc[ind,jnd]=1.-aa[0]
		slmin[ind,jnd]=saa[1]
		sfwhm[ind,jnd]=saa[2]*2.*SQRT(2.*ALOG(2.))
		slinc[ind,jnd]=saa[0]
;weak field restored
                if (probl[ind,jnd] eq 0) then begin
                        yarr=iid[ind,jnd,n_lambda-1,0]*(1.-(aa[0])*exp(-(llwf-aa[1])^2/(2.*aa[2]^2)))
                        dyarr=deriv(llwf,yarr)
                        idyarr=interpol(dyarr,llwf,lla[ind,jnd,*])
                endif else begin
                        print,'using non-restored derivative'
                endelse
; if bad fit use non restored Stokes I derivative
                co=poly_fit(reform(idyarr),-reform(iid[ind,jnd,*,3]),1)
                blon_wf[ind,jnd]=co(1)/(4.67d-13*3.*5250.^2) ; field in G
	
	endfor
endfor

; gaussian fit is not 100 % robuts an creates spikes
; filtering them with a median filter

lmin=median(lmin,3) 
fwhm=median(fwhm,3)
linc=median(linc,3)

lminnr=median(lminnr,3) 
fwhmnr=median(fwhmnr,3)
lincnr=median(lincnr,3)

lmin=lmin-mean(lmin)
lminnr=lminnr-mean(lminnr)

cspeed=299792.458d0 ; c km/s

; we substract some blueshift

lmin=lmin*cspeed/llo-0.2
lminnr=lminnr*cspeed/llo-0.2  ; -0.2 is granulation blueshift (very arbitrary, but not bad)
slmin=slmin*cspeed/llo
slminnr=slminnr*cspeed/llo  ; 

; note: with this correction velocities are in an absolute scale !!!!!

window,0,xsize=900,ysize=800
tvframe,lmin<3>(-3),/bar,/aspect,xrange=[0,936*.055],yrange=[0,936*.055],$          
title='IMaX/Sunrise 09/06 01:54:12 Velocity (km/s)',charsize=1.5  
write_png,'velocity_cal.png',tvrd(/true)

window,1,xsize=900,ysize=800
tvframe,(fwhm*1e3)<400>(0),/bar,/aspect,xrange=[0,936*.055],yrange=[0,936*.055],$          
title='IMaX/Sunrise 09/06 01:54:12 FWHM (mA)',charsize=1.5  
write_png,'fwhm_cal.png',tvrd(/true)

window,2,xsize=900,ysize=800
tvframe,(linc)<1.>(.4),/bar,/aspect,xrange=[0,936*.055],yrange=[0,936*.055],$          
title='IMaX/Sunrise 09/06 01:54:12 Line depth',charsize=1.5  
write_png,'line_depth_cal.png',tvrd(/true)

window,0,xsize=900,ysize=800
tvframe,lminnr<3>(-3),/bar,/aspect,xrange=[0,936*.055],yrange=[0,936*.055],$          
title='IMaX/Sunrise 09/06 01:54:12 Velocity (km/s)',charsize=1.5  
write_png,'velocity_cal_nr.png',tvrd(/true)

window,1,xsize=900,ysize=800
tvframe,(fwhmnr*1e3)<200>(0),/bar,/aspect,xrange=[0,936*.055],yrange=[0,936*.055],$          
title='IMaX/Sunrise 09/06 01:54:12 FWHM (mA)',charsize=1.5  
write_png,'fwhm_cal_nr.png',tvrd(/true)

window,2,xsize=900,ysize=800
tvframe,(lincnr)<1.>(.4),/bar,/aspect,xrange=[0,936*.055],yrange=[0,936*.055],$          
title='IMaX/Sunrise 09/06 01:54:12 Line depth',charsize=1.5  
write_png,'line_depth_cal_nr.png',tvrd(/true)

;  generate linear polarization arrays averaged over the line

linpol=0.25*(sqrt(iid[*,*,0,1]^2+iid[*,*,0,2]^2)+ $
sqrt(iid[*,*,1,1]^2+iid[*,*,1,2]^2)+ $
sqrt(iid[*,*,2,1]^2+iid[*,*,2,2]^2)+ $
sqrt(iid[*,*,3,1]^2+iid[*,*,3,2]^2)) $
/(iid[*,*,4,0])

linpolnr=0.25*(sqrt(iidn[*,*,0,1]^2+iidn[*,*,0,2]^2)+ $
sqrt(iidn[*,*,1,1]^2+iidn[*,*,1,2]^2)+ $
sqrt(iidn[*,*,2,1]^2+iidn[*,*,2,2]^2)+ $
sqrt(iidn[*,*,3,1]^2+iidn[*,*,3,2]^2)) $
/(iidn[*,*,4,0])

; generate circular polarization arrays averaged over the line

vvv=0.25*(iid[*,*,0,3]+iid[*,*,1,3]       $
          -iid[*,*,2,3]-iid[*,*,3,3])     $
	  /(iid[*,*,4,0])

vvvnr=0.25*(iidn[*,*,0,3]+iidn[*,*,1,3]       $
          -iidn[*,*,2,3]-iidn[*,*,3,3])     $
	  /(iidn[*,*,4,0])

; generate continuum frame

cct=iid[*,*,4,0]/mean(iid[*,*,4,0])
cctnr=iidn[*,*,4,0]/mean(iidn[*,*,4,0])

; generate Stokes V continuum frame

cctv=iid[*,*,4,3]/iid[*,*,4,0]
cctvnr=iidn[*,*,4,3]/iidn[*,*,4,0]

; generate longitudinal and transverse magnetograms with
; correction for blueshift !!!

lm=lm*1.e3 ; in mA

; note glon_ref and gtra_ref are arbiitrary values. If longitudinal
; fields are always at 0 degrees and transverse fields always at
; 90 degrees they are right. For the QS signals we analyze, probably
; we should be using 45 degrees for both. To be discussed among competent
; scientists

glon_ref=0.
;glon_ref=45.
clon=4762./cos(glon_ref*!pi/180.)-8.97250*lm+3.03625*lm*lm

gtra_ref=90.
;gtra_ref=45.
ctra=2527.00/sin(gtra_ref*!pi/180.)+1.72167*lm+0.322500*lm*lm 

blon=clon*vvv
btra=ctra*sqrt(linpol)

blonnr=clon*vvvnr
btranr=ctra*sqrt(linpolnr)

ofile=pathso+'parameters_reduc_rnr_wf_'+strcompress(string(fix(obsen)),/remove_all)+ $
'_'+strcompress(string(fix(cyclen)),/remove_all)+'.save'

save, $
lmin,fwhm,linc,slmin,sfwhm,slinc,cct,cctv,linpol,vvv,blon,btra,probl, $
lminnr,fwhmnr,lincnr,slminnr,sfwhmnr,slincnr, $
cctnr,cctvnr,linpolnr,vvvnr,blonnr,btranr,problnr,blon_wf,blonnr_wf,  $
filename=ofile


end
