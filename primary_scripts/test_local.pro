function test_local,pol,xta,iii,subfr

; if there is locally in the polarization images any
; residual correlation with I or motion induced xtalk
; it removes it. It helps Q and U a lot.

sz1=sizeof(pol,1)
sz2=sizeof(pol,2)

; zero out signal in pol

dd=where(abs(pol) gt 3.*stdev(pol),nns)
pols=pol
if (nns gt 0) then pol[dd]=0.

num_local=sz1/subfr

cte=fltarr(num_local,num_local)
co1=fltarr(num_local,num_local)
co2=fltarr(num_local,num_local)
szlocal=long(subfr)*long(subfr)

for ind=0,sz1-1,subfr do begin
	for jnd=0,sz2-1,subfr do begin
		pola=pol[ind:ind+subfr-1,jnd:jnd+subfr-1]
		iiia=iii[ind:ind+subfr-1,jnd:jnd+subfr-1]
		xtaa=xta[ind:ind+subfr-1,jnd:jnd+subfr-1]
		xa=transpose([[reform(xtaa,szlocal)],[reform(iiia,szlocal)]])
		res=regress(xa,reform(pola,szlocal),const=ccc)
		dd=(finite(res,/NAN))
		if (dd[0] eq 1) then begin
			res[*]=0
			ccc=0
		endif
		if (dd[1] eq 1) then begin
			res[*]=0
			ccc=0
		endif
		cte(ind/subfr,jnd/subfr)=ccc
		co1(ind/subfr,jnd/subfr)=res[0]
		co2(ind/subfr,jnd/subfr)=res[1]
	endfor
endfor

ctee=smooth(rebin(cte,sz1,sz2),10)
co1e=smooth(rebin(co1,sz1,sz2),10)
co2e=smooth(rebin(co2,sz1,sz2),10)

return,pols-(ctee+co1e*xta+co2e*iii)

end


