pro IMAX_SUFI397_ALIGN

; Get data and info
;imax_list = ['fringe_norm_May8th_07.fits','fringe_norm_May8th_12.fits', 'fringe_norm_May8th_21.fits', 'fringe_norm_May8th_24.fits']
;sufi_list = ['sufi_resamp_300_161240.fits', 'sufi_resamp_300_161526.fits', 'sufi_resamp_300_162022.fits', 'sufi_resamp_300_162206.fits']
imax_list = 'fringe_norm_May8th_07.fits'
sufi_list = 'sufi_resamp_397_161239.fits'


for i = 0,0 do begin
	imax = readfits(imax_list[i])
	sufi = readfits(sufi_list[i])
	full_sufi_dim = size(sufi, /dim)
	imax_sub = imax[300+full_sufi_dim[0]/2:299+full_sufi_dim[0], 80:80+full_sufi_dim[0]/2, 1, 0]
	sufi_sub = sufi[full_sufi_dim[0]/2:full_sufi_dim[0]-1, 0:full_sufi_dim[0]/2]
	imax_dim = size(imax_sub, /dim)
	sufi_dim = size(sufi_sub, /dim)
	imax_sub = imax_sub / max(imax_sub)
	sufi_sub = sufi_sub / max(sufi_sub)
	
	off = shc(imax_sub, sufi_sub)
	x = off[0]
	y = off[1]
	
	print, x, y
	
	;both = imax_sub
	;;both[0:sufi_dim(0)-17, 6:sufi_dim(1)-1] = sufi_sub[16:sufi_dim(0)-1, 0:sufi_dim(1)-7]
	test  = imax[*, *, 1, 0]
	test2 = imax[*, *, 1, 0]
	test[300+full_sufi_dim[0]/2+x:300+full_sufi_dim[0]/2+sufi_dim[0]+x-1, 80+y:80+y+sufi_dim[0]-1] = sufi_sub
	test2[300+x : 299+x+full_sufi_dim[0], 80+y : 79+y+full_sufi_dim[1]] = sufi
	
	print, 'imax aligned to sufi: (', 300+x, ':', 299+x+full_sufi_dim[0], ', ', 80+y, ': ',  79+y+full_sufi_dim[1], ')'
	
	x = 0
	for x = 0, 15 do begin
		tvscl,imax[*,*,1,0] 
		wait,0.2
		tvscl,test2
		wait, 0.2
		x = x+1 
	endfor
	tvscl, test
endfor
end
	

