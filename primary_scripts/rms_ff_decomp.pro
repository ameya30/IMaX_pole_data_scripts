;; decomposition of flat-fied image into a low frequency smooth field
;; (ff_low), interference fringes (ff_fringes), and a high-frequency
;; component (ff_high) including dust, etc.. 
;; opt is the degree of the polinomial to fit the ff_low component
;; smooth: how many pixels to smooth for fringes fitting

pro rms_ff_decomp, file_in

save_name = 'cam2_flat_full_no_med2.fits'

deg   = 4
smuth = 10

full_cycle = readfits(file_in)

sz = size(full_cycle, /dim)

nx = sz[0]
ny = sz[1]

no_fringe_full = make_array(sz)

for wvln = 0, 4 do begin
	for stok = 0, 3 do begin
			
		ff = full_cycle[*, *, wvln, stok]

		ff_lo = fltarr(nx, ny)
		ff_fr = fltarr(nx, ny)
		ff_hi = fltarr(nx, ny)
		
		;; ff_lo:
		
		aux = fltarr(nx, ny)
		
		coeffs = fltarr(deg+1, ny)
		scoeffs = coeffs
		
		for i = 0, ny-1 do coeffs(*, i) = poly_fit(findgen(nx), ff[*, i], deg, yfit=yfit)
		   
		for i = 0, deg do begin
		   res = poly_fit(findgen(ny), coeffs[i, *], 4, yfit=yfit)
		   scoeffs[i, *] = yfit
		endfor
		
		for j = 0, ny-1 do ff_lo[*, j] = scoeffs[deg, j] 
		for i = deg-1, 0, -1 do for j = 0, ny-1 do ff_lo[*, j] = ff_lo[*, j] * findgen(nx) + scoeffs[i, j]
		
		for j = 0, ny-1 do aux[*, j] = coeffs[deg, j] 
		for i = deg-1, 0, -1 do for j = 0, ny-1 do aux[*, j] = aux[*, j] * findgen(nx) + coeffs[i, j]
		
		;; ff_fringes:
		
		dff = ff - ff_lo
		
		for i = 0, ny-1 do ff_fr[*, i] = smooth(smooth(dff[*, i], smuth, /edge_mirror), smuth, /edge_mirror)
		
		;; ff_hi
		
		ff_hi = dff - ff_fr
		
		no_fringe_full[*, *, wvln, stok] = ff_hi
	endfor
endfor

writefits, save_name, no_fringe_full

end
