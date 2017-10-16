@lib_wiener
@mean_power2
@masking_power
@filtering_fringes
@rubbish

pro rem_med_flat, file_in, save_out

all = readfits(file_in)

size_all = size(all, /dim)

xdim = size_all[0]
ydim = size_all[1]

out_flat = make_array(size_all)

for stok = 0, 3 do begin
	for wvln = 0, 4 do begin
	     	
		ima = all[*, *, wvln, stok]

		; Create apopodisation window

		apo_win = win_cos(xdim, ydim, 10)
	
		; Compute fft of image
	
		fft_ima = shift(fft(ima * apo_win, -1), xdim/2, ydim/2)
	
		; Compute power spectrum
	
		pow_ima = (abs(fft_ima)) ^ 2
	
		; filter fringes
	
		masking_power, pow_ima, pow_corr, mask, weight
	
		filtering_fringes, 0, 0, 5, weight, ima, flat_no_med			
	 		
		out_flat[*, *, wvln, stok] = flat_no_med	
		
		undefine, pow_corr
		undefine, mask	
                undefine, weight
	endfor
endfor

writefits, save_out, out_flat 

end
