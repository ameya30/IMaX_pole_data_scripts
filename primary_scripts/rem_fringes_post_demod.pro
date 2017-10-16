;**********************************************************************

@lib_wiener
@mean_power2
@masking_power
@filtering_fringes
@rubbish

;*********************************************************************************

pro rem_fringes_post_demod, stok = stok, wvln = wvln;, save_out = save_out;file_in
files = file_search('saves_May31st/mk_mag*.sav')

for i = 0, n_elements(files)-1 do begin	
        print, files[i]
	; Read file
	restore, files[i]
	full_cycle = iid
	
	split_name = strsplit(files[i], '_', /extract)
	split_2    = strsplit(split_name[6], '.', /extract)

        save_out = 'saves_May31st/post_demod_output_' + split_2[0]
	
	print, save_out
	; Take single image
        if i EQ 0 then begin	
		ima = full_cycle[*, *, wvln, stok]
		
		xdim = (size(ima))[1]
		ydim = (size(ima))[2]
		
		; Create apopodisation window
		
		apo_win = win_cos(xdim, ydim, 10)
		
		; Compute fft of image
		
		fft_ima = shift(fft(ima * apo_win, -1), xdim/2, ydim/2)
		
		; Compute power spectrum
		
		pow_ima = (abs(fft_ima)) ^ 2
		
		tvscl, pow_ima
		
		; Based on the power spectrum, construct the mask defining the areas including
		; spurious signal (fringes). Also construct the weighting function to filter out this 
		; spurious signal.
		
		print, 'RUNNING PROGRAM masking_power.pro'
		print, 'INTERACTIVE PROCESS WITH THE MOUSE: "rubber tool"'
		masking_power, pow_ima, pow_corr, mask, weight
	endif
	; Remove fringes from the image
	
	print, 'RUNNING PROGRAM filtering_fringes.pro'
	
	for wvln = 0, 4 do begin
		for stok = 0, 3 do begin
			
			ima = full_cycle[*, *, wvln, stok]
			
			filtering_fringes, 0, 0, 5, weight, ima, ima_no_fringes
			
			full_cycle[*, *, wvln, stok] = ima_no_fringes
	
		endfor
	endfor	
	
	save, weight, ima, ima_no_fringes, pow_ima, pow_corr, mask, filename = save_out + '_params.sav'
	
	writefits, save_out + '.fits', full_cycle

endfor
end
