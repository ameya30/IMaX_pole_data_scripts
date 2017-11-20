pro imax_find_sigma_run_dial
	
	; Read in input image	
	
	files = file_search('../imax_lp_max_/imax_lp_max_08*')
	
	for i = 0, n_elements(files)-1 do begin
	        
		im = readfits(files[i])

		seg_out = make_array(size(im, /dim))	

		; Create output name

		split_name = strsplit(files[i], '_', /extract)
		out_name   = '../imax_max_LP_seg_/imax_max_LP_seg_' + split_name[6] 	
		print, out_name

		for wv = 0, 4 do begin

			; Find noise level
		
			noise_lev = stddev(im[100:699, 240:639, wv])
		
			; Run dialation routine to find segmented image with centre 3 sigma and outer 2 sigma

			mlt, abs(im), seg, levels = [3*noise_lev, 2*noise_lev]
			
			seg_out[*, *, wv] = seg

		endfor
	
		writefits, out_name, seg_out
	endfor
end
