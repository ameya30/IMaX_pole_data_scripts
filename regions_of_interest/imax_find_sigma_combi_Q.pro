pro imax_find_sigma_combi_Q,filename,noise_lev=noise_lev

  ; Read in input image

  

    im = readfits(filename)

    seg_out = make_array(size(im, /dim))

    ; Create output name

    split_name = strsplit(filename, '_', /extract)
    s = size(split_name)
    if s[1] eq 10 then begin
      out_name   = '/scratch/prabhu/backup_workstation/sunrise_holly/imax_maxlp_combi_wv_seg_/imax_maxlp_combi_wv_seg_' + split_name[9]
      print, out_name
    endif else begin
      out_name   = '/scratch/prabhu/backup_workstation/sunrise_holly/imax_maxlp_combi_wv_seg_/imax_maxlp_combi_wv_seg_' + split_name[9]+'_'+split_name[10]
      print, out_name
    endelse
     
    print,noise_lev

    ; Run dialation routine to find segmented image with centre 3 sigma and outer 2 sigma
    
    mlt, abs(im), seg, levels = [3*noise_lev, 2*noise_lev]

    seg_out = seg

    

    writefits, out_name, seg_out
    
end
