pro check_align

imax_file = 'fringe_norm_May8th_07.fits'
sufi_file = 'sufi_resamp_397_161239.fits'

imax = readfits(imax_file)
sufi = readfits(sufi_file)

sufi_dim = size(sufi, /dim)
print, size(imax, /dim)

x = 0

while x LE 20 do begin
	tvscl, imax[310:598, 114:810, 1, 0]
	wait, 0.2
	tvscl, sufi 
	wait, 0.2
	x += 1
endwhile


end
