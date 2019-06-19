PRO FTS_tool_convolve_fits_image_idl, fits_file, from_beam = from_beam, to_beam = to_beam, output_file = output_file
    
    !EXCEPT=2 ; check math error
    
    ; 
    fits_image = mrdfits(fits_file, 0, fits_header)
    fits_image_nonan = fits_image
    fits_image_nonan[where(~finite(fits_image),/null)] = 0.0
    
    ; 
    fits_key_naxis1 = fxpar(fits_header, 'NAXIS1')
    if fits_key_naxis1 ne '' then fits_key_naxis1 = fix(fits_key_naxis1) else fits_key_naxis1 = float('nan') ; pixel
    fits_key_naxis2 = fxpar(fits_header, 'NAXIS2')
    if fits_key_naxis2 ne '' then fits_key_naxis2 = fix(fits_key_naxis2) else fits_key_naxis2 = float('nan') ; pixel
    ; 
    fits_key_bmaj = fxpar(fits_header, 'BMAJ')
    if fits_key_bmaj ne '' then fits_key_bmaj = double(fits_key_bmaj) * 3600.0 else fits_key_bmaj = float('nan') ; arcsec
    fits_key_bmin = fxpar(fits_header, 'BMIN')
    if fits_key_bmin ne '' then fits_key_bmin = double(fits_key_bmin) * 3600.0 else fits_key_bmin = float('nan') ; arcsec
    fits_key_bpa = fxpar(fits_header, 'BPA')
    if fits_key_bpa ne '' then fits_key_bpa = double(fits_key_bpa) else fits_key_bpa = float('nan') ; degree
    
    ; 
    if n_elements(from_beam) eq 1 then begin
        if from_beam le 0 then begin
            if finite(fits_key_bmaj) and finite(fits_key_bmin) then begin
                if finite(fits_key_bpa) then begin
                    from_beam = [fits_key_bmaj, fits_key_bmin, fits_key_bpa]
                endif else begin
                    from_beam = [fits_key_bmaj, fits_key_bmin, 0.0]
                endelse
            endif else if finite(fits_key_bmaj) then begin
                from_beam = [fits_key_bmaj, fits_key_bmaj, 0.0]
            endif else begin
                message, 'Error! No BMAJ BMIN BPA keyword found! Please input "from_beam"!'
                return
            endelse
        endif
    endif
    
    ; 
    if n_elements(from_beam) eq 2 then begin
        from_beam = [ from_beam[0], from_beam[1], 0.0 ]
    endif
    
    ; 
    if n_elements(from_beam) eq 1 then begin
        from_beam = [ from_beam, from_beam, 0.0 ]
    endif
    
    ; 
    if n_elements(to_beam) eq 2 then begin
        to_beam = [ to_beam[0], to_beam[1], 0.0 ]
    endif
    
    ; 
    if n_elements(to_beam) eq 1 then begin
        to_beam = [ to_beam, to_beam, 0.0 ]
    endif
    
    ; 
    from_beam_str = 'beam ' + string(format='(F0.3)',from_beam[0]) + ' x ' + string(format='(F0.3)',from_beam[1]) + ' arcsec, PA ' + string(format='(F0.3)',from_beam[2]) + ' degree'
    to_beam_str = 'beam ' + string(format='(F0.3)',to_beam[0]) + ' x ' + string(format='(F0.3)',to_beam[1]) + ' arcsec, PA ' + string(format='(F0.3)',to_beam[2]) + ' degree'
    
    
    ; 
    message, 'Convolving "'+strtrim(fits_file,2)+'" from '+from_beam_str+' to '+to_beam_str, /continue
    
    ; 
    pix_scale = double(fxpar(fits_header, 'CDELT2'))
    if pix_scale le 0 then pix_scale = double(fxpar(fits_header, 'CD2_2'))
    if pix_scale le 0 then message, 'Error! Failed to get pix_scale!'
    pix_scale = pix_scale * 3600.0
    message, 'pix_scale = '+strtrim(string(pix_scale),2)+' arcsec', /continue
    
    ; 
    NPIXEL = max([from_beam[0], from_beam[1], to_beam[0], to_beam[1]]) / pix_scale * 10
    NPIXEL = long(NPIXEL/2)*2+1
    if NPIXEL gt fits_key_naxis1 then NPIXEL = fits_key_naxis1
    if NPIXEL gt fits_key_naxis2 then NPIXEL = fits_key_naxis2
    if NPIXEL lt 255 then NPIXEL = 255
    if NPIXEL lt 511 then NPIXEL = 511
    ;if NPIXEL lt 511 then NPIXEL = 511
    message, 'NPIXEL = '+strtrim(string(NPIXEL),2), /continue
    
    ; 
    ; Generate to 2D Gaussian
    ; 
    ;Gauss1 = PSF_Gaussian(NPIXEL=NPIXEL, FWHM=from_beam, /double) ; use XY_CORREL to turn the angle
    ;Gauss2 = PSF_Gaussian(NPIXEL=NPIXEL, FWHM=to_beam, /double) ; use XY_CORREL to turn the angle
    Gauss1 = CrabImageGaussian2D(NPIXEL=NPIXEL, FWHM=from_beam[0:1], PA=from_beam[2], LOG10GAUSSIAN=LOG10Gauss1)
    Gauss2 = CrabImageGaussian2D(NPIXEL=NPIXEL, FWHM=to_beam[0:1], PA=to_beam[2], LOG10GAUSSIAN=LOG10Gauss2)
    ;message, 'CrabImageGaussian2D done', /continue
    
    ; 
    ; Convert 2D Guassian to Fourier space
    ; 
    FTGauss1 = FFT(Gauss1, /CENTER, /DOUBLE) ; forward
    FTGauss2 = FFT(Gauss2, /CENTER, /DOUBLE) ; forward
    ;FTGaussKernel = DOUBLE(FTGauss2)/DOUBLE(FTGauss1) ; here reports warning
    LOG10FTGauss1 = ALOG10(SQRT(real_part(FTGauss1)^2 + imaginary(FTGauss1)^2))
    LOG10FTGauss2 = ALOG10(SQRT(real_part(FTGauss2)^2 + imaginary(FTGauss2)^2))
    FTGaussKernel = 10.0D^(LOG10FTGauss2-LOG10FTGauss1)
    FTGaussKernel = DOUBLE(FTGaussKernel)/TOTAL(ABS(FTGaussKernel),/DOUBLE)
    GaussKernel = FFT(FTGaussKernel, /CENTER, /DOUBLE, /INVERSE)
    
    GaussKernel = SHIFT(GaussKernel, -([NPIXEL,NPIXEL]+1)/2)
    
    ;plot_var = Image(Gauss1, LAYOUT=[2,2,1]) ; top left
    ;plot_var = Image(FTGauss1, LAYOUT=[2,2,2], /CURRENT) ; top right
    ;plot_var = Image(Gauss2, LAYOUT=[2,2,3], /CURRENT) ; bottom left
    ;plot_var = Image(FTGauss2, LAYOUT=[2,2,4], /CURRENT) ; bottom right
    
    ;print, size(SQRT(real_part(FTGauss1)^2 + imaginary(FTGauss1)^2), /tname)
    mwrfits, SQRT(real_part(FTGauss1)^2 + imaginary(FTGauss1)^2), 'Gauss1FT.fits', /create
    mwrfits, SQRT(real_part(FTGauss2)^2 + imaginary(FTGauss2)^2), 'Gauss2FT.fits', /create
    mwrfits, SQRT(real_part(FTGaussKernel)^2 + imaginary(FTGaussKernel)^2), 'GaussKernelFT.fits', /create
    mwrfits, ALOG10(SQRT(real_part(FTGaussKernel)^2 + imaginary(FTGaussKernel)^2)), 'GaussKernelFTLOG10.fits', /create
    mwrfits, Gauss1, 'Gauss1.fits', /create
    mwrfits, Gauss2, 'Gauss2.fits', /create
    mwrfits, LOG10Gauss1, 'Gauss1LOG10.fits', /create
    mwrfits, LOG10Gauss2, 'Gauss2LOG10.fits', /create
    mwrfits, LOG10FTGauss1, 'Gauss1FTLOG10.fits', /create
    mwrfits, LOG10FTGauss2, 'Gauss2FTLOG10.fits', /create
    mwrfits, SQRT(real_part(GaussKernel)^2 + imaginary(GaussKernel)^2), 'GaussKernel.fits', /create
    ; -- we can visually check these images with ds9 softwares:
    ; -- ds9 -tile grid layout 2 3 -lock frame image Gauss1.fits Gauss1FT.fits Gauss2.fits Gauss2FT.fits GaussKernel.fits GaussKernelFT.fits -frame last -scale log
    ; -- ds9 -tile grid layout 2 3 -lock frame image Gauss1LOG10.fits Gauss1FTLOG10.fits Gauss2LOG10.fits Gauss2FTLOG10.fits GaussKernel.fits GaussKernelFTLOG10.fits -frame last -scale log
    
    convolved_image = convol_fft(fits_image_nonan, GaussKernel, KERNEL_FFT=FTGaussKernel)
    ;convolved_image = convolve(fits_image_nonan, GaussKernel, FT_PSF=FTGaussKernel)
    
    convolved_image = convolved_image / total(convolved_image,/NAN,/DOUBLE) * total(fits_image_nonan,/DOUBLE)
    
    output_header = fits_header
    sxaddpar, output_header, 'BMAJ', to_beam[0]/3600.0, 'degree'
    sxaddpar, output_header, 'BMIN', to_beam[1]/3600.0, 'degree'
    sxaddpar, output_header, 'BPA', to_beam[2], 'degree'
    sxaddhist, '', output_header
    sxaddhist, SYSTIME(), output_header
    sxaddhist, 'Convolved from '+from_beam_str+' ', output_header
    sxaddhist, +'to '+to_beam_str+' ', output_header
    sxaddhist, 'with SpireLines script FTS_tool_convolve_fits_image.', output_header
    sxaddhist, '', output_header
    
    convolved_image_a = SQRT(real_part(convolved_image)^2 + imaginary(convolved_image)^2)
    convolved_image_a[where(~finite(fits_image),/null)] = double('nan')
    ;mwrfits, convolved_image_a, 'convolved_image.fits', /create
    if n_elements(output_file) eq 0 then begin
        output_file = strmid(fits_file,0,strpos(fits_file,'.fits'))+'_Convolved_to_'+string(format='(F0.3)',to_beam[0])+'_arcsec_beam.fits'
    endif else begin
        if output_file eq '' then begin
            output_file = strmid(fits_file,0,strpos(fits_file,'.fits'))+'_Convolved_to_'+string(format='(F0.3)',to_beam[0])+'_arcsec_beam.fits'
        endif
    endelse
    mwrfits, convolved_image_a, output_file, output_header, /create
    print, 'Output to "'+output_file+'"!'
    
    ;Plot_device = image(Gauss1) & Plot_device.save, 'Plot_Gauss1.pdf'
    ;Plot_device = image(Gauss2) & Plot_device.save, 'Plot_Gauss2.pdf'
    ;Plot_device = image(GaussKernel) & Plot_device.save, 'Plot_GaussKernel.pdf'
    
    
END