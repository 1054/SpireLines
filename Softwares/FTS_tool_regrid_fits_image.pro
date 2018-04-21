PRO FTS_tool_regrid_fits_image, fits_file, x = to_dimension_x, y = to_dimension_y
    
    ; 
    fits_image = mrdfits(fits_file, 0, fits_header)
    fits_image_nonan = fits_image
    fits_image_nonan[where(~finite(fits_image),/null)] = 0.0
    
    ; 
    from_dimension_x = fxpar(fits_header, 'NAXIS1')
    from_dimension_y = fxpar(fits_header, 'NAXIS2')
    
    ; 
    if to_dimension_x le 0 and to_dimension_y ne 0 then begin
        to_dimension_x = fix(to_dimension_y * float(from_dimension_x)/float(from_dimension_y) + 0.5)
    endif else if to_dimension_y le 0 and to_dimension_x ne 0 then begin
        to_dimension_y = fix(to_dimension_x * float(from_dimension_y)/float(from_dimension_x) + 0.5)
    endif else if to_dimension_x le 0 and to_dimension_y le 0 then begin
        message, 'Error! Please input positive dimensions!'
        return
    endif
        
    ; 
    message, 'Regridding "'+strtrim(fits_file,2)+'" to '+strtrim(string(to_dimension_x),2)+'x'+strtrim(string(to_dimension_y),2)+' dimension.', /continue
    
    ; 
    pix_scale = double(fxpar(fits_header, 'CDELT2'))
    if pix_scale le 0 then pix_scale = double(fxpar(fits_header, 'CD2_2'))
    if pix_scale le 0 then message, 'Error! Failed to get pix_scale!'
    pix_scale = pix_scale * 3600.0
    message, 'pix_scale = '+strtrim(string(pix_scale),2)+' arcsec', /continue
    
    ; 
    regridded_image = frebin(fits_image, to_dimension_x, to_dimension_y)
    
    output_header = fits_header
    sxaddpar, output_header, 'NAXIS1', to_dimension_x
    sxaddpar, output_header, 'NAXIS2', to_dimension_y
    sxaddhist, '', output_header
    sxaddhist, SYSTIME(), output_header
    sxaddhist, 'Regridded from '+string(format='(I0)',from_dimension_x)+'x'+string(format='(I0)',from_dimension_y)+' to '+string(format='(I0)',to_dimension_x)+'x'+string(format='(I0)',to_dimension_y), output_header
    sxaddhist, 'with SpireLines script FTS_tool_regrid_fits_image.', output_header
    sxaddhist, '', output_header
    
    if fxpar(fits_header, 'CDELT2') ne '' then begin
        sxaddpar, output_header, 'CDELT2', double(fxpar(fits_header, 'CDELT2')) * float(from_dimension_y) / float(to_dimension_y)
    endif
    if fxpar(fits_header, 'CDELT1') ne '' then begin
        sxaddpar, output_header, 'CDELT1', double(fxpar(fits_header, 'CDELT1')) * float(from_dimension_x) / float(to_dimension_x)
    endif
    if fxpar(fits_header, 'CD2_2') ne '' then begin
        sxaddpar, output_header, 'CD2_2', double(fxpar(fits_header, 'CD2_2')) * float(from_dimension_y) / float(to_dimension_y)
    endif
    if fxpar(fits_header, 'CD1_1') ne '' then begin
        sxaddpar, output_header, 'CD1_1', double(fxpar(fits_header, 'CD1_1')) * float(from_dimension_x) / float(to_dimension_x)
    endif
    if fxpar(fits_header, 'CRPIX2') ne '' then begin
        sxaddpar, output_header, 'CRPIX2', (double(fxpar(fits_header, 'CRPIX2')) - 0.5) / float(from_dimension_y) * float(to_dimension_y) + 0.5
    endif
    if fxpar(fits_header, 'CRPIX1') ne '' then begin
        sxaddpar, output_header, 'CRPIX1', (double(fxpar(fits_header, 'CRPIX1')) - 0.5) / float(from_dimension_x) * float(to_dimension_x) + 0.5
    endif
    
    mwrfits, regridded_image, strmid(fits_file,0,strpos(fits_file,'.fits'))+'_Regridded.fits', output_header, /create
    
END