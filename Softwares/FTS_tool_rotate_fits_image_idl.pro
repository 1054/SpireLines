PRO FTS_tool_rotate_fits_image_idl, fits_file, angle_in_degrees, output_file = output_file
    
    ; 
    fits_image = mrdfits(fits_file, 0, fits_header)
    fits_image_nonan = fits_image
    fits_image_nonan[where(~finite(fits_image),/null)] = 0.0
    
    ; 
    ;dimension_x = fxpar(fits_header, 'NAXIS1')
    ;dimension_y = fxpar(fits_header, 'NAXIS2')
    
    ; 
    
        
    ; 
    message, 'Rotating "'+strtrim(fits_file,2)+'" by '+strtrim(string(angle_in_degrees),2)+' degrees.', /continue
    
    ; 
    xc = -1
    yc = -1
    interp = 1 ; 0 for nearest neighbor, 1 for bilinear interpolation, 2 for cubic interpolation.
    hrot, fits_image, fits_header, new_image, new_header, angle_in_degrees, xc, yc, interp
    
    if n_elements(output_file) eq 0 then output_file = strmid(fits_file,0,strpos(fits_file,'.fits'))+'_Rotated.fits'
    
    mwrfits, new_image, output_file, new_header, /create
    print, 'Output to "'+output_file+'"!'
    
END