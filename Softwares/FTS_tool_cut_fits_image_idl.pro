PRO FTS_tool_cut_fits_image_idl, fits_file, RA, Dec, FoV_X, FoV_Y, output_file = output_file
    
    ; 
    fits_image = mrdfits(fits_file, 0, fits_header)
    fits_image_nonan = fits_image
    fits_image_nonan[where(~finite(fits_image),/null)] = 0.0
    
    ; 
    dimension_x = fxpar(fits_header, 'NAXIS1')
    dimension_y = fxpar(fits_header, 'NAXIS2')
    
    ; 
    extast, fits_header, fits_wcs
    ad2xy, RA, Dec, fits_wcs, X, Y
    
    ; compute pixscale
    xy2ad, X+10.0, Y, fits_wcs, RA2, Dec2
    pixscale_X = SQRT((Dec2-Dec)^2 + ((RA2-RA)*COS(Dec/180.0*!PI))^2) / 10.0 * 3600.0 ; arcsec
    xy2ad, X, Y+10.0, fits_wcs, RA2, Dec2
    pixscale_Y = SQRT((Dec2-Dec)^2 + ((RA2-RA)*COS(Dec/180.0*!PI))^2) / 10.0 * 3600.0 ; arcsec
    
    ; 
    print, 'pixscale_X =', pixscale_X
    print, 'pixscale_Y =', pixscale_Y
    side_half_length_x = FoV_X / pixscale_X / 2.0
    side_half_length_y = FoV_Y / pixscale_Y / 2.0
    
    ; 
    IF side_half_length_x LT 0 THEN BEGIN
        message, 'Cutting "'+strtrim(fits_file,2)+'" Failed! FoV_X smaller than one pixel!'
    ENDIF
    IF side_half_length_y LT 0 THEN BEGIN
        message, 'Cutting "'+strtrim(fits_file,2)+'" Failed! FoV_Y smaller than one pixel!'
    ENDIF
    
    ; 
    X1 = FIX(ROUND(X - (side_half_length_x-0.5))) ;  - 0.5
    Y1 = FIX(ROUND(Y - (side_half_length_y-0.5))) ;  - 0.5
    X2 = FIX(ROUND(X + (side_half_length_x-0.5))) ;  + 0.5
    Y2 = FIX(ROUND(Y + (side_half_length_y-0.5))) ;  + 0.5
    
    ; 
    IF X1 LT 0 OR X2 GE dimension_x OR Y1 LT 0 OR Y2 GE dimension_y THEN BEGIN
        message, 'Cutting "'+strtrim(fits_file,2)+'" Failed! FoV out of the input image ([X1,Y1,X2,Y2]='+'['+strtrim(string(X1),2)+','+strtrim(string(Y1),2)+','+strtrim(string(X2),2)+','+strtrim(string(Y2),2)+']'+' vs '+'['+strtrim(string(0),2)+','+strtrim(string(0),2)+strtrim(string(dimension_x-1),2)+','+strtrim(string(dimension_y-1),2)+']'+')!'
    ENDIF
    
        
    ; 
    message, 'Cutting "'+strtrim(fits_file,2)+'" to [X1,Y1,X2,Y2]='+'['+strtrim(string(X1),2)+','+strtrim(string(Y1),2)+','+strtrim(string(X2),2)+','+strtrim(string(Y2),2)+'].', /continue
    
    ; 
    hextract, fits_image, fits_header, new_image, new_header, X1, X2, Y1, Y2
    
    if n_elements(output_file) eq 0 then output_file = strmid(fits_file,0,strpos(fits_file,'.fits'))+'_Cutout.fits'
    
    mwrfits, new_image, output_file, new_header, /create
    print, 'Output to "'+output_file+'"!'
    
END