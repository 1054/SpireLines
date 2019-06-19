#!/usr/bin/env python3
# 

import os, sys, re
import astropy
import numpy as np
from numpy import isfinite, sqrt, log, log10, pi; finite = isfinite
from numpy.fft import fft2, ifft2, rfft2, irfft2, fftshift, ifftshift
from numpy import real, imag, abs, sum, nansum, convolve
from astropy.io import fits
from datetime import datetime
from astropy.modeling.models import Gaussian2D
from scipy.signal import convolve2d, fftconvolve

if sys.version_info.major>=3:
    long = int
    double = np.float64
fix = long
string = str
n_elements = len
ALOG10 = log10
SQRT = sqrt
real_part = real
imaginary = imag
total = nansum
TOTAL = nansum
ABS = abs
DOUBLE = double


def fxpar(fits_header, key_name):
    if key_name in fits_header:
        return str(fits_header[key_name])
    return ''


def mwrfits(image_data, fits_file, header_data = None, overwrite=True):
    hdu = fits.PrimaryHDU(image_data, header=header_data)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(fits_file, overwrite=overwrite)


def SYSTIME():
    return datetime.today().strftime('%Y-%m-%d %Hh%Mm%Ss %Z')


def FTS_tool_convolve_fits_image(fits_file, from_beam = [], to_beam = [], output_file = ''):
    
    # 
    fits_image = fits.getdata(fits_file)
    fits_header = fits.getheader(fits_file)
    fits_image_nonan = fits_image
    fits_image_nonan[~finite(fits_image)] = 0.0
    
    # 
    fits_key_naxis1 = fxpar(fits_header, 'NAXIS1')
    if fits_key_naxis1 != '':
        fits_key_naxis1 = fix(fits_key_naxis1) 
    else:
        fits_key_naxis1 = float('nan') # pixel
    # 
    fits_key_naxis2 = fxpar(fits_header, 'NAXIS2')
    if fits_key_naxis2 != '':
        fits_key_naxis2 = fix(fits_key_naxis2) 
    else:
        fits_key_naxis2 = float('nan') # pixel
    # 
    fits_key_bmaj = fxpar(fits_header, 'BMAJ')
    if fits_key_bmaj != '':
        fits_key_bmaj = double(fits_key_bmaj) * 3600.0 
    else:
        fits_key_bmaj = float('nan') # arcsec
    fits_key_bmin = fxpar(fits_header, 'BMIN')
    if fits_key_bmin != '':
        fits_key_bmin = double(fits_key_bmin) * 3600.0 
    else:
        fits_key_bmin = float('nan') # arcsec
    fits_key_bpa = fxpar(fits_header, 'BPA')
    if fits_key_bpa != '':
        fits_key_bpa = double(fits_key_bpa) 
    else:
        fits_key_bpa = float('nan') # degree
    # 
    if n_elements(from_beam) == 1:
        if from_beam[0] <= 0:
            if finite(fits_key_bmaj) and finite(fits_key_bmin):
                if finite(fits_key_bpa):
                    from_beam = [fits_key_bmaj, fits_key_bmin, fits_key_bpa]
                else:
                    from_beam = [fits_key_bmaj, fits_key_bmin, 0.0]
            elif finite(fits_key_bmaj):
                from_beam = [fits_key_bmaj, fits_key_bmaj, 0.0]
            else:
                print('Error! No BMAJ BMIN BPA keyword found! Please input "from_beam"!')
                return
    
    # 
    if n_elements(from_beam) == 2:
        from_beam = [ from_beam[0], from_beam[1], 0.0 ]
    
    # 
    if n_elements(from_beam) == 1:
        from_beam = [ from_beam[0], from_beam[0], 0.0 ]
    
    # 
    if n_elements(to_beam) == 2:
        to_beam = [ to_beam[0], to_beam[1], 0.0 ]
    
    # 
    if n_elements(to_beam) == 1:
        to_beam = [ to_beam[0], to_beam[0], 0.0 ]
    
    # 
    from_beam_str = 'beam %0.3f x %0.3f arcsec, PA %0.3f degree'%(from_beam[0],from_beam[1],from_beam[2])
    to_beam_str = 'beam %0.3f x %0.3f arcsec, PA %0.3f degree'%(to_beam[0],to_beam[1],to_beam[2])
    
    
    # 
    print('Convolving "'+fits_file+'" from '+from_beam_str+' to '+to_beam_str )
    
    # 
    pix_scale = double(fxpar(fits_header, 'CDELT2'))
    if pix_scale <= 0: pix_scale = double(fxpar(fits_header, 'CD2_2'))
    if pix_scale <= 0: raise Exception('Error! Failed to get pix_scale!')
    pix_scale = pix_scale * 3600.0
    print('pix_scale = %0.6f arcsec'%(pix_scale) )
    
    # 
    NPIXEL = max([from_beam[0], from_beam[1], to_beam[0], to_beam[1]]) / pix_scale * 10.0
    NPIXEL = long(NPIXEL/2)*2+1
    if NPIXEL > fits_key_naxis1: NPIXEL = fits_key_naxis1
    if NPIXEL > fits_key_naxis2: NPIXEL = fits_key_naxis2
    if NPIXEL < 255: NPIXEL = 255
    if NPIXEL < 511: NPIXEL = 511
    if NPIXEL < 1023: NPIXEL = 1023
    #if NPIXEL < 2047: NPIXEL = 2047
    #if NPIXEL < 3095: NPIXEL = 3095
    print('NPIXEL = %s'%(NPIXEL) )
    
    # 
    # Generate to 2D Gaussian
    # 
    xc = (NPIXEL+1.0)/2.0
    yc = (NPIXEL+1.0)/2.0
    xv, yv = np.meshgrid(np.arange(1.0,NPIXEL+1.0,1.0), np.arange(1.0,NPIXEL+1.0,1.0))
    xv = xv - xc
    yv = yv - yc
    print('xv.shape', xv.shape)
    print('yv.shape', yv.shape)
    Gauss1Func = Gaussian2D(amplitude=1.0, x_mean=0.0, y_mean=0.0, x_stddev=from_beam[0]/(2.0*sqrt(2.0*log(2.0))), y_stddev=from_beam[1]/(2.0*sqrt(2.0*log(2.0))), theta=(from_beam[2]+90.0)/180.0*pi)
    Gauss2Func = Gaussian2D(amplitude=1.0, x_mean=0.0, y_mean=0.0, x_stddev=to_beam[0]/(2.0*sqrt(2.0*log(2.0))), y_stddev=to_beam[1]/(2.0*sqrt(2.0*log(2.0))), theta=(to_beam[2]+90.0)/180.0*pi)
    Gauss1 = Gauss1Func(xv, yv)
    Gauss2 = Gauss2Func(xv, yv)
    Gauss1 = Gauss1 / total(Gauss1)
    Gauss2 = Gauss2 / total(Gauss2)
    LOG10Gauss1 = ALOG10(Gauss1)
    LOG10Gauss2 = ALOG10(Gauss2)
    
    # 
    # Convert 2D Guassian to Fourier space
    # 
    FTGauss1 = fftshift(fft2(Gauss1)) # forward
    FTGauss2 = fftshift(fft2(Gauss2)) # forward
    FTGauss1 = ABS(FTGauss1) # ABS, because Gaussian FFT is also a Gaussian with zero imaginary
    FTGauss2 = ABS(FTGauss2) # ABS, because Gaussian FFT is also a Gaussian with zero imaginary
    FTGauss1[(FTGauss1<3E-16)] = 3E-16
    FTGauss2[(FTGauss2<3E-16)] = 3E-16
    LOG10FTGauss1 = ALOG10(ABS(FTGauss1))
    LOG10FTGauss2 = ALOG10(ABS(FTGauss2))
    #FTGaussKernel = DOUBLE(FTGauss2)/DOUBLE(FTGauss1) # here reports warning
    FTGaussKernel = 10.0**(LOG10FTGauss2-LOG10FTGauss1)
    FTGaussKernel[np.logical_and((FTGauss1<3E-16),(FTGauss2<3E-16))] = 1.0
    FTGaussKernel = FTGaussKernel/TOTAL(ABS(FTGaussKernel))
    GaussKernel = fftshift(ifft2(ifftshift(FTGaussKernel)))
    GaussKernel = ABS(GaussKernel)
    print('Gauss1.shape', Gauss1.shape)
    print('FTGauss1.shape', FTGauss1.shape)
    print('Gauss2.shape', Gauss2.shape)
    print('FTGauss2.shape', FTGauss2.shape)
    print('GaussKernel.shape', GaussKernel.shape)
    print('FTGaussKernel.shape', FTGaussKernel.shape)
    
    #GaussKernel = SHIFT(GaussKernel, -([NPIXEL,NPIXEL]+1)/2)
    
    #plot_var = Image(Gauss1, LAYOUT=[2,2,1]) # top left
    #plot_var = Image(FTGauss1, LAYOUT=[2,2,2], /CURRENT) # top right
    #plot_var = Image(Gauss2, LAYOUT=[2,2,3], /CURRENT) # bottom left
    #plot_var = Image(FTGauss2, LAYOUT=[2,2,4], /CURRENT) # bottom right
    
    #print, size(SQRT(real_part(FTGauss1)**2 + imaginary(FTGauss1)**2), /tname)
    mwrfits(SQRT(real_part(FTGauss1)**2 + imaginary(FTGauss1)**2), 'Gauss1FT.fits' )
    mwrfits(SQRT(real_part(FTGauss2)**2 + imaginary(FTGauss2)**2), 'Gauss2FT.fits' )
    mwrfits(SQRT(real_part(FTGaussKernel)**2 + imaginary(FTGaussKernel)**2), 'GaussKernelFT.fits' )
    mwrfits(ALOG10(SQRT(real_part(FTGaussKernel)**2 + imaginary(FTGaussKernel)**2)), 'GaussKernelFTLOG10.fits' )
    mwrfits(Gauss1, 'Gauss1.fits' )
    mwrfits(Gauss2, 'Gauss2.fits' )
    mwrfits(LOG10Gauss1, 'Gauss1LOG10.fits' )
    mwrfits(LOG10Gauss2, 'Gauss2LOG10.fits' )
    mwrfits(LOG10FTGauss1, 'Gauss1FTLOG10.fits' )
    mwrfits(LOG10FTGauss2, 'Gauss2FTLOG10.fits' )
    mwrfits(SQRT(real_part(GaussKernel)**2 + imaginary(GaussKernel)**2), 'GaussKernel.fits' )
    # -- we can visually check these images with ds9 softwares:
    # -- ds9 -tile grid layout 2 3 -lock frame image Gauss1.fits Gauss1FT.fits Gauss2.fits Gauss2FT.fits GaussKernel.fits GaussKernelFT.fits -frame last -scale log
    # -- ds9 -tile grid layout 2 3 -lock frame image Gauss1LOG10.fits Gauss1FTLOG10.fits Gauss2LOG10.fits Gauss2FTLOG10.fits GaussKernel.fits GaussKernelFTLOG10.fits -frame last -scale log
    
    #convolved_image = convol_fft(fits_image_nonan, GaussKernel, KERNEL_FFT=FTGaussKernel)
    #convolved_image = convolve(fits_image_nonan, GaussKernel, FT_PSF=FTGaussKernel)
    #FT_fits_image_nonan = fftshift(rfft2(fits_image_nonan)) # forward
    #convolved_image = fftshift(irfft2(ifftshift(FT_fits_image_nonan * FTGaussKernel)))
    #convolved_image = convolve2d(fits_image_nonan, GaussKernel, mode='full', boundary='fill', fillvalue=0)
    print('type(fits_image_nonan)', type(fits_image_nonan))
    print('type(GaussKernel)', type(GaussKernel))
    print('fits_image_nonan.dtype', fits_image_nonan.dtype)
    print('GaussKernel.dtype', GaussKernel.dtype)
    print('fits_image_nonan.shape', fits_image_nonan.shape)
    print('GaussKernel.shape', GaussKernel.shape)
    # fftconvolve
    # -- https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.fftconvolve.html
    convolved_image = fftconvolve(fits_image_nonan.astype(double), GaussKernel, mode='same') 
    print('fits_image_nonan.shape', fits_image_nonan.shape)
    print('GaussKernel.shape', GaussKernel.shape)
    print('convolved_image.shape', convolved_image.shape)
    
    convolved_image = convolved_image / total(convolved_image) * total(fits_image_nonan)
    
    output_header = fits_header
    output_header['BMAJ'] = (to_beam[0]/3600.0, 'degree')
    output_header['BMIN'] = (to_beam[1]/3600.0, 'degree')
    output_header['BPA'] = (to_beam[2], 'degree')
    output_header['HISTORY'] = ''
    output_header['HISTORY'] = SYSTIME()
    output_header['HISTORY'] = 'Convolved from '+from_beam_str+' '
    output_header['HISTORY'] = 'to '+to_beam_str+' '
    output_header['HISTORY'] = 'with SpireLines script FTS_tool_convolve_fits_image.'
    output_header['HISTORY'] = ''
    
    convolved_image_a = SQRT(real_part(convolved_image)**2 + imaginary(convolved_image)**2)
    convolved_image_a[~finite(fits_image)] = double('nan')
    if n_elements(output_file) == 0:
        output_file = re.sub(r'\.fits$',r'',fits_file)+'_Convolved_to_%0.3f'%(to_beam[0])+'_arcsec_beam.fits'
    else: 
        if output_file == '':
            output_file = re.sub(r'\.fits$',r'',fits_file)+'_Convolved_to_%0.3f'%(to_beam[0])+'_arcsec_beam.fits'
    mwrfits(convolved_image_a, output_file, output_header)
    print('Output to "'+output_file+'"!')
    
    #Plot_device = image(Gauss1) & Plot_device.save, 'Plot_Gauss1.pdf'
    #Plot_device = image(Gauss2) & Plot_device.save, 'Plot_Gauss2.pdf'
    #Plot_device = image(GaussKernel) & Plot_device.save, 'Plot_GaussKernel.pdf'










# 
# main
# 
if __name__ == '__main__':
    
    if len(sys.argv) <= 1:
        print('Usage: FTS_tool_convolve_fits_image')
        sys.exit()
    
    fits_file = ''
    from_beam = []
    to_beam = []
    output_file = ''
    argmode = ''
    iarg = 1
    while iarg < len(sys.argv):
        if sys.argv[iarg] == '-out':
            argmode = 'out'
            iarg += 1
            continue
        elif sys.argv[iarg] == '-from':
            argmode = 'from'
            iarg += 1
            continue
        elif sys.argv[iarg] == '-to':
            argmode = 'to'
            iarg += 1
            continue
        # 
        if argmode == 'out':
            output_file = sys.argv[iarg]
            argmode = ''
        elif argmode == 'from':
            from_beam.append(float(sys.argv[iarg]))
        elif argmode == 'to':
            to_beam.append(float(sys.argv[iarg]))
        else:
            if fits_file == '':
                fits_file = sys.argv[iarg]
        iarg += 1
    
    print('fits_file', fits_file)
    print('from_beam', from_beam)
    print('to_beam', to_beam)
    print('output_file', output_file)
    
    FTS_tool_convolve_fits_image(fits_file, from_beam, to_beam, output_file)
    
    







    
