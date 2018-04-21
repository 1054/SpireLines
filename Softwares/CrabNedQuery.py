#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# use astropy query 
# to check if the input 
# galaxies are AGN or not
# 
# 

def getInfoAGN(InputList, OutputInfoAGN=[], OutputInfoMetallicity=[]):
    # 
    # check input type
    if type(InputList) is not list:
        InputList = [InputList]
    # 
    # prepare output array
    OutputTypeAGN = ["Unknown"]*len(InputList)
    # 
    # loop each Input 
    for oio in range(len(InputList)):
        # 
        InputItem = InputList[oio]
        # 
        QueryObject = Ned.query_object(InputItem)
        # print(QueryObject) # an astropy.table.Table
        QueryPhotometry = Ned.get_table(InputItem, table='photometry') # Must be one of [‘photometry’|’positions’|’diameters’|’redshifts’|’references’|’object_notes’]
        # print QueryPhotometry[['Observed Passband','Units']] # ('No.', 'Observed Passband', 'Photometry Measurement', 'Uncertainty', 'Units', 'Frequency', 'NED Photometry Measurement', 'NED Uncertainty', 'NED Units', 'Refcode', 'Significance', 'Published frequency', 'Frequency Mode', 'Coordinates Targeted', 'Spatial Mode', 'Qualifiers', 'Comments')
        
        # 
        # extract OIII5007 Hbeta4861 Halpha6563 NII6583
        LineFlux_OII3727 = 0.0
        LineFlux_OIII4945 = 0.0
        LineFlux_OIII5007 = 0.0
        LineFlux_Hbeta4861 = 0.0
        LineFlux_Halpha6563 = 0.0
        LineFlux_NII6583 = 0.0
        for i in range(len(QueryPhotometry)):
            if LineFlux_OII3727 <= 0 and QueryPhotometry['Observed Passband'][i].startswith("[O II] 3727"):
                LineFlux_OII3727 = QueryPhotometry['Photometry Measurement'][i]
                print "LineFlux", QueryPhotometry['Observed Passband'][i], QueryPhotometry['Photometry Measurement'][i], QueryPhotometry['Units'][i]
            if LineFlux_OIII4945 <= 0 and QueryPhotometry['Observed Passband'][i].startswith("[O III] 4945"):
                LineFlux_OIII4945 = QueryPhotometry['Photometry Measurement'][i]
                print "LineFlux", QueryPhotometry['Observed Passband'][i], QueryPhotometry['Photometry Measurement'][i], QueryPhotometry['Units'][i]
            if LineFlux_OIII5007 <= 0 and QueryPhotometry['Observed Passband'][i].startswith("[O III] 5007"):
                LineFlux_OIII5007 = QueryPhotometry['Photometry Measurement'][i]
                print "LineFlux", QueryPhotometry['Observed Passband'][i], QueryPhotometry['Photometry Measurement'][i], QueryPhotometry['Units'][i]
            if LineFlux_Hbeta4861 <= 0 and QueryPhotometry['Observed Passband'][i].startswith("H{beta}"):
                LineFlux_Hbeta4861 = QueryPhotometry['Photometry Measurement'][i]
                print "LineFlux", QueryPhotometry['Observed Passband'][i], QueryPhotometry['Photometry Measurement'][i], QueryPhotometry['Units'][i]
            if LineFlux_Halpha6563 <= 0 and QueryPhotometry['Observed Passband'][i].startswith("H{alpha}"):
                LineFlux_Halpha6563 = QueryPhotometry['Photometry Measurement'][i]
                print "LineFlux", QueryPhotometry['Observed Passband'][i], QueryPhotometry['Photometry Measurement'][i], QueryPhotometry['Units'][i]
            if LineFlux_NII6583 <= 0 and QueryPhotometry['Observed Passband'][i].startswith("[N II] 6584"):
                LineFlux_NII6583 = QueryPhotometry['Photometry Measurement'][i]
                print "LineFlux", QueryPhotometry['Observed Passband'][i], QueryPhotometry['Photometry Measurement'][i], QueryPhotometry['Units'][i]
            if LineFlux_NII6583 <= 0 and QueryPhotometry['Observed Passband'][i].startswith("[N II] 6583"):
                LineFlux_NII6583 = QueryPhotometry['Photometry Measurement'][i]
                print "LineFlux", QueryPhotometry['Observed Passband'][i], QueryPhotometry['Photometry Measurement'][i], QueryPhotometry['Units'][i]
        
        # 
        # <TODO> fix NED bugs
        if InputItem.startswith('NGC1068'):
            if LineFlux_OIII5007>1000.0:
                LineFlux_OIII5007 = LineFlux_OIII5007/1e14 # [O III] 5007 (ESO) 	35580.0E-1  <--> [O III] 5007 (ESO) 	35580.0E-15
        
        # 
        # compute OIII/Hbeta
        LineRatio_OIIIHbeta = 0.0
        if LineFlux_Hbeta4861 > 0.0 and LineFlux_OIII5007 > 0.0:
            LineRatio_OIIIHbeta = (LineFlux_OIII5007 / LineFlux_Hbeta4861)
            print "LineRatio_OIII_Hbeta = %0.3f"%(LineRatio_OIIIHbeta)
            LineRatio_OIIIHbeta_log10 = math.log10(LineFlux_OIII5007 / LineFlux_Hbeta4861)
            print "LineRatio_OIII_Hbeta_log10 = %0.3f dex"%(LineRatio_OIIIHbeta_log10)
        else:
            print "Warning! No Hbeta4861 found!"
            #return OutputTypeAGN
        
        # 
        # compute NII/Halpha
        LineRatio_NIIHalpha = 0.0
        if LineFlux_Halpha6563 > 0.0 and LineFlux_NII6583 > 0.0:
            LineRatio_NIIHalpha = (LineFlux_NII6583 / LineFlux_Halpha6563)
            print "LineRatio_NII_Halpha = %0.3f"%(LineRatio_NIIHalpha)
            LineRatio_NIIHalpha_log10 = math.log10(LineFlux_NII6583 / LineFlux_Halpha6563)
            print "LineRatio_NII_Halpha_log10 = %0.3f dex"%(LineRatio_NIIHalpha_log10)
        else:
            print "Warning! No Halpha6563 found!"
            #return OutputTypeAGN
        
        # 
        # check Seyfert or LINER -- Kauffmann et al. 2013 MNRAS 346 105 Fig.1 Caption -- http://adsabs.harvard.edu/abs/2003MNRAS.346.1055K
        if LineRatio_OIIIHbeta >= (3.0) and LineRatio_NIIHalpha >= (0.6):
            print "%s is %s"%(InputItem,"Seyfert")
            OutputTypeAGN[oio] = "Seyfert"
        if LineRatio_OIIIHbeta < (3.0) and LineRatio_NIIHalpha >= (0.6):
            print "%s is %s"%(InputItem,"LINER")
            OutputTypeAGN[oio] = "LINER"
        
        # 
        # write to file
        print OutputTypeAGN
        if LineFlux_Hbeta4861 > 0.0 and LineFlux_Halpha6563 > 0.0 and LineFlux_NII6583 > 0.0 and LineFlux_OIII5007 > 0.0:
            if type(OutputInfoAGN) is str:
                OutputFileUnit = open(OutputInfoAGN,'w')
                OutputFileUnit.write("LineRatio_OIII_Hbeta = %f\n"%(LineRatio_OIIIHbeta))
                OutputFileUnit.write("LineRatio_OIII_Hbeta_log10 = %f # dex\n"%(LineRatio_OIIIHbeta_log10))
                OutputFileUnit.write("LineRatio_NII_Halpha = %f\n"%(LineRatio_NIIHalpha))
                OutputFileUnit.write("LineRatio_NII_Halpha_log10 = %f # dex\n"%(LineRatio_NIIHalpha_log10))
                OutputFileUnit.write("AGN_Type = \"%s\"\n"%(OutputTypeAGN[oio]))
                OutputFileUnit.write("# %s type is %s\n"%(InputItem,OutputTypeAGN[oio]))
                OutputFileUnit.close()
        
        # 
        # compute 12+log(O/H)
        Metallicity = 0.0
        
        if Metallicity == 0.0 and LineRatio_OIIIHbeta > 0.0 and LineRatio_NIIHalpha > 0.0:
            if LineRatio_OIIIHbeta/LineRatio_NIIHalpha >= 2.0:
                Metallicity = 8.73 - 0.32 * math.log10(LineRatio_OIIIHbeta/LineRatio_NIIHalpha)
                if type(OutputInfoMetallicity) is str:
                    OutputFileUnit = open(OutputInfoMetallicity,'w')
                    OutputFileUnit.write("Metallicity = %f\n"%(Metallicity))
                    OutputFileUnit.write("Metallicity_method = \"12+log(O/H) = Pettini & Pagel (2004) - PP04 O3N2\" # http://www.ifa.hawaii.edu/~kewley/Metallicity/ms.pdf\n")
                    OutputFileUnit.close()
                print "Metallicity = %f"%(Metallicity)
                print "Metallicity_method = \"12+log(O/H) = Pettini & Pagel (2004) - PP04 O3N2\" # http://www.ifa.hawaii.edu/~kewley/Metallicity/ms.pdf"
        
        if Metallicity == 0.0 and LineRatio_NIIHalpha > 0.0:
            Metallicity_lgN2 = math.log10(LineRatio_NIIHalpha)
            if Metallicity_lgN2 > -2.5 and Metallicity_lgN2 < -0.3:
                Metallicity = 9.37 + 2.03*Metallicity_lgN2 + 1.26*Metallicity_lgN2**2 + 0.32*Metallicity_lgN2**3
                if type(OutputInfoMetallicity) is str:
                    OutputFileUnit = open(OutputInfoMetallicity,'w')
                    OutputFileUnit.write("Metallicity = %f\n"%(Metallicity))
                    OutputFileUnit.write("Metallicity_method = \"12+log(O/H) = Pettini & Pagel (2004) - PP04 N2\" # http://www.ifa.hawaii.edu/~kewley/Metallicity/ms.pdf\n")
                    OutputFileUnit.close()
                print "Metallicity = %f"%(Metallicity)
                print "Metallicity_method = \"12+log(O/H) = Pettini & Pagel (2004) - PP04 N2\" # http://www.ifa.hawaii.edu/~kewley/Metallicity/ms.pdf"
        
        if Metallicity == 0.0 and LineFlux_OII3727 > 0.0 and LineFlux_OIII4945 > 0.0 and LineFlux_OIII5007 > 0.0 and LineFlux_Hbeta4861 > 0.0:
            Metallicity_lgR23 = math.log10((LineFlux_OII3727+LineFlux_OIII4945+LineFlux_OIII5007)/LineFlux_Hbeta4861)
            Metallicity_lgO32 = math.log10((LineFlux_OIII4945+LineFlux_OIII5007)/LineFlux_OII3727)
            Metallicity_lower = 12-4.944+0.767*Metallicity_lgR23+0.602*Metallicity_lgR23**2-(Metallicity_lgO32)*(0.29+0.332*Metallicity_lgR23-0.331*Metallicity_lgR23**2)
            Metallicity_upper = 12-2.939-0.200*Metallicity_lgR23-0.237*Metallicity_lgR23**2-0.305*Metallicity_lgR23**3-0.0283*Metallicity_lgR23**4-(Metallicity_lgO32)*(0.0047-0.0221*Metallicity_lgR23-0.102*Metallicity_lgR23**2-0.0817*Metallicity_lgR23**3-0.00717*Metallicity_lgR23**4)
            Metallicity = (Metallicity_lower+Metallicity_upper)/2.0
            if type(OutputInfoMetallicity) is str:
                OutputFileUnit = open(OutputInfoMetallicity,'w')
                OutputFileUnit.write("Metallicity = %f\n"%(Metallicity))
                OutputFileUnit.write("Metallicity_method = \"12+log(O/H) = McGaugh (1991) - M91\" # http://www.ifa.hawaii.edu/~kewley/Metallicity/ms.pdf\n")
                OutputFileUnit.close()
            print "Metallicity = %f"%(Metallicity)
            print "Metallicity_method = \"12+log(O/H) = McGaugh (1991) - M91\" # http://www.ifa.hawaii.edu/~kewley/Metallicity/ms.pdf"
        
        if Metallicity == 0.0 and LineRatio_NIIHalpha > 0.0:
            Metallicity_lgN2 = math.log10(LineRatio_NIIHalpha)
            Metallicity = 9.12 + 0.73*Metallicity_lgN2
            if type(OutputInfoMetallicity) is str:
                OutputFileUnit = open(OutputInfoMetallicity,'w')
                OutputFileUnit.write("Metallicity = %f\n"%(Metallicity))
                OutputFileUnit.write("Metallicity_error = %f\n"%(Metallicity*0.2))
                OutputFileUnit.write("Metallicity_method = \"12+log(O/H) = Denicolo, Terlevich & Terlevich (2002) - D02\" # http://www.ifa.hawaii.edu/~kewley/Metallicity/ms.pdf\n")
                OutputFileUnit.close()
            print "Metallicity = %f"%(Metallicity)
            print "Metallicity_error = %f"%(Metallicity*0.2) # 0.2 dex large error
            print "Metallicity_method = \"12+log(O/H) = Denicolo, Terlevich & Terlevich (2002) - D02\" # http://www.ifa.hawaii.edu/~kewley/Metallicity/ms.pdf"
        
    # 
    # return
    return OutputTypeAGN



def getInfoOptical(InputList):
    # 
    # check input type
    if type(InputList) is not list:
        InputList = [InputList]
    # 
    # prepare output array
    OutputInfoOptical = numpy.recarray((0,), \
                                  dtype=[('wave', float), ('flux', float), ('error', float), ('passband', '|S255'), ('reference', '|S255')])
    # 
    # loop each Input 
    for InputItem in InputList:
        QueryObject = Ned.query_object(InputItem)
        # print(QueryObject) # an astropy.table.Table
        QueryPhotometry = Ned.get_table(InputItem, table='photometry') # Must be one of [‘photometry’|’positions’|’diameters’|’redshifts’|’references’|’object_notes’]
        # print QueryPhotometry[['Observed Passband','Units']] # ('No.', 'Observed Passband', 'Photometry Measurement', 'Uncertainty', 'Units', 'Frequency', 'NED Photometry Measurement', 'NED Uncertainty', 'NED Units', 'Refcode', 'Significance', 'Published frequency', 'Frequency Mode', 'Coordinates Targeted', 'Spatial Mode', 'Qualifiers', 'Comments')
        
        # 
        # extract Optical bands
        iOptical = []
        for i in range(len(QueryPhotometry)):
            if QueryPhotometry['Frequency Mode'][i].startswith("Broad-band "):
                if QueryPhotometry['Observed Passband'][i].startswith("u ") or \
                   QueryPhotometry['Observed Passband'][i].startswith("U ") or \
                   QueryPhotometry['Observed Passband'][i].startswith("B ") or \
                   QueryPhotometry['Observed Passband'][i].startswith("g ") or \
                   QueryPhotometry['Observed Passband'][i].startswith("V ") or \
                   QueryPhotometry['Observed Passband'][i].startswith("r ") or \
                   QueryPhotometry['Observed Passband'][i].startswith("R ") or \
                   QueryPhotometry['Observed Passband'][i].startswith("i ") or \
                   QueryPhotometry['Observed Passband'][i].startswith("I ") or \
                   QueryPhotometry['Observed Passband'][i].startswith("z ") or \
                   QueryPhotometry['Observed Passband'][i].startswith("j ") or \
                   QueryPhotometry['Observed Passband'][i].startswith("J ") or \
                   QueryPhotometry['Observed Passband'][i].startswith("J_") or \
                   QueryPhotometry['Observed Passband'][i].startswith("H ") or \
                   QueryPhotometry['Observed Passband'][i].startswith("H_") or \
                   QueryPhotometry['Observed Passband'][i].startswith("K ") or \
                   QueryPhotometry['Observed Passband'][i].startswith("K_") or \
                   QueryPhotometry['Observed Passband'][i].startswith("Ks ") or \
                   QueryPhotometry['Observed Passband'][i].startswith("Ks_"):
                    iOptical.extend([i])
                    print "Optical Flux", QueryPhotometry['Observed Passband'][i], QueryPhotometry['Photometry Measurement'][i], QueryPhotometry['Units'][i]
        
        # 
        # get Optical fluxes
        fOptical_Passband = (QueryPhotometry['Observed Passband'][iOptical])
        fOptical_FluxValue = (QueryPhotometry['NED Photometry Measurement'][iOptical])
        fOptical_FluxError = (QueryPhotometry['NED Uncertainty'][iOptical])
        fOptical_Frequency = (QueryPhotometry['Frequency'][iOptical])
        fOptical_Reference = QueryPhotometry['Refcode'][iOptical]
        
        # 
        # print
        for i in range(len(iOptical)):
            fOptical_FluxError[i] = fOptical_FluxError[i].replace('+/-','').replace('...','0.0')
            if len(fOptical_FluxError[i])==0: fOptical_FluxError[i]=0.0
            if fOptical_FluxError[i]=='...': fOptical_FluxError[i]=0.0 #<20180201>#
            fOptical_FluxError[i] = float(fOptical_FluxError[i])*1e3 # mJy
            fOptical_FluxValue[i] = float(fOptical_FluxValue[i])*1e3 # mJy
            fOptical_Wavelengthi = 2.99792458e5/(float(fOptical_Frequency[i])/1e9) # um
            print fOptical_Wavelengthi,fOptical_FluxValue[i],fOptical_FluxError[i],fOptical_Passband[i],fOptical_Reference[i]
            TemporaryOptical = numpy.array([(fOptical_Wavelengthi,fOptical_FluxValue[i],fOptical_FluxError[i],fOptical_Passband[i],fOptical_Reference[i])], dtype=OutputInfoOptical.dtype)
            OutputInfoOptical = rec.stack_arrays( (OutputInfoOptical, TemporaryOptical), usemask=False, asrecarray=True, autoconvert=True )
        
    # 
    # return
    return OutputInfoOptical



def getInfoInfrared(InputList):
    # 
    # check input type
    if type(InputList) is not list:
        InputList = [InputList]
    # 
    # prepare output array
    OutputInfoInfrared = numpy.recarray((0,), \
                                  dtype=[('wave', float), ('flux', float), ('error', float), ('passband', '|S255'), ('reference', '|S255')])
    # 
    # loop each Input 
    for InputItem in InputList:
        QueryObject = Ned.query_object(InputItem)
        # print(QueryObject) # an astropy.table.Table
        QueryPhotometry = Ned.get_table(InputItem, table='photometry') # Must be one of [‘photometry’|’positions’|’diameters’|’redshifts’|’references’|’object_notes’]
        # print QueryPhotometry[['Observed Passband','Units']] # ('No.', 'Observed Passband', 'Photometry Measurement', 'Uncertainty', 'Units', 'Frequency', 'NED Photometry Measurement', 'NED Uncertainty', 'NED Units', 'Refcode', 'Significance', 'Published frequency', 'Frequency Mode', 'Coordinates Targeted', 'Spatial Mode', 'Qualifiers', 'Comments')
        
        # 
        # extract Infrared bands
        iInfrared = []
        for i in range(len(QueryPhotometry)):
            if QueryPhotometry['Frequency Mode'][i].startswith("Broad-band "):
                if QueryPhotometry['Observed Passband'][i].startswith("7.9 microns Spitzer") or \
                   QueryPhotometry['Observed Passband'][i].startswith("8 microns (Spitzer)") or \
                   QueryPhotometry['Observed Passband'][i].startswith("10 microns (Spitzer)") or \
                   QueryPhotometry['Observed Passband'][i].startswith("IRAS 12 microns") or \
                   QueryPhotometry['Observed Passband'][i].startswith("14 microns (IRS)") or \
                   QueryPhotometry['Observed Passband'][i].startswith("15 microns (Spitzer)") or \
                   QueryPhotometry['Observed Passband'][i].startswith("20 microns (Spitzer)") or \
                   QueryPhotometry['Observed Passband'][i].startswith("25 microns (Spitzer)") or \
                   QueryPhotometry['Observed Passband'][i].startswith("IRAS 25 microns") or \
                   QueryPhotometry['Observed Passband'][i].startswith("30 microns (Spitzer)") or \
                   QueryPhotometry['Observed Passband'][i].startswith("IRAS 60 microns") or \
                   QueryPhotometry['Observed Passband'][i].startswith("IRAS 100 microns") or \
                   QueryPhotometry['Observed Passband'][i].startswith("350 microns (CSO)") or \
                   QueryPhotometry['Observed Passband'][i].startswith("850 microns (SCUBA)") or \
                   QueryPhotometry['Observed Passband'][i].startswith("Herschel"):
                    iInfrared.extend([i])
                    print "Infrared Flux", QueryPhotometry['Observed Passband'][i], QueryPhotometry['Photometry Measurement'][i], QueryPhotometry['Units'][i]
        
        # 
        # get Infrared fluxes
        fInfrared_Passband = (QueryPhotometry['Observed Passband'][iInfrared])
        fInfrared_FluxValue = (QueryPhotometry['NED Photometry Measurement'][iInfrared])
        fInfrared_FluxError = (QueryPhotometry['NED Uncertainty'][iInfrared])
        fInfrared_Frequency = (QueryPhotometry['Frequency'][iInfrared])
        fInfrared_Reference = QueryPhotometry['Refcode'][iInfrared]
        
        # 
        # print
        for i in range(len(iInfrared)):
            fInfrared_FluxError[i] = fInfrared_FluxError[i].replace('+/-','').replace('...','0.0')
            if len(fInfrared_FluxError[i])==0: fInfrared_FluxError[i]=0.0
            if fInfrared_FluxError[i]=='...': fInfrared_FluxError[i]=0.0 #<20180201>#
            fInfrared_FluxError[i] = float(fInfrared_FluxError[i])*1e3 # mJy
            fInfrared_FluxValue[i] = float(fInfrared_FluxValue[i])*1e3 # mJy
            fInfrared_Wavelengthi = 2.99792458e5/(float(fInfrared_Frequency[i])/1e9) # um
            print fInfrared_Wavelengthi,fInfrared_FluxValue[i],fInfrared_FluxError[i],fInfrared_Passband[i],fInfrared_Reference[i]
            TemporaryInfrared = numpy.array([(fInfrared_Wavelengthi,fInfrared_FluxValue[i],fInfrared_FluxError[i],fInfrared_Passband[i],fInfrared_Reference[i])], dtype=OutputInfoInfrared.dtype)
            OutputInfoInfrared = rec.stack_arrays( (OutputInfoInfrared, TemporaryInfrared), usemask=False, asrecarray=True, autoconvert=True )
        
    # 
    # return
    return OutputInfoInfrared



def getInfoIRAS(InputList):
    # 
    # check input type
    if type(InputList) is not list:
        InputList = [InputList]
    # 
    # prepare output array
    OutputInfoIRAS = numpy.recarray((0,), \
                                  dtype=[('f12', float), ('df12', float), ('ref12', '|S200'), \
                                         ('f25', float), ('df25', float), ('ref25', '|S200'), \
                                         ('f60', float), ('df60', float), ('ref60', '|S200'), \
                                         ('f100', float), ('df100', float), ('ref100', '|S200')])
    # 
    # loop each Input 
    for InputItem in InputList:
        QueryObject = Ned.query_object(InputItem)
        QueryPhotometry = Ned.get_table(InputItem, table='photometry') # Must be one of [‘photometry’|’positions’|’diameters’|’redshifts’|’references’|’object_notes’]
        
        # 
        # print InputItem
        print "----------------------------------------------------------------------------------------------------------------"
        print InputItem
        
        # 
        # prepare output array
        OutputFluxIRAS12 = 0.0
        OutputFluxIRAS25 = 0.0
        OutputFluxIRAS60 = 0.0
        OutputFluxIRAS100 = 0.0
        OutputFErrIRAS12 = 0.0
        OutputFErrIRAS25 = 0.0
        OutputFErrIRAS60 = 0.0
        OutputFErrIRAS100 = 0.0
        OutRefcodeIRAS12 = ""
        OutRefcodeIRAS25 = ""
        OutRefcodeIRAS60 = ""
        OutRefcodeIRAS100 = ""
        
        # 
        # extract IRAS fluxes
        for i in range(len(QueryPhotometry)):
            if QueryPhotometry['Observed Passband'][i].startswith("IRAS "):
                print "LineFlux", QueryPhotometry['Observed Passband'][i], QueryPhotometry['Photometry Measurement'][i], QueryPhotometry['Uncertainty'][i], QueryPhotometry['Units'][i], QueryPhotometry['Refcode'][i]
                # OutputInfoIRAS.append(QueryPhotometry['Photometry Measurement'][i])
                # print ""
                # print str(QueryPhotometry['NED Photometry Measurement'][i])
                if str(QueryPhotometry['NED Photometry Measurement'][i])=="--":
                    # if the IRAS flux is an upper limit
                    if QueryPhotometry['Observed Passband'][i].startswith("IRAS 12 microns"):
                        if (OutputFluxIRAS12 == 0.0 and OutputFErrIRAS12 == 0.0) or QueryPhotometry['Refcode'][i].startswith("2003AJ....126.1607S"):
                            OutputFluxIRAS12 = 0.0                                                              # mJy
                            OutputFErrIRAS12 = float(QueryPhotometry['NED Uncertainty'][i].replace('<',''))/3.0 # mJy
                            OutRefcodeIRAS12 = QueryPhotometry['Refcode'][i]
                    if QueryPhotometry['Observed Passband'][i].startswith("IRAS 25 microns"):
                        if (OutputFluxIRAS12 == 0.0 and OutputFErrIRAS12 == 0.0) or QueryPhotometry['Refcode'][i].startswith("2003AJ....126.1607S"):
                            OutputFluxIRAS25 = 0.0                                                              # mJy
                            OutputFErrIRAS25 = float(QueryPhotometry['NED Uncertainty'][i].replace('<',''))/3.0 # mJy
                            OutRefcodeIRAS25 = QueryPhotometry['Refcode'][i]
                    if QueryPhotometry['Observed Passband'][i].startswith("IRAS 60 microns"):
                        if (OutputFluxIRAS12 == 0.0 and OutputFErrIRAS12 == 0.0) or QueryPhotometry['Refcode'][i].startswith("2003AJ....126.1607S"):
                            OutputFluxIRAS60 = 0.0                                                              # mJy
                            OutputFErrIRAS60 = float(QueryPhotometry['NED Uncertainty'][i].replace('<',''))/3.0 # mJy
                            OutRefcodeIRAS60 = QueryPhotometry['Refcode'][i]
                    if QueryPhotometry['Observed Passband'][i].startswith("IRAS 100 microns"):
                        if (OutputFluxIRAS12 == 0.0 and OutputFErrIRAS12 == 0.0) or QueryPhotometry['Refcode'][i].startswith("2003AJ....126.1607S"):
                            OutputFluxIRAS100 = 0.0                                                              # mJy
                            OutputFErrIRAS100 = float(QueryPhotometry['NED Uncertainty'][i].replace('<',''))/3.0 # mJy
                            OutRefcodeIRAS100 = QueryPhotometry['Refcode'][i]
                else:
                    if QueryPhotometry['Observed Passband'][i].startswith("IRAS 12 microns"):
                        if (OutputFluxIRAS12 == 0.0 and OutputFErrIRAS12 == 0.0) or QueryPhotometry['Refcode'][i].startswith("2003AJ....126.1607S"):
                            OutputFluxIRAS12 = float(QueryPhotometry['NED Photometry Measurement'][i])         # mJy
                            OutputFErrIRAS12 = float(QueryPhotometry['NED Uncertainty'][i].replace('+/-','').replace('...','0.0'))  # mJy
                            OutRefcodeIRAS12 = QueryPhotometry['Refcode'][i]
                    if QueryPhotometry['Observed Passband'][i].startswith("IRAS 25 microns"):
                        if (OutputFluxIRAS12 == 0.0 and OutputFErrIRAS12 == 0.0) or QueryPhotometry['Refcode'][i].startswith("2003AJ....126.1607S"):
                            OutputFluxIRAS25 = float(QueryPhotometry['NED Photometry Measurement'][i])         # mJy
                            OutputFErrIRAS25 = float(QueryPhotometry['NED Uncertainty'][i].replace('+/-','').replace('...','0.0'))  # mJy
                            OutRefcodeIRAS25 = QueryPhotometry['Refcode'][i]
                    if QueryPhotometry['Observed Passband'][i].startswith("IRAS 60 microns"):
                        if (OutputFluxIRAS12 == 0.0 and OutputFErrIRAS12 == 0.0) or QueryPhotometry['Refcode'][i].startswith("2003AJ....126.1607S"):
                            OutputFluxIRAS60 = float(QueryPhotometry['NED Photometry Measurement'][i])         # mJy
                            OutputFErrIRAS60 = float(QueryPhotometry['NED Uncertainty'][i].replace('+/-','').replace('...','0.0'))  # mJy
                            OutRefcodeIRAS60 = QueryPhotometry['Refcode'][i]
                    if QueryPhotometry['Observed Passband'][i].startswith("IRAS 100 microns"):
                        if (OutputFluxIRAS12 == 0.0 and OutputFErrIRAS12 == 0.0) or QueryPhotometry['Refcode'][i].startswith("2003AJ....126.1607S"):
                            OutputFluxIRAS100 = float(QueryPhotometry['NED Photometry Measurement'][i])        # mJy
                            OutputFErrIRAS100 = float(QueryPhotometry['NED Uncertainty'][i].replace('+/-','').replace('...','0.0')) # mJy
                            OutRefcodeIRAS100 = QueryPhotometry['Refcode'][i]
        
        # 
        # Combine for output
        TemporaryIRAS = numpy.array([(OutputFluxIRAS12,OutputFErrIRAS12,OutRefcodeIRAS12,OutputFluxIRAS25,OutputFErrIRAS25,OutRefcodeIRAS25,OutputFluxIRAS60,OutputFErrIRAS60,OutRefcodeIRAS60,OutputFluxIRAS100,OutputFErrIRAS100,OutRefcodeIRAS100)], dtype=OutputInfoIRAS.dtype)
        OutputInfoIRAS = rec.stack_arrays( (OutputInfoIRAS, TemporaryIRAS), usemask=False, asrecarray=True, autoconvert=True )
    
    # 
    # return
    return OutputInfoIRAS



def getInfoPHOT(InputList):
    # 
    # check input type
    if type(InputList) is not list:
        InputList = [InputList]
    # 
    # prepare output array
    OutputInfoPHOT = []
    # 
    # loop each Input 
    for InputItem in InputList:
        QueryObject = Ned.query_object(InputItem)
        QueryPhotometry = Ned.get_table(InputItem, table='photometry') # Must be one of [‘photometry’|’positions’|’diameters’|’redshifts’|’references’|’object_notes’]
        
        # 
        # print InputItem
        print "--------------------------------------------------------"
        print InputItem
        
        # 
        # extract PHOT fluxes
        print "#", "%13s"%("Wave"), "%15s"%("Flux"), "%15s"%("FErr"), "%30s"%("Passband"), "%30s"%("Refcode")
        print "#", "%13s"%("um"), "%15s"%("Jy"), "%15s"%("Jy"), "%30s"%(" "), "%30s"%(" ")
        for i in range(len(QueryPhotometry)):
            if QueryPhotometry['Frequency Mode'][i].startswith("Broad-band "):
                
                # remove uncertainty "+/-" string
                # remove flux "--" string upper limit
                if len(QueryPhotometry['NED Uncertainty'][i]) == 0:
                    QueryPhotometry['NED Uncertainty'][i] = "0"
                elif QueryPhotometry['NED Uncertainty'][i].startswith('+/-'):
                    QueryPhotometry['NED Uncertainty'][i] = QueryPhotometry['NED Uncertainty'][i].replace('+/-','').replace('...','0.0')
                else:
                    QueryPhotometry['NED Photometry Measurement'][i] = 0.0
                
                # 
                QueryFlux = float(QueryPhotometry['NED Photometry Measurement'][i])
                QueryFErr = float(QueryPhotometry['NED Uncertainty'][i].replace('+/-','').replace('...','0.0'))
                QueryWave = 2.99792458e5/(float(QueryPhotometry['Frequency'][i])/1e9)
                QueryBand = QueryPhotometry['Observed Passband'][i].replace(' ','')
                
                # convert uncertainty "<" upper limit to 1-sigma
                if QueryPhotometry['NED Photometry Measurement'][i] == 0.0:
                    QueryFErr = QueryFErr/3.0
                
                # print results
                # print "%15g"%(QueryWave), "%15g"%(QueryFlux), "%15g"%(QueryFErr), "%15s"%(QueryPhotometry['Photometry Measurement'][i]), "%-15s"%(QueryPhotometry['Uncertainty'][i]), "%-15s"%(QueryPhotometry['Units'][i]), "%30s"%(QueryPhotometry['Observed Passband'][i])
                print "%15g"%(QueryWave), "%15g"%(QueryFlux), "%15g"%(QueryFErr), "%30s"%(QueryBand), "%30s"%(QueryPhotometry['Refcode'][i])
                OutputInfoPHOT.append(QueryPhotometry['Photometry Measurement'][i])
        
    # 
    # return
    return OutputInfoPHOT

















# Examples of astroquery NED class
# https://keflavich-astroquery.readthedocs.org/en/doc_frontpage/_generated/astroquery.ned.core.query_ned_photometry.html
# --------------------------------------------------------
# |                       Name | Unit |    Type | Format |
# --------------------------------------------------------
# |                        No. | None |   int32 |    12i |
# |          Observed Passband | None |    |S20 |    20s |
# |     Photometry Measurement | None | float64 | 25.17e |
# |                Uncertainty | None |    |S11 |    11s |
# |                      Units | None |    |S20 |    20s |
# |                  Frequency |   Hz | float64 | 25.17e |
# | NED Photometry Measurement |   Jy | float64 | 25.17e |
# |            NED Uncertainty | None |    |S11 |    11s |
# |                  NED Units | None |     |S2 |     2s |
# |                    Refcode | None |    |S19 |    19s |
# |               Significance | None |    |S23 |    23s |
# |        Published frequency | None |    |S17 |    17s |
# |             Frequency Mode | None |    |S71 |    71s |
# |       Coordinates Targeted | None |    |S31 |    31s |
# |               Spatial Mode | None |    |S24 |    24s |
# |                 Qualifiers | None |    |S40 |    40s |
# |                   Comments | None |   |S161 |   161s |
# --------------------------------------------------------



















############################################################################
# 
############################################################################

try:
    import os,sys
except ImportError: 
    print "Error! Could not import os,sys!"

try:
    import glob
except ImportError: 
    print "Error! Could not import glob!"

try:
    import re
except ImportError: 
    print "Error! Could not import re!"

try:
    import string
except ImportError: 
    print "Error! Could not import string!"

try:
    import math
except ImportError: 
    print "Error! Could not import math!"

try:
    import numpy
except ImportError: 
    print "Error! Could not import numpy!"

try:
    import numpy.lib.recfunctions as rec
except ImportError: 
    print "Error! Could not import numpy.lib.recfunctions!"

try:
    import astropy
except ImportError: 
    print "Error! Could not import astropy!"

try:
    import astropy.io.ascii as asciitable
except ImportError: 
    print "Error! Could not import astropy.io.ascii as asciitable!"

try:
    import astropy.io.fits as fitsfile
except ImportError: 
    print "Error! Could not import astropy.io.fits as fitsfile!"

try:
    import astroquery
except ImportError: 
    print "Error! Could not import astroquery!"

try:
    from astroquery.ned import Ned
except ImportError: 
    print "Error! Could not from astroquery.ned import Ned!"











