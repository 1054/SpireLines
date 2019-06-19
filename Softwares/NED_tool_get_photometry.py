#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 

import os, sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import CrabNedQuery


if len(sys.argv) <= 1:
    print('Usage: ')
    print('  NED_tool_get_photometry.py UGC5101 > extracted_flux_from_NED_for_UGC05101.txt')
    print('')
    sys.exit()


f = CrabNedQuery.getInfoPHOT(sys.argv[1:])



