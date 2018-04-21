#!/usr/bin/env python
# 
# Query the astroquery for galaxy radius R25





# Import python packages

import os, sys, time, numpy, astropy, astroquery

from astroquery.vizier import Vizier




# Read user input

input_sources = []

i = 1
while i < len(sys.argv):
    str_argv = sys.argv[i].lower().replace('--','-').replace('_','-')
    if str_argv.startswith('-'):
        pass
    else:
        input_sources.append(sys.argv[i])
    i = i+1

if len(input_sources) == 0:
    print('Usage: ')
    print('    Galaxy_info_tool_query_radius_R25.py')
    print('')



# Query astroquery

catalog_list = Vizier.find_catalogs('de Vaucouleurs+ 1991')
print({k:v.description for k,v in catalog_list.items()})

for i in range(len(input_sources)):
    result = Vizier.query_object(input_sources[i], catalog='de Vaucouleurs+ 1991')
    print(result)




