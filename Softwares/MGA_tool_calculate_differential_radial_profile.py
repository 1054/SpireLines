#!/usr/bin/env python
# 

import os, sys, astropy
import astropy.io.ascii as asciitable


data_table = asciitable.read(sys.argv[1])
#print(data_table)


radius_pix = data_table['col3']
num_pix = data_table['col4']
sum_pix = data_table['col5']
diff_pix = sum_pix

print('# %13s %25s'%('radius', 'surface_brightness'))

for i in range(len(data_table)):
    if i == 0:
        diff_pix[i] = (sum_pix[i]) / (num_pix[i])
    else:
        diff_pix[i] = (sum_pix[i] - sum_pix[i-1]) / (num_pix[i] - num_pix[i-1])
    # 
    print('%15.5f %25g'%(radius_pix[i], diff_pix[i]))






