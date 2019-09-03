#!/usr/bin/env fish
# 


set BIN_SETUP_SCRIPT (dirname (status --current-filename))/Software/bin_setup.bash

#echo 
#echo "PATH = $PATH"
#echo 

set -x PATH (string split ":" (bash -c "source $BIN_SETUP_SCRIPT -debug -check pdbi-uvt-go-average pdbi-uvt-go-uvfit -print" | tail -n 1))

type FTS_tool_convolve_fits_image
type FTS_tool_draw_pointing_pattern

#echo 
#echo "PATH = $PATH"
#echo 


