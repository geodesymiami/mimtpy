#!/usr/bin/env python3
############################################################
# Program is part of MimtPy                                #
# Copyright(c) 2013-2019, Lv Xiaoran                       #
# Author:  Lv Xiaoran                                      #
############################################################

import os
import argparse
import numpy as np
from mintpy.utils import readfile, writefile
from mintpy.objects import RAMP_LIST, deramp

##############################################################################
EXAMPLE = """example:
   subtract_h5.py geo_20180610_20190430_ramp.h5 geo_20181125_20181207_ramp.h5 -b 34 34.6 45.5 46 --output geo_cumulative_des.h5
   subtract_h5.py geo_20180610_20190430_ramp.h5 geo_20181125_20181207_ramp.h5 --output geo_cumulative_des.h5
   subtract_h5.py geo_20180610_20190430_ramp.h5 geo_20181125_20181207_ramp.h5 --output geo_cumulative_des.h5 --ramp -s linear -m mask.h5 --ramp_file ramp.h5
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Subtract the second data from the first data.And calculate ramp from the subtracted data',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs=2, help='HDF5 file to be calculated.')

    parser.add_argument('-b', '--bbox', dest='SNWE', type=float, nargs=4, metavar=('S', 'N', 'W', 'E'),
                        help='Bounding box of area to be geocoded.\n' +
                        'Include the uppler left corner of the first pixel' +
                        '    and the lower right corner of the last pixel')    

    parser.add_argument('-o', '--output', dest='outfile', help='output file name.')

    parser.add_argument('--outdir', dest='outdir', type=str, help='output file dir')

    ramp_option = parser.add_argument_group(title='options for calculate ramps')

    ramp_option.add_argument('--ramp', action='store_true', default=False,
                             help='whether calculate the ramp of subtracted data')
    
    ramp_option.add_argument('-s', dest='surface_type', default='linear', choices=RAMP_LIST,
                             help='type of surface/ramp to remove, linear by default')

    ramp_option.add_argument('-m', dest='mask', type=str, help='mask data')
    
    ramp_option.add_argument('--ramp_file', type=str, help='output name of calculated ramp')
    
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    
    return inps

def latlontoyx(inps):
    '''change SNWD to xy'''
    atr = readfile.read_attribute(inps.file[0])
    row_max = int(np.rint((inps.SNWE[0] - float(atr['Y_FIRST'])) / float(atr['Y_STEP'])))
    row_min = int(np.rint((inps.SNWE[1] - float(atr['Y_FIRST'])) / float(atr['Y_STEP'])))
    
    clom_min = int(np.rint((inps.SNWE[2] - float(atr['X_FIRST'])) / float(atr['X_STEP'])))
    clom_max = int(np.rint((inps.SNWE[3] - float(atr['X_FIRST'])) / float(atr['X_STEP'])))
    
    return row_min, row_max, clom_min, clom_max

def subtract_data(inps):
    # metadata
    atr = readfile.read_attribute(inps.file[0])
    data1 = readfile.read(inps.file[0])[0]
    data2 = readfile.read(inps.file[1])[0]

    #if 'WAVELENGTH' in atr.keys():
    #    phase2range = float(atr['WAVELENGTH'])/(-4 * np.pi)
    
    #print('convert phase to range')
    #data1 *= phase2range
    #data2 *= phase2range
    if inps.SNWE is not None:
        row_min, row_max, clom_min, clom_max = latlontoyx(inps)

        print('extract displacement')
        data1[row_min:row_max,clom_min:clom_max] -= data2[row_min:row_max,clom_min:clom_max]
    else:
        data1 -= data2
    
    # get rid of starting . if output as hdf5 file
    if inps.outfile.endswith('.h5'):
        if atr['FILE_TYPE'].startswith('.'):
            atr['FILE_TYPE'] = atr['FILE_TYPE']

    atr['PROCESSOR'] = 'roipac'

    return data1, atr, inps.outfile

def calculate_ramp(data, atr, inps):
    """calcualte ramp from data file"""
    
    ramp_type = inps.surface_type
    mask = readfile.read(inps.mask)[0]
    ramp = deramp(data, mask, ramp_type=ramp_type, metadata=atr)[1]
    ramp[mask == False] = np.nan

    return ramp
##############################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    data, atr, out_file = subtract_data(inps)

    if inps.ramp:
        ramp = calculate_ramp(data, atr, inps)
        writefile.write(ramp, out_file=inps.outdir + inps.ramp_file, metadata=atr)

    writefile.write(data, out_file=inps.outdir + out_file, metadata=atr)

    return  
##########################################################################
if __name__ == '__main__':
    main()
