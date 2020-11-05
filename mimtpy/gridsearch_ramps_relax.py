#!/usr/bin/env python3
#################################################################
# Program is used for preparing data for kite software          #
# Author: Lv Xiaoran                                            #
# Created: June 2020                                            #
#################################################################

import os
import argparse
import numpy as np

import mintpy
from mintpy.objects import RAMP_LIST
from mintpy.utils import readfile, writefile

import mimtpy
import mimtpy.workflow
from mimtpy.utils import multitrack_utilities as mu

######################################################################################
EXAMPLE = """example:
  
  gridsearch_ramps_relax.py KokoxiliModelData_AT143_ramped.h5 KokoxiliModel_AT143.h5 mask.h5 --ramp_estimate -s linear -o ramp_noise.h5 --ramp_file ramps.h5 --outdir ./
  
  gridsearch_ramps_relax.py KokoxiliModelData_AT143_ramped.h5 KokoxiliModel_AT143.h5 mask.h5 --outdir ./
  
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Prepare data for Kite software',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('Obs_file', nargs=1, type=str, help='Observed data')

    parser.add_argument('Model_file', nargs=1, type=str, help='Modeled data')

    parser.add_argument('mask', nargs=1, type=str, help='Mask file')

    parser.add_argument('--ramp_estimate', dest='e_ramp', action='store_true', help='estimate ramp or not')

    parser.add_argument('-s', dest='surface_type', default='linear', choices=RAMP_LIST,
                        help='type of surface/ramp to calculated, linear by default')

    parser.add_argument('-o', dest='output', type=str, nargs=1,
                        help='output name for ramp+noise data')
    
    parser.add_argument('--ramp_file',dest='ramp_file',nargs=1,
                        help='estimated ramp data')
    
    parser.add_argument('--outdir', dest='outdir', type=str, nargs=1,
                        help='output file dir')

    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps

def calculate_ramp_noise(inps):
    """subtract forward model data from observed data"""
    observed_file = inps.Obs_file[0]
    mask_file = inps.mask[0]
    ramp_type = inps.surface_type
    ramp_file = inps.ramp_file[0]
    ramp_noise_file = inps.output[0]
    model_file = inps.Model_file[0]

    scp_args = [observed_file, model_file, '-o', ramp_noise_file, '--ramp', '-s', ramp_type, '-m', mask_file, '--ramp_file', ramp_file, '--outdir', inps.outdir[0]]
    scp_args = mu.seperate_str_byspace(scp_args)

    print('subtract_h5.py', scp_args)
    print('ramp_noise file is: %s' % ramp_noise_file)
    print('estimated ramp file is: %s' % ramp_file)
    mimtpy.subtract_h5.main(scp_args.split())

    return ramp_file, ramp_noise_file

def calculate_residual(inps, ramp_noise_file, ramp_file):
    """calculate residual data which is observed data - model data - estimated ramp"""
    residual_file = inps.Model_file[0][:-3] +  '_residual.h5'
    ramp_noise_data = readfile.read(inps.outdir[0] + ramp_noise_file)[0]
    ramp_data, atr = readfile.read(inps.outdir[0] + ramp_file)

    residual = ramp_noise_data - ramp_data

    print('the residual file after subtracting model data and estimated ramp is: %s' % residual_file)
    writefile.write(residual, out_file=inps.outdir[0] + residual_file, metadata=atr)

    return 

    
######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)   
   
    print('--------------------------------------------------\n') 
    
    if inps.e_ramp:
        print('subtract forward model data from observed data')
        ramp_file, ramp_noise_file = calculate_ramp_noise(inps)

        print('--------------------------------------------------\n') 
        print('calculate residual data with estimating ramps')
        calculate_residual(inps, ramp_noise_file, ramp_file)

    else:
        print('calculate residual data without estimating ramp')
        observed_file = inps.Obs_file[0]
        mask_file = inps.mask[0]
        residual_file = inps.Model_file[0][:-3] +  '_residual_noramp.h5'
        model_file = inps.Model_file[0]
        scp_args = [observed_file, model_file, '-o', residual_file, '--outdir', inps.outdir[0]]
        scp_args = mu.seperate_str_byspace(scp_args)

        print('subtract_h5.py', scp_args)
        mimtpy.subtract_h5.main(scp_args.split())

    print('finish!')
######################################################################################
if __name__ == '__main__':
    main()
