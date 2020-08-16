#!/usr/bin/env python3
#################################################################
# Program is used for generate synthetic data using S1 geometry #
# Author: Lv Xiaoran                                            #
# Created: Aug 2020                                             #
#################################################################

import os
import argparse
import numpy as np
import copy

import mintpy
from mintpy.utils import readfile, writefile, utils as ut
from mimtpy.utils import multitrack_utilities as mut
######################################################################################
EXAMPLE = """example:
  
  save_kite.py geo_mask.h5 --lls1 43.6 99.7 46.3 100.5 --lls2 43.7 100.9 46.4 101.7 -g ./inputs/geometryRadar.h5 --tramp 1e-08 1e-08 --tiramp 1e-07 1e-07 5e-07 5e-07 8e-07 8e-07 -o BogdSenDT106_synthetic_linearramp
  
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Prepare data for Kite software',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs=1, type=str, help='geocoded mask file\n')

    parser.add_argument('--lls1', dest='latlons1', type=float, nargs=4, metavar=('LowLat', 'LowLon', 'HighLat', 'HighLon'),
                        help='lat and lon for the two points located at the intersect corner of swath 1 and swath 2.')
    
    parser.add_argument('--lls2', dest='latlons2', type=float, nargs=4, metavar=('LowLat', 'LowLon', 'HighLat', 'HighLon'),
                        help='lat and lon for the two points located at the intersect corner of swath 2 and swath 3.')

    parser.add_argument('-g', '--geometryRadar', dest='geometry', type=str, nargs=1,
                        help='geometry file')
    
    parser.add_argument('--tramp', dest='tramp', type=float, nargs=2, metavar=('txramp','tyramp'),
                         help='tectonic linear ramp in X direction and Y direction.')

    parser.add_argument('--tiramp', dest='tiramp', type=float, nargs=6, metavar=('s1xramp', 's1yramp', 's2xramp', 's2yramp', 's3xramp', 's3yramp'),
                         help='residual troposphere / ionosphere linear ramp for swath1, swath2 and swath3 in X and Y direction.')

    parser.add_argument('-o','--outfile',dest='outfile',nargs=1, type=str,
                        help='outfile name')
    parser.add_argument('--outdir', nargs=1, type=str, help='output dir')    

    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps

def ll2xy(inps):
    """transfer lat/lon to local coordination"""
    inps.metadata = readfile.read_attribute(inps.file[0])

    # read geometry
    inps.lat, inps.lon = ut.get_lat_lon(inps.metadata)
    inps.inc_angle = readfile.read(inps.geometry[0], datasetName='incidenceAngle')[0]
    inps.head_angle = np.ones(inps.inc_angle.shape, dtype=np.float32) * float(inps.metadata['HEADING'])
    inps.height = readfile.read(inps.geometry[0], datasetName='height')[0]

    # read mask file
    inps.mask = readfile.read(inps.file[0])[0]
    # mask data
    #inps.lat[inps.mask==0] = np.nan
    #inps.lon[inps.mask==0] = np.nan
    #inps.inc_angle[inps.mask==0] = np.nan
    #inps.head_angle[inps.mask==0] = np.nan
    #inps.height[inps.mask==0] = np.nan

    # conver latlon to xy
    # origin point
    origin_lat = (inps.lat[0,0] + inps.lat[-1,0]) / 2
    origin_lon = (inps.lon[0,0] + inps.lon[0,-1]) / 2

    lat = np.transpose(inps.lat.reshape(-1,1))
    lon = np.transpose(inps.lon.reshape(-1,1))

    llh = np.vstack((lon, lat))
    origin = np.array([origin_lon, origin_lat], dtype=float)
    XY = np.transpose(mut.llh2xy(llh,origin)) * 1000 # unit of X/Y is meter and is a [N,2] matrix with [N,0] is X; [N,1] is Y
    X = XY[:,0]
    Y = XY[:,1]

    return X, Y, origin

def generate_swaths_mask(inps, original_mask, X, Y, origin):
    """generate mask to distinguish the 3 swaths based on the intersection cover points"""
    # convert lat/lon to local coordination
    X, Y, origin = ll2xy(inps)
    
    # conver lat/lon of cover points to local value
    LowLat_s1 = float(inps.latlons1[0])
    LowLon_s1 = float(inps.latlons1[1])
    HighLat_s1 = float(inps.latlons1[2])
    HighLon_s1 = float(inps.latlons1[3])

    s1 = np.array([[LowLon_s1, HighLon_s1], [LowLat_s1, HighLat_s1]], dtype=float)

    s1_XY = np.transpose(mut.llh2xy(s1, origin)) * 1000
    LowX_s1 = s1_XY[0,0]
    LowY_s1 = s1_XY[0,1]
    HighX_s1 = s1_XY[1,0]
    HighY_s1 = s1_XY[1,1]

    LowLat_s2 = float(inps.latlons2[0])
    LowLon_s2 = float(inps.latlons2[1])
    HighLat_s2 = float(inps.latlons2[2])
    HighLon_s2 = float(inps.latlons2[3])

    s2 = np.array([[LowLon_s2, HighLon_s2], [LowLat_s2, HighLat_s2]], dtype=float)
  
    s2_XY = np.transpose(mut.llh2xy(s2, origin)) * 1000
    LowX_s2 = s2_XY[0,0]
    LowY_s2 = s2_XY[0,1]
    HighX_s2 = s2_XY[1,0]
    HighY_s2 = s2_XY[1,1]

    # generate the indicator to judge swath1 swath2 swath3
    row, colm = original_mask.shape
    mask_s123 = np.zeros((row, colm), dtype=int) * np.nan
    mask_s123_tmp = mask_s123.reshape(-1,1)
  
    swath1_indi = (HighY_s1 - LowY_s1) * X + (LowX_s1 - HighX_s1) * Y + HighX_s1 * LowY_s1 - LowX_s1 * HighY_s1
    swath2_indi = (HighY_s2 - LowY_s2) * X + (LowX_s2 - HighX_s2) * Y + HighX_s2 * LowY_s2 - LowX_s2 * HighY_s2
    
    # for swath1
    # swath1_indi >= 0 
    mask_s123_tmp[swath1_indi < 0] = 1
    # for swath2
    # swath1_indi <0 && swath2_indi >=0
    idx = ((swath1_indi >= 0) * (swath2_indi < 0))
    mask_s123_tmp[idx] = 2
    # for swath3
    # swath2_indi < 0    
    mask_s123_tmp[swath2_indi > 0] = 3

    mask_s123 = mask_s123_tmp.reshape(row, colm)
    mask_swaths = mask_s123 * original_mask

    writefile.write(mask_swaths, out_file=inps.outdir[0] + 'mask_swaths.h5', metadata=inps.metadata)

    return mask_swaths

def generate_linear_ramp(mask_swaths, X, Y, inps):
    """generate linear ramp for the synthetic data"""
    # ramp data
    row, colm = mask_swaths.shape
    data_ramp = np.zeros((row, colm), dtype=float)
    data_ramp = data_ramp.reshape(-1,1)
    
    # flatten mask_swaths
    mask_swaths_flat = mask_swaths.reshape(-1,1)
     
    # for tectonic linear ramp
    txramp = float(inps.tramp[0])
    tyramp = float(inps.tramp[1])
  
    # calculate tectonic linear ramp
    data_tramp = (X * (mask_swaths_flat != 0)) * txramp + (Y * (mask_swaths_flat != 0)) * tyramp
    
    data_ramp += data_tramp

    # for residual troposphere / ionosphere linear ramp
    s1xramp = float(inps.tiramp[0])
    s1yramp = float(inps.tiramp[1])
    s2xramp = float(inps.tiramp[2])
    s2yramp = float(inps.tiramp[3])
    s3xramp = float(inps.tiramp[4])
    s3yramp = float(inps.tiramp[5])

    # calculate residual troposphere / ionosphere linear ramp for swath 1/2/3
    data_tiramp = np.zeros((row, colm), dtype=float).reshape(-1,1)

    data_tiramp += ((X * (mask_swaths_flat == 1)) * s1xramp + (Y * (mask_swaths_flat == 1)) * s1yramp)
    data_tiramp += ((X * (mask_swaths_flat == 2)) * s2xramp + (Y * (mask_swaths_flat == 2)) * s2yramp)
    data_tiramp += ((X * (mask_swaths_flat == 3)) * s3xramp + (Y * (mask_swaths_flat == 3)) * s3yramp)

    data_ramp += data_tiramp

    # wrapped
    vmin = -float(inps.metadata['WAVELENGTH']) / 4
    vmax = float(inps.metadata['WAVELENGTH']) / 4
    data_ramp = vmin + np.mod(data_ramp - vmin, vmax - vmin)
    data_tramp = vmin + np.mod(data_tramp - vmin, vmax - vmin)
    data_tiramp = vmin + np.mod(data_tiramp - vmin, vmax - vmin)
    #data_ramp = (data_ramp / float(inps.metadata['WAVELENGTH'])) * (4 * np.pi)
    #data_tramp = (data_tramp / float(inps.metadata['WAVELENGTH'])) * (4 * np.pi)
    #data_tiramp = (data_tiramp / float(inps.metadata['WAVELENGTH'])) * (4 * np.pi)
    
    # write ramp and ramped data
    inps.metadata['UNIT'] = '.unw'
    inps.metadata['FILE_TYPE'] = 'unw'
    data_tramp = data_tramp.reshape(row, colm)
    data_tiramp = data_tiramp.reshape(row, colm)
    data_ramp = data_ramp.reshape(row, colm)
    
    # change 0 to np.nan
    data_ramp[inps.mask == False] = np.nan
    data_tramp[inps.mask == False] = np.nan
    data_tiramp[inps.mask == False] = np.nan
    
    writefile.write(data_tramp, out_file=inps.outdir[0] + 'tectonic_ramp.h5', metadata=inps.metadata)
    writefile.write(data_tiramp, out_file=inps.outdir[0] + 'troposphere_ionosphere_ramp.h5', metadata=inps.metadata)
    writefile.write(data_ramp, out_file=inps.outdir[0] + 'data_ramped.h5', metadata=inps.metadata)
######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)   
    
    X, Y, origin = ll2xy(inps)
    original_mask = readfile.read(inps.file[0])[0]
    X = np.transpose(np.array([X]))
    Y = np.transpose(np.array([Y]))
    mask_swaths = generate_swaths_mask(inps, original_mask, X, Y, origin)
    generate_linear_ramp(mask_swaths, X, Y, inps)
######################################################################################
if __name__ == '__main__':
    main()
