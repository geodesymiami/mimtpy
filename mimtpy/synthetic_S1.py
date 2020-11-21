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
import h5py

import mintpy
from mintpy.utils import readfile, writefile, utils as ut
from mimtpy.utils import multitrack_utilities as mut
######################################################################################
EXAMPLE = """example:
    
    synthetic_S1.py mask.h5 BogdModelData_33.h5 --ramp --subswath --lls1 43.60 101.86 46.38 102.47 --lls2 43.53 102.94 46.47 103.67 --tramp 1e-06 1e-06 0.003 --tiramp 4e-07 4e-07 0.002 6e-07 6e-07 0.005 8e-07 8e-07 0.008 --wrap --outdir ./4ramps/

   synthetic_S1.py mask.h5 BogdModelData_106.h5 --ramp --subswath --lls1 43.56 99.80 46.34 100.42 --lls2 43.65 100.91 46.42 101.59 --tramp 3e-07 3e-07 0.002 --tiramp 1e-07 1e-07 0.001 3e-07 3e-07 0.004 6e-07 6e-07 0.008--wrap --outdir ./4ramps/

    synthetic_S1.py mask.h5 BogdModelData_4.h5 --ramp --subswath --lls1 43.62 97.77 46.38 98.37 --lls2 43.54 98.85 46.48 99.57 --tramp 3e-08 3e-08 0.001 --tiramp 1e-08 1e-08 0.001 2e-08 2e-08 0.001 5e-08 5e-08 0.001 --wrap --outdir ./4ramps/

    synthetic_S1.py mask.h5 KokoxiliModelData_AT143.h5 --ramp KokoxiliTramp_AT143.h5 --subswath --lls1 43.62 97.77 46.38 98.37 --lls2 43.54 98.85 46.48 99.57 --tramp 3e-08 3e-08 --tiramp 1e-08 1e-08 2e-08 2e-08 5e-08 5e-08 --outdir ./4ramps/ 
    
    synthetic_S1.py mask.h5 --ramp --subswath --lls1 43.62 97.77 46.38 98.37 --lls2 43.54 98.85 46.48 99.57 --tramp 3e-08 3e-08 --tiramp 1e-08 1e-08 2e-08 2e-08 5e-08 5e-08 --outdir ./4ramps/ 

    synthetic_S1.py mask.h5 KokoxiliModelData_AT143.h5 --outdir ./
    
    synthetic_S1.py mask.h5 KokoxiliModelData_AT143.h5 --ramp --tramp 1e-06 1e-06 0.003 --tiramp_whole 1e-07 1e-08 0.001 --outdir ./
    
    synthetic_S1.py mask.h5 KokoxiliModelData_AT143.h5 --ramp --tramp 1e-06 1e-06 0.003 --tiramp_whole 1e-07 1e-08 0.001 --atmo --atmofile atmospheric_noise.h5 --outdir ./
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Prepare data for Kite software',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs=1, type=str, help='geocoded mask file\n')

    parser.add_argument('modelfile', nargs='?', type=str, help='coseismic model field.')

    parser.add_argument('--ramp', action='store_true', default=False, help='whether simulate tectonic and atmospheric ramp')

    parser.add_argument('trampfile', nargs='?', type=str, help='tectonic ramp file.')

    parser.add_argument('--tramp', dest='tramp', type=float, nargs=3, metavar=('txramp','tyramp','toffset'),
                         help='tectonic linear ramp in X direction and Y direction.')
   
    parser.add_argument('--tiramp_whole', dest='tiramp_whole', type=float, nargs=3, metavar=('tixramp','tiyramp','tioffset'),
                         help='troposphere and ionosphere linear ramp in X direction and Y direction.')
 
    parser.add_argument('--subswath', action='store_true', default=False, help='whether simulate troposphere/ionosphere ramp for each subswath')
    
    parser.add_argument('--lls1', dest='latlons1', type=float, nargs=4, metavar=('LowLat', 'LowLon', 'HighLat', 'HighLon'),
                        help='lat and lon for the two points located at the intersect corner of swath 1 and swath 2.')
    
    parser.add_argument('--lls2', dest='latlons2', type=float, nargs=4, metavar=('LowLat', 'LowLon', 'HighLat', 'HighLon'),
                        help='lat and lon for the two points located at the intersect corner of swath 2 and swath 3.')

    #parser.add_argument('-g', '--geometryRadar', dest='geometry', type=str, nargs=1,
    #                    help='geometry file')
    
    parser.add_argument('--tiramp', dest='tiramp', type=float, nargs=9, metavar=('s1xramp', 's1yramp', 's1offset', 's2xramp', 's2yramp', 's2offset', 's3xramp', 's3yramp', 's3offset'),
                         help='residual troposphere / ionosphere linear ramp for swath1, swath2 and swath3 in X and Y direction.')

    parser.add_argument('--wrap', action='store_true', default=False, help='whether wrap the synthetic data.\n')

    parser.add_argument('--atmo', action='store_true', default=False, help='whether add atmo noise data.')
  
    parser.add_argument('--atmofile', dest='atmofile', nargs=1, type=str, help='atmospheric noise')

    #parser.add_argument('-o','--outfile',dest='outfile',nargs=1, type=str,
    #                    help='outfile name')
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
    #inps.inc_angle = readfile.read(inps.geometry[0], datasetName='incidenceAngle')[0]
    #inps.head_angle = np.ones(inps.inc_angle.shape, dtype=np.float32) * float(inps.metadata['HEADING'])
    #inps.height = readfile.read(inps.geometry[0], datasetName='height')[0]

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
    row, colm =  original_mask.shape
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

def read_model_h5(files):
    """read h5 file generated by matlab"""
    model_file = files
    model_idx = h5py.File(model_file,'r')
    data = model_idx['/model'][()]
    model_data = np.transpose(data)

    return model_data

def generate_linear_ramp(mask_swaths, X, Y, inps):
    """generate linear ramp for the synthetic data"""
    # ramp data
    inps.row, inps.colm = mask_swaths.shape
    data_ramp = np.zeros((inps.row, inps.colm), dtype=float)
    data_ramp = data_ramp.reshape(-1,1)
    
    # flatten mask_swaths
    mask_swaths_flat = mask_swaths.reshape(-1,1)
     
    # for tectonic linear ramp
    if inps.trampfile:
        data_tramp = read_model_h5(inps.trampfile)
        data_tramp = data_tramp.reshape(-1, 1)
    else:
        txramp = float(inps.tramp[0])
        tyramp = float(inps.tramp[1])
        toffset = float(inps.tramp[2]) 
 
        # calculate tectonic linear ramp
        data_tramp = (X * (mask_swaths_flat != 0)) * txramp + (Y * (mask_swaths_flat != 0)) * tyramp + toffset
    
    data_ramp += data_tramp

    # for residual troposphere / ionosphere linear ramp
    if inps.subswath:
        s1xramp = float(inps.tiramp[0])
        s1yramp = float(inps.tiramp[1])
        s1offset = float(inps.tiramp[2])
        s2xramp = float(inps.tiramp[3])
        s2yramp = float(inps.tiramp[4])
        s2offset = float(inps.tiramp[5])
        s3xramp = float(inps.tiramp[6])
        s3yramp = float(inps.tiramp[7])
        s3offset = float(inps.tiramp[8])

        # calculate residual troposphere / ionosphere linear ramp for subswath 1/2/3
        data_tiramp = np.zeros((inps.row, inps.colm), dtype=float).reshape(-1,1)

        data_tiramp += ((X * (mask_swaths_flat == 1)) * s1xramp + (Y * (mask_swaths_flat == 1)) * s1yramp + s1offset)
        data_tiramp += ((X * (mask_swaths_flat == 2)) * s2xramp + (Y * (mask_swaths_flat == 2)) * s2yramp + s2offset)
        data_tiramp += ((X * (mask_swaths_flat == 3)) * s3xramp + (Y * (mask_swaths_flat == 3)) * s3yramp + s3offset)

    else:
        tixramp = float(inps.tiramp_whole[0])
        tiyramp = float(inps.tiramp_whole[1])
        tioffset = float(inps.tiramp_whole[2]) 
 
        # calculate tectonic linear ramp
        data_tiramp = (X * (mask_swaths_flat != 0)) * tixramp + (Y * (mask_swaths_flat != 0)) * tiyramp + tioffset
        
    data_ramp += data_tiramp

    # read the atmospheric data
    if inps.atmo:
        inps.atmo_data = read_model_h5(inps.atmofile[0])
        inps.atmo_data = inps.atmo_data[0:inps.row, 0:inps.colm]
        write_atmo_data(inps.atmo_data, inps)

    # calculated the synthetic model data + ramp
    if inps.modelfile:
        inps.model_data = read_model_h5(inps.modelfile)
        model_data_copy = copy.deepcopy(inps.model_data)
        write_model_data(model_data_copy, inps)
        row_model, colm_model = inps.model_data.shape
        if row_model != inps.row or colm_model != inps.colm:
            raise Exception('Error! The dimension of modeled data and simulated ramp is disagree!')
        inps.model_data = inps.model_data.reshape(-1,1)
        inps.model_data += data_ramp
    
    # wrapped
    if inps.wrap:
        vmin = -float(inps.metadata['WAVELENGTH']) / 4
        vmax = float(inps.metadata['WAVELENGTH']) / 4
        data_ramp_wrap = vmin + np.mod(data_ramp - vmin, vmax - vmin)
        data_tramp_wrap = vmin + np.mod(data_tramp - vmin, vmax - vmin)
        data_tiramp_wrap = vmin + np.mod(data_tiramp - vmin, vmax - vmin)
        if inps.modelfile:
            inps.model_data_wrap = vmin + np.mod(inps.model_data - vmin, vmax - vmin)
            write_data(inps, data_tramp_wrap, data_tiramp_wrap, data_ramp_wrap, model_data_ramp='yes', wrap_flag='yes')
        else:
            write_data(inps, data_tramp_wrap, data_tiramp_wrap, data_ramp_wrap, wrap_flag='yes')       

    # range to phase 
    #data_ramp = (data_ramp / float(inps.metadata['WAVELENGTH'])) * (4 * np.pi)
    #data_tramp = (data_tramp / float(inps.metadata['WAVELENGTH'])) * (4 * np.pi)
    #data_tiramp = (data_tiramp / float(inps.metadata['WAVELENGTH'])) * (4 * np.pi)
    if inps.modelfile:
        write_data(inps, data_tramp, data_tiramp, data_ramp, model_data_ramp='yes')
    else:
        write_data(inps, data_tramp, data_tiramp, data_ramp)   

    return 

def write_data(inps, data_tramp, data_tiramp, data_ramp, model_data_ramp=None, wrap_flag=None):
    """write data"""
    # write ramp and ramped data
    inps.metadata['UNIT'] = 'm'
    if wrap_flag:
        inps.metadata['FILE_TYPE'] = 'wrap'
    else:
        inps.metadata['FILE_TYPE'] = '.unw'
 
    data_tramp = data_tramp.reshape(inps.row, inps.colm)
    data_tiramp = data_tiramp.reshape(inps.row, inps.colm)
    data_ramp = data_ramp.reshape(inps.row, inps.colm)
    if model_data_ramp:
        if wrap_flag:
            inps.model_data_wrap = inps.model_data_wrap.reshape(inps.row, inps.colm)
        else:
            inps.model_data = inps.model_data.reshape(inps.row, inps.colm)

    # mask data
    data_ramp[inps.mask == False] = np.nan
    data_tramp[inps.mask == False] = np.nan
    data_tiramp[inps.mask == False] = np.nan
    if model_data_ramp:
        if wrap_flag:
            inps.model_data_wrap[inps.mask == False] = np.nan
        else:
            inps.model_data[inps.mask == False] = np.nan

    if wrap_flag:
        tramp_name = inps.outdir[0] + 'tectonic_ramp_wrap.h5'
        tiramp_name = inps.outdir[0] + 'troposphere_ionosphere_ramp_wrap.h5'
        ramp_name = inps.outdir[0] + 'data_ramped_wrap.h5'
        if model_data_ramp:
            model_ramp_name = inps.outdir[0] + inps.modelfile.split('.')[0] + '_ramped_wrap.h5' 
    else:
        tramp_name = inps.outdir[0] + 'tectonic_ramp.h5'
        tiramp_name = inps.outdir[0] + 'troposphere_ionosphere_ramp.h5'
        ramp_name = inps.outdir[0] + 'data_ramped.h5'
        if model_data_ramp:
            model_ramp_name = inps.outdir[0] + inps.modelfile.split('.')[0] + '_ramped.h5' 
        
    writefile.write(data_tramp, out_file=tramp_name, metadata=inps.metadata)
    writefile.write(data_tiramp, out_file=tiramp_name, metadata=inps.metadata)
    writefile.write(data_ramp, out_file=ramp_name, metadata=inps.metadata)
    if model_data_ramp:
        if wrap_flag:
            writefile.write(inps.model_data_wrap, out_file=model_ramp_name, metadata=inps.metadata)
        else:
            writefile.write(inps.model_data, out_file=model_ramp_name, metadata=inps.metadata)

def write_model_data(model_data_copy, inps):
    """change simulated model data to mintpy format"""
    model_data_copy[inps.mask == False] = np.nan
    outfile = inps.outdir[0] + inps.modelfile.split('.')[0] + '_mintpy.h5'
    writefile.write(model_data_copy, out_file=outfile, metadata=inps.metadata)
    # wrap 
    vmin = -float(inps.metadata['WAVELENGTH']) / 4
    vmax = float(inps.metadata['WAVELENGTH']) / 4
    model_data_copy_wrap = vmin + np.mod(model_data_copy - vmin, vmax - vmin)
    outfile = inps.outdir[0] + inps.modelfile.split('.')[0] + '_mintpy_wrap.h5'
    writefile.write(model_data_copy_wrap, out_file=outfile, metadata=inps.metadata)

def write_atmo_data(atmo_data, inps):
    """change atmospheric noise data to mintpy format"""
    atmo_data[inps.mask == False] = np.nan
    outfile = inps.outdir[0] + inps.atmofile[0].split('.')[0] + '_mintpy.h5'
    writefile.write(atmo_data, out_file=outfile, metadata=inps.metadata)
######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)   
    
    X, Y, origin = ll2xy(inps)
    original_mask = readfile.read(inps.file[0])[0]
    X = np.transpose(np.array([X]))
    Y = np.transpose(np.array([Y]))
    
    if inps.ramp:
        # generate model data/ramp data/model+ramp data
        if inps.tiramp_whole:
            row, colm =  original_mask.shape
            mask_swaths = np.ones((row, colm), dtype=float)
            mask_swaths[inps.mask == False] = np.nan
        else:
            mask_swaths = generate_swaths_mask(inps, original_mask, X, Y, origin)
    
        generate_linear_ramp(mask_swaths, X, Y, inps)
    else:
        # only generate model data
        if inps.modelfile:
            model_data = read_model_h5(inps.modelfile)
            model_data_copy = copy.deepcopy(model_data)
            write_model_data(model_data_copy, inps)
######################################################################################
if __name__ == '__main__':
    main()
