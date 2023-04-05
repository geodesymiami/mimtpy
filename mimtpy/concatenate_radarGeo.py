#!/usr/bin/env python3
#################################################################
# Program is used for concatenating data under radar geometry   #
# Author: Lv Xiaoran                                            #
# Created: Mar 2023                                             #
#################################################################

import os
import argparse
import numpy as np
import h5py

import mintpy
import mintpy.workflow
from mintpy.utils import readfile, ptime, writefile, utils as ut
######################################################################################
NOTE = """
Please Note:
1. This script concatenates two datasets together. The input files are displacement (velocity/timeseries) datasets, their corresponding geometryRadar.h5 files
2. Regarding to the reference point, now the script only support the specified point. Given the specified point, the script adjust the displacement by using the value of a point which is the nearest point to the specified point.
3. If a batch concatenation needed, please use the concatnate_patches.py script.

"""

EXAMPLE = """example:

    concatenate_radarGeo.py miaplpy_NE/velocity_msk.h5  miaplpy_NNE/velocity_msk.h5 --geo_file1 miaplpy_NE/inputs/geometryRadar.h5 --geo_file2 miaplpy_NNE/intpus/geometryRadar.h5 --ref_ll 39.04 118.45 --datatype velocity --output velocity_NE_NNE --outdir miaplpyBig

"""

def create_parser():
    parser = argparse.ArgumentParser(description='Concatenate miaplpy patches',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=NOTE + '\n' + EXAMPLE)

    parser.add_argument('patch_files', nargs='+', type=str, help='two displacement datasets to be concatenated \n')

    parser.add_argument('--geo_file1', nargs=1, type=str, help='geometryRadar file of dataset1. \n')
    
    parser.add_argument('--geo_file2', nargs=1, type=str, help='geometryRadar file of dataset2. \n')
    
    parser.add_argument('--ref_ll', nargs=2, type=float, help='reference point. \n')
    
    parser.add_argument('--datatype', nargs=1, type=str, help='data type. \n')
    
    parser.add_argument('--output', nargs=1, type=str, help='output name of concatenated data. \n')

    parser.add_argument('--outdir',dest='outdir',nargs=1,
                        help='output directory to store the concatenated results')
    
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps

def find_minelement_index(lat_flatten, lon_flatten, ref_lat, ref_lon):
    abslat = np.abs(lat_flatten - ref_lat)
    abslon = np.abs(lon_flatten - ref_lon)
    
    # distance
    dis = (abslat ** 2 + abslon ** 2)
    loc = np.where(dis == np.nanmin(dis))

    return loc    

def concatenate_vel(inps, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten):
    data1 = inps.patch_files[0]
    data2 = inps.patch_files[1]
    
    ref_lat = inps.ref_ll[0]
    ref_lon = inps.ref_ll[1]

    vel1, vel1_atr = readfile.read(data1, datasetName='velocity')
    vel2, vel2_atr = readfile.read(data2, datasetName='velocity')

    vel1_flatten = vel1.flatten('F')
    loc = find_minelement_index(lat1_flatten, lon1_flatten, ref_lat, ref_lon)
    # change vel value
    if np.isnan(vel1_flatten[loc[0][0]]):
        raise ValueError('The value at the reference point is nan. Please change reference point!')
    else:
        vel1_flatten -= vel1_flatten[loc[0][0]]

    vel2_flatten = vel2.flatten('F')
    loc = find_minelement_index(lat2_flatten, lon2_flatten, ref_lat, ref_lon)
    # change vel value
    if np.isnan(vel2_flatten[loc[0][0]]):
        raise ValueError('The value at the reference point is nan. Please change reference point!')
    else:
        vel2_flatten -= vel2_flatten[loc[0][0]]

    vel_joined = np.transpose(np.array([np.hstack((vel1_flatten, vel2_flatten))]))

    return vel_joined

def concatenate_ts(inps, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten):
    data1 = inps.patch_files[0]
    data2 = inps.patch_files[1]

    ref_lat = inps.ref_ll[0]
    ref_lon = inps.ref_ll[1]

    ts1, ts1_atr = readfile.read(data1, datasetName='timeseries')
    ts2, ts2_atr = readfile.read(data2, datasetName='timeseries')

    bperp_date = h5py.File(data1,'r')
    bperp = bperp_date['/bperp']
    date = np.array(bperp_date['/date']).tolist()

    ndis = ts1.shape[0]
    pnum = len(lat1_flatten) + len(lat2_flatten)
    ts_join = np.empty(shape=(ndis, pnum, 1), dtype=np.float)

    loc1 = find_minelement_index(lat1_flatten, lon1_flatten, ref_lat, ref_lon)
    loc2 = find_minelement_index(lat2_flatten, lon2_flatten, ref_lat, ref_lon)

    for i in np.arange(ndis):
        ts1_tmp = ts1[i, :, :].flatten()
        ts2_tmp = ts2[i, :, :].flatten()

        # change vel value
        if np.isnan(ts1_tmp[loc1[0][0]]):
            raise ValueError('The value at the reference point is nan. Please change reference point!')
        else:
            ts1_tmp -= ts1_tmp[loc1[0][0]]

        # change vel value
        if np.isnan(ts2_tmp[loc2[0][0]]):
            raise ValueError('The value at the reference point is nan. Please change reference point!')
        else:
            ts2_tmp -= ts2_tmp[loc2[0][0]]
        ts_tmp_joined = np.transpose(np.array([np.hstack((ts1_tmp, ts2_tmp))]))

        ts_join[i, :, :] = ts_tmp_joined

    return ts_join, ts1_atr, bperp, date
        
def concatenate_geo(inps):
    """concatenate geometry data"""
    data_geo1 = inps.geo_file1[0]
    data_geo2 = inps.geo_file2[0]

    lat1, lat_atr1 = readfile.read(data_geo1, datasetName='latitude')
    lon1, lon_atr1 = readfile.read(data_geo1, datasetName='longitude')
    inc1, inc_atr1 = readfile.read(data_geo1, datasetName='incidenceAngle')
    azi1, azi_atr1 = readfile.read(data_geo1, datasetName='azimuthAngle')
    hgt1, hgt_atr1 = readfile.read(data_geo1, datasetName='height')

    lat2, lat_atr2 = readfile.read(data_geo2, datasetName='latitude')
    lon2, lon_atr2 = readfile.read(data_geo2, datasetName='longitude')
    inc2, inc_atr2 = readfile.read(data_geo2, datasetName='incidenceAngle')
    azi2, azi_atr2 = readfile.read(data_geo2, datasetName='azimuthAngle')
    hgt2, hgt_atr2 = readfile.read(data_geo2, datasetName='height')
    
    lat1_flatten = lat1.flatten('F') # flatten matrix according colms
    lon1_flatten = lon1.flatten('F')
    inc1_flatten = inc1.flatten('F')
    azi1_flatten = azi1.flatten('F')
    hgt1_flatten = hgt1.flatten('F')

    lat2_flatten = lat2.flatten('F')
    lon2_flatten = lon2.flatten('F')
    inc2_flatten = inc2.flatten('F')
    azi2_flatten = azi2.flatten('F')
    hgt2_flatten = hgt2.flatten('F')
    
    # do concatenation
    lat_joined = np.transpose(np.array([np.hstack((lat1_flatten, lat2_flatten))]))
    lon_joined = np.transpose(np.array([np.hstack((lon1_flatten, lon2_flatten))]))
    inc_joined = np.transpose(np.array([np.hstack((inc1_flatten, inc2_flatten))]))
    azi_joined = np.transpose(np.array([np.hstack((azi1_flatten, azi2_flatten))]))
    hgt_joined = np.transpose(np.array([np.hstack((hgt1_flatten, hgt2_flatten))]))

    return lat_joined, lon_joined, inc_joined, azi_joined, hgt_joined, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten

def write_vel(vel_joined, inps):
    
    row, colm = vel_joined.shape
    ref_lat = inps.ref_ll[0]
    ref_lon = inps.ref_ll[1]

    atr_vel = dict()
    atr_vel['WIDTH'] = str(colm)
    atr_vel['LENGTH'] = str(row)
    atr_vel['REF_LAT'] = str(ref_lat)
    atr_vel['REF_LON'] = str(ref_lon)
    atr_vel['FILE_TYPE'] = 'velocity'
    
    vel_data = dict()
    vel_data['velocity'] = vel_joined

    output_dir = inps.outdir[0]
    data_outname = inps.output[0]
    vel_filename = output_dir + data_outname + '.h5'
    writefile.write(datasetDict=vel_data, out_file=vel_filename, metadata=atr_vel)

    return 

def write_ts(ts_joined, ts_atr, bperp, date, inps):
    row, colm = ts_joined.shape[1: ]
    ref_lat = inps.ref_ll[0]
    ref_lon = inps.ref_ll[1]

    atr_ts = ts_atr
    atr_ts['WIDTH'] = str(colm)
    atr_ts['LENGTH'] = str(row)
    atr_ts['REF_LAT'] = str(ref_lat)
    atr_ts['REF_LON'] = str(ref_lon)
    atr_ts['FILE_TYPE'] = 'timeseries'
    
    ts_data = dict()
    ts_data['timeseries'] = ts_joined
    
    ts_bp = np.array(bperp[:],dtype=np.float64)
    ts_date= np.array(date, dtype=np.string_)
    ts_data['bperp'] = ts_bp
    ts_data['date'] = ts_date

    output_dir = inps.outdir[0]
    data_outname = inps.output[0]
    ts_filename = output_dir + data_outname + '.h5'
    writefile.write(datasetDict=ts_data, out_file=ts_filename, metadata=atr_ts)
    
    return

def write_geo(lat_joined, lon_joined, inc_joined, azi_joined, hgt_joined, inps):

    row, colm = lat_joined.shape
    ref_lat = inps.ref_ll[0]
    ref_lon = inps.ref_ll[1]

    # write simple attribution
    atr_geo = dict()
    atr_geo['WIDTH'] = str(colm)
    atr_geo['LENGTH'] = str(row)
    atr_geo['FILE_TYPE'] = 'geometry'


    lat_data = dict()
    lat_data['latitude'] = lat_joined

    lon_data = dict()
    lon_data['longitude'] = lon_joined

    geo_data = dict()
    geo_data['azimuthAngle'] = azi_joined
    geo_data['incidenceAngle'] = inc_joined
    geo_data['height'] = hgt_joined
    geo_data['latitude'] = lat_joined
    geo_data['longitude'] = lon_joined

    output_dir = inps.outdir[0]
    data_outname = inps.output[0]
    geo_outname = data_outname.replace(data_outname.split('_')[1], 'geometry')
    suffix = "_".join(data_outname.split('_')[2:]) 
    
    lat_filename = output_dir + '/latitude_' + suffix + '.h5'
    lon_filename = output_dir + '/longitude_' + suffix + '.h5'

    geo_filename = output_dir + geo_outname + '.h5'

    # write h5 file
    writefile.write(datasetDict=lat_data, out_file=lat_filename, metadata=atr_geo)
    writefile.write(datasetDict=lon_data, out_file=lon_filename, metadata=atr_geo)
    writefile.write(datasetDict=geo_data, out_file=geo_filename, metadata=atr_geo)

    print('Finish!')

def main(iargs=None):
    inps = cmd_line_parse(iargs)   
    
    print('process the geometry info')
    lat_joined, lon_joined, inc_joined, azi_joined, hgt_joined, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten = concatenate_geo(inps)
    write_geo(lat_joined, lon_joined, inc_joined, azi_joined, hgt_joined, inps)

    print('process %s data' % inps.datatype[0])
    if inps.datatype[0] == 'velocity':
        vel_joined = concatenate_vel(inps, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten)
        write_vel(vel_joined, inps)
    elif inps.datatype[0] == 'timeseries':
        ts_joined, ts_atr, bperp, date = concatenate_ts(inps, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten)
        write_ts(ts_joined, ts_atr, bperp, date, inps)

######################################################################################
if __name__ == '__main__':
    main()
