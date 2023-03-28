#!/usr/bin/env python3
#################################################################
# Program is used for concatenating data under radar geometry   #
# Author: Lv Xiaoran                                            #
# Created: Mar 2023                                             #
#################################################################

import os
import argparse
import numpy as np

import mintpy
import mintpy.workflow
from mintpy.utils import readfile, ptime, writefile, utils as ut
######################################################################################
EXAMPLE = """example:

    concatenate_radarGeo.py $SCRATCHDIR/TangshanSenAT69/miaplpy_NE_201410_202212/network_delaunay_4/rdr_velocity_msk.h5  $SCRATCHDIR/TangshanSenAT69/miaplpy_NNE_201410_202212/network_delaunay_4/rdr_velocity_msk.h5 --geo_file1 $SCRATCHDIR/TangshanSenAT69/miaplpy_NE_201410_202212/network_delaunay_4/rdr_geometry_msk.h5 --geo_file2 $SCRATCHDIR/TangshanSenAT69/miaplpy_NNE_201410_202212/network_delaunay_4/rdr_geometry_msk.h5 --ref_ll 39.04 118.45 --output rdr_velocity_NE_NNE --outdir miaplpyBig

"""

def create_parser():
    parser = argparse.ArgumentParser(description='Concatenate miaplpy patches',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('patch_files', nargs='+', type=str, help='dataset to be concatenated \n')

    parser.add_argument('--geo_file1', nargs=1, type=str, help='geometryRadar file of dataset1. \n')
    
    parser.add_argument('--geo_file2', nargs=1, type=str, help='geometryRadar file of dataset2. \n')
    
    parser.add_argument('--ref_ll', nargs=2, type=float, help='reference point. \n')
    
    parser.add_argument('--output', nargs=1, type=str, help='output name of concatenated data. \n')

    parser.add_argument('--outdir',dest='outdir',nargs=1,
                        help='output directory to store the concatenated results')
    
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps

def find_minelement_index(lat_flatten, lon_flatten, vel_flatten, ref_lat, ref_lon):
    abslat = np.abs(lat_flatten - ref_lat)
    abslon = np.abs(lon_flatten - ref_lon)
    
    # distance
    dis = (abslat ** 2 + abslon ** 2)
    loc = np.where(dis == np.nanmin(dis))
    
    # change vel value
    if np.isnan(vel_flatten[loc[0][0]]):
        raise ValueError('The value at the reference point is nan. Please change reference point!')
    else:
        vel_flatten -= vel_flatten[loc[0][0]]
    
    return vel_flatten
    #row, colm = array.shape
    #index_1D = (np.absolute(array - value)).argmin()

    #row_index = np.int(index_1D / colm)
    #colm_index = np.int(index_1D % colm)

    #return row_index, colm_index

def concatenate(inps):
    data1 = inps.patch_files[0]
    data2 = inps.patch_files[1]

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
    
    ref_lat = inps.ref_ll[0]
    ref_lon = inps.ref_ll[1]
    
    vel1, vel1_atr = readfile.read(data1, datasetName='velocity')
    vel2, vel2_atr = readfile.read(data2, datasetName='velocity')

    lat1_flatten = lat1.flatten('F') # flatten matrix according colms
    lon1_flatten = lon1.flatten('F')
    inc1_flatten = inc1.flatten('F')
    azi1_flatten = azi1.flatten('F')
    hgt1_flatten = hgt1.flatten('F')
    vel1_flatten = vel1.flatten('F')
    vel1_flatten = find_minelement_index(lat1_flatten, lon1_flatten, vel1_flatten, ref_lat, ref_lon)

    lat2_flatten = lat2.flatten('F')
    lon2_flatten = lon2.flatten('F')
    inc2_flatten = inc2.flatten('F')
    azi2_flatten = azi2.flatten('F')
    hgt2_flatten = hgt2.flatten('F')
    vel2_flatten = vel2.flatten('F')
    vel2_flatten = find_minelement_index(lat2_flatten, lon2_flatten, vel2_flatten, ref_lat, ref_lon)
    
    # do concatenation
    lat_joined = np.transpose(np.array([np.hstack((lat1_flatten, lat2_flatten))]))
    lon_joined = np.transpose(np.array([np.hstack((lon1_flatten, lon2_flatten))]))
    inc_joined = np.transpose(np.array([np.hstack((inc1_flatten, inc2_flatten))]))
    azi_joined = np.transpose(np.array([np.hstack((azi1_flatten, azi2_flatten))]))
    hgt_joined = np.transpose(np.array([np.hstack((hgt1_flatten, hgt2_flatten))]))
    vel_joined = np.transpose(np.array([np.hstack((vel1_flatten, vel2_flatten))]))

    return lat_joined, lon_joined, inc_joined, azi_joined, hgt_joined, vel_joined

def write(lat_joined, lon_joined, inc_joined, azi_joined, hgt_joined, vel_joined, inps):

    row, colm = lat_joined.shape
    ref_lat = inps.ref_ll[0]
    ref_lon = inps.ref_ll[1]

    # write simple attribution
    atr_geo = dict()
    atr_geo['WIDTH'] = str(colm)
    atr_geo['LENGTH'] = str(row)
    atr_geo['FILE_TYPE'] = 'geometry'

    atr_vel = dict()
    atr_vel['WIDTH'] = str(colm)
    atr_vel['LENGTH'] = str(row)
    atr_vel['REF_LAT'] = str(ref_lat)
    atr_vel['REF_LON'] = str(ref_lon)
    atr_vel['FILE_TYPE'] = 'velocity'

    lat_data = dict()
    lat_data['latitude'] = lat_joined

    lon_data = dict()
    lon_data['longitude'] = lon_joined

    vel_data = dict()
    vel_data['velocity'] = vel_joined

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
    vel_filename = output_dir + data_outname + '.h5'

    # write h5 file
    writefile.write(datasetDict=lat_data, out_file=lat_filename, metadata=atr_geo)
    writefile.write(datasetDict=lon_data, out_file=lon_filename, metadata=atr_geo)
    writefile.write(datasetDict=geo_data, out_file=geo_filename, metadata=atr_geo)
    writefile.write(datasetDict=vel_data, out_file=vel_filename, metadata=atr_vel)

    print('Finish!')

def main(iargs=None):
    inps = cmd_line_parse(iargs)   
    
    lat_joined, lon_joined, inc_joined, azi_joined, hgt_joined, vel_joined = concatenate(inps)
    write(lat_joined, lon_joined, inc_joined, azi_joined, hgt_joined, vel_joined, inps)
######################################################################################
if __name__ == '__main__':
    main()
