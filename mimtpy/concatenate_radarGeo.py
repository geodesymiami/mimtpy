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
import copy
import pandas as pd
import math
from sklearn.neighbors import KNeighborsClassifier

import mintpy
import mintpy.workflow
from mintpy.utils import readfile, ptime, writefile, utils as ut
from mintpy.objects import timeseries

import mimtpy.workflow
######################################################################################
NOTE = """
Please Note:
1. Four types of data are supported: velocity, timeseries, maskPS and maskTempCoh
2. This script concatenates two datasets together. The input files are object datasets, their corresponding geometryRadar.h5 files
3. If a batch concatenation needed, please use the concatnate_patches.py script.
4. For timeseries, datasets should have same reference date

"""

EXAMPLE = """example:

    concatenate_radarGeo.py miaplpy_NE/velocity_msk.h5  miaplpy_NNE/velocity_msk.h5 --geo_file1 miaplpy_NE/inputs/geometryRadar.h5 --geo_file2 miaplpy_NNE/intpus/geometryRadar.h5 --datatype vel --geo_write --output NE_NNE --outdyyir miaplpyBig

    concatenate_radarGeo.py miaplpy_NE/velocity_msk.h5  miaplpy_NNE/velocity_msk.h5 --geo_file1 miaplpy_NE/inputs/geometryRadar.h5 --geo_file2 miaplpy_NNE/intpus/geometryRadar.h5 --datatype vel --output NE_NNE --outdyyir miaplpyBig
    
    concatenate_radarGeo.py miaplpy_NE/timeseries_msk.h5 miaplpy_NNE/timeseries_msk.h5 --geo_file1 miaplpy_NE/inputs/geometryRadar.h5 --geo_file2 miaplpy_NNE/inputs/geometryRadar.h5 --datatype ts --output NE_NNE --outdir ./miaplpyBig/

    concatenate_radarGeo.py miaplpy_NE/maskPS.h5  miaplpy_NNE/maskPS.h5 --geo_file1 miaplpy_NE/inputs/geometryRadar.h5 --geo_file2 miaplpy_NNE/intpus/geometryRadar.h5 --datatype maskPS --output NE_NNE --outdir miaplpyBig

    concatenate_radarGeo.py miaplpy_NE/maskTempCoh.h5  miaplpy_NNE/maskTempCoh.h5 --geo_file1 miaplpy_NE/inputs/geometryRadar.h5 --geo_file2 miaplpy_NNE/intpus/geometryRadar.h5 --datatype maskTC --output NE_NNE --outdyyir miaplpyBig
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Concatenate miaplpy patches',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=NOTE + '\n' + EXAMPLE)

    parser.add_argument('patch_files', nargs='+', type=str, help='two displacement datasets to be concatenated \n')

    parser.add_argument('--geo_file1', nargs=1, type=str, help='geometryRadar file of dataset1. \n')
    
    parser.add_argument('--geo_file2', nargs=1, type=str, help='geometryRadar file of dataset2. \n')
    
    parser.add_argument('--datatype', nargs=1, type=str, help='data type: vel, ts, maskPS, maskTC\n')
   
    parser.add_argument('--geo_write',action='store_true', default=False, help='whether write the concatenated geometryRadar.h5 results. \n')
 
    parser.add_argument('--output', nargs=1, type=str, help='output name of concatenated data. \n')

    parser.add_argument('--outdir',dest='outdir',nargs=1,
                        help='output directory to store the concatenated results')
    
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps

def haversin(theta):
    v = np.sin(theta / 2)
    return v * v

def distance2points(lat1, lon1, lat2, lon2):
    radius = 6370

    lat1 = math.radians(lat1)
    lon1 = math.radians(lon1)
    lat2 = math.radians(lat2)
    lon2 = math.radians(lon2)
    
    dlon = lon2 - lon1
    dlat = lat2 - lat1

    h = haversin(dlat) + np.cos(lat1) * np.cos(lat2) * haversin(dlon)

    dis = 2 * radius * np.sin(np.sqrt(h))

    return dis

def overlay_lalo(lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten):
    # calculate the overlay region between two images
    lat1_min = np.min(lat1_flatten)
    lat1_max = np.max(lat1_flatten)
    lon1_min = np.min(lon1_flatten)
    lon1_max = np.max(lon1_flatten)

    lat2_min = np.min(lat2_flatten)
    lat2_max = np.max(lat2_flatten)
    lon2_min = np.min(lon2_flatten)
    lon2_max = np.max(lon2_flatten)

    # calculate the overlay lat and lon
    over_lat_min = max(lat1_min,lat2_min)
    over_lon_min = max(lon1_min,lon2_min)
    over_lat_max = min(lat1_max,lat2_max)
    over_lon_max = min(lon1_max,lon2_max)

    return over_lat_min, over_lon_min, over_lat_max, over_lon_max

def PS_overlay(latlon, over_lat_min, over_lon_min, over_lat_max, over_lon_max):
    # extract the PS points located in the overlay region
    flag_lat = np.where((latlon[:,0]<over_lat_max) & (latlon[:,0]>over_lat_min)) 
    flag_lon = np.where((latlon[:,1]<over_lon_max) & (latlon[:,1]>over_lon_min))

    flag = np.intersect1d(flag_lat, flag_lon)
    PS_num = len(flag)
    print('The total number of PS located at the overlay region is %d' % PS_num)

    return flag

def calculate_offset_matrix(vel_ref, lat_ref, lon_ref, vel_aff, lat_aff, lon_aff):
    # calculate the offset between reference and afflicate image overlay region
    # constructure PD frame
    find_PS = {'lon':lon_ref, 'lat':lat_ref,'vel':vel_ref}
    find = pd.DataFrame(find_PS)
    data_PS = {'lon':lon_aff, 'lat':lat_aff,'vel':vel_aff}
    data = pd.DataFrame(data_PS)
    
    data_fit = data.iloc[:, [0, 1]]
    y = [1] * len(data_fit)

    find_x = find.iloc[:, [0, 1]]

    knn = KNeighborsClassifier(n_neighbors=1,
                               algorithm='ball_tree',
                               metric=lambda s1, s2: distance2points(*s1, *s2))

    # train the knn model
    knn.fit(data_fit, y)
    # calculate the nearest point
    distance, point = knn.kneighbors(find_x, n_neighbors=1, return_distance=True)

    # calculate the median of difference between reference image and affilicate image
    #offset_overlay = m_overlay - s_overlay
    #offset = np.nanmedian(offset_overlay) 
    offset = np.array([[1]])
    for i, row in find.iterrows():
        tmp = data.iloc[point[i]]
        if distance[i][0] < 0.005: 
            find_s = pd.DataFrame(row).T
            vel_ref_value = find_s.loc[i, 'vel']
            vel_aff_value = tmp['vel'].get(point[i][0])
            vel_delta = vel_ref_value - vel_aff_value
            offset = np.append(offset, vel_delta)
    
    return offset[1:]

def concatenate_process(data1_flatten, data2_flatten, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten):
    # for two datasets, do concatenation
    latlon1 = np.hstack((np.transpose(np.array([lat1_flatten])), np.transpose(np.array([lon1_flatten]))))
    latlon2 = np.hstack((np.transpose(np.array([lat2_flatten])), np.transpose(np.array([lon2_flatten]))))
    
    # compare the value of reference between orginal data and concatenated data
    #import pdb
    #pdb.set_trace()
    #print(vel1[int(vel1_atr['REF_Y']), int(vel1_atr['REF_X'])])
    #vel_joined = np.transpose(np.array([np.hstack((vel1_flatten, vel2_flatten))]))
    #pos = (int(vel1_atr['REF_Y']) + 1) + ((int(vel1_atr['REF_X']) + 1) - 1) * int(vel1_atr['LENGTH']) - 1
    #print(vel_joined[pos])

    # calculate the overlay latlon
    over_lat_min, over_lon_min, over_lat_max, over_lon_max = overlay_lalo(lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten)

    # calculate the PS points located in the overlay region
    PS_flag1 = PS_overlay(latlon1, over_lat_min, over_lon_min, over_lat_max, over_lon_max)
    PS_flag2 = PS_overlay(latlon2, over_lat_min, over_lon_min, over_lat_max, over_lon_max)

    # calculate the offset between dataset1(reference image) and dataset2 (afflicate image)
    # extract the PS points whose vel value is not Nan
    data1_tmp = data1_flatten[PS_flag1]
    mask1 = ~np.isnan(data1_tmp)
    data2_tmp = data2_flatten[PS_flag2]
    mask2 = ~np.isnan(data2_tmp)   
 
    data1_overlay_num = len(data1_tmp[mask1])
    print('The Nan PS point of reference image located in the overlay region is %d' % data1_overlay_num)
    data2_overlay_num = len(data2_tmp[mask2])
    print('The Nan PS point of affilicate image located in the overlay region is %d' % data2_overlay_num)

    if data1_overlay_num <= data2_overlay_num:
        data_ref = data1_tmp[mask1]
        lat_ref = latlon1[:, 0][PS_flag1][mask1]
        lon_ref = latlon1[:, 1][PS_flag1][mask1]
        data_aff = data2_tmp[mask2]
        lat_aff = latlon2[:, 0][PS_flag2][mask2]
        lon_aff = latlon2[:, 1][PS_flag2][mask2]
        offset = calculate_offset_matrix(data_ref, lat_ref, lon_ref, data_aff, lat_aff, lon_aff)
    else:
        data_ref = data2_tmp[mask2]
        lat_ref = latlon2[:, 0][PS_flag2][mask2]
        lon_ref = latlon2[:, 1][PS_flag2][mask2]
        data_aff = data1_tmp[mask1]
        lat_aff = latlon1[:, 0][PS_flag1][mask1]
        lon_aff = latlon1[:, 1][PS_flag1][mask1]
        offset = calculate_offset_matrix(data_ref, lat_ref, lon_ref, data_aff, lat_aff, lon_aff)
        offset *= (-1)

    overlay_offset = np.nanmedian(offset)

    # adjust the affiliate image
    data2_flatten_adjust = data2_flatten + overlay_offset

    # do concatenation
    data_joined = np.transpose(np.array([np.hstack((data1_flatten, data2_flatten))]))

    return data_joined

def concatenate_vel(inps, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten):
    data1 = inps.patch_files[0] # reference image
    data2 = inps.patch_files[1] # affiliate image
   
    print('Read the first dataset') 
    vel1, vel1_atr = readfile.read(data1, datasetName='velocity')
    vel1_flatten = vel1.flatten('F')
    
    print('Read the second dataset') 
    vel2, vel2_atr = readfile.read(data2, datasetName='velocity')
    vel2_flatten = vel2.flatten('F')

    vel_joined = concatenate_process(vel1_flatten, vel2_flatten, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten)

    return vel_joined, vel1_atr

def concatenate_ts(inps, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten):
    data1 = inps.patch_files[0]
    data2 = inps.patch_files[1]

    print('Read the first dataset') 
    ts1, ts1_atr = readfile.read(data1, datasetName='timeseries')
    print('Read the second dataset') 
    ts2, ts2_atr = readfile.read(data2, datasetName='timeseries')

    bperp_date1 = h5py.File(data1,'r')
    bperp1 = bperp_date1['/bperp']
    dateList1 = timeseries(data1).get_date_list()
    #date1 = np.array(bperp_date1['/date']).tolist()

    bperp_date2 = h5py.File(data2,'r')
    bperp2 = bperp_date2['/bperp']
    dateList2 = timeseries(data2).get_date_list()
    #date2 = np.array(bperp_date2['/date']).tolist()

    # judging whether dominant and affiliate data have same dimension
    dim1 = ts1.shape[0]
    rows1, colms1 = ts1.shape[1:3]
    dim2 = ts2.shape[0]
    rows2, colms2 = ts2.shape[1:3]

    #calculate the intersected date betwee two datasets    
    date_final, Date1, Date2, bperp = mimtpy.concatenate_offset.date_match(dateList1, dateList2, dim1, dim2, bperp1, bperp2)
    
    # prepare to concatenate
    join_dim = len(Date1)
    pnum = len(lat1_flatten) + len(lat2_flatten)

    ts_join_datast = dict()
    ts_join = np.empty(shape=(join_dim, pnum, 1), dtype=float)
    
    # do concatenation
    i = 0
    for date1, date2 in zip(Date1, Date2):
        print('Process displacement data of date %s' % date1)
        dis1 = readfile.read(data1, datasetName=date1)[0]
        dis1_flatten = dis1.flatten('F')
        dis2 = readfile.read(data2, datasetName=date2)[0]
        dis2_flatten = dis2.flatten('F')

        ts_join[i, :, :] = concatenate_process(dis1_flatten, dis2_flatten, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten)
        i += 1

    ts_join_dataset['bperp'] = np.array(bperp, dtype=float)
    ts_join_dataset['date'] = np.array(date_final, dtype=np.string_)
    ts_join_dataset['timeseries'] = ts_join


    return ts_join_dataset, ts1_atr, date_final
        
def concatenate_mask(inps):
    """concantenate maskPS data"""
    data1 = inps.patch_files[0]
    data2 = inps.patch_files[1]

    print('Read the first dataset') 
    msk1, msk1_atr = readfile.read(data1) 
    print('Read the second dataset') 
    msk2, msk2_atr = readfile.read(data2)

    msk1_flatten = msk1.flatten('F') # flatten matrix according colms
    msk2_flatten = msk2.flatten('F') # flatten matrix according colms

    # do concatenation
    msk_joined = np.transpose(np.array([np.hstack((msk1_flatten, msk2_flatten))]))

    return msk_joined, msk1_atr

def concatenate_geo(inps):
    """concatenate geometry data"""
    data_geo1 = inps.geo_file1[0]
    data_geo2 = inps.geo_file2[0]

    print('Read the first dataset') 
    lat1, lat_atr1 = readfile.read(data_geo1, datasetName='latitude')
    lon1, lon_atr1 = readfile.read(data_geo1, datasetName='longitude')
    inc1, inc_atr1 = readfile.read(data_geo1, datasetName='incidenceAngle')
    azi1, azi_atr1 = readfile.read(data_geo1, datasetName='azimuthAngle')
    hgt1, hgt_atr1 = readfile.read(data_geo1, datasetName='height')

    print('Read the second dataset') 
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

def write_vel(vel_joined, vel_atr, inps):
    
    row, colm = vel_joined.shape

    atr_vel = dict()
    atr_vel['WIDTH'] = str(colm)
    atr_vel['LENGTH'] = str(row)
    atr_vel['REF_X'] = str(1 - 1)
    atr_vel['REF_Y'] = str((int(vel_atr['REF_Y']) + 1) + ((int(vel_atr['REF_X']) + 1) - 1) * int(vel_atr['LENGTH']) - 1)
    #print('The value is %d' % vel_joined[int(atr_vel['REF_Y'])])
    atr_vel['FILE_TYPE'] = 'velocity'
    
    vel_data = dict()
    vel_data['velocity'] = vel_joined

    output_dir = inps.outdir[0]
    outname = inps.output[0]
    
    vel_filename = output_dir + inps.datatype[0] + '_' + outname + '.h5'
    writefile.write(datasetDict=vel_data, out_file=vel_filename, metadata=atr_vel)

    return 

def write_ts(ts_joined_dataset, ts_atr, date_final, inps):
    row, colm = ts_joined_dataset['timeseries'].shape[1: ]

    atr_ts = ts_atr
    atr_ts['WIDTH'] = str(colm)
    atr_ts['LENGTH'] = str(row)
    atr_ts['REF_LAT'] = str(1 - 1)
    atr_ts['REF_LON'] = str((int(vel_atr['REF_Y']) + 1) + ((int(vel_atr['REF_X']) + 1) - 1) * int(vel_atr['LENGTH']) - 1)
    atr_ts['FILE_TYPE'] = 'timeseries'
    
    output_dir = inps.outdir[0]
    outname = inps.output[0]

    ts_filename = output_dir + inps.datatype[0] + '_' + outname + '.h5'
    writefile.write(datasetDict=ts_data, out_file=ts_filename, metadata=atr_ts)
    
    return

def write_mask(msk_joined, inps):
    row, colm = msk_joined.shape
    
    # write simple attribution
    atr_msk = dict()
    atr_msk['WIDTH'] = str(colm)
    atr_msk['LENGTH'] = str(row)
    atr_msk['FILE_TYPE'] = 'mask'

    msk_data = dict()
    msk_data['mask'] = msk_joined

    output_dir = inps.outdir[0]
    outname = inps.output[0]
    msk_outname = inps.datatype[0] + '_' + outname
    msk_filename = output_dir + '/' + msk_outname + '.h5'

    writefile.write(datasetDict=msk_data, out_file=msk_filename, metadata=atr_msk)
     
    return 

def write_geo(lat_joined, lon_joined, inc_joined, azi_joined, hgt_joined, inps):

    row, colm = lat_joined.shape

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
    outname = inps.output[0]
    geo_outname = 'geometry_' + outname

    #suffix = "_".join(data_outname.split('_')[2:]) 
    #lat_filename = output_dir + '/latitude_' + suffix + '.h5'
    #lon_filename = output_dir + '/longitude_' + suffix + '.h5'

    geo_filename = output_dir + '/' + geo_outname + '.h5'

    # write h5 file
    #writefile.write(datasetDict=lat_data, out_file=lat_filename, metadata=atr_geo)
    #writefile.write(datasetDict=lon_data, out_file=lon_filename, metadata=atr_geo)
    writefile.write(datasetDict=geo_data, out_file=geo_filename, metadata=atr_geo)

    print('Finish!')

def main(iargs=None):
    inps = cmd_line_parse(iargs)   
    
    print('process the geometry info')
    lat_joined, lon_joined, inc_joined, azi_joined, hgt_joined, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten = concatenate_geo(inps)
    if inps.geo_write:
        write_geo(lat_joined, lon_joined, inc_joined, azi_joined, hgt_joined, inps)

    print('process %s data' % inps.datatype[0])
    if inps.datatype[0] == 'vel':
        vel_joined, vel_atr = concatenate_vel(inps, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten)
        write_vel(vel_joined, vel_atr, inps)
    elif inps.datatype[0] == 'ts':
        ts_join_dataset, ts1_atr, date_final = concatenate_ts(inps, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten)
        write_ts(ts_joined_dataset, ts_atr, date_final, inps)
    elif inps.datatype[0] == 'maskPS':
        msk_joined, msk_atr = concatenate_mask(inps)
        write_mask(msk_joined, inps)
    elif inps.datatype[0] == 'maskTC':
        msk_joined, msk_atr = concatenate_mask(inps)
        write_mask(msk_joined, inps)
        

######################################################################################
if __name__ == '__main__':
    main()
