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
from functools import partial

import mintpy
import mintpy.workflow
from mintpy.utils import readfile, ptime, writefile, utils as ut
from mintpy.objects import timeseries

import mimtpy.workflow
######################################################################################
NOTE = """
Please Note:
1. Four types of data are supported: velocity, timeseries, maskPS and maskTempCoh
2. This script concatenates two datasets together. The input files are object datasets, their corresponding geometryRadar.h5 files, and the whole region geometryRadar.h5 file processed with 1:1 looks
3. If a batch concatenation needed, please use the concatnate_patches.py script.
4. For timeseries, datasets should have same reference date

"""

EXAMPLE = """example:

    concatenate_radarGeo.py miaplpy_NE/velocity_msk.h5  miaplpy_NNE/velocity_msk.h5 --geo_file1 miaplpy_NE/inputs/geometryRadar.h5 --geo_file2 miaplpy_NNE/intpus/geometryRadar.h5 --geo_full ./inputs/geometryRadar.h5 --datatype vel --geo_write --output NE_NNE --outdir miaplpyBig

    concatenate_radarGeo.py miaplpy_NE/velocity_msk.h5  miaplpy_NNE/velocity_msk.h5 --geo_file1 miaplpy_NE/inputs/geometryRadar.h5 --geo_file2 miaplpy_NNE/intpus/geometryRadar.h5 --geo_full ./inputs/geometryRadar.h5 --datatype vel --output NE_NNE --outdir miaplpyBig
    
    concatenate_radarGeo.py miaplpy_NE/timeseries_msk.h5 miaplpy_NNE/timeseries_msk.h5 --geo_file1 miaplpy_NE/inputs/geometryRadar.h5 --geo_file2 miaplpy_NNE/inputs/geometryRadar.h5 --geo_full ./inputs/geometryRadar.h5 --datatype ts --output NE_NNE --outdir ./miaplpyBig/

    concatenate_radarGeo.py miaplpy_NE/maskPS.h5  miaplpy_NNE/maskPS.h5 --geo_file1 miaplpy_NE/inputs/geometryRadar.h5 --geo_file2 miaplpy_NNE/intpus/geometryRadar.h5 --geo_full ./inputs/geometryRadar.h5 --datatype maskPS --output NE_NNE --outdir miaplpyBig

    concatenate_radarGeo.py miaplpy_NE/maskTempCoh.h5  miaplpy_NNE/maskTempCoh.h5 --geo_file1 miaplpy_NE/inputs/geometryRadar.h5 --geo_file2 miaplpy_NNE/intpus/geometryRadar.h5 --geo_full ./inputs/geometryRadar.h5 --datatype maskTC --output NE_NNE --outdir miaplpyBig
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Concatenate miaplpy patches',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=NOTE + '\n' + EXAMPLE)

    parser.add_argument('patch_files', nargs='+', type=str, help='two displacement datasets to be concatenated \n')

    parser.add_argument('--geo_file1', nargs=1, type=str, help='geometryRadar file of dataset1. \n')
    
    parser.add_argument('--geo_file2', nargs=1, type=str, help='geometryRadar file of dataset2. \n')
    
    parser.add_argument('--geo_full', nargs=1, type=str, help='Whole region geometryRadar.h5 file processed with 1:1 looks. \n')
    
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

def flatten_trans(x):
    original_shape = x.shape
    return partial(np.reshape, newshape=original_shape)

def geo_position(lat_sub, lon_sub, lat, lon, ref_flag=None):
    if ref_flag is None:
        lat_corner = lat_sub[0][0] # use the upper left point
        lon_corner = lon_sub[0][0]
    elif ref_flag == 1:
        lat_corner = lat_sub[0][0] # use the upper left point
        lon_corner = lon_sub[0][0]
    elif ref_flag == 2:
        lat_corner = lat_sub[0][-1] # use the upper right point
        lon_corner = lon_sub[0][-1]

    lat_flag = (lat == lat_corner)
    lon_flag = (lon == lon_corner)

    flag = lat_flag * lon_flag
    pos_flag = np.where(flag == True)

    return pos_flag[0][0], pos_flag[1][0]

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
        if distance[i][0] < 0.006: 
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
    print('The overlay offset is %f' % overlay_offset)

    # adjust the affiliate image
    data2_flatten_adjust = data2_flatten + overlay_offset

    return data2_flatten_adjust

def concatenate_2D(geo_ref, geo_aff, row_ref, col_ref, row_aff, col_aff, row_a_r, col_a_r, ref_flag, data_type, geometry=None):
    if ref_flag == 1:
        # join the data
        row_sum = geo_aff.shape[0] + row_a_r
        col_sum = geo_aff.shape[1] + col_a_r
        # join geo
        print('The concatenated data size is {} rows and {} colms'.format(row_sum, col_sum))
        geo_joined = np.ones((row_sum, col_sum)) * np.nan
        geo_joined[0: geo_ref.shape[0], 0: geo_ref.shape[1]] = geo_ref
        geo_joined[row_a_r: , col_a_r: ] = geo_aff
        if data_type == 'geo':
            geo_joined[0: row_a_r, geo_ref.shape[1]: ] = geometry[row_ref: row_aff, col_ref + geo_ref.shape[1]: col_aff + geo_aff.shape[1]]
            geo_joined[geo_ref.shape[0]: , 0: col_a_r] = geometry[row_ref + geo_ref.shape[0] : row_aff + geo_aff.shape[0], col_ref: col_aff]
        elif data_type == 'msk':
            geo_joined[0: row_a_r, geo_ref.shape[1]: ] = 0
            geo_joined[geo_ref.shape[0]: , 0: col_a_r] = 0

    elif ref_flag == 2:
        # join the data
        row_sum = geo_aff.shape[0] + row_a_r
        col_sum = geo_aff.shape[1] + geo_ref.shape[1] - col_a_r
        # join geo
        print('The concatenated data size is {} rows and {} colms'.format(row_sum, col_sum))
        geo_joined = np.ones((row_sum, col_sum)) * np.nan
        geo_joined[0: geo_ref.shape[0], geo_aff.shape[1] - col_a_r: ] = geo_ref
        geo_joined[row_a_r: , 0: geo_aff.shape[1]] = geo_aff
        if data_type == 'geo':
            geo_joined[0: row_a_r, 0: geo_aff.shape[1] - col_a_r] = geometry[row_ref: row_aff, col_aff: col_ref]
            geo_joined[geo_ref.shape[0]:, geo_aff.shape[1]: ] = geometry[row_ref + geo_ref.shape[0]: row_aff + geo_aff.shape[0], col_aff + geo_aff.shape[1]: col_ref + geo_ref.shape[1]]
        elif datatype == 'msk':
            geo_joined[0: row_a_r, 0: geo_aff.shape[1] - col_a_r] = 0
            geo_joined[geo_ref.shape[0]:, geo_aff.shape[1]: ] = 0

    return geo_joined 

def concatenate_vel(inps, lat_ref_flatten, lon_ref_flatten, lat_aff_flatten, lon_aff_flatten, unflatten_trans_aff, row_ref, col_ref, row_aff, col_aff, row_a_r, col_a_r, ref_flag):
    if ref_flag == 1:
        data_ref = inps.patch_files[0]
        data_aff = inps.patch_files[1]
    else:
        data_ref = inps.patch_files[1]
        data_aff = inps.patch_files[0]   


    print('Read the reference dataset') 
    vel_ref, vel_ref_atr = readfile.read(data_ref, datasetName='velocity')
    vel_ref_flatten = vel_ref.flatten()
    
    print('Read the affilicate dataset') 
    vel_aff, vel_aff_atr = readfile.read(data_aff, datasetName='velocity')
    vel_aff_flatten = vel_aff.flatten()

    vel_aff_flatten_adjust = concatenate_process(vel_ref_flatten, vel_aff_flatten, lat_ref_flatten, lon_ref_flatten, lat_aff_flatten, lon_aff_flatten)
    vel_aff_adjust = unflatten_trans_aff(vel_aff_flatten_adjust)

    # generate 2D concatenation results
    data_type = 'vel'
    vel_joined = concatenate_2D(vel_ref, vel_aff_adjust, row_ref, col_ref, row_aff, col_aff, row_a_r, col_a_r, ref_flag, data_type)

    # adjust the attribute table
    vel_atr = vel_ref_atr
    vel_atr['LENGTH'] = vel_joined.shape[0]
    vel_atr['WIDTH'] = vel_joined.shape[1]
    
    return vel_joined, vel_atr

def concatenate_ts(inps, lat_ref_flatten, lon_ref_flatten, lat_aff_flatten, lon_aff_flatten, unflatten_trans_aff, row_ref, col_ref, row_aff, col_aff, row_a_r, col_a_r, ref_flag):
    if ref_flag == 1:
        data_ref = inps.patch_files[0]
        data_aff = inps.patch_files[1]
    else:
        data_ref = inps.patch_files[1]
        data_aff = inps.patch_files[0]   

    print('Read the reference dataset') 
    ts_ref, ts_ref_atr = readfile.read(data_ref, datasetName='timeseries')
    print('Read the affilite dataset') 
    ts_aff, ts2_affatr = readfile.read(data_aff, datasetName='timeseries')

    bperp_date_ref = h5py.File(data_ref,'r')
    bperp_ref = bperp_date_ref['/bperp']
    dateList_ref = timeseries(data_ref).get_date_list()

    bperp_date_aff = h5py.File(data_aff,'r')
    bperp_aff = bperp_date_aff['/bperp']
    dateList_aff = timeseries(data_aff).get_date_list()

    # judging whether dominant and affiliate data have same dimension
    dim_ref = ts_ref.shape[0]
    rows_ref, colms_ref = ts_ref.shape[1:3]
    dim_aff = ts_aff.shape[0]
    rows_aff, colms_aff = ts_aff.shape[1:3]

    #calculate the intersected date betwee two datasets    
    date_final, Date1, Date2, bperp = mimtpy.concatenate_offset.date_match(dateList_ref, dateList_aff, dim_ref, dim_aff, bperp_ref, bperp_aff)
    
    # prepare to concatenate
    join_dim = len(Date1)

    ts_join_dataset = dict()
    if ref_flag == 1:
        # join the data
        row_sum = ts_aff.shape[1] + row_a_r
        col_sum = ts_aff.shape[2] + col_a_r
    elif ref_flag == 2:
        # join the data
        row_sum = ts_aff.shape[1] + row_a_r
        col_sum = ts_aff.shape[2] + ts_ref.shape[2] - col_a_r
    ts_join = np.empty(shape=(join_dim, row_sum, col_sum), dtype=float)
    # do concatenation
    i = 0
    for date1, date2 in zip(Date1, Date2):
        print('Process displacement data of date %s' % date1)
        dis_ref = readfile.read(data_ref, datasetName=date1)[0]
        dis_ref_flatten = dis_ref.flatten()
        dis_aff = readfile.read(data_aff, datasetName=date2)[0]
        dis_aff_flatten = dis_aff.flatten()

        dis_aff_flatten_adjust = concatenate_process(dis_ref_flatten, dis_aff_flatten, lat_ref_flatten, lon_ref_flatten, lat_aff_flatten, lon_aff_flatten)
        dis_aff_adjust = unflatten_trans_aff(dis_aff_flatten_adjust)
        # generate 2D concatenation results
        data_type = 'ts'
        ts_join[i, :, :] = concatenate_2D(dis_ref, dis_aff_adjust, row_ref, col_ref, row_aff, col_aff, row_a_r, col_a_r, ref_flag, data_type)
        i += 1

    ts_join_dataset['bperp'] = np.array(bperp, dtype=float)
    ts_join_dataset['date'] = np.array(date_final, dtype=np.string_)
    ts_join_dataset['timeseries'] = ts_join

    # adjust the attribute table
    ts_atr = ts_ref_atr
    ts_atr['LENGTH'] = ts_join.shape[1]
    ts_atr['WIDTH'] = ts_join.shape[2]

    return ts_join_dataset, ts_atr, date_final
        
def concatenate_mask(inps, row_ref, col_ref, row_aff, col_aff, row_a_r, col_a_r, ref_flag):
    """concantenate maskPS data"""
    if ref_flag == 1:
        data_ref = inps.patch_files[0]
        data_aff = inps.patch_files[1]
    else:
        data_ref = inps.patch_files[1]
        data_aff = inps.patch_files[0]   

    print('Read the reference dataset') 
    msk_ref, msk_ref_atr = readfile.read(data_ref) 
    print('Read the affilite dataset') 
    msk_aff, msk_aff_atr = readfile.read(data_aff)

    data_type = 'msk'
    msk_joined = concatenate_2D(msk_ref, msk_aff, row_ref, col_ref, row_aff, col_aff, row_a_r, col_a_r, ref_flag, data_type)

    # adjust the attribute table
    msk_atr = msk_ref_atr
    msk_atr['LENGTH'] = msk_joined.shape[0]
    msk_atr['WIDTH'] = msk_joined.shape[1]

    return msk_joined, msk_atr

def concatenate_geo(inps, ref_flag):
    """concatenate geometry data"""
    if ref_flag == 1:
        data_geo_ref = inps.geo_file1[0]
        data_geo_aff = inps.geo_file2[0]
    else:
        data_geo_ref = inps.geo_file2[0]
        data_geo_aff = inps.geo_file1[0]
        
    geo_full = inps.geo_full[0]

    print('Read the geometry data for the full region') 
    lat_full = readfile.read(geo_full, datasetName='latitude')[0]
    lon_full = readfile.read(geo_full, datasetName='longitude')[0]
    inc_full = readfile.read(geo_full, datasetName='incidenceAngle')[0]
    azi_full = readfile.read(geo_full, datasetName='azimuthAngle')[0]
    hgt_full = readfile.read(geo_full, datasetName='height')[0]
    
    print('Read the reference dataset') 
    lat_ref, lat_atr_ref = readfile.read(data_geo_ref, datasetName='latitude')
    lon_ref, lon_atr_ref = readfile.read(data_geo_ref, datasetName='longitude')
    inc_ref, inc_atr_ref = readfile.read(data_geo_ref, datasetName='incidenceAngle')
    azi_ref, azi_atr_ref = readfile.read(data_geo_ref, datasetName='azimuthAngle')
    hgt_ref, hgt_atr_ref = readfile.read(data_geo_ref, datasetName='height')

    print('Read the affilicate dataset') 
    lat_aff, lat_atr_aff = readfile.read(data_geo_aff, datasetName='latitude')
    lon_aff, lon_atr_aff = readfile.read(data_geo_aff, datasetName='longitude')
    inc_aff, inc_atr_aff = readfile.read(data_geo_aff, datasetName='incidenceAngle')
    azi_aff, azi_atr_aff = readfile.read(data_geo_aff, datasetName='azimuthAngle')
    hgt_aff, hgt_atr_aff = readfile.read(data_geo_aff, datasetName='height')
    
    lat_ref_flatten = lat_ref.flatten() # flatten matrix according rows
    lon_ref_flatten = lon_ref.flatten()

    lat_aff_flatten = lat_aff.flatten()
    lon_aff_flatten = lon_aff.flatten()
    
    # calculate the unflatten pattern
    unflatten_trans_ref = flatten_trans(lat_ref)
    unflatten_trans_aff = flatten_trans(lat_aff)

    print('Find the position of reference data using lat/lon')
    row_ref, col_ref = geo_position(lat_ref, lon_ref, lat_full, lon_full)
    print('Find the position of affilite data using lat/lon')
    row_aff, col_aff = geo_position(lat_aff, lon_aff, lat_full, lon_full)
    print('Find the position of affilite data in reference data using lat/lon')
    row_a_r, col_a_r = geo_position(lat_aff, lon_aff, lat_ref, lon_ref, ref_flag)
    
    data_type = 'geo'
    lat_joined = concatenate_2D(lat_ref, lat_aff, row_ref, col_ref, row_aff, col_aff, row_a_r, col_a_r, ref_flag, data_type, lat_full)
    lon_joined = concatenate_2D(lon_ref, lon_aff, row_ref, col_ref, row_aff, col_aff, row_a_r, col_a_r, ref_flag, data_type, lon_full)
    inc_joined = concatenate_2D(inc_ref, inc_aff, row_ref, col_ref, row_aff, col_aff, row_a_r, col_a_r, ref_flag, data_type, inc_full)
    azi_joined = concatenate_2D(azi_ref, azi_aff, row_ref, col_ref, row_aff, col_aff, row_a_r, col_a_r, ref_flag, data_type, azi_full)
    hgt_joined = concatenate_2D(hgt_ref, hgt_aff, row_ref, col_ref, row_aff, col_aff, row_a_r, col_a_r, ref_flag, data_type, hgt_full)

    return lat_joined, lon_joined, inc_joined, azi_joined, hgt_joined, lat_ref_flatten, lon_ref_flatten, lat_aff_flatten, lon_aff_flatten, unflatten_trans_ref, unflatten_trans_aff, row_ref, col_ref, row_aff, col_aff, row_a_r, col_a_r

def write_vel(vel_joined, vel_atr, inps):
    
    row, colm = vel_joined.shape

    atr_vel = dict()
    atr_vel['WIDTH'] = str(colm)
    atr_vel['LENGTH'] = str(row)
    atr_vel['FILE_TYPE'] = 'velocity'
    
    vel_data = dict()
    vel_data['velocity'] = vel_joined

    output_dir = inps.outdir[0]
    outname = inps.output[0]
    
    vel_filename = output_dir + '/' + inps.datatype[0] + '_' + outname + '.h5'
    writefile.write(datasetDict=vel_data, out_file=vel_filename, metadata=atr_vel)

    return 

def write_ts(ts_joined_dataset, ts_atr, date_final, inps):
    row, colm = ts_joined_dataset['timeseries'].shape[1: ]

    atr_ts = ts_atr
    atr_ts['WIDTH'] = str(colm)
    atr_ts['LENGTH'] = str(row)
    atr_ts['FILE_TYPE'] = 'timeseries'
    
    output_dir = inps.outdir[0]
    outname = inps.output[0]

    ts_filename = output_dir + '/' + inps.datatype[0] + '_' + outname + '.h5'
    writefile.write(datasetDict=ts_joined_dataset, out_file=ts_filename, metadata=atr_ts)
    
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

    geo_filename = output_dir + '/' + geo_outname + '.h5'

    # write h5 file
    writefile.write(datasetDict=geo_data, out_file=geo_filename, metadata=atr_geo)

    print('Finish!')

def find_the_reference(inps):
    """Find the reference data based on geo info"""
    geo1_lat = readfile.read(inps.geo_file1[0], datasetName='latitude')[0]
    geo1_lon = readfile.read(inps.geo_file1[0], datasetName='longitude')[0]
    geo2_lat = readfile.read(inps.geo_file2[0], datasetName='latitude')[0]
    geo2_lon = readfile.read(inps.geo_file2[0], datasetName='longitude')[0]

    lat1_ul = geo1_lat[0][0]
    lon1_ul = geo1_lon[0][0]

    lat2_ul = geo2_lat[0][0]
    lon2_ul = geo2_lon[0][0]

    if lat1_ul <= lat2_ul and lon1_ul <= lon2_ul:
        return 1
    else:
        return 2
    #elif lat1_ul < lat2_ul and lon1_ul < lon2_ul:
    #    return 2
    #elif lat1_ul > lat2_ul and lon1_ul > lon2_ul:
    #    return 3
    #elif lat1_ul < lat2_ul and lon1_ul > lon2_ul:
    #    return 4

def main(iargs=None):
    inps = cmd_line_parse(iargs)   
    
    print('Find which dataset is reference dataset')
    ref_flag = find_the_reference(inps)

    print('process the geometry info')
    lat_joined, lon_joined, inc_joined, azi_joined, hgt_joined, lat_ref_flatten, lon_ref_flatten, lat_aff_flatten, lon_aff_flatten, unflatten_trans_ref, unflatten_trans_aff, row_ref, col_ref, row_aff, col_aff, row_a_r, col_a_r = concatenate_geo(inps, ref_flag)
    if inps.geo_write:
        write_geo(lat_joined, lon_joined, inc_joined, azi_joined, hgt_joined, inps)

    print('process %s data' % inps.datatype[0])
    if inps.datatype[0] == 'vel':
        vel_joined, vel_atr = concatenate_vel(inps, lat_ref_flatten, lon_ref_flatten, lat_aff_flatten, lon_aff_flatten, unflatten_trans_aff, row_ref, col_ref, row_aff, col_aff, row_a_r, col_a_r, ref_flag)
        write_vel(vel_joined, vel_atr, inps)
    elif inps.datatype[0] == 'ts':
        ts_join_dataset, ts_atr, date_final = concatenate_ts(inps, lat_ref_flatten, lon_ref_flatten, lat_aff_flatten, lon_aff_flatten, unflatten_trans_aff, row_ref, col_ref, row_aff, col_aff, row_a_r, col_a_r, ref_flag)
        write_ts(ts_join_dataset, ts_atr, date_final, inps)
    elif inps.datatype[0] == 'maskPS':
        msk_joined, msk_atr = concatenate_mask(inps, row_ref, col_ref, row_aff, col_aff, row_a_r, col_a_r, ref_flag)
        write_mask(msk_joined, inps)
    elif inps.datatype[0] == 'maskTC':
        msk_joined, msk_atr = concatenate_mask(inps, row_ref, col_ref, row_aff, col_aff, row_a_r, col_a_r, ref_flag)
        write_mask(msk_joined, inps)
        

######################################################################################
if __name__ == '__main__':
    main()
