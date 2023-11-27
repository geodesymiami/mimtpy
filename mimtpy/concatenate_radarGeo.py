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

    concatenate_radarGeo.py miaplpy_NE/velocity_msk.h5  miaplpy_NNE/velocity_msk.h5 --geo-file1 miaplpy_NE/inputs/geometryRadar.h5 --geo-file2 miaplpy_NNE/intpus/geometryRadar.h5 --geo-full ./inputs/geometryRadar.h5 --geo-write --outdir miaplpyBig

    concatenate_radarGeo.py miaplpy_NE/velocity_msk.h5  miaplpy_NNE/velocity_msk.h5 --geo-full ./inputs/geometryRadar.h5 --outdir miaplpyBig
    
    concatenate_radarGeo.py miaplpy_NE/timeseries_msk.h5 miaplpy_NNE/timeseries_msk.h5 --geo-file1 miaplpy_NE/inputs/geometryRadar.h5 --geo-file2 miaplpy_NNE/inputs/geometryRadar.h5 --geo-full ./inputs/geometryRadar.h5  --out-suffix NE_NNE --outdir ./miaplpyBig/

    concatenate_radarGeo.py miaplpy_NE/maskPS.h5  miaplpy_NNE/maskPS.h5 --geo-file1 miaplpy_NE/inputs/geometryRadar.h5 --geo-file2 miaplpy_NNE/intpus/geometryRadar.h5 --geo-full ./inputs/geometryRadar.h5  --out-suffix NE_NNE --outdir miaplpyBig

    concatenate_radarGeo.py miaplpy_NE/maskTempCoh.h5  miaplpy_NNE/maskTempCoh.h5 --geo-file1 miaplpy_NE/inputs/geometryRadar.h5 --geo-file2 miaplpy_NNE/intpus/geometryRadar.h5 --geo-full ./inputs/geometryRadar.h5  --out-suffix NE_NNE --outdir miaplpyBig
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Concatenate miaplpy patches',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=NOTE + '\n' + EXAMPLE)

    parser.add_argument('patch_files', nargs='+', type=str, help='two displacement datasets to be concatenated \n')

    parser.add_argument('--geo-file1', nargs=1, type=str, default=['default'], 
                        help='geometryRadar file of dataset1.' + 
                        'The default means using geometryRadar.h5 existing in the corresponding patch_files inputs folder\n')
    
    parser.add_argument('--geo-file2', nargs=1, type=str, default=['default'], 
                        help='geometryRadar file of dataset2.' + 
                        'The default means using geometryRadar.h5 existing in the corresponding patch_files inputs folder\n')
    
    parser.add_argument('--geo-full', nargs=1, type=str, help='Whole region geometryRadar.h5 file processed with 1:1 looks. \n')
    
    parser.add_argument('--geo-write',action='store_true', default=False, help='whether write the concatenated geometryRadar.h5 results. \n')
 
    parser.add_argument('--out-suffix', nargs=1,  default=[''], help='suffix of output name of concatenated data. \n')

    parser.add_argument('--outdir',dest='outdir',nargs=1, default=[],
                        help='output directory to store the concatenated results')
    
    parser.add_argument('--list',action='store_true', default=False, help='list the files used for concatenation and check the order of latitude. \n')
    
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps

def flatten_trans(x):
    original_shape = x.shape
    return partial(np.reshape, newshape=original_shape)

def geo_position_backup(lat_sub, lon_sub, lat, lon, ref_flag=None):
    if ref_flag is None:
        lat_corner = lat_sub[0][0] # use the upper left point
        lon_corner = lon_sub[0][0]
    elif ref_flag == 1 or ref_flag == 4 or ref_flag == 6 or ref_flag == 7:
        lat_corner = lat_sub[0][0] # use the upper left point
        lon_corner = lon_sub[0][0]
    elif ref_flag == 2 or ref_flag == 3 or ref_flag == 8:
        lat_corner = lat_sub[0][-1] # use the upper right point
        lon_corner = lon_sub[0][-1]
    elif ref_flag == 5:
        lat_corner = lat_sub[-1][0] # use the lower left point
        lon_corner = lon_sub[-1][0]

    lat_flag = (lat == lat_corner)
    lon_flag = (lon == lon_corner)

    flag = lat_flag * lon_flag
    pos_flag = np.where(flag == True)

    return pos_flag[0][0], pos_flag[1][0]

def get_position(lat_sub, lon_sub, lat, lon):
    row_list = [0, -1]
    col_list = [0, -1]

    row_col_matrix = np.ones((4, 2), dtype=int) * np.nan
    num = 0

    for row in row_list:
        for col in col_list:
            lat_corner = lat_sub[row][col]
            lon_corner = lon_sub[row][col]
            lat_flag = (lat == lat_corner)
            lon_flag = (lon == lon_corner)
            flag = lat_flag * lon_flag
            pos = np.where(flag == True)
            row_col_matrix[num][0] = pos[0][0]
            row_col_matrix[num][1] = pos[1][0]
            num += 1

    return row_col_matrix

def design_joined_matrix(rc_ref, rc_aff):
    row_join_upper = int(np.min([rc_ref[0, 0], rc_aff[0, 0]]))
    row_join_lower = int(np.max([rc_ref[-1, 0], rc_aff[-1, 0]]))
    col_join_left = int(np.min([rc_ref[0, 1], rc_aff[0, 1]]))
    col_join_right = int(np.max([rc_ref[-1, 1], rc_aff[-1, 1]]))
    return row_join_upper, row_join_lower, col_join_left, col_join_right

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

    # calculate the median of difference between reference image and affiliate image
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
    print('The Nan PS point of affiliate image located in the overlay region is %d' % data2_overlay_num)

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

def concatenate_2D(val_ref, val_aff, rc_ref, rc_aff, ref_flag, data_type):
    def common_cal(overlap_ref, overlap_aff, offset):
        temp = (overlap_ref + overlap_aff + offset) / 2
        temp[np.isnan(overlap_ref)] = overlap_aff[np.isnan(overlap_ref)] + offset
        temp[np.isnan(overlap_aff)] = overlap_ref[np.isnan(overlap_aff)]
        return temp

    row_join_start, row_join_end, col_join_start, col_join_end = design_joined_matrix(rc_ref, rc_aff)
    val_join = np.ones((row_join_end - row_join_start + 1, col_join_end - col_join_start + 1)) * np.nan
    row_a_r = int(np.abs(rc_ref[0, 0] - rc_aff[0, 0]))
    col_a_r = int(np.abs(rc_ref[0, 1] - rc_aff[0, 1]))
    print('The ref_flag is %d' % ref_flag)
    if ref_flag == 1 or ref_flag == 4:
        # join geo
        print('Full the concatenated data: {}, {}'.format(val_join.shape[0], val_join.shape[1]))
        overlap_ref = val_ref[row_a_r:, col_a_r:]
        overlap_aff = val_aff[0: val_ref.shape[0] - row_a_r, 0: val_ref.shape[1] - col_a_r]
        offset = np.nanmedian(overlap_ref - overlap_aff)
        print('The offset is %f' % offset)
        val_join[0: val_ref.shape[0], 0: val_ref.shape[1]] = val_ref
        val_join[row_a_r: , col_a_r: ] = val_aff + offset
        val_join[row_a_r: val_ref.shape[0], col_a_r: val_ref.shape[1]] = common_cal(overlap_ref, overlap_aff, offset)
        if data_type == 'msk':
            val_join[np.where(val_join == np.nan)] = 0

    elif ref_flag == 2 or ref_flag == 3:
        # join geo
        print('Full the concatenated data: {}, {}'.format(val_join.shape[0], val_join.shape[1]))
        overlap_ref = val_ref[row_a_r:, 0: val_aff.shape[1] - col_a_r]
        overlap_aff = val_aff[0: val_ref.shape[0] - row_a_r, col_a_r: ]
        offset = np.nanmedian(overlap_ref - overlap_aff)
        print('The offset is %f' % offset)
        val_join[0: val_ref.shape[0], col_a_r: ] = val_ref
        val_join[row_a_r: , 0: val_aff.shape[1]] = val_aff + offset
        val_join[row_a_r: val_ref.shape[0], col_a_r: val_aff.shape[1]] = common_cal(overlap_ref, overlap_aff, offset)
        if data_type == 'msk':
            val_join[np.where(val_join == np.nan)] = 0

    elif ref_flag == 5:
        print('Full the concatenated data: {}, {}'.format(val_join.shape[0], val_join.shape[1]))
        overlap_ref = val_ref[0: val_aff.shape[0] - row_a_r, col_a_r: col_a_r + val_aff.shape[1]]
        overlap_aff = val_aff[row_a_r:, :]
        offset = np.nanmedian(overlap_ref - overlap_aff)
        print('The offset is %f' % offset)
        val_join[row_a_r: , 0: val_ref.shape[1]] = val_ref
        val_join[0: val_aff.shape[0], col_a_r: col_a_r + val_aff.shape[1]] = val_aff + offset
        val_join[row_a_r: val_aff.shape[0], col_a_r: col_a_r + val_aff.shape[1]] = common_cal(overlap_ref, overlap_aff, offset)
        if data_type == 'msk':
            val_join[np.where(val_join == np.nan)] = 0

    elif ref_flag == 6:
        print('Full the concatenated data: {}, {}'.format(val_join.shape[0], val_join.shape[1]))
        overlap_ref = val_ref[row_a_r: , col_a_r: col_a_r + val_aff.shape[1]]
        overlap_aff = val_aff[0: val_ref.shape[0] - row_a_r, :]
        offset = np.nanmedian(overlap_ref - overlap_aff)
        print('The offset is %f' % offset)
        val_join[0: val_ref.shape[0], 0: val_ref.shape[1]] = val_ref
        val_join[row_a_r: , col_a_r: col_a_r + val_aff.shape[1]] = val_aff + offset
        val_join[row_a_r: val_ref.shape[0], col_a_r: col_a_r + val_aff.shape[1]] = common_cal(overlap_ref, overlap_aff, offset)
        if data_type == 'msk':
            val_join[np.where(val_join == np.nan)] = 0

    elif ref_flag == 8:
        print('Full the concatenated data: {}, {}'.format(val_join.shape[0], val_join.shape[1]))
        overlap_ref = val_ref[row_a_r: row_a_r + val_aff.shape[0], 0: val_aff.shape[1] - col_a_r]
        overlap_aff = val_aff[:, col_a_r:]
        offset = np.nanmedian(overlap_ref - overlap_aff)
        print('The offset is %f' % offset)
        val_join[0: val_ref.shape[0], col_a_r: ] = val_ref
        val_join[row_a_r: row_a_r + val_aff.shape[0], 0: val_aff.shape[1]] = val_aff + offset
        val_join[row_a_r: row_a_r + val_aff.shape[0], col_a_r: val_aff.shape[1]] = common_cal(overlap_ref, overlap_aff, offset)
        if data_type == 'msk':
            val_join[np.where(val_join == np.nan)] = 0

    elif ref_flag == 7:
        print('Full the concatenated data: {}, {}'.format(val_join.shape[0], val_join.shape[1]))
        overlap_ref = val_ref[row_a_r: row_a_r + val_aff.shape[0], col_a_r:]
        overlap_aff = val_aff[:, 0: val_ref.shape[1] - col_a_r]
        offset = np.nanmedian(overlap_ref - overlap_aff)
        print('The offset is %f' % offset)
        val_join[0: val_ref.shape[0], 0: val_ref.shape[1]] = val_ref
        val_join[row_a_r: row_a_r + val_aff.shape[0], col_a_r: col_a_r + val_aff.shape[1]] = val_aff + offset
        val_join[row_a_r: row_a_r + val_aff.shape[0], col_a_r: val_ref.shape[1]] = common_cal(overlap_ref, overlap_aff, offset)
        if data_type == 'msk':
            val_join[np.where(val_join == np.nan)] = 0

    return val_join

def concatenate_vel(inps, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten, unflatten_trans1, unflatten_trans2, rc1, rc2, ref_flag):
    ref_No, aff_No, rc_ref, rc_aff,lat_ref_flatten, lon_ref_flatten, lat_aff_flatten, lon_aff_flatten = appoint_ref(rc1, rc2, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten, ref_flag)
    data_ref = inps.patch_files[ref_No]
    data_aff = inps.patch_files[aff_No]

    print('Read the reference dataset') 
    vel_ref, vel_ref_atr = readfile.read(data_ref, datasetName='velocity')
    vel_ref_flatten = vel_ref.flatten()
    
    print('Read the affiliate dataset') 
    vel_aff, vel_aff_atr = readfile.read(data_aff, datasetName='velocity')
    vel_aff_flatten = vel_aff.flatten()

    # using the KNN method to calculate the offset
    #vel_aff_flatten_adjust = concatenate_process(vel_ref_flatten, vel_aff_flatten, lat_ref_flatten, lon_ref_flatten, lat_aff_flatten, lon_aff_flatten)
    #if aff_No == 0:
    #    vel_aff_adjust = unflatten_trans1(vel_aff_flatten_adjust)
    #else:
    #    vel_aff_adjust = unflatten_trans2(vel_aff_flatten_adjust)

    # generate 2D concatenation results
    data_type = 'vel'
    #vel_joined = concatenate_2D(vel_ref, vel_aff_adjust, rc_ref, rc_aff, ref_flag, data_type) # line for KNN method
    vel_joined = concatenate_2D(vel_ref, vel_aff, rc_ref, rc_aff, ref_flag, data_type)

    # adjust the attribute table
    vel_atr = vel_ref_atr
    vel_atr['LENGTH'] = vel_joined.shape[0]
    vel_atr['WIDTH'] = vel_joined.shape[1]
    
    return vel_joined, vel_atr

def concatenate_ts(inps, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten, unflatten_trans1, unflatten_trans2, rc1, rc2, ref_flag):
    ref_No, aff_No, rc_ref, rc_aff,lat_ref_flatten, lon_ref_flatten, lat_aff_flatten, lon_aff_flatten = appoint_ref(rc1, rc2, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten, ref_flag)
    data_ref = inps.patch_files[ref_No]
    data_aff = inps.patch_files[aff_No]

    print('Read the reference dataset') 
    ts_ref, ts_ref_atr = readfile.read(data_ref, datasetName='timeseries')
    print('Read the affiliate dataset') 
    ts_aff, ts_aff_atr = readfile.read(data_aff, datasetName='timeseries')

    # check whether two timeseries have same ref_date
    if ts_ref_atr['REF_DATE'] != ts_aff_atr['REF_DATE']:
        raise ValueError('Two timeseries datasets should have same REF_DATE. PLease change!')

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
    row_join_start, row_join_end, col_join_start, col_join_end = design_joined_matrix(rc_ref, rc_aff)
    row_sum = row_join_end - row_join_start + 1
    col_sum = col_join_end - col_join_start + 1

    ts_join_dataset = dict()
    ts_join = np.empty(shape=(join_dim, row_sum, col_sum), dtype=float)
    # do concatenation
    i = 0
    for date1, date2 in zip(Date1, Date2):
        print('Process displacement data of date %s' % date1)
        dis_ref = readfile.read(data_ref, datasetName=date1)[0]
        dis_ref_flatten = dis_ref.flatten()
        dis_aff = readfile.read(data_aff, datasetName=date2)[0]
        dis_aff_flatten = dis_aff.flatten()

        # using the KNN method to calculate the offset
        #dis_aff_flatten_adjust = concatenate_process(dis_ref_flatten, dis_aff_flatten, lat_ref_flatten, lon_ref_flatten, lat_aff_flatten, lon_aff_flatten)
        #if aff_No == 0:
        #    dis_aff_adjust = unflatten_trans1(dis_aff_flatten_adjust)
        #else:
        #    dis_aff_adjust = unflatten_trans2(dis_aff_flatten_adjust)

        # generate 2D concatenation results
        data_type = 'ts'
        #ts_join[i, :, :] = concatenate_2D(dis_ref, dis_aff_adjust, rc_ref, rc_aff, ref_flag, data_type)
        ts_join[i, :, :] = concatenate_2D(dis_ref, dis_aff, rc_ref, rc_aff, ref_flag, data_type)
        i += 1

    ts_join_dataset['bperp'] = np.array(bperp, dtype=float)
    ts_join_dataset['date'] = np.array(date_final, dtype=np.string_)
    ts_join_dataset['timeseries'] = ts_join

    # adjust the attribute table
    ts_atr = ts_ref_atr
    ts_atr['LENGTH'] = ts_join.shape[1]
    ts_atr['WIDTH'] = ts_join.shape[2]

    return ts_join_dataset, ts_atr, date_final
        
def concatenate_mask(inps, rc1, rc2, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten, ref_flag):
    """concantenate maskPS data"""
    ref_No, aff_No, rc_ref, rc_aff,lat_ref_flatten, lon_ref_flatten, lat_aff_flatten, lon_aff_flatten = appoint_ref(rc1, rc2, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten, ref_flag)
    data_ref = inps.patch_files[ref_No]
    data_aff = inps.patch_files[aff_No]

    print('Read the reference dataset') 
    msk_ref, msk_ref_atr = readfile.read(data_ref) 
    print('Read the affiliate dataset') 
    msk_aff, msk_aff_atr = readfile.read(data_aff)

    #get the type of mask
    if msk_ref_atr['DATA_TYPE'] == 'bool': # the mask is maskTempCoh.h5
        msk_ref = msk_ref + 0
        msk_aff = msk_aff + 0 
 
    data_type = 'msk'
    msk_joined = concatenate_2D(msk_ref, msk_aff, rc_ref, rc_aff, ref_flag, data_type)

    # adjust the attribute table
    msk_atr = msk_ref_atr
    msk_atr['LENGTH'] = msk_joined.shape[0]
    msk_atr['WIDTH'] = msk_joined.shape[1]

    return msk_joined, msk_atr

def concatenate_geo(inps):
    """concatenate geometry data"""
    data_geo1 = inps.geo_file1[0]
    data_geo2 = inps.geo_file2[0]

    if data_geo1 == 'default':
        data_geo1 = os.path.dirname(inps.patch_files[0]) + '/inputs/' + 'geometryRadar.h5'
    if data_geo2 == 'default':
        data_geo2 = os.path.dirname(inps.patch_files[1]) + '/inputs/' + 'geometryRadar.h5'
        
    geo_full = inps.geo_full[0]

    print('Read the geometry data for the whole region') 
    lat_full = readfile.read(geo_full, datasetName='latitude')[0]
    lon_full = readfile.read(geo_full, datasetName='longitude')[0]
    inc_full = readfile.read(geo_full, datasetName='incidenceAngle')[0]
    azi_full = readfile.read(geo_full, datasetName='azimuthAngle')[0]
    hgt_full = readfile.read(geo_full, datasetName='height')[0]
    
    print('Read the geometry data of the first given dataset') 
    lat1, lat_atr1 = readfile.read(data_geo1, datasetName='latitude')
    lon1, lon_atr1 = readfile.read(data_geo1, datasetName='longitude')
    inc1, inc_atr1 = readfile.read(data_geo1, datasetName='incidenceAngle')
    azi1, azi_atr1 = readfile.read(data_geo1, datasetName='azimuthAngle')
    hgt1, hgt_atr1 = readfile.read(data_geo1, datasetName='height')

    print('Read the geometry data of the second given dataset') 
    lat2, lat_atr2 = readfile.read(data_geo2, datasetName='latitude')
    lon2, lon_atr2 = readfile.read(data_geo2, datasetName='longitude')
    inc2, inc_atr2 = readfile.read(data_geo2, datasetName='incidenceAngle')
    azi2, azi_atr2 = readfile.read(data_geo2, datasetName='azimuthAngle')
    hgt2, hgt_atr2 = readfile.read(data_geo2, datasetName='height')
    
    lat1_flatten = lat1.flatten() # flatten matrix according rows
    lon1_flatten = lon1.flatten()

    lat2_flatten = lat2.flatten()
    lon2_flatten = lon2.flatten()
    
    # calculate the unflatten pattern
    unflatten_trans1 = flatten_trans(lat1)
    unflatten_trans2 = flatten_trans(lat2)

    print('Convert to X-Y coordinate based on the geometry info')
    rc1 = get_position(lat1, lon1, lat_full, lon_full)
    rc2 = get_position(lat2, lon2, lat_full, lon_full)

    # extract the geometry info for the joined region
    row_join_start, row_join_end, col_join_start, col_join_end = design_joined_matrix(rc1, rc2)

    lat_joined = lat_full[row_join_start: row_join_end + 1, col_join_start: col_join_end + 1]
    lon_joined = lon_full[row_join_start: row_join_end + 1, col_join_start: col_join_end + 1]
    inc_joined = inc_full[row_join_start: row_join_end + 1, col_join_start: col_join_end + 1]
    azi_joined = azi_full[row_join_start: row_join_end + 1, col_join_start: col_join_end + 1]
    hgt_joined = hgt_full[row_join_start: row_join_end + 1, col_join_start: col_join_end + 1]
    
    return lat_joined, lon_joined, inc_joined, azi_joined, hgt_joined, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten, unflatten_trans1, unflatten_trans2, rc1, rc2

def write_vel(vel_joined, vel_atr, datatype, inps):
    
    row, colm = vel_joined.shape
    atr_vel = dict()
    atr_vel['WIDTH'] = str(colm)
    atr_vel['LENGTH'] = str(row)
    atr_vel['FILE_TYPE'] = 'velocity'

    #vel_atr['FILE_TYPE'] = 'velocity'
    
    vel_data = dict()
    vel_data['velocity'] = vel_joined

    output_dir = inps.outdir[0]
    outname = inps.out_suffix[0]
    if len(outname) is 0:
        vel_filename = output_dir + '/' + datatype +  '.h5'
    else:
        vel_filename = output_dir + '/' + datatype + '_' + outname + '.h5'

    writefile.write(datasetDict=vel_data, out_file=vel_filename, metadata=atr_vel)

    return 

def write_ts(ts_joined_dataset, ts_atr, date_final, datatype, inps):
    row, colm = ts_joined_dataset['timeseries'].shape[1: ]

    atr_ts = ts_atr
    atr_ts['WIDTH'] = str(colm)
    atr_ts['LENGTH'] = str(row)
    atr_ts['FILE_TYPE'] = 'timeseries'
    #ts_atr['FILE_TYPE'] = 'timeseries'
    
    output_dir = inps.outdir[0]
    outname = inps.out_suffix[0]
    file_name = os.path.basename(inps.patch_files[0])
    if len(outname) is 0:
        #ts_filename = output_dir + '/' + file_name.split('.')[0] +  '.h5'
        ts_filename = output_dir + '/' + datatype +  '.h5'
    else:
        ts_filename = output_dir + '/' + datatype + '_' + outname + '.h5'

    writefile.write(datasetDict=ts_joined_dataset, out_file=ts_filename, metadata=atr_ts)
    
    return

def write_mask(msk_joined, msk_atr, datatype, inps):
    row, colm = msk_joined.shape
    
    # write simple attribution
    atr_msk = dict()
    atr_msk['WIDTH'] = str(colm)
    atr_msk['LENGTH'] = str(row)
    atr_msk['FILE_TYPE'] = 'mask'
    #msk_atr['FILE_TYPE'] = 'mask'
    atr_msk['DATA_TYPE'] = msk_atr['DATA_TYPE']
    if msk_atr['DATA_TYPE'] == 'bool': # the mask is maskTempCoh.h5
        msk_joined = msk_joined > 0

    msk_data = dict()
    msk_data['mask'] = msk_joined

    output_dir = inps.outdir[0]
    outname = inps.out_suffix[0]
    file_name = os.path.basename(inps.patch_files[0])
    if len(outname) is 0:
        msk_filename = output_dir + '/' + datatype +  '.h5'
    else:
        #msk_filename = output_dir + '/' + file_name.split('.')[0] + '_' + outname + '.h5'
        msk_filename = output_dir + '/' + datatype + '_' + outname + '.h5'

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
    outname = inps.out_suffix[0]
    if len(outname) is 0:
        geo_outname = 'geometryRadar'
    else:
        geo_outname = 'geometryRadar_' + outname

    geo_filename = output_dir + '/' + geo_outname + '.h5'

    # write h5 file
    writefile.write(datasetDict=geo_data, out_file=geo_filename, metadata=atr_geo)

    print('Finish!')

def find_the_reference(rc1, rc2):
    """Find the reference data based on geo info"""
    row1_ul = rc1[0][0] # upper left point
    col1_ul = rc1[0][1]
    row1_lr = rc1[3][0] # lower right point
    col1_lr = rc1[3][1]

    row2_ul = rc2[0][0] # upper left point
    col2_ul = rc2[0][1]
    row2_lr = rc2[3][0] # lower right point
    col2_lr = rc2[3][1]

    if row1_ul <= row2_ul and col1_ul <= col2_ul and row1_lr <= row2_lr and col1_lr <= col2_lr:
        return 1
    elif row1_ul > row2_ul and col1_ul < col2_ul and row1_lr > row2_lr and col1_lr < col2_lr:
        return 2
    elif row1_ul < row2_ul and col1_ul > col2_ul and row1_lr < row2_lr and col1_lr > col2_lr:
        return 3
    elif row1_ul > row2_ul and col1_ul > col2_ul and row1_lr > row2_lr and col1_lr > col2_lr:
        return 4
    elif row1_ul > row2_ul and col1_ul < col2_ul and row1_lr > row2_lr and col1_lr > col2_lr:
        return 5
    elif row1_ul < row2_ul and col1_ul < col2_ul and row1_lr < row2_lr and col1_lr > col2_lr:
        return 6
    elif row1_ul < row2_ul and col1_ul < col2_ul and row1_lr > row2_lr and col1_lr < col2_lr:
        return 7
    elif row1_ul < row2_ul and col1_ul > col2_ul and row1_lr > row2_lr and col1_lr > col2_lr:
        return 8

def appoint_ref(rc1, rc2, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten, ref_flag):
    if ref_flag == 2 or ref_flag == 4:
        data_ref = 1
        data_aff = 0
        rc_ref = rc2
        rc_aff = rc1
        lat_ref_flatten = lat2_flatten
        lon_ref_flatten = lon2_flatten
        lat_aff_flatten = lat1_flatten
        lon_aff_flatten = lon1_flatten
        print('The second dataset is the reference image')
    else:
        data_ref = 0
        data_aff = 1
        rc_ref = rc1
        rc_aff = rc2
        lat_ref_flatten = lat1_flatten
        lon_ref_flatten = lon1_flatten
        lat_aff_flatten = lat2_flatten
        lon_aff_flatten = lon2_flatten
        print('The first dataset is the reference image')

    return data_ref, data_aff, rc_ref, rc_aff, lat_ref_flatten, lon_ref_flatten, lat_aff_flatten, lon_aff_flatten

def determine_datatype(inps):
    data = inps.patch_files[0]
    data_atr = readfile.read_attribute(data)
    datatype = data_atr['FILE_TYPE'] 

    return datatype

def list_and_check(inps):
    def state_judge(dataset):
        if os.path.exists(dataset):
            return 'True'
        else:
            return 'False'
    def check_ordering(geo_data):
        lat = readfile.read(geo_data, datasetName='latitude')[0]
        lon = readfile.read(geo_data, datasetName ='longitude')[0]
        if lat[0][0] < lat[-1][0] and lon[0][0] < lon[0][-1]:
            print('{}: Correct lat and lon order'.format(geo_data)) 
        else:
            raise ValueError('Wrong lat/lon order! The lat should increase from north to south. The lon should increase from west to east')
        return
    print('Check and list the attribute dataset to be concatenated:')
    dataset1 = inps.patch_files[0]
    dataset2 = inps.patch_files[1]
    print('---The first dataset is {} and the state of being is {}'.format(dataset1, state_judge(dataset1)))
    print('---The second dataset is {} and the state of being is {}'.format(dataset2, state_judge(dataset2)))

    print('Check and list the geometry dataset used:')
    data_geo1 = inps.geo_file1[0]
    if data_geo1 == 'default':
        data_geo1 = os.path.dirname(inps.patch_file[0]) + './inputs/' + 'geometryRadar.h5'
    data_geo2 = inps.geo_file2[0]
    if data_geo2 == 'default':
        data_geo2 = os.path.dirname(inps.patch_file[1]) + './inputs/' + 'geometryRadar.h5'
    geo_full = inps.geo_full[0]
    print('---The first geometry dataset is {} and the state of being is {}'.format(data_geo1, state_judge(data_geo1)))
    print('---The second geometry dataset is {} and the state of being is {}'.format(data_geo2, state_judge(data_geo2)))
    print('---The geometry dataset of whole region is {} and the state of being is {}'.format(geo_full, state_judge(geo_full)))

    print('Check the latitude and longitude ordering')
    check_ordering(data_geo1)   
    check_ordering(data_geo2)   

    return

def main(iargs=None):
    inps = cmd_line_parse(iargs)   

    if inps.list:
        list_and_check(inps)
        exit()

    datatype = determine_datatype(inps)
    print('Data type is: ', datatype)

    print('Process the geometry info')
    lat_joined, lon_joined, inc_joined, azi_joined, hgt_joined, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten, unflatten_trans1, unflatten_trans2, rc1, rc2 = concatenate_geo(inps)
    if inps.geo_write:
        write_geo(lat_joined, lon_joined, inc_joined, azi_joined, hgt_joined, inps)

    print('Find the relative position between the two datasets')
    ref_flag = find_the_reference(rc1, rc2)

    print('process %s data' % datatype)
    if datatype == 'velocity':
        vel_joined, vel_atr = concatenate_vel(inps, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten, unflatten_trans1, unflatten_trans2, rc1, rc2, ref_flag)
        write_vel(vel_joined, vel_atr, datatype, inps)
    elif datatype == 'timeseries':
        ts_join_dataset, ts_atr, date_final = concatenate_ts(inps, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten, unflatten_trans1, unflatten_trans2, rc1, rc2, ref_flag)
        write_ts(ts_join_dataset, ts_atr, date_final, datatype, inps)
    elif datatype == 'mask':
        msk_joined, msk_atr = concatenate_mask(inps, rc1, rc2, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten, ref_flag)
        write_mask(msk_joined, msk_atr, datatype, inps)
    #elif datatype == 'mask':
    #    msk_joined, msk_atr = concatenate_mask(inps, row_ref, col_ref, row_aff, col_aff, row_a_r, col_a_r, ref_flag)
    #    write_mask(msk_joined, inps)
######################################################################################
if __name__ == '__main__':
    main()
