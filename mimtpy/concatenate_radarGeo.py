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
2. This script concatenates two datasets together. The input files are object datasets, their corresponding geometryRadar.h5 files
3. If a batch concatenation needed, please use the concatnate_patches.py script.
4. For timeseries, datasets should have same reference date

"""

EXAMPLE = """example:

    concatenate_radarGeo.py miaplpy_NE/velocity_msk.h5  miaplpy_NNE/velocity_msk.h5 --geo_file1 miaplpy_NE/inputs/geometryRadar.h5 --geo_file2 miaplpy_NNE/intpus/geometryRadar.h5 --geo_full ./inputs/geometryRadar.h5 --datatype vel --geo_write --output NE_NNE --outdyyir miaplpyBig

    concatenate_radarGeo.py miaplpy_NE/velocity_msk.h5  miaplpy_NNE/velocity_msk.h5 --geo_file1 miaplpy_NE/inputs/geometryRadar.h5 --geo_file2 miaplpy_NNE/intpus/geometryRadar.h5 --geo_full ./inputs/geometryRadar.h5 --datatype vel --output NE_NNE --outdyyir miaplpyBig
    
    concatenate_radarGeo.py miaplpy_NE/timeseries_msk.h5 miaplpy_NNE/timeseries_msk.h5 --geo_file1 miaplpy_NE/inputs/geometryRadar.h5 --geo_file2 miaplpy_NNE/inputs/geometryRadar.h5 --geo_full ./inputs/geometryRadar.h5 --datatype ts --output NE_NNE --outdir ./miaplpyBig/

    concatenate_radarGeo.py miaplpy_NE/maskPS.h5  miaplpy_NNE/maskPS.h5 --geo_file1 miaplpy_NE/inputs/geometryRadar.h5 --geo_file2 miaplpy_NNE/intpus/geometryRadar.h5 --geo_full ./inputs/geometryRadar.h5 --datatype maskPS --output NE_NNE --outdir miaplpyBig

    concatenate_radarGeo.py miaplpy_NE/maskTempCoh.h5  miaplpy_NNE/maskTempCoh.h5 --geo_file1 miaplpy_NE/inputs/geometryRadar.h5 --geo_file2 miaplpy_NNE/intpus/geometryRadar.h5 --geo_full ./inputs/geometryRadar.h5 --datatype maskTC --output NE_NNE --outdyyir miaplpyBig
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Concatenate miaplpy patches',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=NOTE + '\n' + EXAMPLE)

    parser.add_argument('patch_files', nargs='+', type=str, help='two displacement datasets to be concatenated \n')

    parser.add_argument('--geo_file1', nargs=1, type=str, help='geometryRadar file of dataset1. \n')
    
    parser.add_argument('--geo_file2', nargs=1, type=str, help='geometryRadar file of dataset2. \n')
    
    parser.add_argument('--geo_full', nargs=1, type=str, help='geometryRadar file of the full region. \n')
    
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

def geo_position_drop(data_sub, data):
    a2 = np.lib.stride_tricks.as_strided(data, shape=data.shape + data_sub.shape, strides=2*data.strides)
    a2 = a2[:-data_sub.shape[0]+1, :-data_sub.shape[1]+1]
    import pdb
    pdb.set_trace() 
    pos_matrix = (a2 == data_sub).all(axis=(-2, -1))
    #pos_temp = np.array(pos_matrix == True) 
    
    return pos_matrix

def geo_position(lat_sub, lon_sub, lat, lon):
    lat_ul = lat_sub[0][0]
    lon_ul = lon_sub[0][0]

    lat_flag = (lat == lat_ul)
    lon_flag = (lon == lon_ul)

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
    print('The overlay offset is %f' % overlay_offset)

    # adjust the affiliate image
    data2_flatten_adjust = data2_flatten + overlay_offset

    # do concatenation
    #data_joined = np.transpose(np.array([np.hstack((data1_flatten, data2_flatten_adjust))])) # for the flatten data concatenation

    # return data_joined
    return data2_flatten_adjust

def concatenate_2D_geo(geo1_, geo2_, row1, col1, row2, col2, geometry):
    # join the data
    row_sum = geo1_.shape[0] + geo2_.shape[0]
    col_sum = geo1_.shape[1] + geo2_.shape[1]

    # join geo
    geo_joined = np.ones((row_sum, col_sum)) * np.nan
    geo_joined[0: geo1_.shape[0], 0: geo1_.shape[1]] = geo1_
    geo_joined[geo1_.shape[0]: , geo1_.shape[1]: ] = geo2_
    geo_joined[0: geo1_.shape[0], geo1_.shape[1]:] = geometry[row1: row1 + geo1_.shape[0], col1 + geo1_.shape[1]: col1 + geo1_.shape[1] + geo2_.shape[1]]
    geo_joined[geo1_.shape[0]:, 0: geo1_.shape[1]] = geometry[row2: row2 + geo2_.shape[0], col2 - geo1_.shape[1]: col2]

    return geo_joined 

def concatenate_2D(data1_flatten, data2_flatten, unflatten_trans1, unflatten_trans2):
    # transfer to 2D matrix
    data1_ = unflatten_trans1(data1_flatten)
    data2_ = unflatten_trans2(data2_flatten)

    # join the data
    row_sum = data1_.shape[0] + data2_.shape[0]
    col_sum = data1_.shape[1] + data2_.shape[1]

    # join value
    data_joined = np.ones((row_sum, col_sum)) * np.nan
    data_joined[0: data1_.shape[0], 0: data1_.shape[1]] = data1_
    data_joined[data1_.shape[0]: , data1_.shape[1]: ] = data2_

    return data_joined

def concatenate_vel(inps, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten, unflatten_trans1, unflatten_trans2):
    data1 = inps.patch_files[0] # reference image
    data2 = inps.patch_files[1] # affiliate image
   
    print('Read the first dataset') 
    vel1, vel1_atr = readfile.read(data1, datasetName='velocity')
    vel1_flatten = vel1.flatten()
    
    print('Read the second dataset') 
    vel2, vel2_atr = readfile.read(data2, datasetName='velocity')
    vel2_flatten = vel2.flatten()

    #vel_joined = concatenate_process(vel1_flatten, vel2_flatten, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten)
    vel2_flatten_adjust = concatenate_process(vel1_flatten, vel2_flatten, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten)

    # generate 2D concatenation results
    vel_joined = concatenate_2D(vel1_flatten, vel2_flatten_adjust, unflatten_trans1, unflatten_trans2)

    # adjust the attribute table
    vel_atr = vel1_atr
    vel_atr['LENGTH'] = vel_joined.shape[0]
    vel_atr['WIDTH'] = vel_joined.shape[1]
    
    return vel_joined, vel_atr

def concatenate_ts(inps, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten, unflatten_trans1, unflatten_trans2):
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
    #pnum = len(lat1_flatten) + len(lat2_flatten)

    ts_join_dataset = dict()
    row_sum = ts1.shape[1] + ts2.shape[1]
    col_sum = ts1.shape[2] + ts2.shape[2]
    ts_join = np.empty(shape=(join_dim, row_sum, col_sum), dtype=float)
    #ts_join = np.empty(shape=(join_dim, pnum, 1), dtype=float)
    
    # do concatenation
    i = 0
    for date1, date2 in zip(Date1, Date2):
        print('Process displacement data of date %s' % date1)
        dis1 = readfile.read(data1, datasetName=date1)[0]
        dis1_flatten = dis1.flatten()
        dis2 = readfile.read(data2, datasetName=date2)[0]
        dis2_flatten = dis2.flatten()

        #ts_join[i, :, :] = concatenate_process(dis1_flatten, dis2_flatten, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten)i
        dis2_flatten_adjust = concatenate_process(dis1_flatten, dis2_flatten, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten)
        # generate 2D concatenation results
        ts_join[i, :, :] = concatenate_2D(dis1_flatten, dis2_flatten_adjust, unflatten_trans1, unflatten_trans2)
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

    #msk1_flatten = msk1.flatten('F') # flatten matrix according colms
    #msk2_flatten = msk2.flatten('F') # flatten matrix according colms

    # do concatenation
    #msk_joined = np.transpose(np.array([np.hstack((msk1_flatten, msk2_flatten))]))
    # join the data
    row_sum = msk1.shape[0] + msk2.shape[0]
    col_sum = msk1.shape[1] + msk2.shape[1]

    # join value
    msk_joined = np.zeros((row_sum, col_sum))
    msk_joined[0: msk1.shape[0], 0: msk1.shape[1]] = msk1
    msk_joined[msk1.shape[0]: , msk1.shape[1]: ] = msk2

    return msk_joined, msk1_atr

def concatenate_geo(inps):
    """concatenate geometry data"""
    data_geo1 = inps.geo_file1[0]
    data_geo2 = inps.geo_file2[0]
    geo_full = inps.geo_full[0]

    print('Read the geometry data for the full region') 
    lat_full = readfile.read(geo_full, datasetName='latitude')[0]
    lon_full = readfile.read(geo_full, datasetName='longitude')[0]
    inc_full = readfile.read(geo_full, datasetName='incidenceAngle')[0]
    azi_full = readfile.read(geo_full, datasetName='azimuthAngle')[0]
    hgt_full = readfile.read(geo_full, datasetName='height')[0]
    
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
    
    lat1_flatten = lat1.flatten() # flatten matrix according rows
    lon1_flatten = lon1.flatten()

    lat2_flatten = lat2.flatten()
    lon2_flatten = lon2.flatten()
    
    # calculate the unflatten pattern
    unflatten_trans1 = flatten_trans(lat1)
    unflatten_trans2 = flatten_trans(lat2)

    # do concatenation
    #lat_joined = np.transpose(np.array([np.hstack((lat1_flatten, lat2_flatten))]))
    #lon_joined = np.transpose(np.array([np.hstack((lon1_flatten, lon2_flatten))]))
    #inc_joined = np.transpose(np.array([np.hstack((inc1_flatten, inc2_flatten))]))
    #azi_joined = np.transpose(np.array([np.hstack((azi1_flatten, azi2_flatten))]))
    #hgt_joined = np.transpose(np.array([np.hstack((hgt1_flatten, hgt2_flatten))]))

    # find the position of subset lat/lon in the Geometry data
    #pos_la1 = geo_position(lat1, lat_full)
    #pos_lo1 = geo_position(lon1, lon_full)
    #pos1_ = pos_la1 * pos_lo1
    #row1 = np.array(pos1_ == True)[0][0]
    #col1 = np.array(pos1_ == True)[1][0]

    #pos_la2 = geo_position(lat2, lat_full)
    #pos_lo2 = geo_position(lon2, lon_full)
    #pos2_ = pos_la2 * pos_lo2
    #row2 = np.array(pos2_ == True)[0][0]
    #col2 = np.array(pos2_ == True)[1][0]

    row1, col1 = geo_position(lat1, lon1, lat_full, lon_full)
    row2, col2 = geo_position(lat2, lon2, lat_full, lon_full)
    
    lat_joined = concatenate_2D_geo(lat1, lat2, row1, col1, row2, col2, lat_full)
    lon_joined = concatenate_2D_geo(lon1, lon2, row1, col1, row2, col2, lon_full)
    inc_joined = concatenate_2D_geo(inc1, inc2, row1, col1, row2, col2, inc_full)
    azi_joined = concatenate_2D_geo(azi1, azi2, row1, col1, row2, col2, azi_full)
    hgt_joined = concatenate_2D_geo(hgt1, hgt2, row1, col1, row2, col2, hgt_full)

    return lat_joined, lon_joined, inc_joined, azi_joined, hgt_joined, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten, unflatten_trans1, unflatten_trans2

def write_vel(vel_joined, vel_atr, inps):
    
    row, colm = vel_joined.shape

    atr_vel = dict()
    atr_vel['WIDTH'] = str(colm)
    atr_vel['LENGTH'] = str(row)
    #atr_vel['REF_X'] = str(1 - 1)
    #atr_vel['REF_Y'] = str((int(vel_atr['REF_Y']) + 1) + ((int(vel_atr['REF_X']) + 1) - 1) * int(vel_atr['LENGTH']) - 1)
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
    #atr_ts['REF_LAT'] = str(1 - 1)
    #atr_ts['REF_LON'] = str((int(ts_atr['REF_Y']) + 1) + ((int(ts_atr['REF_X']) + 1) - 1) * int(ts_atr['LENGTH']) - 1)
    atr_ts['FILE_TYPE'] = 'timeseries'
    
    output_dir = inps.outdir[0]
    outname = inps.output[0]

    ts_filename = output_dir + inps.datatype[0] + '_' + outname + '.h5'
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
    lat_joined, lon_joined, inc_joined, azi_joined, hgt_joined, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten, unflatten_trans1, unflatten_trans2 = concatenate_geo(inps)
    if inps.geo_write:
        write_geo(lat_joined, lon_joined, inc_joined, azi_joined, hgt_joined, inps)

    print('process %s data' % inps.datatype[0])
    if inps.datatype[0] == 'vel':
        vel_joined, vel_atr = concatenate_vel(inps, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten, unflatten_trans1, unflatten_trans2)
        write_vel(vel_joined, vel_atr, inps)
    elif inps.datatype[0] == 'ts':
        ts_join_dataset, ts_atr, date_final = concatenate_ts(inps, lat1_flatten, lon1_flatten, lat2_flatten, lon2_flatten, unflatten_trans1, unflatten_trans2)
        write_ts(ts_join_dataset, ts_atr, date_final, inps)
    elif inps.datatype[0] == 'maskPS':
        msk_joined, msk_atr = concatenate_mask(inps)
        write_mask(msk_joined, inps)
    elif inps.datatype[0] == 'maskTC':
        msk_joined, msk_atr = concatenate_mask(inps)
        write_mask(msk_joined, inps)
        

######################################################################################
if __name__ == '__main__':
    main()
