#!/usr/bin/env python3
######################################################################################################
# Program is used for concatenating two tracks based on the overlay region                           #
# Author: Lv Xiaoran                                                                                 #
# Created: January  2020                                                                             #
######################################################################################################

import os
import sys
import argparse
import copy
import numpy as np
import json
import math
import re
import h5py

import matplotlib.pyplot as plt
from mintpy.utils import readfile, writefile, ptime
from mintpy.objects import (
    giantIfgramStack, 
    giantTimeseries, 
    ifgramStack, 
    timeseries, 
    HDFEOS,
)

from mimtpy.objects.profiles import Profile
import mimtpy.objects.profiles as profiles

######################################################################################
EXAMPLE = """example:
    Note:
    if you want to make sure the REF is also correct, please make the dominant image is at the east of the affiliate image!
    concatenate_offset.py $SCRATCHDIR/BogdSenDT33/mintpy/geo_velocity.h5 $SCRATCHDIR/BogdSenDT106/mintpy/geo_velocity.h5  --rewrite_affiliate --outdir /data/lxr/insarlab/SCRATCHDIR/BalochistanSenDT/cumulative/

    concatenate_offset.py $SCRATCHDIR/BogdSenDT33/mintpy/geo_velocity.h5 $SCRATCHDIR/BogdSenDT106/mintpy/geo_velocity.h5 --output mosaic --outdir /data/lxr/insarlab/SCRATCHDIR/BalochistanSenDT/cumulative/

    concatenate_offset.py $SCRATCHDIR/BogdSenDT33/mintpy/geo_velocity.h5 $SCRATCHDIR/BogdSenDT106/mintpy/geo_velocity.h5 --output Bogd_mosaic --plotpair --pro_num 1 --ll 99.5 45 -g $SCRATCHDIR/BogdSenDT33/mintpy/inputs/geo_geometry.h5 --outdir ./
    
    concatenate_offset.py $SCRATCHDIR/BogdSenDT33/mintpy/geo_velocity.h5 $SCRATCHDIR/BogdSenDT106/mintpy/geo_velocity.h5 --output Bogd_mosaic --plotpair --pro_num 3 -g $SCRATCHDIR/BogdSenDT33/mintpy/inputs/geo_geometry.h5 --outdir ./

    concatenate_offset.py $SCRATCHDIR/BogdSenDT33/mintpy/inputs/geo_geometryRadar.h5 $SCRATCHDIR/BogdSenDT106/mintpy/inputs/geo_geometryRadar.h5 -ot ave --output Bogd_mosaic --outdir ./

    concatenate_offset.py ../../SpainSenAT147/mintpy/timeseries/timeseries_demErr_ERA5.h5 ../../SpainSenAT74/mintpy/timeseries/timeseries_demErr_ERA5.h5 --output timeseries_AT147_74 --outdir ./

    concatenate_offset.py ../../KokoxiliBigSenAT114/velocity/velocity_lat30_41.h5 ../../KokoxiliBigSenAT41/velocity/velocity_lat30_41.h5 --output steeks --outdir ./ --plotpair -g ../../KokoxiliBigSenAT114/geometry/geo_lat30_41.h5 --pro_num 5

"""

def create_parser():
    parser = argparse.ArgumentParser(description='Concatenating two different tracks based on the overlay region',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('dominant', nargs=1, help='dominant track. \n')

    parser.add_argument('affiliate', nargs=1, help='affiliate track. \n ')

    parser.add_argument('-b', '--bbox', dest='SNWE', type=float, nargs=4, metavar=('S', 'N', 'W', 'E'),
                        help='Bounding box of area to be geocoded.\n' +
                        'Include the uppler left corner of the first pixel' +
                        '    and the lower right corner of the last pixel')
   
    parser.add_argument('--rewrite_affiliate', action='store_true', default=False, help='whether rewrite affiliate *.h5 file after add the offset.\n')

    plotpair_opt = parser.add_argument_group(title='whether plot profiles for overlay regions of both tracks')

    plotpair_opt.add_argument('--plotpair', action='store_true', default=False, help='whether plot profiles for overlay regions of both tracks. \n')

    plotpair_opt.add_argument('--pro_num', dest='pro_num', type=int, help='profile numbers\n')

    plotpair_opt.add_argument('--ll', dest='LONLAT', type=float, nargs=2,
                              help='longitude and latitude of point that you want to profile pass through.'
                                    'use this option when pro_num = 1.')
    plotpair_opt.add_argument('-g', '--geometry', dest='geometry', type=str, help='geoemtry data of dominant track.')

    mosaic = parser.add_argument_group(title='mosaic options')
    #mosaic.add_argument('--mosaic', action='store_true', default=False, help='whether mosaic two track data.') 
    mosaic.add_argument('-ot', '--overlapTreatment', dest="overlapTreatment", nargs='?', type=str, 
                         help='process method of incidence and azimuth angle for overlapping region,'
                         'average: using average value of dominant and affiliate track for the overlapping region;'
                         'dominant: using dominant swath value for the overlapping region;'
                         'affiliate: using affiliate swath value for the overlapping region')   
    
    mosaic.add_argument('--output',dest='output',nargs=1,help='output name')
 
    parser.add_argument('--outdir',dest='outdir',type=str,nargs=1,help='outdir')    

    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps    

def get_bounding_box(atr):
    """get lon/lat range
    lat0,lon0- starting latitude/longitude (first row/column)
    lat1,lon1- ending latitude/longitude (last row/column)
    """
    length,width = int(atr['LENGTH']),int(atr['WIDTH'])
    
    if 'Y_FIRST' in atr.keys():
        # geo coordinates
        lat0 = float(atr['Y_FIRST'])
        lon0 = float(atr['X_FIRST'])
        lat_step = float(atr['Y_STEP'])
        lon_step = float(atr['X_STEP'])
        lat1 = lat0 + (length - 1) * lat_step
        lon1 = lon0 + (width - 1) * lon_step
    else:
        # radar coordinates
        lats = [float(atr['LAT_REF{}'.format(i)]) for i in [1,2,3,4]]
        lons = [float(atr['LON_REF{}'.format(i)]) for i in [1,2,3,4]]
        lat0 = np.mean(lats[0:2])
        lat1 = np.mean(lats[2:4])
        lon0 = np.mean(lons[0:3:2])
        lon1 = np.mean(lons[1:4:2])
    return lat0, lon0, lat1, lon1

def calculate_rc(lat0,lon0,lat_n,lon_n,lat_step,lon_step):
    row = int((lat_n - lat0) / lat_step + 0.5)
    colm = int((lon_n - lon0) / lon_step + 0.5)
     
    return row, colm

def calculate_roco_overlay(lat0,lon0,lat1,lon1,lat_step,lon_step):
    rows = int((lat1 - lat0) / lat_step + 0.5) + 1
    colms = int((lon1 - lon0) / lon_step + 0.5) + 1

    return rows, colms

def max(a,b):
    if a > b:
        return a
    else:
        return b

def min(a,b):
    if a < b:
        return a
    else:
        return b

def calculate_overlay(inps,m_atr,m_data,s_atr,s_data,typeflag):
    parser = create_parser()
    
    # calculate overlay region
    # lat lon range of dominant and affiliate image
    m_lat0, m_lon0, m_lat1, m_lon1 = get_bounding_box(m_atr)
    
    s_lat0, s_lon0, s_lat1, s_lon1 = get_bounding_box(s_atr)

    #get lalo of overlay region
    over_lat0_tmp = min(m_lat0,s_lat0)
    over_lon0_tmp = max(m_lon0,s_lon0)
    over_lat1_tmp = max(m_lat1,s_lat1)
    over_lon1_tmp = min(m_lon1,s_lon1)

    #intersection of overlay region and user bounding box
    if inps.SNWE == None:
        over_lat0 = over_lat0_tmp
        over_lon0 = over_lon0_tmp
        over_lat1 = over_lat1_tmp
        over_lon1 = over_lon1_tmp
    else:
        user_lat0 = float(inps.SNWE[1])
        user_lon0 = float(inps.SNWE[2])
        user_lat1 = float(inps.SNWE[0])
        user_lon1 = float(inps.SNWE[3])
        if user_lat0 < user_lat1:
            parser.print_usage()
            raise Exception('input bounding box error! Wrong latitude order!')
        elif user_lon0 > user_lon1:
            parser.print_usage()
            raise Exception('input bounding box error! Wrong longitude order!')
        else:
            over_lat0 = min(over_lat0_tmp,user_lat0)
            over_lon0 = max(over_lon0_tmp,user_lon0)
            over_lat1 = max(over_lat1_tmp,user_lat1)
            over_lon1 = min(over_lon1_tmp,user_lon1)
    
    # get row/colm number of overlay region
    overlay_rows,overlay_colms = calculate_roco_overlay(over_lat0,over_lon0,over_lat1,over_lon1,float(m_atr['Y_STEP']),float(m_atr['X_STEP']))
    
    #get the row/column position of overlay region
    m_row0, m_colm0 = calculate_rc(m_lat0,m_lon0,over_lat0,over_lon0,float(m_atr['Y_STEP']),float(m_atr['X_STEP']))
    #m_row1, m_colm1 = calculate_rc(m_lat0,m_lon0,over_lat1,over_lon1,float(m_atr['Y_STEP']),float(m_atr['X_STEP']))
    s_row0, s_colm0 = calculate_rc(s_lat0,s_lon0,over_lat0,over_lon0,float(s_atr['Y_STEP']),float(s_atr['X_STEP']))
    #s_row1, s_colm1 = calculate_rc(s_lat0,s_lon0,over_lat1,over_lon1,float(s_atr['Y_STEP']),float(s_atr['X_STEP']))

    if typeflag == 0 or typeflag == 2:
        # calculate overlay region in dominant and affiliate
        m_overlay = m_data[m_row0 : m_row0 + overlay_rows,m_colm0 : m_colm0 + overlay_colms]
        s_overlay = s_data[s_row0 : s_row0 + overlay_rows,s_colm0 : s_colm0 + overlay_colms]

        # calculate offset between dominant and affiliate
        offset_overlay = m_overlay - s_overlay
        offset = np.nanmedian(offset_overlay)
        if 'Y_FIRST' in m_atr.keys():
            ref_lat = m_atr['REF_LAT']
            ref_lon = m_atr['REF_LON']
            ref_x = m_atr['REF_X']
            ref_y = m_atr['REF_Y']
        else:
            ref_x = m_atr['REF_X']
            ref_y = m_atr['REF_Y']
        
        #print('The average offset is : %f \n' % offset)
       
        #if inps.plotpair:
        #    out_dir = inps.outdir[0]
        #    angle = inps.azimuth
        #    overlay_profile(angle, m_overlay, s_overlay, m_name, s_name, out_dir)

        return offset, m_row0, m_colm0, s_row0, s_colm0, over_lat0, over_lon0, overlay_rows, overlay_colms
    else:
        return m_row0, m_colm0, s_row0, s_colm0, over_lat0, over_lon0, overlay_rows, overlay_colms

def rewrite_affiliate(inps,offset,s_atr,s_data):

    s_data_offset = s_data + offset
    #if inps.output == None: 
    out_file = '{}{}{}'.format(os.path.split("".join(inps.affiliate))[1].split('.')[0], '_offset.', os.path.split("".join(inps.affiliate))[1].split('.')[1])
    #else:
    #   out_file = '{}{}{}'.format("".join(inps.output),'.',os.path.split("".join(inps.affiliate))[1].split('.')[1])

    if inps.outdir != None:
        out_dir = inps.outdir[0]

    out_dir_file = out_dir + out_file
    
    if inps.rewrite_affiliate:
        #print('Writing affiliate data %s' % "".join(inps.affiliate))
        writefile.write(s_data_offset, out_file=out_dir_file, metadata=s_atr)

    return s_data_offset

def sum_matrix_nan(matrix1,matrix2):
    rows,colms = matrix1.shape
    matrix_sum = np.zeros((rows,colms),dtype=np.float32) * np.nan
    for row in range(rows):
        for colm in range(colms):
            if matrix1[row,colm] == np.nan and matrix2[row,colm] == np.nan:
                matrix_sum[row,colm] = np.sum([matrix1[row,colm],matrix2[row,colm]])
            else:
                matrix_sum[row,colm] = np.nansum([matrix1[row,colm],matrix2[row,colm]])
    return matrix_sum
 
def mosaic_tracks(inps,m_atr,m_data,s_atr,s_data_offset,m_row0,m_colm0,s_row0,s_colm0,over_lat0,over_lon0,overlay_rows,overlay_colms,dslice=None):
    """mosaic two tracks"""
    mm_atr = copy.deepcopy(m_atr) 
    m_rows,m_colms = m_data.shape
    
    s_data = s_data_offset
    s_rows,s_colms = s_data.shape
    
    # get dominant and affiliate images overlay region data
    m_overlay = m_data[m_row0 : m_row0 + overlay_rows,m_colm0 : m_colm0 + overlay_colms]
    s_overlay = s_data[s_row0 : s_row0 + overlay_rows,s_colm0 : s_colm0 + overlay_colms]

    # calculate the mean difference between dominant and affiliate images overlay region
    ms_difference = np.nanmean(np.abs(m_overlay - s_overlay))
    #print('the difference between dominant and affiliate overlay region is %f' % ms_difference)
    
    # calculate values in overlay region
    mosaic_freq = np.ones((overlay_rows,overlay_colms),dtype = np.float32)
    m_over_pos = np.isnan(m_overlay)
    s_over_pos = np.isnan(s_overlay)
    # change np.nan to zero
    m_overlay[m_over_pos] = 0
    s_overlay[s_over_pos] = 0

    # record the common np.nan in m and s(np.nan in mosaic_freq); one np.nan in m or s(1 in mosaic_freq);none np.nan in m and s(2 in mosaic_freq)
    both_nan = m_over_pos & s_over_pos
    none_nan = m_over_pos | s_over_pos
    mosaic_freq[both_nan] = np.nan
    mosaic_freq[~none_nan] = 2
    
    mosaic_overlay_tmp = (m_overlay + s_overlay)
   
    if dslice:
        if inps.overlapTreatment == 'dominant':
            mosaic_overlay_tmp[~none_nan] = m_overlay[~none_nan]
            mosaic_freq[~none_nan] = 1        
        elif inps.overlapTreatment == 'affiliate':
            mosaic_overlay_tmp[~none_nan] = s_overlay[~none_nan]
            mosaic_freq[~none_nan] = 1        
    
    # calculated final mosaic_overlay
    mosaic_overlay = mosaic_overlay_tmp / mosaic_freq
    
    # generate mosaic dataset
    mosaic_rows = m_rows + s_rows - overlay_rows
    mosaic_colms = m_colms + s_colms - overlay_colms

    # store mosaic
    mosaic_data = np.zeros((mosaic_rows,mosaic_colms),dtype=np.float32) * np.nan 
    #mosaic_freq = np.zeros((mosaic_rows,mosaic_colms),dtype=np.float32)

    # lat lon range of dominant and affiliate image
    m_lat0, m_lon0, m_lat1, m_lon1 = get_bounding_box(m_atr)
    
    s_lat0, s_lon0, s_lat1, s_lon1 = get_bounding_box(s_atr)
    
    # get the lat/lon range of mosaic region
    mosaic_lat0 = max(m_lat0,s_lat0)
    mosaic_lon0 = min(m_lon0,s_lon0)
    mosaic_lat1 = min(m_lat1,s_lat1)
    mosaic_lon1 = max(m_lon1,s_lon1)

    mosaic_lat_step = float(m_atr['Y_STEP'])
    mosaic_lon_step = float(m_atr['X_STEP'])
  
    # calculate m_data offset to mosaic
    m_row_loc0, m_colm_loc0 = calculate_rc(mosaic_lat0,mosaic_lon0,m_lat0,m_lon0,mosaic_lat_step,mosaic_lon_step)
    
    # caculate s_data offset to mosaic
    s_row_loc0,s_colm_loc0 = calculate_rc(mosaic_lat0,mosaic_lon0,s_lat0,s_lon0,mosaic_lat_step,mosaic_lon_step)
    
    # calcualte overlay offset to mosaic
    overlay_row_loc0,overlay_colm_loc0 = calculate_rc(mosaic_lat0,mosaic_lon0,over_lat0,over_lon0,mosaic_lat_step,mosaic_lon_step)
    
    # mosaic data
    mosaic_data[m_row_loc0:m_row_loc0+m_rows,m_colm_loc0:m_colm_loc0+m_colms] = m_data
    mosaic_data[s_row_loc0:s_row_loc0+s_rows,s_colm_loc0:s_colm_loc0+s_colms] = s_data
    mosaic_data[overlay_row_loc0:overlay_row_loc0+overlay_rows,overlay_colm_loc0:overlay_colm_loc0+overlay_colms] = mosaic_overlay 
    # sum freq
    #mosaic_freq[m_row_loc0:m_row_loc0+m_rows,m_colm_loc0:m_colm_loc0+m_colms] += 1
    #mosaic_freq[s_row_loc0:s_row_loc0+s_rows,s_colm_loc0:s_colm_loc0+s_colms] += 1
    #mosaic_data = mosaic_data / mosaic_freq

    # write the mosaic data
    mosaic_atr = mm_atr
    mosaic_atr['LENGTH'] = mosaic_rows
    mosaic_atr['WIDTH'] = mosaic_colms
    mosaic_atr['X_FIRST'] = mosaic_lon0
    mosaic_atr['Y_FIRST'] = mosaic_lat0
    mosaic_atr['Y_FIRST'] = mosaic_lat0
    mosaic_atr['Y_FIRST'] = mosaic_lat0
    mosaic_atr['Y_FIRST'] = mosaic_lat0
    if m_atr['ORBIT_DIRECTION'] == 'DESCENDING':
        mosaic_atr['LAT_REF1'] = m_atr['LAT_REF1']
        mosaic_atr['LON_REF1'] = m_atr['LON_REF1']
        mosaic_atr['LAT_REF3'] = m_atr['LAT_REF3']
        mosaic_atr['LON_REF3'] = m_atr['LON_REF3']
        mosaic_atr['LAT_REF2'] = s_atr['LAT_REF2']
        mosaic_atr['LON_REF2'] = s_atr['LON_REF2']
        mosaic_atr['LAT_REF4'] = s_atr['LAT_REF4']
        mosaic_atr['LON_REF4'] = s_atr['LON_REF4']
    else:
        mosaic_atr['LAT_REF1'] = s_atr['LAT_REF1']
        mosaic_atr['LON_REF1'] = s_atr['LON_REF1']
        mosaic_atr['LAT_REF3'] = s_atr['LAT_REF3']
        mosaic_atr['LON_REF3'] = s_atr['LON_REF3']
        mosaic_atr['LAT_REF2'] = m_atr['LAT_REF2']
        mosaic_atr['LON_REF2'] = m_atr['LON_REF2']
        mosaic_atr['LAT_REF4'] = m_atr['LAT_REF4']
        mosaic_atr['LON_REF4'] = m_atr['LON_REF4']
        
    if inps.plotpair:
        outdir = inps.outdir[0]
        # read geometry data for elevation matrix
        geometry_file = inps.geometry
        topography = readfile.read(geometry_file, datasetName='height')[0]
        overlap_dem = topography[m_row0 : m_row0 + overlay_rows,m_colm0 : m_colm0 + overlay_colms]

        # calculate over_lat1 and over_lon1
        over_lat1 = over_lat0 + (overlay_rows - 1) * mosaic_lat_step
        over_lon1 = over_lon0 + (overlay_colms - 1) * mosaic_lon_step
        # calculate angle for the track
        polygon = m_atr['scene_footprint']
        lonlats = re.findall(r'([\d+\.]+)',polygon)
        # lat/lon of the upper right point
        lon_ur = float(lonlats[0])
        lat_ur = float(lonlats[1])
        # lat/lon of the lower right point
        lon_lr = float(lonlats[2])
        lat_lr = float(lonlats[3])
        # azimuth angle for the track (clockwise from the North)
        angle_rad = math.atan((lon_lr - lon_ur) / (lat_lr - lat_ur))
        if inps.pro_num == 1:
            point_lon = float(inps.LONLAT[0])
            point_lat = float(inps.LONLAT[1])
            
            pro_obj = Profile(1, angle_rad, point_lon, point_lat, m_overlay, s_overlay, over_lat0, over_lon0, m_atr, s_atr, outdir)
            pro_obj.profile_extract()
            # dem value along the profile
            dem_profile = overlap_dem[pro_obj.row_no, pro_obj.colm_no]
            profiles.profile_plot(pro_obj.lon_start, pro_obj.lat_start, pro_obj.lon_end, pro_obj.lat_end, pro_obj.m_profile, pro_obj.s_profile, pro_obj.m_name, pro_obj.s_name, dem_profile, outdir)
        else:
            pro_catalog = profiles.search_profiles(inps.pro_num, over_lat0, over_lon0, over_lat1, over_lon1, m_atr, s_atr)         
            profile_dict_list = []
            profile_dem_list = []
            for pro_NO in pro_catalog[:,0]:
                profile_dict = dict()
                pro_obj = Profile(int(pro_NO), angle_rad, pro_catalog[int(pro_NO)-1,1], pro_catalog[int(pro_NO)-1,2], m_overlay, s_overlay, over_lat0, over_lon0, m_atr, s_atr, outdir)        
                pro_obj.profile_extract()
                profile_dict['NO'] = int(pro_NO)
                profile_dict['p_start'] = np.array([pro_obj.lon_start, pro_obj.lat_start])
                profile_dict['p_end'] = np.array([pro_obj.lon_end, pro_obj.lat_end])
                profile_dict['m_data'] = pro_obj.m_profile
                profile_dict['s_data'] = pro_obj.s_profile
                profile_dict_list.append(profile_dict)
                
                profile_dem = dict()
                profile_dem['NO'] = int(pro_NO)
                profile_dem['value'] = overlap_dem[pro_obj.row_no, pro_obj.colm_no]
                profile_dem_list.append(profile_dem)

                m_name = pro_obj.m_name
                s_name = pro_obj.s_name

            # process multiprofiles
            profile_dict_final, profile_dem_final = profiles.profile_average(inps.pro_num, profile_dict_list, profile_dem_list)
            if mm_atr['FILE_TYPE'] == 'velocity':
                filetype = 'velocity'
            else:
                filetype = 'displacement'
            profiles.profiles_plot(profile_dict_final, profile_dem_final, m_name, s_name, filetype, outdir)             
            #print('the mean difference between dominant and affiliate data is %f' % profile_dict_final[-1]['data'])

    return mosaic_data, mosaic_atr

def write_mosaic(inps, mosaic_data, mosaic_atr):
    """write mosaic data"""
    if inps.output == None: 
        out_file = '{}{}{}'.format(os.path.split("".join(inps.affiliate))[1].split('.')[0], '_mosaic.', os.path.split("".join(inps.affiliate))[1].split('.')[1])
        #out_file = '{}{}{}'.format(os.path.split("".join(inps.affiliate))[1].split('.')[0], '_mosaic.', 'unw')
    else:
        out_file = '{}{}{}'.format("".join(inps.output),'.',os.path.split("".join(inps.affiliate))[1].split('.')[1])

    if inps.outdir != None:
        out_dir = inps.outdir[0]

    out_dir_file = out_dir + out_file

    #print('writing mosaic files:\n')
    writefile.write(mosaic_data, out_file=out_dir_file, metadata=mosaic_atr)

def judge_data_datasets(m_atr):
    """judage dominant/affiliate image is data(velocity or .unw) or dataset (geometry)"""
    file_type = m_atr['FILE_TYPE']
    if file_type == 'geometry':
        typeflag = 1
    elif file_type == 'timeseries':
        typeflag = 2
    elif file_type.find('coherence') != -1:
        typeflag = 3
    else:
        typeflag = 0

    return typeflag

def date_match(m_dateList, s_dateList, m_dim, s_dim, m_bperp, s_bperp):
    """match the date in dominant and affiliate timeseries dataset"""
    m_datevector = ptime.date_list2vector(m_dateList)[1]
    s_datevector = ptime.date_list2vector(s_dateList)[1]

    m_Date = []
    s_Date = []
    bperp = []

    if m_dim <= s_dim:
        for i in np.arange(m_dim):
            date = m_datevector[i]
            date_sub = np.abs(np.array(s_datevector) - date)
            # select object date based on the difference between two dates is less than 7 days [0.02 in vector]
            date_obj = date_sub[date_sub < 0.02]
            if date_obj.size != 0:
                m_Date.append(m_dateList[i])
                bperp.append(m_bperp[i])
                s_Date.append(np.array(s_dateList)[date_sub < 0.02][0])
        date_final = m_Date
    else:
        for i in np.arange(s_dim):
            date = s_datevector[i]
            date_sub = np.abs(np.array(m_datevector) - date)
            # select object date based on the difference between two dates is less than 7 days [0.02 in vector]
            date_obj = date_sub[date_sub < 0.02]
            if date_obj.size != 0:
                s_Date.append(s_dateList[i])
                bperp.append(s_bperp[i])
                m_Date.append(np.array(m_dateList)[date_sub < 0.02][0])
        date_final = s_Date
    return date_final, m_Date, s_Date, bperp
        
######################################################################################
def main(iargs=None):
    """simulate LOS displacement"""
    inps = cmd_line_parse(iargs)

    # dominant and affiliate data attribute
    m_atr = readfile.read_attribute("".join(inps.dominant))
    s_atr = readfile.read_attribute("".join(inps.affiliate))

    # judge file type
    typeflag = judge_data_datasets(m_atr)
    
    if typeflag == 0:
        # dominant data
        m_data = readfile.read("".join(inps.dominant))[0]
        # affiliate data
        s_data = readfile.read("".join(inps.affiliate))[0]
         
        # calculated offset and mosaic two tracks
        offset,m_row0,m_colm0,s_row0,s_colm0,over_lat0,over_lon0,overlay_rows, overlay_colms = calculate_overlay(inps,m_atr,m_data,s_atr,s_data,typeflag)
        s_data_offset = rewrite_affiliate(inps,offset,s_atr,s_data)
            
        #print('prepare mosaicing:\n')
        mosaic_data, mosaic_atr = mosaic_tracks(inps,m_atr,m_data,s_atr,s_data_offset,m_row0,m_colm0,s_row0,s_colm0,over_lat0,over_lon0,overlay_rows,overlay_colms)
        write_mosaic(inps,mosaic_data, mosaic_atr)
            
    elif typeflag == 1:
        #dataslice = ['azimuthAngle', 'height', 'incidenceAngle', 'latitude', 'longitude', 'shadowMask', 'slantRangeDistance']
        dataslice = ['azimuthAngle', 'height', 'incidenceAngle', 'latitude', 'longitude', 'slantRangeDistance']
        m_file = "".join(inps.dominant)
        s_file = "".join(inps.affiliate)
        
        mosaic_dataset = dict()
        
        for dslice in dataslice:
            m_data = readfile.read(m_file, datasetName=dslice)[0]
            s_data = readfile.read(s_file, datasetName=dslice)[0]
            # calculated overlay region
            m_row0,m_colm0,s_row0,s_colm0,over_lat0,over_lon0,overlay_rows, overlay_colms = calculate_overlay(inps,m_atr,m_data,s_atr,s_data,typeflag)
                
            #print('prepare mosaicing for: %s\n' % dslice)
            if dslice == 'incidenceAngle' or dslice == 'azimuthAngle':
                mosaic_data, mosaic_atr = mosaic_tracks(inps,m_atr,m_data,s_atr,s_data,m_row0,m_colm0,s_row0,s_colm0,over_lat0,over_lon0,overlay_rows,overlay_colms,dslice=dslice)
            else:
                mosaic_data, mosaic_atr = mosaic_tracks(inps,m_atr,m_data,s_atr,s_data,m_row0,m_colm0,s_row0,s_colm0,over_lat0,over_lon0,overlay_rows,overlay_colms)
            mosaic_dataset[dslice] = mosaic_data 
            #print('finish mosaic %s' % dslice)       
        write_mosaic(inps,mosaic_dataset, mosaic_atr)
    elif typeflag == 2:
        m_file = "".join(inps.dominant)
        s_file = "".join(inps.affiliate)

        m_data = readfile.read(m_file, datasetName='timeseries')[0]
        s_data = readfile.read(s_file, datasetName='timeseries')[0]
        
        m_bperp_date = h5py.File(m_file,'r')
        m_bperp = m_bperp_date['bperp']
        #m_date = m_bperp_date['date']
        m_dateList = timeseries(m_file).get_date_list()

        s_bperp_date = h5py.File(s_file,'r')
        s_bperp = s_bperp_date['bperp']
        #s_date = s_bperp_date['date']
        s_dateList = timeseries(s_file).get_date_list()
        
        # judging whether dominant and affiliate data have same dimension
        m_dim = m_data.shape[0]
        m_rows, m_colms = m_data.shape[1:3]
        s_dim = s_data.shape[0]
        s_rows, s_colms = s_data.shape[1:3]
        mosaic_dataset = dict()

        date_final, m_Date, s_Date, bperp = date_match(m_dateList, s_dateList, m_dim, s_dim, m_bperp, s_bperp)
        mosaic_dim = len(m_Date)

        offset,m_row0,m_colm0,s_row0,s_colm0,over_lat0,over_lon0,overlay_rows, overlay_colms = calculate_overlay(inps,m_atr,m_data[0,:,:],s_atr,s_data[0,:,:],typeflag)
        mosaic_rows = m_rows + s_rows - overlay_rows
        mosaic_colms = m_colms + s_colms - overlay_colms
        mosaic_timeseries = np.empty(shape=(mosaic_dim, mosaic_rows, mosaic_colms))
        
        i = 0
        for m_date,s_date in zip(m_Date, s_Date):
            # calculated offset and mosaic two tracks
            dominant_data = readfile.read(m_file, datasetName=m_date)[0]
            affiliate_data = readfile.read(s_file, datasetName=s_date)[0]
            offset,m_row0,m_colm0,s_row0,s_colm0,over_lat0,over_lon0,overlay_rows, overlay_colms = calculate_overlay(inps,m_atr,dominant_data,s_atr,affiliate_data,typeflag)
            s_data_offset = rewrite_affiliate(inps,offset,s_atr,affiliate_data)
                 
            #print('prepare mosaicing:\n')
            mosaic_timeseries[i,:,:], mosaic_atr = mosaic_tracks(inps,m_atr,dominant_data,s_atr,s_data_offset,m_row0,m_colm0,s_row0,s_colm0,over_lat0,over_lon0,overlay_rows,overlay_colms)
            i = i + 1
    
        mosaic_dataset['bperp'] = np.array(bperp, dtype=np.float32)
        mosaic_dataset['date'] = np.array(date_final, dtype=np.string_)
        mosaic_dataset['timeseries'] = mosaic_timeseries

        mosaic_atr['REF_DATE'] = str(date_final[0])
        write_mosaic(inps,mosaic_dataset, mosaic_atr)
    elif typeflag == 3:
        m_file = "".join(inps.dominant)
        s_file = "".join(inps.affiliate)
        
        mosaic_dataset = dict()
        
        m_data = readfile.read(m_file)[0]
        s_data = readfile.read(s_file)[0]
        # calculated overlay region
        m_row0,m_colm0,s_row0,s_colm0,over_lat0,over_lon0,overlay_rows, overlay_colms = calculate_overlay(inps,m_atr,m_data,s_atr,s_data,typeflag)
                
        #print('prepare mosaicing for: %s\n' % dslice)
        mosaic_data, mosaic_atr = mosaic_tracks(inps,m_atr,m_data,s_atr,s_data,m_row0,m_colm0,s_row0,s_colm0,over_lat0,over_lon0,overlay_rows,overlay_colms)
        #print('finish mosaic %s' % dslice)       
        write_mosaic(inps,mosaic_dataset, mosaic_atr)
        
                  
######################################################################################
if __name__ == '__main__':
    main()
