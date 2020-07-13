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
import matplotlib.pyplot as plt
from mintpy.utils import readfile, writefile

######################################################################################
EXAMPLE = """example:
    track_offset.py $SCRATCHDIR/BogdSenDT33/mintpy/geo_velocity.h5 $SCRATCHDIR/BogdSenDT106/mintpy/geo_velocity.h5  --rewrite_slave --outdir /data/lxr/insarlab/SCRATCHDIR/BalochistanSenDT/cumulative/

    track_offset.py $SCRATCHDIR/BogdSenDT33/mintpy/geo_velocity.h5 $SCRATCHDIR/BogdSenDT106/mintpy/geo_velocity.h5 --output mosaic --outdir /data/lxr/insarlab/SCRATCHDIR/BalochistanSenDT/cumulative/

    track_offset.py $SCRATCHDIR/BogdSenDT33/mintpy/geo_velocity.h5 $SCRATCHDIR/BogdSenDT106/mintpy/geo_velocity.h5 --output Bogd_mosaic --plotpair --azi_angle 11 --outdir ./

"""

def create_parser():
    parser = argparse.ArgumentParser(description='Concatenating two different tracks based on the overlay region',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('master', nargs=1, help='master track. \n')

    parser.add_argument('slave', nargs=1, help='slave track. \n ')

    parser.add_argument('-b', '--bbox', dest='SNWE', type=float, nargs=4, metavar=('S', 'N', 'W', 'E'),
                        help='Bounding box of area to be geocoded.\n' +
                        'Include the uppler left corner of the first pixel' +
                        '    and the lower right corner of the last pixel')
   
    parser.add_argument('--rewrite_slave', action='store_true', default=False, help='whether rewrite slave *.h5 file after add the offset.\n')
    
    plotpair_opt = parser.add_argument_group(title='whether plot profiles for overlay regions of both tracks')
    plotpair_opt.add_argument('--plotpair', action='store_true', default=False, help='whether plot profiles for overlay regions of both tracks. \n')
    plotpair_opt.add_argument('--azi_angle', dest='azimuth', type=float, help='profile direction relative to North in clockwise direction.'
                                                                              ' Range=[0,pi); 0 degree: N-S; 90 degree: E-W.\n')

    mosaic = parser.add_argument_group(title='mosaic options')
    #mosaic.add_argument('--mosaic', action='store_true', default=False, help='whether mosaic two track data.') 
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
    # lat lon range of master and slave
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

    if typeflag == 0:
        # calculate overlay region in master and slave
        m_overlay = m_data[m_row0 : m_row0 + overlay_rows,m_colm0 : m_colm0 + overlay_colms]
        s_overlay = s_data[s_row0 : s_row0 + overlay_rows,s_colm0 : s_colm0 + overlay_colms]

        # calculate offset between master and slave
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
        
        print('The average offset is : %f \n' % offset)
       
        #if inps.plotpair:
        #    out_dir = inps.outdir[0]
        #    angle = inps.azimuth
        #    overlay_profile(angle, m_overlay, s_overlay, m_name, s_name, out_dir)

        return offset, m_row0, m_colm0, s_row0, s_colm0, over_lat0, over_lon0, overlay_rows, overlay_colms
    else:
        return m_row0, m_colm0, s_row0, s_colm0, over_lat0, over_lon0, overlay_rows, overlay_colms

def profile_latlon(row_overlay, colm_overlay, row0, colm0, atr):
    """transfer pixel coordinate to lat/lon"""
    row_wholeMatrix = row_overlay + row0
    colm_wholeMatrix = colm_overlay + colm0
    
    x_step = float(atr["X_STEP"])
    y_step = float(atr["Y_STEP"])
    lat0 = float(atr["Y_FIRST"])
    lon0 = float(atr["X_FIRST"]) 

    lon = lon0 + x_step * colm_wholeMatrix
    lat = lat0 + y_step * row_wholeMatrix

    return lat, lon

def profile_write(lat_start, lon_start, lat_end, lon_end, name, outdir):
    """write lat/lon of two endpoints into gmt file"""
    gmt_file = outdir + 'profile_latlon_' + name + '.gmt'
    
    f = open(gmt_file, mode='w')
    f.write('# @VGMT1.0 @GLINESTRING \n')
    f.writelines(['# @R',str(min(lon_start,lon_end)),'/',str(max(lon_start,lon_end)),'/',str(min(lat_start,lat_end)),'/', str(max(lat_start,lat_end)),'\n'])
    f.write('# @Je4326 \n')
    f.write('# @Jp"+proj=longlat +datum=WGS84 +no_defs" \n')
    f.write('# @Jw"GEOGCS[\\"WGS 84\\",DATUM[\\"WGS_1984\\",SPHEROID[\\"WGS 84\\",6378137,298.257223563,AUTHORITY[\\"EPSG\\",\\"7030\\"]],AUTHORITY[\\"EPSG\\",\\"6326\\"]],PRIMEM[\\"Greenwich\\",0,AUTHORITY[\\"EPSG\\",\\"8901\\"]],UNIT[\\"degree\\",0.0174532925199433,AUTHORITY[\\"EPSG\\",\\"9122\\"]],AXIS[\\"Latitude\\",NORTH],AXIS[\\"Longitude\\",EAST],AUTHORITY[\\"EPSG\\",\\"4326\\"]]" \n')
    f.write('# @NId \n')
    f.write('# @Tinteger \n')
    f.write('# FEATURE_DATA \n')
    f.write('>')
    f.write('# @D0 \n')
    f.writelines([str(lon_start), ' ', str(lat_start), '\n'])
    f.writelines([str(lon_end), ' ' , str(lat_end), '\n'])
    f.close()
 
    return

def overlay_profile(angle, m_overlay, s_overlay, m_row0, m_colm0, s_row0, s_colm0, m_atr, s_atr, m_name, s_name, outdir):
    """plot profiles for overlay region of both tracks"""
    # get the size of overlay region
    rows, colms = np.shape(m_overlay)
    # get the origin position of the overlay region.
    colm_x = colms / 2
    row_y = rows / 2

    # calculat the intersect pixels between overlay region and profile
    if angle >= 45 and angle <= 135:
        # use colm to calculate row
        colm_no = np.arange(colms)
        if angle != 90:
            tan_value = -1 * math.tan(angle * np.pi / 180)
            row_no = np.ceil(((colm_no - colm_x) / tan_value) + row_y)
        elif anlge == 90:
            row_no = np.ceil(np.ones(colms) * row_y)
    else:
        # use row to calculate colm
        row_no = np.arange(rows)
        if angle != 0:
            tan_value = -1 * math.tan(angle * np.pi / 180)
            colm_no = np.ceil((row_no - row_y) * tan_value + colm_x)
        elif anlge == 0:
            colm_no = np.ceil(np.ones(rows) * colm_x)
   
    row_no = row_no.astype(dtype=np.int)
    colm_no = colm_no.astype(dtype=np.int) 
    m_profile = m_overlay[row_no, colm_no]    
    s_profile = s_overlay[row_no, colm_no] 

    # change zero value to np.nan
    m_profile[(m_profile == 0)] = np.nan
    s_profile[(s_profile == 0)] = np.nan

    # calaculate lat/lon for profiles of two tracks
    row_start = row_no[0]
    row_end = row_no[-1]
    colm_start = colm_no[0]
    colm_end = colm_no[-1]
    lat_start_m, lon_start_m = profile_latlon(row_start, colm_start, m_row0, m_colm0, m_atr)
    lat_end_m, lon_end_m = profile_latlon(row_end, colm_end, m_row0, m_colm0, m_atr)
    lat_start_s, lon_start_s = profile_latlon(row_start, colm_start, s_row0, s_colm0, s_atr)
    lat_end_s, lon_end_s = profile_latlon(row_end, colm_end, s_row0, s_colm0, s_atr) 
    
    # save lat/lon files in gmt format
    profile_write(lat_start_m, lon_start_m, lat_end_m, lon_end_m, m_name, outdir)
    profile_write(lat_start_s, lon_start_s, lat_end_s, lon_end_s, s_name, outdir)

    # plot two profiles
    figure_size = [10,8]
    fig,axes = plt.subplots(1,1,figsize = figure_size)
    ax1 = axes
    print('*****************************ploting profile************************')       
    x_axis = np.arange(1,len(m_profile)+1)
    ax1.plot(x_axis, m_profile, color='black', linestyle='-', label=m_name)
    ax1.plot(x_axis, s_profile, color='blue', linestyle='-', label=s_name)

    ax1.tick_params(which='both', direction='in', labelsize=18, bottom=True, top=True, left=True, right=True)
    font1 = {'family' : 'serif',
             'weight': 'normal',
             'size' : 18.}
    ax1.set_xlabel('Distance [km]',font1)
    ax1.set_ylabel('LOS Displacement [cm]',font1)
    labels = ax1.get_xticklabels() + ax1.get_yticklabels()
    [label.set_fontname('serif') for label in labels]
   
    ax1.legend(loc='upper left', prop=font1)
     
    #save figure
    fig_name = 'Profiles.png'
    fig_output = outdir + fig_name
    fig.savefig(fig_output, dpi=300, bbox_inches='tight')

def rewrite_slave(inps,offset,s_atr,s_data):

    s_data_offset = s_data + offset
    #if inps.output == None: 
    out_file = '{}{}{}'.format(os.path.split("".join(inps.slave))[1].split('.')[0], '_offset.', os.path.split("".join(inps.slave))[1].split('.')[1])
    #else:
    #   out_file = '{}{}{}'.format("".join(inps.output),'.',os.path.split("".join(inps.slave))[1].split('.')[1])

    if inps.outdir != None:
        out_dir = inps.outdir[0]

    out_dir_file = out_dir + out_file
    
    if inps.rewrite_slave:
        print('Writing slave data %s' % "".join(inps.slave))
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
 
def mosaic_tracks(inps,m_atr,m_data,s_atr,s_data_offset,m_row0,m_colm0,s_row0,s_colm0,over_lat0,over_lon0,overlay_rows,overlay_colms):
    """mosaic two tracks"""
    mm_atr = copy.deepcopy(m_atr) 
    m_rows,m_colms = m_data.shape
    
    s_data = s_data_offset
    s_rows,s_colms = s_data.shape
    
    # get master and slave overlay region data
    m_overlay = m_data[m_row0 : m_row0 + overlay_rows,m_colm0 : m_colm0 + overlay_colms]
    s_overlay = s_data[s_row0 : s_row0 + overlay_rows,s_colm0 : s_colm0 + overlay_colms]

    # calculate values in overlay region
    mosaic_freq = np.ones((overlay_rows,overlay_colms),dtype = np.float32)
    m_over_pos = np.isnan(m_overlay)
    s_over_pos = np.isnan(s_overlay)
    # change np.nan to zero
    m_overlay[m_over_pos] = 0
    s_overlay[s_over_pos] = 0

    mosaic_overlay_tmp = (m_overlay + s_overlay)
    # record the common np.nan in m and s(np.nan in mosaic_freq); one np.nan in m or s(1 in mosaic_freq);none np.nan in m and s(2 in mosaic_freq)
    both_nan = m_over_pos & s_over_pos
    none_nan = m_over_pos | s_over_pos
    mosaic_freq[both_nan] = np.nan
    mosaic_freq[~none_nan] = 2
    
    # calculated final mosaic_overlay
    mosaic_overlay = mosaic_overlay_tmp / mosaic_freq

    # generate mosaic dataset
    mosaic_rows = m_rows + s_rows - overlay_rows
    mosaic_colms = m_colms + s_colms - overlay_colms

    # store mosaic
    mosaic_data = np.zeros((mosaic_rows,mosaic_colms),dtype=np.float32) * np.nan 
    #mosaic_freq = np.zeros((mosaic_rows,mosaic_colms),dtype=np.float32)

    # lat lon range of master and slave
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
    
    if inps.plotpair:
        # master and slave data name
        m_name_tmp = os.path.split("".join(inps.master))[0]
        s_name_tmp = os.path.split("".join(inps.slave))[0]
        m_name= re.search('Sen([^/]+)/', m_name_tmp)[1]
        s_name = re.search('Sen([^/]+)/', s_name_tmp)[1]
        out_dir = inps.outdir[0]
        angle = inps.azimuth
        overlay_profile(angle, m_overlay, s_overlay, m_row0, m_colm0, s_row0, s_colm0, m_atr, s_atr, m_name, s_name, out_dir)

    return mosaic_data, mosaic_atr

def write_mosaic(inps, mosaic_data, mosaic_atr):
    """write mosaic data"""
    if inps.output == None: 
        out_file = '{}{}{}'.format(os.path.split("".join(inps.slave))[1].split('.')[0], '_mosaic.', os.path.split("".join(inps.slave))[1].split('.')[1])
        #out_file = '{}{}{}'.format(os.path.split("".join(inps.slave))[1].split('.')[0], '_mosaic.', 'unw')
    else:
        out_file = '{}{}{}'.format("".join(inps.output),'.',os.path.split("".join(inps.slave))[1].split('.')[1])

    if inps.outdir != None:
        out_dir = inps.outdir[0]

    out_dir_file = out_dir + out_file

    print('writing mosaic files:\n')
    writefile.write(mosaic_data, out_file=out_dir_file, metadata=mosaic_atr)

def judge_data_datasets(m_atr):
    """judage master/slave is data(velocity or .unw) or dataset (geometry)"""
    file_type = m_atr['FILE_TYPE']
    if file_type == 'geometry':
        typeflag = 1
    else:
        typeflag = 0

    return typeflag
        
######################################################################################
def main(iargs=None):
    """simulate LOS displacement"""
    inps = cmd_line_parse(iargs)

    # master and slave data attribute
    m_atr = readfile.read_attribute("".join(inps.master))
    s_atr = readfile.read_attribute("".join(inps.slave))
    
    # judge file type
    typeflag = judge_data_datasets(m_atr)
    
    if typeflag == 0:
        # master data
        m_data = readfile.read("".join(inps.master))[0]
        # slave data
        s_data = readfile.read("".join(inps.slave))[0]
         
        # calculated offset and mosaic two tracks
        offset,m_row0,m_colm0,s_row0,s_colm0,over_lat0,over_lon0,overlay_rows, overlay_colms = calculate_overlay(inps,m_atr,m_data,s_atr,s_data,typeflag)
        s_data_offset = rewrite_slave(inps,offset,s_atr,s_data)
            
        print('prepare mosaicing:\n')
        mosaic_data, mosaic_atr = mosaic_tracks(inps,m_atr,m_data,s_atr,s_data_offset,m_row0,m_colm0,s_row0,s_colm0,over_lat0,over_lon0,overlay_rows,overlay_colms)
        write_mosaic(inps,mosaic_data, mosaic_atr)
            
    else:
        dataslice = ['azimuthAngle', 'height', 'incidenceAngle', 'latitude', 'longitude', 'shadowMask', 'slantRangeDistance']
        m_file = "".join(inps.master)
        s_file = "".join(inps.slave)
        
        mosaic_dataset = dict()
        
        for dslice in dataslice:
            m_data = readfile.read(m_file, datasetName=dslice)[0]
            s_data = readfile.read(s_file, datasetName=dslice)[0]
            # calculated overlay region
            m_row0,m_colm0,s_row0,s_colm0,over_lat0,over_lon0,overlay_rows, overlay_colms = calculate_overlay(inps,m_atr,m_data,s_atr,s_data,typeflag)
                
            print('prepare mosaicing for: %s\n' % dslice)
            mosaic_data, mosaic_atr = mosaic_tracks(inps,m_atr,m_data,s_atr,s_data,m_row0,m_colm0,s_row0,s_colm0,over_lat0,over_lon0,overlay_rows,overlay_colms)
            mosaic_dataset[dslice] = mosaic_data 
            print('finish mosaic %s' % dslice)       
        write_mosaic(inps,mosaic_dataset, mosaic_atr)

          
######################################################################################
if __name__ == '__main__':
    main()
