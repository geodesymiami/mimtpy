#!/usr/bin/env python3
######################################################################################################
# Program is used for concatenating two tracks based on the overlay region                           #
# Author: Lv Xiaoran                                                                                 #
# Created: January  2020                                                                             #
######################################################################################################

import os
import sys
import argparse
import numpy as np
import json
from mintpy.utils import readfile, writefile
######################################################################################
EXAMPLE = """example:
    
    track_offset.py $SCRATCHDIR/BogdSenDT33/full_mintpy/mintpy/geo_velocity.h5 $SCRATCHDIR/BogdSenDT106/full_mintpy/mintpy/geo_velocity.h5 -b 45 46 101 102 
    track_offset.py $SCRATCHDIR/BogdSenDT33/full_mintpy/mintpy/geo_velocity.h5 $SCRATCHDIR/BogdSenDT106/full_mintpy/mintpy/geo_velocity.h5 -output velocity_off.h5  
    track_offset.py $SCRATCHDIR/BogdSenDT33/full_mintpy/mintpy/geo_velocity.h5 $SCRATCHDIR/BogdSenDT106/full_mintpy/mintpy/geo_velocity.h5 --outdir /data/lxr/insarlab/SCRATCHDIR/BalochistanSenAT/  
    track_offset.py $SCRATCHDIR/BogdSenDT33/full_mintpy/mintpy/geo_velocity.h5 $SCRATCHDIR/BogdSenDT106/full_mintpy/mintpy/geo_velocity.h5  --rewrite_slave --mosaic --outdir /data/lxr/insarlab/SCRATCHDIR/BalochistanSenDT/cumulative/
    track_offset.py $SCRATCHDIR/BogdSenDT33/full_mintpy/mintpy/geo_velocity.h5 $SCRATCHDIR/BogdSenDT106/full_mintpy/mintpy/geo_velocity.h5 --mosaic --output mosaic --outdir /data/lxr/insarlab/SCRATCHDIR/BalochistanSenDT/cumulative/
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
    parser.add_argument('--mosaic', action='store_true', default=False, help='whether mosaic two track data.') 
    parser.add_argument('--output',dest='output',nargs=1,help='output name')
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

def calculate_overlay(inps):
    parser = create_parser()
    
    m_atr = readfile.read_attribute("".join(inps.master))
    m_data = readfile.read("".join(inps.master))[0]
    s_atr = readfile.read_attribute("".join(inps.slave))
    s_data = readfile.read("".join(inps.slave))[0]

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
   
    return offset, m_row0, m_colm0, s_row0, s_colm0, over_lat0, over_lon0, overlay_rows, overlay_colms

def rewrite_slave(inps,offset):
    s_atr = readfile.read_attribute("".join(inps.slave))
    s_data = readfile.read("".join(inps.slave))[0] 

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
 
def mosaic_tracks(inps,s_data_offset,m_row0,m_colm0,s_row0,s_colm0,over_lat0,over_lon0,overlay_rows,overlay_colms):
    """mosaic two tracks"""
    m_atr = readfile.read_attribute("".join(inps.master))
    m_data = readfile.read("".join(inps.master))[0]
    m_rows,m_colms = m_data.shape
    
    s_atr = readfile.read_attribute("".join(inps.slave))
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
    mosaic_atr = m_atr
    mosaic_atr['LENGTH'] = mosaic_rows
    mosaic_atr['WIDTH'] = mosaic_colms
    mosaic_atr['X_FIRST'] = mosaic_lon0
    mosaic_atr['Y_FIRST'] = mosaic_lat0
    
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
        
######################################################################################
def main(iargs=None):
    """simulate LOS displacement"""
    inps = cmd_line_parse(iargs)

    offset,m_row0,m_colm0,s_row0,s_colm0,over_lat0,over_lon0,overlay_rows, overlay_colms = calculate_overlay(inps)
    s_data_offset = rewrite_slave(inps,offset)
    if inps.mosaic:
        print('prepare mosaicing:\n')
        mosaic_tracks(inps,s_data_offset,m_row0,m_colm0,s_row0,s_colm0,over_lat0,over_lon0,overlay_rows,overlay_colms)
######################################################################################
if __name__ == '__main__':
    main()
