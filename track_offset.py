#!/usr/bin/env python3
######################################################################################################
# Program is used for mosaic two tracks based on the  overlay region                                 #
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
    
    track_offset.py $SCRATCHDIR/BogdSenDT33/full_mintpy/mintpy/velocity.h5 $SCRATCHDIR/BogdSenDT106/full_mintpy/mintpy/velocity.h5 -b 45 46 101 102 
    track_offset.py $SCRATCHDIR/BogdSenDT33/full_mintpy/mintpy/velocity.h5 $SCRATCHDIR/BogdSenDT106/full_mintpy/mintpy/velocity.h5 -output velocity_off.h5  
    track_offset.py $SCRATCHDIR/BogdSenDT33/full_mintpy/mintpy/velocity.h5 $SCRATCHDIR/BogdSenDT106/full_mintpy/mintpy/velocity.h5  
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Mosaic two different tracks based on the overlay region',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('master', nargs=1, help='master track. \n')

    parser.add_argument('slave', nargs=1, help='slave track. \n ')

    parser.add_argument('-b', '--bbox', dest='SNWE', type=float, nargs=4, metavar=('S', 'N', 'W', 'E'),
                        help='Bounding box of area to be geocoded.\n' +
                        'Include the uppler left corner of the first pixel' +
                        '    and the lower right corner of the last pixel')
    
    parser.add_argument('--output',dest='output',nargs=1,help='output name')

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
    row = int(np.floor((lat_n - lat0) / lat_step))
    colm = int(np.floor((lon_n - lon0) / lon_step))
     
    return row, colm

def calculate_roco_overlay(lat0,lon0,lat1,lon1,lat_step,lon_step):
    rows = int(np.ceil((lat1 - lat0) / lat_step))
    colms = int(np.ceil((lon1 - lon0) / lon_step))

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
    
    m_overlay = m_data[m_row0 : m_row0 + overlay_rows,m_colm0 : m_colm0 + overlay_colms]
    s_overlay = s_data[s_row0 : s_row0 + overlay_rows,s_colm0 : s_colm0 + overlay_colms]

    offset_overlay = m_overlay - s_overlay
    offset = np.nanmean(offset_overlay)
    if 'Y_FIRST' in m_atr.keys():
        ref_lat = m_atr['REF_LAT']
        ref_lon = m_atr['REF_LON']
        ref_x = m_atr['REF_X']
        ref_y = m_atr['REF_Y']
    else:
        ref_x = m_atr['REF_X']
        ref_y = m_atr['REF_Y']
    
    print('The average offset is :')
    print(offset) 
    return offset

def rewrite_slave(inps,offset):
    s_atr = readfile.read_attribute("".join(inps.slave))
    s_data = readfile.read("".join(inps.slave))[0] 

    s_data_offset = s_data - offset
    if inps.output == None: 
        out_file = '{}{}{}'.format(os.path.splitext("".join(inps.slave))[0], '_offset', os.path.splitext("".join(inps.slave))[1])
    else:
        out_file = "".join(inps.output)

    writefile.write(s_data_offset, out_file=out_file, metadata=s_atr)

######################################################################################
def main(iargs=None):
    """simulate LOS displacement"""
    inps = cmd_line_parse(iargs)

    offset = calculate_overlay(inps)
    rewrite_slave(inps,offset)
    
######################################################################################
if __name__ == '__main__':
    main()
