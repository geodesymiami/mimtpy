#!/usr/bin/env python3
#################################################################
# Program is used for preparing data for kite software          #
# Author: Lv Xiaoran                                            #
# Created: Jan 2021                                             #
#################################################################

import os
import argparse
import numpy as np
from kite import Scene

import mintpy
import mintpy.workflow
from mintpy.utils import readfile, ptime, utils as ut
from mintpy import view
######################################################################################
EXAMPLE = """example:
  
  save_kite.py S1****.he5 --date 20141016_20171221 -o SenAT115_201410_201712_kite
  
  save_kite.py S1****.he5 --date 20141016_20171221 --bbox 26 27 64 65 -o SenAT115_201410_201712_kite
  
  save_kite.py S1****.he5 --date 20141016_20171221 --velocity -o SenAT115_201410_201712_kite
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Prepare data for Kite software',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('input_HDFEOS', nargs=1, type=str, help='HDFEOS file. \n')

    parser.add_argument('--date', dest='date',nargs=1, help='formate: date1_date2. \n')

    parser.add_argument('--velocity', action='store_true', default=False, help='whether calculate velocity. \n')

    parser.add_argument('--bbox', dest='SNWE', type=float, nargs=4, metavar=('S', 'N','W','E'),
                        help='Bounding box of area to be geocoded.\n'+
                        'Include the uppler left corner of the first pixel' +
                        'and the lower right corner of the last pixel')
    
    parser.add_argument('-o','--outfile',dest='outfile',nargs=1,
                        help='outfile name.')
    
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps

def extract_data_based_bbox(inps):
    """return the row_no,sample_no and rows and samples"""
    # metadata
    atr = readfile.read_attribute("".join(inps.input_HDFEOS))
    ul_lat = float(atr['Y_FIRST'])
    ul_lon = float(atr['X_FIRST'])
    lat_step = float(atr["Y_STEP"])
    lon_step = float(atr["X_STEP"])
    length = int(atr['LENGTH'])
    width = int(atr['WIDTH'])
    lr_lat = ul_lat + lat_step * length
    lr_lon = ul_lon + lon_step * width
    # bbox
    user_lat0 = float(inps.SNWE[1])
    user_lon0 = float(inps.SNWE[2])
    user_lat1 = float(inps.SNWE[0])
    user_lon1 = float(inps.SNWE[3])

    if user_lat0 > ul_lat or user_lat1 < lr_lat or user_lon0 < ul_lon or user_lon1 > lr_lon:
        print('data range is {} {} {} {}'.format(lr_lat, ul_lat, ul_lon, lr_lon))
        raise Exception('input bounding box is beyound the data range!\n')

    if user_lat0 < user_lat1:
        parser.print_usage()
        raise Exception('input bounding box error! Wrong latitude order!')
    elif user_lon0 > user_lon1:
        parser.print_usage()
        raise Exception('input bounding box error! Wrong longitude order!')

    row = int((user_lat0 - ul_lat) / lat_step + 0.5)
    sample = int((user_lon0 - ul_lon) / lon_step + 0.5)
    rows = int((user_lat1 - user_lat0) / lat_step + 0.5) + 1
    samples = int((user_lon1 - user_lon0) / lon_step + 0.5) + 1
    
    return user_lat1, user_lon0, row, sample, rows,samples

def read_HDFEOS(inps):
    """read displacement from HDFEOS"""
    print('read displacement, incidence and azimuth information')
    # read metadata
    HDFEOS_file = inps.input_HDFEOS[0]
    atr = readfile.read_attribute(HDFEOS_file)

    if inps.date == None:
        date1 = atr['START_DATE']
        date2 = atr['END_DATE']
    else:
        # date1 and date2
        if '_' in "".join(inps.date):
            date1, date2 = ptime.yyyymmdd("".join(inps.date).split('_'))
        else:
            date1 = atr['START_DATE']
            date2 = ptime.yyyymmdd("".join(inps.date))
    # read angle infomation
    azimuth = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/geometry/azimuthAngle')[0]
    incidence = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/geometry/incidenceAngle')[0]
    
    if inps.velocity:
        vel_file = 'velocity.h5'

        iargs = [HDFEOS_file, '--start-date', date1, '--end-date', date2, '-o', vel_file, '--update']
        print('\ntimeseries2velocity.py', ' '.join(iargs))
        mintpy.timeseries2velocity.main(iargs)

        data = readfile.read(vel_file, datasetName='velocity')[0]
        os.remove(vel_file) 
    else: 
        # read / prepare data
        slice_list = readfile.get_slice_list(HDFEOS_file)
        # read/prepare data
        dname = 'displacement'
        slice_name1 = view.search_dataset_input(slice_list, '{}-{}'.format(dname, date1))[0][0]
        slice_name2 = view.search_dataset_input(slice_list, '{}-{}'.format(dname, date2))[0][1]
        data = readfile.read("".join(inps.input_HDFEOS), datasetName=slice_name2)[0]
        data -= readfile.read("".join(inps.input_HDFEOS), datasetName=slice_name1)[0]
    
    print("mask file")
    maskfile = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/quality/mask')[0]
    data[maskfile == 0] = np.nan
    azimuth[maskfile == 0] = np.nan
    incidence[maskfile == 0] = np.nan
    
    return data, atr, incidence, azimuth

def save_kite(inps):
    """save kite"""
    # read mintpy data
    disp, disp_atr, incidence, azimuth = read_HDFEOS(inps)

    # subset data based on bbox
    if inps.SNWE:
        print('Subset data based on bbox')
        lat_user, lon_user, row, sample, rows,samples = extract_data_based_bbox(inps)
        disp = disp[row: row+rows, sample: sample+samples]
        incidence = incidence[row: row+rows, sample: sample+samples]
        azimuth = azimuth[row: row+rows, sample: sample+samples]
    
    # convert to head angle
    print('convert azimuth angle to head angle')
    head = ut.azimuth2heading_angle(azimuth)
    phi = -head + 180

    # convert degree to radian
    incidence *= np.pi/180
    phi *= np.pi/180

    sc = Scene()
    # flip up-down of displacement and angle matrix
    sc.displacement = np.flipud(disp)
    
    sc.theta = np.flipud(incidence)
    sc.phi = np.flipud(phi)
    # calculate the scene's frame lower left corner, in geographical coordinate
    lon_ul = float(disp_atr['X_FIRST'])
    lat_ul = float(disp_atr['Y_FIRST'])

    lat_step = float(disp_atr['Y_STEP'])
    length = int(disp_atr['LENGTH'])

    lat_ll = lat_ul + lat_step * length
    lon_ll = lon_ul
 
    if inps.SNWE:
        lat_ll = lat_user 
        lon_ll = lon_user

    sc.frame.llLat = lat_ll
    sc.frame.llLon = lon_ll

    # the pixel spacing can be either 'meter' or 'degree'
    sc.frame.spacing = disp_atr['X_UNIT'][:-1]
    sc.frame.dN = float(disp_atr['Y_STEP']) * (-1)
    sc.frame.dE = float(disp_atr['X_STEP'])

    # Saving the scene
    print('write Kite scene')
    sc.save(inps.outfile[0])
######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)   
    
    save_kite(inps)
######################################################################################
if __name__ == '__main__':
    main()
