#!/usr/bin/env python3
#############################################################
# Program is part of MimtPy                                 #
# Author: Lv Xiaoran Feb 2020                               #
#############################################################


import os
import argparse
import numpy as np
from mintpy.utils import readfile
import json
###########################################################################################
EXAMPLE = """example:
  relax_grd2point.py 063-058-relax-geo --tfile $SCRATCHDIR/BogdSenDT4/mintpy/displacement/quadtree/lola.txt --outdir ../../result_coarse/ --outfile vs_1717
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Extract displacement and incidence and azimuth angle from InSAR',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)
    parser.add_argument('file', nargs='+',
                        help='file to extract value\n')

    parser.add_argument('--tfile',type=str, dest='tfile', nargs=1,
                       help='txt file stored lalo information of points')
    
    parser.add_argument('--outdir',type=str,dest='outdir',nargs=1,
                        help='out put dir.')
    
    parser.add_argument('--outfile',type=str,dest='outfile',nargs=1,
                        help='out put dir.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    return inps

###########################################################################################
def read_txt(inps):
    """read txt file."""
    txt_file = "".join(inps.tfile)
    lon, lat = np.loadtxt(txt_file,skiprows=0,dtype=np.float,usecols =(0, 1),unpack=True) 
    return lat,lon

def get_location(file_name,lat,lon):
    """get the row/clom in file according the lat/lon.""" 
    #get lon0/lat0
    #lat0,lon0- starting latitude/longitude (first row/column)
    #lat1,lon1- ending latitude/longitude (last row/column)
    atr = readfile.read(file_name)[1]
    length,width = int(atr['LENGTH']),int(atr['WIDTH'])
   
    if atr['Y_FIRST']: 
        # geo coordinates
        lat0 = float(atr['Y_FIRST'])
        lon0 = float(atr['X_FIRST'])
        lat_step = float(atr['Y_STEP'])
        lon_step = float(atr['X_STEP'])
        #lat1 = lat0 + (length - 1) * lat_step
        #lon1 = lon0 + (width - 1) * lon_step
    else:
        parser.print_usage()
        raise Exception('input bounding box error! Wrong latitude order!')
    row = int((lat - lat0) / lat_step + 0.5)
    colm = int((lon - lon0) / lon_step + 0.5)
    #extract value
    #value = data[row,colm]
    #coord = ut.coordinate(atr)
    #y, x = coord.geo2radar(lat, lon)[0:2]
    
    return row, colm

def extract_value(file_name,row,colm):
    """extract value"""
    # read *.unw file
    data, atr = readfile.read(file_name) 
    
    value = data[row,colm]
    
    return value

def extract_points_value(inps):
    """extract north east down three component cumulative displacement"""
    
    # obtain the name of relax grd file
    east_name = inps.file[0] + '-east.h5'
    north_name = inps.file[0] + '-north.h5'
    up_name = inps.file[0] + '-up.h5'

    # read meta data
    atr = readfile.read(east_name)[1]
    
    # read lalo txt file
    lat,lon = read_txt(inps)
    row_number = lat.shape[0]
    print('%d points' %row_number)
    
    # extrac three components
    # the content of values is [lon,lat,north,east,down]
    points_value = np.empty(shape=[0,5],dtype=float)
    
    for lat_single,lon_single in zip(lat,lon):
        row, colm = get_location(east_name,lat_single,lon_single) 
        
        value_east = extract_value(east_name, row, colm)
        value_north = extract_value(north_name, row, colm)
        value_up = extract_value(up_name, row, colm)
        values_tmp = np.array([lon_single, lat_single, value_north, value_east, -1*value_up]).reshape(1,5)
        points_value = np.append(points_value,values_tmp, axis=0)

 
    #write in json file
    text_file = inps.outdir[0] + inps.outfile[0] 
    write_files(points_value,text_file)

def write_files(data,file_name):
    """write lalo displacement data into *.json"""

    displacement = {"displacement": data.tolist()}
    open(file_name + '.json', "w").write(json.dumps(displacement))

###########################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    extract_points_value(inps)

#########################################################################################
if __name__ == '__main__':
    main()
