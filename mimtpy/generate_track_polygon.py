#!/usr/bin/env python3
################################################################################
# Program is used for generating scene footprint of each track with gmt format #
# Author: Lv Xiaoran                                                           #
# Created: July 2020                                                           #
################################################################################

import os
import argparse
import re

from mintpy.utils import readfile
######################################################################################
EXAMPLE = """example:
  
  generate_track_polygon.py S1***.he5 --option scfoot --output BogdSenDT4.gmt --outdir ./
  generate_track_polygon.py S1***.he5 --option cbox --output BogdSenDT4.gmt --outdir ./
  
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Generate scene footprint of each track',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs=1, type=str, help='HDFEOS file for the project\n')

    parser.add_argument('--option', nargs=1, type=str, help='two options: generate screen_footprint or coverage box'
                                                            'sfoot: screen_footprint; cbox: coverage box')

    parser.add_argument('--output', dest='output', type=str, nargs=1,
                        help='output file name')
    
    parser.add_argument('--outdir',dest='outdir',nargs=1,
                        help='output dir')
    
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps

def generate_polygon_sfoot(inps):
    """generate polygon for track using gmt formate"""
    file_name = inps.file[0]
  
    atr = readfile.read(file_name)[1]

    # extract lon/lat of four points
    polygon = atr['scene_footprint']
    lonlats = re.findall(r'([\d+\.]+)',polygon)
    
    # lon/lat of four points: upperright; lowerright; lowerleft; upperleft 
    lon_ur = float(lonlats[0])
    lat_ur = float(lonlats[1])

    lon_lr = float(lonlats[2])
    lat_lr = float(lonlats[3])

    lon_ll = float(lonlats[4])
    lat_ll = float(lonlats[5])

    lon_ul = float(lonlats[6])
    lat_ul = float(lonlats[7])
    
    return lon_ur, lat_ur, lon_lr, lat_lr, lon_ll, lat_ll, lon_ul, lat_ul
    
def write_gmt(lon_ur, lat_ur, lon_lr, lat_lr, lon_ll, lat_ll, lon_ul, lat_ul, inps): 
    '''write gmt file'''
    gmt_file = inps.outdir[0] + '/' + inps.output[0]
    
    f = open(gmt_file, mode='w')
    f.write('# @VGMT1.0 @GLINESTRING \n')
    f.writelines(['# @R',str(min(lon_ll, lon_ul)),'/',str(max(lon_ur, lon_lr)),'/',str(min(lat_ll, lat_lr)),'/', str(max(lat_ur, lat_ul)),'\n'])
    f.write('# @Je4326 \n')
    f.write('# @Jp"+proj=longlat +datum=WGS84 +no_defs" \n')
    f.write('# @Jw"GEOGCS[\\"WGS 84\\",DATUM[\\"WGS_1984\\",SPHEROID[\\"WGS 84\\",6378137,298.257223563,AUTHORITY[\\"EPSG\\",\\"7030\\"]],AUTHORITY[\\"EPSG\\",\\"6326\\"]],PRIMEM[\\"Greenwich\\",0,AUTHORITY[\\"EPSG\\",\\"8901\\"]],UNIT[\\"degree\\",0.0174532925199433,AUTHORITY[\\"EPSG\\",\\"9122\\"]],AXIS[\\"Latitude\\",NORTH],AXIS[\\"Longitude\\",EAST],AUTHORITY[\\"EPSG\\",\\"4326\\"]]" \n')
    f.write('# @NId \n')
    f.write('# @Tinteger \n')
    f.write('# FEATURE_DATA \n')
    f.write('>')
    f.write('# @D0 \n')
    f.writelines([str(lon_ur), ' ', str(lat_ur), '\n'])
    f.writelines([str(lon_lr), ' ' , str(lat_lr), '\n'])
    f.writelines([str(lon_ll), ' ' , str(lat_ll), '\n'])
    f.writelines([str(lon_ul), ' ' , str(lat_ul), '\n'])
    f.writelines([str(lon_ur), ' ', str(lat_ur), '\n'])
    f.close() 

def generate_polygon_cbox(inps):
    """generate covarage box
    lat0,lon0- starting latitude/longitude (first row/column)
    lat1,lon1- ending latitude/longitude (last row/column)
    """
    file_name = inps.file[0] 
    atr = readfile.read(file_name)[1]
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
    
    # lon/lat of four points: upperright; lowerright; lowerleft; upperleft 
    lon_ur = lon1
    lat_ur = lat0
 
    lon_lr = lon1
    lat_lr = lat1
 
    lon_ll = lon0
    lat_ll = lat1

    lon_ul = lon0
    lat_ul = lat0

    return lon_ur, lat_ur, lon_lr, lat_lr, lon_ll, lat_ll, lon_ul, lat_ul
######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)   
    
    option = inps.option
    if option == 'sfoot':
        lon_ur, lat_ur, lon_lr, lat_lr, lon_ll, lat_ll, lon_ul, lat_ul = generate_polygon_sfoot(inps)
    else:
        lon_ur, lat_ur, lon_lr, lat_lr, lon_ll, lat_ll, lon_ul, lat_ul = generate_polygon_cbox(inps)

    write_gmt(lon_ur, lat_ur, lon_lr, lat_lr, lon_ll, lat_ll, lon_ul, lat_ul, inps)
######################################################################################
if __name__ == '__main__':
    main()
