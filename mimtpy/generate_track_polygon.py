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
  
  generate_track_polygon.py S1***.he5 --output BogdSenDT4.gmt --outdir ./
  
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Generate scene footprint of each track',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs=1, type=str, help='HDFEOS file for the project\n')

    parser.add_argument('--output', dest='output', type=str, nargs=1,
                        help='output file name')
    
    parser.add_argument('--outdir',dest='outdir',nargs=1,
                        help='output dir')
    
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps

def generate_polygon(inps):
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
    
    # write gmt file
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

######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)   
    
    generate_polygon(inps)
######################################################################################
if __name__ == '__main__':
    main()
