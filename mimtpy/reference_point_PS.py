#!/usr/bin/env python3
#########################################################################################################
# Program is used for changing the reference point for data under radar coordinate using lat/lon        #
# Author: Lv Xiaoran                                                                                    #
# Created: June 2024                                                                                    #
#########################################################################################################

import os
import argparse
import numpy as np

import mintpy
from mintpy.utils import readfile, writefile, utils as ut
######################################################################################
EXAMPLE = """example:
  
    reference_point_PS.py velcity.h5 --geo-fle ./inputs/geoemtryRadar.h5 --lalo 37.0 118.2 --output velocity_refpoi.h5
"""

def create_parser():
    parser = argparse.ArgumentParser(description='generating horz and vert data for PS points',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs=1, type=str, help='ascending and descending data, respectively\n')

    parser.add_argument('--geo-file', dest='geoFile', type=str, nargs=1, help='geometry file')
    
    parser.add_argument('--lalo', dest='lalo', type=float, nargs=2, help='Selected reference point')
    
    parser.add_argument('--output', dest='output', nargs=1, type=str,
                        help='output name')
    
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps

def change_refpoi(lat_data, lon_data, nanMask, ref_lat, ref_lon):
    def haversin(theta):
        v = np.sin(theta / 2)
        return v * v

    def distance2points(lat1, lon1, lat2, lon2):
        radius = 6370
    
        lat1 = np.radians(lat1)
        lon1 = np.radians(lon1)
        lat2 = np.radians(lat2)
        lon2 = np.radians(lon2)
        
        dlon = lon2 - lon1
        dlat = lat2 - lat1
    
        h = haversin(dlat) + np.cos(lat1) * np.cos(lat2) * haversin(dlon)
    
        dis = 2 * radius * np.sin(np.sqrt(h))
        return dis

    # calculate the distance between points
    dis_matrix = distance2points(ref_lat, ref_lon, lat_data, lon_data)
    dis_matrix += nanMask
    row, col = divmod(np.nanargmin(dis_matrix), np.shape(dis_matrix)[1])

    return row, col

def reference_point(data, data_atr, lat_data, lon_data, ref_lat, ref_lon):
    # generate nan mask
    nanMask = np.zeros((data.shape[0], data.shape[1]))
    nanMask[np.isnan(data)] = np.nan
    
    # find the row and col of reference point
    row, col = change_refpoi(lat_data, lon_data, nanMask, ref_lat, ref_lon)

    refLat = lat_data[row, col]
    refLon = lon_data[row, col]
    print('The reference lat and lon is {} and {}'.format(refLat, refLon))

    # adjust the data
    data -= data[row, col]
    data_atr['REF_X'] = str(col)
    data_atr['REF_Y'] = str(row)

    return data, data_atr
######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)   

    data, data_atr = readfile.read(inps.file[0])
    lat_data = readfile.read(inps.geoFile[0], datasetName='latitude')[0]
    lon_data = readfile.read(inps.geoFile[0], datasetName='longitude')[0]

    ref_lat = inps.lalo[0]
    ref_lon = inps.lalo[1]

    data, data_atr = reference_point(data, data_atr, lat_data, lon_data, ref_lat, ref_lon)

    outfile = inps.output[0]
    writefile.write(datasetDict=data, out_file=outfile, metadata=data_atr)
######################################################################################
if __name__ == '__main__':
    main()
