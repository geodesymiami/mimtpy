#!/usr/bin/env python3
############################################################
# Program is part of MimtPy                                #
# Author: Zhang Yunjun                                     #
############################################################



import os
import sys
import argparse
import numpy as np
import scipy.io as sio

EXAMPLE = """example:
  convert_Kite2GBIS.py KiteScene_name 0.056 -o AlosDT73_20081012_20100302.mat
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Convert Kite Scene to GBIS.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('KiteScene', help='The name of KiteScene without suffix.')

    parser.add_argument('wavelength', type= float, help='the wavelength of data. Unit is meter')
    
    parser.add_argument('-o', '--output', dest='outfile', help='output file name.')
    
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    return inps

def read_yml(inps):
    """read yaml file"""
    Kite_config = inps.KiteScene + '.yml'
    with open(Kite_config,'r') as f:
        lines = []
        for line in f.readlines():
            if line != '\n':
                lines.append(line)
    f.close() 
   
    # check unit
    unit = lines[16].split(':')[1].strip().replace('\n', '').replace('\r', '')
    if unit == 'degree':
        llat = float(lines[12].split(':')[1])
        llon = float(lines[13].split(':')[1])
        lat_step = float(lines[14].split(':')[1])
        lon_step = float(lines[15].split(':')[1])

    return llat, llon, lat_step, lon_step

def read_data(inps):
    """
    Returns: defo: 2D np.array with in-valid/masked-out pixel in NaN
    """
    # metadata
    ll_lat, ll_lon, lat_step, lon_step = read_yml(inps) 

    # data
    kite_dict = np.load(inps.KiteScene + '.npz')

    # the disp data should be phase with unit is radian
    inps.disp = np.flipud(kite_dict['arr_0'])
    range2phase =  -4. * np.pi / inps.wavelength
    inps.disp *= range2phase

    # obtain incidence angle in decimal degrees 
    r2d = 180 / np.pi
    inps.incidence = np.flipud(kite_dict['arr_1'])
    inps.incidence *= r2d
  
    # obtain head angle in decimal degrees
    phi = np.flipud(kite_dict['arr_2'])
    phi *= r2d
    inps.head = 180 - phi

    # obtain latitude in decimal degrees
    row, colm = inps.disp.shape
    lat_max = ll_lat + row * lat_step
    lat_list = np.arange(ll_lat + lat_step, lat_max + lat_step, lat_step)
    lat_colm = np.flipud(np.transpose(np.array([lat_list])))
    inps.lat = np.tile(lat_colm, colm)

    lon_max = ll_lon + colm * lon_step
    lon_list = np.arange(ll_lon + lon_step, lon_max + lon_step, lon_step)
    lon_row = np.array([lon_list])
    inps.lon = np.tile(lon_row, (row, 1))

    # remove the nan pixel
    inps.mask = ~np.isnan(inps.disp)
      
    return

def save2mat(inps):
    """write mat file"""
    mdict = {}
    # required by GBIS
    mdict['Heading'] = inps.head[inps.mask].reshape(-1,1)
    mdict['Inc'] = inps.incidence[inps.mask].reshape(-1,1)
    mdict['Lat'] = inps.lat[inps.mask].reshape(-1,1)
    mdict['Lon'] = inps.lon[inps.mask].reshape(-1,1)
    mdict['Phase'] = inps.disp[inps.mask].reshape(-1,1)
    # optional
    #mdict['Height'] = inps.height[inps.mask].reshape(-1,1)
    #mdict['Mask'] = inps.mask
    #mdict['Metadata'] = inps.metadata
    # save to mat file
    sio.savemat(inps.outfile, mdict, long_field_names=True)
    print('save to file: {}'.format(os.path.abspath(inps.outfile)))

    return

##############################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    read_data(inps)

    save2mat(inps)

    return inps.outfile

##########################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
