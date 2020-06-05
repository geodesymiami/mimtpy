#!/usr/bin/env python3

######################################################################################################
# Program is used for transfering *.unw and *.h5 files to geotiff file with WGS projection           #
# Author: Lv Xiaoran                                                                                 #
# Created: February  2020                                                                            #
######################################################################################################

import sys
import argparse
import os
import re
import numpy as np
import scipy.io as io
from osgeo import gdal, gdal_array, osr
from mintpy.utils import readfile, writefile

######################################################################################
EXAMPLE = """example:
    H5UNW_to_geotiff.py /data/lxr/insarlab/SCRATCHDIR/BalochistanSenAT/ERA/geo_velocity_msk_mosaic.h5 --outdir /data/lxr/insarlab/SCRATCHDIR/BalochistanSenAT/ERA/ --unit radian --output geo_velocity_msk_mosaic.tif
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Generate *.tiff file with WGS projection based *.h5 / *.unw files',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('input_HDFEOS', nargs=1, type=str, help='directory stored *.he5 files. \n')
    
    parser.add_argument('--unit', dest = 'unit', nargs='?', type=str, help='the unit of input data.Make transefer for m and radian.\n')   
 
    parser.add_argument('--outdir',dest='outdir',nargs=1, type = str, help='output directory.\n')

    parser.add_argument('--output',dest='output', nargs=1, type=str, help='output file name.\n')
    
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps  

def hdf_to_geotif(inps):
    """transfer"""
    fname = inps.input_HDFEOS[0] # for asc
    output_path = inps.outdir[0]
    output_tif = output_path + inps.output[0]
    print('output tiff file name : %s' % output_tif)
    atr = readfile.read_attribute(fname)
    data = readfile.read(fname)[0]

    if not inps.unit == 'None':
        print('using the orgianl unit of input')
        disp = data
    elif inps.unit[0] == 'm':
        # the unit of velocity is m/year
        disp = data * 100 # change unit to cm/year
    elif inps.unit[0] == 'radian':
        # the unit of data is radian
        wavelength = float(atr['WAVELENGTH'])
        disp = (((-1) * (data * wavelength)) / (4 * np.pi)) * 100 # change unit to cm

    xmin = float(atr['X_FIRST'])  # corresponding to X_FIRST
    ymax = float(atr['Y_FIRST'])  # corresponding to Y_FIRST
    nrows, ncols = np.shape(disp)
    
    xres = float(atr['X_STEP']) # corresponding to X_STEP
    yres = (-1) * float(atr['Y_STEP'])  # corresponding to X_STEP

    # xmin, ymin, xmax, ymax = [np.nanmin(lon), np.nanmin(lat), np.nanmax(lon), np.nanmax(lat)]
    # nrows, ncols = np.shape(disp)
    # xres = (xmax - xmin) / float(ncols)
    # yres = (ymax - ymin) / float(nrows)
    geotransform = [xmin, xres, 0, ymax, 0, -yres]
    raster = gdal.GetDriverByName('GTiff').Create(output_tif, ncols, nrows, 1, gdal.GDT_Float32)
    raster.SetGeoTransform(geotransform)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)  # WGS 84 - WGS84 - World Geodetic System 1984, used in GPS
    raster.SetProjection(srs.ExportToWkt())
    raster.GetRasterBand(1).WriteArray(disp)
    raster.FlushCache()

    print('finish conversion!')
######################################################################################
def main(iagrs=None):
    inps = cmd_line_parse(iagrs)
    hdf_to_geotif(inps)
######################################################################################
if __name__ == '__main__':
    main()
