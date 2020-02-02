#!/usr/bin/env python3
######################################################################################################
# Program is used for generating simulated LOS displacement based on MintPy and PSGRN/PSCMP results  #
# Author: Lv Xiaoran                                                                                 #
# Created: January  2020                                                                             #
######################################################################################################

import os
import sys
import argparse
import numpy as np
import json
from osgeo import gdal, osr
######################################################################################
EXAMPLE = """example:
    simulate_LOS_displacement.py $MODELOUT/pscmp/Bogd_test2/model4/poseis-10y_20y/subtract/  --outdir $MODELOUT/pscmp/Bogd_test2/model4/poseis-10y_20y/subtract/ 
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Simulate PSGRN/PSCMP LOS deformation',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('inputdir', nargs=1, help='directory stored *.json files. \n')
    
    parser.add_argument('--outdir',dest='outdir',nargs='?',default=os.getenv('MODELOUT'), help='output directory')

    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps    

def geocode(disp_data,disp_name,disp_sim,outdir):
    """geocode step"""
    # calculate geo information
    #xmin, ymin, xmax, ymax = [np.nanmin(lon), np.nanmin(lat), np.nanmax(lon), np.nanmax(lat)]
    #nrows, ncols = np.shape(disp_data)
    #xres = (xmax - xmin) / float(ncols)
    #yres = (ymax - ymin) / float(nrows)
    ymin = disp_sim["ul_lat"] + disp_sim["lat_step"] * disp_sim["rows"]
    ymax = disp_sim["ul_lat"]
    xmin = disp_sim["ul_lon"]
    xmax = disp_sim["ul_lon"] + disp_sim["lon_step"] * disp_sim["samples"]
    nrows = disp_sim["rows"]
    ncols = disp_sim["samples"]
    yres = abs(disp_sim["lat_step"])
    xres = disp_sim["lon_step"]
    # output
    output_tif = outdir +'/' + disp_name + '.tiff'
    geotransform = [xmin,xres,0,ymax,0,-yres]
    raster = gdal.GetDriverByName('GTiff').Create(output_tif,ncols,nrows,1,gdal.GDT_Float32)
    raster.SetGeoTransform(geotransform)
    srs=osr.SpatialReference()
    #srs.ImportFromEPSG(32638) #wgs84 UTM 38N
    srs.ImportFromEPSG(4326) #WGS 84 - WGS84 - World Geodetic System 1984, used in GPS
    raster.SetProjection(srs.ExportToWkt())
    raster.GetRasterBand(1).WriteArray(disp_data)
    raster.FlushCache()

def write_json(disp_data,disp_name,disp_sim_d,outdir):
    """write json file"""
    ul_lon = disp_sim_d["ul_lon"]
    ul_lat = disp_sim_d["ul_lat"]
    lat_step = disp_sim_d["lat_step"]
    lon_step = disp_sim_d["lon_step"]
    samples = disp_sim_d["samples"]
    rows = disp_sim_d["rows"]
    displacement = {"ul_lon": ul_lon, "ul_lat": ul_lat, "lat_step": lat_step, "lon_step": lon_step, "samples": samples, "rows": rows, "data": disp_data.tolist()}
    outname = disp_name + '.json'
    open(outdir + '/' + outname, "w").write(json.dumps(displacement))
    
######################################################################################
def main(iargs=None):
    """simulate LOS displacement"""
    inps = cmd_line_parse(iargs)
    
    inputdir = "".join(inps.inputdir)
    path_list = os.listdir(inputdir)
    disp_sim_name = []
    for file in path_list:
        if os.path.splitext(file)[0].split("_")[-1] == 'subtract' and os.path.splitext(file)[1] == '.json':
            disp_sim_name.append(file)
    disp_sim_name.sort()
    # read 2d simulated displacement data
    disp_sim_e = json.loads(open(inputdir + disp_sim_name[0]).read())
    disp_sim_n = json.loads(open(inputdir + disp_sim_name[1]).read())
    disp_sim_d = json.loads(open(inputdir + disp_sim_name[2]).read())
    # attribuate
    ul_lon_sim = disp_sim_d["ul_lon"]
    ul_lat_sim = disp_sim_d["ul_lat"]
    lat_step = disp_sim_d["lat_step"]
    lon_step = disp_sim_d["lon_step"]
    samples = disp_sim_d["samples"]
    rows = disp_sim_d["rows"]
    
    # displacement data
    disp_data_d = np.array(disp_sim_d["data"])
    disp_data_e = np.array(disp_sim_e["data"])
    disp_data_n = np.array(disp_sim_n["data"])
    
    #read incidence and azimuth file    
    angle_name = []
    for file in path_list:
        if os.path.splitext(file)[0].split("_")[-1] == 'resam':
            angle_name.append(file)
    angle_name.sort()
    # data
    incidence = json.loads(open(inputdir + angle_name[-1]).read())
    inc_data = np.array(incidence["data"])
    azimuth = json.loads(open(inputdir + angle_name[0]).read())
    az_data = np.array(azimuth["data"])
    
    #Calculate LOS simulated displacement
    inc_data *= np.pi/180.
    # heading angle
    head_angle = az_data
    head_angle *= np.pi/180.
    
    az_angle = 270 * np.pi/180.
    # construct design matrix
    coefficient_d = np.cos(inc_data)
    coefficient_e = np.sin(inc_data) * np.sin(head_angle - az_angle)
    coefficient_n = np.sin(inc_data) * np.cos(head_angle - az_angle)
    # los displacement
    disp_los = disp_data_d * coefficient_d - disp_data_n * coefficient_e - disp_data_e * coefficient_n
    # write
    disp_name = "_".join(disp_sim_name[0].split('_')[0:2]) + '_LOS'
    outdir = "".join(inps.outdir)
    print(disp_name)
    print(outdir)
    write_json(disp_los,disp_name,disp_sim_d,outdir)
    #geocode
    geocode(disp_los,disp_name,disp_sim_d,outdir)

    
######################################################################################
if __name__ == '__main__':
    main()
