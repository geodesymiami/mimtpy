#!/usr/bin/env python3
######################################################################################################################################
# Program is used for extract three directions psgrn/pscmp results to have same area with MintPy result's area                       #
#    and generate incidence and azimuth *.json files with same spatial resolution and area for calculating simulated LOS deformation #
# Author: Lv Xiaoran                                                                                                                 #
# Created: January  2020                                                                                                             #
######################################################################################################################################

import os
import sys
import argparse
import numpy as np
import json
from mintpy.utils import readfile
from osgeo import gdal, osr
######################################################################################
EXAMPLE = """example:
    extract_simresults_to_Insar_region.py $MODELOUT/pscmp/Bogd_test2/model4/poseis-10y_20y/ $SCRATCHDIR/BogdSenDT106/mintpy/S1_IW123_106_0441_0447_20161013_20191215.he5 --outdir $MODELOUT/pscmp/Bogd_test2/model4/poseis-10y_20y/subtract/ 
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Extract three directions *.json files and incidence/azimuth *.json for calculating Los deformation',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('input_json', nargs=1, help='directory stored *.json files. \n')

    parser.add_argument('input_HDFEOS', nargs=1, help='directory stored *.he5 files. \n')
    
    parser.add_argument('--outdir',dest='outdir',nargs='?',default=os.getenv('MODELOUT'), help='output directory')

    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps    

def resample(angle,times,row_resample,sample_resample):
    """resample incidence/azimuth to simulated data resolution"""
    angle_resam = np.zeros((row_resample,sample_resample),dtype=np.float32) * np.nan
    rows_h5,samples_h5 = angle.shape
    for i in range(0,rows_h5,times):
        for j in range(0,samples_h5,times):
            ii = int(i/times)
            jj = int(j/times)
            #print(ii)
            #print(jj)
            if i!= np.floor(rows_h5/times)*times and j!= np.floor(samples_h5/times)*times:
                angle_resam[ii, jj] = np.nanmean(angle[i:i + times,j:j + times])
            elif i!= np.floor(rows_h5/times)*times and j == np.floor(samples_h5/times)*times:
                angle_resam[ii, jj] = np.nanmean(angle[i:i + times, j:])
            elif i == np.floor(rows_h5/times)*times and j!= np.floor(samples_h5/times)*times:
                angle_resam[ii, jj] = np.nanmean(angle[i:,j:j + times])
            elif i == np.floor(rows_h5/times)*times and j== np.floor(samples_h5/times)*times:
                angle_resam[ii, jj] = np.nanmean(angle[i:,j:])       
    return angle_resam
    
def read_hdfeos(inps):
    """read HDFEOS files"""
    #read HDFEOS file
    # metadata
    atr = readfile.read_attribute("".join(inps.input_HDFEOS))
    # data
    incidence = readfile.read("".join(inps.input_HDFEOS), datasetName='/HDFEOS/GRIDS/timeseries/geometry/incidenceAngle')[0]
    azimuth = readfile.read("".join(inps.input_HDFEOS), datasetName='/HDFEOS/GRIDS/timeseries/geometry/azimuthAngle')[0]
    return incidence,azimuth,atr

def geocode(disp_data,displacement,disp_name,outdir):
    """geocode step"""
    ymin = displacement["ul_lat"] + displacement["lat_step"] * displacement["rows"]
    ymax = displacement["ul_lat"]
    xmin = displacement["ul_lon"]
    xmax = displacement["ul_lon"] + displacement["lon_step"] * displacement["samples"]
    nrows = displacement["rows"]
    ncols = displacement["samples"]
    yres = abs(displacement["lat_step"])
    xres = displacement["lon_step"]
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

def subtract_data(inps,disp_sim_name,incidence,azimuth,atr):
    """subtract corresponding region from simulated displacement accroding InSAR results"""
    # read 2d simulated displacement data
    disp_data_name = "".join(inps.input_json) + '/' + disp_sim_name
    disp_sim = json.loads(open(disp_data_name).read())
    # attribuate
    ul_lon_sim = disp_sim["ul_lon"]
    ul_lat_sim = disp_sim["ul_lat"]
    
    lat_step_sim = disp_sim["lat_step"]
    lon_step_sim = disp_sim["lon_step"]
    # displacement data
    disp_data = np.array(disp_sim["disp_data"])
    
    # read attribute of HDFEOS
    ul_lat_h5 = float(atr['Y_FIRST'])
    ul_lon_h5 = float(atr['X_FIRST'])
    samples_h5 = int(atr["WIDTH"])
    rows_h5 = int(atr["LENGTH"])
    
    lat_step_h5 = float(atr["Y_STEP"][1:])
    lon_Step_h5 = float(atr["X_STEP"])

    #locate Insar position in Simulated field and subtract corresponding region from simulated displacement accroding InSAR results 
    row_locate = int((ul_lat_h5 - ul_lat_sim) / disp_sim["lat_step"] + 0.5)
    sample_locate = int((ul_lon_h5 - ul_lon_sim) / disp_sim["lon_step"] + 0.5)
    
    # judge whether simulated data resolution equal to InSAR data resolution 
    if lat_step_h5 == lat_step_sim:
        print("equal resolution!")
        rows_offset = rows_h5
        samples_offset = samples_h5
        disp_sub = disp_data[row_locate:row_locate + rows_offset, sample_locate:sample_locate + samples_offset]
        position = np.isnan(incidence)
        disp_sub[position] = np.nan  
    else:
        print("not equal resolution")
        # generally lat_step_sim is larger than lat_step_h5
        times = int(np.abs(lat_step_sim / lat_step_h5))
        rows_offset = int(np.ceil(rows_h5 / times))
        samples_offset = int(np.ceil(samples_h5 / times))
        disp_sub = disp_data[row_locate:row_locate + rows_offset, sample_locate:sample_locate + samples_offset]
        
        # resample incidence and azimuth data
        print('resample incidence')
        incidence = resample(incidence,times,rows_offset,samples_offset)
        
        print('resample azimuth')
        azimuth = resample(azimuth,times,rows_offset,samples_offset)
        position = np.isnan(incidence)
        disp_sub[position] = np.nan
        
    # write subtracted simulated displacement
    outdir = "".join(inps.outdir)
    displacement = {"ul_lon": ul_lon_h5, "ul_lat": ul_lat_h5, "lat_step": lat_step_sim, "lon_step": lon_step_sim, "samples": samples_offset, "rows": rows_offset, "data": disp_sub.tolist()}
    outname = disp_sim_name.split('.')[0] +'_subtract.json'
    open(outdir + '/' + outname, "w").write(json.dumps(displacement))
    
    # geocode data
    disp_name = disp_sim_name.split('.')[0]+'_subtract'
    geocode(disp_sub,displacement,disp_name,outdir)
    
    return incidence, azimuth, atr, disp_sim, rows_offset, samples_offset

def write_inci_az(inps,incidence,azimuth,atr,disp_sim,rows_offset,samples_offset):
    """write incidence file and azimuth file""" 
    lat_step_sim = disp_sim["lat_step"]
    lon_step_sim = disp_sim["lon_step"]
    ul_lat_h5 = float(atr['Y_FIRST'])
    ul_lon_h5 = float(atr['X_FIRST'])
    samples_h5 = int(atr["WIDTH"])
    rows_h5 = int(atr["LENGTH"])
    outdir = "".join(inps.outdir)

    inc = {"ul_lon": ul_lon_h5, "ul_lat": ul_lat_h5, "lat_step": lat_step_sim, "lon_step": lon_step_sim, "samples": samples_offset, "rows": rows_offset, "data": incidence.tolist()}
    outname = 'incidence_resam.json'
    open(outdir + '/' + outname, "w").write(json.dumps(inc))
    
    az = {"ul_lon": ul_lon_h5, "ul_lat": ul_lat_h5, "lat_step": lat_step_sim, "lon_step": lon_step_sim, "samples": samples_offset, "rows": rows_offset, "data": azimuth.tolist()}
    outname = 'azimuth_resam.json'
    open(outdir + '/' + outname, "w").write(json.dumps(az))    

######################################################################################
def main(iagrs=None):
    inps = cmd_line_parse(iagrs)
    # read HDFEOS file
    incidence_or,azimuth_or,atr = read_hdfeos(inps)
        
    # find displacement with *eu.json,*ns.json,*up.json files
    path_list = os.listdir("".join(inps.input_json))
    json_file = []
    for filename in path_list:
        if os.path.splitext(filename)[1] == '.json':
            json_file.append(filename)
    json_file.sort()
    print(json_file)
    for json in json_file:
         incidence,azimuth,atr,disp_sim,rows_offset,samples_offset = subtract_data(inps,json,incidence_or,azimuth_or,atr)
    write_inci_az(inps,incidence,azimuth,atr,disp_sim,rows_offset,samples_offset)
######################################################################################
if __name__ == '__main__':
    main()
