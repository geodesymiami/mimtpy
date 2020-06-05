#!/usr/bin/env python3
######################################################################################################
# Program is used for generate h5 and geotiff(WGS projection) files from mintpy HDFEOS file          #
# Author: Lv Xiaoran                                                                                 #
# Created: February  2020                                                                            #
######################################################################################################

import os
import argparse
import h5py
import numpy as np
from osgeo import gdal, osr
from mintpy.utils import readfile, writefile, ptime
from mintpy import view
from mimtpy.utils import multitrack_utilities
######################################################################################
EXAMPLE = """example:
    HDFEOS_to_geotif.py $SCRATCHDIR/BogdSenDT106/mintpy/S1_IW123_106_0441_0447_20161013_20191215.he5 mask
    HDFEOS_to_geotif.py $SCRATCHDIR/BogdSenDT106/mintpy/S1_IW123_106_0441_0447_20161013_20191215.he5 temporalCoherence --outdir $SCRATCHDIR/BogdSenDT106/
    HDFEOS_to_geotif.py $SCRATCHDIR/BogdSenDT106/mintpy/S1_IW123_106_0441_0447_20161013_20191215.he5 displacement --date 20191215 --mask
    
    HDFEOS_to_geotif.py $SCRATCHDIR/BogdSenDT106/mintpy/S1_IW123_106_0441_0447_20161013_20191215.he5 velocity --outdir $SCRATCHDIR/BalochistanSenDT151/test
     
    HDFEOS_to_geotif.py $SCRATCHDIR/BogdSenDT106/mintpy/S1_IW123_106_0441_0447_20161013_20191215.he5 displacement --date 20161013_20191215 --mask
    
    HDFEOS_to_geotif.py $SCRATCHDIR/BogdSenDT106/mintpy/S1_IW123_106_0441_0447_20161013_20191215.he5 displacement --date 20161013_20191215 --mask --bbox 45 46 101 102
    
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Generate *.h5 and *.tiff file with WGS projection based HDFEOS file',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('input_HDFEOS', nargs=1, help='directory stored *.he5 files. \n')
    
    parser.add_argument('attribute', nargs=1, help='chosen attribute. \n')
    
    parser.add_argument('--date', dest='date',nargs=1, help='date1 or date1_2 for displacement. \n')
    
    parser.add_argument('--mask',action='store_true', default=False, help='whether mask outputfile. \n')

    parser.add_argument('--bbox', dest='SNWE', type=float, nargs=4, metavar=('S', 'N','W','E'),
                        help='Bounding box of area to be geocoded.\n'+
                        'Include the uppler left corner of the first pixel' + 
                        'and the lower right corner of the last pixel')

    parser.add_argument('--outdir',dest='outdir',nargs='?',default=os.getenv('SCRATCHDIR'), help='output directory')

    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps    

def located_dataset_based_attribut(inps):
    """find the dataset based on the given attribution"""
    attr = "".join(inps.attribute)
    geometry_attr = ['azimuthAngle','height','incidenceAngle','latitude','longitude','shadowMask','slantRangeDistance']
    observation_attr = ['bperp','date','displacement','velocity']
    quality_attr = ['mask','temporalCoherence']
    
    #define the dataset
    if attr in geometry_attr:
        dataset='/HDFEOS/GRIDS/timeseries/geometry/'
    if attr in quality_attr:
        dataset = '/HDFEOS/GRIDS/timeseries/quality/'
    if attr in observation_attr:
        dataset = '/HDFEOS/GRIDS/timeseries/observation/'
        if attr == 'bperp' or attr == 'date':
            raise Exception("ERROR! This attributaion is 1D array!")

    return dataset
    #azimuth = readfile.read("".join(inps.input_HDFEOS), datasetName='/HDFEOS/GRIDS/timeseries/geometry/azimuthAngle')[0]

def extract_data_based_bbox(inps):
    """return the row_no,sample_no and rows and samples"""
    # metadata
    atr = readfile.read_attribute("".join(inps.input_HDFEOS))
    ul_lat = float(atr['Y_FIRST'])
    ul_lon = float(atr['X_FIRST'])
    lat_step = float(atr["Y_STEP"])
    lon_step = float(atr["X_STEP"])
    # bbox
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
    
    row = int((user_lat0 - ul_lat) / lat_step + 0.5)
    sample = int((user_lon0 - ul_lon) / lon_step + 0.5)    
    rows = int((user_lat1 - user_lat0) / lat_step + 0.5) + 1
    samples = int((user_lon1 - user_lon0) / lon_step + 0.5) + 1 
    return row, sample, rows,samples

def extract_data(inps,dataset,outdir):
    """extract data from HDFEOS file based on the given attribute"""
    # read HDFEOS file
    # metadata
    atr = readfile.read_attribute("".join(inps.input_HDFEOS))
    
    attr = "".join(inps.attribute)
    # read 2d data
    if attr == 'displacement' or attr == 'velocity':
        if inps.date == None:
            date1 = atr['START_DATE']
            date2 = atr['END_DATE']
            #raise Exception("ERROR! Date for displacement must be given!")
        else:
            # date1 and date2
            if '_' in "".join(inps.date):
                date1, date2 = ptime.yyyymmdd("".join(inps.date).split('_'))
            else:
                date1 = atr['START_DATE']
                date2 = ptime.yyyymmdd("".join(inps.date))
            #date12 = '{}_{}'.format(date1, date2)
        if attr == 'displacement':            
            # read / prepare data
            slice_list = readfile.get_slice_list("".join(inps.input_HDFEOS))
            # read/prepare data
            dname = 'displacement'
            slice_name1 = view.check_dataset_input(slice_list, '{}-{}'.format(dname, date1))[0][0]
            slice_name2 = view.check_dataset_input(slice_list, '{}-{}'.format(dname, date2))[0][1]
            data = readfile.read("".join(inps.input_HDFEOS), datasetName=slice_name1)[0]
            data -= readfile.read("".join(inps.input_HDFEOS), datasetName=slice_name2)[0]
            data_name = '{}_{}_{}'.format(attr,date1,date2)
            if inps.mask:
                print("mask file")
                maskfile = readfile.read("".join(inps.input_HDFEOS),datasetName='/HDFEOS/GRIDS/timeseries/quality/mask')[0]
                data[maskfile == 0] = np.nan
            outfile = outdir + '/' + data_name + '.h5'
            atr['FILE_TYPE'] = '.unw'
            writefile.write(data, out_file=outfile, metadata=atr)
            #print('converting range to phase')
            #data *= range2phase
            #if inps.ref_yx:
            #    data -= data[inps.ref_yx[0], inps.ref_yx[1]]
        elif attr == 'velocity':
            dname = 'displacement'
            data_timeseries = readfile.read("".join(inps.input_HDFEOS),datasetName = dataset+dname)[0]
            if inps.mask:
                print("mask file")
                maskfile = readfile.read("".join(inps.input_HDFEOS),datasetName='/HDFEOS/GRIDS/timeseries/quality/mask')[0]
                data_timeseries[:,maskfile == 0] = np.nan
            bperp_date = h5py.File("".join(inps.input_HDFEOS),'r')
            data_bperp = bperp_date[(dataset+'bperp')]
            data_date = bperp_date[(dataset+'date')]
            ts_data_name = 'timeseries_ERA_demErr.h5'
            outfile = outdir + '/' + ts_data_name
            # write to HDF5 file
            dsDict = dict()
            dsDict['bperp'] = data_bperp
            dsDict['date'] = data_date
            dsDict['timeseries'] = data_timeseries
            atr['FILE_TYPE'] = attr 
            writefile.write(dsDict, out_file=outfile, metadata=atr)    
            # run timeseries2velocity.py
            output_vel = 'velocity_' + date1 + '_' + date2 + '.h5'
            data_name = output_vel.split('.')[0]
            os.chdir(outdir)
            completion_status = os.system(multitrack_utilities.seperate_str_byspace(['timeseries2velocity.py', outfile, '--start-date', date1, '--end-date', date2, '--exclude', '../exclude_date.txt', '--output', output_vel]))
            if completion_status == 1:
                raise Exception('error when generate velocity!')  
            date_vel = readfile.read(output_vel, datasetName='velocity')[0]
            data = date_vel
            os.chdir('../')
    else:
        data = readfile.read("".join(inps.input_HDFEOS),datasetName=dataset+attr)[0]
        outfile = outdir + '/' + attr + '.h5'
        data_name = attr
        atr['FILE_TYPE'] = attr 
        if inps.mask:
            print("mask file")
            maskfile = readfile.read("".join(inps.input_HDFEOS),datasetName='/HDFEOS/GRIDS/timeseries/quality/mask')[0]
            data[maskfile == 0] = np.nan
        writefile.write(data, out_file=outfile, metadata=atr) 
    # mask data
    #if inps.mask:
    #    print("mask file")
    #    maskfile = readfile.read("".join(inps.input_HDFEOS),datasetName='/HDFEOS/GRIDS/timeseries/quality/mask')[0]
    #    data[maskfile == 0] = np.nan
        
    return data, data_name, atr

def geocode(data,data_name,atr,outdir):
    """geocode step"""
    # calculate geo information
    # read attribute of HDFEOS
    ul_lat = float(atr['Y_FIRST'])
    ul_lon = float(atr['X_FIRST'])
    samples = int(atr["WIDTH"])
    rows = int(atr["LENGTH"])
    lat_step = float(atr["Y_STEP"])
    lon_step = float(atr["X_STEP"])
    
    # prepare geometry information for gdal
    ymax = ul_lat
    ymin = ul_lat + lat_step * (rows - 1)
    xmin = ul_lon
    xmax = ul_lon + lon_step * (samples - 1)
    nrows = rows
    ncols = samples
    xres = lon_step
    yres = abs(lat_step)
    
    # output
    output_tif = outdir +'/' + data_name + '.tiff'
    print('Start calculating geotiff data {}.\n'.format(output_tif))
    geotransform = [xmin,xres,0,ymax,0,-yres]
    raster = gdal.GetDriverByName('GTiff').Create(output_tif,ncols,nrows,1,gdal.GDT_Float32)
    raster.SetGeoTransform(geotransform)
    srs=osr.SpatialReference()
    #srs.ImportFromEPSG(32638) #wgs84 UTM 38N
    srs.ImportFromEPSG(4326) #WGS 84 - WGS84 - World Geodetic System 1984, used in GPS
    raster.SetProjection(srs.ExportToWkt())
    raster.GetRasterBand(1).WriteArray(data)
    raster.FlushCache()
    print('Finishing calculating geotiff data {}.\n'.format(output_tif))

def geocode_bbox(data,data_name,inps,atr,outdir):
    """geocode step"""
    # calculate geo information
    # read attribute of HDFEOS
    lat_step = float(atr["Y_STEP"])
    lon_step = float(atr["X_STEP"])
    
    # prepare geometry information for gdal
    ymax = float(inps.SNWE[1])
    xmin = float(inps.SNWE[2])
    nrows,ncols = data.shape
    xres = lon_step
    yres = abs(lat_step)
    
    # output
    output_tif = outdir +'/' + data_name + '.tiff'
    geotransform = [xmin,xres,0,ymax,0,-yres]
    raster = gdal.GetDriverByName('GTiff').Create(output_tif,ncols,nrows,1,gdal.GDT_Float32)
    raster.SetGeoTransform(geotransform)
    srs=osr.SpatialReference()
    #srs.ImportFromEPSG(32638) #wgs84 UTM 38N
    srs.ImportFromEPSG(4326) #WGS 84 - WGS84 - World Geodetic System 1984, used in GPS
    raster.SetProjection(srs.ExportToWkt())
    raster.GetRasterBand(1).WriteArray(data)
    raster.FlushCache()

######################################################################################
def main(iagrs=None):
    inps = cmd_line_parse(iagrs)
    # read HDFEOS file
    dataset = located_dataset_based_attribut(inps)
   
    outdir = "".join(inps.outdir) 
    data,data_name,atr = extract_data(inps,dataset,outdir)
    # bbox
    if inps.SNWE == None:
        print('\n generate geotiff')
        #geocode data
        geocode(data,data_name,atr,outdir)    
    else:
        row_No,sample_No,rows,samples = extract_data_based_bbox(inps)
        data_bbox = data[row_No : row_No + rows,sample_No : sample_No + samples]
        data_name_bbox = data_name + '_subarea'
        #geocode data
        geocode_bbox(data_bbox,data_name_bbox,inps,atr,outdir)    
######################################################################################
if __name__ == '__main__':
    main()
