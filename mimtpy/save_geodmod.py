#!/usr/bin/env python3
#################################################################
# Program is used for preparing data for geodmod software       #
# Author: Lv Xiaoran                                            #
# Created: June 2020                                            #
#################################################################

import os
import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
import shutil
import numpy as np
from PIL import Image

import mintpy
import mintpy.workflow  #dynamic import for modules used by pysarApp workflow
from mintpy.utils import readfile, writefile,utils as ut
from mimtpy.utils import multitrack_utilities as mu
######################################################################################
EXAMPLE = """example:
  save_geodmod.py S1_IW23_026_0108_0113_20171117_XXXXXXXX.he5 -s 20171128 -e 20181210 
  
  save_geodmod.py S1_IW23_026_0108_0113_20171117_XXXXXXXX.he5 -b 26.0 27.5 65.0 66.0 -s 20171128 -e 20181210 
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Prepare data for Geodmod software',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs='?', help='ascending or descending files\n')

    parser.add_argument('-b', '--bbox', dest='SNWE', type=float, nargs=4, metavar=('S', 'N', 'W', 'E'),
                        help='Bounding box of area to be geocoded.\n' +
                        'Include the uppler left corner of the first pixel' +
                        '    and the lower right corner of the last pixel')
    parser.add_argument('-s','--startDate',dest='startDate',nargs='?',
                        help='date1 of timeseires to be converted.The default is the StartDate')
    parser.add_argument('-e','--endDate',dest='endDate',nargs='?',
                        help='date2 of timeseries to be converted.The default is the EndDate')
    
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps
    
def write_rsc_file(inps,in_file,out_file):
    """ write rsc file for Geodmod just estract several properities from rsc file"""
    # read file
    meta = readfile.read_roipac_rsc(in_file)
    # initiate dict
    rsc = dict()
    rsc['FILE_DIR'] = os.getenv('pwd')
    rsc['FILE_LENGTH'] = meta["FILE_LENGTH"]
    rsc['WIDTH'] = meta["WIDTH"]
    rsc['XMIN'] = 0
    rsc['XMAX'] = int(meta["WIDTH"]) - 1
    rsc['YMIN'] = 0
    rsc['YMAX'] = int(meta["FILE_LENGTH"]) - 1
    rsc['X_FIRST'] = float(meta["X_FIRST"])
    rsc['Y_FIRST'] = float(meta["Y_FIRST"])
    rsc['X_STEP'] = float(meta["X_STEP"])
    rsc['Y_STEP'] = float(meta["Y_STEP"])
    rsc['X_UNIT'] = 'degrees'
    rsc['Y_UNIT'] = 'degrees'
    rsc['RLOOKS'] = meta["RLOOKS"]
    rsc['ALOOKS'] = meta["ALOOKS"]
    rsc['Z_OFFSET'] = 0
    rsc['Z_SCALE'] = 1
    rsc['PROJECTION'] = 'LATLON'
    rsc['DATE12'] = '111111-222222'
    # write rsc file
    writefile.write_roipac_rsc(rsc, out_file, print_msg=True)
    return out_file 

def dem_jpeg(dem_file):
    """generate dem.jepg file based on Yunjun's code"""
    out_file = dem_file+'.jpeg'
    rsc_file = out_file+'.rsc'
    shutil.copy2(dem_file+'.rsc', rsc_file)
    # read data
    dem = readfile.read(dem_file)[0]
    print('dem.shape:',dem.shape)
    # figure size
    ds_shape = tuple(reversed(dem.shape))
    fig_dpi = 300
    fig_size = [i / fig_dpi for i in ds_shape]
    print('fig_size:',fig_size)
    # color range
    disp_min = np.nanmin(dem) - 4000
    disp_max = np.nanmax(dem) + 2000
    # prepare shaded relief
    ls = LightSource(azdeg=315, altdeg=45)
    dem_shade = ls.shade(dem, vert_exag=0.3, cmap=plt.get_cmap('gray'), vmin=disp_min, vmax=disp_max)
    dem_shade[np.isnan(dem_shade[:, :, 0])] = np.nan
    print('dem_shade.shape:',dem_shade.shape)
    # plot
    fig, ax = plt.subplots(figsize=fig_size)
    ax.imshow(dem_shade, interpolation='spline16', origin='upper')
    # get rid of whitespace on the side
    ax.axis('off')
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    fig.subplots_adjust(left=0,right=1,bottom=0,top=1)
    # output
    print('save figure to file {}'.format(out_file))
    plt.savefig(out_file, transparent=True, dpi=300, pad_inches=0.0)
    
    #resize to desired size(FA 8/19, unclear why size is wrong)
    im = Image.open(out_file)
    im_out = im.resize(dem.shape, Image.NEAREST)
    im_out.save(out_file)
    
    #plt.show()

def subset_data_based_bbox(inps,dataset):
    """return the row_no,sample_no and rows and samples"""
    # metadata
    atr = readfile.read_attribute("".join(inps.file))
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
    
    # subset data
    data,atr = readfile.read(dataset)
    atr['LENGTH'] = str(rows)
    atr['WIDTH'] = str(samples)
    writefile.write(data, out_file=dataset, metadata=atr)

    return

def mask_filter(inps,dataset):
    """mask data"""
    print("mask {} file".format(dataset))
    maskfile = readfile.read("".join(inps.file),datasetName='/HDFEOS/GRIDS/timeseries/quality/mask')[0]
    data,atr = readfile.read(dataset)
    data[maskfile == 0] = np.nan
    writefile.write(data, out_file=dataset, metadata=atr)

def process_HDFEOS(inps):
    """process *.he5 file"""
    atr_asc = inps.file
        
    #save dataset of unw cor and dem
    filename, extension = mu.seprate_filename_extension("".join(atr_asc))[1:3]
    output_unw = "".join(['geo_',inps.startDate,'_',inps.endDate,'.unw'])
    cmd_args = ['../'+filename+extension, "".join(['displacement-',inps.startDate,'_',inps.endDate]), '-o', output_unw]
    print("save_roipac.py", cmd_args)
    asct_str = mu.seperate_str_byspace(cmd_args)
    os.system(mu.seperate_str_byspace(['save_roipac.py', asct_str.split()]))   

    output_cor = "".join(['geo_',inps.startDate,'_',inps.endDate,'.cor'])
    cmd_args = ['../'+filename+extension, 'temporalCoherence', '-o', output_cor]
    print("save_roipac.py", cmd_args)
    asct_str = mu.seperate_str_byspace(cmd_args)
    completion_status=os.system(mu.seperate_str_byspace(['save_roipac.py', asct_str.split()])) 
    
    output_dem = 'srtm.dem'
    cmd_args = ['../'+filename+extension, 'height', '-o', output_dem]
    print("save_roipac.py", cmd_args)
    asct_str = mu.seperate_str_byspace(cmd_args)
    os.system(mu.seperate_str_byspace(['save_roipac.py', asct_str.split()]))
    
    mask_filter(inps,output_unw)
    mask_filter(inps,output_cor)
    #mask_filter(inps,output_dem)

    if inps.SNWE:
        print('Subset data based on bbox')
        subset_data_based_bbox(inps,output_unw)
        subset_data_based_bbox(inps,output_cor)
        subset_data_based_bbox(inps,output_dem)
   
    
######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)   
    print(inps) 
    # process HDFEOS data
    process_HDFEOS(inps)

    # rename *.rsc1 to *.rsc
    outfile = mu.seperate_str_byspace(['srtm.dem' + '.rsc'])
    write_rsc_file(inps,outfile,mu.seperate_str_byspace(['srtm.dem' +'.rsc1']))
    os.remove(outfile)
    print('rename *.rsc1 to *.rsc')
    os.rename(mu.seperate_str_byspace(['srtm.dem' +'.rsc1']),outfile)
    
    # generate dem.jpeg
    dem_jpeg('srtm.dem')

######################################################################################
if __name__ == '__main__':
    main()
