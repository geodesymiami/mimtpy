#!/usr/bin/env python3
#################################################################
# Program is used for concatenat data velocity loop             #
# Author: Lv Xiaoran                                            #
# Created: October 2021                                         #
#################################################################

import os
import argparse
import numpy as np
import shutil

import mintpy
from mintpy.utils import readfile, writefile

import mimtpy
import mimtpy.workflow
from mimtpy.utils import multitrack_utilities as mu
######################################################################################
EXAMPLE = """example:
    loop_concatenation.py --project KokoxiliBig --tracks SenDT19 --lat 30 41
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Concatenate data',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('--project', nargs=1, type=str, help='project name\n')

    parser.add_argument('--tracks', type=str, nargs=1, help='track name\n')

    parser.add_argument('--lat', nargs=2, type=int, help='start and end lat of chunks\n')
    
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps

######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)   
    
    pro_name = inps.project[0]
    track_name = inps.tracks[0]

    lat_min = int(inps.lat[0])
    lat_max = int(inps.lat[1]) + 1

    # make sure the existence for project folders
    dataset_dirs = []
    for lat in np.arange(lat_min, lat_max):
        projects = pro_name + str(lat) + track_name
        dataset_dir = os.path.abspath(os.path.join(os.getenv('SCRATCHDIR'),projects,'mintpy'))
        if not os.path.isdir(dataset_dir):
            raise Exception('Error! No such dir : {}'.format(dataset_dir))
        else:
            dataset_dirs.append(dataset_dir)
    
    # extract velocity for each dataset
    for dataset in dataset_dirs:
        os.chdir(dataset)
        print('\nGo to project dir:', dataset)
        # creat velocity dir
        outdir = os.path.abspath(os.path.join(dataset, 'velocity'))
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        print('\nthe output dir for velocity is {}.\n'.format(outdir))
 
        # find HDFEOS file
        HDFEOS_file = mu.find_HDFEOS_fullname(dataset)
        print('\nThe HDFEOS file is {}'.format(HDFEOS_file))
    
        scp_args = [HDFEOS_file]
        scp_args = mu.seperate_str_byspace(scp_args)
        
        # extract mask
        maskfile = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/quality/mask')[0] 
        
        # run HDFEOS_to_geotiff.py
        print('timeseries2velocity.py',scp_args)
        mintpy.timeseries2velocity.main(scp_args.split())
        shutil.move(os.path.join(dataset, 'velocity.h5'), os.path.join(outdir, 'velocity.h5'))
        # mask the velocity
        os.chdir(outdir)
        vel_data, vel_atr = readfile.read('velocity.h5')
        vel_data[maskfile == 0] = np.nan
        writefile.write(vel_data, out_file='velocity.h5', metadata=vel_atr)

    # concatenate chunks for one track
    # generate output dir 
    track_dir = os.path.abspath(os.path.join(os.getenv('SCRATCHDIR') + '/' +pro_name + track_name + '/mimtpy/' + 'velocity/'))
    if not os.path.isdir(track_dir):
        os.makedirs(track_dir)
    print('\nthe output dir for concatenation is {}.\n'.format(track_dir))
    # start the loop
    temp_files = []
    pro_num = len(dataset_dirs)
    lat_range = np.arange(lat_min, lat_max)
    for i in np.arange(pro_num - 1):
        if i == 0:
            temp_name = 'velocity_lat_' + str(lat_range[0]) + '_' + str(lat_range[1])
            scp_args = [dataset_dirs[0] + '/velocity/' + 'velocity.h5', dataset_dirs[1]  + '/velocity/' + 'velocity.h5', '--output', temp_name,  '--outdir', track_dir + '/']
            
        else:
            temp_name = temp_files[i-1] + '_' + str(lat_range[i+1])
            scp_args = [track_dir + '/' + temp_files[i-1] +'.h5', dataset_dirs[i+1]  + '/velocity/' + 'velocity.h5', '--output', temp_name, '--outdir', track_dir + '/']

        temp_files.append(temp_name)
        scp_args = mu.seperate_str_byspace(scp_args)
        
        # run track_offset.py
        print('track_offset.py',scp_args)
        mimtpy.track_offset.main(scp_args.split())   
       
######################################################################################
if __name__ == '__main__':
    main()
