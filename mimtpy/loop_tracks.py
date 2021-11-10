#!/usr/bin/env python3
#################################################################
# Program is used for concatenat data velocity loop             #
# Author: Lv Xiaoran                                            #
# Created: October 2021                                         #
#################################################################

import os
import re
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
    loop_tracks.py --project KokoxiliBig --tracks SenDT19 SenDT150 SenDT48
    loop_tracks.py --project KokoxiliBig --tracks SenDT19 SenDT150 SenDT48 --chunk
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Concatenate data',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('--project', nargs=1, type=str, help='project name\n')

    parser.add_argument('--tracks', type=str, nargs='+', help='tracks name\n')
    
    parser.add_argument('--chunk', action='store_true', default=False, help='whether concatenate chunks for each track \n')

    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps

def search_chunks_number(file_dir, project, track_name):   
    L=[]
    chunk_number = []
    for root, dirs, files in os.walk(file_dir):  
        for dir in dirs:  
            ret = re.match(project + r'(\d{2})' + track_name, dir)
            if ret:
                #print(ret, ret[1])
                L.append(os.path.join(root, dir, 'mintpy'))
                chunk_number.append(dir.replace(project, '').replace(track_name,''))
    chunk_number = [int(x) for x in chunk_number]
    return L, chunk_number

def concatenation_chunks(pro_name, track_name):
    #dataset_dirs = []
    # search for chunks number and name
    scratch_dir = os.getenv('SCRATCHDIR')
    print(scratch_dir)
    print(pro_name)
    print(track_name)
    dataset_dirs, chunk_number = search_chunks_number(scratch_dir, pro_name, track_name)
    dataset_dirs.sort()
    chunk_number.sort()   
    
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
    track_dir = os.path.abspath(os.path.join(os.getenv('SCRATCHDIR') + '/' +pro_name[0: -5] + 'Big' + track_name + '/mimtpy/' + 'velocity/'))
    if not os.path.isdir(track_dir):
        os.makedirs(track_dir)
    print('\nthe output dir for concatenation is {}.\n'.format(track_dir))
    
    # start the loop
    temp_files = []
    pro_num = len(dataset_dirs)
    #lat_range = np.arange(lat_min, lat_max)
    lat_range = chunk_number
    for i in np.arange(pro_num - 1):
        if i == 0:
            temp_name = 'velocity_lat_' + str(lat_range[0]) + '_' + str(lat_range[1])
            scp_args = [dataset_dirs[0] + '/velocity/' + 'velocity.h5', dataset_dirs[1]  + '/velocity/' + 'velocity.h5', '--output', temp_name,  '--outdir', track_dir + '/']
            
        elif i > 0 and i < (pro_num - 1 - 1):
            temp_name = temp_files[i-1] + '_' + str(lat_range[i+1])
            scp_args = [track_dir + '/' + temp_files[i-1] +'.h5', dataset_dirs[i+1]  + '/velocity/' + 'velocity.h5', '--output', temp_name, '--outdir', track_dir + '/']
        else:
            temp_name = 'velocity_track'
            scp_args = [track_dir + '/' + temp_files[i-1] +'.h5', dataset_dirs[i+1]  + '/velocity/' + 'velocity.h5', '--output', temp_name, '--outdir', track_dir + '/']

        temp_files.append(temp_name)
        scp_args = mu.seperate_str_byspace(scp_args)
        # run track_offset.py
        print('track_offset.py',scp_args)
        mimtpy.track_offset.main(scp_args.split())   

    for dataset in dataset_dirs:
        os.chdir(dataset)
        print('\nGo to project dir:', dataset)
        # delete velocity folder
        print('\n Delete velocity folder from %s' % dataset)
        outdir = os.path.abspath(os.path.join(dataset, 'velocity'))
        shutil.rmtree(outdir)
 
def concatenation_tracks(track_names, pro_name):
    velocity_track = 'velocity_track.h5'
 
    # generate output dir
    track_name = track_names[0]
    pro_dir = os.path.abspath(os.path.join(os.getenv('SCRATCHDIR') + '/' + pro_name[0: -5] + 'Big' +  track_name[0:5] + '/mimtpy/' + 'velocity/'))
    if not os.path.isdir(pro_dir):
        os.makedirs(pro_dir)
    print('\nthe output dir for concatenation is {}.\n'.format(pro_dir))
    
    track_number = len(track_names)
    temp_Files = []
    pro_dir_tmp = os.path.abspath(os.path.join(os.getenv('SCRATCHDIR') + '/' +pro_name[0: -1] + 'Big'))
    for i in np.arange(track_number - 1):  
        if i == 0:
            temp_name = 'velocity_' + track_names[0] + '_' + track_names[1]
            scp_args = [pro_dir_tmp + track_names[0] + '/mimtpy/velocity/velocity_track.h5', pro_dir_tmp + track_names[1] + '/mimtpy/velocity/velocity_track.h5', '--output', temp_name, '--outdir', pro_dir + '/']
        else:
            temp_name = temp_Files[i-1] + '_' + track_names[i+1]
            scp_args = [pro_dir + '/' + temp_Files[i-1] + '.h5', pro_dir_tmp + track_names[i+1] + '/mimtpy/velocity/velocity_track.h5', '--output', temp_name, '--outdir', pro_dir + '/']

        temp_Files.append(temp_name) 
        scp_args = mu.seperate_str_byspace(scp_args)

        # run track_offset.py
        print('track_offset.py',scp_args)
        mimtpy.track_offset.main(scp_args.split())        
######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)   
    
    pro_name = inps.project[0]
    track_names = inps.tracks

    # firstly concatenation chunks of each track
    # the output dir is project_name/mintpy/velocity/, for example: KokoxiliBigSenDT19/mintpy/velocity/
    # make sure the existence for project folders
    if inps.chunk:
        for track_name in track_names:
            concatenation_chunks(pro_name, track_name)

    # secondly concatenation tracks 
    concatenation_tracks(track_names, pro_name)
    
######################################################################################
if __name__ == '__main__':
    main()
