#!/usr/bin/env python3
#################################################################
# Program is used for concatenate patches of miaplpy            #
# Author: Lv Xiaoran                                            #
# Created: Jan 2021                                             #
#################################################################

import os
import argparse
import numpy as np
import shutil
import copy
import h5py
import sys
import pandas as pd

import mintpy
import mintpy.workflow
from mintpy.utils import readfile, writefile, ptime, utils as ut
from mintpy import view
from mintpy import timeseries2velocity as ts2vel

import mimtpy
import mimtpy.workflow
from mimtpy.utils import multitrack_utilities as mu
######################################################################################
NOTE = """
Please Note
1. this script concatenate all patches together
2. this script calls the concatenate_rdrGeo.py to do the concatenation of two datasets. 
"""

EXAMPLE = """example:

    concatenate_patches.py --project BJSenAT142 --subproject Mentougou --network network_delaunay_4 --DS
  
    concatenate_patches.py --project TangshanSenAT69 --network network_delaunay_4 --PS 
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Concatenate patches of miaplpy under radar coordinate',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=NOTE+'\n'+EXAMPLE)

    parser.add_argument('--project', dest='project',nargs=1, help='Project to be processed.\n')
    
    parser.add_argument('--subproject', dest='subproject',nargs='?', help='Sub project to be processed.\n')
    
    parser.add_argument('--network', dest='network',nargs=1, help='Network name of miaplpy.\n')
    
    parser.add_argument('--DS', action='store_true', default=False, help='whether use S1*_*PSDS.he5 data\n')    
    
    parser.add_argument('--PS', action='store_true', default=False, help='whether use S1*_*PS.he5 data\n')    
    
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps

def read_patches(patches_file):
    patches_name = []
    file = open(patches_file, 'r')
    temps = file.readlines()
    file.close()
    for line in temps:     
        line = line.strip('\n')
        patches_name.append(line)
    
    return patches_name

def find_S1(inps, network_dir):
    # read S1**.he5 file
    if inps.PS:
        HDFEOS_file = mu.find_PS_HDFEOS_fullname(network_dir)
    elif inps.DS:
        HDFEOS_file = mu.find_DS_HDFEOS_fullname(network_dir)
    print('The HDFEOS file is {}'.format(HDFEOS_file))
    return HDFEOS_file

#def concatenate_patches(project, pthList, network, datatype, inps):
def concatenate_patches(project, pthList, network, inps):
    project_dir = os.getenv('SCRATCHDIR') + '/' + project + '/'
    print('change directory to %s' % project_dir)
    os.chdir(project_dir)

    # find the mintpy/inputs/goeometryRadar.h5 processed with 1:1 looks
    geometry_Big = project_dir +  '/mintpy_11/inputs/geometryRadar.h5'
    print('Try to find the geometry data processed with 1:1 looks')
    if os.path.isfile(geometry_Big):
        print('read the geometry data processed with 1:1 looks')
    else:
        raise ValueError('No geometryRadar.h5 processed with 1:1 looks')
 
    # patches directory
    patches_dir = pthList

    # preprocess the S1**_Del4PSDS.he5 and geometryRadar.h5 file
    patches_ID = []
    date_patches = []

    for Pnum, PID in enumerate(patches_dir):
        print('Check whether the patch folder exists.')
        if inps.subproject is not None:
            patch_ID = PID.split('/')[-1].split('_')[2]
        else:
            patch_ID = PID.split('/')[-1].split('_')[1]
        patches_ID.append(patch_ID)

        patch_dir = patches_dir[Pnum]
        if not os.path.isdir(patch_dir):
            raise ValueError('%s directory does not exist. Please check!' % s)
            os.chdir(patch_dir)

    # create the output dir
    if inps.subproject is not None:
        output_dir = os.path.abspath(os.path.join(project_dir + '/' + 'miaplpy_' + inps.subproject + '_Big/'))
    else:
        output_dir = os.path.abspath(os.path.join(project_dir + '/' + 'miaplpyBig/' + '/'))
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    print('************************************************************************************************')
    print('the output dir for concatenation is {}.'.format(output_dir))

    os.chdir(output_dir)
    
    # start the loop
    temp_files = []
    temp_geos = []
    pro_num = len(patches_dir)
    #lat_range = np.arange(lat_min, lat_max)
    patch_range = patches_ID
    if pro_num == 2:
        temp_name = patch_range[0] + '_' + patch_range[1]
        #network_dir = os.path.abspath(os.path.join(patch_dir[0] + '/' + network + '/'))
        S1file_1 = find_S1(inps, os.path.abspath(os.path.join(patches_dir[0] + '/' + network + '/')))
        S1file_2 = find_S1(inps, os.path.abspath(os.path.join(patches_dir[1] + '/' + network + '/')))
        scp_args = [os.path.abspath(os.path.join(patches_dir[0] + '/' + network + '/')) + '/' + S1file_1, os.path.abspath(os.path.join(patches_dir[1] + '/' + network + '/')) + '/' + S1file_2, '--geo-full', geometry_Big, '--out-suffix', temp_name,  '--outdir', output_dir + '/']
        scp_args = mu.separate_string_by_space(scp_args)
        # run concatenate_radarGeo.py
        print('concatenate_radarGeo.py',scp_args)
        mimtpy.concatenate_radarGeo.main(scp_args.split())
    else:
        for i in np.arange(pro_num - 1):
            if i == 0:
                temp_name = patch_range[0] + '_' + patch_range[1]
                S1file_1 = find_S1(inps, os.path.abspath(os.path.join(patches_dir[0] + '/' + network + '/')))
                S1file_2 = find_S1(inps, os.path.abspath(os.path.join(patches_dir[1] + '/' + network + '/')))
                temp_geo = 'geometryRadar_' + patch_range[0] + '_' + patch_range[1]
                scp_args = [os.path.abspath(os.path.join(patches_dir[0] + '/' + network + '/')) + '/' + S1file_1, os.path.abspath(os.path.join(patches_dir[1] + '/' + network + '/')) + '/' + S1file_2, '--geo-full', geometry_Big, '--out-suffix', temp_name,  '--outdir', output_dir + '/']
            elif i > 0 and i < (pro_num - 1 - 1):
                S1file_2 = find_S1(inps, os.path.abspath(os.path.join(patches_dir[i+1] + '/' + network + '/')))
                temp_name = '_'.join(temp_files[i-1].split('_')[1:]) + '_' + patch_range[i+1]
                temp_geo = temp_geos[i-1] + '_' + patch_range[i+1]
                scp_args = [output_dir + '/' + temp_files[i-1] +'.he5', os.path.abspath(os.path.join(patches_dir[i+1] + '/' + network + '/')) + '/' + S1file_2, '--geo-full', geometry_Big, '--out-suffix', temp_name, '--outdir', output_dir + '/']
            else:
                S1file_2 = find_S1(inps, os.path.abspath(os.path.join(patches_dir[i+1] + '/' + network + '/')))
                temp_name = 'full'
                temp_geo = 'geometryRadar_full'
                scp_args = [output_dir + '/' + temp_files[i-1] +'.he5', os.path.abspath(os.path.join(patches_dir[i+1] + '/' + network + '/')) + '/' + S1file_2, '--geo-full', geometry_Big, '--out-suffix', temp_name, '--outdir', output_dir + '/']

            temp_files.append('S1_' + temp_name)
            temp_geos.append(temp_geo)

            scp_args = mu.separate_string_by_space(scp_args)
            # run concatenate_radarGeo.py
            print('concatenate_radarGeo.py',scp_args)
            mimtpy.concatenate_radarGeo.main(scp_args.split())   
    
def getPatches(project, inps):
    """get the patches in the project dir and sorted them by lat/lon range"""
    pro_dir = os.path.realpath(project) if project[0] != '/' else project
    assert(os.path.isdir(pro_dir)), f"error project path `{pro_dir}`"

    # find patch dirs
    patch_dirs = []
    for root, dirs, files in os.walk(pro_dir):
        for dir in dirs:
            if inps.subproject is not None:
                if dir.find('miaplpy_' + inps.subproject) != -1 and dir.find('Big') == -1:
                    patch_dirs.append(os.path.join(root, dir))
            else:
                if dir.find('miaplpy_') != -1 and dir.find('Big') == -1:
                    patch_dirs.append(os.path.join(root, dir))

    
    patch_num = len(patch_dirs)
    patch_obj = []
    patch_lon = []
    patch_lat = []

    for i in np.arange(patch_num):
        patch = patch_dirs[i]
        HDFEOS_file = mu.find_HDFEOS_fullname(patch + '/' + inps.network[0] + '/')
        patch_atr = readfile.read_attribute(patch + '/' + inps.network[0] + '/' + HDFEOS_file)
        maxlo = float(patch_atr['data_footprint'].split(',')[1].split(' ')[0])
        maxla = float(patch_atr['data_footprint'].split(',')[1].split(' ')[1])
        minlo = float(patch_atr['data_footprint'].split(',')[3].split(' ')[0])
        minla = float(patch_atr['data_footprint'].split(',')[3].split(' ')[1])

        patch_obj.append(patch)
        patch_lon.append(minlo)
        patch_lat.append(minla)

    patch_dict = {'patch': patch_obj, 'min_lon':patch_lon, 'min_lat':patch_lat}
    patch_pd = pd.DataFrame(patch_dict)
    patch_sorted = patch_pd.sort_values(by=['min_lon','min_lat'])
    patch_dir_sorted = patch_sorted['patch'].values.tolist()

    return patch_dir_sorted 

######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)   
    
    project = inps.project[0]

    pthList = getPatches(project, inps)

    print('The patches to be concatenated are:')
    print(pthList)

    network = inps.network[0]
    #datatype = inps.datatype[0]

    #concatenate_patches(project, pthList, network, datatype, inps)
    concatenate_patches(project, pthList, network, inps)
     
######################################################################################
if __name__ == '__main__':
    main()
