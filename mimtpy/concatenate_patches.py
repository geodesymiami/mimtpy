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
2. this script first preprocesses the S1**.he5 file to obtain the velocity/timeseries files. Finally it calls the concatenate_rdrGeo.py to do the concatenation of two datasets. 
3. in script, the maskPS means using maskPS.h5. the maskTC means using maskTempCoh.h5 
"""

EXAMPLE = """example:

    concatenate_patches.py --project TangshanSenAT69 --network network_delaunay_4 --datatype vel --maskTC
  
    concatenate_patches.py --project TangshanSenAT69 --network network_delaunay_4 --datatype ts --maskPS 
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Concatenate patches of miaplpy under radar coordinate',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=NOTE+'\n'+EXAMPLE)

    parser.add_argument('--project', dest='project',nargs=1, help='Project to be processed.\n')
    
    parser.add_argument('--network', dest='network',nargs=1, help='Network name of miaplpy.\n')
    
    parser.add_argument('--datatype', dest='datatype',nargs=1, help='datatype to be processed.\n')

    parser.add_argument('--sdate', nargs='?', help='start date.\n')
    
    parser.add_argument('--edate', nargs='?', help='end date.\n')
    
    parser.add_argument('--dryrun', action='store_true', default=False, help='used for timeseries. Only show the removed date for each patch. Don not concatenate patches\n')    

    parser.add_argument('--maskTC', action='store_true', default=False, help='whether use maskTempCoh to mask data\n')    
    
    parser.add_argument('--maskPS', action='store_true', default=False, help='whether use maskPS to mask data\n')    
    
    parser.add_argument('--nogeo', action='store_true', default=False, help='if there is already rdr_geometry_msk file. This option could be used\n')    
    
    parser.add_argument('--nopreprocess', action='store_true', default=False, help='if there is already data in each path. This option could be used\n')    
    
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

def read_ref_point(ref_file):
    data = np.loadtxt(ref_file, dtype=float, delimiter=' ')
    
    return data

def concatenate_patches(project, pthList, network, datatype, inps):
    project_dir = os.getenv('SCRATCHDIR') + '/' + project + '/'
    print('change directory to %s' % project_dir)
    os.chdir(project_dir)
 
    # patches directory
    patches_dir = pthList

    # preprocess the S1**_Del4PSDS.he5 and geometryRadar.h5 file
    patches_ID = []
    date_patches = []

    # start date and end date
    if inps.sdate is not None:
        sdate = inps.sdate[0]
    if inps.edate is not None:
        edate = inps.edate[0]
    
    for Pnum, PID in enumerate(patches_dir):
        patch_ID = PID.split('/')[-1].split('_')[1]
        patches_ID.append(patch_ID)

        if not inps.nopreprocess:
            patch_dir = patches_dir[Pnum]
            if not os.path.isdir(patch_dir):
                raise ValueError('%s directory does not exist. Please check!' % s)
            print('************************************************************************************************')
            os.chdir(patch_dir)
            print('Go to patch dir:', patch_dir)

            # read maskPS.h5
            print('Read maskPS.h5 file')
            maskPS = readfile.read('maskPS.h5')[0]

            print('************************************************************************************************')
            network_dir = os.path.abspath(os.path.join(patch_dir + '/' + network + '/'))
            os.chdir(network_dir)
            print('Go to network dir:', network_dir)

            # read S1**_Del4PSDS.he5 file
            HDFEOS_file = mu.find_Del4PSDS_fullname(network_dir)

            # read mask
            maskfile = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/quality/mask')[0]

            # read and write geometryRadar file
            if not inps.nogeo:
                azi, atr = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/geometry/azimuthAngle')
                inc = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/geometry/incidenceAngle')[0]
                hgt = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/geometry/height')[0]
                lat = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/geometry/latitude')[0]
                lon = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/geometry/longitude')[0]
                sdM = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/geometry/shadowMask')[0]
                srd = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/geometry/slantRangeDistance')[0]

                # write to HDF5 file
                dsDict = dict()
                dsDict['azimuthAngle'] = azi
                dsDict['incidenceAngle'] = inc
                dsDict['height'] = hgt
                dsDict['latitude'] = lat
                dsDict['longitude'] = lon
                dsDict['shadowMask'] = sdM
                dsDict['slantRangeDistance'] = srd 
                atr['FILE_TYPE'] = 'geometry'
                outfile = network_dir + '/' + 'rdr_geometry_msk.h5' 
                writefile.write(dsDict, out_file=outfile, metadata=atr)
            
            if datatype == 'vel':
                if inps.sdate is not None:
                    scp_args = ['timeseries2velocity.py ' + network_dir + '/' + HDFEOS_file + ' --start-date ' + sdate + ' --end-date ' + edate]
                else:
                    scp_args = ['timeseries2velocity.py ' + network_dir + '/' + HDFEOS_file]
                scp_args = mu.separate_string_by_space(scp_args)
                print(scp_args)
                # run timeseries2velocity.py
                os.system(scp_args)
                shutil.move((network_dir + '/velocity.h5'), (network_dir + '/rdr_vel.h5'))
                vel_data, vel_atr = readfile.read('rdr_vel.h5')
                if inps.maskTC:
                    vel_data[maskfile == 0] = np.nan
                if inps.maskPS:
                    vel_data[maskPS == 0] = np.nan
                writefile.write(vel_data, out_file='rdr_vel_msk.h5', metadata=vel_atr)

            elif datatype == 'ts':
                #displacement = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/observation/displacement')[0]
                bperp_date = h5py.File(HDFEOS_file,'r')
                data_date = np.array(bperp_date['/HDFEOS/GRIDS/timeseries/observation/date']).tolist()
                date_patches.append(data_date)

    if not inps.nopreprocess:
        if datatype == 'ts':
            date_intersection = list(set(date_patches[0]).intersection(*date_patches[1:]))
            date_intersection.sort()
            if inps.dryrun:
                date_shared = copy.deepcopy(date_intersection)
                for date_patch, patch_ID in zip(date_patches, patches_ID):
                    print('The removed dates for patch %s is :' % patch_ID)
                    date_diff = list(set(date_patch).difference(set(date_shared)))
                    date_diff.sort()
                    #date_diff = list(set(date_patch)^set(date_shared)).sort()
                    if len(date_diff) == 0:
                        print('NONE')
                    else:
                        date_diff = np.array(date_diff)
                        date_diff_final = [element.decode('utf-8') for element in date_diff]
                        print(date_diff_final)
                sys.exit(0)

        if datatype == 'ts':
            for patch_dir in patches_dir:
                print('************************************************************************************************')
                os.chdir(patch_dir)
                network_dir = os.path.abspath(os.path.join(patch_dir + '/' + network + '/'))
                os.chdir(network_dir)

                # find HDFEOS file
                HDFEOS_file = mu.find_Del4PSDS_fullname(network_dir)
                print('The HDFEOS file is {}'.format(HDFEOS_file))
                maskfile = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/quality/mask')[0]
                # process the date
                displacement, atr = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/observation/displacement')
                bperp_date = h5py.File(HDFEOS_file,'r')
                data_bperp = bperp_date['/HDFEOS/GRIDS/timeseries/observation/bperp']
                data_date = np.array(bperp_date['/HDFEOS/GRIDS/timeseries/observation/date']).tolist()
                #data_date = bperp_date[(dataset+'date')]

                # output info    
                ts_data_name = 'rdr_ts_msk.h5'
                outfile = network_dir + '/' + ts_data_name

                # getting flag
                date_position = []
                for date in date_intersection:
                    date_position.append(data_date.index(date))

                # get the time series displacement
                data_ts = displacement[date_position,:,:]
                data_bp = np.array(data_bperp[date_position],dtype=np.float64)
                date_intersection = np.array(date_intersection, dtype=np.string_)

                # mask data
                if inps.maskTC:
                    data_ts[:, maskfile == 0] = np.nan
                if inps.maskPS:
                    data_ts[:, maskPS == 0] = np.nan

                # write to HDF5 file
                dsDict = dict()
                dsDict['bperp'] = data_bp
                dsDict['date'] = date_intersection
                dsDict['timeseries'] = data_ts
                atr['FILE_TYPE'] = 'timeseries'
                atr['REF_DATE'] = str(date_intersection.astype(np.int)[0]) 
                writefile.write(dsDict, out_file=outfile, metadata=atr)
    
    # create the output dir
    output_dir = os.path.abspath(os.path.join(project_dir + '/' + 'miaplpyBig' + '/rdr/'))
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
        temp_name = 'rdr_' + datatype
        scp_args = [patches_dir[0] + '/' + network + '/rdr_' + datatype + '_msk.h5', patches_dirs[1] + '/' + network + '/rdr_' + datatype + '_msk.h5', '--geo_file1', patches_dir[0] + '/' + network + '/rdr_geometry_msk.h5', '--geo_file2', patches_dir[1] + '/' + network + '/rdr_geometry_msk.h5', '--datatype', datatype, '--geo_write', '--output', temp_name,  '--outdir', outputdir + '/']
        scp_args = mu.separate_string_by_space(scp_args)
        # run concatenate_offset.py
        print('concatenate_offset.py',scp_args)
        mimtpy.concatenate_offset.main(scp_args.split())   
    else:
        for i in np.arange(pro_num - 1):
            if i == 0:
                temp_name = patch_range[0] + '_' + patch_range[1]
                temp_geo = 'geometry_' + patch_range[0] + '_' + patch_range[1]
                scp_args = [patches_dir[0] + '/' + network + '/rdr_' + datatype + '_msk.h5', patches_dir[1] + '/' + network + '/rdr_' + datatype + '_msk.h5', '--geo_file1', patches_dir[0] + '/' + network + '/rdr_geometry_msk.h5', '--geo_file2', patches_dir[1] + '/' + network + '/rdr_geometry_msk.h5', '--datatype', datatype, '--geo_write', '--output', temp_name,  '--outdir', output_dir + '/']
            elif i > 0 and i < (pro_num - 1 - 1):
                temp_name = temp_files[i-1] + '_' + patch_range[i+1]
                temp_geo = temp_geos[i-1] + '_' + patch_range[i+1]
                scp_args = [output_dir + '/' + datatype + '_' + temp_files[i-1] +'.h5', patches_dir[i+1] + '/' + network + '/rdr_' + datatype + '_msk.h5', '--geo_file1', output_dir + '/' + temp_geos[i-1] + '.h5', '--geo_file2', patches_dir[i+1] + '/' + network + '/rdr_geometry_msk.h5', '--datatype', datatype, '--geo_write', '--output', temp_name, '--outdir', output_dir + '/']
            else:
                #temp_name = 'velocity_lat_' + str(lat_range[0]) + '_' + str(lat_range[-1])
                temp_name = 'full'
                temp_geo = 'geometry_full'
                scp_args = [output_dir + '/' + datatype + '_' + temp_files[i-1] +'.h5', patches_dir[i+1]+ '/' + network + '/rdr_' + datatype + '_msk.h5', '--geo_file1', output_dir + '/' + temp_geos[i-1] + '.h5', '--geo_file2', patches_dir[i+1] + '/' + network + '/rdr_geometry_msk.h5', '--datatype', datatype, '--geo_write', '--output', temp_name, '--outdir', output_dir + '/']

            temp_files.append(temp_name)
            temp_geos.append(temp_geo)

            scp_args = mu.separate_string_by_space(scp_args)
            # run concatenate_radarGeo.py
            print('concatenate_radarGeo.py',scp_args)
            mimtpy.concatenate_radarGeo.main(scp_args.split())   
    
def getPatches(project, inps):
    """get the patches in the project dir and sorted them by lat/lon range"""
    pro_dir = os.path.realpath(project) if project[0] != '/' else project
    assert(os.path.isdir(pro_dir)), f"error project path `{pro_dir}`"

    # find all patch dirs
    patch_dirs = []
    for root, dirs, files in os.walk(pro_dir):
        for dir in dirs:
            if dir.find('miaplpy_') != -1:
                patch_dirs.append(os.path.join(root, dir))

    
    patch_num = len(patch_dirs)
    patch_obj = []
    patch_lon = []
    patch_lat = []

    for i in np.arange(patch_num):
        patch = patch_dirs[i]
        HDFEOS_file = mu.find_Del4PSDS_fullname(patch + '/' + inps.network[0] + '/')
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
    datatype = inps.datatype[0]

    concatenate_patches(project, pthList, network, datatype, inps)
     
######################################################################################
if __name__ == '__main__':
    main()
