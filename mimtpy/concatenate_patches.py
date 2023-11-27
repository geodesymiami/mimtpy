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
"""

EXAMPLE = """example:

    concatenate_patches.py --project BJSenAT142 --subproject Mentougou --network network_delaunay_4 --datatype velocity --PSDS --mask
  
    concatenate_patches.py --project BJSenAT142 --subproject Mentougou --network network_delaunay_4 --datatype mask --PSDS --mask
  
    concatenate_patches.py --project TangshanSenAT69 --network network_delaunay_4 --datatype timeseries --PS 
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Concatenate patches of miaplpy under radar coordinate',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=NOTE+'\n'+EXAMPLE)

    parser.add_argument('--project', dest='project',nargs=1, help='Project to be processed.\n')
    
    parser.add_argument('--subproject', dest='subproject',nargs='?', help='Sub project to be processed.\n')
    
    parser.add_argument('--network', dest='network',nargs=1, help='Network name of miaplpy.\n')
    
    parser.add_argument('--datatype', dest='datatype',nargs=1, help='datatype to be processed.\n')

    parser.add_argument('--sdate', nargs='?', help='start date.\n')
    
    parser.add_argument('--edate', nargs='?', help='end date.\n')
    
    parser.add_argument('--dryrun', action='store_true', default=False, help='used for timeseries. Only show the removed date for each patch. Don not concatenate patches\n')    

    parser.add_argument('--PSDS', action='store_true', default=False, help='whether use S1*_*PSDS.he5 data\n')    
    
    parser.add_argument('--PS', action='store_true', default=False, help='whether use S1*_*PS.he5 data\n')    
    
    parser.add_argument('--mask', action='store_true', default=False, help='whether mask data\n')    
    
    parser.add_argument('--nogeo', action='store_true', default=False, help='if there is already rdr_geometry file. This option could be used\n')    
    
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

def concatenate_patches(project, pthList, network, datatype, inps):
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

    # start date and end date
    if inps.sdate is not None:
        sdate = inps.sdate[0]
    if inps.edate is not None:
        edate = inps.edate[0]
    
    for Pnum, PID in enumerate(patches_dir):
        if inps.subproject is not None:
            patch_ID = PID.split('/')[-1].split('_')[2]
        else:
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
            #print('Read maskPS.h5 file')
            #maskPS = readfile.read('maskPS.h5')[0]

            print('************************************************************************************************')
            network_dir = os.path.abspath(os.path.join(patch_dir + '/' + network + '/'))
            os.chdir(network_dir)
            print('Go to network dir:', network_dir)

            if not os.path.isdir(network_dir + '/inputs/'):
                os.makedirs(network_dir + '/inputs/')
            
            # read S1**.he5 file
            if inps.PS:
                HDFEOS_file = mu.find_PS_HDFEOS_fullname(network_dir)
            elif inps.PSDS:
                HDFEOS_file = mu.find_PSDS_HDFEOS_fullname(network_dir)
            print('The HDFEOS file is {}'.format(HDFEOS_file))

            # read mask
            maskfile, msk_atr = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/quality/mask')

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
                outfile = network_dir + '/inputs/' + 'geometryRadar.h5' 
                writefile.write(dsDict, out_file=outfile, metadata=atr)
            
            if datatype == 'velocity':
                if inps.sdate is not None:
                    scp_args = ['timeseries2velocity.py ' + network_dir + '/' + HDFEOS_file + ' --start-date ' + sdate + ' --end-date ' + edate]
                else:
                    scp_args = ['timeseries2velocity.py ' + network_dir + '/' + HDFEOS_file]
                scp_args = mu.separate_string_by_space(scp_args)
                print(scp_args)
                # run timeseries2velocity.py
                os.system(scp_args)
                #shutil.move((network_dir + '/velocity.h5'), (network_dir + '/rdr_vel.h5'))
                vel_data, vel_atr = readfile.read('velocity.h5')
                if inps.mask:
                    vel_data[maskfile == 0] = np.nan
                #writefile.write(vel_data, out_file='rdr_vel_msk.h5', metadata=vel_atr)
                writefile.write(vel_data, out_file='velocity.h5', metadata=vel_atr)

            elif datatype == 'mask':
                msk_atr['FILE_TYPE'] = 'mask'
                mskDict = dict()
                mskDict['mask'] = maskfile
                writefile.write(mskDict, out_file='mask.h5', metadata=msk_atr) 

            elif datatype == 'timeseries':
                #displacement = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/observation/displacement')[0]
                bperp_date = h5py.File(HDFEOS_file,'r')
                data_date = np.array(bperp_date['/HDFEOS/GRIDS/timeseries/observation/date']).tolist()
                date_patches.append(data_date)

    if not inps.nopreprocess:
        if datatype == 'timeseries':
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

        if datatype == 'timeseries':
            for patch_dir in patches_dir:
                print('************************************************************************************************')
                os.chdir(patch_dir)
                network_dir = os.path.abspath(os.path.join(patch_dir + '/' + network + '/'))
                os.chdir(network_dir)

                # read S1**.he5 file
                if inps.PS:
                    HDFEOS_file = mu.find_PS_HDFEOS_fullname(network_dir)
                elif inps.PSDS:
                    HDFEOS_file = mu.find_PSDS_HDFEOS_fullname(network_dir)
                print('The HDFEOS file is {}'.format(HDFEOS_file))
                
                maskfile = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/quality/mask')[0]
                # process the date
                displacement, ts_atr = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/observation/displacement')
                bperp_date = h5py.File(HDFEOS_file,'r')
                data_bperp = bperp_date['/HDFEOS/GRIDS/timeseries/observation/bperp']
                data_date = np.array(bperp_date['/HDFEOS/GRIDS/timeseries/observation/date']).tolist()
                #data_date = bperp_date[(dataset+'date')]

                # output info    
                #ts_data_name = 'rdr_ts_msk.h5'
                ts_data_name = 'timeseries.h5'
                outfile = network_dir + '/' + ts_data_name

                # getting flag
                date_position = []
                for date in date_intersection:
                    date_position.append(data_date.index(date))

                # get the time series displacement
                data_ts = displacement[date_position,:,:]
                data_bp = np.array(data_bperp[date_position],dtype=np.float64)
                date_intersection = np.array(date_intersection, dtype=np.string_)
                # calculating timeseries referencing the first date 
                data_ts = data_ts - data_ts[0,:,:]

                # mask data
                if inps.mask:
                    maskfile = maskfile * 1.0
                    maskfile[np.where(maskfile == 0.0)] = np.nan
                    data_ts = [maskfile] * data_ts.shape[0] * data_ts

                # write to HDF5 file
                dsDict = dict()
                dsDict['bperp'] = data_bp
                dsDict['date'] = date_intersection
                dsDict['timeseries'] = data_ts
                ts_atr['FILE_TYPE'] = 'timeseries'
                ts_atr['REF_DATE'] = str(date_intersection.astype(np.int64)[0]) 
                writefile.write(dsDict, out_file=outfile, metadata=ts_atr)
    
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
        #scp_args = [patches_dir[0] + '/' + network + '/' + datatype + '.h5', patches_dir[1] + '/' + network + '/' + datatype + '.h5', '--geo-file1', patches_dir[0] + '/' + network + '/geometryRadar.h5', '--geo-file2', patches_dir[1] + '/' + network + '/geometryRadar.h5', '--geo-full', geometry_Big, '--geo_write', '--out-suffix', temp_name,  '--outdir', output_dir + '/']
        scp_args = [patches_dir[0] + '/' + network + '/' + datatype + '.h5', patches_dir[1] + '/' + network + '/' + datatype + '.h5', '--geo-full', geometry_Big, '--geo_write', '--out-suffix', temp_name,  '--outdir', output_dir + '/']
        scp_args = mu.separate_string_by_space(scp_args)
        # run concatenate_radarGeo.py
        print('concatenate_radarGeo.py',scp_args)
        mimtpy.concatenate_radarGeo.main(scp_args.split())
    else:
        for i in np.arange(pro_num - 1):
            if i == 0:
                temp_name = patch_range[0] + '_' + patch_range[1]
                temp_geo = 'geometryRadar_' + patch_range[0] + '_' + patch_range[1]
                #scp_args = [patches_dir[0] + '/' + network + '/' + datatype + '.h5', patches_dir[1] + '/' + network + '/' + datatype + '.h5', '--geo-file1', patches_dir[0] + '/' + network + '/geometryRadar.h5', '--geo-file2', patches_dir[1] + '/' + network + '/geometryRadar.h5', '--geo-full', geometry_Big, '--geo_write', '--out-suffix', temp_name,  '--outdir', output_dir + '/']
                scp_args = [patches_dir[0] + '/' + network + '/' + datatype + '.h5', patches_dir[1] + '/' + network + '/' + datatype + '.h5', '--geo-full', geometry_Big, '--geo-write', '--out-suffix', temp_name,  '--outdir', output_dir + '/']
            elif i > 0 and i < (pro_num - 1 - 1):
                temp_name = temp_files[i-1] + '_' + patch_range[i+1]
                temp_geo = temp_geos[i-1] + '_' + patch_range[i+1]
                scp_args = [output_dir + '/' + datatype + '_' + temp_files[i-1] +'.h5', patches_dir[i+1] + '/' + network + '/' + datatype + '.h5', '--geo-file1', output_dir + '/' + temp_geos[i-1] + '.h5', '--geo-file2', patches_dir[i+1] + '/' + network + '/inputs/geometryRadar.h5', '--geo-full', geometry_Big, '--geo-write', '--out-suffix', temp_name, '--outdir', output_dir + '/']
            else:
                #temp_name = 'velocity_lat_' + str(lat_range[0]) + '_' + str(lat_range[-1])
                temp_name = 'full'
                temp_geo = 'geometryRadar_full'
                scp_args = [output_dir + '/' + datatype + '_' + temp_files[i-1] +'.h5', patches_dir[i+1]+ '/' + network + '/' + datatype + '.h5', '--geo-file1', output_dir + '/' + temp_geos[i-1] + '.h5', '--geo-file2', patches_dir[i+1] + '/' + network + '/inputs/geometryRadar.h5', '--geo-full', geometry_Big, '--geo-write', '--out-suffix', temp_name, '--outdir', output_dir + '/']

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
    datatype = inps.datatype[0]

    concatenate_patches(project, pthList, network, datatype, inps)
     
######################################################################################
if __name__ == '__main__':
    main()
