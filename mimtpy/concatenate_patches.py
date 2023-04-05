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
1. this script concatenate all patches of one project together
2. Two files need to be prepared before do concatenation. 
   2.1 concatenate_patch.txt contains the order of patches to be concatenated. Taken TangshanSenAT69 project as example, the content of concatenate_patch.txt is miaplpy_NE_201410_202212 miaplpy_NNE_201410_202212 miaplpy_SE_201410_202212 miaplpy_NW_201410_202212 miaplpy_NNW_201410_202212 miaplpy_SW_201410_202212.
   2.2 reference_point.txt contains the latitude and longitude of the point located at the overlapping region between different patches. The order of these points follows the patches order in the concatenate_patches. Now users need to specify these points. 
3. this script first preprocesses the S1**_Del4PSDS.he5 file to obtain the velocity/timeseries files. Then it calls the concatenate_rdrGeo.py to do the concatenation of two datasets. 
"""

EXAMPLE = """example:

    concatenate_patches.py concatenate_patch.txt --project TangshanSenAT69 --network network_delaunay_4 --datatype velocity --ref_file reference_point.txt 
  
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Concatenate patches of miaplpy under radar coordinate',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=NOTE+'\n'+EXAMPLE)

    parser.add_argument('patches_file', nargs=1, type=str, help='concatenated patches sorted by user. \n')

    parser.add_argument('--project', dest='project',nargs=1, help='Project to be processed.\n')
    
    parser.add_argument('--network', dest='network',nargs=1, help='Network name of miaplpy.\n')
    
    parser.add_argument('--datatype', dest='datatype',nargs=1, help='datatype to be processed.\n')

    parser.add_argument('--start_date', dest='sdate',nargs='?', help='date to start calculate the velocity.\n')
    
    parser.add_argument('--end_date', dest='edate',nargs='?', help='date to finish calculate the velocity.\n')
    
    #parser.add_argument('--velocity', action='store_true', default=False, 
    #                    help='whether calculate velocity data.\n'+
    #                    'The cumulative displacement between start date and end date will be calculated as default.\n')

    parser.add_argument('--dryrun', action='store_true', default=False, help='used for timeseries. Only show the removed date for each patch. Didnot concatenate patches\n')    

    parser.add_argument('--nogeo', action='store_true', default=False, help='if there is already rdr_geometry_msk file. This option could be used\n')    
    
    parser.add_argument('--nopreprocess', action='store_true', default=False, help='if there is already data in each path. This option could be used\n')    
    
    parser.add_argument('--ref_file',dest='ref',nargs=1, 
                        help='reference point file.\n'+
                        'Store all the reference points corresponding to the patches_file.\n')
    
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

def concatenate_patches(project, patches_name, network, datatype, refpoi_list, inps):
    project_dir = os.getenv('SCRATCHDIR') + '/' + project + '/'
    print('change directory to %s' % project_dir) 
    os.chdir(project_dir)
 
    # patches directory
    patches_dir = []
    for patch_name in patches_name:
        patch_dir = os.path.abspath(os.path.join(project_dir + '/' + patch_name + '/'))
        patches_dir.append(patch_dir)

    # preprocess the S1**_Del4PSDS.he5 and geometryRadar.h5 file
    patches_ID = []
    date_patches = []
    
    sdate = inps.sdate
    edate = inps.edate

 
    for Pnum, PID in enumerate(patches_name):
        patch_ID = PID.split('_')[1]
        patches_ID.append(patch_ID)

        if not inps.nopreprocess:
            patch_dir = patches_dir[Pnum]
            if not os.path.isdir(patch_dir):
                raise ValueError('%s directory does not exist. Please check!' % s)
            print('************************************************************************************************')
            os.chdir(patch_dir)
            print('Go to patch dir:', patch_dir)
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

                # mask the data
                #azi[maskfile == 0] = np.nan
                #inc[maskfile == 0] = np.nan
                #hgt[maskfile == 0] = np.nan
                #lat[maskfile == 0] = np.nan
                #lon[maskfile == 0] = np.nan
                #sdM[maskfile == 0] = np.nan
                #srd[maskfile == 0] = np.nan
            
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
            
            if datatype == 'velocity':
                if inps.sdate is not None:
                    scp_args = ['timeseries2velocity.py ' + network_dir + '/' + HDFEOS_file + ' --start-date ' + sdate + ' --end-date ' + edate]
                else:
                    scp_args = ['timeseries2velocity.py ' + network_dir + '/' + HDFEOS_file]
                scp_args = mu.separate_string_by_space(scp_args)
                print(scp_args)
                # run timeseries2velocity.py
                os.system(scp_args)
                shutil.move((network_dir + '/velocity.h5'), (network_dir + '/rdr_velocity.h5'))
                vel_data, vel_atr = readfile.read('rdr_velocity.h5')
                vel_data[maskfile == 0] = np.nan
                writefile.write(vel_data, out_file='rdr_velocity_msk.h5', metadata=vel_atr)

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
                ts_data_name = 'rdr_timeseries_msk.h5'
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
                data_ts[:, maskfile == 0] = np.nan
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
        scp_args = [patches_dir[0] + '/' + network + '/rdr_' + datatype + '_msk.h5', patches_dirs[1] + '/' + network + '/rdr_' + datatype + '_msk.h5', '--geo_file1', patches_dir[0] + '/' + network + '/rdr_geometry_msk.h5', '--geo_file2', patches_dir[1] + '/' + network + '/rdr_geometry_msk.h5', '--ref', str(refpoi_list[0, 0]), str(refpoi_list[0, 1]), '--datatype', datatype, '--output', temp_name,  '--outdir', outputdir + '/']
        scp_args = mu.separate_string_by_space(scp_args)
        # run concatenate_offset.py
        print('concatenate_offset.py',scp_args)
        mimtpy.concatenate_offset.main(scp_args.split())   
    else:
        for i in np.arange(pro_num - 1):
            temp_refpoi = refpoi_list[i, :]
            if i == 0:
                temp_name = 'rdr_' + datatype + '_' + patch_range[0] + '_' + patch_range[1]
                temp_geo = 'rdr_geometry_' + patch_range[0] + '_' + patch_range[1]
                scp_args = [patches_dir[0] + '/' + network + '/rdr_' + datatype + '_msk.h5', patches_dir[1] + '/' + network + '/rdr_' + datatype + '_msk.h5', '--geo_file1', patches_dir[0] + '/' + network + '/rdr_geometry_msk.h5', '--geo_file2', patches_dir[1] + '/' + network + '/rdr_geometry_msk.h5', '--ref', str(temp_refpoi[0]), str(temp_refpoi[1]), '--datatype', datatype, '--output', temp_name,  '--outdir', output_dir + '/']
                
            elif i > 0 and i < (pro_num - 1 - 1):
                temp_name = temp_files[i-1] + '_' + patch_range[i+1]
                temp_geo = temp_geos[i-1] + '_' + patch_range[i+1]
                scp_args = [output_dir + '/' + temp_files[i-1] +'.h5', patches_dir[i+1] + '/' + network + '/rdr_' + datatype + '_msk.h5', '--geo_file1', output_dir + '/' + temp_geos[i-1] + '.h5', '--geo_file2', patches_dir[i+1] + '/' + network + '/rdr_geometry_msk.h5', '--ref', str(temp_refpoi[0]), str(temp_refpoi[1]), '--datatype', datatype, '--output', temp_name, '--outdir', output_dir + '/']
            else:
                #temp_name = 'velocity_lat_' + str(lat_range[0]) + '_' + str(lat_range[-1])
                temp_name = 'rdr_' + datatype + '_full'
                temp_geo = 'rdr_geometry_full'
                scp_args = [output_dir + '/' + temp_files[i-1] +'.h5', patches_dir[i+1]+ '/' + network + '/rdr_' + datatype + '_msk.h5', '--geo_file1', output_dir + '/' + temp_geos[i-1] + '.h5', '--geo_file2', patches_dir[i+1] + '/' + network + '/rdr_geometry_msk.h5', '--ref', str(temp_refpoi[0]), str(temp_refpoi[1]), '--datatype', datatype, '--output', temp_name, '--outdir', output_dir + '/']

            temp_files.append(temp_name)
            temp_geos.append(temp_geo)

            scp_args = mu.separate_string_by_space(scp_args)
            # run concatenate_radarGeo.py
            print('concatenate_radarGeo.py',scp_args)
            mimtpy.concatenate_radarGeo.main(scp_args.split())   
    
######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)   
    
    project = inps.project[0]

    patches_file = inps.patches_file[0]
    patches_name = read_patches(patches_file)
    print('The patches to be concatenated are:')
    print(patches_name)

    refpoi_file = inps.ref[0]
    refpoi_list = read_ref_point(refpoi_file)

    network = inps.network[0]
    datatype = inps.datatype[0]

    concatenate_patches(project, patches_name, network, datatype, refpoi_list, inps)
     
######################################################################################
if __name__ == '__main__':
    main()
