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
import h5py

import mintpy
from mintpy.utils import readfile, writefile
from mintpy.objects import timeseries

import mimtpy
import mimtpy.workflow
from mimtpy.utils import multitrack_utilities as mu
######################################################################################
EXAMPLE = """example:
    concatenate_chunks.py $TE/MakranBigSenDT166.template  --outdir mimtpy --datatype velocity
    concatenate_chunks.py --chunks MakranChunk2*SenDT166 --project MakranBigSenDT166 --outdir mimtpy --datatype velocity
    concatenate_chunks.py --chunks MakranChunk*SenDT166 --project MakranBigSenDT166 --outdir mimtpy --datatype velocity
    concatenate_chunks.py --chunks MakranChunk2{5,6}SenDT166 --project MakranBigSenDT166 --outdir mimtpy --datatype velocity
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Concatenate data',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('template_file', nargs='?', type=str, help='template file\n')

    parser.add_argument('--chunks', nargs='+', type=str, help='chunk name\n')
    
    parser.add_argument('--project', nargs=1, type=str, help='project name\n')

    parser.add_argument('--outdir', type=str, nargs='+', help='output dir of concatenated file\n')
    
    parser.add_argument('--datatype', type=str, nargs=1, help='could be velocity, temporalCoherence, avgSpatialCoherence and timeseries\n')

    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps

def chunks_number(chunks,file_dir):
    
    chunk_number = []
    L = []
    
    for chunk in chunks:
        ret2 = re.findall(r'[A-Za-z]+|[\d.]+', chunk)
        chunk_number.append(ret2[1]) 
        L.append(os.path.join(file_dir, chunk, 'mintpy'))
    
    return L, chunk_number

def concatenation_chunks(chunks, project, datatype, pro_outdir):
    #dataset_dirs = []
    # search for chunks number and name
    scratch_dir = os.getenv('SCRATCHDIR')
   
    dataset_dirs, chunk_number = chunks_number(chunks, scratch_dir)
    dataset_dirs.sort()
    chunk_number.sort()   
   
    date_chunks = []
    # extract velocity/coherence for each dataset
    for ii, dataset in enumerate(dataset_dirs):
        os.chdir(dataset)
        print('\nGo to project dir:', dataset)
        # creat output dir
        outdir = os.path.abspath(os.path.join(dataset, datatype))
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        print('\nthe output dir for {} is {}.\n'.format(datatype, outdir))
 
        # find HDFEOS file
        HDFEOS_file = mu.find_HDFEOS_fullname(dataset)
        print('\nThe HDFEOS file is {}'.format(HDFEOS_file))
        # extract mask
        maskfile = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/quality/mask')[0] 
        
        if datatype == 'velocity':
            scp_args = [HDFEOS_file]
            scp_args = mu.seperate_str_byspace(scp_args)
        
            # run timeseries2velocity.py
            print('timeseries2velocity.py',scp_args)
            mintpy.timeseries2velocity.main(scp_args.split())
            shutil.move(os.path.join(dataset, 'velocity.h5'), os.path.join(outdir, 'velocity.h5'))
            # mask the velocity
            os.chdir(outdir)
            vel_data, vel_atr = readfile.read('velocity.h5')
            vel_data[maskfile == 0] = np.nan
            writefile.write(vel_data, out_file='velocity.h5', metadata=vel_atr)
        elif datatype == 'coherence':
            scp_args = [HDFEOS_file, datatype]
            scp_args = mu.seperate_str_byspace(scp_args)
            
            # run HDFEOS_to_geotiff.py
            print('HDFEOS_to_geotiff.py', scp_args)
            mimtpy.HDFEOS_to_geotiff.main(scp_args.split())
            shutil.move(os.path.join(dataset, datatype+'.h5'), os.path.join(outdir, datatype+'.h5'))
            # mask the velocity
            os.chdir(outdir)
            coh_data, coh_atr = readfile.read(datatype + '.h5')
            coh_data[maskfile == 0] = np.nan
            writefile.write(coh_data, out_file=datatype + '.h5', metadata=coh_atr)
        elif datatype == 'timeseries':
            displacement = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/observation/displacement')[0]
            bperp_date = h5py.File(HDFEOS_file,'r')
            data_date = np.array(bperp_date['/HDFEOS/GRIDS/timeseries/observation/date']).tolist()
            date_chunks.append(data_date)
    if datatype == 'timeseries':
        #mmin_sceneNo = 10000
        #for dataset in dataset_dirs:
        #    os.chdir(dataset)
        #    # find HDFEOS file
        #    HDFEOS_file = mu.find_HDFEOS_fullname(dataset)
        #    bperp_date = h5py.File(HDFEOS_file,'r')
        #    data_bperp = bperp_date['/HDFEOS/GRIDS/timeseries/observation/bperp']
        #    data_date = np.array(bperp_date['/HDFEOS/GRIDS/timeseries/observation/date']).tolist()
        #    date_intersection = list(set(data_date).intersection(set(date_record2)))
        #    if mmin_sceneNo > len(date_intersection):
        #        mmin_sceneNo = len(date_intersection)
        #        date_record = date_intersection
        #date_record.sort()
        date_intersection = list(set(date_chunks[0]).intersection(*date_chunks[1:]))
        date_intersection.sort()
        for dataset in dataset_dirs:
            os.chdir(dataset)
            # creat output dir
            outdir = os.path.abspath(os.path.join(dataset, datatype))
            if not os.path.isdir(outdir):
                os.makedirs(outdir)
            print('\nthe output dir for {} is {}.\n'.format(datatype, outdir))

            # find HDFEOS file
            HDFEOS_file = mu.find_HDFEOS_fullname(dataset)
            print('\nThe HDFEOS file is {}'.format(HDFEOS_file))
            maskfile = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/quality/mask')[0]
            # process the date
            displacement, atr = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/observation/displacement')
            bperp_date = h5py.File(HDFEOS_file,'r')
            data_bperp = bperp_date['/HDFEOS/GRIDS/timeseries/observation/bperp']
            data_date = np.array(bperp_date['/HDFEOS/GRIDS/timeseries/observation/date']).tolist()
            #data_date = bperp_date[(dataset+'date')]
            # output info    
            ts_data_name = 'timeseries.h5'
            outfile = outdir + '/' + ts_data_name
            # getting flag
            # get the intersection of two list
            #date_intersection = list(set(data_date).intersection(set(date_record)))
            #date_intersection.sort()
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
            writefile.write(dsDict, out_file=outfile, metadata=atr)     
                
    # concatenate chunks for one track
    # generate output dir 
    track_dir = os.path.abspath(os.path.join(os.getenv('SCRATCHDIR') + '/' + project + '/' + pro_outdir + '/' + datatype + '/'))
    if not os.path.isdir(track_dir):
        os.makedirs(track_dir)
    print('\nthe output dir for concatenation is {}.\n'.format(track_dir))
    
    # start the loop
    temp_files = []
    pro_num = len(dataset_dirs)
    #lat_range = np.arange(lat_min, lat_max)
    lat_range = chunk_number
    if pro_num == 2:
        temp_name = datatype + '_track'
        scp_args = [dataset_dirs[0] + '/' + datatype + '/' + datatype + '.h5', dataset_dirs[1]  + '/' + datatype + '/' + datatype + '.h5', '--output', temp_name,  '--outdir', track_dir + '/']
        scp_args = mu.seperate_str_byspace(scp_args)
        # run track_offset.py
        print('track_offset.py',scp_args)
        mimtpy.track_offset.main(scp_args.split())   
    else:
        for i in np.arange(pro_num - 1):
            if i == 0:
                temp_name = datatype + '_lat_' + str(lat_range[0]) + '_' + str(lat_range[1])
                scp_args = [dataset_dirs[0] + '/' + datatype + '/' + datatype + '.h5', dataset_dirs[1]  + '/' + datatype + '/' + datatype + '.h5', '--output', temp_name,  '--outdir', track_dir + '/']
                
            elif i > 0 and i < (pro_num - 1 - 1):
                temp_name = temp_files[i-1] + '_' + str(lat_range[i+1])
                scp_args = [track_dir + '/' + temp_files[i-1] +'.h5', dataset_dirs[i+1]  + '/' + datatype + '/' + datatype + '.h5', '--output', temp_name, '--outdir', track_dir + '/']
            else:
                #temp_name = 'velocity_lat_' + str(lat_range[0]) + '_' + str(lat_range[-1])
                temp_name = datatype + '_track'
                scp_args = [track_dir + '/' + temp_files[i-1] +'.h5', dataset_dirs[i+1]  + '/' + datatype + '/' + datatype + '.h5', '--output', temp_name, '--outdir', track_dir + '/']

            temp_files.append(temp_name)
            scp_args = mu.seperate_str_byspace(scp_args)
            # run track_offset.py
            print('track_offset.py',scp_args)
            mimtpy.track_offset.main(scp_args.split())   
    if datatype != 'timeseries':
        for dataset in dataset_dirs:
            os.chdir(dataset)
            print('\nGo to project dir:', dataset)
            # delete velocity folder
            print('\n Delete {} folder from {}'.format(datatype,dataset))
            outdir = os.path.abspath(os.path.join(dataset, datatype))
            shutil.rmtree(outdir)
 
######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)   
    
    outdir = inps.outdir[0] 
 
    if inps.template_file is not None:
        inpsdict =  mu.read_template(inps.template_file)
        if inpsdict['mimtpy.chunk'] == 'yes':
            chunks = inpsdict['mimtpy.chunk.chunks']
            project = inpsdict['mimtpy.chunk.project']
            datatype = inpsdict['mimtpy.chunk.datatype']
            
            scp_args = ['--chunks', chunks, '--project', project, '--datatype', datatype, '--outdir', outdir]
            scp_args = mu.seperate_str_byspace(scp_args)

            # run concatenate_chunks.py
            print('concatenate_chunks.py', scp_args)

            os.system(mu.seperate_str_byspace(['concatenate_chunks.py', scp_args.split()]))
        else:
            raise ValueError('Please set parameters of mimtpy.chunk part')
    else: 
        chunks = inps.chunks
        project = inps.project[0]
        datatype = inps.datatype[0] 

    # firstly concatenation chunks of each track
    # the output dir is project_name/mintpy/velocity/, for example: KokoxiliBigSenDT19/mintpy/velocity/
    # make sure the existence for project folders
    concatenation_chunks(chunks, project, datatype, outdir)

######################################################################################
if __name__ == '__main__':
    main()
