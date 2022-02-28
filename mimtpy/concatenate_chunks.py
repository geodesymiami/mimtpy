#!/usr/bin/env python3
#################################################################
# Program is used for concatenate chunks for one track          #
# Author: Lv Xiaoran                                            #
# Created: October 2021                                         #
#################################################################

import os
import re
import sys
import argparse
import numpy as np
import shutil
import h5py
import copy

import mintpy
from mintpy.utils import ptime, readfile, writefile
from mintpy.objects import timeseries, geometry

import mimtpy
import mimtpy.workflow
from mimtpy.utils import multitrack_utilities as mu

BOOL_ZERO = np.bool_(0)
INT_ZERO = np.int16(0)
FLOAT_ZERO = np.float32(0.0)
CPX_ZERO = np.complex64(0.0)
COMPRESSION = 'lzf'
######################################################################################
NOTE = """ PLEASE NOTE:
1. The concatenated file of several chunks has the name as *_track.h5 
2. All data are extracted from S1*.he5 file.
3. The timeseries2velocity.py of Mintpy was used to calculate the velocity.
4. The velocity_std data was not concatenated
5. The reference point of the chunk have the smallest latitude is used as the reference poitn for the concatenated data.
   If you want to change the reference point, please use the Mintpy to change the reference point of concatenated data
6. The concatenated data are masked data based on the mask of each chunk
7. If you set datatype as HDFEOS, the concatenated data are S1*.he5 file containing displacement, velocity, temporalCoh, avgSpatialCoh and geoemtry.
"""

EXAMPLE = """example:
    concatenate_chunks.py $TE/MakranBigSenDT166.template  --datatype velocity
    concatenate_chunks.py --chunks MakranChunk2*SenDT166 --project MakranBigSenDT166 --dir mimtpy
    concatenate_chunks.py --chunks MakranChunk*SenDT166 --project MakranBigSenDT166
    concatenate_chunks.py --chunks MakranChunk2{5,6}SenDT166 --project MakranBigSenDT166 --dir mimtpy --datatype velocity
    concatenate_chunks.py --chunks $SCRATCHDIR/MakranChunk*SenDT166 --project $SCRATCHDIR/MakranBigSenDT166 --dir mimtpy --datatype velociy
    concatenate_chunks.py GujaratChunk21SenDT107/mintpy/S1_IW123_107_0520_0525_20160926_XXXXXXXX.he5  GujaratChunk22SenDT107/mintpy/S1_IW123_107_0517_0522_20160926_XXXXXXXX.he5 --project GujaratBigSenDT107 --datatype velocity --dir mimtpy
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Concatenate data',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=NOTE+'\n'+EXAMPLE)

    parser.add_argument('template_file', nargs='?', type=str, help='template file\n')

    parser.add_argument('--chunks', nargs='+', type=str, help='chunk name\n')
    
    parser.add_argument('--project', nargs=1, type=str, help='project name\n')

    parser.add_argument('--dir', type=str, nargs='+', default='mimtpy', help='output dir of concatenated file\n')
    
    parser.add_argument('--dryrun', action='store_true', default=False, help='used for timeseries. Only show the removed date for each chunk. Didnot concatenate chunks\n')
    
    parser.add_argument('--datatype', type=str, nargs=1, default='velocity', help='could be velocity, temporalCoherence, avgSpatialCoherence, timeseries, geometry and HDFEOS\n')

    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps


def concatenation_chunks(chunks, project, datatype, pro_outdir, inps):
    dataset_dirs = [os.path.realpath(chunk) + "/mintpy" for chunk in chunks]
    dataset_dirs.sort()
    chunk_number = [re.findall(r'[A-Za-z]+|[\d.]+', os.path.basename(chunk))[1] for chunk in chunks]
   
    date_chunks = []
    # extract velocity/coherence for each dataset
    for ii, dataset in enumerate(dataset_dirs):
        os.chdir(dataset)
        print('Go to project dir:', dataset)
        # creat output dir
        outdir = os.path.abspath(os.path.join(dataset, datatype))
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        print('the output dir for {} is {}.'.format(datatype, outdir))
 
        # find HDFEOS file
        HDFEOS_file = mu.find_HDFEOS_fullname(dataset)
        print('The HDFEOS file is {}'.format(HDFEOS_file))
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
        elif datatype == 'temporalCoherence':
            scp_args = [HDFEOS_file, datatype]
            scp_args = mu.seperate_str_byspace(scp_args)
            
            # run HDFEOS_to_geotiff.py
            print('HDFEOS_to_geotiff.py', scp_args)
            mimtpy.HDFEOS_to_geotiff.main(scp_args.split())
            shutil.move(os.path.join(dataset, datatype+'.h5'), os.path.join(outdir, datatype+'.h5'))
            # mask the coherence
            os.chdir(outdir)
            coh_data, coh_atr = readfile.read(datatype + '.h5')
            coh_data[maskfile == 0] = np.nan
            writefile.write(coh_data, out_file=datatype + '.h5', metadata=coh_atr)
        elif datatype == 'avgSpatialCoherence':
            scp_args = [HDFEOS_file, datatype]
            scp_args = mu.seperate_str_byspace(scp_args)
            
            # run HDFEOS_to_geotiff.py
            print('HDFEOS_to_geotiff.py', scp_args)
            mimtpy.HDFEOS_to_geotiff.main(scp_args.split())
            shutil.move(os.path.join(dataset, datatype+'.h5'), os.path.join(outdir, datatype+'.h5'))
            # mask the coherence
            os.chdir(outdir)
            coh_data, coh_atr = readfile.read(datatype + '.h5')
            coh_data[maskfile == 0] = np.nan
            writefile.write(coh_data, out_file=datatype + '.h5', metadata=coh_atr)
        elif datatype == 'timeseries':
            #displacement = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/observation/displacement')[0]
            bperp_date = h5py.File(HDFEOS_file,'r')
            data_date = np.array(bperp_date['/HDFEOS/GRIDS/timeseries/observation/date']).tolist()
            date_chunks.append(data_date)
        elif datatype == 'geometry':
            azi, atr = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/geometry/azimuthAngle')
            inc = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/geometry/incidenceAngle')[0]
            hgt = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/geometry/height')[0]
            lat = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/geometry/latitude')[0]
            lon = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/geometry/longitude')[0]
            sdM = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/geometry/shadowMask')[0]
            srd = readfile.read(HDFEOS_file, datasetName='/HDFEOS/GRIDS/timeseries/geometry/slantRangeDistance')[0]

            # mask the data
            azi[maskfile == 0] = np.nan
            inc[maskfile == 0] = np.nan
            hgt[maskfile == 0] = np.nan
            lat[maskfile == 0] = np.nan
            lon[maskfile == 0] = np.nan
            sdM[maskfile == 0] = np.nan
            srd[maskfile == 0] = np.nan
            
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
            outfile = outdir + '/' + datatype + '.h5' 
            writefile.write(dsDict, out_file=outfile, metadata=atr)     
           
    if datatype == 'timeseries':
        date_intersection = list(set(date_chunks[0]).intersection(*date_chunks[1:]))
        date_intersection.sort()
        if inps.dryrun:
            date_shared = copy.deepcopy(date_intersection)
            for date_chunk, dataset in zip(date_chunks, dataset_dirs):
                print('The removed date for project %s:' % dataset)
                date_diff = list(set(date_chunk).difference(set(date_shared)))
                #date_diff = list(set(date_chunk)^set(date_shared)).sort()
                if len(date_diff) == 0:
                    print('NONE')
                else:
                    print(date_diff)
            sys.exit(0)

    if datatype == 'timeseries':
        for dataset in dataset_dirs:
            os.chdir(dataset)
            # creat output dir
            outdir = os.path.abspath(os.path.join(dataset, datatype))
            if not os.path.isdir(outdir):
                os.makedirs(outdir)
            print('the output dir for {} is {}.'.format(datatype, outdir))

            # find HDFEOS file
            HDFEOS_file = mu.find_HDFEOS_fullname(dataset)
            print('The HDFEOS file is {}'.format(HDFEOS_file))
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
                
    # concatenate chunks for one track
    # generate output dir 
    track_dir = os.path.abspath(os.path.join(project + '/' + pro_outdir + '/' + datatype + '/'))
    if not os.path.isdir(track_dir):
        os.makedirs(track_dir)
    print('the output dir for concatenation is {}.'.format(track_dir))
    
    # start the loop
    temp_files = []
    pro_num = len(dataset_dirs)
    #lat_range = np.arange(lat_min, lat_max)
    lat_range = chunk_number
    if pro_num == 2:
        temp_name = datatype + '_track'
        scp_args = [dataset_dirs[0] + '/' + datatype + '/' + datatype + '.h5', dataset_dirs[1]  + '/' + datatype + '/' + datatype + '.h5', '--output', temp_name,  '--outdir', track_dir + '/']
        scp_args = mu.seperate_str_byspace(scp_args)
        # run concatenate_offset.py
        print('concatenate_offset.py',scp_args)
        mimtpy.concatenate_offset.main(scp_args.split())   
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
            # run concatenate_offset.py
            print('concatenate_offset.py',scp_args)
            mimtpy.concatenate_offset.main(scp_args.split())   
    # plot the concatenate_data
    os.chdir(track_dir)
    scp_args = [track_dir + '/' + datatype + '_track.h5', '--nodisplay', '--save']
    scp_args = mu.seperate_str_byspace(scp_args)
    mintpy.view.main(scp_args.split())

    #if datatype != 'timeseries':
    for dataset in dataset_dirs:
        os.chdir(dataset)
        print('Go to project dir:', dataset)
        # delete velocity folder
        print('Delete {} folder from {}'.format(datatype,dataset))
        outdir = os.path.abspath(os.path.join(dataset, datatype))
        shutil.rmtree(outdir)
 
def create_hdf5_dataset(group, dsName, data, max_digit=55, compression=COMPRESSION):
    """Create HDF5 dataset and print out message."""
    msg = 'create dataset {d:<{w}}'.format(d='{}/{}'.format(group.name, dsName), w=max_digit)
    msg += ' of {t:<10} in size of {s}'.format(t=str(data.dtype), s=data.shape)
    msg += ' with compression={c}'.format(c=compression)
    print(msg)

    if data.ndim == 1:
        dset = group.create_dataset(dsName,
                                    data=data,
                                    compression=compression)
    elif data.ndim == 2:
        dset = group.create_dataset(dsName,
                                    data=data,
                                    chunks=True,
                                    compression=compression)
    return dset

def write_hdf5_file(metadata, out_file, ts_file, tcoh_file, scoh_file, vel_file, geom_file):
    """Write HDF5 file in HDF-EOS5 format."""
    ts_obj = timeseries(ts_file)
    ts_obj.open(print_msg=False)
    dateList = ts_obj.dateList
    numDate = len(dateList)

    # Open HDF5 File
    f = h5py.File(out_file, 'w')
    print('create HDF5 file: {} with w mode'.format(out_file))
    max_digit = 55

    ##### Group - Observation
    gName = 'HDFEOS/GRIDS/timeseries/observation'
    print('create group   /{}'.format(gName))
    group = f.create_group(gName)

    ## O1 - displacement
    dsName = 'displacement'
    dsShape = (numDate, ts_obj.length, ts_obj.width)
    dsDataType = np.float32
    print(('create dataset /{d:<{w}} of {t:<10} in size of {s}'
           ' with compression={c}').format(d='{}/{}'.format(gName, dsName),
                                           w=max_digit,
                                           t='float32',
                                           s=dsShape,
                                           c=COMPRESSION))
    dset = group.create_dataset(dsName,
                                shape=dsShape,
                                maxshape=(None, dsShape[1], dsShape[2]),
                                dtype=dsDataType,
                                chunks=True,
                                compression=COMPRESSION)

    print('write data acquition by acquition ...')
    prog_bar = ptime.progressBar(maxValue=numDate)
    for i in range(numDate):
        dset[i, :, :] = readfile.read(ts_file, datasetName=dateList[i])[0]
        prog_bar.update(i+1, suffix='{}/{} {}'.format(i+1, numDate, dateList[i]))
    prog_bar.close()

    # attributes
    dset.attrs['Title'] = dsName
    dset.attrs['MissingValue'] = FLOAT_ZERO
    dset.attrs['_FillValue'] = FLOAT_ZERO
    dset.attrs['Units'] = 'meters'

    ## O2 - date
    dsName = 'date'
    data = np.array(dateList, dtype=np.string_)
    dset = create_hdf5_dataset(group, dsName, data)

    ## O3 - perp baseline
    dsName = 'bperp'
    data = np.array(ts_obj.pbase, dtype=np.float32)
    dset = create_hdf5_dataset(group, dsName, data)

    ## O4 - velocity
    dsName = 'velocity'
    data = readfile.read(vel_file)[0]
    dset = create_hdf5_dataset(group, dsName, data)
    # attributes
    dset.attrs['Title'] = dsName
    dset.attrs['MissingValue'] = FLOAT_ZERO
    dset.attrs['_FillValue'] = FLOAT_ZERO
    dset.attrs['Units'] = 'm/yr'

    ##### Group - Quality
    gName = 'HDFEOS/GRIDS/timeseries/quality'
    print('create group   /{}'.format(gName))
    group = f.create_group(gName)

    ## Q1 - temporalCoherence
    dsName = 'temporalCoherence'
    # read
    data = readfile.read(tcoh_file)[0]
    # write
    dset = create_hdf5_dataset(group, dsName, data)
    # attributes
    dset.attrs['Title'] = dsName
    dset.attrs['MissingValue'] = FLOAT_ZERO
    dset.attrs['_FillValue'] = FLOAT_ZERO
    dset.attrs['Units'] = '1'

    ## Q2 - avgSpatialCoherence
    dsName = 'avgSpatialCoherence'
    # read
    data = readfile.read(scoh_file)[0]
    # write
    dset = create_hdf5_dataset(group, dsName, data)
    # attributes
    dset.attrs['Title'] = dsName
    dset.attrs['MissingValue'] = FLOAT_ZERO
    dset.attrs['_FillValue'] = FLOAT_ZERO
    dset.attrs['Units'] = '1'

    ## Q3 - mask
    #dsName = 'mask'
    # read
    #data = readfile.read(mask_file, datasetName='mask')[0]
    # write
    #dset = create_hdf5_dataset(group, dsName, data)
    # attributes
    #dset.attrs['Title'] = dsName
    #dset.attrs['MissingValue'] = BOOL_ZERO
    #dset.attrs['_FillValue'] = BOOL_ZERO
    #dset.attrs['Units'] = '1'
    

    ##### Group - Write Geometry
    # Required: height, incidenceAngle
    # Optional: rangeCoord, azimuthCoord, azimuthAngle, slantRangeDistance, waterMask, shadowMask
    gName = 'HDFEOS/GRIDS/timeseries/geometry'
    print('create group   /{}'.format(gName))
    group = f.create_group(gName)

    geom_obj = geometry(geom_file)
    geom_obj.open(print_msg=False)
    for dsName in geom_obj.datasetNames:
        # read
        data = geom_obj.read(datasetName=dsName, print_msg=False)
        # write
        dset = create_hdf5_dataset(group, dsName, data)

        # attributes
        dset.attrs['Title'] = dsName
        if dsName in ['height', 'slantRangeDistance', 'bperp']:
            dset.attrs['MissingValue'] = FLOAT_ZERO
            dset.attrs['_FillValue'] = FLOAT_ZERO
            dset.attrs['Units'] = 'meters'

        elif dsName in ['incidenceAngle', 'azimuthAngle', 'latitude', 'longitude']:
            dset.attrs['MissingValue'] = FLOAT_ZERO
            dset.attrs['_FillValue'] = FLOAT_ZERO
            dset.attrs['Units'] = 'degrees'

        elif dsName in ['rangeCoord', 'azimuthCoord']:
            dset.attrs['MissingValue'] = FLOAT_ZERO
            dset.attrs['_FillValue'] = FLOAT_ZERO
            dset.attrs['Units'] = '1'

        elif dsName in ['waterMask', 'shadowMask']:
            dset.attrs['MissingValue'] = BOOL_ZERO
            dset.attrs['_FillValue'] = BOOL_ZERO
            dset.attrs['Units'] = '1'

    # Write Attributes to the HDF File
    print('write metadata to root level')
    for key, value in iter(metadata.items()):
        f.attrs[key] = value
    f.close()
    print('finished writing to {}'.format(out_file))

def pre_meta(ts_file):
    # read metadata from ts_file
    ts_obj = timeseries(ts_file)
    ts_obj.open(print_msg=False)
    meta = dict(ts_obj.metadata)

    meta['FILE_TYPE'] = 'HDFEOS'

    return meta

######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)   
    
    outdir = inps.dir[0] 
    datatype = inps.datatype[0] 
 
    if inps.template_file is not None:
        inpsdict =  mu.read_template(inps.template_file)
        if inpsdict['mimtpy.chunk'] == 'auto':
            identifiers = inpsdict['mimtpy.chunk.identifier']
            if identifiers == 'auto':
                chunks = inpsdict['mimtpy.chunk.chunks']
                project = inpsdict['mimtpy.chunk.project']
            else:
                chunks = identifiers[1]
                project = identifiers[0]
             
            scp_args = ['--chunks', chunks, '--project', project, '--datatype', datatype, '--dir', outdir]
            scp_args = mu.seperate_str_byspace(scp_args)

            # run concatenate_chunks.py
            print('concatenate_chunks.py', scp_args)

            os.system(mu.seperate_str_byspace(['concatenate_chunks.py', scp_args.split()]))
        else:
            raise ValueError('Please set parameters of mimtpy.chunk part')
        sys.exit(0)
    else: 
        chunks = inps.chunks
        project = inps.project[0]
    # all trans to full path
    chunks = [os.path.realpath(c) if c[0] != '/' else c for c in chunks]
    project = os.path.realpath(project) if project[0] != '/' else project
    # path validation
    for c in chunks:
        assert(os.path.isdir(c)), f"error chunk path `{c}`"
    
    # firstly concatenation chunks of each track
    # the output dir is project_name/mintpy/velocity/, for example: KokoxiliBigSenDT19/mintpy/velocity/
    # make sure the existence for project folders
    if datatype != 'HDFEOS':
        concatenation_chunks(chunks, project, datatype, outdir, inps)
    elif datatype == 'HDFEOS':
        for dtype in ['geometry', 'temporalCoherence', 'avgSpatialCoherence','velocity','timeseries']:
            print('process %s' % dtype)
            scp_args = ['--chunks', chunks, '--project', project, '--datatype', dtype, '--dir', outdir]
            scp_args = mu.seperate_str_byspace(scp_args)

            # run concatenate_chunks.py
            print('concatenate_chunks.py', scp_args)

            os_output = os.system(mu.seperate_str_byspace(['concatenate_chunks.py', scp_args.split()]))
            if os_output != 0:
                print('Error for datatype of %s' % dtype)
        
        # write concatenated HDFEOS file
        pro_dir = os.path.abspath(project)
        ts_data = pro_dir + '/mimtpy/timeseries/timeseries_track.h5'
        Tcoh_data = pro_dir + '/mimtpy/temporalCoherence/temporalCoherence_track.h5'
        Scoh_data = pro_dir + '/mimtpy/avgSpatialCoherence/avgSpatialCoherence_track.h5'
        vel_data = pro_dir + '/mimtpy/velocity/velocity_track.h5'
        geo_data = pro_dir + '/mimtpy/geometry/geometry_track.h5'

        meta = pre_meta(ts_data)
        bperp_date = h5py.File(ts_data,'r')
        data_date = np.array(bperp_date['/date']).tolist()
        start_date = data_date[0].decode('utf-8')
        end_date = data_date[-1].decode('utf-8')
        out_file = pro_dir + '/mimtpy/S1_' + start_date + '_' + end_date + '_track.he5'
        write_hdf5_file(metadata=meta,
                        out_file=out_file,
                        ts_file=ts_data,
                        tcoh_file=Tcoh_data,
                        scoh_file=Scoh_data,
                        vel_file=vel_data,
                        geom_file=geo_data)
 
######################################################################################
if __name__ == '__main__':
    main()
