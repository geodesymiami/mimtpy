#!/usr/bin/env python3
#################################################################
# Program is used for extract accumulated displacement of period#
# Author: Lv Xiaoran                                            #
# Created: August 2019                                          #
#################################################################

import os
import argparse
import string
import shutil
import numpy as np
import re
import scipy.io as sio

import mintpy
import mimtpy
import mimtpy.workflow
from mintpy.utils import readfile, writefile,utils as ut
from mimtpy.utils import  multitrack_utilities as mu

######################################################################################
EXAMPLE = """example:
  generate_horzvert.py $SCRATCHDIR/BalochistanSenAT/velocity/BalochistanSenAT.h5 $SCRATCHDIR/BalochistanSenDT/velocity/BalochistanSenDT.h5 --bbox 26.0 27.5 63.5 66.0 --outdir ./ 
  generate_horzvert.py $SCRATCHDIR/BalochistanSenAT/velocity/BalochistanSenAT.h5 $SCRATCHDIR/BalochistanSenDT/velocity/BalochistanSenDT.h5 --bbox 26.0 27.5 63.5 66.0 --reference_point 27.5 64.8 --azimuth 15 --outname hz.h5 up.h5 --outdir ./
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Generate horizontal and vertical datafile',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('DataSet', nargs=2, help='ascending or descending files\n')
                        
    parser.add_argument('--bbox', dest='SNWE', type=float, nargs=4, metavar=('S', 'N', 'W', 'E'),
                        help='Bounding box of area to be geocoded.\n' +
                        'Include the uppler left corner of the first pixel' +
                        '    and the lower right corner of the last pixel')
    parser.add_argument('--reference_point',dest='ref_poi', type=float, nargs='*',
                        help='reset reference point')
    parser.add_argument('--azimuth', dest='azimuth', type=float, default=90.0,
                        help='azimuth angle in degree (clockwise) of the direction of the horizontal movement\n' +
                             'default is 90.0 for E-W component, assuming no N-S displacement.\n' +
                             'i.e. azimuth angle of strike-slip fault\n\n')
    
    parser.add_argument('--outname',dest='outname', type=str, nargs='*', default=['horz.h5','vert.h5'],
                        help='output file name for vertical and horizontal components')
    parser.add_argument('--outdir',dest='outdir', type=str, nargs=1,
                        help='output directory')

    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  

    return inps

def horzvert(inps):
    """generate horzontal and vertical files."""
    # get ascending data and descending data
    dataset_asc = inps.DataSet[0]
    dataset_des = inps.DataSet[1]
   
    # get ascending and descending data name
    file_asc = os.path.split(dataset_asc)[1].split('.')[0]
    file_des = os.path.split(dataset_des)[1].split('.')[0]
     
    # output dir
    outdir = inps.outdir[0]

    # judge whether we should change the reference point
    refpoi_lalo = inps.ref_poi
    if refpoi_lalo:
        refpoi_lat = refpoi_lalo[0]
        refpoi_lon = refpoi_lalo[1]

        print('\nchanging ascending and descending data reference point')
        # change ascending data refpoi
        completion_status = os.system(mu.seperate_str_byspace(['reference_point.py', dataset_asc, '-l', refpoi_lat, '-L', refpoi_lon]))
        if completion_status == 1:
            raise Exception('error when runing reference_point.py for %s' % dataset_asc)
        asc_data, asc_atr = readfile.read(dataset_asc)
        asc_atr['REF_LAT'] = refpoi_lat
        asc_atr['REF_LON'] = refpoi_lon
        writefile.write(asc_data, out_file=dataset_asc, metadata=asc_atr)
        # change descending data refpoi 
        completion_status = os.system(mu.seperate_str_byspace(['reference_point.py', dataset_des, '-l', refpoi_lat, '-L', refpoi_lon]))
        if completion_status == 1:
            raise Exception('error when runing reference_point.py for %s' % dataset_des)
        des_data, des_atr = readfile.read(dataset_des)
        des_atr['REF_LAT'] = refpoi_lat
        des_atr['REF_LON'] = refpoi_lon
        writefile.write(des_data, out_file=dataset_des, metadata=des_atr)

    print('\ngo to the output dir {}'.format(outdir))
    os.chdir(outdir)
    
    # spatial range in lat/lon format
    SNWE = inps.SNWE

    # subset ascending and descending to the same spatial region
    print('\nsubset ascending data')
    sub_file_asc = 'subset_' + file_asc + '.h5'
    print(mu.seperate_str_byspace(['subset.py', dataset_asc, '-l', SNWE[0:2], '-L', SNWE[2:4], '-o', sub_file_asc]))
    completion_status = os.system(mu.seperate_str_byspace(['subset.py', dataset_asc, '-l', SNWE[0:2], '-L', SNWE[2:4], '-o', sub_file_asc]))
    if completion_status == 1:
        raise Exception('error when subset ascending data!')
 
    print('\nsubset descending data')
    sub_file_des = 'subset_' + file_des + '.h5'
    print(mu.seperate_str_byspace(['subset.py', dataset_des, '-l', SNWE[0:2], '-L', SNWE[2:4], '-o', sub_file_des]))
    completion_status = os.system(mu.seperate_str_byspace(['subset.py', dataset_des, '-l', SNWE[0:2], '-L', SNWE[2:4], '-o', sub_file_des]))
    if completion_status == 1:
        raise Exception('error when subset descending data!')

    # resolve ascending and descending to horz and vert data.
    azimuth = inps.azimuth
    outname = inps.outname
    horz_name = outname[0]
    vert_name = outname[1]

    print('\n runing asc_desc2horz_vert.py')
    print(mu.seperate_str_byspace(['asc_desc2horz_vert.py', sub_file_asc, sub_file_des, '--az', str(azimuth), '-o', horz_name, vert_name]))
    completion_status = os.system(mu.seperate_str_byspace(['asc_desc2horz_vert.py', sub_file_asc, sub_file_des, '--az', str(azimuth), '-o', horz_name, vert_name]))
    if completion_status == 1:
        raise Exception('error when running asc_desc2horz_vert.py!')
        
    # run H5UNW_to_geotiff.py
    # horizontal data
    scp_args2 = [horz_name, '--outdir', outdir, '--output', horz_name.split('.')[0] + '.tiff'] 
    scp_args2 = mu.seperate_str_byspace(scp_args2)
    print('\nH5UNW_to_geotiff.py',scp_args2)
    mimtpy.H5UNW_to_geotiff.main(scp_args2.split())

    # vertical data
    scp_args2 = [vert_name, '--outdir', outdir, '--output', vert_name.split('.')[0] + '.tiff']
    scp_args2 = mu.seperate_str_byspace(scp_args2)
    print('\nH5UNW_to_geotiff.py',scp_args2)
    mimtpy.H5UNW_to_geotiff.main(scp_args2.split())

    return

######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    print(inps)    

    horzvert(inps)
    
######################################################################################
if __name__ == '__main__':
    main()
