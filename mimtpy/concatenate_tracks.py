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
    concatenate_tracks.py MakranBig*SenDT --datatype velocity
    concatenate_tracks.py MakranBig*SenDT --datatype velocity --reverse
    concatenate_tracks.py MakranBig*SenDT --outdir try --datatype velocity
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Concatenate data',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('track', nargs=1, type=str, help='data to be processed\n')

    parser.add_argument('--outdir', type=str, nargs='?', help='output dir\n')
    
    parser.add_argument('--datatype', type=str, nargs=1, help='could be velocity, temporalCoherence, avgSpatialCoherence\n')
    
    parser.add_argument('--reverse', action='store_true', default=False, help='whether reverse the track order. The default order is west to east\n')
    
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps

def search_tracks(file_dir, track, datatype):   
    L=[]
    for root, dirs, files in os.walk(file_dir):  
        for dir in dirs: 
            if dir.find(track.replace('*','')) != -1:
                if len(re.findall(r'[\d.]+', dir)) != 0:
                    L.append(os.path.join(root, dir, 'mimtpy', datatype))
    return L

def sort_tracks(track_names, datatype, inps):
    velocity_track = datatype + '_track.h5'
    
    longitudes = []

    for track in track_names:
        data_atr = readfile.read_attribute(track + '/' + velocity_track)
        X_FIRST = data_atr['X_FIRST']
        longitudes.append(X_FIRST)

    zipped = zip(track_names, longitudes)
    sort_zipped = sorted(zipped, key=lambda x:(x[1]))
    result = zip(*sort_zipped)
    track_names, longtidues = [list(x) for x in result]
    
    if inps.reverse:
        track_names.reverse()

    return track_names
    
def track_number_count(track_names, track):
    tmp = track.replace('*', '')
    
    track_numbers = []
    for track in track_names:
        for x in track.split('/'):
            if x.find(tmp) != -1:
                track_numbers.append(x.replace(tmp, ''))

    return track_numbers

def concatenation_tracks(track, file_dir, outdir, datatype, inps):
    datatype_track = datatype + '_track.h5'
 
    # search tracks meeting the requirements
    track_names_v1 = search_tracks(file_dir, track, datatype)

    # sort tracks according to their longitude order
    track_names = sort_tracks(track_names_v1, datatype, inps)

    track_number = len(track_names)
    temp_Files = []
    
    # get sorted track number
    track_numbers = track_number_count(track_names, track)
    outdir = outdir + '/mimtpy/' + datatype
   
    for i in np.arange(track_number - 1):  
        if i == 0:
            temp_name = datatype + '_' + track_numbers[0] + '_' + track_numbers[1]
            scp_args = [track_names[0] + '/' + datatype + '_track.h5', track_names[1] + '/' + datatype + '_track.h5', '--output', temp_name, '--outdir', outdir + '/']
        else:
            temp_name = temp_Files[i-1] + '_' + track_numbers[i+1]
            scp_args = [outdir + '/' + temp_Files[i-1] + '.h5', track_names[i+1] + '/' + datatype + '_track.h5', '--output', temp_name, '--outdir', outdir + '/']

        temp_Files.append(temp_name) 
        scp_args = mu.seperate_str_byspace(scp_args)

        # run track_offset.py
        print('track_offset.py',scp_args)
        mimtpy.track_offset.main(scp_args.split())        
######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)   
    
   
    track = inps.track[0]
    datatype = inps.datatype[0]
    if inps.outdir:
        outdir = inps.outdir
    else:
        outdir = track.replace('Big*', '')

    # generate output dir
    pro_dir = os.path.abspath(os.path.join(os.getenv('SCRATCHDIR') + '/' + outdir + '/mimtpy/' + datatype + '/'))
    if not os.path.isdir(pro_dir):
        os.makedirs(pro_dir)
    print('\nthe output dir for concatenation is {}.\n'.format(pro_dir))
    
    # concatenation tracks 
    file_dir = os.getenv('SCRATCHDIR')
    concatenation_tracks(track, file_dir, outdir, datatype, inps)
    
######################################################################################
if __name__ == '__main__':
    main()
