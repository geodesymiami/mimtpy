#!/usr/bin/env python3
#################################################################
# Program is part of MimtPy                                     #
# Author: Lv Xiaoran                                            #
# Created: September   2019                                     #
#################################################################
import os
import argparse
import string
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib.colors import LightSource
import shutil
import numpy as np
import re

import mintpy
import mintpy.workflow  #dynamic import for modules used by pysarApp workflow
from mintpy.objects import sensor
from mintpy.utils import ptime, readfile, writefile,utils as ut
from mintpy.objects import timeseries
from mimtpy.utils import multitrack_utilities

######################################################################################
EXAMPLE = """example:
  for singletrack:
  Note: startDate, endDate and outdir have default values, any of them can be not given:
    startDate default value = atr['START_DATE']
    endDate default value = atr['END_DATE']
    outdir default value = '$MODELDIR/project/Sen**/geodmod_startDate_endDate/'
    DataSet can be not given in template:
  save_for_modelling.py --template $MIMTFILES/Darbandikhan.txt --modelSoftware geodmod --outdir $MODELDIR
  save_for_modelling.py --templat $MIMTFILES/Darbandikhan.txt --modelSoftware gbis
  save_for_modelling.py --templat $MIMTFILES/Darbandikhan.txt --modelSoftware horzvert
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Prepare data for different software',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('--template', dest='templateFile',
                        help="Template file with geocoding options.")
                        
    parser.add_argument('--modelSoftware', dest='ModelSoftware',
                        help="name of software to be prepared for ")
                        
    parser.add_argument('--outdir', dest='outdir', nargs='?',default=os.getenv('SCRATCHDIR'),
                        help="output directory")

    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps    
   
def run_software(inpsdict,inps):
    """run save_geodmod/gbis.py in proper directory"""

    if inpsdict['DataSet'][0] == 'None':
        tempfilename = inps.templateFile
        folders =  multitrack_utilities.find_folder(multitrack_utilities.seprate_filename_extension(tempfilename)[1])
        print(folders)
    else:
        folders = inpsdict['DataSet']
        print(folders)
    # if datatype is HDFEOS,first judge whether they have same lat_step and lon_step
    if inpsdict['DataType'] == 'HDFEOS':
        flag =  multitrack_utilities.check_X_Y_step(folders)
        if flag:
            multitrack_run_software(inpsdict,folders,inps)
    else:
        multitrack_run_software(inpsdict,folders,inps)
        
def multitrack_run_software(inpsdict,folders,inps):
    """run save_geodmod/gbis for each track"""
    for project in folders:        
        os.chdir("".join([os.getenv('SCRATCHDIR')+'/'+project+'/mintpy/']))
        if inpsdict['DataType'] == 'HDFEOS':
            datafile =  multitrack_utilities.find_HDFEOS_fullname("".join([os.getenv('SCRATCHDIR')+'/'+project+'/mintpy/']))
        elif inpsdict['DataType'] == 'timeseries':
            datafile =  multitrack_utilities.find_timeseries("".join([os.getenv('SCRATCHDIR')+'/'+project+'/mintpy/']))
        elif inpsdict['DataType'] == 'ifgramStack':
            datafile = "".join([str(inpsdict['DataType'])+'.h5'])
        elif inpsdict['DataType'] == 'velocity':
            datafile = "".join([str(inpsdict['DataType'])+'.h5'])
        
        if inps.outdir == os.getenv('SCRATCHDIR'):
            inpsdict['outdir'] = inps.outdir + '/' + project + '/mintpy'
        else:
            inpsdict['outdir'] = inps.outdir + '/' + project
            
        if inps.ModelSoftware == 'geodmod':
            print(multitrack_utilities.seperate_str_byspace(['save_geodmod.py', datafile, '-b', inpsdict['SNWE'], '-y', inpsdict['latStep'], '-x', inpsdict['lonStep'], '-s', inpsdict['startDate'], '-e', inpsdict['endDate'], '-outdir', inpsdict['outdir']]))
            completion_status = os.system(multitrack_utilities.seperate_str_byspace(['save_geodmod.py', datafile, '-b', inpsdict['SNWE'], '-y', inpsdict['latStep'], '-x', inpsdict['lonStep'], '-s', inpsdict['startDate'], '-e', inpsdict['endDate'], '-outdir', inpsdict['outdir']]))
            if completion_status == 1:
                print('error when runing save_geodmod.py')
                exit(0)
        elif inps.ModelSoftware == 'gbis':
            if inpsdict['ref_lalo'] =='None':
                if inpsdict['mask_file'] =='None':
                    print('here!')
                    print(multitrack_utilities.seperate_str_byspace(['save_gbis_mimt.py', datafile, '-b', inpsdict['SNWE'], '-y', inpsdict['latStep'], '-x', inpsdict['lonStep'], '-s', inpsdict['startDate'], '-e', inpsdict['endDate'], '-outdir', inpsdict['outdir']]))
                    completion_status = os.system(multitrack_utilities.seperate_str_byspace(['save_gbis_mimt.py', datafile, '-b', inpsdict['SNWE'], '-y', inpsdict['latStep'], '-x', inpsdict['lonStep'], '-s', inpsdict['startDate'], '-e', inpsdict['endDate'], '-outdir', inpsdict['outdir']]))
                    if completion_status == 1:
                        raise Exception('error when runing save_gbis.py')  
                else:
                    print(multitrack_utilities.seperate_str_byspace(['save_gbis_mimt.py', datafile, '-b', inpsdict['SNWE'], '-y', inpsdict['latStep'], '-x', inpsdict['lonStep'], '-s', inpsdict['startDate'], '-e', inpsdict['endDate'], '-m', inpsdict['mask_file'], '-outdir', inpsdict['outdir']]))
                    completion_status = os.system(multitrack_utilities.seperate_str_byspace(['save_gbis_mimt.py', datafile, '-b', inpsdict['SNWE'], '-y', inpsdict['latStep'], '-x', inpsdict['lonStep'], '-s', inpsdict['startDate'], '-e', inpsdict['endDate'], '-m', inpsdict['mask_file'], '-outdir', inpsdict['outdir']]))
                    if completion_status == 1:
                        raise Exception('error when runing save_gbis.py')
            else:
                if inpsdict['mask_file'] =='None':
                    print(multitrack_utilities.seperate_str_byspace(['save_gbis_mimt.py', datafile, '-b', inpsdict['SNWE'], '-y', inpsdict['latStep'], '-x', inpsdict['lonStep'], '-s', inpsdict['startDate'], '-e', inpsdict['endDate'], '--ref-lalo', inpsdict['ref_lalo'], '-outdir', inpsdict['outdir']]))
                    completion_status = os.system(multitrack_utilities.seperate_str_byspace(['save_gbis_mimt.py', datafile, '-b', inpsdict['SNWE'], '-y', inpsdict['latStep'], '-x', inpsdict['lonStep'], '-s', inpsdict['startDate'], '-e', inpsdict['endDate'], '--ref-lalo', inpsdict['ref_lalo'], '-outdir', inpsdict['outdir']]))
                    if completion_status == 1:
                        raise Exception('error when runing save_gbis.py')    
                else:
                    print(multitrack_utilities.seperate_str_byspace(['save_gbis_mimt.py', datafile, '-b', inpsdict['SNWE'], '-y', inpsdict['latStep'], '-x', inpsdict['lonStep'], '-s', inpsdict['startDate'], '-e', inpsdict['endDate'], '-m', inpsdict['mask_file'], '--ref-lalo', inpsdict['ref_lalo'], '-outdir', inpsdict['outdir']]))
                    completion_status = os.system(multitrack_utilities.seperate_str_byspace(['save_gbis_mimt.py', datafile, '-b', inpsdict['SNWE'], '-y', inpsdict['latStep'], '-x', inpsdict['lonStep'], '-s', inpsdict['startDate'], '-e', inpsdict['endDate'], '-m', inpsdict['mask_file'], '--ref-lalo', inpsdict['ref_lalo'], '-outdir', inpsdict['outdir']]))
                    if completion_status == 1:
                        raise Exception('error when runing save_gbis.py')

######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    inpsdict =  multitrack_utilities.read_template2inps(inps)
    print(inpsdict)
    #import pdb
    #pdb.set_trace()
    run_software(inpsdict,inps)


######################################################################################
if __name__ == '__main__':
    main()