#!/usr/bin/env python3
#################################################################
# Program is part of Mimtpy                                     #
# Author: Lv Xiaoran                                            #
# Created: September 2019                                       #
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

def seperate_str_byspace(arr, dst=""):
    """ parse list array(list item or values) to string """
    for k in arr:
        if isinstance(k,list):
            dst = dst + " " + seperate_str_byspace(k)
        else:
            dst = dst + " " + str(k)
    return dst.strip()
    
def seprate_filename_extension(path):
    """ return directory(filepath), filename(shotname), extension """
    (filepath, tempfilename) = os.path.split(os.path.abspath(path))
    (filename, extension) = os.path.splitext(tempfilename)
    return filepath, filename, extension

def set_outdir(inps,modelsoftware):
    """set output directory"""
    dirname = inps.outdir+'/'+modelsoftware+'_'+inps.startDate+'_'+inps.endDate
    return dirname
    
def find_folder(tempfilename):
    """find the project folder and sort the folder in [*AT *DT]"""
    dir = os.getenv('SCRATCHDIR')
    folders = ["".join([dir +'/'])]
    project_folder = []
    for folder in folders:
        for x in os.listdir(folder):
            if os.path.isdir(os.path.join(folder,x)) and str.find(x,tempfilename)!= -1:
                project_folder.append(x)
    project_folder_sort = sorted(project_folder)
    return project_folder_sort
 
def find_timeseries(datadir):
    """find timeseries***.h5 file. The best results is timeseries_***_demErr.h5"""
    datafiles = []
    key1 = 'timeseries'
    key2 = 'Residual'
    for file in os.listdir(datadir):
        if os.path.splitext(file)[1] =='.h5':
            if str.find(file,key1) != -1 and str.find(file,key2) == -1:
                datafiles.append(file)
    datafile = []
    for file in datafiles:
        if len(file)>len(datafile):
            datafile = file
    return datafile

def find_HDFEOS_fullname(datadir):
    """find full name for S1 datatype"""
    datafiles = []
    for file in os.listdir(datadir):
        if os.path.splitext(file)[1] =='.he5':
                 datafiles.append(file)
    return datafiles[0]
    
def find_nearest_date(datadir,date):
    """get the date close to the given date"""
    datafile = find_timeseries(datadir)
    completion_status = os.system(seperate_str_byspace(['info.py', datafile, '--date', '>', 'date_list.txt']))
    if completion_status == 1:
        print('error when runing info.py')
        exit(0)
    if os.path.isfile('date_list.txt'):
        f = open('date_list.txt', 'r')
        lines = f.readlines() 
        f.close()        
    date_part = "".join(date)[0:6]
    date2 = []
    sub = 31
    for dates in lines:
        if str.find(dates,date_part) != -1:
            if abs(int(dates[6:8])-int(date[6:8]))<sub:
                sub = abs(int(dates[6:8])-int(date[6:8]))
                date2 = dates
    return date2.strip()

def find_start_end_date(datadir,inps):
    """find the startdate and enddate of each track"""   
    if not inps.startDate:
        startdate2 = inps.startDate
    if not inps.endDate:
        enddate2 = inps.endDate
    if inps.startDate:
        startdate2 = find_nearest_date(datadir,inps.startDate)
    if inps.endDate:
        enddate2 = find_nearest_date(datadir,inps.endDate)
    return startdate2,enddate2

def check_X_Y_step(folders):
    """for S1*h5 file, check whether the lat_step and lon_step are same for different projects"""
    x_step = []
    y_step = []
    for project in folders:
       os.chdir("".join([os.getenv('SCRATCHDIR')+'/'+project+'/mintpy/']))
       # find S1*.h5 whole name
       datafile = find_HDFEOS_fullname("".join([os.getenv('SCRATCHDIR')+'/'+project+'/mintpy/']))
       atr = readfile.read_attribute(datafile)
       x_step.append(atr['X_STEP'])
       y_step.append(atr['Y_STEP'])
    if len(set(x_step))!= 1:
        raise Exception("error! lon_Step between different tracks is not same!")
    elif len(set(y_step))!= 1:
        raise Exception("error! lat_Step between different tracks is not same!")
    else:
        return True

def read_template2inps(inps):
    """Read input template options into Namespace inps"""
    print('read input option from template file: ' + inps.templateFile)
    inps_dict = {}

    template = readfile.read_template(inps.templateFile)    
    template = ut.check_template_auto_value(template)
    
    prefix = inps.ModelSoftware+'.'
    key_list_tmp = [i for i in template.keys() if str.find(i,prefix) != -1]
    key_list = []
    for i in key_list_tmp:
        key_word = i.split('.')[-1]
        key_list.append(key_word)
    print(key_list)
    for key in key_list:
        value = template[prefix + key]
        if value:
            if key == 'DataSet':
                inps_dict[key] = list(tuple([i for i in value.split(',')]))
            elif key == 'DataType':
                inps_dict[key] = value
            elif key == 'SNWE':
                inps_dict[key] = list(tuple([float(i) for i in value.split(',')]))
            elif key in ['latStep', 'lonStep']:
                inps_dict[key] = float(value)
            elif key in ['startDate','endDate']:
                inps_dict[key] = value
            elif key in ['mask_file']:
                inps_dict[key] = value
            elif key in ['ref_lalo']:
                inps_dict[key] = list(tuple([float(i) for i in value.split(',')]))
            elif key in ['maskfile']:
                inps_dict[key] = value
            elif key in ['outname']:
                inps_dict[key] = value

    return inps_dict
    
def velo_disp(inps):
    """calculated displacement during startDate_endDate period based on linear assumption and velocity.h5"""
    data, atr = readfile.read('geo_velocity.h5')
    # calculate disp
    dt1, dt2 = ptime.date_list2vector([inps.startDate, inps.endDate])[0]
    data *= (dt2 - dt1).days / 365.25
    # displacement to phase
    range2phase =  -4. * np.pi / float(atr['WAVELENGTH'])
    data *= range2phase
    # write atr
    atr['PROCESSOR'] = 'roipac'
    atr['FILE_TYPE'] = '.unw'
    atr['UNIT'] = 'radian'
    out_file = 'geo_'+'{}_{}.unw'.format(inps.startDate, inps.endDate)
    writefile.write(data, out_file=out_file, metadata=atr)
    
def delete_tmpgeo(datadir, key1, key2):
    """delete all geo_*.h5 files in $MODLEDIR/project/SenAT(DT)/geodmod_startdate_enddate/"""
    for file in os.listdir(datadir):
        if os.path.splitext(file)[1] ==key2:
            if str.find(file,key1) != -1:
                os.remove(datadir+'/'+file)
                
def find_intersection_part(s1, s2): 
	m=[[0 for i in range(len(s2)+1)]  for j in range(len(s1)+1)]  
	mmax=0   
	p=0  
	for i in range(len(s1)):
		for j in range(len(s2)):
			if s1[i]==s2[j]:
				m[i+1][j+1]=m[i][j]+1
				if m[i+1][j+1]>mmax:
					mmax=m[i+1][j+1]
					p=i+1
	return s1[p-mmax:p],mmax   