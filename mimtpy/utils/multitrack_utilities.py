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
    datafile = find_HDFEOS_fullname(datadir)
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

def find_start_end_date(datadir,startDate,endDate):
    """find the startdate and enddate of each track"""   
    if startDate == 'None':
        startdate2 = startDate
    else:
        startdate2 = find_nearest_date(datadir,startDate)
    
    if endDate == 'None':
        enddate2 = endDate
    else:
        enddate2 = find_nearest_date(datadir,endDate)
    
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
                #inps_dict[key] = list(tuple([float(i) for i in value.split(',')]))
                inps_dict[key] = list(tuple([i for i in value.split(',')]))
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

# the following def are used for horzvert script                
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
    
def find_folder_horzvert(tempfilename,dir):
    """find the project folder and sort the folder in [*AT *DT]"""
    folders = ["".join([dir +'/'])]
    project_folder = []
    for folder in folders:
        for x in os.listdir(folder):
            if os.path.isdir(os.path.join(folder,x)) and str.find(x,tempfilename)!= -1:
                project_folder.append(x)
    project_folder_sort = sorted(project_folder)
    return project_folder_sort
    
def find_timeseries_horzvert(datadir):
    """find timeseries***.h5 file. The best results is timeseries_***_demErr.h5"""
    datafiles = []
    key1 = 'timeseries'
    key2 = 'Residual'
    key3 = 'msk'
    for file in os.listdir(datadir):
        if os.path.splitext(file)[1] =='.h5':
            if str.find(file,key1) != -1 and str.find(file,key2) == -1 and str.find(file,key3) == -1:
                datafiles.append(file)
    datafile = []
    for file in datafiles:
        if len(file)>len(datafile):
            datafile = file
    return datafile

def llh2xy(llh,origin):
    """Converts from longitude and latitude to local coorindates
      given an origin.  llh (lon; lat; height) and origin should
      be in decimal degrees. Note that heights are ignored and
      that xy is in km"""
    # llh is 2*N matrix, N is the number of points. the first line is longitude and the second line is latitude
    # orgion is 1*2 matrix,the order is longitude and latitude
    #Set ellipsoid constants (WGS84)
    a=6378137.0
    e=0.08209443794970

    #Convert to radians
    llh=llh * np.pi/180
    origin=origin * np.pi/180

    #Do the projection

    z=llh[1,:]!=0

    dlambda=llh[0,z]-origin[0]

    M=a*((1-e**2/4-3*(e**4)/64-5*(e**6)/256)*llh[1,z] - (3*(e**2)/8+3*(e**4)/32+45*(e**6)/1024)*np.sin(2*llh[1,z]) + (15*(e**4)/256 +45*(e**6)/1024)*np.sin(4*llh[1,z]) - (35*(e**6)/3072)*np.sin(6*llh[1,z]))

    M0=a*((1-e**2/4-3*(e**4)/64-5*(e**6)/256)*origin[1] - (3*(e**2)/8+3*(e**4)/32+45*(e**6)/1024)*np.sin(2*origin[1]) + (15*(e**4)/256 +45*(e**6)/1024)*np.sin(4*origin[1]) - (35*(e**6)/3072)*np.sin(6*origin[1]))
   
    N=a/np.sqrt(1-e**2*np.sin(llh[1,z])**2)
    E=dlambda*np.sin(llh[1,z])

    xy = np.zeros((2,len(dlambda)),dtype=float)
    xy[0,z]=N*(1 / np.tan(llh[1,z]))*np.sin(E)
    xy[1,z]=M-M0+N*(1 / np.tan(llh[1,z]))*(1-np.cos(E))

    #Handle special case of latitude = 0

    dlambda=llh[0,~z]-origin[0]
    xy[0,~z]=a*dlambda
    xy[1,~z]=-M0

    #Convert to km
   
    xy=xy/1000
 
    return xy

def calculate_LOS_value(inc,head,north_disp,east_disp,up_disp):
    # calculate LOS
    #Project displacement from LOS to Horizontal and Vertical components
    #    math for 3D: cos(theta)*Uz - cos(alpha)*sin(theta)*Ux + sin(alpha)*sin(theta)*Uy = Ulos
    #    math for 2D: cos(theta)*Uv - sin(alpha-az)*sin(theta)*Uh = Ulos   #Uh_perp = 0.0
    """inc is 1*1 value, head is 1*1 value,north_disp/east_disp/up_disp is n*1 matrix
        positive value means north/east/up"""
    inc_angle = inc
    head_angle = head
    
    az_angle = 90 
    inc_angle *= np.pi/180.

    # heading angle
    if head_angle < 0.:
        head_angle += 360.
    head_angle *= np.pi/180.

    # construct design matrix
    A_up = np.cos(inc_angle)
    A_north = np.sin(inc_angle) * np.sin(head_angle)
    A_east = np.sin(inc_angle) * np.cos(head_angle)
    
    # LOS
    los_sim = up_disp * A_up + north_disp * A_north - east_disp * A_east

    return los_sim

def read_template(fname, delimiter='=', print_msg=True):
    """Reads the template file into a python dictionary structure.
    Parameters: fname : str
                    full path to the template file
                delimiter : str
                    string to separate the key and value
                print_msg : bool
                    print message or not
    Returns:    template_dict : dict
                    file content
    Examples:
        tmpl = read_template(KyushuT424F610_640AlosA.template)
        tmpl = read_template(R1_54014_ST5_L0_F898.000.pi, ':')
        from mintpy.defaults.auto_path import isceAutoPath
        tmpl = read_template(isceAutoPath, print_msg=False)
    """
    template_dict = {}
    plotAttributeDict = {}
    insidePlotObject = False
    plotAttributes = []
    # the below logic for plotattributes object can be made much more simple
    # if we assume that any plot attribute coming after a > belongs to the
    # same object. Must Ask Falk and Yunjung if we can assume this to eliminate
    # all these conditionals

    if os.path.isfile(fname):
        f = open(fname, 'r')
        lines = f.readlines()
    elif isinstance(fname, str):
        lines = fname.split('\n')

    for line in lines:
        line = line.strip()
        # split on the 1st occurrence of delimiter
        c = [i.strip() for i in line.split(delimiter, 1)]
        if len(c) < 2 or line.startswith(('%', '#')):
            if line.startswith(">"):
                plotAttributeDict = {}
                insidePlotObject = True
            # otherwise, if previously inside attributes object, we are now outside
            # unless the line is a comment
            elif insidePlotObject and not line.startswith('%') and not line.startswith('#'):
                # just came from being inside plot object, but now we are outside
                insidePlotObject = False
                plotAttributes.append(plotAttributeDict)
            next  # ignore commented lines or those without variables
        else:
            atrName = c[0]
            atrValue = str.replace(c[1], '\n', '').split("#")[0].strip()
            atrValue = os.path.expanduser(atrValue)
            atrValue = os.path.expandvars(atrValue)

            if insidePlotObject:
                if is_plot_attribute(atrName):
                    plotAttributeDict[atrName] = atrValue
                else:
                    # just came from being inside plot object, but now we are outside
                    insidePlotObject = False
                    plotAttributes.append(plotAttributeDict)
                    template_dict[atrName] = atrValue

            elif atrValue != '':
                template_dict[atrName] = atrValue
    if os.path.isfile(fname):
        f.close()

    # what if no \n at end of file? write out last plot attributes dict
    if insidePlotObject:
        plotAttributes.append(plotAttributeDict)

    if len(plotAttributes) > 0:
        template_dict["plotAttributes"] = json.dumps(plotAttributes)

    return template_dict
