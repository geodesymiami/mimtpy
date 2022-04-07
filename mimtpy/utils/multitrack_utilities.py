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
import subprocess as subp
import matplotlib.path as path

import mintpy
import mintpy.workflow  #dynamic import for modules used by pysarApp workflow
from mintpy.objects import sensor
from mintpy.utils import ptime, readfile, writefile,utils as ut
from mintpy.objects import timeseries

def separate_string_by_space(arr, dst=""):
    """ parse list array(list item or values) to string """
    for k in arr:
        if isinstance(k,list):
            dst = dst + " " + separate_string_by_space(k)
        else:
            dst = dst + " " + str(k)
    return dst.strip()
    
def separate_filename_extension(path):
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
    completion_status = os.system(separate_string_by_space(['info.py', datafile, '--date', '>', 'date_list.txt']))
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
    
def velocity_displacement(inps):
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
    
def delete_temporalGeofile(datadir, key1, key2):
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

# the following part is used for screw dislocation forward and inversion model
def screw_disc(x, s, d, c, xc = 0):
    '''
    Function to calculate displacements/velocities due to slip on a deep 
    screw dislocation (infinitely long strike slip fault. 
    After Savage and Burford (1973).
    v = (s/pi)*arctan((x+xc)/d)
    INPUTS
        x = vector of distances from fault
        s = slip or slip rate on deep dislocation
        d = locking depth [same units as x]
        c = scalar offset in y [same unit as s]
        xc = offset to fault location [same units as x]
    OUTPUTS
        v = vector of displacements or velocities at locations defined by x
          [same units as s]
          
    USEAGE:
        v = deepdisloc(x, s, d)
    '''
    
    v = (s/np.pi) * np.arctan(x/d) + c
    
    return v

def fault_creep(x, s1, s2, d1, d2, c, xc=0):
    '''
    Model for interseismic strain accumulation and fault creep.
    Taken from Hussain et al. (2016).
    INPUTS
        x = vector of distances from fault
        s1 = slip or slip rate on deep dislocation
        s2 = creep or creep rate
        d1 = locking depth
        d2 = depth to bottom of creeping section
        c = scalar offset in y [same unit as s]
        xc = offset of fault location
    OUTPUTS
        v = vector of displacements or velocities at locations defined by x
          [same units as s]
    '''
    
    
    v = (s1/np.pi)*np.arctan((x+xc)/d1) + c - s2*((1/np.pi)*np.arctan((x+xc)/d2) + (x<=0)*0.5 - (x>0)*0.5)
    
    
    
    #(m(1)/pi)*atan(x./m(2)) + m(3) - m(4)*((1/pi)*atan(x./m(5)) + (x<=0)*1/2 - (x>0)*1/2);
    
    return v

def loglike(x, v, m, W):
    '''
    INPUTS
        x = vector of distances from fault
        v = velocities at locations defined by x
        m = model parameters, [0] = slip (mm/yr), [1] = locking depth (km), [2] = scalar offset (mm/yr)
        W = weight matrix (inverse of the VCM)
    OUTPUTS
        ll = value of the loglikelihood function
    '''
    
    #ll = np.sum((np.transpose(v-screw_disc(x, m[0], m[1], m[2]))*W*(v-screw_disc(x, m[0], m[1], m[2]))));
    #ll = np.sum(-0.01*(np.transpose(v-screw_disc(x, m[0], m[1], m[2]))*W*(v-screw_disc(x, m[0], m[1], m[2]))));
    diff = v - screw_disc(x, m[0], m[1], m[2])
    ll = np.nansum(-0.01 * np.dot(np.dot(diff,W),np.transpose(diff)))

    return ll

def logprior(m,m_min,m_max):
    '''
    INPUTS
        m = model values
        m_min = lower model limits
        m_max = upper model limits
    OUTPUTS
        lp = true if within limits, false if any aren't
    '''
    
    lp = np.all(np.all(m>=m_min) & np.all(m<=m_max))
    
    return lp

def rms_misfit(a,b):
    '''
    INPUTS
        a,b = two arrays of same length
    OUTPUTS
        rms = rms misfit between a and b (a-b)
    '''
    
    rms = np.sqrt(np.nanmean((a-b)**2))
    
    return rms

#-------------------------------------------------------------------------------

def run_inversion(n_iterations=10000):
    '''
    Bayesian monte carlo inversion taken from the main notebook.
    Re-packaged as a function here to minimise code repeats.
    '''
    # run inversion
    for ii in range(n_iterations):
        
        # propose model using different step sizes for each parameter
        m_trial = m_current.copy()
        m_trial[0] = m_trial[0] + np.random.normal(loc=0, scale=2.5, size=1)/1000 # slip rate
        m_trial[1] = m_trial[1] + np.random.normal(loc=0, scale=2.5, size=1)*1000 # locking depth
        m_trial[2] = m_trial[2] + np.random.normal(loc=0, scale=1, size=1)/1000 # offset
        
        # check limits and skip the rest of the loop if any parameter is invalid
        if not(lib.logprior(m_trial,m_min,m_max)):
            n_reject += 1
            models_saved[ii,:] = m_current
            continue
        
        # calculate likelihood for the current and trial models
        ll_current = lib.loglike(x, v, m_current, W)
        ll_trial = lib.loglike(x, v, m_trial, W)
        
        #print(np.exp(ll_trial-ll_current))
        
        # test whether to keep trial model
        #if np.exp(ll_current-ll_trial) > np.random.uniform(low=0, high=1,size=1):
        if np.exp(ll_trial-ll_current) > np.random.uniform(low=0, high=1,size=1):
        #if (ll_trial < ll_current) or (np.random.uniform(low=0, high=1,size=1) > 0.75):
            m_current = m_trial
            ll_current = ll_trial
            n_accept += 1
        else:
            n_reject += 1
        models_saved[ii,:] = m_current
        ll_saved[ii] = ll_current
        
    # convert back from metres to mm/yr and km
    models_saved[:,0] = models_saved[:,0] * 1000 # slip rate
    models_saved[:,1] = models_saved[:,1] / 1000 # locking depth
    models_saved[:,2] = models_saved[:,2] * 1000 # offset
    # find best fit model using min of likelihood function
    best_model = models_saved[np.nanargmax(ll_saved),:]

#-------------------------------------------------------------------------------

def profile_data(x,y,data,prof_start,prof_end,params):
    
    '''
    Generates a profile through gridded data.
    
    INPUTS:
    data = numpy array of values to profile
    x = vector of coords for the x axis
    y = vector of coords for the y axis
    prof_start = (x, y) pair for the start of the profile line
    prof_end = (x, y) pair for the end of the profile line
    params = dictionary of parameters for the profiler (currently nbins and width)
    
    '''
    
    xx,yy = np.meshgrid(x,y)
    
    prof_start = np.array(prof_start)
    prof_end = np.array(prof_end)
    
    # Profile dimensions relative to profile itself
    prof_dist = np.sqrt((prof_start[1]-prof_end[1])**2 + (prof_start[0]-prof_end[0])**2)
    prof_bin_edges = np.linspace(0, prof_dist ,params["nbins"]+1)    
    prof_bin_mids = (prof_bin_edges[:-1] + prof_bin_edges[1:]) / 2
    
    # Profile points in lat long space
    bin_mids = np.linspace(0,1,params["nbins"]+1)
    bin_grad = prof_end - prof_start
    x_mids = prof_start[0] + (bin_mids * bin_grad[0])
    y_mids = prof_start[1] + (bin_mids * bin_grad[1])
    
    # Gradient of line perpendicular to profile
    bin_grad_norm = (params["width"]/2) * bin_grad / np.linalg.norm(bin_grad)
    
    # Corner points of bins
    bin_x1 = x_mids + bin_grad_norm[1]
    bin_x2 = x_mids - bin_grad_norm[1]
    bin_y1 = y_mids - bin_grad_norm[0]
    bin_y2 = y_mids + bin_grad_norm[0]
    
    # Pre-allocate outputs
    bin_val = np.zeros_like((bin_x1[:-1]))
    bin_std = np.zeros_like(bin_val)
    
    # Trim data set to points inside any bin (improves run time)
    full_poly = path.Path([(bin_x1[0], bin_y1[0]), (bin_x1[-1], bin_y1[-1]), (bin_x2[-1], bin_y2[-1]), (bin_x2[0], bin_y2[0])])
    poly_points = full_poly.contains_points(np.transpose([xx.flatten(),yy.flatten()]))
    poly_points = poly_points.reshape(data.shape)
    trim_data = data[poly_points]
    trim_xx = xx[poly_points]
    trim_yy = yy[poly_points]
    
    # Loop through each bin identifying the points that they contain
    for ii in range(0,params["nbins"]):
                            
        poly_x = np.array([bin_x1[ii], bin_x1[ii+1], bin_x2[ii+1], bin_x2[ii]]);
        poly_y = np.array([bin_y1[ii], bin_y1[ii+1], bin_y2[ii+1], bin_y2[ii]]);
        
        poly = path.Path([(poly_x[0], poly_y[0]), (poly_x[1], poly_y[1]), (poly_x[2], poly_y[2]), (poly_x[3], poly_y[3])])
        
        poly_points = poly.contains_points(np.transpose([trim_xx,trim_yy]))
                            
        in_poly_vals = trim_data[poly_points]

        bin_val[ii] = np.nanmean(in_poly_vals)
    
    # get point cloud
    poly_x = np.array([bin_x1[0], bin_x1[-1], bin_x2[-1], bin_x2[0]])
    poly_y = np.array([bin_y1[0], bin_y1[-1], bin_y2[-1], bin_y2[0]])
    points_poly = np.vstack((poly_x,poly_y)).T
    points_poly = np.vstack((points_poly,np.array([points_poly[0,0],points_poly[0,1]])))
    
    poly = path.Path([(poly_x[0], poly_y[0]), (poly_x[1], poly_y[1]), (poly_x[2], poly_y[2]), (poly_x[3], poly_y[3])])
    poly_points = poly.contains_points(np.transpose([trim_xx,trim_yy]))
    points_val = trim_data[poly_points]
    points_x = trim_xx[poly_points]
    points_y = trim_yy[poly_points]
    
    prof_m = (prof_start[1] - prof_end[1]) / (prof_start[0] - prof_end[0])
    points_m = (prof_start[1] - points_y) / (prof_start[0] - points_x)
    points_prof_angle = np.arctan((points_m - prof_m) / (1 + prof_m * points_m))
    points2prof_start = np.sqrt((prof_start[1] - points_y)**2 + (prof_start[0] - points_x)**2)
    points_dist = points2prof_start * np.cos(points_prof_angle)
    
    return bin_val, prof_bin_mids, points_val, points_dist, points_poly

#-------------------------------------------------------------------------------

def profile_fault_intersection(prof_start,prof_end,fault_trace):
    '''
    Calculates the distance along a profile at which it intersects a fault, and 
    the angle between the two at this point.
    
    INPUTS:
        prof_start = (x,y) coord
        prof_end = (x,y) coords
        fault_trace = x and y coords of fault trace
        
    OUTPUTS:
        intersection_distance
        intersection_angle
    '''
    
    # loop through all fault segments to find intersection
    for ind in range(fault_trace.shape[0]-1):
        if intersect(prof_start,prof_end,fault_trace[ind,:],fault_trace[ind+1,:]):
            inter_ind = ind
            break
    
    # get coords for either side of intersection segment
    fault_inter_coords = fault_trace[inter_ind:inter_ind+2,:]
    
    # calc gradient of of profile line and fault segment
    prof_m = (prof_start[1] - prof_end[1]) / (prof_start[0] - prof_end[0])
    fault_m = (fault_inter_coords[0,1] - fault_inter_coords[1,1]) / (fault_inter_coords[0,0] - fault_inter_coords[1,0])
    
    # calculate intersection angle
    intersection_angle = np.arctan((fault_m - prof_m) / (1 + prof_m * fault_m));
    
    # calculate distance to intersection point (from https://stackoverflow.com/questions/3252194/numpy-and-line-intersections)
    s = np.vstack([prof_start,prof_end,fault_inter_coords[0,:],fault_inter_coords[1,:]])        # s for stacked
    h = np.hstack((s, np.ones((4, 1)))) # h for homogeneous
    l1 = np.cross(h[0], h[1])           # get first line
    l2 = np.cross(h[2], h[3])           # get second line
    x, y, z = np.cross(l1, l2)          # point of intersection
    intersection_point = np.array([x/z, y/z])
    intersection_distance = np.sqrt((intersection_point[0]-prof_start[0])**2 + (intersection_point[1]-prof_start[1])**2);    
    
    return intersection_distance, np.rad2deg(intersection_angle)


def ccw(A,B,C):
    return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])

def intersect(A,B,C,D):
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)
