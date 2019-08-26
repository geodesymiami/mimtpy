#!/usr/bin/env python3
#################################################################
# Program is used for extract accumulated displacement of period#
# Author: Lv Xiaoran                                            #
# Created: August 2019                                          #
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
from PIL import Image

import mintpy
import mintpy.workflow  #dynamic import for modules used by pysarApp workflow
from mintpy.objects import sensor
from mintpy.utils import ptime, readfile, writefile,utils as ut
from mintpy.objects import timeseries
from minsar.objects import message_rsmas

######################################################################################
EXAMPLE = """example:
  for singletrack:
  Note: startDate, endDate and outdir have default values, any of them can be not given:
    startDate default value = atr['START_DATE']
    endDate default value = atr['END_DATE']
    outdir default value = '$MODELDIR/project/Sen**/geodmod_startDate_endDate/'
  save_geodmod.py timeseries_ECMWF_demErr.h5 -b 34.2 35.2 45.0 46.3 -y 0.001 -x 0.001  
  save_geodmod.py ifgramStack.h5  -b 34.2 35.2 45.0 46.3 -y 0.001 -x 0.001 -s 20171117 -e 20171129 -outdir $MODELDIR/Darbandikhan/SenAT73/
  save_geodmod.py velocity.h5  -b 34.2 35.2 45.0 46.3 -y 0.001 -x 0.001 -s 20171117 -e 20171129 -outdir $MODELDIR/Darbandikhan/SenAT73/
  save_geodmod.py S1_IW23_026_0108_0113_20171117_XXXXXXXX.he5 -s 20171128 -e 20181210 
  for multitrack:
  Note: startDate, endDate and DataSet can be not given in template:
  save_geodmod.py -t $MIMTFILES/Darbandikhan.txt
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Prepare data for Geodmod software',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs='?', help='ascending or descending files\n')
    parser.add_argument('-t', '--template', dest='templateFile',
                        help="Template file with geocoding options.")
                        
    parser.add_argument('-ds', '--dataset', dest='DataSet',nargs='?',
                        help="name of dataset.Seperating by ','. ")
    parser.add_argument('-dt', '--datatype',dest='DataType',nargs='?',
                        help='clarify the type of data.[velocity, ifgramStack, timeseries, S1*.he5].Only used in template file')
    parser.add_argument('-b', '--bbox', dest='SNWE', type=float, nargs=4, metavar=('S', 'N', 'W', 'E'),
                        help='Bounding box of area to be geocoded.\n' +
                        'Include the uppler left corner of the first pixel' +
                        '    and the lower right corner of the last pixel')
    parser.add_argument('-y', '--latstep', dest='latStep', type=float,
                        help='output pixel size in degree in latitude.')
    parser.add_argument('-x', '--lonstep', dest='lonStep', type=float,
                        help='output pixel size in degree in longitude.')
    parser.add_argument('-s','--startDate',dest='startDate',nargs='?',
                        help='date1 of timeseires to be converted.The default is the StartDate')
    parser.add_argument('-e','--endDate',dest='endDate',nargs='?',
                        help='date2 of timeseries to be converted.The default is the EndDate')
    
    parser.add_argument('-outdir','--outdir',dest='outdir',nargs=1,
                        help='output directory')

    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    # read templatefile
    if inps.templateFile:
        inps = read_template2inps("".join(inps.templateFile), inps)
    else:    
        # default startDate and endDate
        if not os.path.isfile("".join(inps.file)):
            file = find_timeseries(os.getcwd())
        else:
            file = "".join(inps.file)    
        atr = readfile.read_attribute(file)
        if not inps.startDate or inps.startDate=='None':
            inps.startDate = atr['START_DATE']
        if not inps.endDate or inps.endDate=='None':
            inps.endDate = atr['END_DATE']
    return inps    


def read_template2inps(templatefile, inps):
    """Read input template options into Namespace inps"""
    print('read input option from template file: ' + templatefile)
    if not inps:
        inps = cmd_line_parse()
    inps_dict = vars(inps)

    template = readfile.read_template(templatefile)    
    template = ut.check_template_auto_value(template)
    
    prefix = 'geodmod.'
    key_list = [i for i in list(inps_dict.keys()) if prefix + i in template.keys()]
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

    inps.laloStep = [inps.latStep, inps.lonStep]
    if None in inps.laloStep:
        inps.laloStep = None
    return inps

def set_outdir(inps,pwdDir):
    """set output directory"""
    if not inps.outdir or inps.outdir[0] == 'None':
        projectTrack = pwdDir.split('/')[-2] #get projectSenAT***
        ret = re.findall(r"^(.+)(Sen[AD]T\d+)$", projectTrack)
        project = ret[0][0]
        Track = ret[0][1]
        dirname = "".join([os.getenv('MODELDIR')+'/'+project+'/'+Track+'/geodmod_'+inps.startDate+'_'+inps.endDate])
    else:
        dirname = inps.outdir[0]+'/geodmod_'+inps.startDate+'_'+inps.endDate
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

def find_S1_fullname(datadir):
    """find full name for S1 datatype"""

    datafiles = []
    key1 = 'S1'
    for file in os.listdir(datadir):
        if os.path.splitext(file)[1] =='.he5':
            if str.find(file,key1) != -1 :
                 datafiles.append(file)
    return datafiles[0]
    
def track_date(datadir,date):
    """get the date close to the given date"""
    datafile = find_timeseries(datadir)
    completion_status = os.system(format_args(['info.py', datafile, '--date', '>', 'date_list.txt']))
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

def find_date(datadir,inps):
    """find the startdate and enddate of each track"""   
    if not inps.startDate:
        startdate2 = inps.startDate
    if not inps.endDate:
        enddate2 = inps.endDate
    if inps.startDate:
        startdate2 = track_date(datadir,inps.startDate)
    if inps.endDate:
        enddate2 = track_date(datadir,inps.endDate)
    return startdate2,enddate2

def check_step(folders):
    """for S1*h5 file, check whether the lat_step and lon_step are same for different projects"""
    x_step = []
    y_step = []
    for project in folders:
       os.chdir("".join([os.getenv('SCRATCHDIR')+'/'+project+'/PYSARTEST/']))
       # find S1*.h5 whole name
       datafile = find_S1_fullname("".join([os.getenv('SCRATCHDIR')+'/'+project+'/PYSARTEST/']))
       atr = readfile.read_attribute(datafile)
       x_step.append(atr['X_STEP'])
       y_step.append(atr['Y_STEP'])
    if len(set(x_step))!= 1:
        raise Exception("error! lon_Step between different tracks is not same!")
    elif len(set(y_step))!= 1:
        raise Exception("error! lat_Step between different tracks is not same!")
    else:
        return True

def multitrack_run_save_geodmod(inps,folders):
    """run save_geodmod for each track"""
    for project in folders:        
        os.chdir("".join([os.getenv('SCRATCHDIR')+'/'+project+'/PYSARTEST/']))
        if inps.DataType == 'S1':
            datafile = find_S1_fullname("".join([os.getenv('SCRATCHDIR')+'/'+project+'/PYSARTEST/']))
        elif inps.DataType == 'timeseries':
            datafile = find_timeseries("".join([os.getenv('SCRATCHDIR')+'/'+project+'/PYSARTEST/']))
        elif inps.DataType == 'ifgramStack':
            datafile = "".join([str(inps.DataType)+'.h5'])
        elif inps.DataType == 'velocity':
            datafile = "".join([str(inps.DataType)+'.h5'])
        
        print(format_args(['save_geodmod.py', datafile, '-b', inps.SNWE, '-y', inps.latStep, '-x', inps.lonStep, '-s', inps.startDate, '-e', inps.endDate, '-outdir', inps.outdir]))
        completion_status = os.system(format_args(['save_geodmod.py', datafile, '-b', inps.SNWE, '-y', inps.latStep, '-x', inps.lonStep, '-s', inps.startDate, '-e', inps.endDate, '-outdir', inps.outdir]))
        if completion_status == 1:
            print('error when runing save_geodmod.py')
            exit(0)

def run_save_geodmod(inps):
    """run save_geodmod.py in proper directory"""
    if not inps.DataSet:
        tempfilename = inps.templateFile
        folders = find_folder(seprate_filename_exten(tempfilename)[1])
        print(folders)
    else:
        folders = inps.DataSet
        print(folders)
    # if datatype is S1,first judge whether they have same lat_step and lon_step
    if inps.DataType == 'S1':
        flag = check_step(folders)
        if flag:
            multitrack_run_save_geodmod(inps,folders)
    else:
        multitrack_run_save_geodmod(inps,folders)
    
def format_args(arr, dst=""):
    """ parse list array(list item or values) to string """
    for k in arr:
        if isinstance(k,list):
            dst = dst + " " + format_args(k)
        else:
            dst = dst + " " + str(k)
    return dst.strip()
   
def seprate_filename_exten(path):
    """ return directory(filepath), filename(shotname), extension """
    (filepath, tempfilename) = os.path.split(os.path.abspath(path))
    (filename, extension) = os.path.splitext(tempfilename)
    return filepath, filename, extension

def delete_tmpgeo(datadir, key1, key2):
    """delete all geo_*.h5 files in $MODLEDIR/project/SenAT(DT)/geodmod_startdate_enddate/"""
    for file in os.listdir(datadir):
        if os.path.splitext(file)[1] ==key2:
            if str.find(file,key1) != -1:
                os.remove(datadir+'/'+file)
    
def write_rsc_file(inps,in_file,out_file):
    """ write rsc file for Geodmod just estract several properities from rsc file"""
    # read file
    meta = readfile.read_roipac_rsc(in_file)
    # initiate dict
    rsc = dict()
    rsc['FILE_DIR'] = "".join(inps.outdir)
    rsc['FILE_LENGTH'] = meta["FILE_LENGTH"]
    rsc['WIDTH'] = meta["WIDTH"]
    rsc['XMIN'] = 0
    rsc['XMAX'] = int(meta["WIDTH"]) - 1
    rsc['YMIN'] = 0
    rsc['YMAX'] = int(meta["FILE_LENGTH"]) - 1
    rsc['X_FIRST'] = float(meta["X_FIRST"])
    rsc['Y_FIRST'] = float(meta["Y_FIRST"])
    rsc['X_STEP'] = float(meta["X_STEP"])
    rsc['Y_STEP'] = float(meta["Y_STEP"])
    rsc['X_UNIT'] = 'degrees'
    rsc['Y_UNIT'] = 'degrees'
    rsc['RLOOKS'] = meta["RLOOKS"]
    rsc['ALOOKS'] = meta["ALOOKS"]
    rsc['Z_OFFSET'] = 0
    rsc['Z_SCALE'] = 1
    rsc['PROJECTION'] = 'LATLON'
    rsc['DATE12'] = '111111-222222'
    # write rsc file
    writefile.write_roipac_rsc(rsc, out_file, print_msg=True)
    return out_file 

def dem_jpeg(dem_file):
    """generate dem.jepg file based on Yunjun's code"""
    out_file = dem_file+'.jpeg'
    rsc_file = out_file+'.rsc'
    shutil.copy2(dem_file+'.rsc', rsc_file)
    # read data
    dem = readfile.read(dem_file)[0]
    print('dem.shape:',dem.shape)
    # figure size
    ds_shape = tuple(reversed(dem.shape))
    fig_dpi = 300
    fig_size = [i / fig_dpi for i in ds_shape]
    print('fig size:', fig_size)
    # color range
    disp_min = np.nanmin(dem) - 4000
    disp_max = np.nanmax(dem) + 2000

    # prepare shaded relief
    ls = LightSource(azdeg=315, altdeg=45)
    dem_shade = ls.shade(dem, vert_exag=0.3, cmap=plt.get_cmap('gray'), vmin=disp_min, vmax=disp_max)
    dem_shade[np.isnan(dem_shade[:, :, 0])] = np.nan
    print('dem_shade.shape:', dem_shade.shape)
    # plot
    fig, ax = plt.subplots(figsize=fig_size)
    ax.imshow(dem_shade, interpolation='spline16', origin='upper')
    # get rid of whitespace on the side
    ax.axis('off')
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    fig.subplots_adjust(left=0,right=1,bottom=0,top=1)
    # output
    print('save figure to file {}'.format(out_file))
    plt.savefig(out_file, transparent=True, dpi=300, pad_inches=0.0)
    # resize to desired size  (FA 8/19, unclear why size is wrong))
    im = Image.open(out_file)
    im_out = im.resize(dem.shape, Image.NEAREST) 
    im_out.save(out_file)

    plt.show()

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

def process_geocode(inps):
    """process temporalCoherence.h5 and geometryRadar.h5 file"""
    # process cor and dem dataset
    if os.path.exists("".join(inps.outdir))==False:
        os.makedirs("".join(inps.outdir))
        
    # geocode coherence
    fname = 'temporalCoherence.h5'
    cmd_args = [fname, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    args_str = format_args(cmd_args)
    print("geocode.py", cmd_args)
    message_rsmas.log(os.getcwd(),'geocode.py ' + args_str)
    mintpy.geocode.main(args_str.split())
    
    # geocode dem
    fname = 'geometryRadar.h5'
    if not os.path.isfile(fname):
        fname = './inputs/geometryRadar.h5'
    else:
        fname = 'geometryRadar.h5'
    cmd_args = [fname, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    args_str = format_args(cmd_args)
    print("geocode.py", cmd_args)
    message_rsmas.log(os.getcwd(),'geocode.py ' + args_str)
    mintpy.geocode.main(args_str.split())

def process_saveroi(inps):    
    #save_roipac
    cmd_args = ['geo_temporalCoherence.h5', '-o', "".join(['geo_',inps.startDate,'_',inps.endDate,'.cor'])]    
    print("save_roipac.py", cmd_args)
    args_str = format_args(cmd_args)
    os.system(format_args(['save_roipac.py', args_str.split()]))
    
    cmd_args = ['geo_geometryRadar.h5', 'height', '-o', 'srtm.dem']
    args_str = format_args(cmd_args)
    print("save_roipac.py", cmd_args)
    message_rsmas.log(os.getcwd(),'save_roipac.py ' + args_str)
    os.system(format_args(['save_roipac.py', args_str.split()]))

def process_time(inps):
    """geocode timeseries**.h5 file and get the deformation field of two time periods"""
    atr_asc = inps.file
   
    #unw file
    cmd_args = [atr_asc, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    args_str = format_args(cmd_args)
    print("geocode.py", cmd_args)
    message_rsmas.log(os.getcwd(),'geocode.py ' + args_str)
    mintpy.geocode.main(args_str.split())
        
    #save dataset of unw cor and dem
    os.chdir("".join(inps.outdir))
    filename, extension = seprate_filename_exten("".join(atr_asc))[1:3]
    
    cmd_args = ['geo_'+filename+extension, "".join([inps.startDate,'_',inps.endDate])]
    args_str = format_args(cmd_args)
    print("save_roipac.py", cmd_args)
    message_rsmas.log(os.getcwd(),'save_roipac.py ' + args_str)
    os.system(format_args(['save_roipac.py', args_str.split()]))

    process_saveroi(inps)
    delete_tmpgeo(inps.outdir,'geo_','.h5')

def process_ifgS(inps):
    """process ifgramStack.h5 file"""
    if os.path.exists("".join(inps.outdir))==False:
        os.makedirs("".join(inps.outdir))
    
    # dem file
    demname='geometryRadar.h5'
    if not os.path.isfile(demname):
        demname_f = './inputs/geometryRadar.h5'
    else:
        demname_f = 'geometryRadar.h5'
    cmd_args = [demname_f, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    args_str = format_args(cmd_args)
    print("geocode.py", cmd_args)
    message_rsmas.log(os.getcwd(),'geocode.py ' + args_str)
    mintpy.geocode.main(args_str.split())
    
    #ifgramStack file
    atr_asc = ['./inputs/'+inps.file]
    cmd_args = [atr_asc, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    args_str = format_args(cmd_args)
    print("geocode.py", cmd_args)
    message_rsmas.log(os.getcwd(),'geocode.py ' + args_str)
    mintpy.geocode.main(args_str.split())
        
    #save dataset of unw cor and dem
    os.chdir("".join(inps.outdir))
    filename, extension = seprate_filename_exten("".join(atr_asc))[1:3]
    
    cmd_args = ['geo_'+filename+extension, "".join(['unwrapPhase-',inps.startDate,'_',inps.endDate])]
    args_str = format_args(cmd_args)
    print("save_roipac.py", cmd_args)
    message_rsmas.log(os.getcwd(),'save_roipac.py ' + args_str)
    os.system(format_args(['save_roipac.py', args_str.split()]))   

    cmd_args = ['geo_'+filename+extension, "".join(['coherence-',inps.startDate,'_',inps.endDate])]
    args_str = format_args(cmd_args)
    print("save_roipac.py", cmd_args)
    message_rsmas.log(os.getcwd(),'save_roipac.py ' + args_str)
    completion_status = os.system(format_args(['save_roipac.py', args_str.split()])) 
    
    cmd_args = ['geo_geometryRadar.h5', 'height', '-o', 'srtm.dem']
    args_str = format_args(cmd_args)
    print("save_roipac.py", cmd_args)
    message_rsmas.log(os.getcwd(),'save_roipac.py ' + args_str)
    os.system(format_args(['save_roipac.py', args_str.split()]))
    
    delete_tmpgeo(inps.outdir,'geo_','.h5')

def process_vel(inps):
    """process velocity.h5 file"""
    atr_asc = inps.file
   
    #velocity file
    cmd_args = [atr_asc, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    args_str = format_args(cmd_args)
    print("geocode.py", cmd_args)
    message_rsmas.log(os.getcwd(),'geocode.py ' + args_str)
    mintpy.geocode.main(args_str.split())
    
    os.chdir("".join(inps.outdir))
    process_saveroi(inps)
    print('save unw file')
    velo_disp(inps)
    
    delete_tmpgeo(inps.outdir,'geo_','.h5')

def process_S1(inps):
    """process S1*.h5 file"""

    atr_asc = inps.file
        
    #save dataset of unw cor and dem
    filename, extension = seprate_filename_exten("".join(atr_asc))[1:3]
    cmd_args = [filename+extension, "".join(['displacement-',inps.startDate,'_',inps.endDate]), '-o', "".join(['geo_',inps.startDate,'_',inps.endDate,'.unw'])]
    args_str = format_args(cmd_args)
    print("save_roipac.py", cmd_args)
    message_rsmas.log(os.getcwd(),'save_roipac.py ' + args_str)
    os.system(format_args(['save_roipac.py', args_str.split()]))   

    cmd_args = [filename+extension, 'temporalCoherence', '-o', "".join(['geo_',inps.startDate,'_',inps.endDate,'.cor'])]
    args_str = format_args(cmd_args)
    print("save_roipac.py", cmd_args)
    message_rsmas.log(os.getcwd(),'save_roipac.py ' + args_str)
    completion_status=os.system(format_args(['save_roipac.py', args_str.split()])) 
    
    cmd_args = [filename+extension, 'height', '-o', 'srtm.dem']
    args_str = format_args(cmd_args)
    print("save_roipac.py", cmd_args)
    message_rsmas.log(os.getcwd(),'save_roipac.py ' + args_str)
    os.system(format_args(['save_roipac.py', args_str.split()]))
    
    # mv 
    if os.path.exists(inps.outdir)==False:
        os.makedirs(inps.outdir)
    else:
        shutil.rmtree(inps.outdir)
        os.makedirs(inps.outdir)
    key1 = 'geo_'
    key2 = 'srtm'
    for file in os.listdir(os.getcwd()):
        if str.find(file,key1) != -1 or str.find(file,key2) != -1:
            shutil.move(file,inps.outdir) 
    os.chdir("".join(inps.outdir))  

    delete_tmpgeo(inps.outdir,'geo_','.h5')    
    
######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    print(inps)    
    if not inps.templateFile:
        print('single track!')
        inps.startDate,inps.endDate = find_date(os.getcwd(),inps)
        inps.outdir = set_outdir(inps,os.getcwd())
        print(inps)
        if str.find(inps.file,'ifgramStack.h5') != -1:
            process_ifgS(inps)
        elif str.find(inps.file,'velocity.h5') != -1:
            process_geocode(inps)
            process_vel(inps)
        elif str.find(inps.file,'timeseries') != -1:
            process_geocode(inps)
            process_time(inps)
        else:
            process_S1(inps)
        # rename *.rsc1 to *.rsc
        outfile = format_args(['srtm.dem' + '.rsc'])
        write_rsc_file(inps,outfile,format_args(['srtm.dem' +'.rsc1']))
        os.remove(outfile)
        print('rename *.rsc1 to *.rsc')
        os.rename(format_args(['srtm.dem' +'.rsc1']),outfile)
        # generate dem.jpeg
        dem_jpeg('srtm.dem')
        
    else:
        print('multi track!')
        run_save_geodmod(inps)

######################################################################################
if __name__ == '__main__':
    main()
