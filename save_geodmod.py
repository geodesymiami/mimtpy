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

import mintpy
import mintpy.workflow  #dynamic import for modules used by pysarApp workflow
from mintpy.objects import sensor
from mintpy.utils import readfile, writefile,utils as ut
from mintpy.objects import timeseries

######################################################################################
EXAMPLE = """example:
  for singletrack:
  save_geodmod.py timeseries_ECMWF_demErr.h5 -b 34.2 35.2 45.0 46.3 -y 0.001 -x 0.001 -s 20171117 -e 20180603 
  for multitrack:
  save_geodmod.py -t $MODELDIR/Darbandikhan.txt
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Prepare data for Geodmod software',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs='?', help='ascending and descending timeseries files\n')
    parser.add_argument('-t', '--template', dest='templateFile',
                        help="Template file with geocoding options.")
                        
    parser.add_argument('-ds', '--dataset', dest='DataSet',nargs='?',
                        help="name of dataset.Seperating by ','. ")                        
    parser.add_argument('-b', '--bbox', dest='SNWE', type=float, nargs=4, metavar=('S', 'N', 'W', 'E'),
                        help='Bounding box of area to be geocoded.\n' +
                        'Include the uppler left corner of the first pixel' +
                        '    and the lower right corner of the last pixel')
    parser.add_argument('-y', '--latstep', dest='latStep', type=float,
                        help='output pixel size in degree in latitude.')
    parser.add_argument('-x', '--lonstep', dest='lonStep', type=float,
                        help='output pixel size in degree in longitude.')
    parser.add_argument('-s','--startDate',dest='StartDate',nargs='?',
                        help='date1 of timeseires to be converted.The default is the StartDate')
    parser.add_argument('-e','--endDate',dest='EndDate',nargs='?',
                        help='date2 of timeseries to be converted.The default is the EndDate')
    
    parser.add_argument('-outdir','--outdir',dest='outdir',nargs=1,
                        help='output directory',default='geodmod')

    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    # read templatefile
    if inps.templateFile:
        inps = read_template2inps("".join(inps.templateFile), inps)
    else:    
        # default startDate and endDate
        print('come here')
        atr = readfile.read_attribute("".join(inps.file))
        if not inps.StartDate or inps.StartDate=='None':
            inps.StartDate=atr['START_DATE']
        if not inps.EndDate or inps.EndDate=='None':
            inps.EndDate=atr['END_DATE']
    print(inps)
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
            elif key == 'SNWE':
                inps_dict[key] = list(tuple([float(i) for i in value.split(',')]))
            elif key in ['latStep', 'lonStep']:
                inps_dict[key] = float(value)
            elif key in ['StartDate','EndDate']:
                inps_dict[key] = value

    inps.laloStep = [inps.latStep, inps.lonStep]
    if None in inps.laloStep:
        inps.laloStep = None
    return inps

    
def find_folder(tempfilename):
    """find the project folder and sort the folder in [*AT *DT]"""
    dir=os.getenv('SCRATCHDIR')
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
    datafiles=[]
    key1 = 'timeseries'
    key2 = 'Residual'
    for file in os.listdir(datadir):
        if os.path.splitext(file)[1] =='.h5':
            if str.find(file,key1) != -1 and str.find(file,key2) == -1:
                datafiles.append(file)
    datafile=[]
    for file in datafiles:
        if len(file)>len(datafile):
            datafile=file
    return datafile

def track_date(datafile,date):
    """get the date close to the given date"""
    completion_status=os.system(format_args(['info.py', datafile, '--date', '>', 'date_list.txt']))
    if completion_status == 1:
        print('error when runing info.py')
        exit(0)
    if os.path.isfile('date_list.txt'):
        f = open('date_list.txt', 'r')
        lines = f.readlines() 
        f.close()        
    date_part = "".join(date)[0:6]
    date2=[]
    sub=31
    for dates in lines:
        if str.find(dates,date_part) != -1:
            if abs(int(dates[6:8])-int(date[6:8]))<sub:
                sub=abs(int(dates[6:8])-int(date[6:8]))
                date2=dates
    return date2.strip()

def find_date(datafile,inps):
    """find the startdate and enddate of each track"""   
    if not inps.StartDate:
        startdate2=inps.StartDate
    if not inps.EndDate:
        enddate2=inps.EndDate
    if inps.StartDate:
        startdate2=track_date(datafile,inps.StartDate)
    if inps.EndDate:
        enddate2=track_date(datafile,inps.EndDate)
    return startdate2,enddate2

def run_save_geodmod(inps):
    """run save_geodmod.py in proper directory"""
    if not inps.DataSet:
        tempfilename=inps.templateFile
        folders = find_folder(seprate_filename_exten(tempfilename)[1])
        print(folders)
    else:
        folders = inps.DataSet
        print(folders)
    for project in folders:        
        os.chdir("".join([os.getenv('SCRATCHDIR')+'/'+project+'/PYSAR/']))
        datafile = find_timeseries("".join([os.getenv('SCRATCHDIR')+'/'+project+'/PYSAR/']))
        StartDate,EndDate = find_date(datafile,inps)
        print(format_args(['save_geodmod.py', datafile, '-b', inps.SNWE, '-y', inps.latStep, '-x', inps.lonStep, '-s', StartDate, '-e', EndDate, '-outdir', inps.outdir]))
        completion_status = os.system(format_args(['save_geodmod.py', datafile, '-b', inps.SNWE, '-y', inps.latStep, '-x', inps.lonStep, '-s', StartDate, '-e', EndDate, '-outdir', inps.outdir]))
        if completion_status == 1:
            print('error when runing save_geodmod.py')
            exit(0)

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
    print(dem.shape)
    # figure size
    ds_shape = tuple(reversed(dem.shape))
    fig_dpi = 300
    fig_size = [i / fig_dpi for i in ds_shape]
    print(fig_size)
    # color range
    disp_min = np.nanmin(dem) - 4000
    disp_max = np.nanmax(dem) + 2000
    # prepare shaded relief
    ls = LightSource(azdeg=315, altdeg=45)
    dem_shade = ls.shade(dem, vert_exag=0.3, cmap=plt.get_cmap('gray'), vmin=disp_min, vmax=disp_max)
    dem_shade[np.isnan(dem_shade[:, :, 0])] = np.nan
    print(dem_shade.shape)
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
    plt.show()


def processdata(inps):
    #use geocode.py and save_roipac.py to process data"
    #  geocode timeseries**.h5 file and get the deformation field of two time periods
    #  and geocode ifgramStack.h5 file and get the coherence of two time periods
    #  and geocode geometryRadar.h5 file and get the dem  
    #atr_asc = inps.file[0]
    atr_asc = inps.file
    
    if os.path.exists("".join(inps.outdir))=='False':
        os.mkdir("".join(inps.outdir))
   
    # process cor and dem dataset
    corname='temporalCoherence.h5'
    cmd_args = [corname, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    print("geocode.py", cmd_args)
    args_str = format_args(cmd_args)
    mintpy.geocode.main(args_str.split())
    
    demname='geometryRadar.h5'
    if not os.path.isfile(demname):
        demname_f='./INPUTS/geometryRadar.h5'
    else:
        demname_f='geometryRadar.h5'
    cmd_args = [demname_f, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    print("geocode.py", cmd_args)
    args_str = format_args(cmd_args)
    mintpy.geocode.main(args_str.split())
   
    #unw file
    cmd_args = [atr_asc, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    print("geocode.py", cmd_args)
    args_str = format_args(cmd_args)
    mintpy.geocode.main(args_str.split())
        
    #save dataset of unw cor and dem
    os.chdir("".join(inps.outdir))
    filename, extension = seprate_filename_exten("".join(atr_asc))[1:3]
    
    cmd_args = ['geo_'+filename+extension, "".join([inps.StartDate,'_',inps.EndDate])]
    print("save_roipac.py", cmd_args)
    asct_str = format_args(cmd_args)
    os.system(format_args(['save_roipac.py', asct_str.split()]))
    #mintpy.save_roipac.main(asct_str.split())     

    cmd_args = ['geo_temporalCoherence.h5', '-o', "".join(['geo_',inps.StartDate,'_',inps.EndDate,'.cor'])]    
    print("save_roipac.py", cmd_args)
    asct_str = format_args(cmd_args)
    os.system(format_args(['save_roipac.py', asct_str.split()]))
    #mintpy.save_roipac.main(asct_str.split())
    
    cmd_args = ['geo_geometryRadar.h5', 'height', '-o', 'srtm.dem']
    print("save_roipac.py", cmd_args)
    asct_str = format_args(cmd_args)
    os.system(format_args(['save_roipac.py', asct_str.split()]))
    #mintpy.save_roipac.main(asct_str.split())

######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
        
    if not inps.templateFile:
        print('single track!')
        processdata(inps)
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
