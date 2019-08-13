#!/usr/bin/env python3
#################################################################
# Program is used for extract accumulated displacement of period#
# Author: Lv Xiaoran                                            #
# Created: August 2019                                          #
#################################################################

import os
import argparse
import string

import mintpy
import mintpy.workflow  #dynamic import for modules used by pysarApp workflow
from mintpy.objects import sensor
from mintpy.utils import readfile, writefile,utils as ut
from mintpy.objects import timeseries

######################################################################################
EXAMPLE = """example:
  for singletrack:
  save_geodmod.py -f timeseries_ECMWF_demErr.h5 -b 34.2 35.2 45.0 46.3 -y 0.001 -x 0.001 -startDate 20171117 -endDate 20180603 -outdir $MODELDIR 
  for multitrack:
  save_geodmod.py -t $MODELDIR/Darbandikhan.txt -outdir $MODELDIR
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Prepare data for Geodmod software',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    #parser.add_argument('file', nargs=1, help='ascending and descending timeseries files\n')
    parser.add_argument('-f','--file', dest='File', nargs=1, help='ascending and descending timeseries files\n')
    parser.add_argument('-t', '--template', dest='templateFile',
                        help="Template file with geocoding options.")
                        
    parser.add_argument('-m', '--multitrack', dest='multiTrack',
                        help="flag about whether use multitrack data.")
    parser.add_argument('-ds1', '--dataset1', dest='DataSet1',
                        help="name of dataset1.")
    parser.add_argument('-ds2', '--dataset2', dest='DataSet2',
                        help="name of dataset2.")                          
    parser.add_argument('-b', '--bbox', dest='SNWE', type=float, nargs=4, metavar=('S', 'N', 'W', 'E'),
                        help='Bounding box of area to be geocoded.\n' +
                        'Include the uppler left corner of the first pixel' +
                        '    and the lower right corner of the last pixel')
    parser.add_argument('-y', '--latstep', dest='latStep', type=float,
                        help='output pixel size in degree in latitude.')
    parser.add_argument('-x', '--lonstep', dest='lonStep', type=float,
                        help='output pixel size in degree in longitude.')
    parser.add_argument('-startDate','--startDate',dest='StartDate',nargs='?',help='date1 of timeseires to be converted.The default is the StartDate')
    parser.add_argument('-endDate','--endDate',dest='EndDate',nargs='?',help='date2 of timeseries to be converted.The default is the EndDate')
    
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
        print('come here')
        atr = readfile.read_attribute(format_args(inps.File))
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
            if key == 'multiTrack':
                inps_dict[key] = value
            elif key == 'DataSet1':
                inps_dict[key] = value
            elif key == 'DataSet2':
                inps_dict[key] = value
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
    """find the project folder"""
    dir=os.getenv('SCRATCHDIR')
    folders = ["".join([dir +'/'])]
    
    for folder in folders:
        folders += [os.path.join(folder,x) for x in os.listdir(folder)\
                   if os.path.isdir(os.path.join(folder,x) and tempfilename in x)]
    return folders
 
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

def des_date(datafile,date):
    """get the descending date close to ascending date"""
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
    """find the startdate and enddate of descending track"""   
    if not inps.StartDate:
        startdate2=inps.StartDate
    if not inps.EndDate:
        enddate2=inps.EndDate
    if inps.StartDate:
        startdate2=des_date(datafile,inps.StartDate)
    if inps.EndDate:
        enddate2=des_date(datafile,inps.EndDate)
    return startdate2,enddate2

def run_save_geodmod(inps):
    """run save_geodmod.py in proper directory"""
    if not inps.DataSet1 or not inps.DataSet2:
        tempfilename=inps.templateFile
        folders = find_folder(tempfilename)
        project1 = folders[0]
        project2 = folders[1]
    else:
        project1 = inps.DataSet1
        project2 = inps.DataSet2
    
    os.chdir("".join([os.getenv('SCRATCHDIR')+'/'+project1+'/PYSARTEST/']))
    datafile1 = find_timeseries("".join([os.getenv('SCRATCHDIR')+'/'+project1+'/PYSARTEST/']))
    print(format_args(['save_geodmod.py', '-f', datafile1, '-b', inps.SNWE, '-y', inps.latStep, '-x', inps.lonStep, '-startDate', inps.StartDate, '-endDate', inps.EndDate, '-outdir', inps.outdir]))
    completion_status = os.system(format_args(['save_geodmod.py', '-f', datafile1, '-b', inps.SNWE, '-y', inps.latStep, '-x', inps.lonStep, '-startDate', inps.StartDate, '-endDate', inps.EndDate, '-outdir', inps.outdir]))
    if completion_status == 1:
        print('error when runing save_geodmod.py')
        exit(0)
        
    os.chdir("".join([os.getenv('SCRATCHDIR')+'/'+project2+'/PYSARTEST/']))
    datafile2 = find_timeseries("".join([os.getenv('SCRATCHDIR')+'/'+project2+'/PYSARTEST/']))
    StartDate2,EndDate2 = find_date(datafile2,inps)
    print(format_args(['save_geodmod.py', '-f', datafile2, '-b', inps.SNWE, '-y', inps.latStep, '-x', inps.lonStep, '-startDate', StartDate2, '-endDate', EndDate2, '-outdir', inps.outdir]))
    completion_status = os.system(format_args(['save_geodmod.py', '-f', datafile2, '-b', inps.SNWE, '-y', inps.latStep, '-x', inps.lonStep, '-startDate', StartDate2, '-endDate', EndDate2, '-outdir', inps.outdir]))
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
   
def get_path_com(path):
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

def processdata(inps):
    #use geocode.py and save_roipac.py to process data"
    #  geocode timeseries**.h5 file and get the deformation field of two time periods
    #  and geocode ifgramStack.h5 file and get the coherence of two time periods
    #  and geocode geometryRadar.h5 file and get the dem  
    #atr_asc = inps.file[0]
    atr_asc = inps.File
    
    if os.path.exists("".join(inps.outdir))=='False':
        os.mkdir("".join(inps.outdir))
   
    # process cor and dem dataset
    corname='temporalCoherence.h5'
    cor_args = [corname, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    print("geocode.py", cor_args)
    args_str = format_args(cor_args)
    mintpy.geocode.main(args_str.split())
    
    demname='geometryRadar.h5'
    dem_args = [demname, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    print("geocode.py", dem_args)
    args_str = format_args(dem_args)
    mintpy.geocode.main(args_str.split())
   
    #unw file
    asc_args = [atr_asc, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    print("geocode.py", asc_args)
    args_str = format_args(asc_args)
    mintpy.geocode.main(args_str.split())
        
    #save dataset of unw cor and dem
    os.chdir("".join(inps.outdir))
    filename, extension = get_path_com("".join(atr_asc))[1:3]
    
    asct_args = ['geo_'+filename+extension, "".join([inps.StartDate,'_',inps.EndDate])]
    print("save_roipac.py", asct_args)
    asct_str = format_args(asct_args)
    mintpy.save_roipac.main(asct_str.split())     

    asct_args = ['geo_temporalCoherence.h5', '-o', "".join(['geo_',inps.StartDate,'_',inps.EndDate,'.cor'])]    
    print("save_roipac.py", asct_args)
    asct_str = format_args(asct_args)
    mintpy.save_roipac.main(asct_str.split())
    
    asct_args = ['geo_geometryRadar.h5', 'height', '-o', 'srtm.dem']
    print("save_roipac.py", asct_args)
    asct_str = format_args(asct_args)
    mintpy.save_roipac.main(asct_str.split())

######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    
    if inps.multiTrack and inps.multiTrack=='on':
        run_save_geodmod(inps)
    else:
        processdata(inps)
        # rename *.rsc1 to *.rsc
        outfile = format_args(['srtm.dem' + '.rsc'])
        write_rsc_file(inps,outfile,format_args(['srtm.dem' +'.rsc1']))
        os.remove(outfile)
        print('rename *.rsc1 to *.rsc')
        os.rename(format_args(['srtm.dem' +'.rsc1']),outfile)

######################################################################################
if __name__ == '__main__':
    main()
