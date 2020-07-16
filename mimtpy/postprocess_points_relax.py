#!/usr/bin/env python3
######################################################################
# Generate displacement for certain time period based on RELAX       #
# Author: Lv Xiaoran Mar 2020                                        #
######################################################################

import numpy as np
import argparse
import os
import json
import math
from mimtpy.utils import multitrack_utilities
##############################################################
EXAMPLE = """examples:
   postprocess_points_relax_density.py ./ --lalo ../radar/pos.lalo --minvisco 1700 --start 0.6 --end 2.3 --Tstep 0.1 --outname vs_17002100 --outdir ../result_density/
"""  

def create_parser():
    parser = argparse.ArgumentParser(description='postprocess displacement for certain time period based on output of Relax',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('txtdir', nargs=1,
                        help='observation file\n')
    parser.add_argument('--lalo', dest='lalo', nargs=1,
                        help='lalo information for observation points.\n')
    parser.add_argument('--minvisco', nargs=1, dest='minvisco', type=float,
                        help='the smallest value of the viscosity (value = viscosity * 100)\n')
    parser.add_argument('--start', dest='stime', type=float, nargs=1,
                        help='start time of the studied postseismic time.\n')
    parser.add_argument('--end', dest='etime', type=float, nargs=1,
                        help='end time of the studied postseismic time.\n')
    parser.add_argument('--Tstep',dest='tstep',type=float, nargs=1, metavar=('Interval'),
                        help='time step.\n')
    parser.add_argument('--outname', dest='outname', nargs=1,help='outname.\n')
    parser.add_argument('--outdir', dest='outdir', nargs=1,
                        help='output dir\n')
    return parser

def cmd_line_parser(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    return inps

def lookup_table(minvisco, time_interval, maxwell_time):
    """locate the lines according to the stime and etime, based on the mintvisco value"""
    lines_number = time_interval / (maxwell_time / 10)

    # this condition need more instances to support, only now is correct.
    if str(lines_number)[2] != 0:
        lines_number = np.ceil(lines_number)
    else:
        lines_number = np.floor(lines_number)
    # Instances:when time_step = 0.1yr and minvisco is ~10^17
    #if minvisco == 1700: (time_step/(maxwell_time)=9.59)
    #     lines_number = 10
    #elif minvisco == 1725: (time_step/(maxwell_time)=5.41)
    #    lines_number = 6
    #elif minvisco == 1750: (time_step/(maxwell_time)=3.03)
    #    lines_number = 3
    #elif minvisco == 1775: (time_step/(maxwell_time)=1.71)
    #    lines_number = 2
    
    return lines_number

def integrate_points(inps):
    """load txt file"""
    # read lalo file
    lalo_file = "".join(inps.lalo)
    lon, lat = np.loadtxt(lalo_file, skiprows=0, dtype=np.float, usecols=(0,1),unpack=True)
    # located time period
    # shear modulu
    mu = 3.00E+10
    sec_to_yr = 3.20E+7
    maxwell_time = ( math.pow(10,inps.minvisco[0]/100) / mu ) / sec_to_yr
    print('maxwell time is %f' % maxwell_time)
    time_interval = inps.tstep[0]
    
    if time_interval > maxwell_time / 10:
        print('using lookup table')
        lines_number = lookup_table(inps.minvisco[0], time_interval, maxwell_time)
        stime_location = ( inps.stime[0] / time_interval ) * lines_number
        etime_location = ( inps.etime[0] / time_interval ) * lines_number 
        
    else:
        interval = time_interval
        stime_location = inps.stime[0] / interval
        etime_location = inps.etime[0] / interval
    
    # integrate displacement files
    file_dir = "".join(inps.txtdir)
    
    #search G**.txt files
    Gfiles = []
    path_list = os.listdir(file_dir)
    for Gfile in path_list:
        if os.path.splitext(Gfile)[1] == '.txt' and str.find(os.path.split(Gfile)[-1],'G') != -1:
            Gfiles.append(Gfile)
    Gfiles.sort()
    
    obser_points = np.empty(shape=[0,5],dtype=float)
    for Gfile in Gfiles:
        print('process %s file' % Gfile)
        cumulate_disp = load_txt(file_dir + Gfile,stime_location,etime_location)
        points_order = int(Gfile.split('.')[0][1:])
        lalo = np.array([lon[points_order-1],lat[points_order-1]])
        obser_point = np.array([lalo[0],lalo[1],cumulate_disp[0],cumulate_disp[1],cumulate_disp[2]])
        obser_points = np.append(obser_points,[obser_point],axis=0)
 
    outdir = "".join(inps.outdir)
    outname = "".join(inps.outname)
    write_json(obser_points,outname,outdir)
 
def load_txt(Gfile,stime_location,etime_location):
    """load txt file"""
    time, north, east, down = np.loadtxt(Gfile, skiprows=1, dtype=np.float, usecols=(0,1,2,3), unpack=True, encoding='utf-8')
    print('start line is %.0f' % stime_location)
    print('start time is %.2f' % time[int(stime_location+1)])
    print('end line is %.0f' % etime_location)
    print('end time is %.2f' % time[int(etime_location+1)])   
    
    # located time period 
    cumulate_disp_stime = np.array([north[int(stime_location+1) + 1],east[int(stime_location+1) + 1],down[int(stime_location+1) + 1]])
    cumulate_disp_etime = np.array([north[int(etime_location+1) + 1],east[int(etime_location+1) + 1],down[int(etime_location+1) + 1]])
 
    # displacement between this time period
    cumulate_disp = cumulate_disp_etime - cumulate_disp_stime
    
    return cumulate_disp

def write_json(obser_points,json_name,outdir):
    """write into json"""
    json_file = json_name + '.json'
    json_data = {"displacement":obser_points.tolist()}
    open(outdir + json_file,"w").write(json.dumps(json_data))

######################################################################
def main(iargs=None):
    inps = cmd_line_parser(iargs)
    integrate_points(inps)
#################################################################################################
if __name__ == '__main__':
    main()
