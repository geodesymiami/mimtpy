#!/usr/bin/env python3
#################################################################
# Program is used for preparing data for kite software          #
# Author: Lv Xiaoran                                            #
# Created: Oct 2020                                             #
#################################################################

import os
import argparse
import numpy as np
import json
import matplotlib.pyplot as plt

import mintpy
from mintpy.utils import readfile
from mintpy.objects import RAMP_LIST

import mimtpy
import mimtpy.workflow
from mimtpy.utils import multitrack_utilities as mu

######################################################################################
EXAMPLE = """example:
  
  gridsearch_ramps_batch_relax.py KokoxiliModelData_AT143_ramped.h5 ./model_data mask.h5 --ramp -s linear --outdir ./
  
  gridsearch_ramps_batch_relax.py KokoxiliModelData_AT143_ramped.h5 ./model_data mask.h5 --outdir ./
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Prepare data for Kite software',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('Obs_file', nargs=1, type=str, help='Observed data')

    parser.add_argument('Model_dir', nargs=1, type=str, help='Modeled data dir')

    parser.add_argument('mask', nargs=1, type=str, help='Mask file')

    parser.add_argument('--ramp', dest='ramp', action='store_true', help='estimate ramps or not')

    parser.add_argument('-s', dest='surface_type', default='linear', choices=RAMP_LIST,
                        help='type of surface/ramp to calculated, linear by default')

    parser.add_argument('--outdir', dest='outdir', type=str, nargs=1,
                        help='output file dir')

    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps

def run_gridsearch_ramps(inps, model_file):
    """run_gridsearch_ramps"""
    observed_file = inps.Obs_file[0]
    ramp_type = inps.surface_type
    mask_file = inps.mask[0]


    if inps.ramp:
        scp_args = ['../' + observed_file, model_file, '../' + mask_file, '-s', ramp_type, '-o', model_file[:-3] + '_ramp_noise.h5', '--ramp_file', model_file[:-3] + '_ramps.h5', '--outdir', inps.outdir[0] + 'residual/']
        scp_args = mu.seperate_str_byspace(scp_args)

        print('gridsearch_ramps_relax.py', scp_args)
        print('run gridsearch_ramps_relax.py for model data: %s' % model_file)
        mimtpy.gridsearch_ramps_relax.main(scp_args.split())

    else:
        scp_args = ['../' + observed_file, model_file, '../' + mask_file, '--outdir', inps.outdir[0] + 'residual/']
        scp_args = mu.seperate_str_byspace(scp_args)

        print('gridsearch_ramps_relax.py', scp_args)
        print('run gridsearch_ramps_relax.py for model data: %s' % model_file)
        mimtpy.gridsearch_ramps_relax.main(scp_args.split())

    return 

def calculate_residual_RMSE(residual_file):
    """calculate RMSE for each residual data"""
    residual_data = readfile.read(residual_file)[0]

    # calcualte rms of the residual data
    rms = np.sqrt((np.nansum(residual_data * residual_data)) / (np.nansum(~np.isnan(residual_data))))
    return rms

def write_json(data, name):
    """write json file"""
    residual = {"rmse": data.tolist()}
    open(name + '.json', 'w').write(json.dumps(residual))

    return 

def plot_residual(residual):
    """plot residual"""
    print("****************************************\n")
    print("Start drawing results\n")
    x_num = len(np.unique(residual[:,0]))
    y_num = len(np.unique(residual[:,1]))
    x = np.linspace(np.min(residual[:,0]),np.max(residual[:,0]),x_num)
    y = np.linspace(np.min(residual[:,1]),np.max(residual[:,1]),y_num)
    value = residual[:,2].reshape(x_num,y_num)
    X,Y = np.meshgrid(x,y)
    im = plt.contourf(X,Y,value,alpha=0.75,cmap=plt.cm.jet)
    cbar=plt.colorbar(im)
    contour = plt.contour(X, Y, value, 8, colors = 'black')
    plt.clabel(contour, fontsize=6, colors='k',fmt='%.1f')
    #cbar.ax.tick_params(labelsize=10)
    cbar.set_label('rms')
    #cbar.set_ticks(np.linspace(0,1,10))
    # set the range of legend
    #cbar.set_ticks(np.linspace(0,1,50))
    plt.title("InSAR grid search") 
    plt.xlabel('para1')
    plt.ylabel('para2')
    # C = plt.contour(X, Y, value, 8, colors = 'black', linewidth = 0.5)
    #plt.clabel(C, inline = True, fontsize = 10)
    plt.savefig('gird_search_ramp.png', dpi=300, bbox_inches='tight')

######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)   
   
    model_dir = inps.Model_dir[0]

    model_files = []
    path_list = os.listdir(model_dir)
    for model_file in path_list:
        if model_file.find('.h5') != -1:
            model_files.append(model_file)

    model_files.sort()

    rms_summary = np.empty(shape=[0,3], dtype=float)

    os.chdir(model_dir)
    os.mkdir('residual')
    for model_file in model_files:
        print('\n--------------------------------------------------')
        print('process %s forward model data' % model_file)
        run_gridsearch_ramps(inps, model_file)
        if inps.ramp:
            residual_file = model_file[:-3] + '_residual.h5'
        else:
            residual_file = model_file[:-3] + '_residual_noramp.h5'
        
        rms = calculate_residual_RMSE('./residual/' + residual_file)

        # for synthetic data
        #para1 = float(model_file.split('_')[1][0:2])
        #para2 = float(model_file.split('_')[1][4:])

        # for RELAX
        para1 = float(model_file.split('_')[0])
        para2 = float(model_file.split('_')[1])
 
        rms_summary = np.append(rms_summary, [np.array([para1, para2, rms])], axis=0)
        print('------------------------------------------------------------')
    print(rms_summary)

    # store rms file
    file_name = inps.outdir[0] + 'residual/' + 'rmse'    
    write_json(rms_summary, file_name)

    # plot 
    os.chdir(inps.outdir[0] + 'residual/')
    plot_residual(rms_summary)

    print('finish!')
######################################################################################
if __name__ == '__main__':
    main()
