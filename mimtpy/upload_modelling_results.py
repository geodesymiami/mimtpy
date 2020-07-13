#!/usr/bin/env python3
########################
# Author:  Falk Amelung
#######################

import os
import subprocess
import sys
import glob
import time
import shutil
import argparse

#######################################################################################
EXAMPLE = """example:
    upload_modelling_product.py --model relax --project Bogd
"""
model_dir = {'relax': os.getenv('RELAX_DIR')}
#######################################################################################
def create_parser():
    parser = argparse.ArgumentParser(description='upload modelling results to Jetstream',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    model_option = parser.add_argument_group(title='model options.')
    
    model_option.add_argument('--model', dest='model', nargs=1, type=str, help='Model software name\n')
    
    model_option.add_argument('--project', dest='project', nargs=1, type=str, help='Project name. \n')
 
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps

##############################################################################
def main(iargs=None):

    inps = cmd_line_parse(iargs)

    model = inps.model[0]  
  
    project_name = inps.project[0]
 
    work_dir = model_dir[model] + '/examples/' + project_name + '/'  

    os.chdir(work_dir)

    # get DATA_SERVER and return if it does not exist

    DATA_SERVER = 'centos@129.114.104.223'

    REMOTE_DIR = '/data/HDF5EOS/'
    
    destination = DATA_SERVER + ':' + REMOTE_DIR + 'RELAX/'

    pattern = '/' + project_name + '_vs_*'
    
    #rsync_list = [
    #        '/mintpy/pic',
    #        '/mintpy/*.he5',
    #        '/mintpy/inputs'
    #        ]

    command = 'ssh ' + DATA_SERVER + ' mkdir -p ' + REMOTE_DIR + 'RELAX/' + project_name + '/'
    print (command)
    status = subprocess.Popen(command, shell=True).wait()
    if status is not 0:
         raise Exception('ERROR in upload_modelling_products.py')

    #for pattern in rsync_list:
        
    command = 'scp -r ' + work_dir + pattern + ' ' + destination + project_name
    print (command)
    status = subprocess.Popen(command, shell=True).wait()
    if status is not 0:
        raise Exception('ERROR in upload_modelling_products.py')

    print ('\nAdjusting permissions:')
    command = 'ssh ' + DATA_SERVER + ' chmod -R u=rwX,go=rX ' + REMOTE_DIR + 'RELAX/' + project_name 
    print (command)
    status = subprocess.Popen(command, shell=True).wait()
    if status is not 0:
         raise Exception('ERROR in upload_modelling_products.py')

    return None

if __name__ == "__main__":
    main()
