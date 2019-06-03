#! /usr/bin/env python3
"""This script downloads SAR data and deals with various errors produced by download clients
   Author: Falk Amelung
   Created:12/2018
"""
###############################################################################

import os
import sys
import argparse
import subprocess

from minsar.objects import message_rsmas
import mimt.utils.process_utilities as putils
from minsar.objects.dataset_template import Template
from mimt.utils.process_utilities import get_work_directory, get_project_name, create_or_update_template
import minsar.job_submission as js

###############################################################################
EXAMPLE = '''example:
  mimt.py $SAMPLESDIR/GalapagosSenDT128.template
'''


def command_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps

def create_parser():
    """ Creates command line argument parser object. """
    parser = argparse.ArgumentParser(description='Downloads SAR data using a variety of scripts',
                                     formatter_class=argparse.RawTextHelpFormatter, epilog=EXAMPLE)
    parser.add_argument('customTemplateFile', help='template file containing ssaraopt field')
    parser.add_argument( '--submit', dest='submit_flag', action='store_true', help='submits job')

    # parser.add_argument('customTemplateFile', help='template file containing ssaraopt field', nargs='?')
    return parser

###############################################################################

def main(iargs=None):
    """Calculate horz and vert."""

    inps = command_line_parse(iargs)

    inps.project_name = get_project_name(inps.customTemplateFile)
    inps.work_dir = get_work_directory(None, inps.project_name)
    if not os.path.isdir(inps.work_dir):
       os.mkdir(inps.work_dir)
    os.chdir(inps.work_dir)

    command = os.path.basename(__file__) + ' ' + iargs[0]
    message_rsmas.log(command)

    #########################################
    # Submit job
    #########################################
    if inps.submit_flag:
        job_file_name = 'download_rsmas'
        work_dir = os.getcwd()
        job_name = inps.customTemplateFile.split(os.sep)[-1].split('.')[0]
        wall_time = '24:00'

        js.submit_script(job_name, job_file_name, sys.argv[:], work_dir, wall_time)
        sys.exit(0)

    inps = putils.create_or_update_template(inps)

    if not os.path.isdir(inps.work_dir):
        os.makedirs(inps.work_dir)

    import pdb; pdb.set_trace()
    if bool(inps.template['horzvert']):
       generate_verthorz(inps.customTemplateFile)
       
    #if not os.path.isdir(slc_dir):
    #    os.makedirs(slc_dir)

    os.chdir(inps.work_dir)

###########################################################################################

if __name__ == '__main__':
    main(sys.argv[1:])
