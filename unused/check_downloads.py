#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright(c) 2019, Lv Xiaoran                            #
# Author:  Lv Xiaoran                                      #
############################################################



import os
import glob
import sys
import argparse
import zipfile


EXAMPLE = """example:
  check_downloads.py  $TESTDATA_ISCE/project/SLC/
  check_downloads.py  $TESTDATA_ISCE/project/SLC/ --delete
"""

def create_parser():
    parser = argparse.ArgumentParser(description='delete broken zipfiles in $TESTDATA_ISCE/project/SLC/.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('inputdir', nargs=1, help='directory for download zipfiles')
    parser.add_argument('--delete', action='store_true', default=False, help='whether delete data.')

    return parser


def cmd_line_parse(iargs=None):

    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    
    return inps


def check_zipfiles(inps):
    """
    check download zipfiles and get the reports
    """
    inputdir = "".join(inps.inputdir)
    os.chdir(inputdir)
    filelist = glob.glob('*.zip')
    broken_list = []
    for file in filelist:
        try:
            zf = zipfile.ZipFile(file,'r')
        except zipfile.BadZipFile:
            broken_list.append(file)
    print('The broken zipfiles are:')
    for filename in broken_list:
        print(filename)
    return broken_list


def delete_zipfiles(inps,broken_list):
    """delete data"""
    inputdir = "".join(inps.inputdir)
    os.chdir(inputdir)
    for file in broken_list:
        real_path = os.path.realpath(file)
        os.remove(real_path)
    return


##############################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    broken_files=check_zipfiles(inps)

    if inps.delete:
        delete_zipfiles(inps,broken_files)

##########################################################################
if __name__ == '__main__':
    main(sys.argv[1:])

