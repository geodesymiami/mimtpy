#!/usr/bin/env python3
###########################################################################
# Program is used for calculating [X Y] of lower edge middle point        #
# Author: Lv Xiaoran                                                      #
# Created: Feb  2021                                                      #
###########################################################################

import os
import argparse
import numpy as np

######################################################################################
EXAMPLE = """example:
    dislocation_conversion.py --para_upper -58.92 -18.67 8 56 257 32 
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Prepare data for insamp software based on ROIPAC format',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('--para_upper', dest='upper', type=float, nargs=6, 
                        metavar=('x_up', 'y_up','depth_up','dip', 'strike', 'width'),
                        help='information of upper edge middle point and fault geometry')

    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps

######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)   


    d2r = np.pi / 180

    # input parameters
    x_up = float(inps.upper[0]) 
    y_up = float(inps.upper[1])
    depth_up = float(inps.upper[2])
    dip = float(inps.upper[3])
    strike = float(inps.upper[4])
    width = float(inps.upper[5])

    # calculate 3-D dip vector:(1row, 3colm)
    dip_vector = np.array([(np.cos(dip * d2r) * np.cos(strike * d2r)),(-np.cos(dip * d2r) * np.sin(strike * d2r)), np.sin(dip * d2r)])

    # calculate 3-D strike vector:(1row, 3colm)
    str_vector = np.array([np.sin(strike * d2r),np.cos(strike * d2r),0])

    # calculate center point of the lower edge
    position = np.array([x_up, y_up, depth_up]) + (width * dip_vector)

    print('x coordinate of lower edge middle point is {}km'.format(position[0]))
    print('y coordinate of lower edge middle point is {}km'.format(position[1]))
    print('depth of lower edge {}km'.format(position[2]))
    
    return
######################################################################################
if __name__ == '__main__':
    main()
