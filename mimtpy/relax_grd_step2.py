#!/usr/bin/env python3
##############################################################################
# Program is used for masking Relax grd file based on Mintpy product         #
# Author: Lv Xiaoran                                                         #
# Created: July 2020                                                         #
##############################################################################

import os
import argparse
import numpy as np
import netCDF4 as nc
import copy

import mintpy
from mintpy.utils import readfile, writefile, utils as ut
######################################################################################
EXAMPLE = """example:
  Note:
  The script has two functions:
  1.It can be used for masking relax grd file (three dirctions) based on mintpy product covering the same region.
  2.It can calculat LOS direction simulation data and generate residual between model and InSAR results.
  
  Prerequisites of using this scripts:
  1. using xyz2grd.sh to extract same coverage and same resolution (normally 0.002) of simulated data. 
  2. using grdsub.sh to calculate the cumulative simulated displacement between two time periods.

  Output:
  simulated displacement in three direction;
  simulated LOS displacement;
  residual between simulated data and InSAR data;  

  relax_grd_step2.py 063-058-relax-geo ../Mintpy/dis_2015_2020.h5 --outdir ./HDFEOS/
  
  relax_grd_step2.py 063-058-relax-geo ../Mintpy/dis_2015_2020.h5 --LOS -g $SCRATCHDIR/BogdSenDT/geometry/geo_geometry.h5 --output 063-058-relax --outdir ./HDFEOS/
  
  relax_grd_step2.py 063-058-relax-geo ../Mintpy/dis_2015_2020.h5 --LOS --inc_angle 34 --azi_angle -102 --output 063-058-relax --outdir ./HDFEOS/
  
  relax_grd_step2.py 063-058-relax-geo ../Mintpy/dis_2015_2020.h5 --residual -g $SCRATCHDIR/BogdSenDT/geometry/geo_geometry.h5 --output 063-058-relax --outdir ./HDFEOS/
"""
Direction=['east',
           'north',
           'up']

#######################################################################################
def create_parser():
    parser = argparse.ArgumentParser(description='masking relax grd file based on mintpy product covering the same region',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('grdfile', nargs=1, type=str, help='geocoded unw or h5 files to be converted\n')

    parser.add_argument('mintfile', type=str, nargs=1, help='mintpy file to provide mask\n')

    parser.add_argument('--los', action='store_true', default=False, help='whether generate LOS simulation data\n')
 
    parser.add_argument('--residual', action='store_true', default=False, help='whether calculat resicual between LOS simualtion and InSAR\n')

    parser.add_argument('-g', '--geometry', dest='geometry', nargs='?', type=str, help='geometry data containing incidence and azimuth\n')

    parser.add_argument('--inc_angle', dest='inc_angle', nargs='?', type=float, help='incidence angle\n')

    parser.add_argument('--azi_angle', dest='azi_angle', nargs='?', type=float, help='azimuth angle\n')    

    parser.add_argument('--output', dest='output', nargs=1, type=str, help='outfile name\n')
 
    parser.add_argument('--outdir', dest='outdir', nargs=1, type=str, help='output dir\n')
    
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps

#####################################################################################################
def relax_mask(inps):
    """mask relax grd file using mintpy product"""
    # mintpy data
    mint_file = inps.mintfile[0]
    print('processing mintpy data {}'.format(mint_file))
    mint_data, atr = readfile.read(mint_file)
    
    # process relax grd file 
    grd_dict = dict()
    for direction in Direction:
        # read relax grd file
        relax_file = inps.grdfile[0] + '-' + direction + '.grd'
        print("processing relax grd file: {}".format(relax_file))
        grd_file = nc.Dataset(relax_file, 'r')
        grd_dataset = grd_file.variables['z'][:,:]
        grd_data = np.flipud(grd_dataset.data)
  
        # mask grd file
        grd_data[np.isnan(mint_data)] = np.nan  
        grd_dict[direction] = grd_data

        # write masked grd file attribute table
        grd_lon_dataset = grd_file.variables['x'][:]
        grd_lon = grd_lon_dataset.data
        grd_lat_dataset = grd_file.variables['y'][:]
        grd_lat = grd_lat_dataset.data

        grd_atr = dict()
        grd_atr["X_FIRST"] = str(grd_lon[0])
        grd_atr["X_STEP"] = atr["X_STEP"]
        grd_atr["X_UNIT"] = atr["X_UNIT"]
        grd_atr["Y_STEP"] = atr["Y_STEP"]
        grd_atr["Y_FIRST"] = str(grd_lat[-1])
        grd_atr["Y_UNIT"] = atr["Y_UNIT"]
        grd_atr["LENGTH"] = atr["LENGTH"]
        grd_atr["WIDTH"]  = atr["WIDTH"]
        grd_atr["PROCESSOR"] = 'RELAX'
        grd_atr["FILE_TYPE"] = '.unw'
        grd_atr["unit"] = 'meter'

        # write masked relax grd file to h5 file
        outfile_h5 = inps.outdir[0] + inps.grdfile[0] + '-' + direction + '.h5'
        writefile.write(datasetDict=grd_data, out_file=outfile_h5, metadata=grd_atr)   

    return grd_dict, grd_atr, mint_data

def LOS_calculation(grd_dict, grd_atr, inps):
    """calcualte simulated LOS deformation"""
    if not inps.geometry:
        inc_angle = inps.inc_angle

        head_angle = ut.azimuth2heading_angle(inps.azi_angle)
        if head_angle < 0.:
            head_angle += 360.
    else:
        # geometry data 
        geometryRadar = inps.geometry[0]
        print('processing geometry data {}'.format(geometryRadar))
        inc_angle = readfile.read(geometryRadar, datasetName='incidenceAngle')[0]
        azi_angle = readfile.read(geometryRadar, datasetName='azimuthAngle')[0]

        # heading angle
        head_angle = ut.azimuth2heading_angle(azi_angle) 
        head_angle[head_angle<0.]+= 360.
    
    head_angle *= np.pi/180.

    # incidence angle
    inc = copy.deepcopy(inc_angle)
    inc *= np.pi/180.

    # construct design matrix
    A_up = np.cos(inc)
    A_east = - np.sin(inc) * np.cos(head_angle)
    A_north = np.sin(inc) * np.sin(head_angle)
    
    # three component
    east_disp = grd_dict['east']
    north_disp = grd_dict['north']
    up_disp = grd_dict['up']

    # LOS simulated results. Note:the unit of RELAX displacement is meter.
    #Project displacement from LOS to Horizontal and Vertical components
    #    math for 3D: cos(theta)*Uz - cos(alpha)*sin(theta)*Ux + sin(alpha)*sin(theta)*Uy = Ulos
    #    math for 2D: cos(theta)*Uv - sin(alpha-az)*sin(theta)*Uh = Ulos   #Uh_perp = 0.0
    los_sim = up_disp * A_up + north_disp * A_north + east_disp * A_east
    
    # write masked relax grd file to h5 file
    outfile_h5 = inps.outdir[0] + inps.grdfile[0] + '-' + 'LOS.h5'
    writefile.write(datasetDict=los_sim, out_file=outfile_h5, metadata=grd_atr)   

    return los_sim

def residual_cal(inps, los_sim, los_insar, grd_atr):
    """calculate residual"""
    res = los_insar - los_sim
    
    # write masked relax grd file to h5 file
    outfile_h5 = inps.outdir[0] + inps.grdfile[0] + '-' + 'insar_sim_res.h5'
    writefile.write(datasetDict=res, out_file=outfile_h5, metadata=grd_atr)

    return   

######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)   
    # generate three component of relax simulation data 
    grd_dict, grd_atr, los_insar = relax_mask(inps)

    if inps.residual:
        los_sim = LOS_calculation(grd_dict, grd_atr, inps)
        residual_cal(inps, los_sim, los_insar, grd_atr)
    
    if inps.los:
        los_sim = LOS_calculation(grd_dict, grd_atr, inps)
        
######################################################################################
if __name__ == '__main__':
    main()
