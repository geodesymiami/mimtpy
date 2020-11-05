#!/usr/bin/env python3
##############################################################################################################
# Program is used for calculating horizontal and vertical component for overlapping area of two tracks       #
# Author: Lv Xiaoran                                                                                         #
# Created: July 2020                                                                                         #
############################################################################################################## 

import os
import argparse
import shutil
import numpy as np

import mintpy
from mintpy.utils import readfile, writefile, utils as ut

import mimtpy
import mimtpy.workflow 
from mimtpy.utils import multitrack_utilities as mu
######################################################################################
EXAMPLE = """example:
  
   overlap2horz_vert.py -m $SCRATCHDIR/BogdSenDT106/mintpy/velocity/velocity.h5 -mg $SCRATCHDIR/BogdSenDT106/mintpy/inputs/geo_geometry.h5 -s $SCRATCHDIR/BogdSenDT4/mintpy/velocity/velocity.h5 -sg $SCRATCHDIR/BogdSenDT4/mintpy/inputs/geo_geometry.h5 -o horizontal.h5 vertical.h5  

Reference paper:
[1] Gourmelen N , Amelung F , Casu F , et al. Mining-related ground deformation in Crescent Valley, Nevada: Implications for sparse GPS networks[J]. Geophysical Research Letters, 2007, 34(9):252-254
[2] Wright, T. J., B. E. Parsons, and Z. Lu (2004), Toward mapping surface deformation in three dimensions using InSAR, Geophys. Res. Lett., 31, L01607, doi:10.1029/2003GL018827.

Note:
Horizonatl component: Negative means west; Positive means east
Vertical component: Positive means up; Negative means down
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Calculate horizontal and vertical component for overlapping region',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('-m', '--master', nargs=1, dest='master', type=str, help='master track\n')

    parser.add_argument('-mg', '--mgeometry', dest='mgeometry', type=str, nargs=1,
                        help='geometry file of master data\n')
    
    parser.add_argument('-s', '--slave', nargs=1, dest='slave', type=str, help='slave track\n')

    parser.add_argument('-sg', '--sgeometry', dest='sgeometry', type=str, nargs=1,
                        help='geometry file of slave data\n')
    
    parser.add_argument('-o', '--outfile', dest='outfile', nargs=1, type=str,
                        help='outfile name')
    
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    m_atr = readfile.read_attribute(inps.master[0])
    s_atr = readfile.read_attribute(inps.slave[0])

    # check coordinates
    if any ('X_FIRST' not in i for i in [m_atr, s_atr]):
        raise Exception('Not all input files are geocoded!')

    # check spatial resolution
    if any (m_atr[i] != s_atr[i] for i in ['X_STEP', 'Y_STEP']):
        raise Exception('input files do not have the same spatial resolution')

    # check reference point
    # round to 3 decimal digits (~100m)
    ref_lalo_m = ['{:.1f}'.format(float(m_atr[i])) for i in ['REF_LAT', 'REF_LON']]
    ref_lalo_s = ['{:.1f}'.format(float(s_atr[i])) for i in ['REF_LAT', 'REF_LON']]

    #if ref_lalo_m != ref_lalo_s:
    #    raise Exception('input fils do not have the same reference point from REF_LAT/REF_LON value')

    return inps

def A_design(m_inc,s_inc):
    """parameters to calculate horizontal and vertical component"""
    A_horz_denominator = np.sin(m_inc) - (np.sin(s_inc) * (np.cos(m_inc) / np.cos(s_inc)))
    A_horz_numerator = np.cos(m_inc) / np.cos(s_inc)

    A_vert_denominator = np.cos(m_inc) - (np.cos(s_inc) * (np.sin(m_inc) / np.sin(s_inc)))
    A_vert_numerator = np.sin(m_inc) / np.sin(s_inc)

    return A_horz_denominator, A_horz_numerator, A_vert_denominator, A_vert_numerator

    

def latlon_overlay(m_atr, s_atr):
    """get lat/lon range of overlapping area
    lat0,lon0- starting latitude/longitude (first row/column)
    lat1,lon1- ending latitude/longitude (last row/column)
    """
    # calculate overlay region
    # lat lon range of master and slave
    m_lat0, m_lon0, m_lat1, m_lon1 = mimtpy.track_offset.get_bounding_box(atr=m_atr)
    
    s_lat0, s_lon0, s_lat1, s_lon1 = mimtpy.track_offset.get_bounding_box(atr=s_atr)

    #get lalo of overlay region
    over_lat0 = min(m_lat0,s_lat0)
    over_lon0 = max(m_lon0,s_lon0)
    over_lat1 = max(m_lat1,s_lat1)
    over_lon1 = min(m_lon1,s_lon1)    
   
    # get row/colm number of overlay region
    overlay_rows,overlay_colms = mimtpy.track_offset.calculate_roco_overlay(lat0=over_lat0,lon0=over_lon0,lat1=over_lat1,lon1=over_lon1,lat_step=float(m_atr['Y_STEP']),lon_step=float(m_atr['X_STEP']))
    
    #get the row/column position of overlay region
    m_row0, m_colm0 = mimtpy.track_offset.calculate_rc(lat0=m_lat0,lon0=m_lon0,lat_n=over_lat0,lon_n=over_lon0,lat_step=float(m_atr['Y_STEP']),lon_step=float(m_atr['X_STEP']))
    s_row0, s_colm0 = mimtpy.track_offset.calculate_rc(lat0=s_lat0,lon0=s_lon0,lat_n=over_lat0,lon_n=over_lon0,lat_step=float(s_atr['Y_STEP']),lon_step=float(s_atr['X_STEP']))
    
    return over_lat0, over_lon0, overlay_rows, overlay_colms, m_row0, m_colm0, s_row0, s_colm0
 
def overlap2horz_vert(inps):
    """calculate horizontal and vertical component for overlapping region of two tracks"""
    # master data
    m_data, m_atr = readfile.read(inps.master[0])
    
    # master geometry data 
    m_geometry = inps.mgeometry[0]
    print('processing master geometry data {}'.format(m_geometry))
    m_inc_angle = readfile.read(m_geometry, datasetName='incidenceAngle')[0]
    m_azi_angle = readfile.read(m_geometry, datasetName='azimuthAngle')[0]
    # heading angle
    m_head_angle = ut.azimuth2heading_angle(m_azi_angle)

    # slave data
    s_data, s_atr = readfile.read(inps.slave[0])

    # slave geometry data
    s_geometry = inps.sgeometry[0]
    print('processing slave geometry data {}'.format(s_geometry))
    s_inc_angle = readfile.read(s_geometry, datasetName='incidenceAngle')[0]
    s_azi_angle = readfile.read(s_geometry, datasetName='azimuthAngle')[0] 
    # heading angle
    s_head_angle = ut.azimuth2heading_angle(s_azi_angle)
    
    # calculate position of overlapping region
    over_lat0, over_lon0, overlay_rows, overlay_colms, m_row0, m_colm0, s_row0, s_colm0 = latlon_overlay(m_atr, s_atr)    
 
    # get master and slave overlay region data
    m_data_overlay = m_data[m_row0 : m_row0 + overlay_rows,m_colm0 : m_colm0 + overlay_colms]
    s_data_overlay = s_data[s_row0 : s_row0 + overlay_rows,s_colm0 : s_colm0 + overlay_colms]

    m_inc_angle_overlay = m_inc_angle[m_row0 : m_row0 + overlay_rows,m_colm0 : m_colm0 + overlay_colms] 
    m_inc_angle_overlay *= np.pi/180
    
    s_inc_angle_overlay = s_inc_angle[s_row0 : s_row0 + overlay_rows,s_colm0 : s_colm0 + overlay_colms]
    s_inc_angle_overlay *= np.pi/180

    # calculate horizontal and vertical component
    A_horz_denominator, A_horz_numerator, A_vert_denominator, A_vert_numerator = A_design(m_inc_angle_overlay, s_inc_angle_overlay)
    data_horz = (m_data_overlay - s_data_overlay * A_horz_numerator) / A_horz_denominator
    data_vert = (m_data_overlay - s_data_overlay * A_vert_numerator) / A_vert_denominator

    overlay_atr = dict()
    overlay_atr['LENGTH'] = overlay_rows
    overlay_atr['WIDTH'] = overlay_colms
    overlay_atr['X_FIRST'] = over_lon0
    overlay_atr['Y_FIRST'] = over_lat0 
    overlay_atr['X_UNIT'] = m_atr['X_UNIT']
    overlay_atr['Y_UNIT'] = m_atr['Y_UNIT']
    overlay_atr['X_STEP'] = m_atr['X_STEP']
    overlay_atr['Y_STEP'] = m_atr['Y_STEP']
    overlay_atr['FILE_TYPE'] = m_atr['FILE_TYPE']
    overlay_atr['UNIT'] = m_atr['UNIT']
    overlay_atr['WAVELENGTH'] = m_atr['WAVELENGTH']

    # write horzontal and vertical component
    horz_name = inps.outfile[0] + '_hz.h5' 
    writefile.write(data_horz, out_file=horz_name, metadata=overlay_atr)

    vert_name = inps.outfile[0] + '_up.h5'
    writefile.write(data_vert, out_file=vert_name, metadata=overlay_atr) 

    return
######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)   
    
    overlap2horz_vert(inps)
######################################################################################
if __name__ == '__main__':
    main()
