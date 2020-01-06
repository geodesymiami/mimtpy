#!/usr/bin/env python3
######################################################################################################
# Program is used for generating simulated LOS displacement based on MintPy and PSGRN/PSCMP results  #
# Author: Lv Xiaoran                                                                                 #
# Created: January  2020                                                                             #
######################################################################################################

import os
import sys
import argparse
import numpy as np
import json
from mintpy.utils import readfile
######################################################################################
EXAMPLE = """example:
    subtract_simresults_to_Insar_region.py $MODELOUT/pscmp/Bogd_test2/model4/poseis-10y_20y_east_west.json $SCRATCHDIR/BogdSenDT106/mintpy/S1_IW123_106_0441_0447_20161013_20191215.he5 --outdir $MODELOUT/pscmp/Bogd_test2/model4/ 
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Generate *.tiff files based on psgrn/pscmp *.dat files',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('input_json', nargs=1, help='directory stored *.json files. \n')

    parser.add_argument('input_HDFEOS', nargs=1, help='directory stored *.he5 files. \n')
    
    parser.add_argument('--outdir',dest='outdir',nargs='?',default=os.getenv('MODELOUT'), help='output directory')

    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps    
    
######################################################################################
def main(iargs=None):
    """subtract corresponding region from simulated displacement accroding InSAR results"""
    inps = cmd_line_parse(iargs)
    # read 2d simulated displacement data
    disp_sim = json.loads(open("".join(inps.input_json)).read())
    # attribuate
    ul_lon_sim = disp_sim["ul_lon"]
    ul_lat_sim = disp_sim["ul_lat"]
    ll_lon_sim = ul_lon_sim
    ll_lat_sim = ul_lat_sim + disp_sim["lat_step"] * disp_sim["rows"]
    # displacement data
    disp_data = np.array(disp_sim["disp_data"])
    
    #read HDFEOS file
    # metadata
    atr = readfile.read_attribute("".join(inps.input_HDFEOS))
    # attribute
    ul_lat_h5 = float(atr['Y_FIRST'])
    ul_lon_h5 = float(atr['X_FIRST'])
    samples_h5 = int(atr["WIDTH"])
    rows_h5 = int(atr["LENGTH"])
    ll_lon_h5 = ul_lon_h5
    ll_lat_h5 = ul_lat_h5 + float(atr["Y_STEP"][1:]) * rows_h5
    # data
    incidence = readfile.read("".join(inps.input_HDFEOS), datasetName='/HDFEOS/GRIDS/timeseries/geometry/incidenceAngle')[0]
    azimuth = readfile.read("".join(inps.input_HDFEOS), datasetName='/HDFEOS/GRIDS/timeseries/geometry/azimuthAngle')[0]
    
    #locate Insar position in Simulated field and subtract corresponding region from simulated displacement accroding InSAR results 
    
    row_locate = int(np.floor((ll_lat_h5 - ll_lat_sim) / np.abs(disp_sim["lat_step"])))
    sample_locate = int(np.floor((ll_lon_h5 - ll_lon_sim) / np.abs(disp_sim["lon_step"])))
    
    disp_sub = disp_data[row_locate:row_locate + rows_h5, sample_locate:sample_locate + samples_h5]
    position = np.isnan(incidence)
    disp_sub[position] = np.nan  
    
    #write subtracted simulated displacement
    displacement = {"ul_lon": ul_lon_h5, "ul_lat": ul_lat_h5, "lat_step": atr["Y_STEP"], "lon_step": atr["X_STEP"], "samples": samples_h5, "rows": rows_h5, "disp_data": disp_sub.tolist()}
    outname = os.path.basename(inps.input_json).split('.')[0] + 'subtract.json'
    open(outdir + '/' + outname, "w").write(json.dumps(displacement))
    
    
######################################################################################
if __name__ == '__main__':
    main()
