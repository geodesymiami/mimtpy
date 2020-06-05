#!/usr/bin/env python3

######################################################################################################
# Program is used for plot geotiff figures                                                           #
# Author: Lv Xiaoran                                                                                 #
# Created: February  2020                                                                            #
######################################################################################################
import argparse
import os
import rasterio
import rasterio.plot
import matplotlib.pyplot as plt
import geopandas
import numpy as np

######################################################################################
EXAMPLE = """example:
    plot_geotiff.py geotiff_file --fault /data/lxrtest/Balochistan/shp/multifault.shp --refpoi /data/lxrtest/Balochistan/shp/refpoi_DT.shp --outfile velocity --outdir /data/lxrtest/BalochistanSenAT/shp/ 
    
    plot_geotiff.py geotiff_file --fault /data/lxrtest/Balochistan/shp/multifault.shp --refpoi /data/lxrtest/Balochistan/shp/refpoi_DT.shp --vlim -0.02 0.02 --outfile velocity --outdir /data/lxrtest/BalochistanSenAT/shp/ 
"""

def create_parser():
    parser = argparse.ArgumentParser(description='plot *.tiff',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('input_geotiff', nargs=1, type=str, help='geotiff file. \n')

    parser.add_argument('--fault', nargs=1, type=str, help='shp file of fault. \n')
    
    parser.add_argument('--refpoi', nargs=1, type=str, help='shp file of reference point. \n')
   
    parser.add_argument('--vlim', dest = 'vlim', nargs=2, type=float, help='max and min value of colorbar. \n')   
 
    parser.add_argument('--outdir',dest='outdir',nargs=1, type = str, help='output directory.\n')

    parser.add_argument('--outfile',dest='outfile', nargs=1, type=str, help='output file name.\n')
    
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps

def plot_geotiff(inps):
    """read geotiff data"""
    raster = rasterio.open(inps.input_geotiff[0])
    shp_fault = geopandas.read_file(inps.fault[0])
    shp_refpoi = geopandas.read_file(inps.refpoi[0])

    """plot geotiff"""
    raster_data = raster.read(1)
    if inps.vlim == None: 
        raster_min = np.nanmin(raster_data)
        raster_max = np.nanmax(raster_data)
    else:
        raster_min = inps.vlim[0]
        raster_max = inps.vlim[1]
    
    cmap = plt.cm.jet
    figure_size = [8.0,8.0]
    fig,axes = plt.subplots(1,1,figsize = figure_size)
    ax1 = axes
    print('\n\n***************************************ploting raster data****************************************')
    rasterio.plot.show(raster,1,ax = ax1, cmap=cmap,vmin = raster_min, vmax = raster_max, alpha=0.8)
    print('\n\n***************************************ploting fault**********************************************')
    shp_fault.plot(color='black',ax = ax1,linestyle='dashed',linewidth=1)
    print('\n\n************************************ploting reference point****************************************')
    shp_refpoi.plot(color='black',ax=ax1,marker='s')
    ax1.tick_params(which='both', direction='in', labelsize=18, bottom=True, top=True, left=True, right=True)
    
    cax1 = fig.add_axes([0.92, 0.18, 0.02, 0.6])
    sm1 = plt.cm.ScalarMappable(cmap=cmap)
    sm1.set_array([])
    sm1.set_clim(vmin = raster_min*100, vmax = raster_max*100)
    #fig.colorbar(sm1,cax1, orientation = 'vertical', label = 'LOS velocity(m/year)')
    cb = fig.colorbar(sm1,cax1, orientation = 'vertical',format='%.2f')
    cb.ax.tick_params(labelsize=18)
    #cb.set_ticks(np.linspace(-20,20,5))
    font2 = {'family' : 'serif',
             'weight': 'normal',
             'size' : 20.}
    filetype = os.path.split(inps.input_geotiff[0])[-1].split('_')[0]
    print('File type is %s' % filetype)
    if filetype == 'velocity': 
        cb.set_label('velocity [cm/year]', fontdict=font2)
    elif filetype == 'displacement':
        cb.set_label('displacement [cm/year]', fontdict=font2)

    #设置横纵坐标刻度值大小及字体
    font1 = {'family' : 'serif',
             'weight': 'normal',
             'size' : 18.}
    ax1.set_xlabel('Longitude [degree]',font1)
    ax1.set_ylabel('Latitude [degree]',font1)
    labels = ax1.get_xticklabels() + ax1.get_yticklabels()
    [label.set_fontname('serif') for label in labels]


    #save figure
    fig_name = inps.outfile[0] + '.png'
    fig_output = inps.outdir[0] + fig_name
    fig.savefig(fig_output, dpi=300, bbox_inches='tight')
######################################################################################
def main(iagrs=None):
    inps = cmd_line_parse(iagrs)
    plot_geotiff(inps)    
######################################################################################
if __name__ == '__main__':
    main()
