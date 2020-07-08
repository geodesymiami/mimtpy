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
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import geopandas
import numpy as np

from mintpy.objects import timeseries
from datetime import datetime as dt
import mintpy.objects.gps as gps
from mintpy.objects.gps import GPS
from mintpy.utils import utils as ut


min_figsize_single = 6.0       # default min size in inch, for single plot
max_figsize_single = 10.0      # default min size in inch, for single plot
max_figsize_height = 8.0       # max figure size in vertical direction in inch

######################################################################################
EXAMPLE = """example:
   
    # only plot faults and reference points
    plot_geotiff.py geotiff_file --shpdir /data/lxrtest/Balochistan/shp/ --fault multifault.shp nearbyfault.shp --fcolor b r --fstyle s d --refpoi refpoi_DT.shp --outfile velocity --outdir /data/lxrtest/BalochistanSenAT/shp/ 
    
    plot_geotiff.py geotiff_file --shpdir /data/lxrtest/Balochistan/shp/ --fault multifault.shp nearbyfault.shp --fcolor o m --fstyle s s --refpoi refpoi_DT.shp --vlim -0.02 0.02 --outfile velocity --outdir /data/lxrtest/BalochistanSenAT/shp/ 
   
    # plot shape file and gps velocity vector field
    # plot single gps site with extracting incidence and heading angle from geometreRadar 
    plot_geotiff.py geotiff_file --shpdir /data/lxrtest/Balochistan/shp/ --fault multifault.shp nearbyfault.shp --fcolor o m --fstyle s s --refpoi refpoi_DT.shp --gps_putton --gps_name GV05 --geometry $SCRATCHDIR/BalochistanSenAT115/mintpy/inputs/geometryRadar.h5 --gpsdir /data/lxr/insarlab/GPS/ --vlim -0.02 0.02 --outfile velocity --outdir /data/lxrtest/BalochistanSenAT/shp/ 
    
    # plot single gps site using typical incidence and heading angle for Seitinel1
    plot_geotiff.py geotiff_file --shpdir /data/lxrtest/Balochistan/shp/ --fault multifault.shp nearbyfault.shp --fcolor o m --fstyle s s --refpoi refpoi_DT.shp --gps_putton --gps_name GV05 --gpsdir /data/lxr/insarlab/GPS/ --vlim -0.02 0.02 --outfile velocity --outdir /data/lxrtest/BalochistanSenAT/shp/ 
    
    # plot multiple gps sites with extracting incidence and heading angle from geometreRadar 
    plot_geotiff.py geotiff_file --shpdir /data/lxrtest/Balochistan/shp/ --fault multifault.shp nearbyfault.shp --fcolor o m --fstyle s s --refpoi refpoi_DT.shp --gps_putton --search --bbox 36 40 36 39 --geometry $SCRATCHDIR/BalochistanSenAT115/mintpy/inputs/geometryRadar.h5 --gpsdir /data/lxr/insarlab/GPS/ --vlim -0.02 0.02 --outfile velocity --outdir /data/lxrtest/BalochistanSenAT/shp/ 
    
    # plot multiple gps sites with extracting incidence and heading angle from geometreRadar 
    plot_geotiff.py geotiff_file --shpdir /data/lxrtest/Balochistan/shp/ --fault multifault.shp nearbyfault.shp --fcolor o m --fstyle s s --refpoi refpoi_DT.shp --gps_putton --search --bbox 36 40 36 39 --gpsdir /data/lxr/insarlab/GPS/ --vlim -0.02 0.02 --outfile velocity --outdir /data/lxrtest/BalochistanSenAT/shp/ 
"""

fcolor_table = {'b':'black','r':'red','y':'yellow','o':'orange','m':'megenta'}
flinestyle_table = {'s':'solid','d':'dashed'}
#######################################################################################
def create_parser():
    parser = argparse.ArgumentParser(description='plot *.tiff',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('input_geotiff', nargs=1, type=str, help='geotiff file. \n')

    shp_option = parser.add_argument_group(title='option for shape files (faults and reference point).')   
 
    shp_option.add_argument('--shpdir', dest='shpdir', nargs=1, type=str, help='directory of shp files. \n')

    shp_option.add_argument('--fault', dest='fault', nargs='+', type=str, help='shp files of fault. \n')
    
    shp_option.add_argument('--fcolor', dest='fcolor', nargs='+', type=str, help='define color for faults:'
                                                              'b -> black; r -> red; y -> yellow;'
                                                              'o -> orange; m-> magenta \n')
    shp_option.add_argument('--fstyle', dest='fstyle', nargs='+', type=str, help='define line stype for faults:'
                                                              's -> solid; d -> dashed \n')
    shp_option.add_argument('--refpoi', dest='refpoi', nargs=1, type=str, help='shp file of reference point. \n')
   
    gps_option = parser.add_argument_group(title='options for plot gps velocity vector field')
    
    gps_option.add_argument('--gps_putton', action='store_true', default=False, help='whether plot gps velocity vector field. \n')

    gps_option.add_argument('--search', action='store_true', default=False, help='whether search gps sites based on the given SNWE. \n')
   
    gps_option.add_argument('--bbox', dest='SNWE', type=float, nargs=4, metavar=('S', 'N','W','E'),
                            help='Bounding box of area to be geocoded.\n'+
                            'Include the uppler left corner of the first pixel' + 
                            'and the lower right corner of the last pixel')
    gps_option.add_argument('--gps_name', dest='gps_name', nargs='?', type=str, help='gps site name. \n')

    gps_option.add_argument('--geometry', dest='geometry', nargs='?', type=str, help='geometry data. \n')

    gps_option.add_argument('--gpsdir', dest='gpsdir',nargs=1, type=str, help='directory to store gps data. \n')

    parser.add_argument('--vlim', dest = 'vlim', nargs=2, type=float, help='max and min value for displacement to display. The unit is m. \n')   
 
    parser.add_argument('--outdir',dest='outdir',nargs=1, type = str, help='output directory.\n')

    parser.add_argument('--outfile',dest='outfile', nargs=1, type=str, help='output file name.\n')

    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps

def auto_figure_size(shape, disp_cbar=False, ratio=1.0):
    """Get auto figure size based on input data shape"""
    length, width = shape
    plot_shape = [width*1.25, length]
    if not disp_cbar:
        plot_shape = [width, length]
    fig_scale = min(min_figsize_single/min(plot_shape),
                    max_figsize_single/max(plot_shape),
                    max_figsize_height/plot_shape[1])
    fig_size = [i*fig_scale*ratio for i in plot_shape]
    return fig_size

def fcolor(fcolors, shp_faults):
    """define fault color"""
    fault_colors = dict()
    for fcolor, fault in zip(fcolors,shp_faults):
        color = fcolor_table[fcolor]
        fault_colors[fault] = color

    return fault_colors

def flinestyle(fstyles, shp_faults):
    """defind fault linestyle"""
    fault_linestyles = dict()
    for fstyle, fault in zip(fstyles, shp_faults):
        style = flinestyle_table[fstyle]
        fault_linestyles[fault] = style

    return fault_linestyles

def gps_multisites(gpsname, geometry, gps_dir):
    """obtain multiple gps sites velocity field"""
    gps_obj = GPS(site=gpsname, data_dir=gps_dir)
    gps_obj.open()
    if not geometry:
        dis_los = ut.enu2los(gps_obj.dis_e, gps_obj.dis_n, gps_obj.dis_u)
        dates = gps_obj.dates
        los_vel = dis_velocity(dates,dis_los) 
    else:
        dates, dis_los = gps_obj.read_gps_los_displacement(geometry)
        gps_obj.get_gps_los_velocity(geometry)
        los_vel = gps_obj.velocity
    
    e_vel = dis_velocity(dates, gps_obj.dis_e)
    n_vel = dis_velocity(dates, gps_obj.dis_n)
   
    return np.array([[los_vel, n_vel, e_vel]])

def gps_single(gpsname, geometry, gps_dir):
    """process single gps site""" 
    gps_obj = GPS(site=gpsname, data_dir=gps_dir)
    gps_obj.open()
    if not geometry:
        dis_los = ut.enu2los(gps_obj.dis_e, gps_obj.dis_n, gps_obj.dis_u)
        dates = gps_obj.dates
        los_vel = dis_velocity(dates,dis_los) 
    else:
        dates, dis_los = gps_obj.read_gps_los_displacement(geometry)
        gps_obj.get_gps_los_velocity(geometry)
        los_vel = gps_obj.velocity
    
    e_vel = dis_velocity(dates, gps_obj.dis_e)
    n_vel = dis_velocity(dates, gps_obj.dis_n)
  
    site_inve = np.array([[gps_obj.site, gps_obj.site_lat, gps_obj.site_lon, los_vel, n_vel, e_vel]])
    
    return site_inve

def dis_velocity(dates, dis):
    """calculate velocity"""
    date_list = [dt.strftime(i, '%Y%m%d') for i in dates]
    if len(date_list) > 2:
        A = timeseries.get_design_matrix4average_velocity(date_list)
        vel = np.dot(np.linalg.pinv(A), dis)[0]
    else:
        vel = np.nan
    
    return vel
 
def search_gpssites(inps):
    """search gps sites based on given SNWE"""
    site_names, site_lats, site_lons = gps.search_gps(inps.SNWE)
    sites_baseinfo = np.transpose(np.vstack((np.vstack((np.array([site_names]),np.array([site_lats]))),np.array([site_lons])))) 
    
    return sites_baseinfo

def plot_geotiff(inps, site_infovel=None):
    """read geotiff data"""

    shpdir = inps.shpdir[0]

    # read raster file
    raster = rasterio.open(inps.input_geotiff[0])
    
    # read shape files
    shp_refpoi = geopandas.read_file(shpdir + inps.refpoi[0])
    
    # read fault files and set color/linestyle
    shp_faults = dict()
    for fault in inps.fault:
        shp_faults[fault] = geopandas.read_file(shpdir + fault)
    fault_colors = fcolor(inps.fcolor, shp_faults)
    fault_linestyles = flinestyle(inps.fstyle, shp_faults)

    """plot geotiff"""
    raster_data = raster.read(1)
    shape = np.shape(raster_data)
    if inps.vlim == None: 
        raster_min = np.nanmin(raster_data)
        raster_max = np.nanmax(raster_data)
    else:
        raster_min = inps.vlim[0]
        raster_max = inps.vlim[1]
    
    cmap = plt.cm.jet
    figure_size = auto_figure_size(shape)
    fig,axes = plt.subplots(1,1,figsize = figure_size)
    ax1 = axes
    print('\n\n***************************************ploting raster data****************************************')
    rasterio.plot.show(raster,1,ax = ax1, cmap=cmap,vmin = raster_min, vmax = raster_max, alpha=0.8)
    print('\n\n***************************************ploting fault**********************************************')
    for fault in inps.fault:
        shp_faults[fault].plot(color=fault_colors[fault],ax=ax1, linestyle=fault_linestyles[fault], linewidth=1)
    #shp_fault.plot(color='black',ax = ax1,linestyle='dashed',linewidth=1)
    print('\n\n************************************ploting reference point****************************************')
    shp_refpoi.plot(color='black',ax=ax1,marker='s')
   
    # plot gps velocity vector field
    if site_infovel is not None:
        lon = list(site_infovel[:,2].astype(float))
        lat = list(site_infovel[:,1].astype(float))
        x_component = list(site_infovel[:,5].astype(float))
        y_component = list(site_infovel[:,4].astype(float))
    
        ax1.quiver(lon, lat, x_component, y_component, color='deepskyblue')

    ax1.tick_params(which='both', direction='in', labelsize=18, bottom=True, top=True, left=True, right=True)
   
    # method 1 for drawing colorbar
    cax = make_axes_locatable(ax1).append_axes("right", "3%", pad="3%")
    norm = mpl.colors.Normalize(vmin=raster_min*100, vmax=raster_max*100)
    cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
    cbar.ax.tick_params(labelsize=18)

    # method 2 for drawing colorbar
    #cax1 = fig.add_axes([0.92, 0.18, 0.02, 0.6])
    #sm1 = plt.cm.ScalarMappable(cmap=cmap)
    #sm1.set_array([])
    #sm1.set_clim(vmin = raster_min*100, vmax = raster_max*100)
    ##fig.colorbar(sm1,cax1, orientation = 'vertical', label = 'LOS velocity(m/year)')
    #cbar = fig.colorbar(sm1,cax1, orientation = 'vertical',format='%.2f')
    #cbar.ax.tick_params(labelsize=18)
    ##cbar.set_ticks(np.linspace(-20,20,5))
    font2 = {'family' : 'serif',
             'weight': 'normal',
             'size' : 20.}
    
    filetype = os.path.split(inps.input_geotiff[0])[-1].split('_')[0]
    print('File type is %s' % filetype)
    if filetype == 'velocity': 
        cbar.set_label('velocity [cm/year]', fontdict=font2)
    elif filetype == 'displacement':
        cbar.set_label('displacement [cm]', fontdict=font2)

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
    # whether plot gps velocity vector field
    if inps.gps_putton:
        geodata = inps.geometry
        gpsdir = inps.gpsdir[0]
        if inps.search:
            sites_baseinfo = search_gpssites(inps)
            sites_velocity = np.zeros([sites_baseinfo.shape[0],3])
            for gpsname, i in zip(sites_baseinfo[:,0],np.arange(sites_baseinfo.shape[0])):
                site_ve = gps_multisites(str(gpsname),geodata,gpsdir)
                sites_velocity[i,:] = site_ve
            site_infovel = np.concatenate((sites_baseinfo,sites_velocity), axis=1)
            #print(sites_inve)
        else:
            gpsname = inps.gps_name
            site_infovel = gps_single(gpsname, geodata, gpsdir)
            #print(site_inve)    
        plot_geotiff(inps,site_infovel)
    else:
        plot_geotiff(inps)
######################################################################################
if __name__ == '__main__':
    main()
