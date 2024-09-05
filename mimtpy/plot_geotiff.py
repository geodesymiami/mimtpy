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
import matplotlib.colors as colors
import geopandas
import numpy as np

from mintpy.objects import timeseries
from datetime import datetime as dt
import mintpy.objects.gnss as gnss
from mintpy.objects.gnss import GNSS_UNR
from mintpy.utils import utils as ut

import mimtpy.objects.cgps as cgps
from mimtpy.objects.cgps import CGPS

min_figsize_single = 6.0       # default min size in inch, for single plot
max_figsize_single = 10.0      # default min size in inch, for single plot
max_figsize_height = 8.0       # max figure size in vertical direction in inch

######################################################################################
EXAMPLE = """example:
   
    # only plot faults and reference points
    plot_geotiff.py geotiff_file --shpdir /data/lxrtest/Balochistan/shp/ --fault multifault.shp nearbyfault.shp --fcolor b r --fstyle s d --refpoi refpoi_DT.shp --outfile velocity --outdir /data/lxrtest/BalochistanSenAT/shp/ 
    
    plot_geotiff.py geotiff_file --shpdir /data/lxrtest/Balochistan/shp/ --fault multifault.shp nearbyfault.shp --fcolor o m --fstyle s s --refpoi refpoi_DT.shp --vlim -0.02 0.02 --outfile velocity --outdir /data/lxrtest/BalochistanSenAT/shp/ 
   
    # plot data using nonlinear colorbar
    plot_geotiff.py geotiff_file --shpdir /data/lxrtest/Balochistan/shp/ --fault multifault.shp nearbyfault.shp --fcolor b r --fstyle s d --refpoi refpoi_DT.shp --nonlinear --outfile velocity --outdir /data/lxrtest/BalochistanSenAT/shp/ 
    

    # plot shape file and gps velocity vector field
    # plot single gps site with extracting incidence and heading angle from geometreRadar 
    plot_geotiff.py geotiff_file --shpdir /data/lxrtest/Balochistan/shp/ --fault multifault.shp nearbyfault.shp --fcolor o m --fstyle s s --refpoi refpoi_DT.shp --gps_putton --type NGPS --gps_name GV05 --geometry $SCRATCHDIR/BalochistanSenAT115/mintpy/inputs/geometryRadar.h5 --gpsdir /data/lxr/insarlab/GPS/ --scale 1 --scale_key 0.01 --vlim -0.02 0.02 --outfile velocity --outdir /data/lxrtest/BalochistanSenAT/shp/ 
    
    # plot single gps site using typical incidence and heading angle for Seitinel1
    plot_geotiff.py geotiff_file --shpdir /data/lxrtest/Balochistan/shp/ --fault multifault.shp nearbyfault.shp --fcolor o m --fstyle s s --refpoi refpoi_DT.shp --gps_putton --type CGPS --gps_name AHAQ --scale 1 --scale_key 0.01 --vlim -0.02 0.02 --outfile velocity --outdir /data/lxrtest/BalochistanSenAT/shp/ 
    
    # plot multiple gps sites with extracting incidence and heading angle from geometreRadar 
    plot_geotiff.py geotiff_file --shpdir /data/lxrtest/Balochistan/shp/ --fault multifault.shp nearbyfault.shp --fcolor o m --fstyle s s --refpoi refpoi_DT.shp --gps_putton --type NGPS --search --bbox 36 40 36 39 --geometry $SCRATCHDIR/BalochistanSenAT115/mintpy/inputs/geometryRadar.h5 --gpsdir /data/lxr/insarlab/GPS/ --scale 1 --scale_key 0.01 --vlim -0.02 0.02 --outfile velocity --outdir /data/lxrtest/BalochistanSenAT/shp/ 
    
    # plot multiple gps sites with extracting incidence and heading angle from geometreRadar 
    plot_geotiff.py geotiff_file --shpdir /data/lxrtest/Balochistan/shp/ --fault multifault.shp nearbyfault.shp --fcolor o m --fstyle s s --refpoi refpoi_DT.shp --gps_putton --type CGPS --search --bbox 30 35 101 106 --scale 1 --scale_key 0.01 --vlim -0.02 0.02 --outfile velocity --outdir /data/lxrtest/BalochistanSenAT/shp/ 

    # plot GPS data alone
    plot_geotiff.py velocity_201504_202010.tiff --gps_putton --type NGPS --search --bbox 37 43 -121 -115 --gpsdir /data/lxr/insarlab/GPS/ --scale 1 --scale_key 0.01 --vlim -0.02 0.02 --outfile velocity3 --outdir ./
"""

fcolor_table = {'b':'black','r':'red','y':'yellow','o':'orange','m':'magenta'}
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
    shp_option.add_argument('--refpoi', dest='refpoi', nargs='?', type=str, help='shp file of reference point. \n')
   
    gps_option = parser.add_argument_group(title='options for plot gps velocity vector field')
    
    gps_option.add_argument('--gps_putton', action='store_true', default=False, help='whether plot gps velocity vector field. \n')

    gps_option.add_argument('--type', dest='Gtype', type=str, nargs=1, help='CGPS: China GPS database; NGPS: Nevada GPS database')
    
    gps_option.add_argument('--search', action='store_true', default=False, help='whether search gps sites based on the given SNWE. \n')
   
    gps_option.add_argument('--bbox', dest='SNWE', type=float, nargs=4, metavar=('S', 'N','W','E'),
                            help='Bounding box of area to be geocoded.\n'+
                            'Include the uppler left corner of the first pixel' + 
                            'and the lower right corner of the last pixel')
    gps_option.add_argument('--gps_name', dest='gps_name', nargs='?', type=str, help='gps site name. \n')

    gps_option.add_argument('--geometry', dest='geometry', nargs='?', type=str, help='geometry data. \n')

    gps_option.add_argument('--gpsdir', dest='gpsdir',nargs='?', type=str, help='directory to store Nevada Geodetic Lab gps data. \n')

    gps_option.add_argument('--scale', dest='scale', nargs=1, type=int, help='parameter for quiver.scale')

    gps_option.add_argument('--scale_key', dest='scale_key', nargs=1, type=float, help='parameter for quiverkey.U \n')

    parser.add_argument('--vlim', dest = 'vlim', nargs=2, type=float, help='max and min value for displacement to display. The unit is m. \n')   
 
    colorbar_option = parser.add_argument_group(title='options for the colorbar: linear colorbar or nonlinear colorbar')

    colorbar_option.add_argument('--nonlinear', action='store_true', default=False, help='whether use nonlinear colorbar.')
 
    #colorbar_option.add_argument('--level', dest='level', type=float, nargs=3, metavar=('L', 'H', 'N'),
    #                              help='the lower bound and higher bound of range that display with nonlinear colorbar' + 
    #                                   'N is the number of interval')

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

def cgps_process(gpsname, geometry):
    """process single gps site""" 
    gps_obj = CGPS(site=gpsname)
    gps_obj.open()
    if not geometry:
        vel_los = ut.enu2los(gps_obj.vel_e, gps_obj.vel_n, gps_obj.vel_u)
        los_vel = vel_los 
    else:
        gps_obj.read_gps_los_velocity(geometry)
        los_vel = gps_obj.vel_los
    
    site_inve = np.array([[gps_obj.site, gps_obj.site_lat, gps_obj.site_lon, los_vel, gps_obj.vel_n, gps_obj.vel_e]])
    
    return site_inve

def ngps_process(gpsname, geometry, gps_dir):
    """process ngps sites""" 
    gps_obj = GPS(site=gpsname, data_dir=gps_dir)
    print('process {} site'.format(gpsname))
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
    if inps.Gtype[0] == 'NGPS':
        site_names, site_lats, site_lons = gps.search_gps(inps.SNWE)
    elif inps.Gtype[0] == 'CGPS':
        site_names, site_lats, site_lons = cgps.search_gps(inps.SNWE)

    sites_baseinfo = np.transpose(np.vstack((np.vstack((np.array([site_names]),np.array([site_lats]))),np.array([site_lons])))) 
    
    return sites_baseinfo

#class nlcmap(object):
#    def __init__(self, cmap, levels):
#        self.cmap = cmap
#        self.N = cmap.N
#        self.monochrome = self.cmap.monochrome
#        self.levels = np.asarray(levels, dtype='float64')
#        self._x = self.levels
#        self.levmax = self.levels.max()
#        self.transformed_levels = np.linspace(0.0, self.levmax,  len(self.levels))

#    def __call__(self, xi, alpha=1.0, **kw):
#        yi = np.interp(xi, self._x, self.transformed_levels)
#        return self.cmap(yi / self.levmax, alpha)

def plot_geotiff(inps, site_infovel=None):
    """read geotiff data"""

    # read raster file
    raster = rasterio.open(inps.input_geotiff[0])
    
    """plot geotiff"""
    raster_data = raster.read(1)
    shape = np.shape(raster_data)
    if inps.vlim == None: 
        raster_min = np.nanmin(raster_data)
        raster_max = np.nanmax(raster_data)
    else:
        raster_min = inps.vlim[0]
        raster_max = inps.vlim[1]
    
    # setting colorbar:linear or nonlinear
#    if inps.nonlinear:
#        level_lower = float(inps.level[0])
#        level_higher = float(inps.level[1])
#        level_num = int(inps.level[2])
#        level_tmp = np.linspace(level_lower, level_higher, level_num)
#        raster_mmin = np.array([raster_min])
#        raster_mmax = np.array([raster_max])
#        levels = np.concatenate((raster_mmin, level_tmp, raster_mmax), axis=0)

#       cmap = nlcmap(plt.cm.jet, levels)     
#    else:        
    cmap = plt.cm.jet

    figure_size = auto_figure_size(shape)
    fig, axes = plt.subplots(1,1,figsize = figure_size)
    ax1 = axes
    print('\n\n***************************************ploting raster data****************************************')
    if inps.nonlinear:
        im = rasterio.plot.show(raster, 1, ax=ax1, norm=colors.SymLogNorm(linthresh=0.001, linscale=0.001,
                                                                    base=10, vmin=np.nanmin(raster_data), 
                                                                    vmax=np.nanmax(raster_data)),
                                  cmap=cmap, alpha=0.8)
    else:
        rasterio.plot.show(raster, 1, ax=ax1, cmap=cmap, vmin=raster_min, vmax=raster_max, alpha=0.8)
   
    if inps.shpdir:
        shpdir = inps.shpdir[0]
        # read shape files
        if inps.refpoi:
            shp_refpoi = geopandas.read_file(shpdir + inps.refpoi)
        
            print('\n\n************************************ploting reference point****************************************')
            shp_refpoi.plot(color='black', ax=ax1, marker='s')
        
        # read fault files and set color/linestyle
        if inps.fault:
            shp_faults = dict()
            for fault in inps.fault:
                shp_faults[fault] = geopandas.read_file(shpdir + fault)
            fault_colors = fcolor(inps.fcolor, shp_faults)
            fault_linestyles = flinestyle(inps.fstyle, shp_faults)
    
            print('\n\n***************************************ploting fault**********************************************')
            for fault in inps.fault:
                shp_faults[fault].plot(color=fault_colors[fault], ax=ax1, linestyle=fault_linestyles[fault], linewidth=1)
                #shp_fault.plot(color='black',ax = ax1,linestyle='dashed',linewidth=1)
   
    # plot gps velocity vector field
    if site_infovel is not None:
        lon = list(site_infovel[:,2].astype(float))
        lat = list(site_infovel[:,1].astype(float))
        x_component = list(site_infovel[:,5].astype(float))
        y_component = list(site_infovel[:,4].astype(float))
        
        quiver_scale = inps.scale[0]
        h1 = ax1.quiver(lon, lat, x_component, y_component, color='black', scale=quiver_scale)
        
        U_parameters = inps.scale_key[0]
        u_label = 'v:' + str(U_parameters) + 'm/s'
        ax1.quiverkey(h1,X=0.09, Y = 0.05, U = U_parameters, angle = 0, label=u_label, labelpos='S', color = 'b',labelcolor = 'b')

    ax1.tick_params(which='both', direction='in', labelsize=18, bottom=True, top=True, left=True, right=True)
   
    # method 1 for drawing colorbar
    if inps.nonlinear:
        cax1 = fig.add_axes([0.92, 0.18, 0.02, 0.8])
        sm1 = plt.cm.ScalarMappable(norm=colors.SymLogNorm(linthresh=0.001, linscale=0.001, base=10, 
                                    vmin=np.nanmin(raster_data), vmax=np.nanmax(raster_data)), cmap=cmap)
        sm1.set_array([])
        sm1.set_clim(np.nanmin(raster_data),np.nanmax(raster_data))
        #ig.colorbar(im,ax=ax1, extend='both')
        cbar = fig.colorbar(sm1,cax1, orientation = 'vertical',format='%.2f')
        cbar.ax.tick_params(labelsize=18)
        cbar.set_ticks(np.array([np.nanmin(raster_data), -0.01, 0, 0.01, np.nanmax(raster_data)]))
   
    else:
        cax = make_axes_locatable(ax1).append_axes("right", "3%", pad="3%")
        norm = mpl.colors.Normalize(vmin=raster_min, vmax=raster_max)
        cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
        cbar.ax.tick_params(labelsize=18)

    # method 2 for drawing colorbar
    #cax1 = fig.add_axes([0.92, 0.18, 0.02, 0.6])
    #sm1 = plt.cm.ScalarMappable(norm=colors.SymLogNorm(linthresh=0.001, linscale=0.001, base=10), cmap=cmap)
    #sm1.set_array([])
    #sm1.set_clim(-1,1)
    #fig.colorbar(sm1,cax1, orientation = 'vertical', label = 'LOS velocity(m/year)')
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
        if inps.search:
            sites_baseinfo = search_gpssites(inps)
            sites_vel = np.zeros([sites_baseinfo.shape[0],5])
            for gpsname, i in zip(sites_baseinfo[:,0],np.arange(sites_baseinfo.shape[0])):
                if inps.Gtype[0] == 'NGPS':
                    gpsdir = inps.gpsdir
                    site_vel = ngps_process(str(gpsname),geodata,gpsdir)
                elif inps.Gtype[0] == 'CGPS':
                    site_vel = cgps_process(str(gpsname), geodata)          
                sites_vel[i,:] = site_vel[0,1:]
            site_infovel = np.concatenate((sites_baseinfo[:,0].reshape(sites_baseinfo.shape[0],1),sites_vel), axis=1)
            #print(sites_inve)
        else:
            gpsname = inps.gps_name
            if inps.Gtype[0] == 'NGPS':
                gpsdir = inps.gpsdir
                site_infovel = ngps_process(gpsname, geodata, gpsdir)
            elif inps.Gtype[0] == 'CGPS':
                site_infovel = cgps_process(str(gpsname), geodata)          
                
            #print(site_inve)    
        plot_geotiff(inps,site_infovel)
    else:
        plot_geotiff(inps)
######################################################################################
if __name__ == '__main__':
    main()
