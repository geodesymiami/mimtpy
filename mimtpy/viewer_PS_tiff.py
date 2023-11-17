#!/usr/bin/env python3
#################################################################
# Program is used for viewing the PS point with satellite image #
# Author: Lv Xiaoran                                            #
# Created: Sep 2023                                             #
#################################################################
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import rasterio
from rasterio.mask import mask
from rasterio.plot import show
import geopandas
import pandas as pd
import copy
from shapely.geometry import box
from pyproj import CRS
from matplotlib.backend_bases import MouseButton
from matplotlib import widgets
import scipy
from shapely.geometry import Polygon

from mintpy.utils import readfile, ptime, writefile, utils as ut
from mintpy.objects import timeseries

######################################################################################
NOTE = """
NOTE:
1. Here is some links to download the optical satellite images. But the highest spatial resolution of the open and free dataset is 10m.
https://gisgeography.com/free-satellite-imagery-data-list/
https://skywatch.com/free-sources-of-satellite-data/
2. Planet is a good choose to try. And the Planet's policy is to download 5,000 square kilometers of free video per month
"""

EXAMPLE = """example:

    viewer_PS_tiff.py velocity.h5 --tiff_file tangshan_NE_patch.tiff --geo_file ./inputs/geometryRadar.h5 --two_poi 39.6 118.2 39.7 118.3 --output two_poi.png --outdir ./ 
    
    viewer_PS_tiff.py velocity.h5 --tiff_file tangshan_NE_patch.tiff --geo_file ./inputs/geometryRadar.h5 --subset 39.6 39.7 118.2 118.3 --output subset.png --outdir ./ 

    viewer_PS_tiff.py velocity.h5 --tiff_file tangshan_NE_patch.tiff --geo_file ./inputs/geometryRadar.h5 --subset 39.6 39.7 118.2 118.3 --vlim -0.4 0.4 --output subset.png --outdir ./ 

    viewer_PS_tiff.py velocity.h5 --tiff_file tangshan_NE_patch.tiff --geo_file ./inputs/geometryRadar.h5 --subset 39.6 39.7 118.2 118.3 --vlim -0.4 0.4 --ts_file ./TangshanSenAT69/miaplpy_NE_201410_202212/network_single_reference/timeseries.h5 --interactive 

    viewer_PS_tiff.py ./stamps_results/2y/ps_plot_v.mat --tiff_file ./TangshanSenDT149/Tangshan_NE_patch.tif --shp_file ./shpfile/road.shp --geo_file ./stamps_results/2y/ps2.mat --subset 39.50 39.55 118.30 118.35 --output tiff_Try.png --outdir ./ --vlim -3 3

    viewer_PS_tiff.py ./stamps_results/2y/ps_plot_v.mat --shp_file ./shpfile/road.shp --geo_file ./stamps_results/2y/ps2.mat --subset 39.50 39.55 118.30 118.35 --output tiff_Try.png --outdir ./ --vlim -3 3
    
    viewer_PS_tiff.py ./stamps_results/2y/ps_plot_v.mat --tiff_file ./TangshanSenDT149/Tangshan_NE_patch.tif --geo_file ./stamps_results/2y/ps2.mat --subset 39.50 39.55 118.30 118.35 --output tiff_Try.png --outdir ./ --vlim -3 3  

    viewer_PS_tiff.py ./TangshanSenAT69/miaplpy_NE_201410_202212/network_single_reference/velocity_msk.h5 ./stamps_results/2y/ps_plot_v.mat --tiff_file ./TangshanSenDT149/Tangshan_NE_patch.tif --shp_file ./shpfile/road.shp --geo_file ./TangshanSenAT69/miaplpy_NE_201410_202212/network_single_reference/inputs/geometryRadar.h5 ./stamps_results/2y/ps2.mat --subset 39.50 39.55 118.30 118.35 --output tiff_Try.png --outdir ./ --vlim -3 3
    
"""

def create_parser():
    parser = argparse.ArgumentParser(description='View PS points on satellite image and/or shape file vector',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE+NOTE)

    parser.add_argument('input_file', nargs='+', type=str, help='velocity file under radar coordinate.\n')

    parser.add_argument('--tiff_file', nargs='*', type=str, help='high spatial resolution image with geotiff format.\n')
    
    parser.add_argument('--geo_file', nargs='+', type=str, help='geometryRadar data corresponding to the Wphase file.\n')
    
    parser.add_argument('--shp_file', nargs='*', type=str, 
                        help='shapefile of vector. For example, *.osm files downloaded from OpenStreetMap can be converted to *.shp formate by QGIS'
                        ' and plotted with PS points.\n')
    
    parser.add_argument('--two_poi', nargs='*', type=float, help='lat1 lon1 lat2 lon2 of two points.\n')
    
    parser.add_argument('--subset', nargs='*', type=float, help='lat1 lat2 lon1 lon2 of two points.\n')
    
    parser.add_argument('--vlim', nargs='*', type=float, help='velocity range to display.\n')
   
    parser.add_argument('--output', nargs='?', type=str, help='output name of differential wrap phase timeseries txt file.\n')

    parser.add_argument('--outdir', nargs='?', type=str, help='output directory to store the differential wrap phase timeseries file.\n')
    
    parser.add_argument('--ts_file', nargs='?', type=str, help='time series file under radar coordinate.\n')

    parser.add_argument('--interactive',action='store_true', default=False, help='whether displace interactive map. \n')

    parser.add_argument('--display',action='store_true', default=False, help='whether displace the final plot. \n')
    
    parser.add_argument('--buffer', nargs='*', type=float, help='buffer distance. Unit is meter.\n')
    
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps

def find_row_col(lat_data, lon_data, lat_p1, lon_p1):
    lat_data1 = copy.deepcopy(lat_data)
    lon_data1 = copy.deepcopy(lon_data)

    lat_data1[np.where(lat_data1 != lat_p1)] = 0
    lat_data1[np.where(lat_data1 == lat_p1)] = 1

    lon_data1[np.where(lon_data1 != lon_p1)] = 0
    lon_data1[np.where(lon_data1 == lon_p1)] = 1

    flag1 = lat_data1 * lon_data1

    pos_row1 = np.where(flag1 == 1)[0]
    pos_col1 = np.where(flag1 == 1)[1]

    return pos_row1, pos_col1

def subset_PS_poi(lat_min, lat_max, lon_min, lon_max, lat_data, lon_data, vel_mask):
    lat_data1 = copy.deepcopy(lat_data)
    lon_data1 = copy.deepcopy(lon_data)

    lat_data1[np.where((lat_data1 > lat_min) & (lat_data1 < lat_max))] = 1
    lon_data1[np.where((lon_data1 > lon_min) & (lon_data1 < lon_max))] = 1

    flag1 = lat_data1 * lon_data1

    flag1[np.isnan(vel_mask)] = np.nan

    pos_row = np.where(flag1 == 1)[0]
    pos_col = np.where(flag1 == 1)[1]

    return pos_row, pos_col

def generate_geopandas(lat_sub, lon_sub, vel_sub, crs):
    df = pd.DataFrame(
        {
            "Latitude": lat_sub,
            "Longitude": lon_sub,
            "Value": vel_sub * 100,
        }
    )
    # convert to shapely object
    gdf = geopandas.GeoDataFrame(
        df, geometry=geopandas.points_from_xy(df.Longitude, df.Latitude), crs=crs
    )

    return gdf

def plot_tiff_PS(file_src, gdf_obj, inps, gdf_obj_s=None):
    cmap = plt.cm.jet
    figure_size = [8.0, 8.0]
    fig, axes = plt.subplots(1, 1, figsize=figure_size)
    ax1 = axes

    if inps.vlim is not None:
        vmin = inps.vlim[0]
        vmax = inps.vlim[1]
    else:
        vmin = np.min(gdf_obj['Value'])
        vmax = np.max(gdf_obj['Value'])

    # generate bounding box
    if inps.subset is not None:
        bbox = box(inps.subset[2], inps.subset[0], inps.subset[3], inps.subset[1])
    else:
        bbox = box(np.nanmin(gdf_obj['Longitude']) - 0.005, np.nanmin(gdf_obj['Latitude']) - 0.005, np.nanmax(gdf_obj['Longitude']) + 0.005, np.nanmax(gdf_obj['Latitude']) + 0.005)
    geo_box = geopandas.GeoDataFrame({'geometry': bbox}, index=[0], crs=file_src.crs)

    # start plot
    if inps.tiff_file is not None:
        def getFeatures(gdf):
            import json
            return [json.loads(gdf.to_json())['features'][0]['geometry']]
        coords_box = getFeatures(geo_box)
        tif_clip, transform_clip = mask(file_src, coords_box, crop=True)
        meta_clip = file_src.meta.copy()
        meta_clip.update({"driver": "GTiff",
                          "height": tif_clip.shape[1],
                          "width": tif_clip.shape[2],
                          "transform": transform_clip})
        # write the cut geotif file
        tif_clip_outname = os.path.dirname(inps.tiff_file[0]) + '/' + os.path.basename(inps.tiff_file[0]).split('.')[0] + '_clip.tif'
        with rasterio.open(tif_clip_outname, "w", **meta_clip) as dest:
            dest.write(tif_clip)
        bg_tif = rasterio.open(tif_clip_outname)

        # delete the cut geotif
        os.remove(tif_clip_outname)

        # plot the high spatial resolution optical image
        show(bg_tif, ax=ax1, alpha=0.8)
        # plot shp file
        if inps.shp_file is not None:
            shp_file = geopandas.read_file(inps.shp_file[0])
            shp_file_clip = shp_file.clip(geo_box)
            old_crs = shp_file.crs
            new_crs = CRS('EPSG:3857')
            if inps.buffer is not None:
                shp_file_clip = shp_file_clip.to_crs(new_crs)
                shp_file_clip_buffer = shp_file_clip.buffer(inps.buffer[0])
                shp_file_clip_buffer = shp_file_clip_buffer.to_crs(old_crs)

                # mask the PS points according to the buffer
                gdf_obj_msk = geopandas.clip(gdf_obj, shp_file_clip_buffer)
                if gdf_obj_s is not None:
                    gdf_obj_s_msk = geopandas.clip(gdf_obj_s, shp_file_clip_buffer)
                # plot the masked PS points
                gdf_obj_msk.plot('Value', ax=ax1, cmap=cmap, vmin=vmin, vmax=vmax, markersize=1)
                if gdf_obj_s is not None:
                    gdf_obj_s_msk.plot('Value', ax=ax1, cmap=plt.cm.gray, vmin=vmin, vmax=vmax, markersize=8, alpha=0.7)
            
                shp_file_clip = shp_file_clip.to_crs(old_crs)
                shp_file_clip.plot(color='black', ax=ax1, linestyle='solid', linewidth=0.3)
            else:
                # plot PS points
                gdf_obj.plot('Value', ax=ax1, cmap=cmap, vmin=vmin, vmax=vmax, markersize=1)
                if gdf_obj_s is not None:
                    gdf_obj_s.plot('Value', ax=ax1, cmap=plt.cm.gray, vmin=vmin, vmax=vmax, markersize=8, alpha=0.7)
                
                shp_file_clip.plot(color='black', ax=ax1, linestyle='solid', linewidth=0.3)
                #shp_file_clip.plot(ax=ax1, color='black', marker='s', markersize=10, alpha=0.3)
        else:          
            # plot PS points
            gdf_obj.plot('Value', ax=ax1, cmap=cmap, vmin=vmin, vmax=vmax, markersize=1)
            if gdf_obj_s is not None:
                gdf_obj_s.plot('Value', ax=ax1, cmap=plt.cm.gray, vmin=vmin, vmax=vmax, markersize=8, alpha=0.7)
    
    else:
        shp_file_clip = file_src.clip(geo_box)
        old_crs = file_src.crs
        new_crs = CRS('EPSG:3857')
        if inps.buffer is not None:
            shp_file_clip = shp_file_clip.to_crs(new_crs)
            shp_file_clip_buffer = shp_file_clip.buffer(inps.buffer[0])
            shp_file_clip_buffer = shp_file_clip_buffer.to_crs(old_crs)

            # mask the PS points according to the buffer
            gdf_obj_msk = geopandas.clip(gdf_obj, shp_file_clip_buffer)
            if gdf_obj_s is not None:
                gdf_obj_s_msk = geopandas.clip(gdf_obj_s, shp_file_clip_buffer)
            # plot the masked PS points
            gdf_obj_msk.plot('Value', ax=ax1, cmap=cmap, vmin=vmin, vmax=vmax, markersize=1)
            if gdf_obj_s is not None:
                gdf_obj_s_msk.plot('Value', ax=ax1, cmap=plt.cm.gray, vmin=vmin, vmax=vmax, markersize=8, alpha=0.7)
        
            shp_file_clip = shp_file_clip.to_crs(old_crs)
            shp_file_clip.plot(color='black', ax=ax1, linestyle='solid', linewidth=0.3)
        else:
            # plot PS points
            gdf_obj.plot('Value', ax=ax1, cmap=cmap, vmin=vmin, vmax=vmax, markersize=1)
            if gdf_obj_s is not None:
                gdf_obj_s.plot('Value', ax=ax1, cmap=plt.cm.gray, vmin=vmin, vmax=vmax, markersize=8, alpha=0.7)
            
            shp_file_clip.plot(color='black', ax=ax1, linestyle='solid', linewidth=0.3)
            #shp_file_clip.plot(ax=ax1, color='black', marker='s', markersize=10, alpha=0.3)
         
    ax1.tick_params(which='both', direction='in', labelsize=8, bottom=True, top=True, left=True, right=True)
    cax1 = fig.add_axes([0.92, 0.18, 0.02, 0.6])
    sm1 = plt.cm.ScalarMappable(cmap=cmap)
    sm1.set_array([])
    sm1.set_clim(vmin=vmin, vmax=vmax)
    cb = fig.colorbar(sm1, cax1, orientation='vertical', format='%.2f')
    cb.ax.tick_params(labelsize=10)
    cb.set_ticks(np.linspace(vmin, vmax, 5))
    font2 = {'family': 'serif',
             'weight': 'normal',
             'size': 10.}
    cb.set_label('velocity [cm/year]', fontdict=font2)

    # set axis
    font1 = {'family': 'serif',
             'weight': 'normal',
             'size': 10.} 
    ax1.set_xlabel('Longitude', font1)
    ax1.set_ylabel('Latitude', font1)
    labels = ax1.get_xticklabels() + ax1.get_yticklabels()
    [label.set_fontname('serif') for label in labels]

    if inps.display:
        plt.show()    
    else:
        # save png
        if inps.outdir is not None:
            fig_name = inps.outdir + '/' + inps.output
            fig.savefig(fig_name, dpi=300, bbox_inches='tight')
        else:
            raise ValueError('Please enter output name and dir!')

class interactive_map:
    def __init__(self, *, file_src, gdf_obj, ts_sub, inps):
        self.src = file_src
        self.gdf_obj = gdf_obj
        self.inps = inps
        self.ts_sub = ts_sub

        # figure variables
        self.fig_img = None
        self.ax_img = None
        self.cbar_img = None
        self.img = None

        self.fig_pts = None
        self.ax_pts = None
        self.latlon = [gdf_obj['Latitude'].values[0], gdf_obj['Longitude'].values[0]]
        self.index = 0

        return

    def plot(self):

        # Figure 1 - Axes 1 - Amplitude Map
        self.fig_img = plt.figure()
        self.ax_img = self.fig_img.add_axes([0.05, 0.1, 0.8, 0.8])
        self.plot_tiff_image()

        # Figure 2 - Axes 1 - Time Series Displacement - Point
        self.fig_pts = plt.figure()
        self.ax_pts = self.fig_pts.add_subplot(111)
        self.plot_point_timeseries(self.latlon)

        # Final linking of the canvas to the plots.
        self.fig_img.canvas.mpl_connect('button_press_event', self.update_plot_timeseries)
        self.fig_pts.canvas.mpl_connect('button_press_event', self.update_plot_timeseries)

        plt.show() 

    def update_plot_timeseries(self, event):
        """Event function to get y/x from button press"""
        self.index += 1
        if event.button is MouseButton.RIGHT:
            # get lat/lon
            plat = event.ydata
            plon = event.xdata
            self.latlon = (plat, plon)
            # plot time-series displacement
            self.plot_point_timeseries(self.latlon)

            plt.show()

    def plot_point_timeseries(self, latlon):
        """Plot point displacement time-series at pixel [y, x]
        Parameters: yx     : list of 2 int
        """
        print('Plot the time series of new PS point')
        # get time series data
        poi_lat = latlon[0]
        poi_lon = latlon[1]
        
        lat_index = self.gdf_obj['Latitude'].values.reshape(-1, 1)
        lon_index = self.gdf_obj['Longitude'].values.reshape(-1, 1)

        pos_row, pos_col = find_row_col(lat_index, lon_index, poi_lat, poi_lon)
        if len(pos_row):
            print('Get the PS point')

        ts_dat = self.ts_sub[:, pos_row] * 100 + self.index * 10
        # plot
        if not np.all(np.isnan(ts_dat)):
            self.ax_pts.scatter(np.arange(1,len(ts_dat) + 1), ts_dat, s=10, marker='o')
            #self.ax_pts.plot(ts_dat)
            #self.ax_pts_hist.hist(ts_dat, bins=30, orientation="horizontal")

        self.ax_pts.set_ylabel('Displacement [cm]')
        self.ax_pts.set_xlabel('Date Num')
        self.ax_pts.set_title('Time Series')

        # displacement the geo_info of selectd PS point
        geoinfo_str = 'The latitude is ' + str(poi_lat) + '\nThe longitude is ' + str(poi_lon)
        print(geoinfo_str)
        self.ax_pts.text(1, -5, geoinfo_str, fontsize=10, bbox=dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.7))

        # update figure
        self.fig_pts.canvas.draw()

    def plot_tiff_image(self):
        # Title and Axis Label
        self.ax_img.set_ylabel('Longitude')
        self.ax_img.set_xlabel('Latitude')

        if self.inps.vlim is not None:
            vmin = self.inps.vlim[0]
            vmax = self.inps.vlim[1]
        else:
            vmin = np.min(self.gdf_obj['Value'])
            vmax = np.max(self.gdf_obj['Value'])
        
        bbox = box(self.inps.subset[2], self.inps.subset[0], self.inps.subset[3], self.inps.subset[1])
        geo_box = geopandas.GeoDataFrame({'geometry': bbox}, index=[0], crs=self.src.crs)
        
        if self.inps.tiff_file is not None:
            def getFeatures(gdf):
                import json
                return [json.loads(gdf.to_json())['features'][0]['geometry']]
            coords_box = getFeatures(geo_box)

            tif_clip, transform_clip = mask(self.src, coords_box, crop=True)
            self.transform, self.width, self.height = self.src.transform, self.src.width, self.src.height
            self.crs = self.src.crs
            meta_clip = self.src.meta.copy()
            meta_clip.update({"driver": "GTiff",
                              "height": tif_clip.shape[1],
                              "width": tif_clip.shape[2],
                              "transform": transform_clip})
            # write the cut geotif file
            tif_clip_outname = os.path.dirname(self.inps.tiff_file[0]) + '/' + os.path.basename(self.inps.tiff_file[0]).split('.')[0] + '_clip.tif'
            with rasterio.open(tif_clip_outname, "w", **meta_clip) as dest:
                dest.write(tif_clip)
            bg_tif = rasterio.open(tif_clip_outname)

            # delete the cut geotif
            os.remove(tif_clip_outname)

            # plot the high spatial resolution optical image
            show(bg_tif, ax=self.ax_img, alpha=0.8)
            # plot PS points
            self.gdf_obj.plot('Value', ax=self.ax_img, cmap=plt.cm.jet, vmin=vmin, vmax=vmax, markersize=5)
            # plot shp file
            if self.inps.shp_file is not None:
                shp_file = geopandas.read_file(self.inps.shp_file[0])
                shp_file_clip = shp_file.clip(geo_box)
                shp_file_clip.plot(color='yellow', ax=self.ax_img, linestyle='solid', linewidth=0.3)

        else:
            shpfile_clip = self.src.clip(geo_box)
            # plot the shape file background
            shpfile_clip.plot(ax=self.ax_img)
            # plot the PS points
            self.gdf_obj.plot('Value', ax=self.ax_img, cmap=plt.cm.jet, vmin=vmin, vmax=vmax, markersize=5)
            
        self.ax_img.tick_params(which='both', direction='in', labelsize=8, bottom=True, top=True, left=True, right=True)
        cax1 = self.fig_img.add_axes([0.80, 0.18, 0.02, 0.6])
        sm1 = plt.cm.ScalarMappable(cmap=plt.cm.jet)
        sm1.set_array([])
        sm1.set_clim(vmin=vmin, vmax=vmax)
        cb = self.fig_img.colorbar(sm1, cax1, orientation='vertical', format='%.2f')
        cb.ax.tick_params(labelsize=10)
        cb.set_ticks(np.linspace(vmin, vmax, 5))
        font2 = {'family': 'serif',
                 'weight': 'normal',
                 'size': 10.}
        cb.set_label('velocity [cm/year]', fontdict=font2)

def read_miaplpy_data(input_file, geo_file, ts_file=None):
    vel_data = readfile.read(input_file)[0]
    lat_data = readfile.read(geo_file, datasetName='latitude')[0]
    lon_data = readfile.read(geo_file, datasetName='longitude')[0]

    if ts_file is not None:
        ts_data = readfile.read(ts_file, datasetName='timeseries')[0]
        return vel_data, lat_data, lon_data, ts_data
    else:
        return vel_data, lat_data, lon_data

def read_stamps_data(input_file, geo_file, ts_file=None):
    vel_data = scipy.io.loadmat(input_file)['ph_disp'] / 1000 # the unit of velocity is m/yr
    lon_data = np.float32(scipy.io.loadmat(geo_file)['lonlat'][:,0].reshape(-1, 1))
    lat_data = np.float32(scipy.io.loadmat(geo_file)['lonlat'][:,1].reshape(-1, 1))

    if ts_file is not None:
        ts_data = -1 * np.transpose(scipy.io.loadmat(ts_file)['ph_disp']) *0.056 / (4 * np.pi) # radians to meter
        # the size of matrix is [date_num, PS_num]
        return vel_data, lat_data, lon_data, ts_data
    else:
        return vel_data, lat_data, lon_data

def main():
    inps = cmd_line_parse()
    
    if len(inps.input_file) == 1:
        file_extension = os.path.splitext(inps.input_file[0])[1] 
        if file_extension == '.h5':
            if inps.ts_file is not None:
                vel_data, lat_data, lon_data, ts_data = read_miaplpy_data(inps.input_file[0], inps.geo_file[0], inps.ts_file)
            else:
                vel_data, lat_data, lon_data = read_miaplpy_data(inps.input_file[0], inps.geo_file[0])
            vel_mask = copy.deepcopy(vel_data)
            vel_mask[np.isnan(vel_mask)] = 1
        else:
            if inps.ts_file is not None:
                vel_data, lat_data, lon_data, ts_data = read_stamps_data(inps.input_file[0], inps.geo_file[0], inps.ts_file)
            else:
                vel_data, lat_data, lon_data = read_stamps_data(inps.input_file[0], inps.geo_file[0])
            vel_mask = vel_data

    elif len(inps.input_file) == 2:
        vel_data, lat_data, lon_data = read_miaplpy_data(inps.input_file[0], inps.geo_file[0])
        vel_mask = copy.deepcopy(vel_data)
        vel_mask[np.isnan(vel_mask)] = 1

        vel_data_s, lat_data_s, lon_data_s = read_stamps_data(inps.input_file[1], inps.geo_file[1])
        vel_mask_s = vel_data_s
    else:
        raise ValueError('viewer_PS_tiff.py can only process two datasets at same time!')

    # read the geotiff file
    if inps.tiff_file is not None:
        file_src = rasterio.open(inps.tiff_file[0])
        file_crs = file_src.crs
    # read the shape file
    elif inps.shp_file is not None:
        file_src = geopandas.read_file(inps.shp_file[0])
        file_crs = file_src.crs
    else:
        raise ValueError('Please at least provide one of the Geotiff file and shapefile file!') 

    if inps.two_poi is not None:
        lat_p1 = inps.two_poi[0]
        lon_p1 = inps.two_poi[1]
        lat_p2 = inps.two_poi[2]
        lon_p2 = inps.two_poi[3]
        pos_row1, pos_col1 = find_row_col(lat_data, lon_data, lat_p1, lon_p1)
        pos_row2, pos_col2 = find_row_col(lat_data, lon_data, lat_p2, lon_p2)

        vel_val1 = vel_mask[pos_row1, pos_col1]
        vel_val2 = vel_mask[pos_row2, pos_col2]

        lat_sub = np.array([lat_p1, lat_p2])
        lon_sub = np.array([lon_p1, lon_p2])
        vel_sub = np.array([vel_val1[0], vel_val2[0]])

        gdf_obj = generate_geopandas(lat_sub, lon_sub, vel_sub, file_crs)

    elif inps.subset is not None:
        lat_min = inps.subset[0] 
        lat_max = inps.subset[1] 
        lon_min = inps.subset[2] 
        lon_max = inps.subset[3]

        pos_row, pos_col = subset_PS_poi(lat_min, lat_max, lon_min, lon_max, lat_data, lon_data, vel_mask)
        lat_sub = lat_data[pos_row, pos_col] 
        lon_sub = lon_data[pos_row, pos_col] 
        vel_sub = vel_data[pos_row, pos_col]
        
        if len(inps.input_file) == 2:
            pos_row_s, pos_col_s = subset_PS_poi(lat_min, lat_max, lon_min, lon_max, lat_data_s, lon_data_s, vel_mask_s)
            lat_sub_s = lat_data_s[pos_row_s, pos_col_s] 
            lon_sub_s = lon_data_s[pos_row_s, pos_col_s] 
            vel_sub_s = vel_data_s[pos_row_s, pos_col_s]
            gdf_obj_s = generate_geopandas(lat_sub_s, lon_sub_s, vel_sub_s, file_crs)
             

        print('There are totally %d pixels in subset region' % len(lat_sub))  
        if inps.ts_file is not None:
            if file_extension == '.h5':
                ts_sub = ts_data[:, pos_row, pos_col] # size is [date_num, PS_num]
            else:
                ts_sub = ts_data[:, pos_row]
        gdf_obj = generate_geopandas(lat_sub, lon_sub, vel_sub, file_crs)

    # plot the geo_tiff and PS points together
    if inps.interactive:
        print('Active the interactive Map for the subset region!')
        map_obj = interactive_map(file_src=file_src, gdf_obj=gdf_obj, ts_sub=ts_sub, inps=inps)
        map_obj.plot() 
    else:
        if len(inps.input_file) == 1:
            plot_tiff_PS(file_src, gdf_obj, inps)
        else:
            plot_tiff_PS(file_src, gdf_obj, inps, gdf_obj_s)
######################################################################################
if __name__ == '__main__':
    main()
