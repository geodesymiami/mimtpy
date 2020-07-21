##############################################################
# Program is part of MimtPy                                  #
# Purpose: using China GPS data for plot                     #
# Author: Lv Xiaoran, Jul 2020                               #
# Statement: this script is writen based on mintpy gps.py    #
##############################################################
# Utility scripts for GPS handling
# Recommend import:
#     from mimtpy.objects.cgps import GPS

import os
import codecs
from datetime import datetime as dt
import numpy as np

from mintpy.objects.coord import coordinate
from mintpy.utils import readfile

CGPS_dataset = '/data/lxr/insarlab/GPS/CGPS/IndexForStations.txt'

def search_gps(SNWE, start_date=None, end_date=None, site_list_file=CGPS_dataset, pring_msg=True):
    """Used for search available GPS sites within the geo bounding box from CGPS database
     Parameters: SNWE       : tuple of 4 float, indicating (South, North, West, East) in degrees
                start_date : string in YYYYMMDD format
                end_date   : string in YYYYMMDD format
                site_list_file : string.
                min_num_solution : int, minimum number of solutions available
    Returns:    site_names : 1D np.array of string for GPS station names
                site_lats  : 1D np.array for lat
                site_lons  : 1D np.array for lon
    """
    # read CGPS index file
    txt_data = np.loadtxt(site_list_file, dtype=str, delimiter=',', skiprows=1, usecols=(0,1,2,3,4,5,6), encoding='GBK').astype(str)

    site_names = txt_data[:,0]
    site_lons, site_lats = txt_data[:, 1:3].astype(np.float32).T
    site_lons -= np.round(site_lons / (360.)) * 360.
    t_start = np.array([dt.strptime(i, "%Y%m%d") for i in txt_data[:, 4].astype(str)])
    t_end   = np.array([dt.strptime(i, "%Y%m%d") for i in txt_data[:, 5].astype(str)])

    # limit on space
    idx = ((site_lats >= SNWE[0]) * (site_lats <= SNWE[1]) *
           (site_lons >= SNWE[2]) * (site_lons <= SNWE[3]))

    # limit on time
    if start_date:
        t0 = ptime.date_list2vector([start_date])[0][0]
        idx *= t_end >= t0
    if end_date:
        t1 = ptime.date_list2vector([end_date])[0][0]
        idx *= t_start <= t1

    return site_names[idx], site_lats[idx], site_lons[idx]

#################################### Beginning of CGPS class #################################### 
class CGPS:
    """GPS class for CGPS velocity
    Example:
      import matplotlib.pyplot as plt
      from mimtpy.objects.cgps import CGPS
      from mintpy.utils import utils as ut
      gps_obj = GPS(site='GV05', data_dir='~/insarlab/GPS')
      gps_obj.open()
      vel_los = ut.enu2los(gps_obj.vel_e,
                           gps_obj.vel_n,
                           gps_obj.vel_u)
    """

    def __init__(self, site, data_dir='/data/lxr/insarlab/GPS/CGPS/'):
        self.site = site
        self.data_dir = data_dir
        self.source = 'GNSS data products of China earthquake adiministration'

        self.file = '/data/lxr/insarlab/GPS/CGPS/GPS_vel.txt'
        
    def open(self, print_msg=True):
        if not os.path.isfile(self.file):
            raise Exception('GPS velocity dataset cannot be found. Please check!')
        self.read_velocity(print_msg=print_msg)
       
    def read_velocity(self,print_msg=True):
        """read velocity dataset"""
        if print_msg:
            print('Read velocity dataset')
        if not os.path.isfile(self.file):
            raise Exception('GPS velocity dataset cannot be found. Please check!')

        data = np.loadtxt(self.file, dtype=str, delimiter=',', skiprows=3, usecols=(0,1,2,3,4,5,6,7,8,9), encoding='GBK').astype(str)
        site_names = data[:,0]
        id_position = (site_names == self.site)
        self.site_lat = float(data[:,2][id_position])
        self.site_lon = float(data[:,1][id_position])

        self.vel_e = float(data[:,6][id_position]) / 100
        self.std_e = float(data[:,7][id_position]) / 100
        self.vel_n = float(data[:,4][id_position]) / 100
        self.std_n = float(data[:,5][id_position]) / 100
        self.vel_u = data[:,8][id_position].astype(np.float64)[0] / 100
        self.std_u = data[:,9][id_position].astype(np.float64)[0] / 100
      
        return (self.site_lat, self.site_lon, self.vel_e,
                self.std_e, self.vel_n, self.std_n,
                self.vel_u, self.std_u)

#####################################  Utility Functions ###################################
    def get_los_geometry(self, geom_obj, print_msg=False):
        lat = self.site_lat
        lon = self.site_lon
    
        # get LOS geometry
        if isinstance(geom_obj, str):
            # geometry file
            atr = readfile.read_attribute(geom_obj)
            coord = coordinate(atr, lookup_file=geom_obj)
            y, x = coord.geo2radar(lat, lon, print_msg=print_msg)[0:2]
            box = (x, y, x+1, y+1)
            inc_angle = readfile.read(geom_obj, datasetName='incidenceAngle', box=box, print_msg=print_msg)[0][0,0]
            az_angle  = readfile.read(geom_obj, datasetName='azimuthAngle', box=box, print_msg=print_msg)[0][0,0]
            head_angle = ut.azimuth2heading_angle(az_angle)
        elif isinstance(geom_obj, dict):
            # use mean inc/head_angle from metadata
            inc_angle = ut.incidence_angle(geom_obj, dimension=0, print_msg=print_msg)
            head_angle = float(geom_obj['HEADING'])
            # for old reading of los.rdr band2 data into headingAngle directly
            if (head_angle + 180.) > 45.:
                head_angle = ut.azimuth2heading_angle(head_angle)
        else:
            raise ValueError('input geom_obj is neight str nor dict: {}'.format(geom_obj))
        return inc_angle, head_angle   
    
    def displacement_enu2los(self, inc_angle:float, head_angle:float, gps_comp='enu2los'):
        """Convert displacement in ENU to LOS direction
        Parameters: inc_angle  : float, local incidence angle in degree
                    head_angle : float, satellite orbit heading direction in degree
                        from the north, defined as positive in clock-wise direction
                    gps_comp   : string, GPS components used to convert to LOS direction
        Returns:    dis_los : 1D np.array for displacement in LOS direction
                    std_los : 1D np.array for displacement standard deviation in LOS direction
        """
        # get LOS unit vector
        inc_angle *= np.pi/180.
        head_angle *= np.pi/180.
        unit_vec = [np.sin(inc_angle) * np.cos(head_angle) * -1,
                    np.sin(inc_angle) * np.sin(head_angle),
                    np.cos(inc_angle)]
    
        gps_comp = gps_comp.lower()
        if gps_comp in ['enu2los']:
            pass
        elif gps_comp in ['en2los', 'hz2los']:
            unit_vec[2] = 0.
        elif gps_comp in ['u2los', 'up2los']:
            unit_vec[0] = 0.
            unit_vec[1] = 0.
        else:
            raise ValueError('Un-known input gps components:'+str(gps_comp))
    
        # convert ENU to LOS direction
        self.vel_los = (self.vel_e * unit_vec[0]
                        + self.vel_n * unit_vec[1]
                        + self.vel_u * unit_vec[2])
    
        # assuming ENU component are independent with each other
        self.std_los = ((self.std_e * unit_vec[0])**2
                         + (self.std_n * unit_vec[1])**2
                         + (self.std_u * unit_vec[2])**2)**0.5
        return self.vel_los, self.std_los
    
    def read_gps_los_velocity(self, geom_obj, gps_comp:str='enu2los', print_msg=False):
        """Read GPS velocity in LOS direction
            Parameters: geom_obj : dict / str, metadata of InSAR file, or geometry file path
                        gps_comp   : string, GPS components used to convert to LOS direction
            Retruns:    vel   : 1D np.array of displacement in meters
                        std   : 1D np.array of displacement uncertainty in meters
            """
        # read GPS object
        inc_angle, head_angle = self.get_los_geometry(geom_obj)
        self.displacement_enu2los(inc_angle, head_angle, gps_comp=gps_comp)

#################################### End of GPS-UNR class ####################################
