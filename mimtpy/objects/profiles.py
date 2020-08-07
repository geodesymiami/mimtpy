##############################################################
# Program is part of MimtPy                                  #
# Purpose: processing profiles when doing concatenation      #
# Author: Lv Xiaoran, Jul 2020                               #
##############################################################
# Utility scripts for profiles when doing concatenation
# Recommend import:
#     from mimtpy.objects.profiles import Profile
#     import mimtpy.objects.profiles as profiles

import os
import re
import copy
import json
import math
import matplotlib.pyplot as plt
import numpy as np

def search_profiles(profile_num, over_lat0, over_lon0, over_lat1, over_lon1, m_atr, s_atr):
    """Search lon/lat of one point located on the profiles, respectively
    Parameters: profile_num      : int, number of profiles
                m_atr : dictory, metadata for master dataset
                s_atr : dictory, metadata for slave dataset
                min_num_solution : int, minimum number of solutions available
    Returns:    profile_catalog : 2D np.array of string for number of profile, lon/lat of point that was passed by profile 
    """
    m_polygon = m_atr['scene_footprint']
    m_footprint = re.findall(r'([\d+\.]+)',m_polygon)

    s_polygon = s_atr['scene_footprint']
    s_footprint = re.findall(r'([\d+\.]+)',s_polygon)

    # obtain the overlapping footprint
    # footprint of master and slave track
    m_outline = track_outline_matrix(m_footprint)
    s_outline = track_outline_matrix(s_footprint)

    # identify points that located in the overlapping region 
    m_idx_tmp = ((m_outline[:,1] >= over_lat1) * (m_outline[:,1] <= over_lat0) *
                 (m_outline[:,0] >= over_lon0) * (m_outline[:,0] <= over_lon1))
    s_idx_tmp = ((s_outline[:,1] >= over_lat1) * (s_outline[:,1] <= over_lat0) *
                 (s_outline[:,0] >= over_lon0) * (s_outline[:,0] <= over_lon1))
    
    # get the top points for the footprint overlapping region
    overlap_footprint = np.vstack((m_outline[m_idx_tmp],s_outline[s_idx_tmp])) 
    result_sort = np.sort(overlap_footprint[:,1],axis=0)
    # the second and third latitude after sorting are two points for the overlapping footprint
    target1 = result_sort[1]
    target2 = result_sort[2]

    # find the m_idx and s_idx for target1 and target2
    if (m_outline[:,1] == target1).any(): 
        m_idx = np.argwhere(m_outline[:,1] == target1)[0,0] 
    else:
        s_idx = np.argwhere(s_outline[:,1] == target1)[0,0]

    if (m_outline[:,1] == target2).any(): 
        m_idx = np.argwhere(m_outline[:,1] == target2)[0,0] 
    else:
        s_idx = np.argwhere(s_outline[:,1] == target2)[0,0]
    
    # define upper horizontal line of the overlapping footprint
    if m_outline[m_idx,:][1] > s_outline[s_idx,:][1]:
        if m_idx == 0 or m_idx == 3:
            m_end_idx = 3 - m_idx
        if s_idx == 0 or s_idx == 1:
            s_end_idx = 1 - s_idx
        elif s_idx == 2 or s_idx == 3:
            s_end_idx = 5 - s_idx        
    else:
        if s_idx == 0 or s_idx == 3:
            s_end_idx = 3 - s_idx
        if m_idx == 0 or m_idx == 1:
            m_end_idx = 1 - m_idx
        elif m_idx == 2 or m_idx == 3:
            m_end_idx = 5 - m_idx

    line0_start = m_outline[m_idx,:]
    line0_end = m_outline[m_end_idx,:]
    line1_start = s_outline[s_idx,:]
    line1_end = s_outline[s_end_idx,:]
 
    cross_lon, cross_lat = calculate_cross_lines(line0_start, line0_end, line1_start, line1_end)

    # calculate point for each profile
    if m_outline[m_idx,:][1] > s_outline[s_idx,:][1]:
        profile_start = m_outline[m_idx,:]
    else:
        profile_start = s_outline[s_idx,:]

    # create profile_No profile_lon profile_lat
    profiles_catalog = np.empty((profile_num,3),dtype=float)
    x1 = profile_start[0]
    y1 = profile_start[1]
    x2 = cross_lon
    y2 = cross_lat
    profile_No = np.arange(1,profile_num + 1)
    grade_factor = profile_No / (profile_num + 1)
    profile_lons = x1 + grade_factor * (x2 - x1)
    profile_lats = y1 + grade_factor * (y2 - y1)

    profiles_catalog[:,0] = np.transpose(profile_No)
    profiles_catalog[:,1] = np.transpose(profile_lons)
    profiles_catalog[:,2] = np.transpose(profile_lats)
    
    return profiles_catalog

def track_outline_matrix(footprint):
    """store lon/lat of each points for track as dictionary"""
    # lon/lat of four points: upperright; lowerright; lowerleft; upperleft 
    outline = dict()
    lonlats = footprint
    # the structure of outline matrix is:
    #             longitude   latitude
    # upper right
    # lower right
    # lower left
    # upper left
    outline = np.array([[float(lonlats[0]),float(lonlats[1])],[float(lonlats[2]),float(lonlats[3])],
                       [float(lonlats[4]),float(lonlats[5])],[float(lonlats[6]),float(lonlats[7])]])
    
    return outline

def calculate_cross_lines(line0_pos0, line0_pos1, line1_pos0, line1_pos1):
    """calculate the intersect point for two lines
       param line0_pos0: the start point for the first line
       param line0_pos1: the end point for the first line
       param line1_pos0: the start point for the second line
       param line1_pos1: the end point for the second line 
    """
    line0_a =line0_pos0[1] - line0_pos1[1]
    line0_b = line0_pos1[0] - line0_pos0[0]
    line0_c = line0_pos0[0] * line0_pos1[1] - line0_pos1[0] * line0_pos0[1]
    line1_a =line1_pos0[1] - line1_pos1[1]
    line1_b = line1_pos1[0] - line1_pos0[0]
    line1_c = line1_pos0[0] *line1_pos1[1] - line1_pos1[0] * line1_pos0[1]
    d = line0_a * line1_b - line1_a * line0_b
    if d == 0:
        # if these two lines are overlap, there is no intersect point
        return None

    intersect_lon = (line0_b * line1_c - line1_b * line0_c) * 1.0 / d
    intersect_lat = (line0_c * line1_a - line1_c * line0_a) * 1.0 / d

    return intersect_lon, intersect_lat

def profile_average(profile_num, profile_dict_list):
    """for multiprofiles, calculate the average value
       the input parameters is a list composed by dict
       for example:
       [{'NO':1, 'p_start':array([lon_s,lat_s]), 'p_end':array([lon_e,lat_e]), 'm_data':array([1, 2, 3, 4, 5]),'s_data':array([1, 2, 3, 4, 5])},{}...]
       the output:
       [{'NO':1, 'p_start':array([lon_s,lat_s]), 'p_end':array([lon_e,lat_e]), 'm_data':array([1, 2, 3, 4, 5]),'s_data':array([1, 2, 3, 4, 5])},{}...
        ,{'NO':'average','m_data':array([]),'s_data':([])}]
    """
    # length of profile
    length = len(profile_dict_list[0]['m_data'])
    m_sum = np.zeros(length, dtype=float)
    s_sum = np.zeros(length, dtype=float)
    for pro in profile_dict_list:        
       m_sum += pro['m_data']
       s_sum += pro['s_data']

    m_average = m_sum / profile_num
    s_average = s_sum / profile_num

    profile_average = dict()
    profile_average['NO'] = 'average'
    profile_average['m_data'] = m_average
    profile_average['s_data'] = s_average
    
    profile_dict_final = copy.deepcopy(profile_dict_list)
    profile_dict_final.append(profile_average)

    return profile_dict_final

def profiles_plot(profile_dict_final, m_name, s_name, outdir):
    """plot multi profiles with the average value""" 
    figure_size = [10,8]
    fig,axes = plt.subplots(1,1,figsize = figure_size)
    ax1 = axes
    print('*****************************ploting profile************************')       
    length = len(profile_dict_final[0]['m_data'])
    x_axis = np.arange(1,length + 1)
    for pro in profile_dict_final:
        if pro['NO'] != 'average':
            ax1.scatter(x_axis, pro['m_data'], marker='.', c='lightgray',alpha=0.5)
            ax1.scatter(x_axis, pro['s_data'], marker='.', c='lightskyblue',alpha=0.5)
        else:
            ax1.plot(x_axis, pro['m_data'], color='black', linestyle='-', label=m_name)
            ax1.plot(x_axis, pro['s_data'], color='blue', linestyle='-', label=s_name)

    ax1.tick_params(which='both', direction='in', labelsize=18, bottom=True, top=True, left=True, right=True)
    font1 = {'family' : 'serif',
             'weight': 'normal',
             'size' : 18.}
    ax1.set_xlabel('Distance [km]',font1)
    ax1.set_ylabel('LOS Displacement [cm]',font1)
    
    # set x label to km
    lon_s = copy.deepcopy(profile_dict_final[0]['p_start'][0])
    lat_s = copy.deepcopy(profile_dict_final[0]['p_start'][1])
    lon_e = copy.deepcopy(profile_dict_final[0]['p_end'][0])
    lat_e = copy.deepcopy(profile_dict_final[0]['p_end'][1])
    distance = distance_2points(lon_s, lat_s, lon_e, lat_e)
    xtick = np.linspace(0, int(np.ceil(distance)), num=10, endpoint=True)
    
    ax1.set_xticks(np.linspace(0, length, num=10, endpoint=True))
    ax1.set_xticklabels([str(int(round(i))) for i in xtick])
    
    labels = ax1.get_xticklabels() + ax1.get_yticklabels()
    [label.set_fontname('serif') for label in labels]
    
    ax1.legend(loc='upper left', prop=font1)
    
    #save figure
    fig_name = 'Profiles_' + m_name + '_' + s_name + '.png'
    fig_output = outdir + fig_name
    fig.savefig(fig_output, dpi=300, bbox_inches='tight')

def profile_plot(lon_s, lat_s, lon_e, lat_e, m_profile, s_profile, m_name, s_name, outdir):
#def profile_plot(self, m_name, s_name):
    """plot master and slave data along one profile"""
    # plot two profiles
    figure_size = [10,8]
    fig,axes = plt.subplots(1,1,figsize = figure_size)
    ax1 = axes
    print('*****************************ploting profile************************')       
    x_axis = np.arange(1,len(m_profile)+1)
    ax1.plot(x_axis, m_profile, color='black', linestyle='-', label=m_name)
    ax1.plot(x_axis, s_profile, color='blue', linestyle='-', label=s_name)

    ax1.tick_params(which='both', direction='in', labelsize=18, bottom=True, top=True, left=True, right=True)
    font1 = {'family' : 'serif',
             'weight': 'normal',
             'size' : 18.}
    ax1.set_xlabel('Distance [km]',font1)
    ax1.set_ylabel('LOS Displacement [cm]',font1)
    
    # set x label to km
    distance = distance_2points(lon_s, lat_s, lon_e, lat_e)
    xtick = np.linspace(0, int(np.ceil(distance)), num=10, endpoint=True)
    ax1.set_xticks(np.linspace(0, len(m_profile), num=10, endpoint=True))
    ax1.set_xticklabels([str(int(round(i))) for i in xtick])
    
    labels = ax1.get_xticklabels() + ax1.get_yticklabels()
    [label.set_fontname('serif') for label in labels]
   
    ax1.legend(loc='upper left', prop=font1)
    
    #save figure
    fig_name = 'Profile_' + m_name + '_' + s_name + '.png'
    fig_output = outdir + fig_name
    fig.savefig(fig_output, dpi=300, bbox_inches='tight')

def distance_2points(lon_s, lat_s, lon_e, lat_e):
    """calculate distance[km] between two points based on their lat/lon"""
    # the radius of earth
    R = 6373.0
    lon_s_rad = math.radians(lon_s)
    lat_s_rad = math.radians(lat_s)
    lon_e_rad = math.radians(lon_e)
    lat_e_rad = math.radians(lat_e)

    dlon = lon_e_rad - lon_s_rad
    dlat = lat_e_rad - lat_s_rad

    a = math.sin(dlat / 2)**2 + math.cos(lat_s_rad) * math.cos(lat_e_rad) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    distance = R * c
    return distance
###################################### Beginning of CGPS class #################################### 
class Profile:
    """Profile class for profiles
    Example:
      import matplotlib.pyplot as plt
      from mimtpy.objects.profiles import Profile
      import mimtpy.objects.profiles as profiles
      pro_no, pro_lons, pro_lats = profiles.search_profiles(3,over_lat0, over_lon0, over_lat1, over_lon1, m_atr, s_atr)
      pro_obj = Profile(1, pi/6, pro_lon, pro_lat, m_overlay, s_overlay, over_lat0, over_lon0, m_atr, s_atr, outdir='$SCRATCHDIR/project/'+pro_No)
      pro_obj.profile_extract()
    """

    def __init__(self, NO, pro_angle, pro_lon, pro_lat, m_overlay, s_overlay, over_lat0, over_lon0, m_atr, s_atr, outdir):
        self.NO = NO
        self.pro_angle = pro_angle
        self.pro_lon = pro_lon
        self.pro_lat = pro_lat
        self.m_overlay = m_overlay
        self.s_overlay = s_overlay
        self.over_lat0 = over_lat0
        self.over_lon0 = over_lon0
        self.m_atr = m_atr
        self.s_atr = s_atr
        self.outdir = outdir
        
    def lonlat2rowcolm(self):
        """transfer lon/lat of user selected points to local coordinate"""
        lat_step = float(self.m_atr['Y_STEP'])
        lon_step = float(self.m_atr['X_STEP'])
        self.row_y = int((self.pro_lat - self.over_lat0) / lat_step + 0.5)
        self.colm_x = int((self.pro_lon - self.over_lon0) / lon_step + 0.5)

        return self.row_y, self.colm_x

    def rowcolm2ll(self, row, colm):
        """transfer row/colm to geo coordinate"""
        lat_step = float(self.m_atr['Y_STEP'])
        lon_step = float(self.m_atr['X_STEP'])
   
        lon = self.over_lon0 + lon_step * colm
        lat = self.over_lat0 + lat_step * row

        return lon, lat

    def profile_gmt(self, name):
        """write lon/lat of two ends of profile into gmt format"""
        gmt_file = self.outdir + 'profile_latlon_' + name + '.gmt'
    
        f = open(gmt_file, mode='w')
        f.write('# @VGMT1.0 @GLINESTRING \n')
        f.writelines(['# @R',str(min(self.lon_start,self.lon_end)),'/',str(max(self.lon_start,self.lon_end)),'/',str(min(self.lat_start,self.lat_end)),'/', str(max(self.lat_start,self.lat_end)),'\n'])
        f.write('# @Je4326 \n')
        f.write('# @Jp"+proj=longlat +datum=WGS84 +no_defs" \n')
        f.write('# @Jw"GEOGCS[\\"WGS 84\\",DATUM[\\"WGS_1984\\",SPHEROID[\\"WGS 84\\",6378137,298.257223563,AUTHORITY[\\"EPSG\\",\\"7030\\"]],AUTHORITY[\\"EPSG\\",\\"6326\\"]],PRIMEM[\\"Greenwich\\",0,AUTHORITY[\\"EPSG\\",\\"8901\\"]],UNIT[\\"degree\\",0.0174532925199433,AUTHORITY[\\"EPSG\\",\\"9122\\"]],AXIS[\\"Latitude\\",NORTH],AXIS[\\"Longitude\\",EAST],AUTHORITY[\\"EPSG\\",\\"4326\\"]]" \n')
        f.write('# @NId \n')
        f.write('# @Tinteger \n')
        f.write('# FEATURE_DATA \n')
        f.write('>')
        f.write('# @D0 \n')
        f.writelines([str(self.lon_start), ' ', str(self.lat_start), '\n'])
        f.writelines([str(self.lon_end), ' ' , str(self.lat_end), '\n'])
        f.close()
 
        return        

    def profile_json(self, name):
        """save each profile geometry info and data along value as json file"""
        profile_data = {"lon_start": self.lon_start, "lat_start": self.lat_start, "lon_end": self.lon_end, "lat_end": self.lat_end,
                        "m_data": self.m_profile.tolist(), "s_data": self.s_profile.tolist()}
        open(self.outdir + name + '.json', "w").write(json.dumps(profile_data))

    def profile_extract(self):
        """extract m_overlay and s_overlay data along the profile"""
        # get the size of overlay region
        rows, colms = np.shape(self.m_overlay)
        # get the origin position of the overlay region.
        self.lonlat2rowcolm()
        row_y = self.row_y
        colm_x = self.colm_x

        # calculat the intersect pixels between overlay region and profile
        angle = self.pro_angle # the angle is rad not degree
        if angle >= (45 * np.pi / 180) and angle <= (135 * np.pi / 180):
            # use colm to calculate row
            colm_no = np.arange(colms)
            if angle != (90 * np.pi / 180):
                tan_value = -1 * math.tan(angle)
                row_no = np.ceil(((colm_no - colm_x) / tan_value) + row_y)
            elif angle == (90 * np.pi / 180):
                row_no = np.ceil(np.ones(colms) * row_y)
        else:
            # use row to calculate colm
            row_no = np.arange(rows)
            if angle != 0:
                tan_value = -1 * math.tan(angle)
                colm_no = np.ceil((row_no - row_y) * tan_value + colm_x)
            elif anlge == 0:
                colm_no = np.ceil(np.ones(rows) * colm_x)
       
        row_no = row_no.astype(dtype=np.int)
        colm_no = colm_no.astype(dtype=np.int) 
        self.m_profile = self.m_overlay[row_no, colm_no]    
        self.s_profile = self.s_overlay[row_no, colm_no] 

        # change zero value to np.nan
        self.m_profile[(self.m_profile == 0)] = np.nan
        self.s_profile[(self.s_profile == 0)] = np.nan

        # calaculate lat/lon for profiles of two tracks
        row_start = row_no[0]
        row_end = row_no[-1]
        colm_start = colm_no[0]
        colm_end = colm_no[-1]
        self.lon_start, self.lat_start = self.rowcolm2ll(row_start, colm_start)
        self.lon_end, self.lat_end = self.rowcolm2ll(row_end, colm_end)
        
        # get Sentinel track name for master and slave data
        m_name_tmp = self.m_atr['FILE_PATH']
        self.m_name= re.search('Sen([^/]+)/', m_name_tmp)[1]
        s_name_tmp = self.s_atr['FILE_PATH']
        self.s_name= re.search('Sen([^/]+)/', s_name_tmp)[1]
        
        # save lat/lon files in gmt format
        pro_name = self.m_name + '_' + self.s_name + '_' + str(self.NO)
        self.profile_gmt(pro_name)
        self.profile_json(pro_name)
        
        # plot master and slave data along one profile
        # profile_plot()
  
        return 
#################################### End of GPS-UNR class ####################################
