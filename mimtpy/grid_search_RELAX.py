#!/usr/bin/env python3
#############################################################
# Grid search algorithm for viscosity based on RELAX        #
# Author: Lv Xiaoran Feb 2020                               #
#############################################################


import os
import argparse
import copy
import numpy as np
import matplotlib.pyplot as plt
from mintpy.utils import utils as ut
import json
from mimtpy.utils import multitrack_utilities
###########################################################################################
EXAMPLE = """example:
    For ascending: 
    grid_search_RELAX.py points_disp.json points_angle.json --insar-quality 0.00028177 19166.7435 0.00009112 --modeldata /data/lvxr/MODEL/RELAX/Relax/examples/Iran/result/ --outdir /data/lvxr/MODEL/RELAX/Relax/examples/Iran/result/
    For descending:
    grid_search_RELAX.py points_disp_06.json points_angle.json --insar-quality 0.000080557 15838.3728 0.000021162 --modeldata /data/lxrtest/MODELOUT/RELAX/Relax/examples/Iran/result/ --outdir /data/lxrtest/MODELOUT/RELAX/Relax/examples/Iran/result/
"""


def create_parser():
    parser = argparse.ArgumentParser(description='grid search algorithm for optimal viscosity',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)
    parser.add_argument('file', nargs='+',
                        help='LOS displacement Observations\n')
    parser.add_argument('geometry',nargs=1,
                        help='incidence angle and azimuth angle of observation points\n')
    
    parser.add_argument('--epi-lalo',type=float,dest='epilalo',metavar=('LAT','LON'), nargs=2,
                        help='lat and lon of epicenter\n')

    parser.add_argument('--insar-quality',type=float,dest='quality',metavar=('sillExp','range','nugget'),nargs=3,
                        help='sill, range and nugget of insar data.The unit of range is meter.\n')

    parser.add_argument('--modeldata',type=str, dest='mdata', nargs=1,
                       help='model data dir.\n')

    parser.add_argument('--rms', action='store_true',default=False)
    
    parser.add_argument('--outdir',type=str,dest='outdir',nargs=1, default=os.getenv('PWD'),
                        help='out put dir. The default value is $pwd\n')
    
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    return inps

###########################################################################################
def LOS_calculation(azi_angle,inc_angle,simulation,NotNan):
    """calcualte simulated LOS deformation"""
    # azi and inc are n*1 np.array
    #read output file
    inc = copy.deepcopy(inc_angle)
    east_disp = simulation[:,3][NotNan]
    north_disp = simulation[:,2][NotNan]
    up_disp = -simulation[:,4][NotNan]
    # calculate LOS
    #Project displacement from LOS to Horizontal and Vertical components
    #    math for 3D: cos(theta)*Uz - cos(alpha)*sin(theta)*Ux + sin(alpha)*sin(theta)*Uy = Ulos
    #    math for 2D: cos(theta)*Uv - sin(alpha-az)*sin(theta)*Uh = Ulos   #Uh_perp = 0.0
    inc *= np.pi/180.

    # heading angle
    head_angle = ut.azimuth2heading_angle(azi_angle)
    #if head_angle < 0.:
    #    head_angle += 360.
    head_angle[head_angle<0.]+= 360.
    head_angle *= np.pi/180.

    # construct design matrix
    A_up = np.cos(inc)
    A_east = - np.sin(inc) * np.cos(head_angle)
    A_north = np.sin(inc) * np.sin(head_angle)
    
    # LOS simulated results. Note:the unit of RELAX displacement is meter.
    los_sim = up_disp * A_up + north_disp * A_north + east_disp * A_east
    
    return los_sim
     
#def calculate_residual(observation,azi,inc,simulation,epi_lat,epi_lon): 
def calculate_residual(observation,azimuth,incidence,simulation,inv_covariance,NotNan,Vfile):
    """calculate residual between observation and simulation (LOS direction)"""
    #method 1 for residual: calculated residual based on the distance between points and epicenter
    obser_disp = observation[:,2]
    obser_disp = obser_disp[NotNan]
    sim = LOS_calculation(azimuth,incidence,simulation,NotNan)
    diff_tmp = obser_disp - sim
    diff = np.array([diff_tmp])
    #method 1 for residual: calculated residual based on the distance between points and epicenter
    # residual
    #row = len(diff)
    #delta_lat = (obser_lat - epi_lat) * (obser_lat - epi_lat)
    #delta_lon = (obser_lon - epi_lon) * (obser_lon - epi_lon)
    #delta_total = np.sqrt(delta_lat + delta_lon)
    #delta_weight = (delta_total - np.min(delta_total))/(np.max(delta_total) - np.min(delta_total))
    #multipli = delta_weight * (diff * diff)
    #summ = sum(multipli)
    #error = np.sqrt(summ/(row - 2))

    # method 2 for residual: residual = (observation - simulate)T * inv_covariance * (observation - simulate)
    # the unit for observation and simulate should be meter
    residual = np.dot(np.dot(diff,inv_covariance),np.transpose(diff)) 
    
    return residual

def calculate_rms(residual,sample_number):
    """calculate rms"""
    rms = np.sqrt(residual / sample_number)
    print(rms)
    return rms

def observation_covariance(file_name,NotNan,insar_sill,insar_range,insar_nugget):
    """create inverse of covariance"""
    print('The observation file name is %s'%file_name)
    data_json = json.loads(open(file_name).read())
    data_value = np.array(data_json["lalo_value"])
    lat = data_value[:,0]
    lon = data_value[:,1]
    lat = lat[NotNan]
    lon = lon[NotNan]
    reference_lat = float(data_json["Ref_lat"])
    reference_lon = float(data_json["Ref_lon"])
    reference = np.array([reference_lon,reference_lat],dtype=float)
    
    points_num = len(lat)
    llh = np.zeros((2,points_num),dtype=float)
    llh[0,:] = np.transpose(lon)
    llh[1,:] = np.transpose(lat)
    # the unit of x,y is meter
    xy = np.transpose(multitrack_utilities.llh2xy(llh,reference))*1000
    # calculate the distance between points
    x1,x2 = np.meshgrid(xy[:,0],xy[:,0])
    y1,y2 = np.meshgrid(xy[:,1],xy[:,1]) 
    H = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

    # calculate covariance
    covarianceMatrix = insar_sill * np.exp(-H / insar_range) + insar_nugget * np.identity(len(lat))
    inv_covariance = np.linalg.inv(covarianceMatrix)

    return inv_covariance

def read_files(file_name):
    """ read json data """
    data_json = json.loads(open(file_name).read())
    #print(data_json)
    if str.find(os.path.split(file_name)[1],'vs') != -1:
        data = np.array(data_json["displacement"])
    else: 
        data = np.array(data_json["lalo_value"])
    return data

def plot_residual(residual,obser_number):
    """plot residual"""
    print("****************************************\n")
    print("Start drawing results\n")
    #residual[:,2] = np.sqrt(residual[:,2]/obser_number)
    x_num = len(np.unique(residual[:,0]))
    y_num = len(np.unique(residual[:,1]))
    x = np.linspace(np.min(residual[:,0]),np.max(residual[:,1]),x_num)
    y = np.linspace(np.min(residual[:,1]),np.max(residual[:,1]),y_num)
    value = residual[:,2].reshape(x_num,y_num)
    X,Y = np.meshgrid(x,y)
    im = plt.contourf(X,Y,value,alpha=0.75,cmap=plt.cm.jet)
    cbar=plt.colorbar(im)
    contour = plt.contour(X, Y, value, 8, colors = 'black')
    plt.clabel(contour, fontsize=6, colors='k',fmt='%.1f')
    #cbar.ax.tick_params(labelsize=10)
    cbar.set_label('rms')
    #cbar.set_ticks(np.linspace(0,1,10))
    # set the range of legend
    #cbar.set_ticks(np.linspace(0,1,50))
    plt.title("InSAR grid search") 
    plt.xlabel('Upper mantle log(viscosity)')
    plt.ylabel('Lower crust log(viscosity)')
    # C = plt.contour(X, Y, value, 8, colors = 'black', linewidth = 0.5)
    #plt.clabel(C, inline = True, fontsize = 10)
    plt.savefig('gird_search_RELAX.png', dpi=300, bbox_inches='tight')    

def write_json(data,name):
    """write json file"""
    residual = {"lc_um_residual": data.tolist()}
    open(name + '.json', "w").write(json.dumps(residual))

###########################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    obser_dir = os.getenv('PWD')
    # epi_lat epi_lon
    #epi_lat = float(inps.epilalo[0])
    #epi_lon = float(inps.epilalo[1])

    #read head_angle and azimuth_angle
    geo_file = "".join(inps.geometry)
    obs_file = "".join(inps.file)
    data_angle = read_files(geo_file)
    azimuth = data_angle[:,3]
    incidence = data_angle[:,2]

    # read observation
    obser_data = read_files(obs_file)
    obser_did = obser_data[:,2]
    
    # read one of simulation data
    sim_example = "".join(inps.mdata) + '/' + 'vs_17001700.json'
    sim_example_data = read_files(sim_example)
    sim_example_north = sim_example_data[:,2]
    NotNan_sim = ~np.isnan(sim_example_north)

    # get the position of Nan in observation data 
    NotNan_obs = ~np.isnan(azimuth)
   
    # combine the not nan position index
    NotNan = NotNan_sim * NotNan_obs

    # judge the consistence between geo_angle and displacement
    obser_disp_num = np.count_nonzero(obser_did != obser_did)
    obser_angle_num = np.count_nonzero(azimuth != azimuth)
    if obser_disp_num != obser_angle_num:
        parser.print_usage()
        raise Exception('Error! The mask for geo_angle and displacement is different!')
    
    # get array of observation/lat/lon
    azimuth = azimuth[NotNan]
    incidence = incidence[NotNan]
    obser_number = len(azimuth)
    print("sample number is %d" % obser_number)
    # calculate inv_covariance for observation data
    inv_covariance = observation_covariance(obs_file,NotNan,float(inps.quality[0]),float(inps.quality[1]),float(inps.quality[2]))
    
    #start grid search
    
    # change to RELAX dir
    model_path = "".join(inps.mdata) + '/'
    os.chdir(model_path)
   
    residual = np.empty(shape=[0,3],dtype=float) 
    
    #search vs**.json files
    Vfiles = []
    path_list = os.listdir(model_path)
    for Vfile in path_list:
        Vfiles.append(Vfile)
    Vfiles.sort()
    
    for Vfile in Vfiles:
        print('process %s file' % Vfile)  
        viscosity1 = float(str(Vfile.split('.')).split('_')[1][0:4])
        viscosity2 = float(str(Vfile.split('.')).split('_')[1][4:8])
        post_disp = read_files(Vfile)
        error = calculate_residual(obser_data,azimuth,incidence,post_disp,inv_covariance,NotNan,Vfile)[0,0]

        print('viscosity1 and viscosity2 are %.2f,%.2f' % (viscosity1/100,viscosity2/100))
        if inps.rms:
            print('using rms value:')
            rms = calculate_rms(error,obser_number)
            viscosity_error = np.array([viscosity1/100,viscosity2/100,rms])
        else:    
            viscosity_error = np.array([viscosity1/100,viscosity2/100,error])
        residual = np.append(residual,[viscosity_error],axis=0)
    print(residual)
    # write residual file
    # change to observation dir
    os.chdir(obser_dir) 
    
    # store the residual file
    file_name = 'grid_search_residual_RELAX'
    write_json(residual,file_name)
    #with open(file_name,'w') as f:
    #    for line in residual:
    #        f.write('%s' % line) 
    
    # draw gird_search picture
    plot_residual(residual,obser_number)    
#########################################################################################
if __name__ == '__main__':
    main()
