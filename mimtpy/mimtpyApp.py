#!/usr/bin/env python3
#################################################################
# mimtpyApp generate data from HDFEOS file obtained by mintpy   #
# Author: Lv Xiaoran                                            #
# Created: June 2020                                            #
#################################################################
import os
import argparse
from mintpy.utils import readfile, writefile,utils as ut
import mimtpy

import mimtpy.workflow 
from mimtpy.utils import multitrack_utilities as mu

###############################################################################
STEP_LIST = [
     'velcumu',
     'horzvert',
     'concatenation',
     'geodmod',
     'plot'
]

STEP_HELP = """ Command line options for steps processing with names are chosen from the following list:

{}
""".format(STEP_LIST[0:4])

EXAMPLE = """example:
   mimtpyApp.py <template>  # run with default and custom templates
   
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Generate data based on HDFEOS file and mimtpyApp.template',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('template',help="Template file for mimtpyApp.")
                        
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  

    return inps    

def velcumu(inpsdict):
    """generate velocity.h5 and velocity.tiff file or generate cumudisp_date1_date2.h5 and cumudisp_date1_date2_tiff."""
    if inpsdict['mimtpy.velcumu'] == 'yes':
        Dataset = inpsdict['mimtpy.velcumu.DataSet'] 
        Dataset = list(tuple([i for i in Dataset.split(',')]))
        for dataset in Dataset:
            print('\n***********************Start processing {} dataset***********************'.format(dataset))
            # go to dataset_dir
            dataset_dir = os.path.abspath(os.path.join(os.getenv('SCRATCHDIR'),dataset,'mintpy'))
            if not os.path.isdir(dataset_dir):
                raise Exception('Error! No such dir : {}'.format(dataset_dir))
            os.chdir(dataset_dir)
            print('\nGo to project dir:', dataset_dir)
    
            # create velocity or cumulative dir
            filetype = inpsdict['mimtpy.velcumu.type']
            outdir = os.path.abspath(os.path.join(dataset_dir,filetype))
            if not os.path.isdir(outdir):
                os.makedirs(outdir)
            print('\nthe output dir for {} is {}.\n'.format(filetype,outdir))

            # find the HDFEOS file name
            HDFEOS_file = mu.find_HDFEOS_fullname(dataset_dir)
            print('\nThe HDFEOS file is {}'.format(HDFEOS_file))
            atr = readfile.read_attribute(HDFEOS_file)
    
            # get the startdate and enddate
            startDate,endDate = mu.find_start_end_date(dataset_dir,inpsdict['mimtpy.velcumu.startDate'],inpsdict['mimtpy.velcumu.endDate'])
            if startDate == 'None':
                startDate = atr['START_DATE']
            if endDate == 'None':
                endDate = atr['END_DATE']
            
            # generate parameter list for HDFEOS_to_geotiff.py
            if inpsdict['mimtpy.velcumu.SNWE'] == 'None':
                if inpsdict['mimtpy.velcumu.mask'] == 'y':
                   scp_args = [HDFEOS_file, filetype, '--date', startDate+'_'+endDate, '--mask', '--outdir', outdir]
                else:
                   scp_args = [HDFEOS_file, filetype, '--date', startDate+'_'+endDate, '--outdir', outdir]
            else:
                SNWE = inpsdict['mimtpy.velcumu.SNWE']
                SNWE = list(tuple([float(i) for i in SNWE.split(',')]))
                if inpsdict['mimtpy.velcumu.mask'] == 'y':
                   scp_args = [HDFEOS_file, filetype, '--date', startDate+'_'+endDate, '--mask', '--bbox', SNWE, '--outdir', outdir]
                else:
                   scp_args = [HDFEOS_file, filetype, '--date', startDate+'_'+endDate, '--bbox', SNWE, '--outdir', outdir]

            scp_args = mu.seperate_str_byspace(scp_args)
            # run HDFEOS_to_geotiff.py
            print('HDFEOS_to_geotiff.py',scp_args)
            mimtpy.HDFEOS_to_geotiff.main(scp_args.split())   
    else:
        print('\nSkip velcumu process')
    
    return

def horzvert(inpsdict):
    """generate horzontal and vertical files."""
    if inpsdict['mimtpy.horzvert'] == 'yes':
        Dataset = inpsdict['mimtpy.horzvert.DataSet']
        Dataset = list(tuple([i for i in Dataset.split(',')]))
     
        asc_data = Dataset[0]
        asc_dir = os.path.abspath(os.path.join(os.getenv('SCRATCHDIR'),asc_data,'mintpy'))
        if not os.path.isdir(asc_dir):
            raise Exception('Error! No such dir : {}'.format(asc_dir))
               
        des_data = Dataset[1]
        des_dir = os.path.abspath(os.path.join(os.getenv('SCRATCHDIR'),des_data,'mintpy'))
        if not os.path.isdir(des_dir):
            raise Exception('Error! No such dir : {}'.format(des_dir))

        # create workdir
        workdir = inpsdict['mimtpy.horzvert.outdir']
        print('\nthe work dir is {}.\n '.format(workdir))

        files = inpsdict['mimtpy.horzvert.dataname']
        files = list(tuple([i for i in files.split(',')]))
        file_asc = files[0]
        file_des = files[1]

        file_type = file_asc.split('_')[0]
        # out put dir
        outdir = os.path.abspath(os.path.join(workdir,file_type)) + '/'
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        print('\nthe output dir is %s .\n' % outdir)
        
        # go to output dir
        os.chdir(outdir)

        # get the input ascending and descending dir+name 
        file_asc_data = asc_dir + '/' + file_type + '/' + file_asc + '.h5'
        file_des_data = des_dir + '/' + file_type + '/' + file_des + '.h5'
    
        # reference point
        refpoi_lalo = inpsdict['mimtpy.horzvert.referencepoint']
        if refpoi_lalo != 'None':
            refpoi_lalo = list(tuple(i for i in refpoi_lalo.split(',')))

        # spatial range in lat/lon format
        SNWE = inpsdict['mimtpy.horzvert.SNWE']
        SNWE = list(tuple(float(i) for i in SNWE.split(',')))
        
        # resolve ascending and descending to horz and vert data.
        azimuth = inpsdict['mimtpy.horzvert.azimuth']
        outname = inpsdict['mimtpy.horzvert.outname']
        if outname != 'None':
            outname = list(tuple(i for i in outname.split(',')))
        
        if refpoi_lalo != 'None' and azimuth != 'None' and outname != 'None':
            scp_args2 = [file_asc_data, file_des_data, '--bbox', SNWE, '--reference_point', refpoi_lalo, '--azimuth', azimuth, '--outname', outname, '--outdir', outdir] 
        elif refpoi_lalo != 'None' and azimuth != 'None' and outname == 'None':
            scp_args2 = [file_asc_data, file_des_data, '--bbox', SNWE, '--reference_point', refpoi_lalo, '--azimuth', azimuth, '--outdir', outdir] 
        elif refpoi_lalo != 'None' and azimuth == 'None' and outname != 'None':
            scp_args2 = [file_asc_data, file_des_data, '--bbox', SNWE, '--reference_point', refpoi_lalo, '--outname', outname, '--outdir', outdir] 
        elif refpoi_lalo != 'None' and azimuth == 'None' and outname == 'None':
            scp_args2 = [file_asc_data, file_des_data, '--bbox', SNWE, '--reference_point', refpoi_lalo, '--outdir', outdir] 
        elif refpoi_lalo == 'None' and azimuth != 'None' and outname != 'None':
            scp_args2 = [file_asc_data, file_des_data, '--bbox', SNWE, '--azimuth', azimuth, '--outname', outname, '--outdir', outdir] 
        elif refpoi_lalo == 'None' and azimuth != 'None' and outname == 'None':
            scp_args2 = [file_asc_data, file_des_data, '--bbox', SNWE, '--azimuth', azimuth, '--outdir', outdir]
        elif refpoi_lalo == 'None' and azimuth == 'None' and outname != 'None':
            scp_args2 = [file_asc_data, file_des_data, '--bbox', SNWE, '--outname', outname, '--outdir', outdir]
        elif refpoi_lalo == 'None' and azimuth == 'None' and outname == 'None':
            scp_args2 = [file_asc_data, file_des_data, '--bbox', SNWE, '--outdir', outdir]
        scp_args2 = mu.seperate_str_byspace(scp_args2)
        print('\ngenerate_horzvert.py',scp_args2)
        mimtpy.generate_horzvert.main(scp_args2.split())
    else:
        print('\nSkip horzvert step')
 
    return

def concatenation(inpsdict):
    """concatenation adjencant tracks."""
    if inpsdict['mimtpy.concatenation'] == 'yes':
        Dataset = inpsdict['mimtpy.concatenation.DataSet']
        Dataset = list(tuple([i for i in Dataset.split(',')]))
        
        project_master = Dataset[0]
        promaster_dir = os.path.abspath(os.path.join(os.getenv('SCRATCHDIR'),project_master,'mintpy'))
        if not os.path.isdir(promaster_dir):
            raise Exception('Error! No such dir : {}'.format(promaster_dir))
        
        project_slave = Dataset[1]
        proslave_dir = os.path.abspath(os.path.join(os.getenv('SCRATCHDIR'),project_slave,'mintpy'))
        if not os.path.isdir(proslave_dir):
            raise Exception('Error! No such dir : {}'.format(proslave_dir))
        
        # create workdir
        workdir = inpsdict['mimtpy.concatenation.outdir']
        if not os.path.isdir(workdir):
            os.makedirs(workdir)
        print('\nthe work dir is {}.\n'.format(workdir))
        
        # run track_offset.py
        # file type
        files = inpsdict['mimtpy.concatenation.dataname']
        files = list(tuple([i for i in files.split(',')]))
        file_master = files[0]
        file_slave = files[1]
       
        file_type = file_master.split('_')[0]
        # out put dir
        outdir = os.path.abspath(os.path.join(workdir,file_type)) + '/'
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        print('the output dir is %s .\n' % outdir)
        
        #track_offset.py option
        rewrite_opt = inpsdict['mimtpy.concatenation.rewrite']
        plotpair_opt = inpsdict['mimtpy.concatenation.plotpair']

        # file_dir
        fmaster_dir = os.path.abspath(os.path.join(promaster_dir,file_type)) + '/'
        if not os.path.isdir(fmaster_dir):
            raise Exception('Error! No such dir : {}'.format(fmaster_dir))
        
        fslave_dir = os.path.abspath(os.path.join(proslave_dir,file_type)) + '/'
        if not os.path.isdir(fslave_dir):
            raise Exception('Error! No such dir : {}'.format(fslave_dir))
        
        # master and slave data name
        Dmaster = fmaster_dir + file_master + '.h5'
        Dslave = fslave_dir + file_slave + '.h5'
        
        # mosaic file name
        if inpsdict['mimtpy.concatenation.outname'] != 'None':
            outname = inpsdict['mimtpy.concatenation.outname']
        else:
            outname = file_master + '_mosaic'

        if rewrite_opt == 'y':
            if plotpair_opt == 'y':
                azi_angle = float(inpsdict['mimtpy.concatenation.azimuth'])
                scp_args = [Dmaster, Dslave, '--rewrite_slave', '--output', outname, '--plotpair', '--azi_angle', azi_angle, '--outdir', outdir]
            elif plotpair_opt == 'n':
                scp_args = [Dmaster, Dslave, '--rewrite_slave', '--output', outname, '--outdir', outdir]
        elif rewrite_opt == 'n':
            if plotpair_opt == 'y':
                azi_angle = float(inpsdict['mimtpy.concatenation.azimuth'])
                scp_args = [Dmaster, Dslave, '--output', outname, '--plotpair', '--azi_angle', azi_angle, '--outdir', outdir]
            elif plotpair_opt == 'n':
                scp_args = [Dmaster, Dslave, '--output', outname, '--outdir', outdir]
        scp_args = mu.seperate_str_byspace(scp_args)
        
        # run track_offset.py
        print('track_offset.py',scp_args)
        mimtpy.track_offset.main(scp_args.split())   
        
        # run H5UNW_to_geotiff.py
        # mosaic data
        Dmosaic = outdir + outname + '.h5'
        scp_args2 = [Dmosaic, '--outdir', outdir, '--output', outname + '.tiff'] 
        scp_args2 = mu.seperate_str_byspace(scp_args2)
        
        print('\nH5UNW_to_geotiff.py',scp_args2)
        mimtpy.H5UNW_to_geotiff.main(scp_args2.split())
    else:
        print('\nSkip concatenation process')       
    return

def plot(inpsdict):
    """plot."""
    if inpsdict['mimtpy.plot'] == 'yes':
        print('Start ploting data')
        plot_type = inpsdict['mimtpy.plot.type']
        if plot_type == 'velocity' or plot_type == 'displacement':
            Dataset = inpsdict['mimtpy.velcumu.DataSet'] 
            Dataset = list(tuple([i for i in Dataset.split(',')]))
            for dataset in Dataset:
                # go to dataset_dir
                dataset_dir = os.path.abspath(os.path.join(os.getenv('SCRATCHDIR'),dataset,'mintpy'))
                if not os.path.isdir(dataset_dir):
                    raise Exception('Error! No such dir : {}'.format(dataset_dir))
                os.chdir(dataset_dir)
                # go to data to be plotted dir
                plot_dir = os.path.abspath(os.path.join(dataset_dir,plot_type)) + '/'
                if not os.path.isdir(plot_dir):
                    raise Exception('Error! No such dir : {}'.format(plot_dir))
                print('plotting {} data in {} dir'.format(plot_type,plot_dir))
        
                print('Go to %s.\n' % plot_dir)
                os.chdir(plot_dir)
        
                # plot data
                Vfiles = []
                path_list = os.listdir(plot_dir)
                for Vfile in path_list:
                    Vfiles.append(Vfile)
                Vfiles.sort()
                
                for Vfile in Vfiles:
                    if Vfile.find(plot_type) >=0 and Vfile.find('.tiff')>0:
                        filename = Vfile
                        outfile = filename.split('.')[0] 
                        shpdir = inpsdict['mimtpy.plot.shpdir']
                        if shpdir != 'no':
                            faults = inpsdict['mimtpy.plot.fault']
                            faults = list(tuple(i for i in faults.split(',')))
                            fcolor = inpsdict['mimtpy.plot.fcolor']
                            fcolor = list(tuple(i for i in fcolor.split(',')))
                            fstyle = inpsdict['mimtpy.plot.fstyle']
                            fstyle = list(tuple(i for i in fstyle.split(',')))
                            refpoi = inpsdict['mimtpy.plot.refpoi']
                            scp_args = [filename, '--shpdir', shpdir, '--fault', faults, '--fcolor', fcolor, '--fstyle', fstyle, '--refpoi', refpoi, '--outfile', outfile, '--outdir', plot_dir]
                        else:
                            scp_args = [filename, '--outfile', outfile, '--outdir', plot_dir]

                        scp_args = mu.seperate_str_byspace(scp_args)
                        # run plot_geotiff.py
                        print('\n*********************ploting %s file********************' % Vfile)
                        print('plot_geotiff.py',scp_args)
                        mimtpy.plot_geotiff.main(scp_args.split())   
    else:
        print('\nSkip plot process')        

    return
    
def geodmod(inpsdict):
    """prepare data for geodmod software"""
    if inpsdict['mimtpy.geodmod'] == 'yes':
        Dataset = inpsdict['mimtpy.geodmod.DataSet'] 
        Dataset = list(tuple([i for i in Dataset.split(',')]))
        for dataset in Dataset:
            print('\n***********************Start processing {} dataset***********************'.format(dataset))
            # go to dataset_dir
            dataset_dir = os.path.abspath(os.path.join(os.getenv('SCRATCHDIR'),dataset,'mintpy'))
            if not os.path.isdir(dataset_dir):
                raise Exception('Error! No such dir : {}'.format(dataset_dir))
            print('\nGo to project dir:', dataset_dir)
            os.chdir(dataset_dir)
    
            # create geodmod dir
            software = 'geodmod'
            outdir = os.path.abspath(os.path.join(dataset_dir,software))
            if not os.path.isdir(outdir):
                os.makedirs(outdir)

            # find the HDFEOS file name
            HDFEOS_file = mu.find_HDFEOS_fullname(dataset_dir)
            print('\nThe HDFEOS file is {}'.format(HDFEOS_file))
            atr = readfile.read_attribute(HDFEOS_file)
    
            # get the startdate and enddate
            startDate,endDate = mu.find_start_end_date(dataset_dir,inpsdict['mimtpy.geodmod.startDate'],inpsdict['mimtpy.geodmod.endDate'])
            if startDate == 'None':
                startDate = atr['START_DATE']
            if endDate == 'None':
                endDate = atr['END_DATE']
           
            # go to geodmod dir
            print('\nGo to geodmod dir:', outdir)
            os.chdir(outdir)

            # generate parameter list for HDFEOS_to_geotiff.py
            if inpsdict['mimtpy.geodmod.SNWE'] == 'None':
                scp_args = ['../' + HDFEOS_file, '-s', startDate, '-e', endDate]
            else:
                SNWE = inpsdict['mimtpy.geodmod.SNWE']
                SNWE = list(tuple([float(i) for i in SNWE.split(',')]))
                scp_args = ['../' + HDFEOS_file, '-s', startDate, '-e', endDate, '--bbox', SNWE]

            scp_args = mu.seperate_str_byspace(scp_args)
            # run save_geodmod.py
            print('save_geodmod.py',scp_args)
            mimtpy.save_geodmod.main(scp_args.split())   
    else:
       print('\nSkip geodmod step')

def run(inpsdict, steps=STEP_LIST):
    """run the chosen steps."""
    for sname in steps:
        print('\n**************************step - {}***********************'.format(sname))
        if sname == 'velcumu':
            velcumu(inpsdict)
        elif sname == 'horzvert':
            horzvert(inpsdict)
        elif sname == 'concatenation':
            concatenation(inpsdict)
    # plot results
        elif sname == 'plot': 
            plot(inpsdict)  
        elif sname == 'geodmod':
            geodmod(inpsdict)
    return        

######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    inpsdict =  mu.read_template(inps.template)
    run(inpsdict)

######################################################################################
if __name__ == '__main__':
    main()
