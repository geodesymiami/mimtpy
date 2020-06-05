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
     'plot',
     'geodmod',
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
        
        # get the input ascending and descending dir+name 
        file_asc_data = asc_dir + '/' + file_type + '/' + file_asc + '.h5'
        file_des_data = des_dir + '/' + file_type + '/' + file_des + '.h5'
    
        # judge whether we should change the reference point
        refpoi_lalo = inpsdict['mimtpy.horzvert.referencepoint']
        if refpoi_lalo != 'None':
            refpoi_lalo = list(tuple(float(i) for i in refpoi_lalo.split(',')))
            refpoi_lat = refpoi_lalo[0]
            refpoi_lon = refpoi_lalo[1]

            print('\nchanging ascending and descending data reference point')
            # change ascending data refpoi
            completion_status = os.system(mu.seperate_str_byspace(['reference_point.py', file_asc_data, '-l', refpoi_lat, '-L', refpoi_lon]))
            if completion_status == 1:
                raise Exception('error when runing reference_point.py for %s' % file_asc_data)
            asc_data, asc_atr = readfile.read(file_asc_data)
            asc_atr['REF_LAT'] = refpoi_lat
            asc_atr['REF_LON'] = refpoi_lon
            writefile.write(asc_data, out_file=file_asc_data, metadata=asc_atr)
            # change descending data refpoi 
            completion_status = os.system(mu.seperate_str_byspace(['reference_point.py', file_des_data, '-l', refpoi_lat, '-L', refpoi_lon]))
            if completion_status == 1:
                raise Exception('error when runing reference_point.py for %s' % file_des_data)
            des_data, des_atr = readfile.read(file_des_data)
            des_atr['REF_LAT'] = refpoi_lat
            des_atr['REF_LON'] = refpoi_lon
            writefile.write(des_data, out_file=file_des_data, metadata=des_atr)

        print('\ngo to the output dir {}'.format(outdir))
        os.chdir(outdir)
        
        # spatial range in lat/lon format
        SNWE = inpsdict['mimtpy.horzvert.SNWE']
        SNWE = list(tuple(float(i) for i in SNWE.split(',')))

        # subset ascending and descending to the same spatial region
        print('\nsubset ascending data')
        sub_file_asc = 'subset_' + file_asc + '.h5'
        file_asc_data = asc_dir + '/' + file_type + '/' + file_asc + '.h5'
        print(mu.seperate_str_byspace(['subset.py', file_asc_data, '-l', SNWE[0:2], '-L', SNWE[2:4], '-o', sub_file_asc]))
        completion_status = os.system(mu.seperate_str_byspace(['subset.py', file_asc_data, '-l', SNWE[0:2], '-L', SNWE[2:4], '-o', sub_file_asc]))
        if completion_status == 1:
            raise Exception('error when subset ascending data!')
 
        print('\nsubset descending data')
        sub_file_des = 'subset_' + file_des + '.h5'
        file_des_data = des_dir + '/' + file_type + '/' + file_des + '.h5'
        print(mu.seperate_str_byspace(['subset.py', file_des_data, '-l', SNWE[0:2], '-L', SNWE[2:4], '-o', sub_file_des]))
        completion_status = os.system(mu.seperate_str_byspace(['subset.py', file_des_data, '-l', SNWE[0:2], '-L', SNWE[2:4], '-o', sub_file_des]))
        if completion_status == 1:
            raise Exception('error when subset descending data!')

        # resolve ascending and descending to horz and vert data.
        azimuth = inpsdict['mimtpy.horzvert.azimuth']
        outname = inpsdict['mimtpy.horzvert.outname']
        if outname == 'None':
            horz_name = file_type + '_horz.h5'
            vert_name = file_type + '_vert.h5'
        else:
            outname = list(tuple(i for i in outname.split(',')))
            horz_name = file_type + outname[0]
            vert_name = file_type + outname[1]
        if azimuth == 'None':    
            print(mu.seperate_str_byspace(['asc_desc2horz_vert.py', sub_file_asc, sub_file_des, '-o', horz_name, vert_name]))
            completion_status = os.system(mu.seperate_str_byspace(['asc_desc2horz_vert.py', sub_file_asc, sub_file_des, '-o', horz_name, vert_name]))
            if completion_status == 1:
                raise Exception('error when running asc_desc2horz_vert.py!')
        else:
            print(os.system(mu.seperate_str_byspace(['asc_desc2horz_vert.py', sub_file_asc, sub_file_des, '--az', azimuth, '-o', horz_name, vert_name])))
            completion_status = os.system(mu.seperate_str_byspace(['asc_desc2horz_vert.py', sub_file_asc, sub_file_des, '--az', azimuth, '-o', horz_name, vert_name]))
            if completion_status == 1:
                raise Exception('error when running asc_desc2horz_vert.py!')
            
        # run H5UNW_to_geotiff.py
        # horizontal data
        scp_args2 = [horz_name, '--outdir', outdir, '--output', horz_name.split('.')[0] + '.tiff'] 
        scp_args2 = mu.seperate_str_byspace(scp_args2)
        print('\nH5UNW_to_geotiff.py',scp_args2)
        mimtpy.H5UNW_to_geotiff.main(scp_args2.split())

        # vertical data
        scp_args2 = [vert_name, '--outdir', outdir, '--output', vert_name.split('.')[0] + '.tiff']
        scp_args2 = mu.seperate_str_byspace(scp_args2)
        print('\nH5UNW_to_geotiff.py',scp_args2)
        mimtpy.H5UNW_to_geotiff.main(scp_args2.split())

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
        mosaic_opt = inpsdict['mimtpy.concatenation.mosaic']

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
            if mosaic_opt == 'y':
                scp_args = [Dmaster, Dslave, '--rewrite_slave', '--mosaic', '--output', outname, '--outdir', outdir]
            elif mosaic_opt == 'n':
                scp_args = [Dmaster, Dslave, '--rewrite_slave', '--output', outname, '--outdir', outdir]
        elif rewrite_opt == 'n':
            if mosaic_opt == 'y':
                scp_args = [Dmaster, Dslave, '--mosaic', '--output', outname, '--outdir', outdir]
            elif mosaic_opt == 'n':
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
        #elif plot_type == 'horzvert':
        #    plot_dir = inpsdict['mimtpy.horzvert.outdir']
        #elif plot_type == 'coneatenation':
        #    plot_dir = inpsdict['mimtpy.coneatenation.outdir']
        
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
                        fault = inpsdict['mimtpy.plot.fault']
                        refpoi = inpsdict['mimtpy.plot.refpoi']
                        vlim = inpsdict['mimtpy.plot.vlim']
                        outfile = filename.split('.')[0] 
                        if vlim == 'None':
                            scp_args = [filename, '--fault', fault, '--refpoi', refpoi, '--outfile', outfile, '--outdir', plot_dir]
                        else:
                            vlim = list(tuple([float(i) for i in vlim.split(',')]))
                            scp_args = [filename, '--fault', fault, '--refpoi', refpoi, '--vlim', vlim, '--outfile', outfile, '--outdir', plot_dir]

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
