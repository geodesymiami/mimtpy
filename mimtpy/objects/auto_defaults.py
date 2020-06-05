## ALL path and default strings used in the program
# Author: Sara Mirzaee

import os
import datetime


class PathFind:
    def __init__(self):
        self.logdir = os.getenv('OPERATIONS') + '/LOGS'
        self.scratchdir = os.getenv('SCRATCHDIR')
        #self.required_template_options = ['topsStack.subswath', 'topsStack.boundingBox']
        self.required_template_options = []
        self.defaultdir = os.path.expandvars('${RSMAS_INSAR}/sources/mimt/defaults')
        self.processdir = 'PYSAR'
        self.auto_template = self.defaultdir + '/mimt_template.txt'
        return

    def set_isce_defaults(self, inps):

        inps_dict = vars(inps)

        inps_dict['template']['topsStack.demDir'] = inps.work_dir + '/DEM'

        if 'cleanopt' not in inps.template:
            inps_dict['template']['cleanopt'] = '0'

        inps_dict['template']['topsStack.workingDir'] = inps.work_dir

        return

    @staticmethod
    def grab_cropbox(inps):
        try:
            if inps.template['processingMethod'] == 'smallbaseline':
                subset = inps.template['mintpy.subset.lalo']
            else:
                subset = inps.template['minopy.subset']
            subset = subset.split(':')
            cropbox = '{} {} {} {}'.format(subset[0], subset[1].split(',')[0], subset[1].split(',')[1], subset[2])
        except:
            cropbox = inps.template['topsStack.boundingBox']

        return cropbox

    @staticmethod
    def correct_for_ssara_date_format(template_options):

        inps_dict = template_options

        if 'ssaraopt.startDate' in inps_dict:
            inps_dict['ssaraopt.startDate'] = \
                datetime.datetime.strptime(inps_dict['ssaraopt.startDate'], '%Y%m%d').strftime('%Y-%m-%d')

        if 'ssaraopt.endDate' in inps_dict:
            inps_dict['ssaraopt.endDate'] = \
                datetime.datetime.strptime(inps_dict['ssaraopt.endDate'], '%Y%m%d').strftime('%Y-%m-%d')

        return inps_dict


    @staticmethod
    def isce_clean_list():
        cleanlist = []
        cleanlist.append(['stack',  'misreg', 'orbits', 'coarse_interferograms', 'ESD',
                          'interferograms', 'slaves'])
        cleanlist.append(['merged', 'master', 'coreg_slaves', 'baselines', 'geom_master'])
        cleanlist.append(['SLC'])
        cleanlist.append(['MINTPY', 'run_files', 'configs', 'DEM'])

        return cleanlist

