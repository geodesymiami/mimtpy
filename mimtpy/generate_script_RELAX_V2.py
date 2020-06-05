#!/usr/bin/env python3
#########################################################################
# Program is used for generate template files for RELAX                 #
# Author: Lv Xiaoran                                                    #
# Created: December 2019                                                #
#########################################################################

import os
import argparse
import math
import numpy as np
import copy
######################################################################################
EXAMPLE = """example:
    generate_script_RELAX_V2.py --sample_file /data/lvxr/MODEL/RELAX/Relax/examples/Iran/iran_vs_17251725.sh --parameters_number 2 --parameters_symbol V1,V2 --parameters_range 17,0.25,18,18,0.25,19 --template_dir /data/lvxr/MODEL/RELAX/Relax/examples/Iran/no_PBS_512/
    
    generate_script_RELAX_V2.py --sample_file /data/lvxr/MODEL/RELAX/Relax/examples/Iran/iran_vs_17251725.sh --parameters_number 2 --parameters_symbol V1,V2 --parameters_range 17,0.25,18,18,0.25,19 --template_dir /data/lvxr/MODEL/RELAX/Relax/examples/Iran/no_PBS_512/ --PBS

"""

def create_parser():
    parser = argparse.ArgumentParser(description='Generate template files for RELAX.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('--sample_file', dest='sample_file', nargs=1, type=str, help='template file directory for RELAX\n')

    parser.add_argument('--parameters_number', dest='para_num', nargs=1, type=int, help='the number of parameters\n')

    parser.add_argument('--parameters_symbol',dest='para_symbol',nargs='+', type=str, help='the symbol of parameters.Please using "V**" for viscosity, using "H**" for depth\n')

    parser.add_argument('--parameters_range', dest='para_range', nargs='+', help='the range of all parameters\n')

    parser.add_argument('--template_dir', dest='template_dir',nargs=1,help='template files dir.\n')

    parser.add_argument('--PBS', action='store_true', default=False)

    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    return inps

def vis_to_gommadot0(viscosity):
    """transfer viscosity to gommadot0"""
    # change log(viscosity) to viscosity
    vis = math.pow(10,viscosity)
    # first transfer to Maxwell time
    mu = 3.0E+10
    second_to_year = 3.2E+07
    Maxwell_time = (vis / mu) / second_to_year
    # then transfer to gommadot0
    gommadot0 = 1 / Maxwell_time
    return gommadot0

def judge_keys(line,para_dict):
    """judge whether the parameter is in line"""
    position_flag = -1
    key_flag = 'None'
    for key in para_dict:
        key_judge = key + '='
        position = line.find(key_judge)
        if position >= 0:
            position_flag = position
            key_flag = key

    return position_flag, key_flag


def change_text(inps, text_file, template_dir, para_dict,parameters):
    """change text_before to text_after"""
    with open(text_file,'r') as f:
        lines = []
        for line in f.readlines():
            if line != '\n':
                lines.append(line)
    f.close()

    sample_file_name = os.path.split(inps.sample_file[0])[-1]
    text_file_1 = sample_file_name.split('.')[0].split('_')[0]
    text_file_2 = sample_file_name.split('.')[0].split('_')[1]
    tail = ''
    for ppara in parameters:
        key_value = np.str("{:.2f}".format(ppara))
        if key_value.find('.') >0:
            temp = key_value.split('.')[0] + key_value.split('.')[1] + '_'
        else:
            temp = key_value + '_'
        tail += temp

    tail = tail[0:-1]

    text_file_change = text_file_1 + '_' + text_file_2 + '_' + tail

    line_num = 0
    for line in lines:
        position, key = judge_keys(line,para_dict)
        if position >= 0:
            lines[line_num] = key + '=' + str(para_dict[key]) + '\n'
        if line.find('WDIR=') >= 0:
            lines[line_num] = "WDIR=../" + text_file_change + "\n"
        line_num += 1

    # PBS settings
    if inps.PBS:
        line_numm = 0
        for line in lines:
            if line.find('#PBS -N') >= 0:
                lines[line_numm] = '#PBS -N ' + text_file_change + '\n'
            if line.find('#PBS -o') >= 0:
                lines[line_numm] = '#PBS -o ../' + text_file_change + '/out' + '\n'
            line_numm += 1

    template_file_name = text_file_change + '.sh'

    with open(template_dir + template_file_name,'w') as f:
        for line in lines:
            f.write('%s' % line)
    f.close()

def initial_parameters(inps):
    """initial parameters"""
    # initial parameters
    para_num = inps.para_num[0]
    para_symbol = "".join(inps.para_symbol)
    para_range = "".join(inps.para_range)
    para_range = list(tuple(float(i) for i in para_range.split(',')))
    para_symbol = list(tuple(i for i in para_symbol.split(',')))
    para_dict = {}
    for i in np.arange(para_num):
        psymbol = para_symbol[i]
        prange = para_range[3*i:3*i+3]
        pstart = prange[0]
        pstep = prange[1]
        pend = prange[2]
        pvalue_range = np.arange(pstart, pend+pstep, pstep)
        para_dict[psymbol] = pvalue_range

    return para_dict

def set_parameters(inps,parameters):
    """set parameters' value"""
    para_num = inps.para_num[0]
    para_dict_v2 = {}
    para_symbol = "".join(inps.para_symbol)
    para_symbol = list(tuple(i for i in para_symbol.split(',')))
    for i in np.arange(para_num):
        psymbol = para_symbol[i]
        if psymbol.find('V') >= 0:
        # transfer log(viscosity) into gommadot0
            gommadot0 = round(vis_to_gommadot0(parameters[i]),5)
            para_dict_v2[psymbol] = gommadot0
        else:
            para_dict_v2[psymbol] = parameters[i]

    return para_dict_v2

def generate_template(inps,parameters):
    """generate_template based on the given parameters range"""
    para_dict_v2 = set_parameters(inps,parameters)
    #print(para_dict_v2)
    # generate template
    #sample_file_name = os.path.split(inps.sample_file[0])[-1]

    # change text
    change_text(inps,inps.sample_file[0],inps.template_dir[0],para_dict_v2,parameters)

def generate_runfile(inps, template_dir):
    """generate runfile for job_submission.py"""
    runfile_name = 'run_RELAX'
    # search for ".sh" file
    shfiles = []
    path_list = os.listdir(template_dir)
    for shfile in path_list:
        shfiles.append(shfile)
    shfiles.sort()

    if inps.PBS:
        with open(template_dir + runfile_name, 'w') as h:
            for files in shfiles:
                lines = 'qsub ./' + files + '\n'
                h.write('%s' % lines)
        h.close()
    else:
        with open(template_dir + runfile_name, 'w') as h:
            for files in shfiles:
                lines = './' + files + '\n'
                h.write('%s' % lines)
        h.close()

def gen_n_dem_seq(para_dict, curr = [], result = []):
    """ gen n demision sequence """
    para_dict = copy.deepcopy(para_dict)
    if len(para_dict) == 0:
        if len(curr) == 0:
            raise Exception('empty para_dict')
        result += [curr]
        return result
    # get the first key to current para_dict
    key = list(para_dict.keys())[0]
    seq = para_dict[key]
    del para_dict[key]
    for item in seq:
         new_curr = curr + [item]
         gen_n_dem_seq(para_dict, new_curr, result)
    return result

######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    # template file dir
    if inps.template_dir != None:
        template_dir = inps.template_dir[0]
        isExists = os.path.exists(template_dir)
        if not isExists:
            os.makedirs(template_dir)

    # initial parameters
    para_dict = initial_parameters(inps)
    print('parameter and their values:\n')
    print(para_dict)
    # generate the permutations for parameters
    para_list = gen_n_dem_seq(para_dict)
    # generate template files
    for parameters in para_list:
        print('\nGenerate template file for parameter combination:')
        print(parameters)
        generate_template(inps,parameters)

    # generate run.sh file to run all template
    generate_runfile(inps, inps.template_dir[0])
######################################################################################
if __name__ == '__main__':
    main()
