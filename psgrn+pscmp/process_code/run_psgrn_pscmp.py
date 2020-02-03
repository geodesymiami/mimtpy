#!/usr/bin/env python3
#########################################################################
# Program is used for runing psgrn/pscmp software in linux environment  #
# Author: Lv Xiaoran                                                    #
# Created: December 2019                                                #
#########################################################################

import os
import sys
import argparse
import subprocess
######################################################################################
EXAMPLE = """example:
    run_psgrn_pscmp.py --psgrn_file $MIMTFILES/psgrn+pscmp/psgrn/psgrn08-model1.inp -module 0 --psgrn_dir $MODELOUT/psgrn/model3
    run_psgrn_pscmp.py --pscmp_file $MIMTFILES/psgrn+pscmp/pscmp/pscmp08-model1.inp -module 1 --psgrn_dir $MODELOUT/psgrn/model3/ --pscmp_dir $MODELOUT/pscmp/Bogd/model3
    run_psgrn_pscmp.py --psgrn_file $MIMTFILES/psgrn+pscmp/psgrn/psgrn08-model1.inp --pscmp_file $MIMTFILES/psgrn+pscmp/pscmp/pscmp08_model1.inp -module 2 --psgrn_dir $MODELOUT/psgrn/model3 --pscmp_dir $MODELOUT/pscmp/Bogd/model3
"""
# first please complie psgrn/pscmp fortran code
# gfortran psgrn2019-code/*f -o psgrn2019
# gfortran pscmp2019-code/*f -o pscmp2019
def create_parser():
    parser = argparse.ArgumentParser(description='Run psgrn/pscmp software',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('--psgrn_file', dest='psgrn_file', nargs='?', help='template file directory for psgrn\n')
    
    parser.add_argument('--pscmp_file', dest='pscmp_file', nargs='?', help='template file directory for pscmp\n')
    
    parser.add_argument('-module','--module',dest='module',type = int, help='choose module for run. 0: psgrn;' + 
                                                                '1: pscmp; 2: psgrn+pscmp. \n')
    parser.add_argument('--psgrn_dir', dest='psgrn_dir',nargs=1,help='green function dir. \n')

    parser.add_argument('--pscmp_dir', dest='pscmp_dir',nargs=1,help='displacmemnt files dir. \n')
    #parser.add_argument('-outdir','--outdir',dest='outdir',nargs='?',default=os.getenv('MODELOUT'),
    #                    help='output directory')

    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps    
    
def change_text(text_file, text_before, text_after, diff_lines):
    """change text_before to text_after"""
    with open(text_file,'r') as f:
        lines = []
        for line in f.readlines():
            if line != '\n':
                lines.append(line)
    f.close()
    
    line_num = 0
    for line in lines:
        if line == text_before:
            lines[line_num-diff_lines] = text_after
        line_num += 1
        
    with open(text_file,'w') as f:
        for line in lines:
            f.write('%s' % line)

def set_output_dir(inps):
    """set output dir in *.inp template file"""
    psgrn_dir = "".join(inps.psgrn_dir)
    pscmp_dir = "".join(inps.pscmp_dir)
    if inps.module == 0:
        text_before = " 'uz'  'ur'  'ut'\n"
        text_after = " '" + psgrn_dir + "/'\n"
        change_text(inps.psgrn_file, text_before, text_after, 1)
        
    if inps.module == 1:
        text_before1 = " 'uz'  'ur'  'ut'\n"
        text_before2 = "  'U_north.dat'    'U_east.dat'    'U_down.dat'\n"
        text_after1 = " '" + psgrn_dir + "/'\n"
        text_after2 = " '" + pscmp_dir + "/'\n"      
        change_text(inps.pscmp_file, text_before1, text_after1, 1)
        change_text(inps.pscmp_file, text_before2, text_after2, 2)
    if inps.module == 2:
        text_before1 = " 'uz'  'ur'  'ut'\n"
        text_before2 = "  'U_north.dat'    'U_east.dat'    'U_down.dat'\n"
        text_after1 = " '" + psgrn_dir + "/'\n"
        text_after2 = " '" + pscmp_dir + "/'\n"
        change_text(inps.psgrn_file, text_before1, text_after1, 1)        
        change_text(inps.pscmp_file, text_before1, text_after1, 1)
        change_text(inps.pscmp_file, text_before2, text_after2, 2)

def run_module(programpath, parameterfile):
    """ programpath, prameterfile """
    parameterfile_dir = os.path.dirname(parameterfile)
    print(parameterfile_dir)
    parameterfile_name = os.path.basename(parameterfile)
    print(parameterfile_name)
    os.chdir(parameterfile_dir)
    cmd = "cd %s && echo %s | %s" % (parameterfile_dir, parameterfile_name, programpath)
    return os.system(cmd)
    #cmd = "echo %s | %s" % (parameterfile_name, programpath)
    #process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    #while True:
    #    output = process.stdout.readline()
    #    if output =='' and process.poll() is not None:
    #        break
    #    if output:
    #       print(output.strip())
    #rc = process.poll()
    #return rc
    
######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    psgrn_dir = "".join(inps.psgrn_dir)
    pscmp_dir = "".join(inps.pscmp_dir)
    #whether outdir exists
    isExists = os.path.exists(psgrn_dir)
    if not isExists:
        os.makedirs(psgrn_dir)

    isExists = os.path.exists(pscmp_dir)
    if not isExists:
        os.makedirs(pscmp_dir)
    
    # change output_dir in template file as the given output_dir
    set_output_dir(inps)

    if inps.module == 0:
        module_name = os.getenv('MODELDIR') + '/GFZ/psgrn+pscmp/psgrn2019'
        print(module_name)
        template_file = inps.psgrn_file
        run_module(module_name, template_file)
    if inps.module == 1:
        module_name = os.getenv('MODELDIR') + '/GFZ/psgrn+pscmp/pscmp2019'
        template_file = inps.pscmp_file
        run_module(module_name, template_file)
    if inps.module ==2:
        module_grn_name = os.getenv('MODELDIR') + '/GFZ/psgrn+pscmp/psgrn2019'
        template_grn_file = inps.psgrn_file
        run_module(module_grn_name, template_grn_file)  
        
        module_cmp_name = os.getenv('MODELDIR') + '/GFZ/psgrn+pscmp/pscmp2019'
        template_cmp_file = inps.pscmp_file
        run_module(module_cmp_name, template_cmp_file)        

   
######################################################################################
if __name__ == '__main__':
    main()
