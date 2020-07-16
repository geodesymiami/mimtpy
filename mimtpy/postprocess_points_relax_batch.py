#!/usr/bin/env python3
######################################################################
# Generate displacement for certain time period based on RELAX       #
# Author: Lv Xiaoran Mar 2020                                        #
######################################################################

import numpy as np
import argparse
import os
##############################################################
EXAMPLE = """examples:
   postprocess_points_relax_density_batch.py ./ --lalo ../radar/pos.lalo --minvis 17 --maxvis 21 --stepvis 0.25 --start 0.6 --end 2.3 --Tstep 0.1 --outname vs_ --outdir ../result_density/
"""  

def create_parser():
    parser = argparse.ArgumentParser(description='postprocess displacement for certain time period based on output of Relax',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('workdir', nargs=1, type=str,
                        help='work directory\n')
    parser.add_argument('--lalo', dest='lalo', nargs=1,type=str,
                        help='lalo information for observation points.\n')
    parser.add_argument('--minvis', nargs=1, dest='minvis', type=float,
                        help='the smallest value of the viscosity.\n')
    parser.add_argument('--maxvis', nargs=1, dest='maxvis', type=float,
                        help='the largest value of the viscosity.\n')
    parser.add_argument('--stepvis', nargs=1, dest='stepvis', type=float,
                        help='the step of viscosity.\n')
    parser.add_argument('--start', dest='stime', type=float, nargs=1,
                        help='start time of the studied postseismic time.\n')
    parser.add_argument('--end', dest='etime', type=float, nargs=1,
                        help='end time of the studied postseismic time.\n')
    parser.add_argument('--Tstep',dest='tstep',type=float, nargs=1, metavar=('Interval'),
                        help='time step.\n')
    parser.add_argument('--outname', dest='outname', type=str, nargs=1,help='outname.\n')
    parser.add_argument('--outdir', dest='outdir', nargs=1, type=str,
                        help='output dir\n')
    return parser

def cmd_line_parser(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    return inps

def get_minvalue(a,b):
    """return the relative minimum value of a and b"""
    if a < b:
        return a
    else:
        return b

def generate_command(inps):
    """generate batch command"""
    minvisco = inps.minvis[0]
    maxvisco = inps.maxvis[0]
    stepvisco = inps.stepvis[0]

    run_batch = []

    for lc_vis in np.arange(minvisco, maxvisco + stepvisco, stepvisco):
        for um_vis in np.arange(minvisco, maxvisco + stepvisco, stepvisco):
            minviscosity = int(get_minvalue(lc_vis,um_vis) * 100)
            vis_combination = str(int(lc_vis * 100)) + str(int(um_vis * 100))
            workdir = inps.workdir[0] + 'iran_vs_' + vis_combination + '/'
            command = 'postprocess_points_relax.py ' + workdir + ' --lalo ' + inps.lalo[0] + ' --minvisco ' + str(minviscosity) + ' --start ' + str(inps.stime[0]) + ' --end ' + str(inps.etime[0]) + ' --Tstep ' + str(inps.tstep[0]) + ' --outname ' + inps.outname[0] + vis_combination + ' --outdir ' + inps.outdir[0]
            
            run_batch.append(command)            
    
    with open('postprocess_points_relax.sh', 'w') as h:
        for line in run_batch:
            lines = line + '\n'
            h.write('%s' % lines)
    h.close()

######################################################################
def main(iargs=None):
    inps = cmd_line_parser(iargs)
    generate_command(inps)
#################################################################################################
if __name__ == '__main__':
    main()
