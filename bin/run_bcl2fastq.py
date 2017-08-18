#!/usr/bin/env python 

import os
import sys
import re
import errno
import glob
import time
import pickle
import logging
from joblib import Parallel, delayed
import multiprocessing
import subprocess
#sys.path.insert(0,'/usr/local/bin/')
import functions
import argparse

__version__ = 'v01'
# created on 17/08/2016

if __name__ == '__main__':

    """ This script creates a file which needs to be filled with information required for a methylSeq project.
    - It takes one argument, the 'outfile', which is the name of the output file. The default is 'analysis_info.txt'""" 

    parser=argparse.ArgumentParser(prog='analysis_info.py', description='Creates analysis_info.txt')
    parser.add_argument('-v','--version',action='version',version='%(prog)s-'+__version__)
    parser.add_argument('--analysis_info_file', help='Text file with details of the analysis. Default=analysis_info.txt', default='analysis_info.txt')
    args=parser.parse_args()

    # Collect info from analysis_info_file
    ai=functions.read_analysis_info_file(args.analysis_info_file)
    
    functions.make_sure_path_exists(ai['project_location'] + '/' +ai['bcl2fastq_output']) 
    os.system("bcl2fastq -R " + ai['run_folder'] + " -o " +  ai['project_location'] + '/' + ai['bcl2fastq_output'] + " --no-lane-splitting --sample-sheet " + ai['run_samplesheet'] + '&>' + ai['project_location'] + '/bcl_log.txt')
    


