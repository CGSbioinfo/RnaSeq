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
import functions
import argparse

__version__ = 'v01'
# created on 17/08/2016

if __name__ == '__main__':

    """ This script creates a file which needs to be filled with information required for a methylSeq project.
    - It takes one argument, the 'outfile', which is the name of the output file. The default is 'analysis_info.txt'""" 

    parser=argparse.ArgumentParser(prog='create_sampleNames.py', description='Creates sample_names.txt with sample names')
    parser.add_argument('-v','--version',action='version',version='%(prog)s-'+__version__)
    parser.add_argument('--analysis_info_file', help='Text file with details of the analysis. Default=analysis_info.txt', default='analysis_info.txt')
    parser.add_argument('--in_dir', help='Directory saving the fastq files. Default: The one defined in analysis_info in bcl2fastq_output', default='bcl2fastq_output')
    parser.add_argument('--out_file', help='Name of output file. Default=./sample_names.txt', default='sample_names.txt')
    args=parser.parse_args()

    
    # Collect info from analysis_info_file
    if args.in_dir == 'bcl2fastq_output': 
        ai= functions.read_analysis_info_file(args.analysis_info_file)
        allFiles=functions.get_filepaths(ai['project_location'] + '/' + ai['bcl2fastq_output'])
        fastq=[allFiles[y] for y, x in enumerate(allFiles) if re.findall("fastq.gz", x)]
        fastq=[fastq[y] for y,x in enumerate(fastq) if not re.findall('Undetermined', x)]
        #print fastq
    elif args.in_dir != 'bcl2fastq_output':
        ai= functions.read_analysis_info_file(args.analysis_info_file)
        # allFiles=os.listdir(args.in_dir)
        allFiles= functions.get_filepaths(args.in_dir)
        fastq=[allFiles[y] for y, x in enumerate(allFiles) if re.findall("fastq.gz", x)]
        
    fastq=[re.sub('.*/','', f) for f in fastq]
    fastq=[re.sub('_R[0-9_].*','',f) for f in fastq]
    fastq=list(set(fastq))
    fastq=[f.replace('.fastq.gz', '') for f in fastq]
    #print fastq
    fastq=[f.replace('.','') for f in fastq]
    #print fastq 
    

    outfile=open(ai['project_location']+ '/' + args.out_file, 'w')
    for f in fastq:        
        outfile.write(f.split('/')[-1] + '\n')
    outfile.close()        
    
    


