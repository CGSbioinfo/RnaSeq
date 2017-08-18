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

    parser=argparse.ArgumentParser(prog='analysis_info.py', description='Creates analysis_info.txt')
    parser.add_argument('-v','--version',action='version',version='%(prog)s-'+__version__)
    parser.add_argument('--outfile', help='Name of the output file. Default=analysis_info.txt', default='analysis_info.txt')
    args=parser.parse_args()

    outfile_name=args.outfile
    lines=['project_location = ', 'run_folder = ', 'run_samplesheet = ', 'bcl2fastq_output = fastq/', 'reference_genome = ', 'trimgalore_params = --gzip; --paired; --fastqc; --fastqc_args \'--nogroup --extract\'', 'mapping_params = --runThreadN 4; --outSAMtype BAM SortedByCoordinate; --readFilesCommand zcat ','ncores = 8', 'readType = pairedEnd']
    outfile=open(outfile_name,'w')
    for l in lines:
        outfile.write(l + '\n')
    outfile.close()
    
