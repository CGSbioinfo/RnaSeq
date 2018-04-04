# Script to create the blank analysis info file for the Microarray pipeline - to be filled in for each project.

import argparse

__version__ = "v01"
# Created on 04/04/2018


if __name__ == "__main__": # If this program is being run independently, rather than part of a module, run the script.
    """ This script will create an analysis_info blank template text file which requires the user to fill in the 
    specific paths and parameters necessary for running the Microarray pipeline. It takes one argument,
     the 'output_file', which is the name of the output text file. By default this is 'analysis_info.txt'"""

    parser=argparse.ArgumentParser(prog="create_analysisinfo_file.py", description="Creates analysis_info.txt")
    parser.add_argument("-v", "--version", action="version", version="%(prog)s-"+__version__)
    parser.add_argument("--output_file", help="Name of the output file to be written to. Default is 'analysis_info.txt'", default="analysis_info.txt")
    args=parser.parse_args()

    output_file_name = args.output_file  # Take the user given output file name or use the default.
    lines = ["project_location =", "rawCELs_folder =", "QCC_folder ="]

    output_file = open(output_file_name, "w")

    for i in lines:
        output_file.write(i +"\n")

    output_file.close()