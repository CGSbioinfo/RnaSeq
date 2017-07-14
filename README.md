RnaSeq pipeline
================

-   [Getting Started](#getting-started)
-   [Running the pipeline](#running-the-pipeline)
    -   [1. Analysis info file](#analysis-info-file)
        -   [What is it?](#what-is-it)
        -   [How to create one?](#how-to-create-one)
    -   [2. Obtain FASTQ files](#obtain-fastq-files)
        -   [Using bcl2fastq](#using-bcl2fastq)
        -   [Downloading data from basespace (Currently used for Lexogen projects)](#downloading-data-from-basespace-currently-used-for-lexogen-projects)
    -   [3. Create a sample names file](#create-a-sample-names-file)
    -   [4. Organize working directory](#organize-working-directory)
    -   [5. Quality control of Fastq files](#quality-control-of-fastq-files)
    -   [6. Table and plot of number of reads per sample](#table-and-plot-of-number-of-reads-per-sample)

Getting Started
---------------

1.  Create a new folder in fs3 or any other location with the project name - called the "project folder" in this README.
2.  Create a bin/ folder in the project folder.
3.  Download and copy the scripts to the bin/ folder.

Running the pipeline
--------------------

### 1. Analysis info file

A central part of the pipeline is the **analysis info** file. It has information about the project location, the original run folder, the reference fasta, gtf and bed files, and the parameters used throughout the analysis.

#### What is it?

The analysis info file is a simple .txt file, with each line providing information. Parameters are separated by the semicolon (i.e ";") character.

The following is an example of the analysis info file:

|                                                                                                         |
|:--------------------------------------------------------------------------------------------------------|
| **working\_directory =** /mnt/cgs-fs3/Sequencing/Pipelines/RNASeq\_pipeline/example\_small\_files/      |
| **run\_folder =** /mnt/cgs-fs3/Sequencing/Pipelines/RNASeq\_pipeline/example\_small\_files/xxx          |
| **run\_samplesheet =** /mnt/cgs-fs3/Sequencing/Pipelines/RNASeq\_pipeline/example\_small\_files/xxx     |
| **bcl2fastq\_output =** /mnt/cgs-fs3/Sequencing/Pipelines/RNASeq\_pipeline/example\_small\_files/fastq/ |
| **readType =** pairedEnd                                                                                |
| **reference\_genome =**                                                                                 |
| **trimgalore\_params =** --gzip; --paired; --fastqc; --fastqc\_args '--nogroup --extract'               |
| **ncores =** 8                                                                                          |

The following is the explanation of the analysis info file:

|                                                                                                                                     |
|:------------------------------------------------------------------------------------------------------------------------------------|
| **working\_directory =** *&lt;path to directory of the analysis&gt;*                                                                |
| **run\_folder =** *&lt;path to the run folder&gt;*                                                                                  |
| **run\_samplesheet =** *&lt;sample sheet to be used to generate fastq files. This is created using the Illumina Expert Manager&gt;* |
| **bcl2fastq\_output =** *&lt;path to output of bcl2fastq. The defaults is fastq/ and the folder will be created automatically&gt;*  |
| **readType =** *&lt;either pairedEnd or singleEnd&gt;*                                                                              |
| **reference\_genome =** *&lt;path to the STAR reference genome folder that will be used at the mapping step&gt;*                    |
| **trimgalore\_params =** *&lt;parameters to be passed to trim galore&gt;*                                                           |
| **ncores =** *&lt;Number of cores to use to pararellize analysis&gt;*                                                               |

<br>

#### How to create one?

Use the bin/analysis\_info.txt.

``` bash
$ python bin/analysis_info.py
```

This will create a file analysis\_info.txt, which you can open in a text editor and fill.

<br>

### 2. Obtain FASTQ files

#### Using bcl2fastq

Create fastq files (If you have more than one run, you can run this command for each run by changing the run\_folder and run\_samplesheet fields in the analysis\_info.txt file)

``` bash
$ python bin/run_bcl2fastq.py --analysis_info_file analysis_info.py
```

There is also a log file created called "bcl\_log.txt" to check that the conversion didn't have any problems.

#### Downloading data from basespace (Currently used for Lexogen projects)

For lexogen projects, we currently download the data from Basespace. Normally, the downloaded data will have the following structure:

-   A main folder with the project name
-   Subfolders, which correspond to each of the sequenced samples
-   Within each subfolder, four fastq.gz files from the corresponding sample

The tasks done are:

-   Unzip the fastq files:
    -   Move to the downloaded folder.
    -   Use the following commands:

``` bash
$ ls -d * | parallel -j 4 --no-notice "cd {} ; gzip -d *"      
```

-   Change the folder name of each sample to a string consisting of **the sample name, an underscore, a capital S, and one or more numbers (eg. sample1\_S1)**. You can use the following commands, however **make sure that the sample names consists of the original sample name, an underscore, a capital S, and one or more numbers (eg. sample1\_S1)**:
    -   Check what the outcome would be. This command will print the original folder name, followed by a space and the new folder name:
``` bash
$ for i in $(ls -d *); do sample_name=$(ls ${i} | sed 's/.*\///g' | sed 's/_L00[[:digit:]]\+.*//g'  | sort | uniq); echo ${i} ${sample_name}; done   
```   
    -   If it looks correct, use mv to change the folder names:    
``` bash            
$ for i in $(ls -d *); do sample_name=$(ls ${i} | sed 's/.*\///g' | sed 's/_L00[[:digit:]]\+.*//g'  | sort | uniq); mv ${i} ${sample_name}; done     
```       

1.  For each sample, concatenate the four files into one:

    -   The following command tests whether the outcome of the command is right. It prints the sample name, it's four fasta files, and the name of the output concatenated file.

``` bash
$ ls -d * | parallel -j 4 --no-notice "echo {}; ls {}/{}*L001_R1*  {}/{}*L002_R1* {}/{}*L003_R1*  {}/{}*L004_R1*; echo {}/{}_R1_001.fastq"
```

    + If that looks correct, concatenate the files for R1 reads:   

``` bash
$ ls -d * | parallel -j 4 --no-notice "cat {}/{}*L001_R1*  {}/{}*L002_R1* {}/{}*L003_R1*  {}/{}*L004_R1* > {}/{}_R1_001.fastq"
```

    + And for R2 reads for paired end reads:    

``` bash
$ ls -d * | parallel -j 4 --no-notice "cat {}/{}*L001_R2*  {}/{}*L002_R2* {}/{}*L003_R2*  {}/{}*L004_R2* > {}/{}_R2_001.fastq"
```

1.  Download fastq files from Basespace. Basespace splits the data in 4 files. An example of the file names is:
    -   sample1\_S1\_L001\_R1\_001.fastq
    -   sample1\_S1\_L002\_R1\_001.fastq
    -   sample1\_S1\_L003\_R1\_001.fastq
    -   sample1\_S1\_L004\_R1\_001.fastq
        If the run is paired end, there will also be files for R2:
    -   sample1\_S1\_L001\_R2\_001.fastq
    -   sample1\_S1\_L002\_R2\_001.fastq
    -   sample1\_S1\_L003\_R2\_001.fastq
    -   sample1\_S1\_L004\_R2\_001.fastq

2.  3.  Concatenate the 4 files into one, so that at the end there is one file per sample (or two if its paired end). **The name of the concatenated file is important**. The only difference to the original file names, is that concatenated files do not contain the string *L00\[1,2,3,4\]*.
    -   sample1\_S1\_R1\_001.fastq
    -   sample1\_S1\_R2\_001.fastq

4.  It is possible that the **folder names** of the samples will not be in the right format. The scripts expect that each sample has a separate folder and that the folder names consist of the sample name and the \_S\[\\d\] string. From the example above, the folder name would be:
    -   sample1\_S1
        So the relative path be:
    -   sample1\_S1/sample1\_S1\_R1\_001.fastq
    -   sample1\_S1/sample1\_S1\_R2\_001.fastq

<br>

### 3. Create a sample names file

This should be a text file, with one sample name per line. The sample names consist of the sample name and the \_S\[\\d\] string. From the example above, the sample name would be:
- sample1\_S1

<br>

### 4. Organize working directory

### 5. Quality control of Fastq files

To run fastqc in all the samples, use the script bin/qcReads.py.

Script: bin/qcReads.py

Arguments:

|                                                                                                                 |
|:----------------------------------------------------------------------------------------------------------------|
| -h, --help show this help message and exit                                                                      |
| -v, --version show program's version number and exit                                                            |
| --analysis\_info\_file ANALYSIS\_INFO\_FILE. Text file with details of the analysis. Default=analysis\_info.txt |
| --in\_dir IN\_DIR. Path to folder containing fastq files. Default=rawReads/                                     |
| --out\_dir OUT\_DIR. Path to out put folder. Default=rawReads/                                                  |
| --out\_dir\_report OUT\_DIR\_REPORT. Path to out put folder. Default=Report/figure/data/                        |
| --sample\_names\_file SAMPLE\_NAMES\_FILE. Text file with sample names to run. Default=sample\_names.txt        |

Output:
Fastqc files will be created for each sample in the rawReads/ folder

Command example:

``` bash
$ python bin/qcReads.py --analysis_info_file analysis_info.txt --in_dir rawReads/ --out_dir rawReads/ --out_dir_report Report/figure/data/ --sample_names_file sample_names.txt
```

### 6. Table and plot of number of reads per sample

Script: bin/indexQC.R

Requires: ggplot2 and reshape R libraries

Input: The input directory should have the *fastqc\_data.txt* files of all samples of the project. The 7th line in the fastqc\_data.txt files have the total number of sequences. This is the information that the script collects.

Command example:

``` bash
/usr/bin/Rscript bin/indexQX.R rawReads/ Report/figure/data/
```

Output: The output directory will have a csv file with the Total number of Sequences and a bar plot with the proportion of reads represented by each sample from the total number of reads from the run.
