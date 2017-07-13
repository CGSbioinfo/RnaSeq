RnaSeq pipeline
================

-   [Getting Started](#getting-started)
-   [Running the pipeline](#running-the-pipeline)
    -   [1. Analysis info file](#analysis-info-file)
        -   [What is it?](#what-is-it)
        -   [How to create one?](#how-to-create-one)
    -   [2. Convert BCL to FASTQ files](#convert-bcl-to-fastq-files)
        -   [Lexogen projects](#lexogen-projects)
    -   [3. Create a sample names file](#create-a-sample-names-file)
    -   [4. Organize working directory](#organize-working-directory)
    -   [5. Quality control of Fastq files](#quality-control-of-fastq-files)

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

### 2. Convert BCL to FASTQ files

Create fastq files (If you have more than one run, you can run this command for each run by changing the run\_folder and run\_samplesheet fields in the analysis\_info.txt file)

``` bash
$ python bin/run_bcl2fastq.py --analysis_info_file analysis_info.py
```

There is also a log file created called "bcl\_log.txt" to check that the conversion didn't have any problems.

#### Lexogen projects

For lexogen projects, we currently download the data from Basespace. This involves extra steps:

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

2.  Concatenate the 4 files into one, so that at the end there is one file per sample (or two if its paired end). **The name of the concatenated file is important**. The only difference to the original file names, is that concatenated files do not contain the string *L00\[1,2,3,4\]*.
    -   sample1\_S1\_R1\_001.fastq
    -   sample1\_S1\_R2\_001.fastq

3.  It is possible that the **folder names** of the samples will not be in the right format. The scripts expect that each sample has a separate folder and that the folder names consist of the sample name and the \_S\[\\d\] string. From the example above, the folder name would be:
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
