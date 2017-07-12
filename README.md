RnaSeq pipeline
================

-   [Getting Started](#getting-started)
-   [Running the pipeline](#running-the-pipeline)
    -   [Analysis info file](#analysis-info-file)

Getting Started
---------------

Create a new folder in fs3 or any other location with the project name - called the "project folder" in this README.
Create a bin/ folder in the project folder.
Download and copy the scripts to the bin/ folder.

Running the pipeline
--------------------

### Analysis info file

A central part of the pipeline is the **analysis info** file. It has information about the project location, the original run folder, the reference fasta, gtf and bed files, and the parameters used throughout the analysis.

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

|                                                                                                                                                |
|:-----------------------------------------------------------------------------------------------------------------------------------------------|
| **working\_directory =** *&lt;path to directory of the analysis&gt;*                                                                           |
| **run\_folder =** *&lt;path to the run folder&gt;*                                                                                             |
| **run\_samplesheet =** *&lt;sample sheet to be used to generate fastq files. This is created using the Illumina Expert Manager&gt;*            |
| **bcl2fastq\_output =** *&lt;path to the desired output of bcl2fastq. The defaults is fastq/ and the folder will be created automatically&gt;* |
| **readType =** *&lt;either pairedEnd or singleEnd&gt;*                                                                                         |
| **reference\_genome =** *&lt;path to the STAR reference genome folder that will be used at the mapping step&gt;*                               |
| **trimgalore\_params =** *&lt;parameters to be passed to trim galore&gt;*                                                                      |
| **ncores =** *&lt;Number of cores to use to pararellize analysis&gt;*                                                                          |
