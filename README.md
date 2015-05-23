AlignFreeTools
===========

basic tools for computing alignment-free sequence comparison measures

Requirements
---------------

AlignFreeTools is intended to be used in a Unix-based environment. It has been tested
on Mac OS and Linux.


Usage
---------------

This package provides the c++ programs for computing multiple alignment-free sequence comparison meausres, including d2, d2star and d2shepp, hao, Eu, EuF, Ch, Ma, JS, Willner.

The first step for computing alingment-free measures is to count k-tuple frequencies from the fasta/fastq files. The corresponding program is "countKmer.cpp".
To use the program,

> cd /the/path/to/countKmer.cpp/

> g++ countKmer.cpp -o countKmer.out

> ./countKmer.out -k <k-tuple length> -i <input fasta/fastq file> -s <short name> -o <output directory> -q <?fastq> -z <output words with zero count?> -d <count for double strands?> -l <input file is ref seq?> -p <input file has no annotation?>


The second step is to compute multiple alignment-free sequence comparison meausres. The corresponding program is "computeD2MC_multiStat.cpp".
To use the program,

> cd /the/path/to/computeD2MC_multiStat.cpp/

> g++ computeD2MC_multiStat.cpp -o computeD2MC_multiStat.out

> ./computeD2MC_multiStat.out -a <name of species1> -b <MC order of species1> -c <name of species2> -d <MC order of species2> -k <k-tuple length> -i <directory to kmer count output files of species1> -j <directory to kmer count output files of species2> -o <output directory>


For some measures such as d2star, d2shepp and EuF, they model the genome sequence using Markov chain (MC) model. To fit the data with MC, we first need to estimte the order of MC. For a given long genome sequence, we estimate the order of MC using AIC and BIC criteria. The corresponding program is "computeAICBIC.cpp". 
To use the program,

> cd /the/path/to/computeAICBIC.cpp/

> g++ computeAICBIC.cpp -o computeAICBIC.out

> ./computeAICBIC.out -a <name of species> -b <k-tuple length> -i <directory to kmer count output files of the species> -o <output directory>



Contacts and bug reports
------------------------
Jie Ren
renj@usc.edu

Fengzhu Sun
fsun@usc.edu

If you found a bug or mistake in this project, we would like to know about it.
Before you send us the bug report though, please check the following:

1. Are you using the latest version? The bug you found may already have been
fixed.
2. Check that your input is in the correct format and you have selected the
correct options.
3. Please reduce your input to the smallest possible size that still produces
the bug; we will need your input data to reproduce the problem, and the
smaller you can make it, the easier it will be.


Copyright and License Information
---------------------------------
Copyright (C) 2014 University of Southern California, Jie Ren

Authors: Jie Ren

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see http://www.gnu.org/licenses/.