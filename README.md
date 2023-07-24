# BSF Python Library

## Introduction

The [Biomedical Sequencing Facility](http://www.biomedical-sequencing.at/) (BSF) is part of the
joint [Genomics Core Facility](http://corefacilities.meduniwien.ac.at/genomics/?L=1) of the
[Medical University of Vienna](http://www.meduniwien.ac.at/) (MUW) and the
[Research Center for Molecular Medicine](http://www.cemm.oeaw.ac.at/) (CeMM) of the
[Austrian Academy of Sciences](http://www.oeaw.ac.at/) (OeAW).
The BSF is Austria’s first technology platform dedicated to
[next-generation sequencing](http://en.wikipedia.org/wiki/DNA_sequencing#Next-generation_methods) (NGS)
in biomedicine and expected to play a catalyzing role for the
development of genomic medicine in [Vienna](http://en.wikipedia.org/wiki/Vienna) and
[Austria](http://en.wikipedia.org/wiki/Austria).

This Python library and the accompanying scripts are used for day-to-day analysis of
next-generation sequencing data sets. The library consists of two main functions.

The `Analysis` class and its subclasses implement the logic required for submitting processes on a
cluster login node.

The `Runnable` class and its subclasses implements the logic required to run processes on a
cluster compute node via the common `bsf_run_runnable` script.

## BSF Python General Configuration File

General settings for the BSF Python library are configured via a `${HOME}/.bsfpython.ini` file
in the user's home directory. This file is site-specific and its information allows for automatic
discovery of raw (e.g., Illumina run folders) and pre-processed (e.g., de-multiplexed lanes and samples)
NGS data. A template file [`template_bsfpython.ini`](doc/template_bsfpython.ini) can be found in the doc subdirectory.
The template, which documents the configuration options and provides, as far as possible, sensible default settings,
needs copying to `${HOME}/.bsfpython.ini` before editing accordingly.

## Analysis

The BSF `Analysis` is central to the BSF pipeline infrastructure. It encapsulates both, logic and data
for a multistep analysis procedure. Specific Analysis objects are available, tailored to recurring
tasks.

## Analysis Configuration File

BSF Analysis objects are initialised and configured via INI configuration files.
Templates for these files are again provided in the doc subdirectory, document configuration options and
provide, as far as possible, sensible default settings. Generally, only few configuration
options need filling in. Most importantly, the location of sample annotation sheets and, depending on the
analysis type, sample comparison sheets, need to be specified.

## Sample Annotation Sheet

A sample annotation sheet specifies the file system location of NGS reads. For data pre-processed via
Illumina CASAVA, a hierarchy of run folders, projects, samples, and paired reads can be automatically
discovered. Additional reads can be linked into the system by specifying the exact file system path.
Sample annotation sheets also provide grouping of samples that is available to the analysis.

    - Type (CASAVA or External)
    - ProcessRunFolder
    - Project
    - Sample
    - Reads1
    - File1
    - Reads2
    - File2
    - Group

## Analyses

### ChIPSeq

The ChIPSeq analysis aligns each BSF Sample object to the genome sequence via BWA. Regions of interest
are then defined by means of the MACS2 peak caller.

In the context of the ChIPSeq analysis, BSF Paired Reads objects of BSF Sample objects are aligned as a pool.

### RNAseq DESeq

Please see the [RNAseq DESeq](doc/analysis_rnaseq_deseq.md) document.

### RNASeq Tuxedo

The RNA-Seq pipeline is based on the Tuxedo suite. NGS reads are aligned with
[Tophat2](http://ccb.jhu.edu/software/tophat/index.shtml) an aligner that
implements a splice site model and uses a reference transcriptome as the base.
it is based on the [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) short read aligner.

In the context of the Tuxedo analysis, BSF Sample objects are aligned and assembled into transcriptomes individually.
According to a group_replicates configuration option, each BSF Paired Reads object of a BSF Sample object can be
processed individually by TopHat and Cufflinks or pooled before alignment. The resulting transcriptome assemblies
for each BSF Sample resulting from each BSF Paired Reads object are then merged via Cuffmerge. Cuffdiff then compares
the merged assemblies on the basis of the BAM alignments produced by Tophat2.

### Genetic Variant Calling

## Dependencies

* Python ([The Python Software Foundation](https://www.python.org/))
* R ([The R Foundation for Statistical Computing Platform](https://www.r-project.org/))

### ChIP-seq

* MACS2 ([Tao Liu lab at UB and Xiaole Shirley Liu lab at DFCI](https://github.com/taoliu/MACS))
* DiffBind
* Kent Tools
* convert ([ImageMagick Studio LLC](http://www.imagemagick.org/))

## References

## Licence

Copyright 2013 - 2020 Michael K. Schuster

[Biomedical Sequencing Facility](http://www.biomedical-sequencing.at/) (BSF),
part of the joint
[Genomics Core Facility](http://corefacilities.meduniwien.ac.at/genomics/?L=1) of the
[Medical University of Vienna](http://www.meduniwien.ac.at/) (MUW) and the
[Research Center for Molecular Medicine](http://www.cemm.oeaw.ac.at/) (CeMM) of the
[Austrian Academy of Sciences](http://www.oeaw.ac.at/) (OeAW).

This file is part of BSF Python.

BSF Python is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BSF Python is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with BSF Python. If not, see <http://www.gnu.org/licenses/>.
