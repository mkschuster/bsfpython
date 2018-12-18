# BSF Python Library

## Sample Annotation Sheet

A sample annotation sheet defines a hierarchy of BSF ProcessedRunFolder, BSF Project,
BSF Sample, BSF PairedReads, BSF Reads objects in a comma-separated value (CSV) format.
For Illumina run folders post-processed with CASAVA, all hierarchical objects can be
automatically discovered.

### Columns

The Sample Annotation Sheet format is defined in bsf.ngs.Collection._process_row_dict

 - FileType: BSF Data object file_type (i.e. 'CASAVA', 'External' or 'Automatic'), defaults to 'Automatic'.
             The FileType 'CASAVA' allows for auto-discovery of BSF ProcessedRunFolder and subsequent objects.

 - ProcessedRunFolder: bsf.ngs.ProcessedRunFolder.file_path, can be automatically registered.

 - Project: bsf.ngs.Project.name

 - Sample: bsf.ngs.Sample.name

 - File1: bsf.ngs.Reads.file_name instance variable.
          Subjected to os.path.expanduser and os.path.expandvars.
          If still relative, the BSF Collection file_path is prepended.

 - Reads1: bsf.ngs.Reads.name instance variable

 - File2: Same as File1

 - Reads2: Same as Reads1

 - Groups: BSF Sample objects can be grouped for further analysis in
           e.g. RNA-Seq or ChIP-Seq experiments.

### Class Hierarchy

 - ProcessedRunFolder Name
     - Project Name
         - Sample Name
             - PairedReads ReadGroup
                 - Reads1 Name
                 - Reads1 File
                 - Reads2 Name
                 - Reads2 File

## Sample naming and Replicates

Initially, replicates were implemented as BSF PairedReads objects inside BSF Sample objects.
This strategy was in line with CASAVA processed run folder directory structures, as long as all
BSF Sample objects resulted from the same BSF Processed Run Folder object (i.e flow cell).

Later, when it became clear that samples could result from different flow cells resulting in
BSF Sample objects bearing the same name, the code base needed some adjustments.

The bsf.ngs.Collection.sample_group_dict instance variable has been changed from a Python dict of
BSF Sample name keys and BSF Sample values to a Python list of BSF Sample objects. For the moment,
a check is in place testing, whether the Python list contains already the same BSF Sample object.
The Python 'in' comparison operator uses the memory address in lack of a specific __cmp__ method, but
it may make sense to implement such a more specific method.
See method bsf.ngs.Collection._process_row_dict

The bsf.Analysis.samples instance variable also needed changing from a dict to a list. A similar
check avoiding duplicates is in place in the bsf.Analysis.add_sample method.

BSF Sample objects can now have the same name and be treated as replicates via grouping.


## Sample Grouping

BSF Python sample annotation sheets allow for grouping of BSF Sample objects. The meaning of the
grouping is dependent on the analysis.

Most analyses consist of two steps. An alignment stage, where BSF Sample objects get aligned
individually or as a pool and a subsequent comparison stage.

    RNA-Seq

    Normally, each BSF Sample object gets run through the alignment stage (i.e. TopHat and Cufflinks)
    individually, while replicates (i.e. BSF PairedReads objects) can be processed separately or
    as a pool, depending on the 'replicate_grouping' configuration option.

    At the comparison stage (i.e. Cuffmerge and Cuffdiff) the output from BSF Sample objects in the
    alignment stage gets conceptually pooled.

    ChIP-Seq

    Each BSF Sample object gets run through the alignment stage (i.e. Bowtie2) as a pool so that the
    resulting aligned BAM file contains the super-set of all NGS reads. For the moment, grouping of
    BSF Sample objects is not supported, because it would mean that each sample gets aligned
    separately. Grouping could be implemented via two strategies:
    - BSF Sample objects could be merged before the alignment stage.
    - Aligned BAM file could be merged after the alignment stage.


    For the moment, while reading the comparison file, each sample appearing in a comparison gets
    appended to the Analysis.samples list. For RNA-Seq each sample gets aligned separately, while
    for ChIP-Seq they should be aligned in pools.

    To get around this, it would be good to have a SampleGroups object that gets pushed onto the
    bsf.Analysis.sample list. The SampleGroups object could then provide methods to retrieve
    BSF PairedReads objects as list (i.e. pooled) or separately.


    However, BSF PairedReads objects that result from different flow cells could potentially bear
    the same name, which causes problems with SGE job names that need be unique.
    The bsf.ngs.Collection

    A remaining problem is that BSF PairedReads and BSF Sample objects from different flow cells
    may have the same name. Although the BSF Collection could hold them apart, SGE jobs could still
    end up with the same name.
    Therefore, it may be desirable to make available the complete hierarchy of ProcessRunFolder,
    Project, Sample and PairedReads.


bsf.Analysis

The bsf.Analysis class represents a high-level analysis of NGS data.

    Maybe, some sub-classes could be defined, although this may lead to multiple inheritance
    complexities.

        - bsf.Analysis.Genome
          For analyses depending on a genome assembly, defining an analysis.genome_directory

        - bsf.Analysis.Comparison
          For analyses that involve comparisons of BSF Sample objects.

The bsf.process.Executable class represents an executable program, its options and arguments.

bsf.ngs class hierarchy:

    bsf.ngs.Reads - Represents a single FASTQ file.
    bsf.ngs.PairedReads - Represents paired reads holding one or two Reads objects
    bsf.ngs.Sample - Represents one or more PairedReads objects (replicates)
    bsf.ngs.Project - Represents one or more Sample objects
    bsf.ngs.ProcessedRunFolder - Represents a Processed Run Folder (e.g. after CASAVA run)
    bsf.ngs.Collection - Represents a collection of ProcessedRunFolder objects

The bsf.ngs.Sample configuration still has some problems.

Change Collection._process_row_dict() and eventually Collection.get_sample_from_row_dict()
to deal with wildcards. Also, if a column is empty or does not even exist, register
all hierarchical objects underneath, automatically.

However, for a comparison schema for a BSF Analysis, 'ProcessedRunFolder', 'Project' and 'Sample'
columns are an absolute requirement. ProcessedRunFolder auto-discovery registers all hierarchical
objects (i.e. BSF Project and BSF Sample objects) underneath the BSF ProcessedRunFolder.
Since Sample names are not guaranteed to be unique between BSF Project objects
- they typically are 1, 2, 3, ... - both 'Project' and 'ProcessedRunFolder' need specifying.


bsf.ngs.Pool

Would a new bsf.ngs.Pool object make sense for bsf.ngs.Sample objects that have run
on more than one lane and represent sequencing replicates. For the moment, a solution to
merge two samples is implemented in the bsf.analyses.chipseq.ChIPSeq class.


## RNA-Seq Pipeline:

* Integrate with Doris post-processing R-code
* Add PCA and MDS plots as in the cummeRbund manual.
* Andreas Bergthaler:
  * Add individual FPKM values for each replicate.
    Implemented via bsf_process_cuffdiff.R.
  * Rename downloadable files to include the sample or comparison name.
    Otherwise files with identical names may get overwritten.
    Implemented via bsf_process_cufflinks.R and bsf_process_cuffdiff.R.
  * Maybe create a ZIP or Gzip archive to download and archive locally.


## Object Hierarchy

bsf
  - from bsf import Defaults (no further dependency)
  - from bsf.ngs import Collection (no further dependency)
  - from bsf.argument import * (no further dependency)

bsf.Analysis
  - from bsf import Analysis, Configuration, Default, Defaults, Stage, Executable (no further dependency)
  - from bsf.ngs import Collection, ProcessedRunFolder, Sample, SampleAnnotationSheet (no further dependency)
  - from bsf.executables import Bowtie2, Macs14, Macs2Callpeak, Cuffdiff, Cufflinks, Cuffmerge, TopHat, FastQC

bsf.argument
  - no further dependency

bsf.ngs
  - no further dependency

bsf.defaults
  - no further dependency


Configuration file template for picard.extract_illumina_barcodes

illumina_run_folder =

barcode_file =


# Barcode file:

  - lane
  - barcode_sequence_1
  - barcode_sequence_2
  - library_name
  - sample_name
