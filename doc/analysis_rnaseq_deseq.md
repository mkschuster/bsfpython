# RNA-seq DESeq2 Analysis

The RNA-seq analysis requires a distinct set of configuration files. A __project-specific configuration file__ in INI
format defines the project name, the reference genome assembly, the reference transcriptome annotation and additional
analysis parameters specific to the entire project. A __sample annotation sheet__ in comma-separated value (CSV) format
provides sample names, file locations, as well as metadata annotation that is required for the differential expression
modelling. A __design annotation sheet__ in CSV format can specify one or more generalised linear models (GLMs) for the
differential expression modelling. Aside the full model formula, reduced model formulas for likelihood ratio testing
(LRT) need specifying. Optionally, the order of factor levels can also be specified so that base-levels such as controls
can serve as the intercept and thus as a sensible reference. The design annotation sheet also allows encoding of
variable annotation for principal component analysis (PCA), multi-dimensional scaling (MDS) and heatmap plots. Finally,
a __contrast annotation sheet__ in CSV format specifies the contrasts that are to be extracted from the GLM in a
design-specific manner. The contrast annotation sheet can usually only be generated after the main analysis has run,
since DESeq2 result names, especially once they involve interaction terms, can be quite unpredictable.

## User-specific Configuration File

A __user-specific configuration file__ in INI format specifies analysis configuration options that are specific to both,
the user and the laboratory. While the configuration itself is quite comprehensive, an
annotated [template](template_bsfpython.ini) is available and should be adjusted and stored under
`$HOME/.bsfpython.ini`.

## Project-specific Configuration File

The __project-specific configuration file__ in INI format can hold one or more analysis-specific configuration sections.
Thereby, the concatenation of Python module and class name serves as the configuration section header. A particular
analysis will only read the configuration section matching its own class name. Further, analysis-specific configuration
options that remain constant between projects, but are specific to the user or the laboratory, can be put under the same
section names into the user-specific configuration file. For an RNA-seq DESeq analysis, two configuration sections are
mandatory.

The `bsf.analyses.star.Star` section requires a minimum of a project name, and the reference transcriptome version that
should be used for the alignment. The reference transcriptome refers to a directory under the default genomes directory
path that is configured in the __user-specific configuration file__ and also implies a reference genome assembly. An
annotated [template](template_star_config.ini) is available.

The `bsf.analyses.rnaseq.DESeq` section also requires a minimum of project name and reference transcriptome version.
Further, analysis-specific configuration options that remain constant between projects, but are specific to the user,
or the laboratory can be put into a `bsf.analyses.star.Star` section of the user-specific configuration file. An
annotated [template](template_rnaseq_deseq_config.ini) is available.

## Sample Annotation Sheet

A __sample annotation sheet__ in comma-separated value (CSV) format provides sample names, NGS files and metadata
annotation that is important for the differential expression modelling. The table requires some fixed variables and
permits some additional project-specific ones. The file format is similar to the sample annotation sheet written by
the `bsf.analyses.picard.IlluminaDemultiplexSam` analysis into the sample archive folder after demultiplexing.
Therefore, most of the information can be copied from one or more sample annotation sheets from the demultiplexing
stage.

The sample annotation sheet file name should obey the schema `<project_name>_<genome_version>_rnaseq_samples.csv` to
allow for automatic configuration file discovery.

### Variables

- `ProcessedRunFolder Name`
- `Project Name` The library name this sample was a part of.
- `Project Size` The library mean fragment size determined by (automated) gel electrophoresis.
- `Sample Group` A sample group this sample should be part of. Please note that this variable is used by
  the `bsf.analyses.rnaseq.Tuxedo`, but not the `bsf.analyses.rnaseq.DESeq` analysis.
- `Sample Name` A (meaningful) sample name.
- `PairedReads Exclude` for excluding particular `PairedReads` objects, equivalent to read groups in BAM files from a
  high-level analysis.
    - `FALSE` to include a particular PairedReads object
    - `TRUE`  to exclude a particular PairedReads object
- `PairedReads Flow Cell` The flow cell identifier including the experiment name, and the flow cell barcode
  (e.g., BSF_0000_AAAAACXX).
- `PairedReads Flow Cell Lane` The `PairedReads Flow Cell` identifier from above and the lane number
  (e.g., BSF_0000_AAAAACXX_2).
- `PairedReads Index 1` The index sequence of the i7 adapter.
- `PairedReads Index 2` The index sequence of the i5 adapter.
- `PairedReads Lane` The integer lane number.
- `PairedReads ReadGroup`
- `PairedReads Structure` The read structure in Picard format indicating the number of cycles for the first read, the
  first index read, the second index read, and the second read, if available (e.g., `51T8I` or `75T8I8I75T`).
- `Reads1 Name` The logical name of the first read which will be used for file naming.
- `Reads1 File` The path to the file representing the first read in the case of FASTQ files, or the unaligned BAM file.
  Absolute paths are used verbatim, relative paths are prepended with the BSF sample archive path configured in the
  user-specific configuration file.
- `Reads2 Name` Equivalent to `Reads1 Name`.
- `Reads2 File` Equivalent to `Reads1 File`.
- `Sample DESeq designs` A comma-separated string of design names this sample should be part of.
- `Sample DESeq library_type` Indicating the strand orientation of the experiment, if any, and how reads should be
  assigned during counting.
    - `first` NGS reads directly represent the biological mRNA molecule
    - `second` NGS reads represent the opposite strand of the biological mRNA molecule, such as dUTP-based protocols.
    - `unstranded` The protocol does not preserve strand information, such as the _Smart-seq2_ protocol.
- `Sample DESeq sequencing_type` The
  Bioconductor [GenomicAlignments](https://bioconductor.org/packages/release/bioc/html/GenomicAlignments.html)
  summarizeOverlaps function needs information about read pairing for counting.
    - `SE` for single-end data
    - `PE` for paired-end data
- `Sample DESeq RIN` If present, allows for plotting the distribution of the RNA integrity number (RIN) over all
  samples.
- `Sample DESeq run` If present, indicates that technical replicates should be collapsed according to information in
  the `Sample Name` variable. The `Sample DESeq run` variable provides the original sample name before collapsing
  technical replicates. Normally, the `bsf.analyses.star.Star` analysis aligns and post-processes each read
  group-specific BAM file separately, before merging into a sample-specific aligned BAM file. Hence, annotation of this
  variable is only required under exceptional circumstances.
- `Sample DESeq ...` Further variables meaningful for the differential expression modelling can be annotated with
  a `Sample DESeq ...` prefix and will be propagated to the Bioconductor DESeq2 analysis. Examples could be:
    - `Sample DESeq genotye`
    - `Sample DESeq phenotype`
    - `Sample DESeq treatment`
- `Sample DESeq total_counts` A reserved variable, since the total number of reads per sample is automatically
  calculated by the DESeq2 analysis script.

Please note that the `PairedReads Flow Cell`, `PairedReads Lane`, `PairedReads Flow Cell Lane`
and `PairedReads Structure` variables, as well as the automatically calculated `Sample DESeq total_counts` variable can
serve as aesthetics at the exploratory analysis stage to pin down technical variability or serve as covariates in the
differential expression modelling.

## Design Annotation Sheet

A __design annotation sheet__ in comma-separated value (CSV) format can specify one or more designs, each specifying a
_generalised linear model_ (GLM) formula for the differential expression modelling. Aside the full model formula,
further reduced model formulas for _likelihood ratio testing_ (LRT) need specifying. Optionally, the order of factor
levels can also be specified so that base-levels such as controls can serve as the intercept and thus as a sensible
reference. The design annotation sheet also allows encoding of variable annotation for the exploratory stage, such as
_principal component analysis_ (PCA), _multi-dimensional scaling_ (MDS) and _heatmap_ plots.

### Variables

- `design` A (meaningful) design name that has to correlate with the comma-separated values in
  the `Sample DESeq designs` variable of the __sample annotation sheet__.
- `exclude` Exclude particular designs from reporting.
    - `FALSE` to include this design in the report
    - `TRUE` to exclude this design from the report
- `full_formula` The full model formula for the Bioconductor DESeq2 analysis (e.g., `~ genotype * phenotype`).
- `reduced_formulas` A string representation of one or more named, reduced formulas for likelihood ratio testing (LRT),
  separated by semicolons
  (e.g., `wo_genotype:~ phenotype;wo_phenotype:~ genotype;wo_interaction:~ genotype + phenotype`).
- `factor_levels` A string representation to order the levels for one or more factors to properly set the reference for
  comparisons.
- `plot_aes` A string representation of variables mapped to [ggplot2](https://ggplot2.tidyverse.org/) _Geoms_ and
  their _aesthetics_ for one or more plots. The plots are separated by pipe characters (`|`), the _geoms_ and a
  comma-separated list of aesthetics and variable pairs are separated by a colon
  (e.g., `geom_point:colour=total_counts|geom_point:colour=genotype,shape=phenotype`).

## Contrast Annotation Sheet

A __contrast annotation sheet__ in comma-separated value (CSV) format specifies those contrasts (i.e. meaningful
biological comparisons) that should be extracted from the _generalised linear model_ (GLM) in a design-specific manner.
The __contrast annotation sheet__ can usually only be generated after the main analysis has run, since DESeq2 result
names, especially once they involve interaction terms, can be quite unpredictable.

### Variables

- `Design` The design name correlating to the __design annotation sheet__.
- `Numerator` The numerator for the contrast.
- `Denominator` The denominator for the contrast.
- `Label` A human-readable label for this particular contrast.
- `Exclude` Exclude particular contrasts from reporting.
    - `FALSE` to include this contrast in the report
    - `TRUE` to exclude this contrast from the report

## Configuration File Naming Schema

The following configuration files should obey a naming schema to allow for automatic discovery in the configuration file
directory from which the analysis gets submitted. Since project names are used in directory and file names, Perl word
characters (i.e. `[0-9A-Za-z_]+`) would be particularly safe to use on most operating and file systems.

- Project-specific Configuration File
    - `<project_name>_<genome_version>_rnaseq_config.ini`
- Sample Annotation Sheet
    - `<project_name>_<genome_version>_rnaseq_samples.csv`
- Design Annotation Sheet
    - `<project_name>_<genome_version>_rnaseq_designs.csv`
- Contrast Annotation Sheet
    - `<project_name>_<genome_version>_rnaseq_contrasts.csv`

## Process Submission

Once the above configuration files are complete, enter the configuration directory and call the process submission
script. Help on arguments and options is available via `bsf_submit_rnaseq_deseq.py --help`. The most commonly used
option are `--logging-level` to set logging verbosity and `--dry-run` to __not__ submit any processes to the distributed
resource management system.

By default, the script computes the complete process graph and submits it to the distributed resource management
system (DMRS). Processes from individual stages can be submitted by adding the `--stage-name` option followed by the 
analysis stage name. The DESeq analysis supports the `rnaseq_deseq_analsis` and `rnaseq_deseq_results` stages.

`bsf_submit_rnaseq_deseq.py <project_name>_<genome_version>_rnaseq_config.ini`

After process submission, use the relevant DRMS commands to assess the state of analysis completion. Each process –
a `bsf.process.Runnable` object equivalent to a script – will write an empty `*_complete.txt` status file into the
genome directory upon successful completion and a `*_failed.txt` otherwise. Processes, which corresponding status files
indicate successful completion will not be resubmitted upon rerunning of the submission script. While a process is
running, its individual steps – a `bsf.process.RunnableStep` object equivalent to a command in a script – also leave
similar status files, so that processes can rerun just the incomplete sub-steps upon a resubmission. Processes, and
their dependencies are tracked in a SQLite database, for each process a standard output and a standard error file are
written into the `bsfpython_<DRMS>_output` directory.

Once all processes of the analysis have successfully completed, a report can be generated with the submission script.
Use `bsf_submitrnaseq_deaseq.py --stage-name='report'` in the configuration directory to generate the report. 
The report will also be linked symbolically into the web server export directory.

## Troubleshooting

If processes fail, they can be resubmitted by rerunning the above submission script. Status files allow rerunning of the
missing analysis steps for a particular process. If this is not desired, please delete the corresponding files at the
particular stage so that the entire process will be repeated from the beginning. However, sometimes intermediate files
of previous analysis stages are deleted during the analysis path, so that sometimes processes from upstream stages need
also rerunning. Again, the corresponding files including the status files need deleting.
  