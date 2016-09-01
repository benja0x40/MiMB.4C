Chromosome Conformation Capture on Chip (4C): Data Processing
================================================================================

## A. Introduction and prerequisite ##

This archive contains scripts and initial data necessary to perform a demo
analysis of 4C data. See the [references](#1) below for
further informations.

The following dependencies must be installed prior to running the demo:

  - bowtie short read aligner version 1 (http://bowtie-bio.sourceforge.net)
  - R environment version 3.2 or higher
  - R packages: devtools, stringr, getopt, plotrix
  - Bioconductor packages: Biostrings, GenomicRanges

## B. Quick start ##

Running a demo on the 4C study from [Bantignies et al., 2011](#2)

1. Install the package [MRA.TA](https://github.com/benja0x40/MRA.TA) in the R
environment:

```R
library("devtools")
install_github("benja0x40/MRA.TA")
```

2. Run the demo. In the terminal (bash), with current working directory at the
root of MiMB.4C:

```bash
./dataPreparation.sh
Rscript enrichmentAnalysis.R
```

## C. Content of MiMB.4C ##

### 1. Bash scripts ###

  * `importGenome.sh`
  
    Generic tool for automated download of UCSC genomes
  
  * `importRawData.sh`
  
    Tool to facilitate downloading of the demo data from GEO
  
  * `dataPreparation.sh`
  
    Workflow bash script to be executed first when running the demo analysis.
    This script chains several operations:
    importing genome sequence and raw data, updating micro array design probes
    to the lastest (dm6) genome assembly, filtering out non-experimental probes,
    mapping bait sequence (Fab7) to the genome, and computing the restriction
    map associated to the 4C protocol (DpnII).
  
### 2. R scripts ###

  * `updateDesignData.R`
  
    Tool to update micro array design information for new genome releases. This
    tool uses the bowtie alignment software to map probe sequences on a chosen
    assembly of the genome for a precise update of genomic coordinates
    addressed by the microarray platform.
    
  * `computeRestrictionMap.R`
  
    Tool to compute the genomic coordinates of all restriction sites for a given
    restriction enzyme. Restriction sites are defined by the short DNA motif
    (4 to 6bp in common protocols) specifically targeted by the enzyme.
  
  * `enrichmentAnalysis.R`
  
    Workflow R script to be executed secondly (after dataPreparation.sh) when
    running the demo analysis.
    This script chains raw normalization, probes filtering and multi-resolution
    analysis of the 4C enrichments (see the reference publication in section E
    below for more informations).
  
### 3. Default/demo file organisation (before execution) ###

    4C_Bait_Sequences => Fasta file of 4C bait sequence and its bowtie alignments
    Experiment_Design => Tables defining necessary 4C experiment data annotations

### 4. Default/demo file organisation (after execution) ###

    Genome_Data    => UCSC genome data and corresponding bowtie indexes
    UCSC_dm6       => dm6 assembly of the Drosophila melanogaster genome
    Raw_Data       => Raw 4C data
    GSE23887_RAW   => 4C data from Bantignies et al., 2011
    Processed_Data => Pre-processed 4C data
    RestrictionMap => Coordinates of restriction fragments
    Updated        => Updated microarray design data (new probes coordinates)
    Cleaned        => Updated array design filtered out for non-relevant probes

## D. Getting information about bash and R script parameters ##

Runinng `importGenome.sh` or `importRawData.sh` bash scripts without any
parameters will show information about available parameters.
The standalone R scripts `updateDesignData.R` and `computeRestrictionMap.R` can
also provide help on available parameters by using the option -h.

For example (in the terminal):

```bash
./importGenome.sh
Rscript updateDesignData.R -h
Rscript computeRestrictionMap.R -h
```

## E. References ##

<a name="1"></a>1. Leblanc B., Comet I., Bantignies F., and Cavalli G., *Chromosome Conformation Capture on Chip (4C): data processing.* Book chapter to appear in *Polycomb Group Proteins.* Lanzuolo C., Bodega B. editors, Methods in Molecular Biology (2016).  
links: [publisher](https://www.springer.com/gp/book/9781493963782)

<a name="2"></a>2. Bantignies F., Roure V., Comet I., Leblanc B., Schuettengruber B., Bonnet J., Tixier V., Mas A., Cavalli G. *Polycomb-dependent regulatory contacts between distant Hox loci in Drosophila.* Cell 144, 214–26 (2011).  
links: [publisher](http://dx.doi.org/10.1016/j.cell.2010.12.026) | [pubmed](https://www.ncbi.nlm.nih.gov/pubmed/21241892) | [pdf](https://www.researchgate.net/publication/49762071_Polycomb-Dependent_Regulatory_Contacts_between_Distant_Hox_Loci_in_Drosophila)
