Chromosome Conformation Capture on Chip (4C): Data Processing
================================================================================

The `MiMB.4C` repository contains a collection of bash and R scripts implementing
our data processing protocol for Chromosome Conformation Capture on Chip (4C).
See [Leblanc et al., 2016](#1) for a detailed presentation of this protocol and
underlying methods, which consists in improved versions of the 4C data
analyses procedures we originally used in [Bantignies et al., 2011](#2).

These procedures, including normalization, probre selection and
multi-resolution analysis, are available via the standalone R package
[MRA.TA](https://github.com/benja0x40/MRA.TA), whereas the `MiMB.4C`
repository illustrates a complete workflow for 4C data processing.

Below are two example of results generated using the
[MRA.TA](https://github.com/benja0x40/MRA.TA) package and based
on 4C data in mouse from [Simonis et al., 2006](#3) (top panel) and
[Schoenfelder et al., 2009](#4) (bottom panel). 

![](./images/examples/MiMB.4C_Examples_smallsize.png "")

A complete demo analysis based on our 4C data in *Drosophila* embryos can be
run as indicated in the following sections.

### A. Prerequisites ###

These dependencies need to be installed before using `MiMB.4C`:

  - `bowtie` short read aligner version 1.x (http://bowtie-bio.sourceforge.net)
  - R environment version 3.x
  - R packages: `devtools`, `stringr`, `getopt`, `plotrix`
  - [Bioconductor](http://www.bioconductor.org/) packages: `Biostrings`, `GenomicRanges`

### B. Quick start ###

1. Install the package [MRA.TA](https://github.com/benja0x40/MRA.TA) in the R
environment:

```R
library("devtools")
install_github("benja0x40/MRA.TA")
```

2. Clone or download and decompress the `MiMB.4C` repository


3. Run the demo in a terminal (bash), with current working directory at the
root of the decompressed `MiMB.4C` folder:

```bash
./dataPreparation.sh
Rscript enrichmentAnalysis.R
```

### C. Content of MiMB.4C ###

#### 1. Bash scripts ####

  * `importGenome.sh`
  
    Tool for automated download of UCSC genomes and index creation for bowtie.
  
  * `importRawData.sh`
  
    Tool to facilitate downloading of the [Bantignies et al., 2011](#2) 4C data from GEO.
  
  * `dataPreparation.sh`
  
    Workflow bash script to be executed first when running the demo analysis.
    This script chains several operations:
    importing genome sequence and raw data, updating micro array design probes
    to the lastest (dm6) genome assembly, filtering out non-experimental probes,
    mapping bait sequence (Fab7) to the genome, and computing the restriction
    map associated to the 4C protocol (DpnII).
  
#### 2. R scripts ####

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
  
#### 3. File organisation (before execution) ####

    4C_Bait_Sequences => Fasta file of 4C bait sequence and its bowtie alignments
    Experiment_Design => Tables defining necessary 4C experiment data annotations

#### 4. File organisation (after execution) ####

    Genome_Data    => UCSC genome data and corresponding bowtie indexes
    UCSC_dm6       => dm6 assembly of the Drosophila melanogaster genome
    Raw_Data       => Raw 4C data
    GSE23887_RAW   => 4C data from Bantignies et al., 2011
    Processed_Data => Pre-processed 4C data
    RestrictionMap => Coordinates of restriction fragments
    Updated        => Updated microarray design data (new probes coordinates)
    Cleaned        => Updated array design filtered out for non-relevant probes

### D. Getting information about bash and R script parameters ###

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

### E. References ###

<a name="1"></a>1. Leblanc B., Comet I., Bantignies F., and Cavalli G., *Chromosome Conformation Capture on Chip (4C): data processing.* Book chapter to appear in *Polycomb Group Proteins.* Lanzuolo C., Bodega B. editors, Methods in Molecular Biology (2016).  
[publisher](https://www.springer.com/gp/book/9781493963782)

<a name="2"></a>2. Bantignies F., Roure V., Comet I., Leblanc B., Schuettengruber B., Bonnet J., Tixier V., Mas A., Cavalli G. *Polycomb-dependent regulatory contacts between distant Hox loci in Drosophila.* Cell (2011).  
[publisher](http://dx.doi.org/10.1016/j.cell.2010.12.026) | [pubmed](https://www.ncbi.nlm.nih.gov/pubmed/21241892)

<a name="3"></a>3.Simonis M., Klous P., Splinter E., Moshkin Y., Willemsen R., de Wit E., van Steensel B., de Laat W. *Nuclear organization of active and inactive chromatin domains uncovered by chromosome conformation capture-on-chip (4C).* Nature Genetics (2006).  
[publisher](http://dx.doi.org/10.1038/ng1896) | [pubmed](https://www.ncbi.nlm.nih.gov/pubmed/17033623)

<a name="4"></a>4. Schoenfelder S., Sexton T., Chakalova L., Cope N.F., Horton A., Andrews S., Kurukuti S., Mitchell J.A., Umlauf D., Dimitrova D.S., Eskiw C.H., Luo Y., Wei C.L., Ruan Y., Bieker J.J, Fraser P. *Preferential associations between co-regulated genes reveal a transcriptional interactome in erythroid cells.* Nature Genetics (2009).  
[publisher](http://dx.doi.org/10.1038/ng.496) | [pubmed](https://www.ncbi.nlm.nih.gov/pubmed/20010836)

