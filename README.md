Chromosome Conformation Capture on Chip (4C): Data Processing
================================================================================

This `MiMB.4C` repository contains a collection of bash and R scripts
implementing our data processing protocol for Chromosome Conformation Capture
on Chip (4C). See [Leblanc et al., 2016](#1) for a brief background on the 4C
technique itself and the detailed presentation of this protocol.
For an extensive perspective on Chromosome Conformation Capture technologies,
see for instance the review from [Denker & de Laat, 2016](#2).

Computational methods undelying this protocol consist in improved versions of
the procedures used in the study from [Bantignies et al., 2011](#3).
These methods, including normalization, probe selection and multi-resolution
visualization and segmentation of the 4C profile, are available via the
standalone R package [MRA.TA](https://github.com/benja0x40/MRA.TA),
whereas the scripts provided in the `MiMB.4C` repository illustrate a complete
workflow for 4C data processing.

Below are two example of results generated using these methods based
on 4C data in mouse from [Simonis et al., 2006](#4) (top panel) and
[Schoenfelder et al., 2009](#5) (bottom panel).

![](./images/examples/MiMB.4C_Examples_smallsize.png "")

Both panels represent the mouse chromosome 7 on the horizontal axis and the
resolution of analysis on the vertical axis, in number of microarray probes.
Frequencies of interactions between the 4C bait sequence (which was targeting
the beta globin locus in both studies) and remote sequences along the chromosome
are indicated by colors, from light blue for the weakest levels to dark red for
the strongest ones.  
The whole colormaps represent statistical scores 
reflecting 4C interaction frequencies for each genonic location and
considering resolutions ranging from single probe to approximately 5000 probes.
The 3 tracks below each colormap show alternative segmentations of the
significant interactions. From top to bottom: segmentation at maximal scale or 
at maximal resolution resulting from our protocol, and segmentation reported in
the original studies using former data analysis methods.

A complete demo analysis based on the 4C data in *Drosophila* anterior larval
tissues from [Bantignies et al., 2011](#3) can be run as indicated in the
following sections.


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
  
    Tool to facilitate downloading of the 4C data used for demo from GEO.
  
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
    This script chains normalization, probes filtering and multi-resolution
    analysis of the 4C enrichments (see the first reference in section E
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

<a name="1"></a>1. Leblanc B., Comet I., Bantignies F., and Cavalli G., *Chromosome Conformation Capture on Chip (4C): data processing.* Book chapter in *Polycomb Group Proteins: Methods and Protocols.* Lanzuolo C., Bodega B. editors, Methods in Molecular Biology (2016).  
[publisher](http://dx.doi.org/10.1007/978-1-4939-6380-5_21) | [pubmed](https://www.ncbi.nlm.nih.gov/pubmed/27659990)

<a name="2"></a>2. Denker A. & de Laat W., *The second decade of 3C technologies: detailed insights into nuclear organization.* Genes & Development (2016).  
[publisher](http://dx.doi.org/10.1101/gad.281964.116) | [pubmed](https://www.ncbi.nlm.nih.gov/pubmed/27340173)

<a name="3"></a>3. Bantignies F., Roure V., Comet I., Leblanc B., Schuettengruber B., Bonnet J., Tixier V., Mas A., Cavalli G. *Polycomb-dependent regulatory contacts between distant Hox loci in Drosophila.* Cell (2011).  
[publisher](http://dx.doi.org/10.1016/j.cell.2010.12.026) | [pubmed](https://www.ncbi.nlm.nih.gov/pubmed/21241892)

<a name="4"></a>4. Simonis M., Klous P., Splinter E., Moshkin Y., Willemsen R., de Wit E., van Steensel B., de Laat W. *Nuclear organization of active and inactive chromatin domains uncovered by chromosome conformation capture-on-chip (4C).* Nature Genetics (2006).  
[publisher](http://dx.doi.org/10.1038/ng1896) | [pubmed](https://www.ncbi.nlm.nih.gov/pubmed/17033623)

<a name="5"></a>5. Schoenfelder S., Sexton T., Chakalova L., Cope N.F., Horton A., Andrews S., Kurukuti S., Mitchell J.A., Umlauf D., Dimitrova D.S., Eskiw C.H., Luo Y., Wei C.L., Ruan Y., Bieker J.J, Fraser P. *Preferential associations between co-regulated genes reveal a transcriptional interactome in erythroid cells.* Nature Genetics (2009).  
[publisher](http://dx.doi.org/10.1038/ng.496) | [pubmed](https://www.ncbi.nlm.nih.gov/pubmed/20010836)

