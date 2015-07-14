# LIBRARIES ####################################################################

message("Loading R/Bioconductor packages")
suppressPackageStartupMessages(library("MRA.4C"))

# SCRIPT PARAMETERS ############################################################

# =============================================================================.
script.args = matrix(c(
  # Column 1: long option name
  # Column 2: short option name 
  # Column 3: 0=no argument, 1=required argument, 2=optional argument
  # Column 4: data type (logical, integer, double, complex, character)
  # Column 5: a brief description of the purpose of the option
  # Column 6: default value
  'enzyme.name',           'n', 1, 'character', "name of the restriction enzyme",      "",
  'enzyme.motif',          'm', 1, 'character', "DNA motif of the restriction enzyme", "",
  'genome.sequence.path',  's', 1, 'character', "path to genome sequence (fasta)",     "",
  'genome.id',             'i', 1, 'character', "identifier of the genome assembly",   "",
  'restriction.data.path', 'r', 2, 'character', "path to folder for storage of restriction data, default = Processed_Data/RestrictionMap/", "Processed_Data/RestrictionMap/"
), byrow=TRUE, ncol=6)
Cfg <- processArgs(script.args)
Cfg$script.args <- NULL
# -----------------------------------------------------------------------------.
# DEBUG VALUES                                          # *** WARNING/TODO *** #
if(is.na(Cfg$script.name)) {
  Cfg$genome.id    <- "dm6"
  Cfg$enzyme.name  <- "DpnII"
  Cfg$enzyme.motif <- "GATC"
  Cfg$genome.sequence.path  <- "Genome_Data/UCSC_dm6/genome.fa.gz"
  Cfg$restriction.data.path <- "Processed_Data/RestrictionMap/"
}
# =============================================================================.
# Parameters validation
lbl.lst <- c("genome.id", "enzyme.name", "enzyme.motif", "genome.sequence.path")
for(lbl in lbl.lst) {
  if(Cfg[[lbl]]=="") {
    msg <- paste("Missing", lbl, "get help using: Rscript", Cfg$script.name, "-h")
    stop(msg, "\n\n")
  }
}
# -----------------------------------------------------------------------------.
attach(Cfg)
# -----------------------------------------------------------------------------.
# Verify that genome sequence file exists
verifyInputFiles(genome.sequence.path)
# =============================================================================.
# Parameter updates
enzyme <- list(name=enzyme.name, motif=enzyme.motif)

cmd <- paste('mkdir -p', restriction.data.path)
system(cmd, intern=F)

out.file.template      <- paste(restriction.data.path, "/", genome.id, "_", enzyme$name, sep="")

# PROCESSING ###################################################################

# Load genome sequence
genome.seq <- readDNAStringSet(genome.sequence.path, "fasta")
# Filter out non canonical sequences
seqlist <- names(genome.seq)
seqlist <- seqlist[! grepl("(^chrUn_)|(_random$)", seqlist, ignore.case=T)]
# Compute restriction sites and restriction fragments
res <- computeRestrictionMap(
  genome.seq, enzyme$motif, output.file=out.file.template, seqlist=seqlist
)
# =============================================================================.
# Compress restriction map files
# gzip Processed_Data/RestrictionMap/*.txt
cmd <- paste("gzip", res$sites)
system(cmd, intern=F)
cmd <- paste("gzip", res$fragments)
system(cmd, intern=F)
# =============================================================================.
detach(Cfg)


