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
  'genome.id',            'i', 1, 'character', "identifier of the genome assembly",   "",
  'genome.index.path',    'p', 1, 'character', "path to genome index for bowtie",                    "",
  'original.design.path', 'o', 1, 'character', "path to the original array design file",             "",
  'updated.design.path',  'u', 2, 'character', "path to folder for storage of updated design data, default = Processed_Data/Updated/", "Processed_Data/Updated/",
  'bowtie.mismatches',    'n', 2, 'integer',   "number of allowed mismatches (bowtie mapping), default = 0",     0,
  'bowtie.threads',       't', 2, 'integer',   "number of CPU threads (bowtie mapping), default = 6",            6
), byrow=TRUE, ncol=6)
Cfg <- processArgs(script.args)
Cfg$script.args <- NULL
# -----------------------------------------------------------------------------.
# DEBUG VALUES                                          # *** WARNING/TODO *** #
if(is.na(Cfg$script.name)) {
  Cfg$genome.id            <- "dm6"
  Cfg$genome.index.path    <- "Genome_Data/UCSC_dm6/Indexes/bowtie/genome"
  Cfg$original.design.path <- "Raw_Data/GSE23887_RAW/GPL10867.ndf.gz"

  Cfg$updated.design.path  <- "Processed_Data/Updated/"
  Cfg$bowtie.mismatches <- 0
  Cfg$bowtie.threads    <- 6
}
# =============================================================================.
# Parameters validation
lbl.lst <- c("genome.id", "genome.index.path", "original.design.path")
for(lbl in lbl.lst) {
  if(Cfg[[lbl]]=="") {
    msg <- paste("Missing", lbl, "get help using: Rscript", Cfg$script.name, "-h")
    stop(msg, "\n\n")
  }
}
# -----------------------------------------------------------------------------.
attach(Cfg)
# -----------------------------------------------------------------------------.
# Verify that bowtie index path is correct
verifyInputFiles(paste(genome.index.path, ".1.ebwt", sep=""))
# -----------------------------------------------------------------------------.
# Verify that array design file exists and has required columns
array.design.format <- file.formats$nimblegen.ndf
verifyInputFiles(original.design.path, array.design.format$colnames)
# =============================================================================.
# Parameter updates
cmd <- paste('mkdir -p', updated.design.path)
system(cmd, intern=F)

design.id <- gsub(
  "[.](ndf|txt)([.]gz)?$", "", basename(original.design.path), perl=T, ignore.case=T
)
probes.file.path     <- paste(updated.design.path, "/", design.id, "_probe_sequences.fa", sep="")
mapped.probes.path   <- paste(updated.design.path, "/", design.id, "_", genome.id, "_aligned_probes.bowtie", sep="")
mapped.probes.format <- file.formats$mapped.probe.sequences
updated.design.path  <- paste(updated.design.path, "/", design.id, "_", genome.id, ".txt", sep="")

# PROCESSING ###################################################################

# -----------------------------------------------------------------------------.
# Load microarray design
design <- readData(original.design.path, array.design.format)

# Remove array layout and control probes
discard <- design$CONTAINER %in% c("V_CODE", "H_CODE", "NGS_CONTROLS", "RANDOM")
design <- design[! discard,]
# -----------------------------------------------------------------------------.
# Create a FASTA file with microarray probe sequences
i.probes <- 1:nrow(design)
n.probes <- length(i.probes)
lines <- rep("", 2*n.probes)
lines[2*(1:n.probes) - 1] <- paste(">", i.probes,"\t",design$PROBE_ID[i.probes], sep="")
lines[2*(1:n.probes)]     <- design$PROBE_SEQUENCE[i.probes]
probes.file <- file(probes.file.path, "w")
writeLines(lines, con = probes.file)
close(probes.file)

# -----------------------------------------------------------------------------.
# Realign probe sequences to latest genome release using bowtie
cmd <- paste(
  "bowtie -f -n", bowtie.mismatches, "-m 1 --strata --best -p", bowtie.threads,
  genome.index.path, probes.file.path, mapped.probes.path,
  sep=" "
)
txt <- system(cmd, intern=T)

# =============================================================================.
# Update microarray design with realigned probes
# -----------------------------------------------------------------------------.
# Main columns of the output from bowtie:
# 1: probe index
# 2: PROBE_ID
# 3: strand
# 4: chromosome/seq_id
# 5: position
# 6: PROBE_SEQUENCE
# -----------------------------------------------------------------------------.
probes <- readData(mapped.probes.path, mapped.probes.format)
# Check consistency 
i <- probes[,1]
if(sum(design$PROBE_ID[i]!=probes[,2])!=0) stop("bad probe alignement file")
# -----------------------------------------------------------------------------.
l <-   nchar(design$PROBE_SEQUENCE[i])
chr <- gsub("(^[ \t\n]+)|([ \t\n]+$)", "", (as.character(probes[,4])))
updated.design <- cbind(
  as.character(design$PROBE_ID[i]),
  as.character(design$PROBE_CLASS[i]),
  as.character(design$CONTAINER[i]),
  as.integer(design$X[i]),
  as.integer(design$Y[i]),
  as.character(design$SEQ_ID[i]),  # chromosome (seqname) from original design
  chr,   						               # chromosome (seqname) from bowtie alignment
  as.character(probes[,3]),        # strand
  as.integer(probes[,5] + 1),      # start
  as.integer(probes[,5] + l - 1),  # end
  as.integer(l)                    # length
)
colnames(updated.design) <- file.formats$nimblegen.updated$table$columns
# -----------------------------------------------------------------------------.
write.table(format(updated.design, scientific=FALSE, digits=11, trim=TRUE, justify='none'), file=updated.design.path , sep="\t", row.names=FALSE, quote=FALSE)
# =============================================================================.
# Cleanup alignment files and compress the updated design files
# rm Processed_Data/Updated/*_probe_sequences.fa
cmd <- paste("rm -f", probes.file.path)
system(cmd, intern=F)
# -----------------------------------------------------------------------------.
# rm Processed_Data/Updated/*_aligned_probes.bowtie
cmd <- paste("rm -f", mapped.probes.path)
system(cmd, intern=F)
# -----------------------------------------------------------------------------.
# gzip Processed_Data/Updated/*.txt
cmd <- paste("gzip", updated.design.path)
system(cmd, intern=F)
# =============================================================================.
detach(Cfg)
