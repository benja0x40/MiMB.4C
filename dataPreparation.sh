# =============================================================================.
GENOME="dm6"
ENZYME="DpnII"
CUTDNA="GATC"
# =============================================================================.
#Import data
# -----------------------------------------------------------------------------.
# Import Drosophila genome sequence from UCSC and create indexes for alignment using bowtie
./importGenome.sh "$GENOME"
# Import Bantignies et al., 2011 4C data from GEO
./importRawData.sh "GSE23887"

# =============================================================================.
# Update design data to latest genome release
# -----------------------------------------------------------------------------.
Rscript updateDesignData.R -i "$GENOME" -p "Genome_Data/UCSC_$GENOME/Indexes/bowtie/genome" -o "Raw_Data/GSE23887_RAW/GPL10867.ndf.gz"

# =============================================================================.
# Select only relevant probes within the array design
# -----------------------------------------------------------------------------.
mkdir -p Processed_Data/Cleaned/
FILTER="\texperimental\tISOTHERMAL\t"
gunzip -c "Processed_Data/Updated/GPL10867_$GENOME.txt.gz" | head -n 1 > "Processed_Data/Cleaned/GPL10867_$GENOME.txt" # copy column names
gunzip -c "Processed_Data/Updated/GPL10867_$GENOME.txt.gz" | grep -P "$FILTER" >> "Processed_Data/Cleaned/GPL10867_$GENOME.txt"
# Compress updated design files
gzip Processed_Data/Cleaned/*.txt

# =============================================================================.
# Align 4C bait sequences to latest genome release
# -----------------------------------------------------------------------------.
bowtie -f -n 1 -p 6 "Genome_Data/UCSC_$GENOME/Indexes/bowtie/genome" "4C_Bait_Sequences/4C_Primers.fa" "4C_Bait_Sequences/4C_Primers_Aligned_$GENOME.txt"

# =============================================================================.
# Compute restriction map and match array probes to restriction fragments
# -----------------------------------------------------------------------------.
Rscript computeRestrictionMap.R -i "$GENOME" -n "$ENZYME" -m "$CUTDNA" -s "Genome_Data/UCSC_$GENOME/genome.fa.gz"


