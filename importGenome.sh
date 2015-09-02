################################################################################
# UCSC's public ftp repository for genome data
FTP_REPOSITORY="hgdownload.cse.ucsc.edu"
# ------------------------------------------------------------------------------
IDXF="Indexes" # Folder for permanent storage of bowtie indexes
################################################################################
# FUNCTIONS
# //////////////////////////////////////////////////////////////////////////////
# $1 = FTP repository
# $2 = remote path
# $3 = local path
# $4 = ls file
# $5 = regular expression
# ------------------------------------------------------------------------------
remote_ls () {
# Check list of available downloads
ftp -i -v -n "$1" << EOF
user anonymous benjamin.leblanc@bric.ku.dk
lcd "$3"
cd "$2"
ls -1 "$4"
EOF
  # Check if download target is available
  if ! grep -iE "$5" "$4" &> /dev/null; then
    echo "Unable to find $5 in $1/$2"
    rm -f "$4"
  else
    echo $2 >> "remote.path"
  fi
}
# //////////////////////////////////////////////////////////////////////////////
# $1 = path to NOT compressed fasta file
# $2 = path to folder for alignment indexes
# ------------------------------------------------------------------------------
bowtie_index () {
  mkdir -p $2/bowtie
  bowtie-build -f $1 $2/bowtie/`basename $1 '.fa'`
}
################################################################################
# MAIN
# //////////////////////////////////////////////////////////////////////////////
# $1 = genome name
# $2 = local path for referenced genomes (optional)
# ------------------------------------------------------------------------------
GENOME=$1
if [[ "$GENOME" == "" ]]; then
  echo -e "\n"
  echo -e "\tUsage: importGenome.sh <genome_name> [<local_path>]\n"
  echo -e "\tgenome_name = UCSC genome (hg38, mm10, dm6, etc.)"
  echo -e "\tlocal_path  = local path for downloaded genome data, default =  Genome_Data/UCSC_<genome_name>\n"
  exit
fi
# ------------------------------------------------------------------------------
unset LOCAL_PATH
if [[ $2 != "" ]] && [ -e $2 ]; then
  LOCAL_PATH="$2/UCSC_$GENOME"
else
  LOCAL_PATH="Genome_Data/UCSC_$GENOME"
fi
mkdir -p "$LOCAL_PATH"
# ------------------------------------------------------------------------------
CURRENT_WD=`pwd`
cd "$LOCAL_PATH"
# ==============================================================================
# Try to find genome data into the public ftp repository from UCSC
# ------------------------------------------------------------------------------
REMOTE_PATH="/goldenPath/"$GENOME"/bigZips"
# ------------------------------------------------------------------------------
remote_ls $FTP_REPOSITORY $REMOTE_PATH $LOCAL_PATH "remote.genome.ls1" "(chromFa.(zip|tar.gz))|($GENOME.fa.gz)"
if ! mv -f remote.genome.ls* "remote.genome.ls" &> /dev/null; then
  echo -e "\tGenome data for $GENOME not found\n"
  exit
fi
# ------------------------------------------------------------------------------
rm -f "remote.path"
# ------------------------------------------------------------------------------
FASTA_PACKAGE=$(grep -oiE "((chromFa.(zip|tar.gz))|($GENOME.fa.gz))$" "remote.genome.ls")
if [[ $FASTA_PACKAGE == "" ]]; then
  echo -e "\tWhole genome fasta package not found\n"
  exit
fi
# ==============================================================================
# Save information on the source of genome data
echo "HOSTURL   = "$FTP_REPOSITORY""$REMOTE_PATH > "source.genome.info"
echo "PACKAGE   = "$FASTA_PACKAGE >> "source.genome.info"
echo "DATE_TIME = "$(date +'%F %T' | sed 's/ /_/g') >> "source.genome.info"
# ==============================================================================
# Download the whole genome assembly (repeats indicated as lower case)
# ------------------------------------------------------------------------------
if grep "chromFa.zip" "remote.genome.ls" &> /dev/null; then
  rsync -av "rsync://$FTP_REPOSITORY""$REMOTE_PATH/chromFa.zip" "chromFa.zip"
  FASTA_PACKAGE="chromFa.zip"
fi
if grep "chromFa.tar.gz" "remote.genome.ls" &> /dev/null; then
  rsync -av "rsync://$FTP_REPOSITORY""$REMOTE_PATH/chromFa.tar.gz" "chromFa.tar.gz"
  FASTA_PACKAGE="chromFa.tar.gz"
fi
if grep "$GENOME.fa.gz" "remote.genome.ls" &> /dev/null; then
  rsync -av "rsync://$FTP_REPOSITORY""$REMOTE_PATH/$GENOME.fa.gz" "$GENOME.fa.gz"
  FASTA_PACKAGE="$GENOME.fa.gz"
fi
# ==============================================================================
CHRF="Chromosomes" # Folder for temporary storage of genome sequences
mkdir -p $CHRF
mkdir -p $IDXF # Folder for permanent storage of bowtie indexes
# ------------------------------------------------------------------------------
# Unpack genome assembly
if [ -e "chromFa.zip" ]; then
  unzip "chromFa.zip"
fi
if [ -e "chromFa.tar.gz" ]; then
  tar xvzf "chromFa.tar.gz"
fi
if [ -e "$GENOME.fa.gz" ]; then
  gunzip "$GENOME.fa.gz"
  # grep "^>" "$GENOME.fa" | sed -r 's/^>//' | cut -d ' ' -f 1 | sort | uniq > "chr.names"
  # Split whole genome fasta file
  csplit -z -f "$GENOME.csplit." "$GENOME.fa" '/^>/' '{*}'
  for FILE in $(ls "$GENOME.csplit."*); do
    # Rename splited fasta files 
    FNAME="$(head -n 1 $FILE | sed -r 's/^>//' | cut -d ' ' -f 1)".fa
    mv -f $FILE $FNAME
  done
fi
for FILE in $(find * | grep ".fa$"); do
  mv -f $FILE $CHRF
done
if [ -e "$CHRF/$GENOME.fa" ]; then
  mv -f "$CHRF/$GENOME.fa" "$GENOME.fa"
fi
for FILE in $(find * -type d | grep -Ev "^($IDXF|$CHRF)"); do
  rmdir $FILE
done
# ----------------------------------------------------------------------------
cd $CHRF
# ----------------------------------------------------------------------------
# Delete non-canonical chromosomes
# rm -f *_*.fa
# ----------------------------------------------------------------------------
# Extract chromosome names, sizes and fasta labels
for FILE in $(ls *.fa); do
  CHRNAME=$(head -n 1 $FILE | sed -r 's/^>//' | cut -d ' ' -f 1)
  CHRSIZE=$(grep -v "^>" $FILE | wc -m)
  NBLINES=$(grep -v "^>" $FILE | wc -l)
  echo -e "$CHRNAME\t$CHRSIZE\t$NBLINES" >> "../chrom.sizes"
done
ls *.fa | sed 's/_\.*$//' | sed 's/.fa//' | sort | uniq > "../chr.names"
head -n 1 *.fa | grep -v "<==" | cut -d ">" -f 2 > "fasta.labels"
grep -v "^$" "fasta.labels" > "../fasta.labels"
rm "fasta.labels"

# ----------------------------------------------------------------------------
# Concatenate all chromosome sequences into single fasta file
cat *.fa > "genome.fa"
mv -f "genome.fa" "../genome.fa"

# ----------------------------------------------------------------------------
# # Delete alternative haplotype sequences
# rm -f chr*_hap*.fa

# ----------------------------------------------------------------------------
# Remove chromosome files
rm -f *.fa
cd ..
rmdir $CHRF

# ==============================================================================
# Create bowtie indexes
bowtie_index "genome.fa" $IDXF

# //////////////////////////////////////////////////////////////////////////////
# Compress any uncompressed fasta file and remove downloaded package
gzip -f *.fa
# ==============================================================================
rm -f "$FASTA_PACKAGE"
rm remote.genome.ls
cd "$CURRENT_WD"
# ==============================================================================
