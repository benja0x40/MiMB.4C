#!/bin/bash
################################################################################
# MAIN
# //////////////////////////////////////////////////////////////////////////////
# $1 = GEO accession number
# ------------------------------------------------------------------------------
ACCESSION=$1
if [[ "$ACCESSION" == "" ]]; then
  echo -e "\n"
  echo -e "\tUsage: importRawData.sh <accession> [<local_path>]\n"
  echo -e "\taccession  = GEO series accession number (GSE prefixed identifier)"
  echo -e "\tlocal_path = local path for downloaded data, default = Raw_Data\n"
  exit
fi
# ------------------------------------------------------------------------------
unset LOCAL_PATH
if [[ $2 != "" ]] && [ -e $2 ]; then
  LOCAL_PATH="$2/"
else
  LOCAL_PATH="Raw_Data/"
fi
mkdir -p "$LOCAL_PATH"
# ------------------------------------------------------------------------------
CURRENT_WD=`pwd`
cd "$LOCAL_PATH"
# ==============================================================================
# Try to find genome data into the public ftp repository from UCSC
# ------------------------------------------------------------------------------
REMOTE_PATH="http://www.ncbi.nlm.nih.gov/geo/download/?acc=$ACCESSION&format=file"
wget -c "$REMOTE_PATH" -O "$ACCESSION""_RAW.tar"
# ------------------------------------------------------------------------------
mkdir -p "$ACCESSION""_RAW"
tar xvf "$ACCESSION""_RAW.tar" -C "$ACCESSION""_RAW"
# ==============================================================================
cd "$CURRENT_WD"
# ==============================================================================

exit
