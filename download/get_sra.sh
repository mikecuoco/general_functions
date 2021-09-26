#!/usr/bin/env bash
#
# Authors: Rohini Gadde & Mike Cuoco
# Usage: ./get_sra.sh $1 $2
#   This script iterates over an accession list to download and convert all
#   SRA files whose accession numbers are in said list
# $1 - Accession list of SRA files
# $2 - Number of threads to use for fasterq-dump (dflt is 6)
# TODO: Automate script more for future use (e.g. create file structure)
# TODO: Automatically disable sra caching
# TODO: Check for correct version of sra-tools

set -e  

# check arguments
if [[ $# -ne 2 ]]; then
  echo "Incorrect number of arguments"
  exit 1
fi

# Set variables
DIR="$( cd "$( dirname "$1" )" && pwd )"
ACC_FILE=$1
THREADS=$2
LOG="../failed_acc.txt"

# enable conda environment if needed
env='data'
if ! command -v 'prefetch' &>/dev/null && \
  command -v 'conda' && \
  [ "$CONDA_DEFAULT_ENV" != $env ] && \
  conda info --envs | grep "$CONDA_PREFIX_1/envs/$env" &>/dev/null; then
    printf "\n\e[0;35m Attempting to switch to $env environment \e[0m\n\n"
    eval "$(conda shell.bash hook)"
    conda activate $env
fi

# create dir for fastq files and enter
mkdir -p raw_fastq 
cd raw_fastq

# Clear log before each run
if [[ -e $LOG ]]; then
  rm $LOG
fi

# setup numeric variables
COUNTER=0
TOTAL=$(wc -l ../$ACC_FILE | awk '{print $1}')

while read -r ACC; do
  # increase counter
  let COUNTER=COUNTER+1 

  # check for existance of file
  if [[ -e ${ACC}*.fastq.gz ]]; then 
     printf "Skipping $ACC since $(ls ${ACC}*.fastq.gz) exists...\n"
  else
    # print in purple
    printf "\n\e[0;35mDownloading and unpacking $ACC ... [$COUNTER/$TOTAL] \e[0m"

    # use command group to store exit code  
    { prefetch -p $ACC && vdb-validate $ACC && fasterq-dump $ACC --split-files -f -e $THREADS -p && gzip $ACC*fastq && rm -rf $ACC; }

    if [ $? -eq 0 ]; then
      # print message in green
      printf "\e[0;32m[✔] $ACC downloading and unpacking succeeded [$COUNTER/$TOTAL] \e[0m\n"
    else
      # print message in red
      printf "\e[0;31m[✖] $ACC downloading and unpacking failed \e[0m\n"
      echo "$ACC" >> $LOG
    fi
  fi
done < ../$ACC_FILE
