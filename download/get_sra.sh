#!/usr/bin/env bash
# Download SRA files from list
# Input = list of SRA accessions

set -e  

# enable conda
env='data'
if ! command -v 'prefetch' &>/dev/null && \
  command -v 'conda' && \
  [ "$CONDA_DEFAULT_ENV" != $env ] && \
  conda info --envs | grep "$CONDA_PREFIX_1/envs/$env" &>/dev/null; then
    printf "\n\e[0;35m Attempting to switch to $env environment \e[0m\n\n"
    eval "$(conda shell.bash hook)"
    conda activate $env
fi

DIR="$( cd "$( dirname "$1" )" && pwd )"
ACC_FILE=$1
THREADS=$2

mkdir -p raw_fastq # create dir for fastq files
cd raw_fastq

# setup numeric variables
COUNTER=0
TOTAL=$(wc -l ../$ACC_FILE | awk '{print $1}')

while read -r ACC;
do
  # increase counter
  let COUNTER=COUNTER+1 

  # check for existance of file
  if ls ${ACC}*.fastq.gz 1> /dev/null 2>&1;
  then 
     printf "Skipping $ACC since $(ls ${ACC}*.fastq.gz) exists...\n"
  else
          # print in purple
    printf "\n\e[0;35mDownloading and unpacking $ACC ... [$COUNTER/$TOTAL] \e[0m"

    # use command group to store exit code  
    { prefetch -p $ACC && vdb-validate $ACC && fasterq-dump $ACC --split-files -f -e $THREADS -p && gzip $ACC*fastq && rm -rf $ACC; }

    if [ $? -eq 0 ]
    then
      # print message in green
      printf "\e[0;32m[✔] $ACC downloading and unpacking succeeded [$COUNTER/$TOTAL] \e[0m\n"
    else
      # print message in red
      printf "\e[0;31m[✖] $ACC downloading and unpacking failed \e[0m\n"
      echo "$ACC" >> ../failed_acc.txt
    fi
  fi
done < ../$ACC_FILE
