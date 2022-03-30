#!/bin/bash
# Usage
#   download_salk.sh {address} {filename_prefix}

WEB=$1
PRE=$2

# get list of fastqs
curl -s -l $WEB/ | grep -o \>$PRE.*fastq.gz | cut -d\> -f 2 > fastqs.fofn

# downlaod
cat fastqs.fofn | xargs -I {} curl --progress-bar --remote-name-all $WEB/{}
