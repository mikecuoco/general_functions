#!/usr/bin/env bash
# Download from ENA

set -e  

# enable conda
env='data'
# activate conda environment
if ! command -v 'wget' &>/dev/null && \
  command -v 'conda' && \
  [ "$CONDA_DEFAULT_ENV" != $env ] && \
  conda info --envs | grep "$CONDA_PREFIX_1/envs/$env" &>/dev/null; then
    printf "\n\e[0;35m Attempting to switch to $env environment \e[0m\n\n"
    eval "$(conda shell.bash hook)"
    conda activate $env
fi

while read line  
do  
  FILE=$(basename $line)
  if [ ! -f "$FILE" ]; then
    aws s3 cp $line . --no-sign-request
  fi
done < $1

