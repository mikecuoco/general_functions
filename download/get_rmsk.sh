#!/usr/bin/env bash

env='shell'
# activate conda environment
if ! command -v 'gtf2bed' &>/dev/null && \
  command -v 'conda' && \
  [ "$CONDA_DEFAULT_ENV" != $env ] && \
  conda info --envs | grep "$CONDA_PREFIX_1/envs/$env" &>/dev/null; then
    printf "\n\e[0;35m Attempting to switch to $env environment \e[0m\n\n"
    eval "$(conda shell.bash hook)"
    conda activate $env
fi

# assign global variables
GENOME=$1
CACHE=~/.cache/annotations

#  downlaod rmsk file if needed
URL="https://hgdownload.cse.ucsc.edu/goldenPath/${GENOME}/database/rmsk.txt.gz"
RMSK="${GENOME}_rmsk.txt"

# create cache directory and code
mkdir -p $CACHE
mkdir -p ~/code/bin

if [ ! -f $CACHE/$RMSK ]; then
    echo "Downloading $RMSK from $URL"
    curl -s $URL > $CACHE/$RMSK.gz && gunzip $CACHE/$RMSK.gz || echo "Downloading $RMSK failed"
fi

makeGTF=~/code/bin/makeTEgtf.pl

# Make into GTF file
GTF="${GENOME}_rmsk.gtf"
if [ ! -f $CACHE/$GTF ]; then
    echo "making $GTF from $RMSK"
    perl $makeGTF -c 6 -s 7 -e 8 -o 10 -n ${GENOME}_rmsk -t 11 -C 12 -f 13 -S 2 ${CACHE}/${RMSK} > ${CACHE}/${GTF} || echo "$RMSK GTF conversion failed"
fi

# Make into BED file
BED="${GENOME}_rmsk.bed"
if [ ! -f $CACHE/$BED ]; then
    echo "making $BED from $GTF"
    gtf2bed --input=gtf < $CACHE/$GTF > $CACHE/$BED || echo "$RMSK BED conversion failed"
fi

# Make into GFF file
GFF="${GENOME}_rmsk.gff"
if [ ! -f $CACHE/$GFF ]; then
    echo "making $GFF from $GTF"
    gffread -E $CACHE/$GTF -o- > $CACHE/$GFF  || echo "$RMSK GFF conversion failed"
fi

# Make symlinks
# TODO: make this into a loop
echo "Linking $CACHE/$RMSK to $RMSK"
ln -fs $CACHE/$RMSK $RMSK
echo "Linking $CACHE/$GTF to $GTF"
ln -fs $CACHE/$GTF $GTF 
echo "Linking $CACHE/$BED to $BED"
ln -fs $CACHE/$BED $BED
echo "Linking $CACHE/$GFF to $GFF"
ln -fs $CACHE/$GFF $GFF
