#!/usr/bin/env bash
# Usage:
#  1 = fastq file
#  2 = reference genome

# attempt to switch proper environment
env=mm2
if ! command -v 'minimap2' &>/dev/null && \
  command -v 'conda' &>/dev/null && \
  [ "$CONDA_DEFAULT_ENV" != "$env" ] && \
  conda info --envs | grep "$CONDA_PREFIX_1/envs/$env" &>/dev/null; then
    printf "\e[0;35mAttempting to switch to $env environment \e[0m\n"
    eval "$(conda shell.bash hook)"
    conda activate $env
fi

# setup file
file=$1
base=$(basename $file)

# use this reference
ref=$2
idx="${ref%.*}.mmi"

# make prefix
if [ "$base" == *".gz"* ]; then
  prefix="${base%.*.*}"
else
  prefix="${base%.*}"
fi 

# align
echo "aligning $base to $ref"
minimap2 -zx map-hifi --secondary=yes -N 100 $idx $file | samtools view -h -b -@ 4 | samtools sort -@ 4 > $prefix.bam 2> $prefix.out
samtools index $prefix.bam $prefix.bam.bai

