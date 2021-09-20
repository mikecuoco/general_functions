#!/usr/bin/env bash
# Usage:
#  1 = list of bam files (file of file names)
#  2 = GTF file
#  3 = Threads
#  4 = Outfile name

env=shell
if ! command -v 'featureCounts' &>/dev/null && \
  command -v 'conda' &>/dev/null && \
  [ "$CONDA_DEFAULT_ENV" != "$env" ] && \
  conda info --envs | grep "$CONDA_PREFIX_1/envs/$env" &>/dev/null; then
    printf "\e[0;35mAttempting to switch to $env environment \e[0m\n"
    eval "$(conda shell.bash hook)"
    conda activate $env
fi

# setup inputs
readarray -t BAMFILES < $1
GTF=$2
THREADS=$3
OUTFILE=$4

# run command


featureCounts -T $THREADS -a "$GTF" -t exon -g transcript_id -fLMO --fracOverlap 0.8 --ignoreDup --verbose -o "$OUTFILE" "${BAMFILES[@]}"
echo "featureCounts -T $THREADS -a $GTF -t exon -g transcript_id -fLMO --fracOverlap 0.8 --ignoreDup --verbose -o $OUTFILE ${BAMFILES[@]}"
