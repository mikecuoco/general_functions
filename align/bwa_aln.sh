#!/usr/bin/env bash

# for paired end reads: 
#    if longer than ~70bp, use bwa mem
#    if shorter than ~70bp, use bwa aln, then bwa sampe

# attempt to switch to appropriate conda env
env=bw
if ! command -v 'bwa' &>/dev/null && \
  command -v 'conda' &>/dev/null && \
  [ "$CONDA_DEFAULT_ENV" != "$env" ] && \
  conda info --envs | grep "$CONDA_PREFIX_1/envs/$env" &>/dev/null; then
    printf "\e[0;35mAttempting to switch to $env environment \e[0m\n"
    eval "$(conda shell.bash hook)"
    conda activate $env
fi

# input files
read1_files=($(cat $1 | tr "\n" " "))
read2_files=($(cat $2 | tr "\n" " "))
ref=$3
n_files=$(( $(wc -l $1 | awk '{print $1}')-1 ))

# parameters
NM=3 #Number of mismatches allowed in bwa realignment
bwa_i=20 #bwa i parameter prevents indels near the edges of a read


mkdir -p bwa_aln

for i in $(seq 0 $n_files)
do
	# get R1 and R2 files
	fq1="${read1_files[$i]}"
	fq2="${read2_files[$i]}"

	# get basename and prefix
	base=$(basename $fq1)
	prefix="${base%_[0-9].fastq*}"

        echo "aligning $(basename $fq1) to $ref"
	bwa aln -k $NM -n $NM -i $bwa_i -t 4 $ref $fq1 > bwa_aln/${prefix}_1.sai 2> bwa_aln/${prefix}_1_aln.out
        echo "aligning $(basename $fq2) to $ref"
        bwa aln -k $NM -n $NM -i $bwa_i -t 4 $ref $fq2 > bwa_aln/${prefix}_2.sai 2> bwa_aln/${prefix}_2_aln.out
        echo "combining $(basename $fq1) and $(basename fq2) alignments, converting to bam file, sorting"
        bwa sampe -a 1000 -n 100 -N 0 $ref bwa_aln/${prefix}_1.sai bwa_aln/${prefix}_2.sai $fq1 $fq2 | samtools view -h -b -@ 4 | samtools sort -@ 4 > bwa_aln/${prefix}.bam 2> bwa_aln/${prefix}_sampe.out
        samtools index -@ 4 bwa_aln/$prefix.bam bwa_aln/$prefix.bam.bai
        rm -f bwa_aln/${prefix}_1.sai bwa_aln/${prefix}_2.sai
done
