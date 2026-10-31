#!/bin/bash
#SBATCH --job-name=bwa_filter_rice
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=300-00:00:00
#SBATCH --partition=life-zhanghk
#SBATCH --output=bwa_filter_%j.out
#SBATCH --error=bwa_filter_%j.err

# Environment
source /datapool/home/2023200496/miniconda3/etc/profile.d/conda.sh
conda activate ssDripPipeline   # make sure this environment contains bwa and samtools

# Configuration
cd /datapool/life-zhanghk/niulk/my_project/19_modified_analysis/10_iMab/GSE_download/GSE184783_results
REF_FASTA="ref/rice_chr1_12.fna"
THREADS=16
samples=("SRR16077410" "SRR16077411" "SRR16077412" "SRR16077413" "SRR16077414")

mkdir -p 02_align

# 1. Build the BWA index (if not already built)
if [ ! -f ${REF_FASTA}.bwt ]; then
    echo "Building BWA index..."
    /datapool/home/2023200496/envs/bwa/bin/bwa index ${REF_FASTA}
fi

# 2. Alignment + filtering
for srr in "${samples[@]}"; do
    echo "Processing $srr at $(date)"
    
    # bwa mem (with -M)
    /datapool/home/2023200496/envs/bwa/bin/bwa mem -t $THREADS -M $REF_FASTA \
        01_trim/${srr}_1.trim.fq.gz \
        01_trim/${srr}_2.trim.fq.gz \
        > 02_align/${srr}.sam
    
    # SAM → BAM (compressed)
    samtools view -@ $THREADS -bS -1 02_align/${srr}.sam -o 02_align/${srr}.bam
    rm 02_align/${srr}.sam
    
    # Sort
    samtools sort -@ $THREADS -m 12G -l 9 02_align/${srr}.bam -o 02_align/${srr}.sorted.bam
    rm 02_align/${srr}.bam
    
    # Filter: keep properly paired reads (-f 2) with MAPQ >= 10 (-q 10)
    samtools view -f 2 -q 10 -b 02_align/${srr}.sorted.bam > 02_align/${srr}.filtered.bam
    02_align/${srr}.sorted.bam
    
    # Build index (needed by Picard and downstream steps)
    samtools index 02_align/${srr}.filtered.bam
    
    echo "Finished $srr (filtered BAM saved: 02_align/${srr}.filtered.bam)"
done

echo "All alignment & filtering done at $(date)"
