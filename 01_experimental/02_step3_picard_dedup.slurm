#!/bin/bash
#SBATCH --job-name=picard_rice
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=300-00:00:00
#SBATCH --partition=life-zhanghk
#SBATCH --output=picard_%j.out
#SBATCH --error=picard_%j.err

source /datapool/home/2023200496/miniconda3/etc/profile.d/conda.sh
conda activate ssDripPipeline   # make sure this environment has the picard command (or use base)

cd /datapool/life-zhanghk/niulk/my_project/19_modified_analysis/10_iMab/GSE_download/GSE184783_results

mkdir -p 03_dedup

samples=("SRR16077410" "SRR16077411" "SRR16077412" "SRR16077413" "SRR16077414")

for srr in "${samples[@]}"; do
    echo "Removing duplicates for $srr at $(date)"
    
    picard MarkDuplicates \
        REMOVE_DUPLICATES=true \
        METRICS_FILE=03_dedup/${srr}.dup_metrics.txt \
        INPUT=02_align/${srr}.filtered.bam \
        OUTPUT=03_dedup/${srr}.dedup.bam
    
    # Build index (for downstream MACS2)
    samtools index 03_dedup/${srr}.dedup.bam
    
    echo "Finished $srr (dedup BAM: 03_dedup/${srr}.dedup.bam)"
done

echo "All Picard dedup finished at $(date)"i
