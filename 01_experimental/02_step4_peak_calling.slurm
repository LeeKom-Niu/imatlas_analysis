#!/bin/bash
#SBATCH --job-name=macs2_rice
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=300-00:00:00
#SBATCH --partition=life-zhanghk
#SBATCH --output=macs2_%j.out
#SBATCH --error=macs2_%j.err

source /datapool/home/2023200496/miniconda3/etc/profile.d/conda.sh
conda activate ssDripPipeline   # make sure this environment contains macs2

cd /datapool/life-zhanghk/niulk/my_project/19_modified_analysis/10_iMab/GSE_download/GSE184783_results

mkdir -p 04_peak

# IP samples (3 biological replicates)
ip_samples=("SRR16077410" "SRR16077411" "SRR16077412")
# Control samples
input_sample="SRR16077413"
igg_sample="SRR16077414"

# 1. IP vs Input
for srr in "${ip_samples[@]}"; do
    echo "Running MACS2: $srr vs Input at $(date)"
    macs2 callpeak \
        -t 03_dedup/${srr}.dedup.bam \
        -c 03_dedup/${input_sample}.dedup.bam \
        -f BAMPE \
        -g 373128865 \
        -p 1e-4 \
        --nomodel \
        -n ${srr}_vs_Input \
        --outdir 04_peak
done

# 2. IP vs IgG
for srr in "${ip_samples[@]}"; do
    echo "Running MACS2: $srr vs IgG at $(date)"
    macs2 callpeak \
        -t 03_dedup/${srr}.dedup.bam \
        -c 03_dedup/${igg_sample}.dedup.bam \
        -f BAMPE \
        -g 373128865 \
        -p 1e-4 \
        --nomodel \
        -n ${srr}_vs_IgG \
        --outdir 04_peak
done

echo "All MACS2 peak calling finished at $(date)"
