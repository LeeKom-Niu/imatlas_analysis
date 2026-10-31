#!/bin/bash
#SBATCH --job-name=fastp_rice
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=300-00:00:00
#SBATCH --partition=life-zhanghk
#SBATCH --output=fastp_%j.out
#SBATCH --error=fastp_%j.err

# Activate conda environment (adjust to the environment where you installed fastp)
#source /datapool/home/2023200496/miniconda3/etc/profile.d/conda.sh
#conda activate fastp   # or the name of the environment where you installed fastp

# Enter the working directory
cd /datapool/life-zhanghk/niulk/my_project/19_modified_analysis/10_iMab/GSE_download/GSE184783_results

# Create the output directory
mkdir -p 01_trim

# Sample list (5 samples)
samples=("SRR16077410" "SRR16077411" "SRR16077412" "SRR16077413" "SRR16077414")

# Loop over each sample
for srr in "${samples[@]}"; do
    echo "Processing $srr at $(date)"
    /datapool/home/2023200496/envs/fastp/bin/fastp \
        -i 00_raw/${srr}_1.fastq.gz \
        -I 00_raw/${srr}_2.fastq.gz \
        -o 01_trim/${srr}_1.trim.fq.gz \
        -O 01_trim/${srr}_2.trim.fq.gz \
        --thread 16 \
        --json 01_trim/${srr}.json \
        --html 01_trim/${srr}.html
done

echo "All samples processed at $(date)"
