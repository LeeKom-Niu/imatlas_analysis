#!/bin/bash
#SBATCH --job-name=intersect_rice
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --partition=life-zhanghk
#SBATCH --output=intersect_%j.out
#SBATCH --error=intersect_%j.err

source /datapool/home/2023200496/miniconda3/etc/profile.d/conda.sh
conda activate base   # make sure this environment contains bedtools

cd /datapool/life-zhanghk/niulk/my_project/19_modified_analysis/10_iMab/GSE_download/GSE184783_results

mkdir -p 05_intersect

# Shared peaks with the Input control
echo "Intersecting peaks (IP vs Input) at $(date)"
bedtools intersect \
    -a 04_peak/SRR16077410_vs_Input_peaks.narrowPeak \
    -b 04_peak/SRR16077411_vs_Input_peaks.narrowPeak \
       04_peak/SRR16077412_vs_Input_peaks.narrowPeak \
    -u > 05_intersect/IP_vs_Input_intersect.narrowPeak

# Shared peaks with the IgG control
echo "Intersecting peaks (IP vs IgG) at $(date)"
bedtools intersect \
    -a 04_peak/SRR16077410_vs_IgG_peaks.narrowPeak \
    -b 04_peak/SRR16077411_vs_IgG_peaks.narrowPeak \
       04_peak/SRR16077412_vs_IgG_peaks.narrowPeak \
    -u > 05_intersect/IP_vs_IgG_intersect.narrowPeak

# Summary statistics
echo ""
echo "=== Peak counts ==="
echo "IP_vs_Input_intersect: $(wc -l < 05_intersect/IP_vs_Input_intersect.narrowPeak)"
echo "IP_vs_IgG_intersect:   $(wc -l < 05_intersect/IP_vs_IgG_intersect.narrowPeak)"

echo "All intersection finished at $(date)"
