#!/bin/bash
# Complete iM-CUT&Tag data processing pipeline (GSE220878)
# Function: bigwig → bedgraph → SEACR peak calling → intersect replicates
# Author: reproduced following the methods in Zanin et al. NAR 2023
# Date: 2026-08-31

set -e  # stop immediately on error

echo "========================================"
echo "Starting GSE220878 iM-CUT&Tag data processing"
echo "========================================"

# 1. Enter the data directory
cd ~/niulk/my_project/19_modified_analysis/10_iMab/GSE_download/GSE220878 || exit

# 2. Create the results directory (if it does not exist)
mkdir -p ../GSE220878_results

# 3. Define the results directory path
RESULT_DIR="../GSE220878_results"

echo ""
echo "[Step 1] Batch convert bigwig → bedgraph and run SEACR..."

# 4. Loop over the 6 iM samples (3 HEK293T + 3 WDLPS)
for bw in GSM7507863_HEK_iM_rep1.bigwig \
          GSM7507864_HEK_iM_rep2.bigwig \
          GSM7507865_HEK_iM_rep3.bigwig \
          GSM7507869_WDLPS_iM_rep1.bigwig \
          GSM7507870_WDLPS_iM_rep2.bigwig \
          GSM7507871_WDLPS_iM_rep3.bigwig; do
    
    # Strip the .bigwig suffix to get the base name
    base=${bw%.bigwig}
    
    echo "  -> Processing ${base} ..."
    
    # Convert to bedgraph (intermediate file)
    bigWigToBedGraph $bw ${base}.bedgraph
    
    # SEACR peak calling (parameters: 0.01 non stringent, consistent with the authors)
    ./SEACR_1.3.sh ${base}.bedgraph 0.01 non stringent ${base}_SEACR
    
    # Move the final BED to the results directory (rename, drop the _SEACR suffix)
    mv ${base}_SEACR.stringent.bed ${RESULT_DIR}/${base}.bed
    
    # Clean up the intermediate bedgraph file
    rm ${base}.bedgraph
done

echo "[Step 1] Done! BED files generated for 6 samples."

echo ""
echo "[Step 2] Intersect replicates to obtain high-confidence iM peaks..."

# 5. Enter the results directory
cd ${RESULT_DIR} || exit

# 6. Intersect the three HEK293T replicates
bedtools intersect -a GSM7507863_HEK_iM_rep1.bed \
                   -b GSM7507864_HEK_iM_rep2.bed \
                      GSM7507865_HEK_iM_rep3.bed \
                   -u > HEK293T_iM_highConf_intersect.bed

echo "  -> HEK293T intersection done: $(wc -l < HEK293T_iM_highConf_intersect.bed) peaks"

# 7. Intersect the three WDLPS replicates
bedtools intersect -a GSM7507869_WDLPS_iM_rep1.bed \
                   -b GSM7507870_WDLPS_iM_rep2.bed \
                      GSM7507871_WDLPS_iM_rep3.bed \
                   -u > WDLPS_iM_highConf_intersect.bed

echo "  -> WDLPS intersection done: $(wc -l < WDLPS_iM_highConf_intersect.bed) peaks"

echo ""
echo "[Step 2] Done!"
echo ""
echo "========================================"
echo "All iM data processing complete! Result file list:"
echo "========================================"
ls -lh *.bed | awk '{print "  " $9 " (" $5 ")"}'
echo ""
echo "High-confidence peaks (intersection):"
echo "  - HEK293T_iM_highConf_intersect.bed"
echo "  - WDLPS_iM_highConf_intersect.bed"
echo ""
echo "Output directory: ${RESULT_DIR}"
