#!/bin/bash
# Complete iM-CUT&Tag data processing pipeline (GSE220878)
# Function: bigwig → bedgraph → SEACR peak calling → intersect replicates
# Author: reproduced following the methods in Zanin et al. NAR 2023
# Date: 2026-08-31

set -e  # stop immediately on error

echo "========================================"
echo "开始处理 GSE220878 iM-CUT&Tag 数据"
echo "========================================"

# 1. Enter the data directory
cd ~/niulk/my_project/19_modified_analysis/10_iMab/GSE_download/GSE220878 || exit

# 2. Create the results directory (if it does not exist)
mkdir -p ../GSE220878_results

# 3. Define the results directory path
RESULT_DIR="../GSE220878_results"

echo ""
echo "[步骤1] 批量转换 bigwig → bedgraph，并运行 SEACR..."

# 4. Loop over the 6 iM samples (3 HEK293T + 3 WDLPS)
for bw in GSM7507863_HEK_iM_rep1.bigwig \
          GSM7507864_HEK_iM_rep2.bigwig \
          GSM7507865_HEK_iM_rep3.bigwig \
          GSM7507869_WDLPS_iM_rep1.bigwig \
          GSM7507870_WDLPS_iM_rep2.bigwig \
          GSM7507871_WDLPS_iM_rep3.bigwig; do
    
    # Strip the .bigwig suffix to get the base name
    base=${bw%.bigwig}
    
    echo "  -> 处理 ${base} ..."
    
    # Convert to bedgraph (intermediate file)
    bigWigToBedGraph $bw ${base}.bedgraph
    
    # SEACR peak calling (parameters: 0.01 non stringent, consistent with the authors)
    ./SEACR_1.3.sh ${base}.bedgraph 0.01 non stringent ${base}_SEACR
    
    # Move the final BED to the results directory (rename, drop the _SEACR suffix)
    mv ${base}_SEACR.stringent.bed ${RESULT_DIR}/${base}.bed
    
    # Clean up the intermediate bedgraph file
    rm ${base}.bedgraph
done

echo "[步骤1] 完成！6 个样本的 BED 文件已生成。"

echo ""
echo "[步骤2] 分别取交集，获得高置信度 iM 峰..."

# 5. Enter the results directory
cd ${RESULT_DIR} || exit

# 6. Intersect the three HEK293T replicates
bedtools intersect -a GSM7507863_HEK_iM_rep1.bed \
                   -b GSM7507864_HEK_iM_rep2.bed \
                      GSM7507865_HEK_iM_rep3.bed \
                   -u > HEK293T_iM_highConf_intersect.bed

echo "  -> HEK293T 交集完成: $(wc -l < HEK293T_iM_highConf_intersect.bed) 个峰"

# 7. Intersect the three WDLPS replicates
bedtools intersect -a GSM7507869_WDLPS_iM_rep1.bed \
                   -b GSM7507870_WDLPS_iM_rep2.bed \
                      GSM7507871_WDLPS_iM_rep3.bed \
                   -u > WDLPS_iM_highConf_intersect.bed

echo "  -> WDLPS 交集完成: $(wc -l < WDLPS_iM_highConf_intersect.bed) 个峰"

echo ""
echo "[步骤2] 完成！"
echo ""
echo "========================================"
echo "所有 iM 数据处理完成！结果文件列表："
echo "========================================"
ls -lh *.bed | awk '{print "  " $9 " (" $5 ")"}'
echo ""
echo "高置信度峰（交集）:"
echo "  - HEK293T_iM_highConf_intersect.bed"
echo "  - WDLPS_iM_highConf_intersect.bed"
echo ""
echo "输出目录: ${RESULT_DIR}"
