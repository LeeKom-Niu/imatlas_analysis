#!/bin/bash
# scripts/run_all_species.sh

GENOME_DIR="/datapool/home/2023200496/niulk/my_project/genome"
OUT_DIR="../results"
mkdir -p $OUT_DIR
OUT_FILE="$OUT_DIR/all_species_region_density.tsv"
rm -f $OUT_FILE  # Clear results before each rerun

# 6 eukaryotic categories
CATS=("fungi" "invertebrate" "plant" "protozoa" "vertebrate_mammalian" "vertebrate_other")

for cat in "${CATS[@]}"; do
    echo "Processing category: $cat ..."
    for sp_dir in "$GENOME_DIR/$cat"/*; do
        if [ -d "$sp_dir/output_genomic" ]; then
            # Run in parallel with 8 threads; adjust -P 8 as needed
            echo "$sp_dir"
        fi
    done | xargs -P 8 -I {} python scripts/calc_one_species.py {} $OUT_FILE
done

echo "All done! Results saved to $OUT_FILE"
