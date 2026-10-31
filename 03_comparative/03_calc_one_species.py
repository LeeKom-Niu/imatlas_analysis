#!/usr/bin/env python3
import os
import sys
import subprocess
import pandas as pd
import numpy as np

def run_cmd(cmd, ignore_error=False):
    res = subprocess.run(cmd, shell=True, check=False if ignore_error else True, 
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0 and not ignore_error:
        print(f"\n[ERROR] Command failed: {cmd}\n{res.stderr}", flush=True)
        raise subprocess.CalledProcessError(res.returncode, cmd)
    return res

def process_species(sp_path, out_dir):
    sp_name = os.path.basename(sp_path)
    if sp_name == 'output_genomic': 
        sp_name = os.path.basename(os.path.dirname(sp_path))
    
    gdir = os.path.join(sp_path, 'output_genomic')
    if not os.path.isdir(gdir):
        gdir = sp_path

    fai = os.path.join(gdir, 'genome.fai')
    prom = os.path.join(gdir, 'promoter_1kb.bed')
    gene = os.path.join(gdir, 'gene_with_strand.bed')
    imotif = os.path.join(gdir, 'imotif_with_strand.bed')
    
    if not all(os.path.exists(f) for f in [fai, prom, gene, imotif]):
        return None

    print(f"Processing: {sp_name}", flush=True)
    tmp = os.path.join(out_dir, f'tmp_{sp_name}')
    os.makedirs(tmp, exist_ok=True)

    # 1. Read fai, build chromosome length dictionary and whole-genome BED
    chr_len = {}
    genome_bed = os.path.join(tmp, 'genome.bed')
    with open(fai) as f_in, open(genome_bed, 'w') as f_out:
        for line in f_in:
            parts = line.strip().split('\t')
            chrom, length = parts[0], int(parts[1])
            chr_len[chrom] = length
            f_out.write(f"{chrom}\t0\t{length}\n")

    # 2. Define a function to safely read BED and truncate coordinates
    def get_safe_bed(in_bed, out_bed, mode='gene_or_prom'):
        out_lines = []
        with open(in_bed) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 3: continue
                chrom = parts[0]
                if chrom not in chr_len: continue
                start, end = int(parts[1]), int(parts[2])
                strand = parts[5] if len(parts) >= 6 else '+'
                
                # Truncate raw input (prevent negative or out-of-range coordinates)
                start = max(0, start)
                end = min(chr_len[chrom], end)
                
                if mode == 'downstream':
                    if strand == '+':
                        start = end
                        end = min(chr_len[chrom], end + 1000)
                    else:
                        end = start
                        start = max(0, start - 1000)
                
                if start < end:
                    out_lines.append(f"{chrom}\t{start}\t{end}\n")
        
        with open(out_bed, 'w') as f_out:
            f_out.writelines(out_lines)

    # 3. Generate safe BED files
    prom_safe = os.path.join(tmp, 'prom_safe.bed')
    gene_safe = os.path.join(tmp, 'gene_safe.bed')
    down_safe = os.path.join(tmp, 'down_safe.bed')
    
    get_safe_bed(prom, prom_safe, mode='gene_or_prom')
    get_safe_bed(gene, gene_safe, mode='gene_or_prom')
    get_safe_bed(gene, down_safe, mode='downstream')

    # 4. Consistently use bedtools sort -g {fai} + merge to compute lengths
    def sort_merge(in_bed, out_bed):
        sort_bed = in_bed + '.sorted'
        run_cmd(f"bedtools sort -g {fai} -i {in_bed} > {sort_bed}")
        run_cmd(f"bedtools merge -i {sort_bed} > {out_bed}")
        run_cmd(f"rm -f {sort_bed}")

    prom_m = os.path.join(tmp, 'prom_m.bed')
    gene_m = os.path.join(tmp, 'gene_m.bed')
    down_m = os.path.join(tmp, 'down_m.bed')
    
    sort_merge(prom_safe, prom_m)
    sort_merge(gene_safe, gene_m)
    sort_merge(down_safe, down_m)

    # 5. Compute Distal (whole genome minus the first three regions)
    cat_bed = os.path.join(tmp, 'cat_pgd.bed')
    run_cmd(f"cat {prom_m} {gene_m} {down_m} > {cat_bed}.raw")
    run_cmd(f"bedtools sort -g {fai} -i {cat_bed}.raw > {cat_bed}.sorted")
    run_cmd(f"bedtools merge -i {cat_bed}.sorted > {cat_bed}")
    
    distal_m = os.path.join(tmp, 'distal_m.bed')
    run_cmd(f"bedtools complement -i {cat_bed} -g {fai} > {distal_m}")

    # 6. Safely extract and deduplicate i-motif coordinates
    imotif_dedup = os.path.join(tmp, 'imotif_dedup.bed')
    run_cmd(f"cut -f1,2,3 {imotif} | bedtools sort -g {fai} -i - > {imotif_dedup}.raw")
    run_cmd(f"bedtools merge -i {imotif_dedup}.raw > {imotif_dedup}")

    # 7. Assign by priority
    im_prom = os.path.join(tmp, 'im_prom.bed')
    run_cmd(f"bedtools intersect -a {imotif_dedup} -b {prom_m} -wa -u > {im_prom}")
    
    im_rem1 = os.path.join(tmp, 'im_rem1.bed')
    run_cmd(f"bedtools intersect -a {imotif_dedup} -b {prom_m} -v > {im_rem1}")
    im_gene = os.path.join(tmp, 'im_gene.bed')
    run_cmd(f"bedtools intersect -a {im_rem1} -b {gene_m} -wa -u > {im_gene}")
    
    im_rem2 = os.path.join(tmp, 'im_rem2.bed')
    run_cmd(f"bedtools intersect -a {im_rem1} -b {gene_m} -v > {im_rem2}")
    im_down = os.path.join(tmp, 'im_down.bed')
    run_cmd(f"bedtools intersect -a {im_rem2} -b {down_m} -wa -u > {im_down}")
    
    im_rem3 = os.path.join(tmp, 'im_rem3.bed')
    run_cmd(f"bedtools intersect -a {im_rem2} -b {down_m} -v > {im_rem3}")
    im_distal = im_rem3

    # 8. Statistics
    def count_lines(f): return sum(1 for _ in open(f)) if os.path.exists(f) else 0
    def calc_len(f):
        if not os.path.exists(f): return 0
        return sum(int(l.split('\t')[2]) - int(l.split('\t')[1]) for l in open(f) if l.strip())

    total_im = count_lines(imotif_dedup)
    genome_len = calc_len(genome_bed)
    
    if total_im == 0 or genome_len == 0: return None

    res = {
        'Species': sp_name,
        'Genome_len_Mb': genome_len / 1e6,
        'Total_iM': total_im,
        'Genome_density': total_im / (genome_len / 1e6),
        
        'Prom_len_Mb': calc_len(prom_m) / 1e6,
        'Prom_iM': count_lines(im_prom),
        'Prom_density': count_lines(im_prom) / (calc_len(prom_m) / 1e6) if calc_len(prom_m) > 0 else 0,
        
        'Gene_len_Mb': calc_len(gene_m) / 1e6,
        'Gene_iM': count_lines(im_gene),
        'Gene_density': count_lines(im_gene) / (calc_len(gene_m) / 1e6) if calc_len(gene_m) > 0 else 0,
        
        'Down_len_Mb': calc_len(down_m) / 1e6,
        'Down_iM': count_lines(im_down),
        'Down_density': count_lines(im_down) / (calc_len(down_m) / 1e6) if calc_len(down_m) > 0 else 0,
        
        'Distal_len_Mb': calc_len(distal_m) / 1e6,
        'Distal_iM': count_lines(im_distal),
        'Distal_density': count_lines(im_distal) / (calc_len(distal_m) / 1e6) if calc_len(distal_m) > 0 else 0,
    }
    
    run_cmd(f"rm -rf {tmp}")
    return res

if __name__ == '__main__':
    sp_path = sys.argv[1]
    out_file = sys.argv[2]
    res = process_species(sp_path, os.path.dirname(out_file))
    if res:
        df = pd.DataFrame([res])
        header = not os.path.exists(out_file)
        df.to_csv(out_file, sep='\t', index=False, mode='a', header=header)
