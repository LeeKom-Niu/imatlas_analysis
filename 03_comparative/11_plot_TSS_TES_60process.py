#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
import time
import signal
import psutil

# Optimized parameter settings - keep 1000 bootstrap iterations
window_size = 500
relative_positions = np.arange(-window_size, window_size + 1)
output_dir = "enrichment_results_TSS_TES_1bp_60process_opt"
n_processes = 60
chunksize = 5
bootstrap_n = 1000  # keep 1000 bootstrap iterations
ci_level = 95

# Data directory
base_dir = "/datapool/home/2023200496/niulk/my_project/genome/archaea"
gene_file_suffix = "gene_with_strand.bed"
imotif_file_suffix = "imotif_with_strand.bed"

def init_worker():
    """初始化工作进程"""
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def calculate_enrichment_optimized(args):
    """优化版本 - 减少内存使用"""
    gene_file, imotif_file, genome_name = args
    try:
        # Efficient data types
        dtype = {
            "chrom": "category", 
            "start": "int32", 
            "end": "int32", 
            "strand": "category"
        }
        
        # Read only needed columns
        genes = pd.read_csv(gene_file, sep="\t", header=None, 
                          names=["chrom", "start", "end", "strand"], 
                          usecols=[0,1,2,5], dtype=dtype)
        imotifs = pd.read_csv(imotif_file, sep="\t", header=None,
                            names=["chrom", "start", "end", "strand"],
                            usecols=[0,1,2,5], dtype=dtype)

        # Pre-group by strand
        imotif_plus = imotifs[imotifs["strand"] == "+"]
        imotif_minus = imotifs[imotifs["strand"] == "-"]
        
        # Convert to numpy arrays
        imotif_data = {
            '+': {
                'starts': imotif_plus["start"].values.astype('int32') if len(imotif_plus) > 0 else np.array([], dtype='int32'),
                'ends': imotif_plus["end"].values.astype('int32') if len(imotif_plus) > 0 else np.array([], dtype='int32')
            },
            '-': {
                'starts': imotif_minus["start"].values.astype('int32') if len(imotif_minus) > 0 else np.array([], dtype='int32'),
                'ends': imotif_minus["end"].values.astype('int32') if len(imotif_minus) > 0 else np.array([], dtype='int32')
            }
        }
        
        # Pre-allocate arrays
        results = {
            "TSS": {
                "template": np.zeros(len(relative_positions), dtype=np.int32),
                "non_template": np.zeros(len(relative_positions), dtype=np.int32)
            },
            "TES": {
                "template": np.zeros(len(relative_positions), dtype=np.int32),
                "non_template": np.zeros(len(relative_positions), dtype=np.int32)
            }
        }
        
        # Process each gene
        for _, gene in genes.iterrows():
            strand = gene["strand"]
            tss = gene["start"] if strand == "+" else gene["end"]
            tes = gene["end"] if strand == "+" else gene["start"]
            template_strand = "-" if strand == "+" else "+"
            
            for region, center in [("TSS", tss), ("TES", tes)]:
                # Template strand calculation
                tmpl_data = imotif_data[template_strand]
                if len(tmpl_data['starts']) > 0:
                    for i, rel_pos in enumerate(relative_positions):
                        pos = center + rel_pos
                        count = np.sum((tmpl_data['starts'] <= pos) & (tmpl_data['ends'] >= pos))
                        results[region]["template"][i] += count
                
                # Non-template strand calculation
                nontmpl_data = imotif_data[strand]
                if len(nontmpl_data['starts']) > 0:
                    for i, rel_pos in enumerate(relative_positions):
                        pos = center + rel_pos
                        count = np.sum((nontmpl_data['starts'] <= pos) & (nontmpl_data['ends'] >= pos))
                        results[region]["non_template"][i] += count
        
        # Normalization
        def normalize(arr):
            total = arr.sum()
            return arr / (total/len(relative_positions)) if total > 0 else arr.astype(np.float64)
        
        normalized_results = {
            region: {
                "template": normalize(results[region]["template"]),
                "non_template": normalize(results[region]["non_template"]),
                "template_raw": results[region]["template"].copy(),
                "non_template_raw": results[region]["non_template"].copy()
            } for region in results
        }
        
        return {"genome": genome_name, "data": normalized_results}
    
    except Exception as e:
        print(f"Error processing {genome_name}: {str(e)}")
        return {"genome": genome_name, "error": str(e)}

def bootstrap_confidence_optimized(data, n_bootstrap=1000, ci=95):
    """优化的bootstrap计算，保持1000次但使用矩阵运算加速"""
    n_samples = len(data)
    if n_samples == 0:
        return np.zeros_like(data[0]), np.zeros_like(data[0])
    
    # Pre-allocate arrays
    bootstraps = np.zeros((n_bootstrap, data.shape[1]))
    
    for i in range(n_bootstrap):
        # Vectorized sampling
        sample_indices = np.random.choice(n_samples, size=n_samples, replace=True)
        sample = data[sample_indices]
        bootstraps[i] = np.mean(sample, axis=0)
    
    lower = np.percentile(bootstraps, (100-ci)/2, axis=0)
    upper = np.percentile(bootstraps, 100-(100-ci)/2, axis=0)
    return lower, upper

def plot_enrichment(df, region, name, output_dir, with_ci=True):
    """绘制富集度图表"""
    plt.figure(figsize=(10, 5))
    
    # Template strand
    plt.plot(df["position"], df["template_enrich"], 
            label="Template Strand", color="blue", linewidth=1)
    
    # Non-template strand
    plt.plot(df["position"], df["non_template_enrich"],
            label="Non-template Strand", color="red", linewidth=1)
    
    if with_ci and "template_lower" in df.columns:
        plt.fill_between(df["position"], df["template_lower"], df["template_upper"],
                        color="blue", alpha=0.1, label="Template 95% CI")
        plt.fill_between(df["position"], df["non_template_lower"], df["non_template_upper"],
                        color="red", alpha=0.1, label="Non-template 95% CI")
    
    plt.axvline(x=0, color="black", linestyle="--", linewidth=0.5, label=region)
    title = f"i-Motif Enrichment around {region}"
    if name:
        title = f"{name}\n{title}"
    plt.title(title, fontsize=10)
    plt.xlabel(f"Relative Position to {region} (bp)")
    plt.ylabel("Normalized Enrichment")
    plt.legend(fontsize=8)
    plt.grid(alpha=0.1)
    
    y_max = max(df["template_enrich"].max(), df["non_template_enrich"].max()) * 1.1
    plt.ylim(0, y_max if y_max > 0 else 1)
    
    plt.tight_layout()
    ci_suffix = "_with_ci" if with_ci else ""
    plt.savefig(os.path.join(output_dir, f"{region}_enrichment{ci_suffix}.png"), 
                dpi=150, bbox_inches='tight')
    plt.close()

def process_genomes_parallel(base_dir):
    """并行处理流程"""
    os.makedirs(output_dir, exist_ok=True)
    print("Scanning genome directories...")
    subdirs = glob(os.path.join(base_dir, '*/output_genomic'))
    
    tasks = []
    for subdir in subdirs:
        genome_name = os.path.basename(os.path.dirname(subdir))
        gene_file = os.path.join(subdir, gene_file_suffix)
        imotif_file = os.path.join(subdir, imotif_file_suffix)
        if os.path.isfile(gene_file) and os.path.isfile(imotif_file):
            tasks.append((gene_file, imotif_file, genome_name))
    
    print(f"Found {len(tasks)} valid genomes")
    print(f"Using {n_processes} processes, bootstrap_n={bootstrap_n}")
    
    # Monitor initial resources
    initial_resources = monitor_system_resources()
    print(f"Initial memory: {initial_resources['memory_used_gb']:.1f}GB used, "
          f"{initial_resources['memory_available_gb']:.1f}GB available")
    
    # Parallel processing
    results = []
    successful = 0
    failed = 0
    
    with Pool(n_processes, initializer=init_worker) as pool:
        with tqdm(total=len(tasks), desc="Processing genomes") as pbar:
            for i, result in enumerate(pool.imap_unordered(calculate_enrichment_optimized, tasks, chunksize=chunksize)):
                if "error" not in result:
                    results.append(result)
                    successful += 1
                    save_single_genome_result(result)
                else:
                    failed += 1
                    print(f"Failed: {result['genome']} - {result['error']}")
                
                pbar.update()
                pbar.set_postfix(success=successful, failed=failed, current=result["genome"][:15])
                
                # Monitor resources every 50 genomes
                if i % 50 == 0 and i > 0:
                    resources = monitor_system_resources()
                    print(f"Progress checkpoint - Processed {i}/{len(tasks)}, Memory: {resources['memory_percent']:.1f}% used")
    
    print(f"Processing completed: {successful} successful, {failed} failed")
    return results

def save_single_genome_result(result):
    """保存单个基因组结果"""
    genome_name = result["genome"]
    genome_dir = os.path.join(output_dir, "genome_results", genome_name)
    os.makedirs(genome_dir, exist_ok=True)
    
    for region in ["TSS", "TES"]:
        df = pd.DataFrame({
            "position": relative_positions,
            "template_enrich": result["data"][region]["template"],
            "non_template_enrich": result["data"][region]["non_template"],
            "template_raw": result["data"][region]["template_raw"],
            "non_template_raw": result["data"][region]["non_template_raw"]
        })
        df.to_csv(os.path.join(genome_dir, f"{region}_results.tsv"), sep="\t", index=False)

def monitor_system_resources():
    """监控系统资源使用情况"""
    memory = psutil.virtual_memory()
    return {
        'memory_used_gb': memory.used / (1024**3),
        'memory_available_gb': memory.available / (1024**3),
        'memory_percent': memory.percent
    }

def save_final_results_with_ci(results):
    """保存最终结果，使用1000次bootstrap计算置信区间"""
    if not results:
        print("No valid results to process!")
        return {}
        
    regions = ["TSS", "TES"]
    merged = {
        region: {
            "template": np.array([r["data"][region]["template"] for r in results]),
            "non_template": np.array([r["data"][region]["non_template"] for r in results])
        } for region in regions
    }
    
    def final_enrich(arr):
        total = arr.sum(axis=0)
        return total / (total.sum()/len(relative_positions)) if total.sum() > 0 else total
    
    for region in regions:
        print(f"Processing {region} with bootstrap_n={bootstrap_n}...")
        template_mean = final_enrich(merged[region]["template"])
        nontmpl_mean = final_enrich(merged[region]["non_template"])
        
        # Compute confidence intervals with 1000 bootstrap iterations
        template_lower, template_upper = bootstrap_confidence_optimized(merged[region]["template"], n_bootstrap=1000)
        nontmpl_lower, nontmpl_upper = bootstrap_confidence_optimized(merged[region]["non_template"], n_bootstrap=1000)
        
        # Save results
        df = pd.DataFrame({
            "position": relative_positions,
            "template_enrich": template_mean,
            "non_template_enrich": nontmpl_mean,
            "template_lower": template_lower,
            "template_upper": template_upper,
            "non_template_lower": nontmpl_lower,
            "non_template_upper": nontmpl_upper,
            "template_raw": np.sum(merged[region]["template"], axis=0),
            "non_template_raw": np.sum(merged[region]["non_template"], axis=0)
        })
        df.to_csv(os.path.join(output_dir, f"{region}_combined_results.tsv"), sep="\t", index=False)
        
        # Plotting
        plot_enrichment(df, region, "All Genomes Combined", output_dir, with_ci=True)
    
    # Summary statistics
    stats = {
        region: {
            "template": sum(np.any(merged[region]["template"] > 0, axis=1)),
            "non_template": sum(np.any(merged[region]["non_template"] > 0, axis=1))
        } for region in regions
    }
    
    return stats

if __name__ == "__main__":
    start_time = time.time()
    
    print("Starting OPTIMIZED i-motif enrichment analysis...")
    print(f"Output directory: {output_dir}")
    print(f"Window size: ±{window_size}bp, Resolution: 1bp")
    print(f"Parallel processes: {n_processes}")
    print(f"Bootstrap iterations: {bootstrap_n} (保持1000次确保统计可靠性)")
    
    # Process data
    results = process_genomes_parallel(base_dir)
    
    if results:
        stats = save_final_results_with_ci(results)
        
        total_time = time.time() - start_time
        print(f"\n=== Analysis Completed ===")
        print(f"Total time: {total_time:.2f} seconds ({total_time/60:.2f} minutes)")
        print(f"Successful genomes: {len(results)}")
        
        for region in ["TSS", "TES"]:
            print(f"\n{region} Statistics:")
            print(f"  Genomes with template strand data: {stats[region]['template']}")
            print(f"  Genomes with non-template strand data: {stats[region]['non_template']}")
        
        print(f"\nResults saved to: {os.path.abspath(output_dir)}")
        
        # Final resource usage report
        final_resources = monitor_system_resources()
        print(f"Final memory usage: {final_resources['memory_percent']:.1f}%")
    else:
        print("No valid results were generated!")
