#!/usr/bin/env python3
"""
Wisecondor/WisecondorX Quality Check Tool

Analyzes Wisecondor/WisecondorX NPZ files (and optionally BAM) to determine
analysis quality and feasibility. Considers:
- Read count per bin
- Uniformity (coefficient of variation)
- Density (zero bin fraction)
- Noise level (from z-scores)

Can work with:
1. BAM file only (basic analysis)
2. Wisecondor NPZ file (converted BAM)
3. WisecondorX NPZ file (with z-scores from prediction)

Usage:
    # With BAM
    python3 wisecondor_quality_check.py --bam sample.bam --bed targets.bed
    
    # With Wisecondor NPZ (converted)
    python3 wisecondor_quality_check.py --wc-npz sample.wc.npz --bed targets.bed
    
    # With WisecondorX NPZ (predicted)
    python3 wisecondor_quality_check.py --wcx-npz sample.wcx.npz --bed targets.bed
    
    # With all inputs (most comprehensive)
    python3 wisecondor_quality_check.py --bam sample.bam --wc-npz sample.wc.npz --wcx-npz sample.wcx.npz --bed targets.bed

Author: Ken
Version: 1.0
"""

import argparse
import sys
import os
from pathlib import Path
from collections import defaultdict
import warnings

try:
    import pysam
    import numpy as np
    import pandas as pd
    from scipy import stats
except ImportError as e:
    print(f"ERROR: Required library not found: {e}", file=sys.stderr)
    print("Please install: pip install pysam numpy pandas scipy", file=sys.stderr)
    sys.exit(1)


class WisecondorQualityAnalyzer:
    """
    Comprehensive quality analyzer for Wisecondor/WisecondorX analysis
    """
    
    def __init__(self, bin_size=200000):
        """
        Initialize quality analyzer
        
        Args:
            bin_size: Bin size used in Wisecondor analysis
        """
        self.bin_size = bin_size
        self.bam_data = None
        self.wc_data = None
        self.wcx_data = None
    
    def load_bam(self, bam_path):
        """Load BAM file for analysis"""
        try:
            self.bamfile = pysam.AlignmentFile(bam_path, "rb")
            if not self.bamfile.has_index():
                print(f"WARNING: BAM file is not indexed. Creating index...", file=sys.stderr)
                pysam.index(bam_path)
                self.bamfile = pysam.AlignmentFile(bam_path, "rb")
            
            # Detect chromosome naming convention (chr1 vs 1)
            self.bam_chromosomes = self.bamfile.references
            self.has_chr_prefix = any(chrom.startswith('chr') for chrom in self.bam_chromosomes[:5])
            
            print(f"✓ Loaded BAM: {bam_path}")
            print(f"  Chromosome naming: {'chr1, chr2...' if self.has_chr_prefix else '1, 2...'}")
            return True
        except Exception as e:
            print(f"ERROR: Cannot open BAM file: {e}", file=sys.stderr)
            return False
    
    def load_wisecondor_npz(self, npz_path):
        """Load Wisecondor converted NPZ file"""
        try:
            data = np.load(npz_path, allow_pickle=True)
            self.wc_data = {
                'chromosomes': data.get('chromosomes', None),
                'starts': data.get('starts', None),
                'sizes': data.get('sizes', None),
                'results': data.get('results', None),
            }
            print(f"✓ Loaded Wisecondor NPZ: {npz_path}")
            print(f"  Chromosomes: {len(self.wc_data['chromosomes']) if self.wc_data['chromosomes'] is not None else 'N/A'}")
            if self.wc_data['results'] is not None:
                print(f"  Total bins: {len(self.wc_data['results'])}")
            return True
        except Exception as e:
            print(f"ERROR: Cannot load Wisecondor NPZ: {e}", file=sys.stderr)
            return False
    
    def load_wisecondorx_npz(self, npz_path):
        """Load WisecondorX converted or predicted NPZ file"""
        try:
            data = np.load(npz_path, allow_pickle=True)
            
            # WCX NPZ can have different structures depending on convert vs predict
            self.wcx_data = {}
            
            # Try to load all possible keys
            for key in data.files:
                self.wcx_data[key] = data[key]
            
            print(f"✓ Loaded WisecondorX NPZ: {npz_path}")
            print(f"  Available data: {', '.join(self.wcx_data.keys())}")
            
            # Check if this is a predicted file (has z-scores)
            if 'z_scores' in self.wcx_data or 'zscore' in self.wcx_data:
                print(f"  Type: Predicted (with z-scores)")
            else:
                print(f"  Type: Converted (read counts only)")
            
            return True
        except Exception as e:
            print(f"ERROR: Cannot load WisecondorX NPZ: {e}", file=sys.stderr)
            return False
    
    def normalize_chromosome(self, chrom):
        """
        Normalize chromosome name to match BAM file format
        
        Converts between '1' <-> 'chr1' as needed
        """
        if not hasattr(self, 'has_chr_prefix'):
            # If BAM not loaded, return as-is
            return chrom
        
        # Remove 'chr' prefix if present
        chrom_num = chrom.replace('chr', '')
        
        # Add 'chr' prefix if BAM uses it and chrom doesn't have it
        if self.has_chr_prefix and not chrom.startswith('chr'):
            return f"chr{chrom_num}"
        # Remove 'chr' prefix if BAM doesn't use it but chrom has it
        elif not self.has_chr_prefix and chrom.startswith('chr'):
            return chrom_num
        else:
            return chrom
    
    def read_bed_file(self, bed_path):
        """Read BED file and return regions"""
        regions = []
        try:
            with open(bed_path, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith('#') or line.startswith('track'):
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) < 3:
                        print(f"WARNING: Skipping invalid BED line {line_num}", file=sys.stderr)
                        continue
                    
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    name = parts[3] if len(parts) > 3 else f"{chrom}:{start}-{end}"
                    
                    # Normalize chromosome name to match BAM format
                    chrom_normalized = self.normalize_chromosome(chrom)
                    
                    regions.append((chrom_normalized, start, end, name))
            
            print(f"✓ Loaded {len(regions)} regions from BED file")
            return regions
        except Exception as e:
            print(f"ERROR: Cannot read BED file: {e}", file=sys.stderr)
            sys.exit(1)
    
    def get_bins_in_region(self, chrom, start, end):
        """Get bin indices that overlap with a region"""
        if self.wc_data is None and self.wcx_data is None:
            return None
        
        # Try WC data first
        data = self.wc_data if self.wc_data else self.wcx_data
        
        if data.get('chromosomes') is None or data.get('starts') is None:
            return None
        
        chromosomes = data['chromosomes']
        starts = data['starts']
        
        # Find bins in this region
        bin_indices = []
        for i, (bin_chrom, bin_start) in enumerate(zip(chromosomes, starts)):
            bin_end = bin_start + self.bin_size
            
            # Check if bin overlaps with region
            if bin_chrom == chrom and not (bin_end <= start or bin_start >= end):
                bin_indices.append(i)
        
        return bin_indices
    
    def analyze_region_from_bam(self, chrom, start, end):
        """Analyze region from BAM file"""
        if not hasattr(self, 'bamfile'):
            return None
        
        num_bins = (end - start + self.bin_size - 1) // self.bin_size
        bin_counts = np.zeros(num_bins, dtype=np.int32)
        
        for i in range(num_bins):
            bin_start = start + (i * self.bin_size)
            bin_end = min(start + ((i + 1) * self.bin_size), end)
            count = self.bamfile.count(chrom, bin_start, bin_end)
            bin_counts[i] = count
        
        return self._calculate_quality_metrics(bin_counts, "BAM")
    
    def analyze_region_from_npz(self, chrom, start, end, source="WC"):
        """Analyze region from NPZ file"""
        data = self.wc_data if source == "WC" else self.wcx_data
        
        if data is None:
            return None
        
        bin_indices = self.get_bins_in_region(chrom, start, end)
        if bin_indices is None or len(bin_indices) == 0:
            return None
        
        # Get read counts or ratios
        if source == "WC" and 'results' in data and data['results'] is not None:
            bin_values = data['results'][bin_indices]
        elif 'bins' in data:
            bin_values = data['bins'][bin_indices]
        elif 'reads' in data:
            bin_values = data['reads'][bin_indices]
        else:
            return None
        
        metrics = self._calculate_quality_metrics(bin_values, source)
        
        # Add z-score analysis if available (WCX predicted)
        if source == "WCX":
            if 'z_scores' in data or 'zscore' in data:
                z_key = 'z_scores' if 'z_scores' in data else 'zscore'
                z_scores = data[z_key]
                if z_scores is not None and len(z_scores) > max(bin_indices):
                    region_z_scores = z_scores[bin_indices]
                    metrics['z_score_metrics'] = self._analyze_z_scores(region_z_scores)
        
        return metrics
    
    def _calculate_quality_metrics(self, values, source=""):
        """Calculate quality metrics from bin values"""
        if len(values) == 0:
            return None
        
        # Convert to numpy array
        values = np.array(values, dtype=np.float64)
        
        # Basic statistics
        total = np.sum(values)
        mean = np.mean(values)
        median = np.median(values)
        std = np.std(values)
        min_val = np.min(values)
        max_val = np.max(values)
        
        # Uniformity: Coefficient of Variation (CV)
        cv = std / mean if mean > 0 else np.inf
        
        # Density: fraction of zero bins
        zero_bins = np.sum(values == 0)
        zero_fraction = zero_bins / len(values)
        non_zero_fraction = 1 - zero_fraction
        
        # Low count bins (< 10 reads)
        low_count_bins = np.sum(values < 10)
        low_count_fraction = low_count_bins / len(values)
        
        # Quartiles
        q25, q75 = np.percentile(values, [25, 75])
        iqr = q75 - q25
        
        # MAD (Median Absolute Deviation) - robust measure of variability
        mad = np.median(np.abs(values - median))
        
        return {
            'source': source,
            'num_bins': len(values),
            'total': total,
            'mean': mean,
            'median': median,
            'std': std,
            'min': min_val,
            'max': max_val,
            'cv': cv,  # Coefficient of Variation
            'zero_bins': zero_bins,
            'zero_fraction': zero_fraction,
            'non_zero_fraction': non_zero_fraction,
            'low_count_bins': low_count_bins,
            'low_count_fraction': low_count_fraction,
            'q25': q25,
            'q75': q75,
            'iqr': iqr,
            'mad': mad,
        }
    
    def _analyze_z_scores(self, z_scores):
        """Analyze z-score distribution"""
        if len(z_scores) == 0:
            return None
        
        z_scores = np.array(z_scores, dtype=np.float64)
        
        # Remove NaN/Inf
        z_scores = z_scores[np.isfinite(z_scores)]
        
        if len(z_scores) == 0:
            return None
        
        return {
            'mean_zscore': np.mean(z_scores),
            'median_zscore': np.median(z_scores),
            'std_zscore': np.std(z_scores),
            'mad_zscore': np.median(np.abs(z_scores - np.median(z_scores))),
            'max_abs_zscore': np.max(np.abs(z_scores)),
            'num_high_zscore': np.sum(np.abs(z_scores) > 3),  # |z| > 3
            'fraction_high_zscore': np.sum(np.abs(z_scores) > 3) / len(z_scores),
        }
    
    def assess_quality(self, metrics):
        """
        Assess overall quality and determine if analysis is feasible
        
        Quality criteria:
        1. Sufficient read count (mean >= 10 reads/bin for 0.2x)
        2. Good uniformity (CV < 1.0, preferably < 0.5)
        3. Low zero fraction (< 20%)
        4. Low noise (MAD z-score < 1.0 if available)
        """
        if metrics is None:
            return {
                'is_analyzable': False,
                'quality_score': 0,
                'quality_grade': 'F',
                'reasons': ['No data available']
            }
        
        reasons = []
        quality_score = 100
        
        # Criterion 1: Read count
        mean_reads = metrics.get('mean', 0)
        if mean_reads < 5:
            quality_score -= 40
            reasons.append(f"Very low read count (mean={mean_reads:.1f} < 5)")
        elif mean_reads < 10:
            quality_score -= 20
            reasons.append(f"Low read count (mean={mean_reads:.1f} < 10)")
        elif mean_reads < 20:
            quality_score -= 5
        
        # Criterion 2: Uniformity (CV)
        cv = metrics.get('cv', np.inf)
        if cv > 2.0:
            quality_score -= 30
            reasons.append(f"Very poor uniformity (CV={cv:.2f} > 2.0)")
        elif cv > 1.0:
            quality_score -= 15
            reasons.append(f"Poor uniformity (CV={cv:.2f} > 1.0)")
        elif cv > 0.5:
            quality_score -= 5
        
        # Criterion 3: Density
        zero_fraction = metrics.get('zero_fraction', 1.0)
        if zero_fraction > 0.5:
            quality_score -= 30
            reasons.append(f"Too many zero bins ({zero_fraction:.1%} > 50%)")
        elif zero_fraction > 0.2:
            quality_score -= 15
            reasons.append(f"Many zero bins ({zero_fraction:.1%} > 20%)")
        
        # Criterion 4: Z-score noise (if available)
        if 'z_score_metrics' in metrics:
            z_metrics = metrics['z_score_metrics']
            mad_z = z_metrics.get('mad_zscore', 0)
            
            if mad_z > 2.0:
                quality_score -= 20
                reasons.append(f"High noise (MAD z-score={mad_z:.2f} > 2.0)")
            elif mad_z > 1.0:
                quality_score -= 10
                reasons.append(f"Moderate noise (MAD z-score={mad_z:.2f} > 1.0)")
        
        # Determine grade
        if quality_score >= 90:
            quality_grade = 'A'
        elif quality_score >= 80:
            quality_grade = 'B'
        elif quality_score >= 70:
            quality_grade = 'C'
        elif quality_score >= 60:
            quality_grade = 'D'
        else:
            quality_grade = 'F'
        
        # Analyzable if score >= 60
        is_analyzable = quality_score >= 60
        
        if not reasons:
            reasons = ['Good quality']
        
        return {
            'is_analyzable': is_analyzable,
            'quality_score': quality_score,
            'quality_grade': quality_grade,
            'reasons': reasons
        }
    
    def analyze_region(self, chrom, start, end, name):
        """Comprehensive region analysis using all available data sources"""
        results = {
            'chrom': chrom,
            'start': start,
            'end': end,
            'name': name,
            'size': end - start,
        }
        
        # Analyze from BAM
        if hasattr(self, 'bamfile'):
            bam_metrics = self.analyze_region_from_bam(chrom, start, end)
            if bam_metrics:
                results['bam_metrics'] = bam_metrics
                results['bam_quality'] = self.assess_quality(bam_metrics)
        
        # Analyze from Wisecondor NPZ
        if self.wc_data is not None:
            wc_metrics = self.analyze_region_from_npz(chrom, start, end, "WC")
            if wc_metrics:
                results['wc_metrics'] = wc_metrics
                results['wc_quality'] = self.assess_quality(wc_metrics)
        
        # Analyze from WisecondorX NPZ
        if self.wcx_data is not None:
            wcx_metrics = self.analyze_region_from_npz(chrom, start, end, "WCX")
            if wcx_metrics:
                results['wcx_metrics'] = wcx_metrics
                results['wcx_quality'] = self.assess_quality(wcx_metrics)
        
        # Overall assessment (use best available)
        best_quality = None
        best_source = None
        
        for source in ['wcx_quality', 'wc_quality', 'bam_quality']:
            if source in results:
                quality = results[source]
                if best_quality is None or quality['quality_score'] > best_quality['quality_score']:
                    best_quality = quality
                    best_source = source.replace('_quality', '').upper()
        
        results['overall_quality'] = best_quality
        results['best_source'] = best_source
        
        return results
    
    def analyze_all_regions(self, regions):
        """Analyze all regions"""
        results = []
        
        print(f"\nAnalyzing {len(regions)} regions...")
        for i, (chrom, start, end, name) in enumerate(regions, 1):
            if i % 10 == 0 or i == len(regions):
                print(f"  Progress: {i}/{len(regions)} regions", end='\r')
            
            result = self.analyze_region(chrom, start, end, name)
            results.append(result)
        
        print(f"  Progress: {len(regions)}/{len(regions)} regions - Complete!")
        return results
    
    def generate_tsv_summary(self, results, output_path):
        """Generate TSV summary table"""
        try:
            with open(output_path, 'w') as f:
                # Header
                header_cols = [
                    'Index',
                    'Region_Name',
                    'Chromosome',
                    'Start',
                    'End',
                    'Size_bp',
                    'Num_Bins',
                    # Overall
                    'Best_Source',
                    'Quality_Grade',
                    'Quality_Score',
                    'Is_Analyzable',
                    # BAM metrics
                    'BAM_Analyzable',
                    'BAM_Grade',
                    'BAM_Score',
                    # Wisecondor metrics
                    'WC_Analyzable',
                    'WC_Grade',
                    'WC_Score',
                    # WisecondorX metrics
                    'WCX_Analyzable',
                    'WCX_Grade',
                    'WCX_Score',
                    # Read statistics
                    'Mean_Reads',
                    'Median_Reads',
                    'CV',
                    'Zero_Fraction',
                    'MAD',
                    'Z_Score_MAD',
                    'High_Zscore_Fraction',
                    'Issues'
                ]
                f.write('\t'.join(header_cols) + '\n')
                
                # Data rows
                for idx, result in enumerate(results, 1):
                    overall = result.get('overall_quality', {})
                    best_source = result.get('best_source', 'N/A')
                    
                    # Get BAM quality
                    bam_quality = result.get('bam_quality', {})
                    bam_analyzable = 'Yes' if bam_quality.get('is_analyzable', False) else 'No' if bam_quality else 'N/A'
                    bam_grade = bam_quality.get('quality_grade', 'N/A') if bam_quality else 'N/A'
                    bam_score = str(bam_quality.get('quality_score', 0)) if bam_quality else 'N/A'
                    
                    # Get Wisecondor quality
                    wc_quality = result.get('wc_quality', {})
                    wc_analyzable = 'Yes' if wc_quality.get('is_analyzable', False) else 'No' if wc_quality else 'N/A'
                    wc_grade = wc_quality.get('quality_grade', 'N/A') if wc_quality else 'N/A'
                    wc_score = str(wc_quality.get('quality_score', 0)) if wc_quality else 'N/A'
                    
                    # Get WisecondorX quality
                    wcx_quality = result.get('wcx_quality', {})
                    wcx_analyzable = 'Yes' if wcx_quality.get('is_analyzable', False) else 'No' if wcx_quality else 'N/A'
                    wcx_grade = wcx_quality.get('quality_grade', 'N/A') if wcx_quality else 'N/A'
                    wcx_score = str(wcx_quality.get('quality_score', 0)) if wcx_quality else 'N/A'
                    
                    # Get best metrics for read statistics
                    metrics = None
                    if best_source == 'BAM' and 'bam_metrics' in result:
                        metrics = result['bam_metrics']
                    elif best_source == 'WC' and 'wc_metrics' in result:
                        metrics = result['wc_metrics']
                    elif best_source == 'WCX' and 'wcx_metrics' in result:
                        metrics = result['wcx_metrics']
                    
                    # Extract values
                    num_bins = metrics.get('num_bins', 0) if metrics else 0
                    mean_reads = f"{metrics.get('mean', 0):.2f}" if metrics else '0.00'
                    median_reads = f"{metrics.get('median', 0):.2f}" if metrics else '0.00'
                    cv = f"{metrics.get('cv', 0):.3f}" if metrics else '0.000'
                    zero_frac = f"{metrics.get('zero_fraction', 0):.3f}" if metrics else '0.000'
                    mad = f"{metrics.get('mad', 0):.2f}" if metrics else '0.00'
                    
                    # Z-score metrics (primarily from WCX)
                    z_mad = 'N/A'
                    high_z_frac = 'N/A'
                    if metrics and 'z_score_metrics' in metrics:
                        z_m = metrics['z_score_metrics']
                        z_mad = f"{z_m.get('mad_zscore', 0):.3f}"
                        high_z_frac = f"{z_m.get('fraction_high_zscore', 0):.3f}"
                    
                    # Issues
                    issues = '; '.join(overall.get('reasons', ['None']))
                    
                    row_cols = [
                        str(idx),
                        result['name'],
                        result['chrom'],
                        str(result['start']),
                        str(result['end']),
                        str(result['size']),
                        str(num_bins),
                        # Overall
                        best_source,
                        overall.get('quality_grade', 'F'),
                        str(overall.get('quality_score', 0)),
                        'Yes' if overall.get('is_analyzable', False) else 'No',
                        # BAM
                        bam_analyzable,
                        bam_grade,
                        bam_score,
                        # Wisecondor
                        wc_analyzable,
                        wc_grade,
                        wc_score,
                        # WisecondorX
                        wcx_analyzable,
                        wcx_grade,
                        wcx_score,
                        # Read statistics
                        mean_reads,
                        median_reads,
                        cv,
                        zero_frac,
                        mad,
                        z_mad,
                        high_z_frac,
                        issues
                    ]
                    f.write('\t'.join(row_cols) + '\n')
            
            print(f"✓ TSV summary written to {output_path}")
            
        except Exception as e:
            print(f"ERROR: Cannot generate TSV summary: {e}", file=sys.stderr)
    
    def generate_report(self, results, output_path):
        """Generate comprehensive quality report"""
        try:
            with open(output_path, 'w') as f:
                f.write("=" * 80 + "\n")
                f.write("Wisecondor/WisecondorX Quality Assessment Report\n")
                f.write("=" * 80 + "\n\n")
                
                # Data sources
                f.write("Data Sources:\n")
                if hasattr(self, 'bamfile'):
                    f.write("  ✓ BAM file\n")
                if self.wc_data:
                    f.write("  ✓ Wisecondor NPZ\n")
                if self.wcx_data:
                    f.write("  ✓ WisecondorX NPZ\n")
                    if 'z_scores' in self.wcx_data or 'zscore' in self.wcx_data:
                        f.write("    (with z-scores - predicted)\n")
                f.write("\n")
                
                # Summary
                total_regions = len(results)
                analyzable = sum(1 for r in results if r.get('overall_quality', {}).get('is_analyzable', False))
                
                grade_counts = defaultdict(int)
                for r in results:
                    grade = r.get('overall_quality', {}).get('quality_grade', 'F')
                    grade_counts[grade] += 1
                
                f.write(f"Summary:\n")
                f.write(f"  Total regions: {total_regions}\n")
                f.write(f"  Analyzable: {analyzable} ({analyzable/total_regions*100:.1f}%)\n")
                f.write(f"  Not analyzable: {total_regions-analyzable} ({(total_regions-analyzable)/total_regions*100:.1f}%)\n\n")
                
                f.write(f"Quality Grade Distribution:\n")
                for grade in ['A', 'B', 'C', 'D', 'F']:
                    count = grade_counts[grade]
                    if count > 0:
                        f.write(f"  Grade {grade}: {count} ({count/total_regions*100:.1f}%)\n")
                f.write("\n")
                
                # Per-region details
                f.write("=" * 80 + "\n")
                f.write("Per-Region Quality Assessment\n")
                f.write("=" * 80 + "\n\n")
                
                for result in results:
                    overall = result.get('overall_quality', {})
                    is_ok = overall.get('is_analyzable', False)
                    grade = overall.get('quality_grade', 'F')
                    score = overall.get('quality_score', 0)
                    
                    status = "✓ PASS" if is_ok else "✗ FAIL"
                    f.write(f"{status} [{grade}] {result['name']} (score: {score}/100)\n")
                    f.write(f"  Location: {result['chrom']}:{result['start']}-{result['end']} ({result['size']:,} bp)\n")
                    f.write(f"  Best data source: {result.get('best_source', 'N/A')}\n")
                    
                    # Show metrics from each source
                    for source_key in ['bam_metrics', 'wc_metrics', 'wcx_metrics']:
                        if source_key in result:
                            metrics = result[source_key]
                            source_name = metrics['source']
                            quality = result[source_key.replace('_metrics', '_quality')]
                            
                            f.write(f"\n  {source_name} Analysis:\n")
                            f.write(f"    Quality: {quality['quality_grade']} (score: {quality['quality_score']}/100)\n")
                            f.write(f"    Bins: {metrics['num_bins']}\n")
                            f.write(f"    Mean reads/bin: {metrics['mean']:.2f}\n")
                            f.write(f"    Median reads/bin: {metrics['median']:.2f}\n")
                            f.write(f"    CV (uniformity): {metrics['cv']:.3f}\n")
                            f.write(f"    Zero bins: {metrics['zero_bins']}/{metrics['num_bins']} ({metrics['zero_fraction']:.1%})\n")
                            f.write(f"    MAD: {metrics['mad']:.2f}\n")
                            
                            if 'z_score_metrics' in metrics:
                                z_m = metrics['z_score_metrics']
                                f.write(f"    Z-score MAD: {z_m['mad_zscore']:.3f}\n")
                                f.write(f"    High |z| bins (>3): {z_m['num_high_zscore']}/{metrics['num_bins']} ({z_m['fraction_high_zscore']:.1%})\n")
                    
                    f.write(f"\n  Issues: {', '.join(overall.get('reasons', ['None']))}\n")
                    f.write("\n")
            
            print(f"✓ Quality report written to {output_path}")
            
        except Exception as e:
            print(f"ERROR: Cannot generate report: {e}", file=sys.stderr)
    
    def close(self):
        """Close open files"""
        if hasattr(self, 'bamfile'):
            self.bamfile.close()


def main():
    parser = argparse.ArgumentParser(
        description='Wisecondor/WisecondorX Quality Check Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic: BAM only
  python3 wisecondor_quality_check.py --bam sample.bam --bed targets.bed -o quality_report.txt
  
  # With Wisecondor NPZ
  python3 wisecondor_quality_check.py --wc-npz sample.wc.npz --bed targets.bed -o report.txt
  
  # With WisecondorX NPZ (most informative if predicted)
  python3 wisecondor_quality_check.py --wcx-npz sample.wcx.npz --bed targets.bed -o report.txt
  
  # All sources (comprehensive)
  python3 wisecondor_quality_check.py --bam sample.bam --wc-npz sample.wc.npz \\
      --wcx-npz sample.wcx.npz --bed targets.bed -o report.txt
        """
    )
    
    # Input files
    parser.add_argument('--bam', help='Input BAM file')
    parser.add_argument('--wc-npz', help='Wisecondor NPZ file (converted)')
    parser.add_argument('--wcx-npz', help='WisecondorX NPZ file (converted or predicted)')
    parser.add_argument('--bed', required=True, help='BED file with target regions')
    
    # Output
    parser.add_argument('--output', '-o', required=True, help='Output quality report file')
    
    # Parameters
    parser.add_argument('--bin-size', type=int, default=200000, help='Bin size (default: 200000)')
    
    args = parser.parse_args()
    
    # Validate: at least one input source
    if not any([args.bam, args.wc_npz, args.wcx_npz]):
        print("ERROR: At least one input source required (--bam, --wc-npz, or --wcx-npz)", file=sys.stderr)
        sys.exit(1)
    
    # Validate files exist
    for path, name in [(args.bam, 'BAM'), (args.wc_npz, 'WC NPZ'), (args.wcx_npz, 'WCX NPZ'), (args.bed, 'BED')]:
        if path and not os.path.exists(path):
            print(f"ERROR: {name} file not found: {path}", file=sys.stderr)
            sys.exit(1)
    
    print("=" * 80)
    print("Wisecondor/WisecondorX Quality Check Tool")
    print("=" * 80)
    
    # Initialize analyzer
    analyzer = WisecondorQualityAnalyzer(bin_size=args.bin_size)
    
    # Load data sources
    if args.bam:
        analyzer.load_bam(args.bam)
    
    if args.wc_npz:
        analyzer.load_wisecondor_npz(args.wc_npz)
    
    if args.wcx_npz:
        analyzer.load_wisecondorx_npz(args.wcx_npz)
    
    # Load regions
    regions = analyzer.read_bed_file(args.bed)
    
    # Analyze
    results = analyzer.analyze_all_regions(regions)
    
    # Generate detailed report
    analyzer.generate_report(results, args.output)
    
    # Generate TSV summary
    tsv_output = args.output.replace('.txt', '_summary.tsv')
    if tsv_output == args.output:  # If no .txt extension
        tsv_output = args.output + '_summary.tsv'
    analyzer.generate_tsv_summary(results, tsv_output)
    
    # Summary
    print("\n" + "=" * 80)
    print("Analysis Complete!")
    print("=" * 80)
    
    total = len(results)
    analyzable = sum(1 for r in results if r.get('overall_quality', {}).get('is_analyzable', False))
    
    print(f"Total regions: {total}")
    print(f"  ✓ Analyzable: {analyzable} ({analyzable/total*100:.1f}%)")
    print(f"  ✗ Not analyzable: {total-analyzable} ({(total-analyzable)/total*100:.1f}%)")
    
    # Grade distribution
    grade_counts = defaultdict(int)
    for r in results:
        grade = r.get('overall_quality', {}).get('quality_grade', 'F')
        grade_counts[grade] += 1
    
    print(f"\nQuality Grades:")
    for grade in ['A', 'B', 'C', 'D', 'F']:
        count = grade_counts[grade]
        if count > 0:
            print(f"  {grade}: {count} ({count/total*100:.1f}%)")
    
    print(f"\nOutput files:")
    print(f"  Detailed report: {args.output}")
    print(f"  TSV summary:     {tsv_output}")
    
    analyzer.close()


if __name__ == '__main__':
    main()

