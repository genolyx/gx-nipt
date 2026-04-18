#!/usr/bin/env python3
"""
MD Coverage Check Tool

Analyzes BAM coverage over BED regions to determine if Wisecondor/WisecondorX analysis is possible.
For regions with insufficient coverage, identifies and outputs the coverable portions.

Usage:
    python3 md_coverage_check.py --bam proper_paired.bam --bed target_regions.bed --output coverable_regions.bed

Options:
    --min-coverage INT    Minimum read depth per base (default: 1)
    --min-region-size INT Minimum region size in bp (default: 1000)
    --coverage-threshold FLOAT Minimum fraction of bases covered (0.0-1.0, default: 0.9)
    --bin-size INT        Bin size for Wisecondor analysis (default: 200000)
    --report              Generate detailed coverage report

Author: Ken
Version: 1.0
"""

import argparse
import sys
import os
from pathlib import Path
from collections import defaultdict

try:
    import pysam
    import numpy as np
    import pandas as pd
except ImportError as e:
    print(f"ERROR: Required library not found: {e}", file=sys.stderr)
    print("Please install required packages: pip install pysam numpy pandas", file=sys.stderr)
    sys.exit(1)


class CoverageAnalyzer:
    """Analyzes BAM coverage for MD detection"""
    
    def __init__(self, bam_path, min_coverage=1, min_region_size=1000, 
                 coverage_threshold=0.9, bin_size=200000):
        """
        Initialize coverage analyzer
        
        Args:
            bam_path: Path to BAM file
            min_coverage: Minimum read depth per base
            min_region_size: Minimum region size in bp
            coverage_threshold: Minimum fraction of bases that must be covered
            bin_size: Bin size for Wisecondor analysis
        """
        self.bam_path = bam_path
        self.min_coverage = min_coverage
        self.min_region_size = min_region_size
        self.coverage_threshold = coverage_threshold
        self.bin_size = bin_size
        
        # Open BAM file
        try:
            self.bamfile = pysam.AlignmentFile(bam_path, "rb")
            if not self.bamfile.has_index():
                print(f"WARNING: BAM file is not indexed. Creating index...", file=sys.stderr)
                pysam.index(bam_path)
                self.bamfile = pysam.AlignmentFile(bam_path, "rb")
        except Exception as e:
            print(f"ERROR: Cannot open BAM file {bam_path}: {e}", file=sys.stderr)
            sys.exit(1)
    
    def read_bed_file(self, bed_path):
        """
        Read BED file and return list of regions
        
        Returns:
            List of tuples (chrom, start, end, name)
        """
        regions = []
        try:
            with open(bed_path, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith('#') or line.startswith('track'):
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) < 3:
                        print(f"WARNING: Skipping invalid BED line {line_num}: {line}", file=sys.stderr)
                        continue
                    
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    name = parts[3] if len(parts) > 3 else f"{chrom}:{start}-{end}"
                    
                    regions.append((chrom, start, end, name))
            
            print(f"Loaded {len(regions)} regions from {bed_path}")
            return regions
            
        except Exception as e:
            print(f"ERROR: Cannot read BED file {bed_path}: {e}", file=sys.stderr)
            sys.exit(1)
    
    def get_coverage_array(self, chrom, start, end):
        """
        Get per-base coverage array for a region
        
        Returns:
            numpy array of coverage values
        """
        try:
            # Get coverage using pileup
            coverage = np.zeros(end - start, dtype=np.int32)
            
            for pileupcolumn in self.bamfile.pileup(chrom, start, end, truncate=True):
                pos = pileupcolumn.pos
                if start <= pos < end:
                    coverage[pos - start] = pileupcolumn.n
            
            return coverage
            
        except Exception as e:
            print(f"WARNING: Error getting coverage for {chrom}:{start}-{end}: {e}", file=sys.stderr)
            return np.zeros(end - start, dtype=np.int32)
    
    def analyze_region_coverage(self, chrom, start, end, name):
        """
        Analyze coverage for a single region
        
        Returns:
            dict with coverage statistics and analysis results
        """
        region_size = end - start
        
        # Get coverage array
        coverage = self.get_coverage_array(chrom, start, end)
        
        if len(coverage) == 0:
            return {
                'chrom': chrom,
                'start': start,
                'end': end,
                'name': name,
                'size': region_size,
                'mean_coverage': 0.0,
                'median_coverage': 0.0,
                'min_coverage': 0,
                'max_coverage': 0,
                'bases_covered': 0,
                'coverage_fraction': 0.0,
                'is_analyzable': False,
                'reason': 'No coverage data',
                'coverable_regions': []
            }
        
        # Calculate statistics
        mean_cov = np.mean(coverage)
        median_cov = np.median(coverage)
        min_cov = np.min(coverage)
        max_cov = np.max(coverage)
        
        # Count bases with sufficient coverage
        covered_bases = np.sum(coverage >= self.min_coverage)
        coverage_fraction = covered_bases / len(coverage) if len(coverage) > 0 else 0.0
        
        # Find contiguous covered regions
        coverable_regions = self._find_coverable_regions(coverage, start)
        
        # Determine if region is analyzable
        is_analyzable = (
            coverage_fraction >= self.coverage_threshold and
            region_size >= self.min_region_size and
            mean_cov >= self.min_coverage
        )
        
        # Determine reason if not analyzable
        reason = "OK"
        if not is_analyzable:
            reasons = []
            if coverage_fraction < self.coverage_threshold:
                reasons.append(f"Low coverage fraction ({coverage_fraction:.2%} < {self.coverage_threshold:.2%})")
            if region_size < self.min_region_size:
                reasons.append(f"Region too small ({region_size} bp < {self.min_region_size} bp)")
            if mean_cov < self.min_coverage:
                reasons.append(f"Low mean coverage ({mean_cov:.2f} < {self.min_coverage})")
            reason = "; ".join(reasons)
        
        return {
            'chrom': chrom,
            'start': start,
            'end': end,
            'name': name,
            'size': region_size,
            'mean_coverage': mean_cov,
            'median_coverage': median_cov,
            'min_coverage': min_cov,
            'max_coverage': max_cov,
            'bases_covered': covered_bases,
            'coverage_fraction': coverage_fraction,
            'is_analyzable': is_analyzable,
            'reason': reason,
            'coverable_regions': coverable_regions
        }
    
    def _find_coverable_regions(self, coverage, offset):
        """
        Find contiguous regions with sufficient coverage
        
        Args:
            coverage: numpy array of coverage values
            offset: genomic offset (start position)
        
        Returns:
            List of (start, end) tuples for coverable regions
        """
        coverable_regions = []
        
        # Find stretches with sufficient coverage
        sufficient = coverage >= self.min_coverage
        
        in_region = False
        region_start = 0
        
        for i, is_sufficient in enumerate(sufficient):
            if is_sufficient and not in_region:
                # Start of a coverable region
                region_start = i
                in_region = True
            elif not is_sufficient and in_region:
                # End of a coverable region
                region_end = i
                region_size = region_end - region_start
                
                if region_size >= self.min_region_size:
                    coverable_regions.append((
                        offset + region_start,
                        offset + region_end
                    ))
                
                in_region = False
        
        # Handle case where region extends to the end
        if in_region:
            region_end = len(sufficient)
            region_size = region_end - region_start
            if region_size >= self.min_region_size:
                coverable_regions.append((
                    offset + region_start,
                    offset + region_end
                ))
        
        return coverable_regions
    
    def analyze_all_regions(self, regions):
        """
        Analyze coverage for all regions
        
        Returns:
            List of analysis results
        """
        results = []
        
        print(f"\nAnalyzing {len(regions)} regions...")
        for i, (chrom, start, end, name) in enumerate(regions, 1):
            if i % 10 == 0 or i == len(regions):
                print(f"  Progress: {i}/{len(regions)} regions", end='\r')
            
            result = self.analyze_region_coverage(chrom, start, end, name)
            results.append(result)
        
        print(f"  Progress: {len(regions)}/{len(regions)} regions - Complete!")
        return results
    
    def write_bed_output(self, results, output_path, include_partial=True):
        """
        Write coverable regions to BED file
        
        Args:
            results: List of analysis results
            output_path: Output BED file path
            include_partial: If True, include partially covered regions
        """
        try:
            with open(output_path, 'w') as f:
                # Write header
                f.write("# MD Coverage Analysis - Coverable Regions\n")
                f.write(f"# Min coverage: {self.min_coverage}\n")
                f.write(f"# Min region size: {self.min_region_size} bp\n")
                f.write(f"# Coverage threshold: {self.coverage_threshold:.2%}\n")
                f.write("#chrom\tstart\tend\tname\tscore\tstrand\tmean_cov\tmedian_cov\tcoverage_frac\n")
                
                regions_written = 0
                
                for result in results:
                    if result['is_analyzable']:
                        # Write fully analyzable region
                        f.write(f"{result['chrom']}\t{result['start']}\t{result['end']}\t"
                               f"{result['name']}\t1000\t.\t"
                               f"{result['mean_coverage']:.2f}\t{result['median_coverage']:.2f}\t"
                               f"{result['coverage_fraction']:.4f}\n")
                        regions_written += 1
                    
                    elif include_partial and result['coverable_regions']:
                        # Write partially covered regions
                        for j, (cov_start, cov_end) in enumerate(result['coverable_regions'], 1):
                            cov_size = cov_end - cov_start
                            # Calculate coverage for this subregion
                            subcov = self.get_coverage_array(result['chrom'], cov_start, cov_end)
                            mean_cov = np.mean(subcov) if len(subcov) > 0 else 0
                            median_cov = np.median(subcov) if len(subcov) > 0 else 0
                            covered = np.sum(subcov >= self.min_coverage) if len(subcov) > 0 else 0
                            cov_frac = covered / len(subcov) if len(subcov) > 0 else 0
                            
                            name = f"{result['name']}_partial{j}"
                            f.write(f"{result['chrom']}\t{cov_start}\t{cov_end}\t"
                                   f"{name}\t500\t.\t"
                                   f"{mean_cov:.2f}\t{median_cov:.2f}\t{cov_frac:.4f}\n")
                            regions_written += 1
            
            print(f"\nWrote {regions_written} coverable regions to {output_path}")
            
        except Exception as e:
            print(f"ERROR: Cannot write output BED file: {e}", file=sys.stderr)
            sys.exit(1)
    
    def generate_report(self, results, report_path):
        """
        Generate detailed coverage report
        """
        try:
            # Create DataFrame
            df = pd.DataFrame(results)
            
            # Summary statistics
            total_regions = len(results)
            analyzable_regions = sum(1 for r in results if r['is_analyzable'])
            partial_regions = sum(1 for r in results if not r['is_analyzable'] and r['coverable_regions'])
            non_analyzable = total_regions - analyzable_regions - partial_regions
            
            with open(report_path, 'w') as f:
                f.write("=" * 80 + "\n")
                f.write("MD Coverage Analysis Report\n")
                f.write("=" * 80 + "\n\n")
                
                f.write(f"BAM File: {self.bam_path}\n")
                f.write(f"Analysis Parameters:\n")
                f.write(f"  - Minimum coverage: {self.min_coverage}x\n")
                f.write(f"  - Minimum region size: {self.min_region_size} bp\n")
                f.write(f"  - Coverage threshold: {self.coverage_threshold:.2%}\n")
                f.write(f"  - Bin size (Wisecondor): {self.bin_size} bp\n\n")
                
                f.write(f"Summary:\n")
                f.write(f"  - Total regions: {total_regions}\n")
                f.write(f"  - Fully analyzable: {analyzable_regions} ({analyzable_regions/total_regions*100:.1f}%)\n")
                f.write(f"  - Partially covered: {partial_regions} ({partial_regions/total_regions*100:.1f}%)\n")
                f.write(f"  - Not analyzable: {non_analyzable} ({non_analyzable/total_regions*100:.1f}%)\n\n")
                
                # Coverage statistics
                mean_coverages = [r['mean_coverage'] for r in results]
                f.write(f"Coverage Statistics:\n")
                f.write(f"  - Mean coverage (avg): {np.mean(mean_coverages):.2f}x\n")
                f.write(f"  - Mean coverage (median): {np.median(mean_coverages):.2f}x\n")
                f.write(f"  - Mean coverage (min): {np.min(mean_coverages):.2f}x\n")
                f.write(f"  - Mean coverage (max): {np.max(mean_coverages):.2f}x\n\n")
                
                # Per-region details
                f.write("=" * 80 + "\n")
                f.write("Per-Region Analysis\n")
                f.write("=" * 80 + "\n\n")
                
                for result in results:
                    status = "✓ PASS" if result['is_analyzable'] else "✗ FAIL"
                    f.write(f"{status} {result['name']} ({result['chrom']}:{result['start']}-{result['end']})\n")
                    f.write(f"  Size: {result['size']:,} bp\n")
                    f.write(f"  Mean coverage: {result['mean_coverage']:.2f}x\n")
                    f.write(f"  Median coverage: {result['median_coverage']:.2f}x\n")
                    f.write(f"  Coverage range: {result['min_coverage']}-{result['max_coverage']}x\n")
                    f.write(f"  Bases covered: {result['bases_covered']:,}/{result['size']:,} "
                           f"({result['coverage_fraction']:.2%})\n")
                    
                    if not result['is_analyzable']:
                        f.write(f"  Reason: {result['reason']}\n")
                        
                        if result['coverable_regions']:
                            f.write(f"  Coverable subregions: {len(result['coverable_regions'])}\n")
                            for i, (start, end) in enumerate(result['coverable_regions'], 1):
                                f.write(f"    {i}. {result['chrom']}:{start}-{end} ({end-start:,} bp)\n")
                    
                    f.write("\n")
                
                # Wisecondor bin analysis
                f.write("=" * 80 + "\n")
                f.write("Wisecondor/WisecondorX Bin Analysis\n")
                f.write("=" * 80 + "\n\n")
                f.write(f"Bin size: {self.bin_size:,} bp\n\n")
                
                for result in results:
                    num_bins = result['size'] // self.bin_size
                    f.write(f"{result['name']}: {num_bins} bins "
                           f"({result['size']:,} bp / {self.bin_size:,} bp)\n")
            
            print(f"Detailed report written to {report_path}")
            
        except Exception as e:
            print(f"ERROR: Cannot generate report: {e}", file=sys.stderr)
    
    def close(self):
        """Close BAM file"""
        if hasattr(self, 'bamfile'):
            self.bamfile.close()


def main():
    parser = argparse.ArgumentParser(
        description='MD Coverage Check Tool - Analyze BAM coverage for Wisecondor/WisecondorX analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python3 md_coverage_check.py --bam sample.proper_paired.bam --bed MD_Target_8.bed --output coverable.bed
  
  # With custom thresholds
  python3 md_coverage_check.py --bam sample.bam --bed targets.bed --output output.bed \\
      --min-coverage 5 --coverage-threshold 0.95 --min-region-size 5000
  
  # Generate detailed report
  python3 md_coverage_check.py --bam sample.bam --bed targets.bed --output output.bed --report
        """
    )
    
    # Required arguments
    parser.add_argument('--bam', required=True, help='Input BAM file (proper_paired.bam recommended)')
    parser.add_argument('--bed', required=True, help='Input BED file with target regions')
    parser.add_argument('--output', '-o', required=True, help='Output BED file for coverable regions')
    
    # Optional arguments
    parser.add_argument('--min-coverage', type=int, default=1,
                       help='Minimum read depth per base (default: 1)')
    parser.add_argument('--min-region-size', type=int, default=1000,
                       help='Minimum region size in bp (default: 1000)')
    parser.add_argument('--coverage-threshold', type=float, default=0.9,
                       help='Minimum fraction of bases covered (0.0-1.0, default: 0.9)')
    parser.add_argument('--bin-size', type=int, default=200000,
                       help='Bin size for Wisecondor analysis (default: 200000)')
    parser.add_argument('--report', action='store_true',
                       help='Generate detailed coverage report')
    parser.add_argument('--no-partial', action='store_true',
                       help='Do not include partially covered regions in output')
    
    args = parser.parse_args()
    
    # Validate arguments
    if not os.path.exists(args.bam):
        print(f"ERROR: BAM file not found: {args.bam}", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(args.bed):
        print(f"ERROR: BED file not found: {args.bed}", file=sys.stderr)
        sys.exit(1)
    
    if args.coverage_threshold < 0 or args.coverage_threshold > 1:
        print(f"ERROR: Coverage threshold must be between 0 and 1", file=sys.stderr)
        sys.exit(1)
    
    # Print configuration
    print("=" * 80)
    print("MD Coverage Check Tool")
    print("=" * 80)
    print(f"BAM file: {args.bam}")
    print(f"BED file: {args.bed}")
    print(f"Output file: {args.output}")
    print(f"")
    print(f"Parameters:")
    print(f"  - Minimum coverage: {args.min_coverage}x")
    print(f"  - Minimum region size: {args.min_region_size:,} bp")
    print(f"  - Coverage threshold: {args.coverage_threshold:.2%}")
    print(f"  - Bin size: {args.bin_size:,} bp")
    print(f"  - Include partial regions: {not args.no_partial}")
    print("=" * 80)
    
    # Initialize analyzer
    analyzer = CoverageAnalyzer(
        args.bam,
        min_coverage=args.min_coverage,
        min_region_size=args.min_region_size,
        coverage_threshold=args.coverage_threshold,
        bin_size=args.bin_size
    )
    
    try:
        # Read BED file
        regions = analyzer.read_bed_file(args.bed)
        
        # Analyze coverage
        results = analyzer.analyze_all_regions(regions)
        
        # Write output BED file
        analyzer.write_bed_output(results, args.output, include_partial=not args.no_partial)
        
        # Generate report if requested
        if args.report:
            report_path = args.output.replace('.bed', '_report.txt')
            analyzer.generate_report(results, report_path)
        
        # Print summary
        print("\n" + "=" * 80)
        print("Analysis Complete!")
        print("=" * 80)
        
        total = len(results)
        analyzable = sum(1 for r in results if r['is_analyzable'])
        partial = sum(1 for r in results if not r['is_analyzable'] and r['coverable_regions'])
        
        print(f"Total regions: {total}")
        print(f"  ✓ Fully analyzable: {analyzable} ({analyzable/total*100:.1f}%)")
        print(f"  ⚠ Partially covered: {partial} ({partial/total*100:.1f}%)")
        print(f"  ✗ Not analyzable: {total-analyzable-partial} ({(total-analyzable-partial)/total*100:.1f}%)")
        print(f"\nWisecondor/WisecondorX Analysis:")
        if analyzable == total:
            print("  ✓ All regions are suitable for analysis")
        elif analyzable > 0:
            print(f"  ⚠ {analyzable}/{total} regions are suitable for analysis")
            print(f"  ⚠ Consider using coverable regions from {args.output}")
        else:
            print("  ✗ No regions have sufficient coverage")
            print("  ✗ Check BAM quality and coverage depth")
        
    finally:
        analyzer.close()


if __name__ == '__main__':
    main()

