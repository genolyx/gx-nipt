#!/usr/bin/env python3
"""
detect_yff.py - Y chromosome-based fetal fraction calculation
This script calculates fetal fraction based on Y chromosome coverage in a maternal blood sample.
For male fetuses, Y chromosome reads can be used to estimate fetal fraction.

Usage:
    python3 detect_yff.py --sample SAMPLE_ID --bam BAM_FILE --ybed CHR_Y_BED --abed AUTOSOME_BED [--output OUTPUT_FILE]

Arguments:
    --sample    Sample ID
    --bam       Input BAM file
    --ybed      BED file for Y chromosome target regions
    --abed      BED file for autosomal control regions
    --output    Output file (default: sample_id_yff.txt)
"""

import argparse
import os
import sys
import subprocess
import numpy as np
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('detect_yff')

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Calculate Y-based fetal fraction')
    parser.add_argument('--sample', required=True, help='Sample ID')
    parser.add_argument('--bam', required=True, help='Input BAM file')
    parser.add_argument('--ybed', required=True, help='BED file for Y chromosome target regions')
    parser.add_argument('--abed', required=True, help='BED file for autosomal control regions')
    parser.add_argument('--output', help='Output file')
    return parser.parse_args()

def count_reads_in_regions(bam_file, bed_file):
    """Count reads in the specified regions using samtools"""
    try:
        cmd = ['samtools', 'bedcov', bed_file, bam_file]
        logger.info(f"Running command: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        # Parse the bedcov output
        counts = []
        regions = []
        total_bases = 0
        
        for line in result.stdout.strip().split('\n'):
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                #count = int(parts[3])

                count = int(parts[-1])
                
                region_size = end - start
                total_bases += region_size
                
                regions.append((chrom, start, end, region_size))
                counts.append(count)
        
        return counts, regions, total_bases
    
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running samtools: {e}")
        logger.error(f"STDERR: {e.stderr}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        sys.exit(1)

def calculate_normalized_coverage(counts, regions, total_bases):
    """Calculate normalized coverage for regions"""
    normalized_counts = []
    
    for i, count in enumerate(counts):
        region_size = regions[i][3]
        # Normalize by region size
        if region_size > 0:
            normalized_count = count / region_size
            normalized_counts.append(normalized_count)
    
    return normalized_counts

def calculate_yff(y_counts, y_regions, y_total_bases, a_counts, a_regions, a_total_bases):
    """Calculate Y-based fetal fraction"""
    # Normalize coverage
    y_normalized = calculate_normalized_coverage(y_counts, y_regions, y_total_bases)
    a_normalized = calculate_normalized_coverage(a_counts, a_regions, a_total_bases)
    
    # Calculate median coverage
    y_median_coverage = np.median(y_normalized)
    a_median_coverage = np.median(a_normalized)
    
    logger.info(f"Median Y coverage: {y_median_coverage:.6f}")
    logger.info(f"Median Autosome coverage: {a_median_coverage:.6f}")
    
    # Check if we have reasonable values to proceed
    if a_median_coverage <= 0:
        logger.warning("Autosomal median coverage is zero or negative. Cannot calculate YFF.")
        return 0, 0, "FAILED"
    
    # Calculate the ratio of Y to autosome, adjusting for male genome normalization factor
    # In a normal male genome, Y chromosome is present as a single copy while autosomes are in two copies
    # So we multiply the Y/A ratio by 2 to get the fetal fraction for a male fetus
    y_to_a_ratio = y_median_coverage / a_median_coverage
    
    # Calculate fetal fraction assuming a male fetus
    # For a male fetus in maternal blood, the formula is:
    # FF = 2 * (Y/A ratio) since only the male fetus contributes Y chromosome
    fetal_fraction = 2 * y_to_a_ratio * 100  # Convert to percentage
    
    # Determine gender based on Y coverage
    # If Y/A ratio is very low, it likely indicates a female fetus or noise
    # Typically, a Y/A ratio threshold of 0.01-0.02 is used
    gender_threshold = 0.01
    
    if y_to_a_ratio < gender_threshold:
        gender = "FEMALE"
        fetal_fraction = 0  # Set FF to 0 for female fetuses as YFF is not applicable
        status = "NA"  # Not applicable for female fetuses
    else:
        gender = "MALE"
        status = "OK"
    
    return fetal_fraction, y_to_a_ratio, gender, status

def main():
    """Main function"""
    args = parse_args()
    
    # Set default output file if not provided
    if not args.output:
        args.output = f"{args.sample}_yff.txt"
    
    logger.info(f"Processing sample: {args.sample}")
    logger.info(f"BAM file: {args.bam}")
    logger.info(f"Y chromosome BED: {args.ybed}")
    logger.info(f"Autosome BED: {args.abed}")
    
    # Check input files exist
    for file_path in [args.bam, args.ybed, args.abed]:
        if not os.path.exists(file_path):
            logger.error(f"File not found: {file_path}")
            sys.exit(1)
    
    # Count reads in Y chromosome regions
    logger.info("Counting reads in Y chromosome regions...")
    y_counts, y_regions, y_total_bases = count_reads_in_regions(args.bam, args.ybed)
    
    # Count reads in autosomal regions
    logger.info("Counting reads in autosomal control regions...")
    a_counts, a_regions, a_total_bases = count_reads_in_regions(args.bam, args.abed)
    
    # Calculate fetal fraction
    logger.info("Calculating Y-based fetal fraction...")
    fetal_fraction, y_to_a_ratio, gender, status = calculate_yff(
        y_counts, y_regions, y_total_bases, 
        a_counts, a_regions, a_total_bases
    )
    
    # Write results to output file
    logger.info(f"Writing results to {args.output}")
    with open(args.output, 'w') as f:
        f.write(f"sample_id\ty_fraction\ty_to_a_ratio\tgender\tstatus\n")
        f.write(f"{args.sample}\t{fetal_fraction:.6f}\t{y_to_a_ratio:.6f}\t{gender}\t{status}\n")
    
    logger.info(f"Y-based fetal fraction: {fetal_fraction:.2f}%")
    logger.info(f"Y/A ratio: {y_to_a_ratio:.6f}")
    logger.info(f"gender: {gender}")
    logger.info(f"Status: {status}")
    logger.info("Done.")

if __name__ == "__main__":
    main()
