#!/usr/bin/env python3
"""
Create filtered UCL sample list for 2510+ samples
Filter by month >= 2510 and basic QC thresholds
"""

import pandas as pd
import sys
import os

def filter_ucl_2510_plus():
    """Filter UCL samples for 2510+ with quality thresholds"""
    
    # Read all samples
    input_file = 'data/refs/ucl/reference_make/reference_sample_list_UCL_all.tsv'
    output_file = 'data/refs/ucl/reference_make/reference_sample_list_UCL_2510plus.tsv'
    
    print(f"\n{'='*80}")
    print(f"FILTERING UCL SAMPLES FOR 2510+ REFERENCE")
    print(f"{'='*80}\n")
    
    # Read TSV
    df = pd.read_csv(input_file, sep='\t')
    print(f"Total samples in file: {len(df)}")
    
    # Filter by month >= 2510
    df['month_int'] = df['month'].astype(int)
    df_2510 = df[df['month_int'] >= 2510].copy()
    print(f"Samples from 2510+: {len(df_2510)}")
    
    # Show month distribution
    print(f"\nMonth distribution:")
    for month in sorted(df_2510['month'].unique()):
        count = len(df_2510[df_2510['month'] == month])
        print(f"  {month}: {count} samples")
    
    # Basic QC filtering (same as original create_reference.py defaults)
    print(f"\n{'='*80}")
    print(f"APPLYING BASIC QC FILTERS")
    print(f"{'='*80}\n")
    
    initial_count = len(df_2510)
    
    # Filter criteria
    filters = []
    
    # 1. Mapping rate >= 95%
    if 'mapping_rate(%)' in df_2510.columns:
        before = len(df_2510)
        df_2510 = df_2510[df_2510['mapping_rate(%)'] >= 95.0]
        filtered = before - len(df_2510)
        if filtered > 0:
            print(f"Mapping rate < 95%: removed {filtered} samples")
            filters.append(f"mapping_rate >= 95%: {filtered} removed")
    
    # 2. Duplication rate <= 30% (very permissive, to keep failed libraries for analysis)
    # Note: For production, you might want stricter like <= 15%
    if 'duplication_rate(%)' in df_2510.columns:
        before = len(df_2510)
        df_2510 = df_2510[df_2510['duplication_rate(%)'] <= 30.0]
        filtered = before - len(df_2510)
        if filtered > 0:
            print(f"Duplication rate > 30%: removed {filtered} samples")
            filters.append(f"duplication_rate <= 30%: {filtered} removed")
    
    # 3. GC content within reasonable range (35-50%)
    if 'GC_content(%)' in df_2510.columns:
        before = len(df_2510)
        df_2510 = df_2510[(df_2510['GC_content(%)'] >= 35.0) & (df_2510['GC_content(%)'] <= 50.0)]
        filtered = before - len(df_2510)
        if filtered > 0:
            print(f"GC content out of range (35-50%): removed {filtered} samples")
            filters.append(f"GC_content 35-50%: {filtered} removed")
    
    # 4. SeqFF >= 3% (minimum fetal fraction)
    if 'SeqFF' in df_2510.columns:
        before = len(df_2510)
        df_2510 = df_2510[df_2510['SeqFF'] >= 3.0]
        filtered = before - len(df_2510)
        if filtered > 0:
            print(f"SeqFF < 3%: removed {filtered} samples")
            filters.append(f"SeqFF >= 3%: {filtered} removed")
    
    final_count = len(df_2510)
    removed_count = initial_count - final_count
    
    print(f"\n{'='*80}")
    print(f"FILTERING SUMMARY")
    print(f"{'='*80}\n")
    print(f"Initial samples (2510+): {initial_count}")
    print(f"After QC filtering: {final_count}")
    print(f"Removed: {removed_count} samples ({removed_count/initial_count*100:.1f}%)")
    print(f"Retention rate: {final_count/initial_count*100:.1f}%")
    
    # Quality statistics
    print(f"\n{'='*80}")
    print(f"QUALITY STATISTICS OF FILTERED SAMPLES")
    print(f"{'='*80}\n")
    
    metrics = {
        'duplication_rate(%)': 'Duplication Rate',
        'mapping_rate(%)': 'Mapping Rate',
        'GC_content(%)': 'GC Content',
        'mean_coverageData(X)': 'Coverage',
        'SeqFF': 'Fetal Fraction'
    }
    
    for col, name in metrics.items():
        if col in df_2510.columns:
            data = df_2510[col].dropna()
            if len(data) > 0:
                print(f"{name}:")
                print(f"  Mean: {data.mean():.2f}")
                print(f"  Median: {data.median():.2f}")
                print(f"  Std: {data.std():.2f}")
                print(f"  Range: {data.min():.2f} - {data.max():.2f}")
                print()
    
    # Gender distribution
    if 'fetal_gender(gd_2)' in df_2510.columns:
        print(f"Gender distribution:")
        gender_counts = df_2510['fetal_gender(gd_2)'].value_counts()
        for gender, count in gender_counts.items():
            print(f"  {gender}: {count} samples ({count/len(df_2510)*100:.1f}%)")
        print()
    
    # Save filtered list
    # Drop the temporary month_int column
    df_2510 = df_2510.drop('month_int', axis=1)
    
    df_2510.to_csv(output_file, sep='\t', index=False)
    print(f"{'='*80}")
    print(f"SAVED FILTERED SAMPLE LIST")
    print(f"{'='*80}\n")
    print(f"Output file: {output_file}")
    print(f"Total samples: {len(df_2510)}")
    print(f"\nReady for reference creation!")
    print(f"\nNext step:")
    print(f"  python3 bin/scripts/utils/reference/create_reference.py \\")
    print(f"    --sample-list {output_file} \\")
    print(f"    --labcode ucl \\")
    print(f"    --output-dir data/refs/ucl_2510plus \\")
    print(f"    --reference-source data/refs/ucl \\")
    print(f"    --ref-type wc wcx")
    print()

if __name__ == '__main__':
    filter_ucl_2510_plus()
