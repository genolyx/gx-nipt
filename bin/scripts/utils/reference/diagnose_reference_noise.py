#!/usr/bin/env python3
"""
Diagnose why 2510+ reference might show more noise
Check sample quality distribution in both references
"""

import pandas as pd
import numpy as np

def analyze_reference_samples():
    """Compare sample quality between old and new references"""
    
    print("\n" + "="*80)
    print("REFERENCE SAMPLE QUALITY COMPARISON")
    print("="*80)
    
    # Load sample lists
    old_file = 'data/refs/ucl/reference_make/reference_sample_list_UCL_filtered.tsv'
    new_file = 'data/refs/ucl/reference_make/reference_sample_list_UCL_2510plus_filtered.tsv'
    
    df_old = pd.read_csv(old_file, sep='\t')
    df_new = pd.read_csv(new_file, sep='\t')
    
    print(f"\n1. SAMPLE COUNTS:")
    print(f"   Old UCL: {len(df_old)} samples")
    print(f"   New 2510+: {len(df_new)} samples")
    
    # Month distribution
    print(f"\n2. MONTH DISTRIBUTION:")
    print(f"\n   Old UCL:")
    for month in sorted(df_old['month'].unique()):
        count = len(df_old[df_old['month'] == month])
        print(f"     {month}: {count} samples")
    
    print(f"\n   New 2510+:")
    for month in sorted(df_new['month'].unique()):
        count = len(df_new[df_new['month'] == month])
        print(f"     {month}: {count} samples")
    
    # Quality metrics comparison
    print(f"\n3. QUALITY METRICS:")
    metrics = ['duplication_rate(%)', 'mean_coverageData(X)', 'GC_content(%)', 'SeqFF']
    
    for metric in metrics:
        if metric in df_old.columns and metric in df_new.columns:
            old_data = df_old[metric].dropna()
            new_data = df_new[metric].dropna()
            
            print(f"\n   {metric}:")
            print(f"     Old UCL:   {old_data.mean():.2f} ± {old_data.std():.2f} (range: {old_data.min():.2f}-{old_data.max():.2f})")
            print(f"     New 2510+: {new_data.mean():.2f} ± {new_data.std():.2f} (range: {new_data.min():.2f}-{new_data.max():.2f})")
            print(f"     CV old: {(old_data.std()/old_data.mean()*100):.1f}%")
            print(f"     CV new: {(new_data.std()/new_data.mean()*100):.1f}%")
    
    # Check for problematic samples in new reference
    print(f"\n4. HIGH DUPLICATION SAMPLES:")
    old_high_dup = df_old[df_old['duplication_rate(%)'] > 15]
    new_high_dup = df_new[df_new['duplication_rate(%)'] > 15]
    
    print(f"   Old UCL: {len(old_high_dup)} samples with dup > 15%")
    print(f"   New 2510+: {len(new_high_dup)} samples with dup > 15%")
    
    if len(new_high_dup) > 0:
        print(f"\n   High-dup samples in 2510+:")
        for _, row in new_high_dup.head(10).iterrows():
            print(f"     {row['sample_id']}: {row['duplication_rate(%)']:.1f}% (month: {row['month']})")
    
    # Coverage analysis
    print(f"\n5. COVERAGE ANALYSIS:")
    old_cov = df_old['mean_coverageData(X)'].dropna()
    new_cov = df_new['mean_coverageData(X)'].dropna()
    
    print(f"   Old UCL coverage:")
    print(f"     Mean: {old_cov.mean():.3f}X")
    print(f"     Std: {old_cov.std():.3f}X")
    print(f"     CV: {(old_cov.std()/old_cov.mean()*100):.1f}%")
    
    print(f"   New 2510+ coverage:")
    print(f"     Mean: {new_cov.mean():.3f}X")
    print(f"     Std: {new_cov.std():.3f}X")
    print(f"     CV: {(new_cov.std()/new_cov.mean()*100):.1f}%")
    
    # Critical insight
    print(f"\n{'='*80}")
    print("CRITICAL INSIGHTS:")
    print("="*80)
    
    print(f"\n⚠️  POTENTIAL ISSUES:")
    print(f"   1. Sample count: {len(df_old)} -> {len(df_new)} ({(len(df_new)/len(df_old)*100):.1f}%)")
    print(f"   2. High-dup samples: {len(old_high_dup)} -> {len(new_high_dup)}")
    
    # Check if 2510+ has high variation
    if new_cov.std() > old_cov.std():
        print(f"   3. ⚠️  Coverage variation INCREASED ({old_cov.std():.3f} -> {new_cov.std():.3f})")
    else:
        print(f"   3. ✓ Coverage variation decreased ({old_cov.std():.3f} -> {new_cov.std():.3f})")
    
    new_dup = df_new['duplication_rate(%)'].dropna()
    old_dup = df_old['duplication_rate(%)'].dropna()
    if new_dup.std() > old_dup.std():
        print(f"   4. ⚠️  Duplication variation INCREASED ({old_dup.std():.2f} -> {new_dup.std():.2f})")
    else:
        print(f"   4. ✓ Duplication variation decreased ({old_dup.std():.2f} -> {new_dup.std():.2f})")
    
    print(f"\n💡 POSSIBLE REASONS FOR MORE NOISE:")
    print(f"   - Fewer samples = less statistical power to smooth outliers")
    print(f"   - 2510+ still contains high-dup samples (up to {new_dup.max():.1f}%)")
    print(f"   - Coverage variation may not have improved enough")
    print(f"   - Reference might need stricter filtering")
    
    # Recommendation
    print(f"\n{'='*80}")
    print("RECOMMENDATIONS:")
    print("="*80)
    print(f"1. Try stricter filtering:")
    print(f"   - Duplication rate <= 12% (currently: up to {new_dup.max():.1f}%)")
    print(f"   - This would give ~{len(df_new[df_new['duplication_rate(%)'] <= 12])} samples")
    print(f"\n2. Or use more samples:")
    print(f"   - Include 2509 if quality is acceptable")
    print(f"   - Check 2509 sample quality separately")
    print(f"\n3. Test with actual NIPT samples:")
    print(f"   - Run known positives/negatives")
    print(f"   - Compare Z-scores between references")

if __name__ == '__main__':
    analyze_reference_samples()
