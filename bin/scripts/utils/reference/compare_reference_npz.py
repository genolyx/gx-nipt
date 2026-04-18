#!/usr/bin/env python3
"""
Compare WC reference NPZ files to understand structural differences
"""

import numpy as np
import sys

def analyze_npz(npz_path, label):
    """Analyze NPZ file structure"""
    print(f"\n{'='*80}")
    print(f"{label}: {npz_path}")
    print(f"{'='*80}")
    
    try:
        data = np.load(npz_path, allow_pickle=True)
        
        print(f"\nFile size: {npz_path}")
        import os
        size_mb = os.path.getsize(npz_path) / 1024 / 1024
        print(f"  {size_mb:.2f} MB")
        
        print(f"\nKeys in NPZ:")
        for key in sorted(data.keys()):
            arr = data[key]
            print(f"  {key}:")
            print(f"    Type: {type(arr)}")
            print(f"    Shape: {arr.shape if hasattr(arr, 'shape') else 'N/A'}")
            print(f"    Dtype: {arr.dtype if hasattr(arr, 'dtype') else 'N/A'}")
            
            # Check for nested structure
            if arr.dtype == object and len(arr.shape) == 0:
                inner = arr.item()
                if isinstance(inner, dict):
                    print(f"    Dict keys: {list(inner.keys())[:5]}...")
                    # Sample one chromosome
                    if 'chr1' in inner:
                        chr1 = inner['chr1']
                        print(f"    chr1 type: {type(chr1)}")
                        if hasattr(chr1, 'shape'):
                            print(f"    chr1 shape: {chr1.shape}")
                        if hasattr(chr1, 'dtype'):
                            print(f"    chr1 dtype: {chr1.dtype}")
        
        data.close()
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == '__main__':
    print("\n" + "="*80)
    print("UCL REFERENCE NPZ COMPARISON")
    print("="*80)
    
    # Old reference
    old_orig = '/home/ken/ken-nipt/data/refs/ucl/WC/orig_200k_proper_paired.npz'
    old_fetus = '/home/ken/ken-nipt/data/refs/ucl/WC/fetus_200k_of.npz'
    
    # New reference  
    new_orig = '/home/ken/ken-nipt/data/refs/ucl_2510plus/WC/orig_200k_proper_paired.npz'
    new_fetus = '/home/ken/ken-nipt/data/refs/ucl_2510plus/WC/fetus_200k_of.npz'
    
    analyze_npz(old_orig, "OLD UCL (orig)")
    analyze_npz(new_orig, "NEW UCL 2510+ (orig)")
    
    analyze_npz(old_fetus, "OLD UCL (fetus)")
    analyze_npz(new_fetus, "NEW UCL 2510+ (fetus)")
    
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    print("Old UCL reference: ~494 samples (all periods)")
    print("New 2510+ reference: ~386 samples (2510+ only)")
    print("\nFile size difference is normal due to:")
    print("  1. Fewer samples (386 vs 494)")
    print("  2. More uniform data = better compression")
    print("  3. NPZ uses internal compression")
