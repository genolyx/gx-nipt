#!/usr/bin/env python3
"""
Docker 컨테이너 내에서 WisecondorX Reference 생성
"""

import pandas as pd
import subprocess
import sys
import os

def create_wcx_reference(sample_list_file, group, output_dir):
    """WisecondorX reference 생성"""
    
    print(f"=== WisecondorX Reference 생성 ===")
    print(f"Group: {group}")
    print(f"Output: {output_dir}")
    
    # 출력 디렉토리 생성
    os.makedirs(output_dir, exist_ok=True)
    
    # 샘플 리스트 읽기
    df = pd.read_csv(sample_list_file, sep='\t')
    
    # 성별로 분류
    male_df = df[df['fetal_gender(gd_2)'] == 'XY']
    female_df = df[df['fetal_gender(gd_2)'] == 'XX']
    
    print(f"Male samples: {len(male_df)}")
    print(f"Female samples: {len(female_df)}")
    
    # NPZ 파일 수집
    male_files = []
    female_files = []
    
    for _, row in male_df.iterrows():
        sample_id = row['sample_id']
        sample_dir = row['sample_dir']
        
        # 호스트 경로를 Docker 컨테이너 경로로 변환
        sample_dir = sample_dir.replace('/home/ken/ken-nipt/analysis', '/analysis')
        
        if group == 'orig':
            npz_file = f"{sample_dir}/Output_WCX/{sample_id}.wcx.of_orig.npz"
        elif group == 'fetus':
            npz_file = f"{sample_dir}/Output_WCX/{sample_id}.wcx.of_fetus.npz"
        elif group == 'mom':
            npz_file = f"{sample_dir}/Output_WCX/{sample_id}.wcx.of_mom.npz"
        
        if os.path.exists(npz_file):
            male_files.append(npz_file)
    
    for _, row in female_df.iterrows():
        sample_id = row['sample_id']
        sample_dir = row['sample_dir']
        
        # 호스트 경로를 Docker 컨테이너 경로로 변환
        sample_dir = sample_dir.replace('/home/ken/ken-nipt/analysis', '/analysis')
        
        if group == 'orig':
            npz_file = f"{sample_dir}/Output_WCX/{sample_id}.wcx.of_orig.npz"
        elif group == 'fetus':
            npz_file = f"{sample_dir}/Output_WCX/{sample_id}.wcx.of_fetus.npz"
        elif group == 'mom':
            npz_file = f"{sample_dir}/Output_WCX/{sample_id}.wcx.of_mom.npz"
        
        if os.path.exists(npz_file):
            female_files.append(npz_file)
    
    print(f"Male NPZ files found: {len(male_files)}")
    print(f"Female NPZ files found: {len(female_files)}")
    
    # WisecondorX 실행
    wcx_bin = "/opt/conda/envs/nipt/bin/WisecondorX"
    success = True
    
    # 파일명 suffix 결정
    suffix = "proper_paired" if group == "orig" else "of"
    
    # mom은 combined만 생성하고, orig/fetus는 M, F, combined 모두 생성
    if group != 'mom':
        # 1. Male reference 생성 (orig, fetus만)
        if len(male_files) > 0:
            output_npz = f"{output_dir}/{group}_M_200k_{suffix}.npz"
            print(f"\nCreating male reference...")
            
            # Male-only: --nipt --yfrac 0 함께 사용 (generate_ref_gender.sh 방식)
            cmd = [
                wcx_bin, "newref",
                *male_files,
                output_npz,
                "--binsize", "200000",
                "--nipt",
                "--yfrac", "0"
            ]
            
            print(f"Running: {' '.join(cmd[:3])} ... ({len(male_files)} files)")
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                print(f"✓ Created: {output_npz}")
            else:
                print(f"✗ Error: {result.stderr}")
                success = False
        
        # 2. Female reference 생성 (orig, fetus만)
        if len(female_files) > 0:
            output_npz = f"{output_dir}/{group}_F_200k_{suffix}.npz"
            print(f"\nCreating female reference...")
            
            # Female-only: --nipt --yfrac 0 함께 사용 (generate_ref_gender.sh 방식)
            cmd = [
                wcx_bin, "newref",
                *female_files,
                output_npz,
                "--binsize", "200000",
                "--nipt",
                "--yfrac", "0"
            ]
            
            print(f"Running: {' '.join(cmd[:3])} ... ({len(female_files)} files)")
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                print(f"✓ Created: {output_npz}")
            else:
                print(f"✗ Error: {result.stderr}")
                success = False
    
    # 3. Combined reference 생성 (M+F)
    if len(male_files) > 0 and len(female_files) > 0:
        output_npz = f"{output_dir}/{group}_200k_{suffix}.npz"
        print(f"\nCreating combined reference (M+F)...")
        
        if group == 'mom':
            # mom은 female 중심으로 (--yfrac 0 사용, --nipt 사용 안 함)
            cmd = [
                wcx_bin, "newref",
                *female_files,
                *male_files,
                output_npz,
                "--binsize", "200000",
                "--yfrac", "0"
            ]
        else:
            # orig, fetus는 --nipt만 사용 (--yfrac 0 사용 안 함!)
            cmd = [
                wcx_bin, "newref",
                *female_files,
                *male_files,
                output_npz,
                "--binsize", "200000",
                "--nipt"
            ]
        
        print(f"Running: {' '.join(cmd[:3])} ... ({len(female_files)+len(male_files)} files)")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            print(f"✓ Created: {output_npz}")
        else:
            print(f"✗ Error: {result.stderr}")
            success = False
    
    if not success:
        return False
    
    print("\n=== WCX Reference 생성 완료 ===")
    return True

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python create_wcx_docker.py <sample_list_file> <group> <output_dir>")
        print("  group: orig, fetus, mom")
        sys.exit(1)
    
    sample_list_file = sys.argv[1]
    group = sys.argv[2]
    output_dir = sys.argv[3]
    
    success = create_wcx_reference(sample_list_file, group, output_dir)
    sys.exit(0 if success else 1)
