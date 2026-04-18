#!/usr/bin/env python3
"""
Docker 컨테이너 내에서 Wisecondor Reference 생성
"""

import pandas as pd
import subprocess
import sys
import os

def create_wc_reference(sample_list_file, group, output_dir):
    """Wisecondor reference 생성 (Female 먼저, Male 나중 순서로)"""
    
    print(f"=== Wisecondor Reference 생성 ===")
    print(f"Group: {group}")
    print(f"Output: {output_dir}")
    
    # 출력 디렉토리 생성
    os.makedirs(output_dir, exist_ok=True)
    
    # 샘플 리스트 읽기
    df = pd.read_csv(sample_list_file, sep='\t')
    
    # 성별로 분류 (원본 스크립트와 동일한 순서: Female 먼저, Male 나중)
    female_df = df[df['fetal_gender(gd_2)'] == 'XX']
    male_df = df[df['fetal_gender(gd_2)'] == 'XY']
    
    print(f"Female samples: {len(female_df)}")
    print(f"Male samples: {len(male_df)}")
    
    # NPZ 파일 수집 (Female 먼저, Male 나중)
    female_files = []
    male_files = []
    
    # Female NPZ 파일 수집
    for _, row in female_df.iterrows():
        sample_id = row['sample_id']
        sample_dir = row['sample_dir']
        
        # 호스트 경로를 Docker 컨테이너 경로로 변환
        sample_dir = sample_dir.replace('/home/ken/ken-nipt/analysis', '/analysis')
        
        if group == 'orig':
            npz_file = f"{sample_dir}/Output_WC/{sample_id}.wc.of_orig.npz"
        elif group == 'fetus':
            npz_file = f"{sample_dir}/Output_WC/{sample_id}.wc.of_fetus.npz"
        elif group == 'mom':
            npz_file = f"{sample_dir}/Output_WC/{sample_id}.wc.of_mom.npz"
        
        if os.path.exists(npz_file):
            female_files.append(npz_file)
    
    # Male NPZ 파일 수집
    for _, row in male_df.iterrows():
        sample_id = row['sample_id']
        sample_dir = row['sample_dir']
        
        # 호스트 경로를 Docker 컨테이너 경로로 변환
        sample_dir = sample_dir.replace('/home/ken/ken-nipt/analysis', '/analysis')
        
        if group == 'orig':
            npz_file = f"{sample_dir}/Output_WC/{sample_id}.wc.of_orig.npz"
        elif group == 'fetus':
            npz_file = f"{sample_dir}/Output_WC/{sample_id}.wc.of_fetus.npz"
        elif group == 'mom':
            npz_file = f"{sample_dir}/Output_WC/{sample_id}.wc.of_mom.npz"
        
        if os.path.exists(npz_file):
            male_files.append(npz_file)
    
    print(f"Female NPZ files found: {len(female_files)}")
    print(f"Male NPZ files found: {len(male_files)}")
    
    # Female + Male 순서로 결합 (원본 스크립트와 동일)
    all_files = female_files + male_files
    
    if len(all_files) == 0:
        print("✗ Error: No NPZ files found!")
        return False
    
    # Wisecondor 실행 (combined only)
    wc_path = "/opt/wisecondor/wisecondor.py"
    
    # 파일명 suffix 결정
    suffix = "proper_paired" if group == "orig" else "of"
    output_npz = f"{output_dir}/{group}_200k_{suffix}.npz"
    
    print(f"\nCreating combined reference...")
    
    # Wisecondor uses single dash for binsize: -binsize (not --binsize)
    cmd = [
        "/usr/bin/python2",
        wc_path,
        "newref",
        *all_files,
        output_npz,
        "-binsize", "200000"
    ]
    
    print(f"Running: /usr/bin/python2 {wc_path} newref ... ({len(all_files)} files)")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode == 0:
        print(f"✓ Created: {output_npz}")
    else:
        print(f"✗ Error: {result.stderr}")
        return False
    
    print("\n=== WC Reference 생성 완료 ===")
    return True

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python create_wc_docker.py <sample_list_file> <group> <output_dir>")
        print("  group: orig, fetus, mom")
        sys.exit(1)
    
    sample_list_file = sys.argv[1]
    group = sys.argv[2]
    output_dir = sys.argv[3]
    
    success = create_wc_reference(sample_list_file, group, output_dir)
    sys.exit(0 if success else 1)
