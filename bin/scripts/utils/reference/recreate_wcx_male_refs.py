#!/usr/bin/env python3
"""
Docker 컨테이너 내에서 WisecondorX Male Reference만 재생성
orig_M_200k_proper_paired.npz
fetus_M_200k_of.npz
"""

import pandas as pd
import subprocess
import sys
import os

def create_wcx_male_reference(sample_list_file, groups, output_dir):
    """WisecondorX Male reference만 생성"""
    
    print(f"=== WisecondorX Male Reference 재생성 ===")
    print(f"Groups: {groups}")
    print(f"Output: {output_dir}")
    
    # 출력 디렉토리 생성
    os.makedirs(output_dir, exist_ok=True)
    
    # 샘플 리스트 읽기
    df = pd.read_csv(sample_list_file, sep='\t')
    
    # Male 샘플만 필터링
    male_df = df[df['fetal_gender(gd_2)'] == 'XY']
    
    print(f"Male samples: {len(male_df)}")
    
    # WisecondorX 실행
    wcx_bin = "/opt/conda/envs/nipt/bin/WisecondorX"
    success = True
    
    for group in groups:
        print(f"\n{'='*60}")
        print(f"Processing group: {group}")
        print(f"{'='*60}")
        
        # NPZ 파일 수집
        male_files = []
        
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
        
        print(f"Male NPZ files found: {len(male_files)}")
        
        if len(male_files) == 0:
            print(f"⚠ No male NPZ files found for {group}")
            continue
        
        # 파일명 suffix 결정
        suffix = "proper_paired" if group == "orig" else "of"
        
        # Male reference 생성
        output_npz = f"{output_dir}/{group}_M_200k_{suffix}.npz"
        print(f"\nCreating male reference: {output_npz}")
        
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
            # 파일 크기 확인
            if os.path.exists(output_npz):
                size_mb = os.path.getsize(output_npz) / (1024 * 1024)
                print(f"  Size: {size_mb:.1f} MB")
        else:
            print(f"✗ Error creating {output_npz}")
            print(f"  stderr: {result.stderr}")
            success = False
    
    if not success:
        return False
    
    print("\n" + "="*60)
    print("=== WCX Male Reference 재생성 완료 ===")
    print("="*60)
    return True

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: python recreate_wcx_male_refs.py <sample_list_file> <output_dir> <group1> [group2] ...")
        print("  groups: orig, fetus, mom")
        print("")
        print("Example:")
        print("  python recreate_wcx_male_refs.py samples.tsv /refs/ucl_new/WCX orig fetus")
        sys.exit(1)
    
    sample_list_file = sys.argv[1]
    output_dir = sys.argv[2]
    groups = sys.argv[3:]  # 나머지 모든 인자를 group으로
    
    print(f"Sample list: {sample_list_file}")
    print(f"Output dir: {output_dir}")
    print(f"Groups: {groups}")
    
    success = create_wcx_male_reference(sample_list_file, groups, output_dir)
    sys.exit(0 if success else 1)
