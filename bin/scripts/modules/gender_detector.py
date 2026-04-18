"""
Gender Detection Module for NIPT Analysis Pipeline

This module provides gender detection functionality using multiple methods:
- gd_1: R script-based Y chromosome analysis
- gd_2: Shell script-based samtools analysis  
- gd_3: YFF (Y-chromosome Fetal Fraction) based analysis
- gd_4: YFF2 (wig normalization) based analysis

Author: NIPT Analysis Team
"""

import pysam
import pandas as pd
import subprocess
import json
import os
import numpy as np
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

class GenderDetector:
    """
    4가지 방법을 사용한 Gender Detection 클래스
    """
    
    def __init__(self, config_file="pipeline_config.json"):
        """
        Gender Detection 클래스 초기화
        
        Args:
            config_file (str): pipeline 전체 설정이 담긴 config 파일 경로
        """
        self.config = self.load_config(config_file)
        # FF_Gender_Config 섹션만 추출
        self.gender_config = self.config.get('FF_Gender_Config', {})
        
    def load_config(self, config_file):
        """
        pipeline config 파일에서 설정값들을 로드
        """
        default_config = {
            "FF_Gender_Config": {
                "fragment_size_cutoff": 160,
                "gd_1_threshold": 0.4,
                "gd_2_threshold": 0.4,
                "gd_3_threshold": 0.01,
                "gd_4_threshold": 0.02
            }
        }
        
        if os.path.exists(config_file):
            with open(config_file, 'r') as f:
                config = json.load(f)
                # 기본값으로 누락된 키들을 채움
                if 'FF_Gender_Config' not in config:
                    config['FF_Gender_Config'] = default_config['FF_Gender_Config']
                else:
                    for key, value in default_config['FF_Gender_Config'].items():
                        if key not in config['FF_Gender_Config']:
                            config['FF_Gender_Config'][key] = value
                return config
        else:
            # config 파일이 없으면 기본값으로 생성
            with open(config_file, 'w') as f:
                json.dump(default_config, f, indent=2)
            return default_config
    
    def get_threshold(self, method):
        """threshold 값 가져오기"""
        return self.gender_config.get(f'{method}_threshold', 0.4)
    
    def get_fragment_size_cutoff(self):
        """fragment size cutoff 값 가져오기"""
        return self.gender_config.get('fragment_size_cutoff', 160)
    
    def gd_1_detection(self, bam_file):
        """
        첫 번째 gender detection 방법 (R script 기반)
        Y 염색체 특정 영역의 coverage 비율을 계산하여 성별 판정
        
        Args:
            bam_file (str): proper_paired.bam 파일 경로
            
        Returns:
            tuple: (gender_value, gender_detect, gender_det)
        """
        try:
            # BAM 파일 열기
            bamfile = pysam.AlignmentFile(bam_file, "rb")
            
            # Y 염색체 전체 영역에서 reads 수집 (2650001-59050000)
            chrY_coverage = {}
            
            # Y 염색체의 유효 영역에서 각 위치별 coverage 계산
            for pileupcolumn in bamfile.pileup("chrY", 2650001, 59050000):
                pos = pileupcolumn.pos + 1  # 1-based coordinate
                if 2650001 <= pos <= 59050000:
                    chrY_coverage[pos] = pileupcolumn.n  # coverage depth
            
            bamfile.close()
            
            if not chrY_coverage:
                return 0.0, "FEMALE", "XX"
            
            # 전체 Y 영역의 총 coverage
            total_y_coverage = sum(chrY_coverage.values())
            logger.info(f"total_y_coverage : {total_y_coverage}")
            
            # 특정 14개 위치에서의 coverage 합계
            specific_positions = [7650001, 7750001, 7800001, 8400001, 8450001, 8500001, 
                                8550001, 8600001, 15500001, 18900001, 22250001, 22450001, 
                                22900001, 23600001]
            
            specific_coverage = 0
            for pos in specific_positions:
                if pos in chrY_coverage:
                    specific_coverage += chrY_coverage[pos]
            
            logger.info(f"specific_coverage : {specific_coverage}")
            # Gender value 계산: (특정 위치 coverage / 전체 Y coverage) * 100
            if total_y_coverage == 0:
                gender_value = 0.0
            else:
                gender_value = (specific_coverage / total_y_coverage) * 100
            
            # threshold 기반 성별 판정
            threshold = self.get_threshold('gd_1')
            gender_detect = "MALE" if gender_value > threshold else "FEMALE"
            gender_det = "XY" if gender_value > threshold else "XX"
            
            logger.info(f"{gender_value}, {gender_detect}, {gender_det}")
            return gender_value, gender_detect, gender_det
            
        except Exception as e:
            logger.error(f"Error in gd_1_detection: {e}")
            return 0.0, "FEMALE", "XX"
    
    def gd_2_detection(self, sample_name, bam_file, bed_file):
        """
        두 번째 gender detection 방법 (shell script 기반)
        
        Args:
            sample_name (str): 샘플 ID
            bam_file (str): proper_paired.bam 파일 경로
            bed_file (str): BED 파일 경로
            
        Returns:
            tuple: (gender_value, gender_detect, gender_det)
        """
        try:
            # samtools를 사용해서 특정 영역의 reads 카운트
            cmd_gv = f"samtools view -c -L {bed_file} {bam_file}"
            result_gv = subprocess.run(cmd_gv, shell=True, capture_output=True, text=True)
            gv_count = int(result_gv.stdout.strip())
            
            # Y 염색체 전체 영역의 reads 카운트
            cmd_whole = f"samtools view -c {bam_file} chrY:2650001-59050000"
            result_whole = subprocess.run(cmd_whole, shell=True, capture_output=True, text=True)
            whole_count = int(result_whole.stdout.strip())

            logger.info(f"gv_count : {gv_count}, whole_count : {whole_count}")
            
            if whole_count == 0:
                gender_value = 0.0
            else:
                gender_value = (gv_count / whole_count) * 100
            
            # 결과 파일에 저장
            output_file = f"/Work/NIPT/analysis/{sample_name}/{sample_name}.new_gv.txt"
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            with open(output_file, 'w') as f:
                f.write(f"{gender_value}\n")
            
            # threshold 기반 성별 판정
            threshold = self.get_threshold('gd_2')
            gender_detect = "MALE" if gender_value > threshold else "FEMALE"
            gender_det = "XY" if gender_value > threshold else "XX"
            
            return gender_value, gender_detect, gender_det
            
        except Exception as e:
            logger.error(f"Error in gd_2_detection: {e}")
            return 0.0, "FEMALE", "XX"
    
    def gd_3_from_ff_result(self, ff_result):
        """
        calculate_fetal_fraction 함수의 calculate_yff 결과를 받아서 gender detection
        
        Args:
            ff_result (dict): calculate_yff 함수의 리턴값
                
        Returns:
            tuple: (gender_value, gender_detect, gender_det)
        """
        try:
            if ff_result is None:
                return 0.0, "FEMALE", "XX"
            
            # FF 결과에서 값 추출
            yff_value = ff_result.get('yff_value', 0.0)
            y_to_a_ratio = ff_result.get('y_to_a_ratio', 0.0)
            
            # config에서 threshold 가져오기
            threshold = self.get_threshold('gd_3')

            logger.info(f"{yff_value}, {y_to_a_ratio}, {threshold}")
            
            # Y/A ratio 기반으로 성별 판정
            if y_to_a_ratio > threshold:
                gender_detect = "MALE"
                gender_det = "XY"
                # male인 경우 실제 FF 값 사용
                final_ff = yff_value
            else:
                gender_detect = "FEMALE" 
                gender_det = "XX"
                # female인 경우 FF를 0으로 설정 (YFF가 의미없으므로)
                final_ff = 0.0
            
            return final_ff, gender_detect, gender_det
            
        except Exception as e:
            logger.error(f"Error in gd_3_from_ff_result: {e}")
            return 0.0, "FEMALE", "XX"
    
    def gd_4_from_ff_result(self, ff_result):
        """
        calculate_fetal_fraction 함수의 calculate_yff2 결과를 받아서 gender detection
        
        Args:
            ff_result (dict): calculate_yff2 함수의 리턴값
                
        Returns:
            tuple: (gender_value, gender_detect, gender_det)
        """
        try:
            if ff_result is None or ff_result.get('status') != 'OK':
                return 0.0, "FEMALE", "XX"
            
            # FF 결과에서 값 추출
            ff_value = ff_result.get('FF_chrY_adjusted', 0.0)
            uar_y = ff_result.get('UAR_Y', 0.0)
            
            # config에서 threshold 가져오기
            threshold = self.get_threshold('gd_4')
            
            logger.info(f"{ff_value}, {uar_y}, {threshold}")
            # UAR_Y 기반으로 성별 판정
            if uar_y > threshold:
                gender_detect = "MALE"
                gender_det = "XY"
                # male인 경우 실제 FF 값 사용
                final_ff = ff_value
            else:
                gender_detect = "FEMALE"
                gender_det = "XX"
                # female인 경우 FF를 0으로 설정
                final_ff = 0.0
            
            return final_ff, gender_detect, gender_det
            
        except Exception as e:
            logger.error(f"Error in gd_4_from_ff_result: {e}")
            return 0.0, "FEMALE", "XX"
    
    def gd_3_detection(self, ff_result):
        """
        세 번째 gender detection 방법 (calculate_yff 결과 기반)
        
        Args:
            ff_result (dict): calculate_fetal_fraction에서 받은 calculate_yff 결과
            
        Returns:
            tuple: (gender_value, gender_detect, gender_det)
        """
        return self.gd_3_from_ff_result(ff_result)
    
    def gd_4_detection(self, ff_result):
        """
        네 번째 gender detection 방법 (calculate_yff2 결과 기반)
        
        Args:
            ff_result (dict): calculate_fetal_fraction에서 받은 calculate_yff2 결과
            
        Returns:
            tuple: (gender_value, gender_detect, gender_det)
        """
        return self.gd_4_from_ff_result(ff_result)
    
    def detect_gender_from_ff(self, sample_name, bam_file, ff_yff_result=None, ff_yff2_result=None, bed_file="/Work/NIPT/data/bed/common/gender_y.bed", save_results=True, output_dir="/Work/NIPT/analysis"):
        """
        calculate_fetal_fraction 함수의 결과를 활용한 전체 gender detection 수행
        
        Args:
            sample_name (str): 샘플 ID
            bam_file (str): proper_paired.bam 파일 경로
            ff_yff_result (dict): calculate_yff 결과 (optional)
            ff_yff2_result (dict): calculate_yff2 결과 (optional)
            bed_file (str): BED 파일 경로 (gd_2용)
            save_results (bool): 결과를 파일로 저장할지 여부
            output_dir (str): 출력 디렉토리
            
        Returns:
            pandas.DataFrame: gender detection 결과
        """
        results = {}
        
        # gd_1 실행 
        gd1_value, gd1_gender, gd1_det = self.gd_1_detection(bam_file)
        results['gd_1'] = {
            'value': gd1_value,
            'gender': gd1_det
        }
        
        # gd_2 실행 (기존 shell script 기반)
        gd2_value, gd2_gender, gd2_det = self.gd_2_detection(sample_name, bam_file, bed_file)
        results['gd_2'] = {
            'value': gd2_value,
            'gender': gd2_det
        }
        
        # gd_3 실행 (calculate_yff 결과 기반) - ff_yff_result가 제공된 경우에만
        if ff_yff_result is not None:
            gd3_value, gd3_gender, gd3_det = self.gd_3_detection(ff_yff_result)
            results['gd_3'] = {
                'value': gd3_value,
                'gender': gd3_det
            }
        
        # gd_4 실행 (calculate_yff2 결과 기반) - ff_yff2_result가 제공된 경우에만
        if ff_yff2_result is not None:
            gd4_value, gd4_gender, gd4_det = self.gd_4_detection(ff_yff2_result)
            results['gd_4'] = {
                'value': gd4_value,
                'gender': gd4_det
            }
        
        # DataFrame으로 변환
        df_results = pd.DataFrame.from_dict(results, orient='index')
        
        # 결과 저장 (옵션)
        if save_results:
            self.save_gender_results(sample_name, df_results, ff_yff_result, ff_yff2_result, output_dir)
        
        return df_results
    
    def detect_gender(self, sample_name, bam_file, bed_file, ff_yff_result=None, ff_yff2_result=None):
        """
        전체 gender detection 수행 (기존 방식과 호환성 유지)
        
        Args:
            sample_name (str): 샘플 ID
            bam_file (str): proper_paired.bam 파일 경로
            bed_file (str): BED 파일 경로 (gd_2용)
            ff_yff_result (dict): calculate_yff 결과 (optional)
            ff_yff2_result (dict): calculate_yff2 결과 (optional)
            
        Returns:
            pandas.DataFrame: gender detection 결과
        """
        return self.detect_gender_from_ff(sample_name, bam_file, ff_yff_result, ff_yff2_result, bed_file)
    
    def save_gender_results(self, sample_name, gender_results, ff_yff_result=None, ff_yff2_result=None, output_dir="/Work/NIPT/analysis"):
        """
        Gender detection 결과를 Output_FF 디렉토리에 저장
        
        Args:
            sample_name (str): 샘플 ID
            gender_results (DataFrame): gender detection 결과 (gd_1, gd_2, gd_3, gd_4 rows)
            ff_yff_result (dict): YFF 결과 (optional)
            ff_yff2_result (dict): YFF2 결과 (optional) 
            output_dir (str): 출력 디렉토리
        """
        try:
            # Output_FF 디렉토리 생성
            output_ff_dir = Path(output_dir)
            output_ff_dir.mkdir(parents=True, exist_ok=True)
            
            # 파일명: <sample_name>.gender_value.txt
            output_file = output_ff_dir / f"{sample_name}.gender_value.txt"
            
            # gender_results DataFrame을 그대로 저장
            # Index: gd_1, gd_2, gd_3, gd_4
            # Columns: value, gender
            gender_results.to_csv(output_file, sep='\t', index=True, header=True)
            
            print(f"Gender detection results saved: {output_file}")
            
        except Exception as e:
            print(f"Error saving gender results: {e}")
    
    def load_gender_results(self, sample_name, output_dir="/Work/NIPT/analysis"):
        """
        저장된 gender detection 결과를 로드
        
        Args:
            sample_name (str): 샘플 ID
            output_dir (str): 출력 디렉토리
            
        Returns:
            DataFrame: gender detection 결과
        """
        try:
            output_ff_dir = Path(output_dir)
            input_file = output_ff_dir / f"{sample_name}.gender_value.txt"
            
            if input_file.exists():
                df = pd.read_csv(input_file, sep='\t', index_col=0)
                return df
            else:
                print(f"Gender results file not found: {input_file}")
                return None
                
        except Exception as e:
            print(f"Error loading gender results: {e}")
            return None


# 편의 함수들
def create_gender_detector(config_file=None):
    """
    GenderDetector 인스턴스를 생성하는 편의 함수
    
    Args:
        config_file (str): config 파일 경로 (optional)
        
    Returns:
        GenderDetector: 초기화된 detector 인스턴스
    """
    if config_file:
        return GenderDetector(config_file)
    else:
        return GenderDetector()


def run_gender_detection(sample_name, bam_file, config_file=None, ff_yff_result=None, ff_yff2_result=None):
    """
    간단한 gender detection 실행 함수
    
    Args:
        sample_name (str): 샘플 ID
        bam_file (str): BAM 파일 경로
        config_file (str): config 파일 경로 (optional)
        ff_yff_result (dict): YFF 결과 (optional)
        ff_yff2_result (dict): YFF2 결과 (optional)
        
    Returns:
        DataFrame: gender detection 결과
    """
    detector = create_gender_detector(config_file)
    return detector.detect_gender_from_ff(
        sample_name=sample_name,
        bam_file=bam_file,
        ff_yff_result=ff_yff_result,
        ff_yff2_result=ff_yff2_result
    )


if __name__ == "__main__":
    # 테스트 예시
    detector = GenderDetector()
    print("Gender Detector initialized successfully!")
    print(f"gd_1 threshold: {detector.get_threshold('gd_1')}")
    print(f"gd_2 threshold: {detector.get_threshold('gd_2')}")
    print(f"gd_3 threshold: {detector.get_threshold('gd_3')}")
    print(f"gd_4 threshold: {detector.get_threshold('gd_4')}")
