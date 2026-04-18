#!/usr/bin/env python3
"""
=========================================================
NIPT Reference Creation Pipeline
=========================================================
통합 reference 생성 파이프라인

4가지 reference type 지원:
1. EZD - UAR & Z-score 기반 threshold
2. PRIZM - Mean & SD 기반 Z-score
3. WC - Wisecondor npz
4. WCX - WisecondorX npz

각각 orig, fetus, mom으로 분류하여 생성

Author: Hyukjung Kwon
Date: 2025-01-06
=========================================================
"""

import argparse
import logging
import os
import sys
import json
import subprocess
from pathlib import Path
from typing import List, Dict, Tuple
import pandas as pd
import numpy as np

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Constants
# Auto-detect workspace root (script is in <workspace>/bin/scripts/utils/reference/)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
WORKSPACE_ROOT = os.environ.get(
    "WORKSPACE",
    os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.dirname(SCRIPT_DIR)))
    ),
)
DATA_DIR = os.path.join(WORKSPACE_ROOT, "data")
ANALYSIS_DIR = os.path.join(WORKSPACE_ROOT, "analysis")
BIN_DIR = os.path.join(WORKSPACE_ROOT, "bin")

# External tools
WCX_BIN = "/usr/bin/miniconda3/bin/WisecondorX"
WC_PATH = "/Work/NIPT/bin/wisecondor/wisecondor.py"
PYTHON2 = "/usr/bin/miniconda2/bin/python2"
PYTHON3 = "/usr/bin/miniconda3/bin/python3"


class SampleListGenerator:
    """샘플 리스트 생성 클래스 (QC + FF + Report 정보 수집)"""
    
    # QC parsing
    WANT_QC_KEYS = {
        "number_of_reads",
        "number_of_mapped_reads",
        "mapping_rate",
        "duplication_rate",
        "mean_mapping_quality",
        "mean_coverageData",
        "GC_content",
    }
    
    @staticmethod
    def strip_bom(s: str) -> str:
        return s.replace("\ufeff", "").replace("\u200b", "")
    
    @staticmethod
    def parse_qc_txt(qc_path: str) -> Dict[str, str]:
        """QC 파일 파싱"""
        out = {k: None for k in SampleListGenerator.WANT_QC_KEYS}
        if not os.path.exists(qc_path):
            return out
        
        with open(qc_path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = SampleListGenerator.strip_bom(line).strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t", 1)
                if len(parts) != 2:
                    continue
                key = parts[0].strip().rstrip(":")
                val = parts[1].strip()
                if key in out:
                    out[key] = val
        return out
    
    @staticmethod
    def parse_gender_gd2(path: str) -> str:
        """Gender 파싱"""
        if not os.path.exists(path):
            return ""
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = SampleListGenerator.strip_bom(line).strip()
                if not line or line.lower().startswith("value"):
                    continue
                parts = line.split()
                if len(parts) >= 3 and parts[0] == "gd_2":
                    return parts[2]
        return ""
    
    @staticmethod
    def parse_ff(path: str) -> Dict[str, str]:
        """Fetal Fraction 파싱"""
        out = {"Fragment_FF": "", "YFF_2": "", "SeqFF": "", "M-SeqFF": ""}
        if not os.path.exists(path):
            return out
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = SampleListGenerator.strip_bom(line).strip()
                if not line or line.lower() == "value":
                    continue
                parts = line.split()
                if len(parts) >= 2 and parts[0] in out:
                    out[parts[0]] = parts[1]
        return out
    
    @staticmethod
    def parse_report_json(report_path: str) -> Dict[str, str]:
        """Report JSON에서 Result, MDResult, Disease 추출"""
        out = {"Result": "", "MDResult": "", "Disease": ""}
        if not os.path.exists(report_path):
            return out
        try:
            with open(report_path, "r", encoding="utf-8") as f:
                j = json.load(f)
            out["Result"] = str(j.get("Result", ""))
            out["MDResult"] = str(j.get("MDResult", ""))
            
            # Extract items from trisomy_result
            trisomy_result = j.get("trisomy_result", [])
            if trisomy_result and isinstance(trisomy_result, list):
                items = [item.get("item", "") for item in trisomy_result if isinstance(item, dict) and "item" in item]
                out["Disease"] = ", ".join(items) if items else ""
        except Exception:
            pass
        return out
    
    @staticmethod
    def to_int(s: str) -> int:
        if s is None:
            return None
        try:
            return int(float(s.replace(",", "")))
        except Exception:
            return None
    
    @staticmethod
    def to_float(s: str) -> float:
        if s is None:
            return None
        import re
        x = s.replace(",", "").replace("%", "")
        x = re.sub(r"[Xx]$", "", x)
        try:
            return float(x)
        except Exception:
            return None
    
    @staticmethod
    def fmt(x: float, nd: int) -> str:
        return "" if x is None else f"{x:.{nd}f}"
    
    @staticmethod
    def generate_sample_list(analysis_dirs: List[str], output_file: str, 
                           prefix: str = "GNMF", output_root: str = None):
        """
        분석 디렉토리들을 순회하며 샘플 리스트 TSV 생성
        
        Args:
            analysis_dirs: 분석 디렉토리 목록
            output_file: 출력 TSV 파일 경로
            prefix: 샘플 ID prefix
            output_root: report JSON 위치
        """
        import re
        import glob
        
        if output_root is None:
            output_root = os.path.join(WORKSPACE_ROOT, "output")
        
        logger.info(f"Generating sample list from {len(analysis_dirs)} directories...")
        logger.info(f"Sample ID prefix: {prefix}")
        
        # Sample directory regex pattern
        SAMPLE_DIR_RE = re.compile(rf"^{re.escape(prefix)}\d{{8}}$")
        
        header = [
            "month", "sample_id", "sample_dir",
            "number_of_reads", "number_of_mapped_reads", "mapping_rate(%)",
            "duplication_rate(%)", "mean_mapping_quality", "mean_coverageData(X)",
            "GC_content(%)", "fetal_gender(gd_2)", "Fragment_FF", "YFF_2",
            "SeqFF", "M-SeqFF", "Result", "Disease", "MDResult",
        ]
        
        rows = []
        
        for base in analysis_dirs:
            base = os.path.abspath(base)
            if not os.path.isdir(base):
                logger.warning(f"Directory not found: {base}")
                continue
            
            month = os.path.basename(base.rstrip("/"))
            logger.info(f"Processing {month}...")
            
            sample_dirs = sorted(
                d for d in glob.glob(os.path.join(base, f"{prefix}*"))
                if os.path.isdir(d) and SAMPLE_DIR_RE.match(os.path.basename(d))
            )
            
            for sample_dir in sample_dirs:
                sample_id = os.path.basename(sample_dir)
                
                # QC
                qc_path = os.path.join(sample_dir, "Output_QC", f"{sample_id}.qc.txt")
                qc = SampleListGenerator.parse_qc_txt(qc_path)
                
                reads = SampleListGenerator.to_int(qc.get("number_of_reads"))
                mapped = SampleListGenerator.to_int(qc.get("number_of_mapped_reads"))
                maprate = SampleListGenerator.to_float(qc.get("mapping_rate"))
                duprate = SampleListGenerator.to_float(qc.get("duplication_rate"))
                mapq = SampleListGenerator.to_float(qc.get("mean_mapping_quality"))
                cov = SampleListGenerator.to_float(qc.get("mean_coverageData"))
                gc = SampleListGenerator.to_float(qc.get("GC_content"))
                
                # FF / gender
                ff_dir = os.path.join(sample_dir, "Output_FF")
                gender = SampleListGenerator.parse_gender_gd2(
                    os.path.join(ff_dir, f"{sample_id}.gender.txt"))
                ff = SampleListGenerator.parse_ff(
                    os.path.join(ff_dir, f"{sample_id}.fetal_fraction.txt"))
                
                # Report JSON
                report_path = os.path.join(output_root, month, sample_id, f"{sample_id}_report.json")
                report = SampleListGenerator.parse_report_json(report_path)
                
                rows.append([
                    month, sample_id, sample_dir,
                    "" if reads is None else str(reads),
                    "" if mapped is None else str(mapped),
                    SampleListGenerator.fmt(maprate, 2),
                    SampleListGenerator.fmt(duprate, 2),
                    SampleListGenerator.fmt(mapq, 4),
                    SampleListGenerator.fmt(cov, 4),
                    SampleListGenerator.fmt(gc, 2),
                    gender, ff["Fragment_FF"], ff["YFF_2"], ff["SeqFF"], ff["M-SeqFF"],
                    report["Result"], report["Disease"], report["MDResult"],
                ])
        
        # Write output
        os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
        with open(output_file, "w", encoding="utf-8") as w:
            w.write("\t".join(header) + "\n")
            for r in rows:
                w.write("\t".join(r) + "\n")
        
        logger.info(f"✓ Sample list created: {output_file} (n={len(rows)} samples)")
        return output_file


class ReferenceConfig:
    """Reference 생성을 위한 설정 클래스"""
    
    def __init__(self, labcode: str, output_base_dir: str = None, source_base_dir: str = None):
        self.labcode = labcode
        
        if output_base_dir is None:
            output_base_dir = os.path.join(DATA_DIR, "refs", labcode)
        
        self.output_base_dir = output_base_dir
        
        # Source directory for existing Positive samples (default to labcode dir)
        if source_base_dir is None:
            source_base_dir = os.path.join(DATA_DIR, "refs", labcode)
        self.source_base_dir = source_base_dir
        
        self.groups = ['orig', 'fetus', 'mom']
        
        # Output directories for each reference type
        self.ezd_dir = os.path.join(output_base_dir, "EZD")
        self.prizm_dir = os.path.join(output_base_dir, "PRIZM")
        self.wc_dir = os.path.join(output_base_dir, "WC")
        self.wcx_dir = os.path.join(output_base_dir, "WCX")
        
        # Source directories for existing Positive samples
        self.source_ezd_dir = os.path.join(source_base_dir, "EZD")
        self.source_prizm_dir = os.path.join(source_base_dir, "PRIZM")
        self.source_wc_dir = os.path.join(source_base_dir, "WC")
        self.source_wcx_dir = os.path.join(source_base_dir, "WCX")
        
        # Create directories
        for dir_path in [self.ezd_dir, self.prizm_dir, self.wc_dir, self.wcx_dir]:
            os.makedirs(dir_path, exist_ok=True)
            for group in self.groups:
                os.makedirs(os.path.join(dir_path, group), exist_ok=True)


class SampleSelector:
    """샘플 선별 클래스"""
    
    def __init__(self, sample_list_file: str):
        """
        Args:
            sample_list_file: TSV 파일 경로 (엑셀에서 추출한 샘플 정보)
        """
        self.sample_list_file = sample_list_file
        self.df = None
        self.selected_samples = {'M': [], 'F': []}  # Male, Female
        
    def load_samples(self):
        """샘플 리스트 로드"""
        logger.info(f"Loading sample list from: {self.sample_list_file}")
        
        try:
            self.df = pd.read_csv(self.sample_list_file, sep='\t')
            logger.info(f"Loaded {len(self.df)} samples")
            logger.info(f"Columns: {list(self.df.columns)}")
            return True
        except Exception as e:
            logger.error(f"Failed to load sample list: {e}")
            return False
    
    def filter_samples(self, 
                      exclude_high_risk: bool = True,
                      exclude_no_call: bool = True,
                      outlier_threshold_zscore: float = 5.0,
                      min_seqff: float = 4.0,
                      max_seqff: float = 30.0,
                      min_mapping_rate: float = 95.0,
                      max_duplication_rate: float = None,
                      min_gc_content: float = None,
                      max_gc_content: float = None,
                      min_reads: int = None):
        """
        샘플 필터링
        
        Args:
            exclude_high_risk: High Risk 샘플 제외
            exclude_no_call: No Call 샘플 제외
            outlier_threshold_zscore: Outlier 제거를 위한 z-score threshold
            min_seqff: 최소 SeqFF 값
            max_seqff: 최대 SeqFF 값
            min_mapping_rate: 최소 mapping rate (기본: 95.0%)
            max_duplication_rate: 최대 duplication rate (None이면 제한 없음)
            min_gc_content: 최소 GC content (None이면 제한 없음)
            max_gc_content: 최대 GC content (None이면 제한 없음)
            min_reads: 최소 read 수 (None이면 제한 없음)
        """
        if self.df is None:
            logger.error("Sample list not loaded. Call load_samples() first.")
            return False
        
        original_count = len(self.df)
        logger.info(f"\n{'='*60}")
        logger.info("Sample Filtering Started")
        logger.info(f"{'='*60}")
        logger.info(f"Initial sample count: {original_count}")
        
        # 1. Result 컬럼 필터링 (case-insensitive, exact match)
        if exclude_high_risk and 'Result' in self.df.columns:
            before = len(self.df)
            # Exact match after stripping whitespace and converting to lowercase
            self.df = self.df[self.df['Result'].str.strip().str.lower() != 'high risk']
            logger.info(f"Exclude Result='High Risk' (case-insensitive): {before} -> {len(self.df)} ({before - len(self.df)} removed)")
        
        if exclude_no_call and 'Result' in self.df.columns:
            before = len(self.df)
            # Exact match after stripping whitespace and converting to lowercase
            self.df = self.df[self.df['Result'].str.strip().str.lower() != 'no call']
            logger.info(f"Exclude Result='No Call' (case-insensitive): {before} -> {len(self.df)} ({before - len(self.df)} removed)")
        
        # 2. MDResult 컬럼 필터링 (case-insensitive, exact match)
        if exclude_high_risk and 'MDResult' in self.df.columns:
            before = len(self.df)
            # Exact match, handle NaN/None by filling with empty string first
            mask = self.df['MDResult'].fillna('').str.strip().str.lower() != 'high risk'
            self.df = self.df[mask]
            logger.info(f"Exclude MDResult='High Risk' (case-insensitive): {before} -> {len(self.df)} ({before - len(self.df)} removed)")
        
        if exclude_no_call and 'MDResult' in self.df.columns:
            before = len(self.df)
            # Exact match, handle NaN/None by filling with empty string first
            mask = self.df['MDResult'].fillna('').str.strip().str.lower() != 'no call'
            self.df = self.df[mask]
            logger.info(f"Exclude MDResult='No Call' (case-insensitive): {before} -> {len(self.df)} ({before - len(self.df)} removed)")
        
        # 3. SeqFF 범위 필터링
        if 'SeqFF' in self.df.columns:
            before = len(self.df)
            self.df = self.df[(self.df['SeqFF'] >= min_seqff) & (self.df['SeqFF'] <= max_seqff)]
            logger.info(f"SeqFF range ({min_seqff}-{max_seqff}%): {before} -> {len(self.df)} ({before - len(self.df)} removed)")
        
        # 4. QC metrics 필터링 (mapping rate, duplication rate 등)
        if 'mapping_rate(%)' in self.df.columns:
            before = len(self.df)
            self.df = self.df[self.df['mapping_rate(%)'] >= min_mapping_rate]
            logger.info(f"Mapping rate >= {min_mapping_rate}%: {before} -> {len(self.df)} ({before - len(self.df)} removed)")
        
        # 5. Duplication rate 필터링
        if max_duplication_rate is not None and 'duplication_rate(%)' in self.df.columns:
            before = len(self.df)
            self.df = self.df[self.df['duplication_rate(%)'] <= max_duplication_rate]
            logger.info(f"Duplication rate <= {max_duplication_rate}%: {before} -> {len(self.df)} ({before - len(self.df)} removed)")
        
        # 6. GC content 필터링
        if min_gc_content is not None and max_gc_content is not None and 'GC_content(%)' in self.df.columns:
            before = len(self.df)
            self.df = self.df[(self.df['GC_content(%)'] >= min_gc_content) & (self.df['GC_content(%)'] <= max_gc_content)]
            logger.info(f"GC content range ({min_gc_content}-{max_gc_content}%): {before} -> {len(self.df)} ({before - len(self.df)} removed)")
        
        # 7. Total reads 필터링
        if min_reads is not None and 'number_of_reads' in self.df.columns:
            before = len(self.df)
            self.df = self.df[self.df['number_of_reads'] >= min_reads]
            logger.info(f"Total reads >= {min_reads:,}: {before} -> {len(self.df)} ({before - len(self.df)} removed)")
        
        # 8. Gender 분류
        self._classify_by_gender()
        
        logger.info(f"\n{'='*60}")
        logger.info(f"Filtering completed: {original_count} -> {len(self.df)} samples")
        logger.info(f"{'='*60}\n")
        
        return True
    
    def _classify_by_gender(self):
        """Gender별 샘플 분류"""
        if 'fetal_gender(gd_2)' in self.df.columns:
            male_samples = self.df[self.df['fetal_gender(gd_2)'] == 'XY']
            female_samples = self.df[self.df['fetal_gender(gd_2)'] == 'XX']
            
            logger.info(f"\nGender distribution:")
            logger.info(f"  Male (XY): {len(male_samples)} samples")
            logger.info(f"  Female (XX): {len(female_samples)} samples")
            
            self.selected_samples['M'] = male_samples['sample_id'].tolist()
            self.selected_samples['F'] = female_samples['sample_id'].tolist()
        else:
            logger.warning("fetal_gender(gd_2) column not found. Cannot classify by gender.")
    
    def save_filtered_list(self, output_file: str):
        """필터링된 샘플 리스트 저장"""
        if self.df is None or len(self.df) == 0:
            logger.error("No samples to save")
            return False
        
        try:
            self.df.to_csv(output_file, sep='\t', index=False)
            logger.info(f"Filtered sample list saved to: {output_file}")
            return True
        except Exception as e:
            logger.error(f"Failed to save filtered list: {e}")
            return False
    
    def get_samples_by_group(self, group: str) -> Dict[str, List[str]]:
        """
        그룹별 샘플 목록 반환
        
        Args:
            group: 'orig', 'fetus', 'mom'
        
        Returns:
            Dict with 'M' and 'F' keys containing sample IDs
        """
        # 모든 그룹에서 동일한 샘플 사용
        # 실제로는 orig/fetus/mom에 따라 다른 필터링이 필요할 수 있음
        return self.selected_samples


class EZDReferenceCreator:
    """EZD Reference 생성 클래스"""
    
    def __init__(self, config: ReferenceConfig, samples: Dict[str, List[str]], sample_info_df: pd.DataFrame = None):
        self.config = config
        self.samples = samples
        self.sample_info_df = sample_info_df  # Full sample info with Result, Disease columns
        
    def create_reference(self, group: str):
        """
        EZD reference 생성 (기존 reference 업데이트 방식)
        
        Args:
            group: 'orig', 'fetus', 'mom'
        """
        logger.info(f"\n{'='*60}")
        logger.info(f"Creating EZD Reference for {group}")
        logger.info(f"{'='*60}")
        
        output_dir = os.path.join(self.config.ezd_dir, group)
        os.makedirs(output_dir, exist_ok=True)
        
        # Disease to chromosome mapping
        disease_chr_map = {
            'T13': 'chr13', 'T18': 'chr18', 'T21': 'chr21',
            'T16': 'chr16', 'T22': 'chr22',
        }
        
        # 각 염색체별 데이터 업데이트
        all_chr_data = {}
        
        for chr_num in range(1, 23):
            chr_name = f"chr{chr_num}"
            chr_file = os.path.join(output_dir, f"{chr_name}.txt")
            
            # 기존 파일 로드 (Positive 샘플 유지) - source directory에서 읽기
            existing_positive = []
            source_chr_file = os.path.join(self.config.source_ezd_dir, group, f"{chr_name}.txt")
            if os.path.exists(source_chr_file):
                try:
                    existing_df = pd.read_csv(source_chr_file, sep='\t')
                    # Normalize column name to lowercase 'sample' if needed
                    if 'Sample' in existing_df.columns:
                        existing_df.rename(columns={'Sample': 'sample'}, inplace=True)
                    # P type 샘플만 선택
                    existing_positive = existing_df[existing_df['type'] == 'P'].to_dict('records')
                    logger.info(f"  {chr_name}: Loaded {len(existing_positive)} existing Positive samples from source")
                except Exception as e:
                    logger.debug(f"  {chr_name}: Could not load existing file: {e}")
            else:
                logger.debug(f"  {chr_name}: No source file found at {source_chr_file}")
            
            # Normal 샘플 수집 (Low Risk만)
            normal_rows = []
            
            # Male Low Risk samples
            for sample_id in self.samples['M']:
                uar, z = self._get_sample_chr_data(sample_id, chr_name, group)
                if uar is not None and z is not None:
                    normal_rows.append({
                        'sample': sample_id,
                        'type': 'N',
                        'UAR': uar,
                        'Z': z
                    })
            
            # Female Low Risk samples
            for sample_id in self.samples['F']:
                uar, z = self._get_sample_chr_data(sample_id, chr_name, group)
                if uar is not None and z is not None:
                    normal_rows.append({
                        'sample': sample_id,
                        'type': 'N',
                        'UAR': uar,
                        'Z': z
                    })
            
            # High Risk 샘플 추가 (Disease 매칭)
            if self.sample_info_df is not None:
                high_risk_samples = self._get_high_risk_samples_for_chr(chr_name, disease_chr_map)
                for sample_id, disease in high_risk_samples:
                    uar, z = self._get_sample_chr_data(sample_id, chr_name, group)
                    if uar is not None and z is not None:
                        existing_positive.append({
                            'sample': sample_id,
                            'type': 'P',
                            'UAR': uar,
                            'Z': z
                        })
            
            # Combine Normal + Positive
            all_rows = normal_rows + existing_positive
            
            if all_rows:
                chr_df = pd.DataFrame(all_rows)
                # Reorder columns (all should have 'sample' now)
                chr_df = chr_df[['sample', 'type', 'UAR', 'Z']]
                chr_df.to_csv(chr_file, sep='\t', index=False)
                logger.info(f"  {chr_name}: {len(normal_rows)} Normal + {len(existing_positive)} Positive = {len(all_rows)} total")
                all_chr_data[chr_name] = chr_df
        
        # Sex chromosomes (X, Y)
        self._create_sex_chromosome_files(output_dir, group)
        
        # Threshold 계산
        self._calculate_thresholds(all_chr_data, output_dir, group)
        
        # SCA config 생성
        self._create_sca_config(output_dir, group)
        
        logger.info(f"EZD reference for {group} created successfully")
    
    def _get_sample_chr_data(self, sample_id: str, chr_name: str, group: str) -> Tuple[float, float]:
        """
        샘플의 특정 염색체 UAR & Z-score 값 가져오기
        
        Returns:
            (UAR, Z-score) tuple
        """
        # 샘플 분석 결과 디렉토리 찾기
        sample_dir = self._find_sample_directory(sample_id)
        
        if not sample_dir:
            logger.debug(f"Sample directory not found for {sample_id}")
            return None, None
        
        # EZD 결과 파일 경로
        ezd_result_file = os.path.join(sample_dir, "Output_EZD", group, f"{group}_ezd_results.tsv")
        
        if not os.path.exists(ezd_result_file):
            logger.debug(f"EZD result file not found: {ezd_result_file}")
            return None, None
        
        try:
            df = pd.read_csv(ezd_result_file, sep='\t')
            chr_row = df[df['Chromosome'] == chr_name]
            
            if chr_row.empty:
                return None, None
            
            ur = float(chr_row['UAR'].iloc[0])
            z = float(chr_row['Z'].iloc[0])
            
            return ur, z
        except Exception as e:
            logger.debug(f"Failed to read EZD data for {sample_id} {chr_name}: {e}")
            return None, None
    
    def _find_sample_directory(self, sample_id: str) -> str:
        """
        샘플의 분석 디렉토리 찾기
        
        Returns:
            Full path to sample directory, or None if not found
        """
        # analysis/YYMM/SAMPLE_ID 형태로 찾기
        
        # Cache for performance
        if not hasattr(self, '_sample_dir_cache'):
            self._sample_dir_cache = {}
        
        if sample_id in self._sample_dir_cache:
            return self._sample_dir_cache[sample_id]
        
        # 월별 디렉토리 탐색
        for month_dir in os.listdir(ANALYSIS_DIR):
            month_path = os.path.join(ANALYSIS_DIR, month_dir)
            
            if not os.path.isdir(month_path):
                continue
            
            sample_path = os.path.join(month_path, sample_id)
            
            if os.path.isdir(sample_path):
                self._sample_dir_cache[sample_id] = sample_path
                return sample_path
        
        return None
    
    def _get_high_risk_samples_for_chr(self, chr_name: str, disease_chr_map: dict) -> list:
        """
        특정 염색체에 해당하는 High Risk 샘플 반환
        
        Args:
            chr_name: 'chr13', 'chr18', 'chr21', etc.
            disease_chr_map: Disease to chromosome mapping
            
        Returns:
            List of (sample_id, disease) tuples
        """
        high_risk_samples = []
        
        if self.sample_info_df is None:
            return high_risk_samples
        
        # Result가 "High Risk"인 샘플 필터링
        high_risk_df = self.sample_info_df[
            self.sample_info_df['Result'].str.contains('High Risk', case=False, na=False)
        ]
        
        for _, row in high_risk_df.iterrows():
            sample_id = row.get('sample_id', '')
            disease_str = row.get('Disease', '')
            
            if not sample_id or not disease_str:
                continue
            
            # Disease 컬럼에 여러 개가 있을 수 있음 (comma-separated)
            diseases = [d.strip() for d in disease_str.split(',')]
            
            for disease in diseases:
                # Disease가 이 염색체에 해당하는지 확인
                if disease in disease_chr_map and disease_chr_map[disease] == chr_name:
                    high_risk_samples.append((sample_id, disease))
                    break  # 같은 샘플을 중복으로 추가하지 않음
        
        return high_risk_samples
    
    def _create_sex_chromosome_files(self, output_dir: str, group: str):
        """Create male.txt and female.txt (기존 Positive 유지 + High Risk 추가)"""
        
        # Sex chromosome aneuploidy mapping
        male_diseases = ['XYY', 'XXY', 'XXYY']  # Male sex chromosome aneuploidies
        female_diseases = ['XO', 'XXX', 'X0']  # Female sex chromosome aneuploidies
        
        # === Female (XX) ===
        female_file = os.path.join(output_dir, "female.txt")
        
        # Load existing Positive samples (type != 'XX') - from source directory
        existing_female_positive = []
        source_female_file = os.path.join(self.config.source_ezd_dir, group, "female.txt")
        if os.path.exists(source_female_file):
            try:
                existing_df = pd.read_csv(source_female_file, sep='\t')
                # Normalize column name to lowercase 'sample'
                if 'Sample' in existing_df.columns:
                    existing_df.rename(columns={'Sample': 'sample'}, inplace=True)
                # Keep non-XX samples (these are positive/disease samples)
                existing_female_positive = existing_df[existing_df['type'] != 'XX'].to_dict('records')
                logger.info(f"  female.txt: Loaded {len(existing_female_positive)} existing Positive samples from source")
            except Exception as e:
                logger.debug(f"  female.txt: Could not load existing file: {e}")
        else:
            logger.debug(f"  female.txt: No source file found at {source_female_file}")
        
        # Normal samples (Low Risk females)
        female_rows = []
        for sample_id in self.samples['F']:
            uar, z = self._get_sample_chr_data(sample_id, "chrX", group)
            if uar is not None and z is not None:
                female_rows.append({
                    'sample': sample_id,
                    'type': 'XX',  # Normal female
                    'UAR': uar,
                    'Z': z
                })
        
        # Add High Risk samples with female sex chromosome aneuploidies
        if self.sample_info_df is not None:
            high_risk_df = self.sample_info_df[
                self.sample_info_df['Result'].str.contains('High Risk', case=False, na=False)
            ]
            
            for _, row in high_risk_df.iterrows():
                sample_id = row.get('sample_id', '')
                disease_str = row.get('Disease', '')
                
                if not sample_id or not disease_str:
                    continue
                
                diseases = [d.strip() for d in disease_str.split(',')]
                
                # Check if any disease is a female sex chromosome aneuploidy
                for disease in diseases:
                    if disease in female_diseases:
                        uar, z = self._get_sample_chr_data(sample_id, "chrX", group)
                        if uar is not None and z is not None:
                            existing_female_positive.append({
                                'sample': sample_id,
                                'type': disease,  # Disease type (XO, XXX, etc.)
                                'UAR': uar,
                                'Z': z
                            })
                        break  # Only add once per sample
        
        # Combine Normal + Positive
        all_female_rows = female_rows + existing_female_positive
        if all_female_rows:
            female_df = pd.DataFrame(all_female_rows)
            female_df = female_df[['sample', 'type', 'UAR', 'Z']]
            female_df.to_csv(female_file, sep='\t', index=False)
            logger.info(f"  female.txt: {len(female_rows)} Normal + {len(existing_female_positive)} Positive = {len(all_female_rows)} total")
        
        # === Male (XY) ===
        male_file = os.path.join(output_dir, "male.txt")
        
        # Load existing Positive samples (type != 'XY') - from source directory
        existing_male_positive = []
        source_male_file = os.path.join(self.config.source_ezd_dir, group, "male.txt")
        if os.path.exists(source_male_file):
            try:
                existing_df = pd.read_csv(source_male_file, sep='\t')
                # Normalize column name to lowercase 'sample'
                if 'Sample' in existing_df.columns:
                    existing_df.rename(columns={'Sample': 'sample'}, inplace=True)
                # Keep non-XY samples (these are positive/disease samples)
                existing_male_positive = existing_df[existing_df['type'] != 'XY'].to_dict('records')
                logger.info(f"  male.txt: Loaded {len(existing_male_positive)} existing Positive samples from source")
            except Exception as e:
                logger.debug(f"  male.txt: Could not load existing file: {e}")
        else:
            logger.debug(f"  male.txt: No source file found at {source_male_file}")
        
        # Normal samples (Low Risk males) - need both chrX and chrY
        male_rows = []
        for sample_id in self.samples['M']:
            uar_x, _ = self._get_sample_chr_data(sample_id, "chrX", group)
            uar_y, _ = self._get_sample_chr_data(sample_id, "chrY", group)
            if uar_x is not None and uar_y is not None:
                male_rows.append({
                    'sample': sample_id,
                    'type': 'XY',  # Normal male
                    'UAR.X': uar_x,
                    'UAR.Y': uar_y
                })
        
        # Add High Risk samples with male sex chromosome aneuploidies
        if self.sample_info_df is not None:
            high_risk_df = self.sample_info_df[
                self.sample_info_df['Result'].str.contains('High Risk', case=False, na=False)
            ]
            
            for _, row in high_risk_df.iterrows():
                sample_id = row.get('sample_id', '')
                disease_str = row.get('Disease', '')
                
                if not sample_id or not disease_str:
                    continue
                
                diseases = [d.strip() for d in disease_str.split(',')]
                
                # Check if any disease is a male sex chromosome aneuploidy
                for disease in diseases:
                    if disease in male_diseases:
                        uar_x, _ = self._get_sample_chr_data(sample_id, "chrX", group)
                        uar_y, _ = self._get_sample_chr_data(sample_id, "chrY", group)
                        if uar_x is not None and uar_y is not None:
                            existing_male_positive.append({
                                'sample': sample_id,
                                'type': disease,  # Disease type (XXY, XYY, etc.)
                                'UAR.X': uar_x,
                                'UAR.Y': uar_y
                            })
                        break  # Only add once per sample
        
        # Combine Normal + Positive
        all_male_rows = male_rows + existing_male_positive
        if all_male_rows:
            male_df = pd.DataFrame(all_male_rows)
            male_df = male_df[['sample', 'type', 'UAR.X', 'UAR.Y']]
            male_df.to_csv(male_file, sep='\t', index=False)
            logger.info(f"  male.txt: {len(male_rows)} Normal + {len(existing_male_positive)} Positive = {len(all_male_rows)} total")
    
    def _calculate_thresholds(self, chr_data: Dict[str, pd.DataFrame], output_dir: str, group: str):
        """
        Threshold 계산 및 저장
        
        Args:
            chr_data: 염색체별 데이터프레임 딕셔너리
            output_dir: 출력 디렉토리
            group: 'orig', 'fetus', 'mom'
        """
        threshold_rows = []
        
        # Sort chromosomes numerically (chr1, chr2, ..., chr22)
        def chr_sort_key(chr_name):
            # Extract number from 'chr1', 'chr2', etc.
            num_str = chr_name.replace('chr', '')
            try:
                return int(num_str)
            except ValueError:
                return 999  # Non-numeric chromosomes go to the end
        
        for chr_name in sorted(chr_data.keys(), key=chr_sort_key):
            df = chr_data[chr_name]
            
            # Normal과 Positive 샘플 분리
            normal_df = df[df['type'] == 'N']
            positive_df = df[df['type'] == 'P']
            
            if len(normal_df) == 0:
                continue
            
            # Normal 샘플의 최대값
            uar_max_n = normal_df['UAR'].max()
            z_max_n = normal_df['Z'].max()
            
            # Threshold 계산: N과 P가 있으면 중간값, 없으면 N 최대값 + 0.01
            if len(positive_df) > 0:
                # Positive 샘플의 최소값
                uar_min_p = positive_df['UAR'].min()
                z_min_p = positive_df['Z'].min()
                
                # N 최대값과 P 최소값의 중간값을 threshold로 사용
                uar_min = (uar_max_n + uar_min_p) / 2
                z_min = (z_max_n + z_min_p) / 2
            else:
                # Positive 샘플이 없으면 기존 방식 (N 최대값 + 0.01)
                uar_min = uar_max_n + 0.01
                z_min = z_max_n + 0.01
            
            # 상한값은 기존 값 유지 또는 새로 계산
            # 여기서는 간단히 percentile 사용
            uar_max_threshold = normal_df['UAR'].quantile(0.95) + 0.1
            z_max_threshold = normal_df['Z'].quantile(0.95) + 0.5
            
            threshold_rows.append({
                'chr': chr_name,
                'UAR_min': round(uar_min, 2),
                'UAR_max': round(uar_max_threshold, 2),
                'Z_min': round(z_min, 2),
                'Z_max': round(z_max_threshold, 2)
            })
        
        if threshold_rows:
            threshold_df = pd.DataFrame(threshold_rows)
            threshold_file = os.path.join(output_dir, f"{group}_thresholds_new.tsv")
            threshold_df.to_csv(threshold_file, sep='\t', index=False)
            logger.info(f"  Thresholds saved to: {threshold_file}")
    
    def _create_sca_config(self, output_dir: str, group: str):
        """
        SCA detection config 파일 생성
        
        Args:
            output_dir: 출력 디렉토리
            group: 'orig', 'fetus', 'mom'
        """
        logger.info("\n" + "="*60)
        logger.info(f"Creating SCA Config for {group}")
        logger.info("="*60)
        
        # Load male and female data
        male_file = os.path.join(output_dir, "male.txt")
        female_file = os.path.join(output_dir, "female.txt")
        
        sca_config = {
            "config_info": {
                "name": f"{group.capitalize()} Reference SCA Detection",
                "type": group,
                "description": f"{group.capitalize()} population reference data for SCA detection",
                "version": "1.0",
                "created_date": pd.Timestamp.now().strftime("%Y-%m-%d"),
                "last_updated": pd.Timestamp.now().strftime("%Y-%m-%d")
            },
            "sca_detection": {},
            "detection_rules": {},
            "plot_settings": {},
            "file_paths": {}
        }
        
        # === Male SCA Config ===
        if os.path.exists(male_file):
            male_config = self._calculate_male_sca_config(male_file)
            sca_config["sca_detection"]["male"] = male_config
            sca_config["detection_rules"]["male"] = {
                "logic": "boundary_line_with_threshold",
                "steps": [
                    "1. Calculate boundary_y = slope * ur_x + intercept",
                    "2. If ur_y > boundary_y + margin: Abnormal",
                    "3. If abnormal and ur_x >= ur_x_threshold: XYY Detected",
                    "4. If abnormal and ur_x < ur_x_threshold: XXY Detected",
                    "5. If ur_y > boundary_y (no margin): Suspected",
                    "6. Else: Normal"
                ]
            }
            sca_config["plot_settings"]["male"] = {
                "x_axis_range": [4.0, 6.0],
                "y_axis_range": [0.0, 0.15],
                "title": f"SCA (Male) - {group.capitalize()} Reference"
            }
        
        # === Female SCA Config ===
        if os.path.exists(female_file):
            female_config = self._calculate_female_sca_config(female_file)
            sca_config["sca_detection"]["female"] = female_config
            sca_config["detection_rules"]["female"] = {
                "logic": "z_score_and_ur_x_ranges",
                "steps": [
                    "1. If z_score < xo_z_threshold and ur_x in [xo_ur_x_min, xo_ur_x_max]: XO Detected",
                    "2. If z_score > xxx_z_threshold and ur_x in [xxx_ur_x_min, xxx_ur_x_max]: XXX Detected",
                    "3. If z_score in [z_normal_low, z_normal_high] and ur_x in [ur_x_low, ur_x_high]: Normal",
                    "4. Else: Not Detected"
                ]
            }
            sca_config["plot_settings"]["female"] = {
                "x_axis_range": [4.0, 7.0],
                "y_axis_range": [-40, 40],
                "title": f"SCA (Female) - {group.capitalize()} Reference"
            }
        
        # Plot settings
        sca_config["plot_settings"]["figure_size"] = [12, 8]
        sca_config["plot_settings"]["dpi"] = 100
        sca_config["plot_settings"]["colors"] = {
            "xy_normal": "blue",
            "xxy": "orange",
            "xyy": "green",
            "xx_normal": "gray",
            "xo": "coral",
            "xxx": "gold",
            "boundary_line": "purple",
            "threshold_line": "red",
            "test_sample": "red"
        }
        
        # File paths
        sca_config["file_paths"] = {
            "male_reference": male_file,
            "female_reference": female_file,
            "output_directory": output_dir
        }
        
        # Save config
        config_file = os.path.join(output_dir, "sca_config.json")
        with open(config_file, 'w', encoding='utf-8') as f:
            json.dump(sca_config, f, indent=4, ensure_ascii=False)
        
        logger.info(f"  SCA config saved to: {config_file}")
    
    def _calculate_male_sca_config(self, male_file: str) -> dict:
        """
        Male SCA detection config 계산 (XYY 대각선 라인)
        
        Args:
            male_file: male.txt 파일 경로
            
        Returns:
            Male SCA config dictionary
        """
        try:
            df = pd.read_csv(male_file, sep='\t')
            
            # Normalize column names
            df.columns = [col.replace('.', '_') for col in df.columns]
            
            # Extract XY (normal) data
            xy_data = df[df['type'] == 'XY'].copy()
            
            if len(xy_data) == 0:
                logger.warning("No XY (normal) data found in male.txt")
                return self._get_default_male_config()
            
            # Calculate boundary line using 95th percentile
            ur_x_col = 'UAR_X'
            ur_y_col = 'UAR_Y'
            
            # UAR_X 범위
            x_min, x_max = xy_data[ur_x_col].min(), xy_data[ur_x_col].max()
            
            # Bin으로 나누어 각 구간의 95th percentile 계산
            n_bins = min(10, len(xy_data) // 5)
            x_bins = np.linspace(x_min, x_max, n_bins + 1)
            
            boundary_points = []
            
            for i in range(len(x_bins) - 1):
                x_bin_min, x_bin_max = x_bins[i], x_bins[i + 1]
                x_center = (x_bin_min + x_bin_max) / 2
                
                bin_mask = ((xy_data[ur_x_col] >= x_bin_min) & 
                           (xy_data[ur_x_col] < x_bin_max))
                bin_data = xy_data[bin_mask]
                
                if len(bin_data) >= 3:
                    y_boundary = np.percentile(bin_data[ur_y_col], 95)
                    boundary_points.append((x_center, y_boundary))
            
            if len(boundary_points) < 2:
                logger.warning("Not enough data to calculate male boundary line")
                return self._get_default_male_config()
            
            # Linear regression using numpy
            X = np.array([p[0] for p in boundary_points])
            y = np.array([p[1] for p in boundary_points])
            
            # Calculate slope and intercept using numpy polyfit
            slope, intercept = np.polyfit(X, y, 1)
            
            # Calculate R² score
            y_pred = slope * X + intercept
            ss_res = np.sum((y - y_pred) ** 2)
            ss_tot = np.sum((y - np.mean(y)) ** 2)
            r2_score = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0
            
            slope = float(slope)
            intercept = float(intercept)
            
            logger.info(f"  Male boundary line: y = {slope:.6f}x + {intercept:.6f} (R²={r2_score:.4f})")
            
            # Calculate ur_x_threshold (median of XY data)
            ur_x_threshold = float(xy_data[ur_x_col].median())
            
            return {
                "slope": slope,
                "intercept": intercept,
                "ur_x_threshold": round(ur_x_threshold, 2),
                "margin": 0.005,
                "description": f"{len(xy_data)} XY samples used for boundary calculation",
                "sample_count": len(xy_data),
                "detection_types": ["XYY", "XXY", "Normal"],
                "r2_score": round(r2_score, 4)
            }
            
        except Exception as e:
            logger.error(f"Error calculating male SCA config: {e}")
            return self._get_default_male_config()
    
    def _calculate_female_sca_config(self, female_file: str) -> dict:
        """
        Female SCA detection config 계산 (XO/XXX thresholds)
        
        Args:
            female_file: female.txt 파일 경로
            
        Returns:
            Female SCA config dictionary
        """
        try:
            df = pd.read_csv(female_file, sep='\t')
            
            # XX (normal) data
            xx_data = df[df['type'] == 'XX'].copy()
            
            # XO and XXX positive data
            xo_data = df[df['type'].isin(['XO', 'X0'])].copy()
            xxx_data = df[df['type'] == 'XXX'].copy()
            
            if len(xx_data) == 0:
                logger.warning("No XX (normal) data found in female.txt")
                return self._get_default_female_config()
            
            # Calculate thresholds based on XX data distribution
            uar_median = xx_data['UAR'].median()
            uar_std = xx_data['UAR'].std()
            z_median = xx_data['Z'].median()
            z_std = xx_data['Z'].std()
            
            # UAR ranges
            ur_x_low = float(xx_data['UAR'].quantile(0.05))
            ur_x_high = float(xx_data['UAR'].quantile(0.95))
            
            # Z-score ranges
            z_normal_low = float(xx_data['Z'].quantile(0.05))
            z_normal_high = float(xx_data['Z'].quantile(0.95))
            
            # XO threshold: XX와 XO 사이 (median 사용으로 outlier에 robust)
            if len(xo_data) > 0:
                xo_z_median = xo_data['Z'].median()
                xx_z_min = xx_data['Z'].min()
                # XX 최소값과 XO median 사이
                xo_z_threshold = float((xx_z_min + xo_z_median) / 2)
                logger.info(f"  XO data: {len(xo_data)} samples, Z_median={xo_z_median:.2f}")
            else:
                xo_z_threshold = float(z_normal_low - 2 * z_std)
                logger.info(f"  No XO data, using calculated threshold")
            
            xo_ur_x_min = float(ur_x_low - 0.5)
            xo_ur_x_max = float(ur_x_low)
            
            # XXX threshold: XX와 XXX 사이 (median 사용으로 outlier에 robust)
            if len(xxx_data) > 0:
                xxx_z_median = xxx_data['Z'].median()
                xx_z_max = xx_data['Z'].max()
                # XX 최대값과 XXX median 사이
                xxx_z_threshold = float((xx_z_max + xxx_z_median) / 2)
                logger.info(f"  XXX data: {len(xxx_data)} samples, Z_median={xxx_z_median:.2f}")
            else:
                xxx_z_threshold = float(z_normal_high + 2 * z_std)
                logger.info(f"  No XXX data, using calculated threshold")
            
            xxx_ur_x_min = float(ur_x_high)
            xxx_ur_x_max = float(ur_x_high + 0.5)
            
            logger.info(f"  Female XX data: UAR={uar_median:.2f}±{uar_std:.2f}, Z={z_median:.2f}±{z_std:.2f}")
            logger.info(f"  XO threshold: Z < {xo_z_threshold:.2f}")
            logger.info(f"  XXX threshold: Z > {xxx_z_threshold:.2f}")
            
            return {
                "xo_z_threshold": round(xo_z_threshold, 2),
                "xxx_z_threshold": round(xxx_z_threshold, 2),
                "ur_x_low": round(ur_x_low, 2),
                "ur_x_high": round(ur_x_high, 2),
                "z_normal_low": round(z_normal_low, 2),
                "z_normal_high": round(z_normal_high, 2),
                "xo_ur_x_min": round(xo_ur_x_min, 2),
                "xo_ur_x_max": round(xo_ur_x_max, 2),
                "xxx_ur_x_min": round(xxx_ur_x_min, 2),
                "xxx_ur_x_max": round(xxx_ur_x_max, 2),
                "description": f"{len(xx_data)} XX samples used for threshold calculation",
                "sample_count": len(xx_data),
                "detection_types": ["XO", "XXX", "Normal"]
            }
            
        except Exception as e:
            logger.error(f"Error calculating female SCA config: {e}")
            return self._get_default_female_config()
    
    def _get_default_male_config(self) -> dict:
        """Default male SCA config (fallback)"""
        return {
            "slope": -0.0837,
            "intercept": 0.4877,
            "ur_x_threshold": 5.3,
            "margin": 0.005,
            "description": "Default configuration (no data available)",
            "sample_count": 0,
            "detection_types": ["XYY", "XXY", "Normal"]
        }
    
    def _get_default_female_config(self) -> dict:
        """Default female SCA config (fallback)"""
        return {
            "xo_z_threshold": -5.5,
            "xxx_z_threshold": 4.5,
            "ur_x_low": 5.25,
            "ur_x_high": 5.49,
            "z_normal_low": -3.0,
            "z_normal_high": 1.0,
            "xo_ur_x_min": 4.9,
            "xo_ur_x_max": 5.2,
            "xxx_ur_x_min": 5.52,
            "xxx_ur_x_max": 6.0,
            "description": "Default configuration (no data available)",
            "sample_count": 0,
            "detection_types": ["XO", "XXX", "Normal"]
        }


class PRIZMReferenceCreator:
    """PRIZM Reference 생성 클래스"""
    
    def __init__(self, config: ReferenceConfig, samples: Dict[str, List[str]], sample_df: pd.DataFrame = None):
        self.config = config
        self.samples = samples
        self.sample_df = sample_df
    
    def create_reference(self, group: str):
        """
        PRIZM reference 생성
        
        Args:
            group: 'orig', 'fetus', 'mom'
        """
        logger.info(f"\n{'='*60}")
        logger.info(f"Creating PRIZM Reference for {group}")
        logger.info(f"{'='*60}")
        
        output_dir = os.path.join(self.config.prizm_dir, group)
        
        # 10mb count 파일 수집
        male_count_files = self._collect_count_files(self.samples['M'], group)
        female_count_files = self._collect_count_files(self.samples['F'], group)
        all_count_files = male_count_files + female_count_files
        
        if len(all_count_files) == 0:
            logger.warning(f"No count files found for PRIZM {group}")
            return
        
        logger.info(f"Total count files: {len(all_count_files)} (M:{len(male_count_files)}, F:{len(female_count_files)})")
        
        # Mean & SD 계산 (chromosome level)
        self._calculate_mean_sd(all_count_files, output_dir, 'total', autosomal=True)
        self._calculate_mean_sd(male_count_files, output_dir, 'male', autosomal=False)
        self._calculate_mean_sd(female_count_files, output_dir, 'female', autosomal=False)
        
        # 10mb level
        self._calculate_10mb_mean_sd(all_count_files, output_dir, 'total', autosomal=True)
        self._calculate_10mb_mean_sd(male_count_files, output_dir, 'male', autosomal=False)
        self._calculate_10mb_mean_sd(female_count_files, output_dir, 'female', autosomal=False)
        
        # 10mb_all level (autosomal + sex chromosomes combined)
        self._calculate_10mb_all_combined(all_count_files, male_count_files, female_count_files, output_dir)
        
        logger.info(f"PRIZM reference for {group} created successfully")
    
    def _collect_count_files(self, sample_ids: List[str], group: str) -> List[str]:
        """
        10mb count 파일 수집
        
        Args:
            sample_ids: 샘플 ID 리스트
            group: 'orig', 'fetus', 'mom'
        
        Returns:
            파일 경로 리스트
        """
        count_files = []
        
        # sample_df가 있으면 효율적으로 파일 찾기
        if self.sample_df is not None:
            sample_subset = self.sample_df[self.sample_df['sample_id'].isin(sample_ids)]
            
            for _, row in sample_subset.iterrows():
                sample_id = row['sample_id']
                sample_dir = row['sample_dir']
                
                # group에 따른 파일 이름 결정
                if group == 'orig':
                    filename = f"{sample_id}.of_orig.10mb.wig.Normalization.txt"
                elif group == 'fetus':
                    filename = f"{sample_id}.of_fetus.10mb.wig.Normalization.txt"
                elif group == 'mom':
                    filename = f"{sample_id}.of_mom.10mb.wig.Normalization.txt"
                else:
                    continue
                
                # Output_hmmcopy 디렉토리 안에 파일이 있음
                file_path = os.path.join(sample_dir, "Output_hmmcopy", filename)
                
                if os.path.exists(file_path):
                    count_files.append(file_path)
        else:
            # 기존 방식 (느림)
            for sample_id in sample_ids:
                if group == 'orig':
                    filename = f"{sample_id}.of_orig.10mb.wig.Normalization.txt"
                elif group == 'fetus':
                    filename = f"{sample_id}.of_fetus.10mb.wig.Normalization.txt"
                elif group == 'mom':
                    filename = f"{sample_id}.of_mom.10mb.wig.Normalization.txt"
                else:
                    continue
                
                # 경로 찾기 (여러 월별 디렉토리 확인)
                file_path = self._find_count_file(sample_id, filename)
                
                if file_path and os.path.exists(file_path):
                    count_files.append(file_path)
        
        return count_files
    
    def _find_count_file(self, sample_id: str, filename: str) -> str:
        """샘플의 count 파일 찾기 (느린 fallback 방법)"""
        # 분석 디렉토리에서 샘플 찾기
        # analysis/YYMM/sample_id/Output_hmmcopy/ 형태
        
        for root, dirs, files in os.walk(ANALYSIS_DIR):
            if sample_id in root and filename in files:
                return os.path.join(root, filename)
        
        return None
    
    def _calculate_mean_sd(self, count_files: List[str], output_dir: str, prefix: str, autosomal: bool):
        """Chromosome level Mean & SD 계산"""
        import re
        
        logger.info(f"  Calculating {prefix} chromosome-level mean & sd...")
        
        mean_file = os.path.join(output_dir, f"{prefix}_mean.csv")
        sd_file = os.path.join(output_dir, f"{prefix}_sd.csv")
        
        def natural_keys(text):
            def atoi(text):
                return int(text) if text.isdigit() else text
            return [atoi(c) for c in re.split(r'(\d+)', text)]
        
        # 각 파일별로 normalized matrix 계산
        dfs = {}
        for idx, filename in enumerate(count_files):
            try:
                count_data = pd.read_csv(filename, sep='\t')
                count_data = count_data[['chr', 'start', 'reads']].rename(
                    columns={'start': 'bin', 'reads': 'count'}
                )
                
                count_data['bin'] = count_data['bin'].apply(lambda x: x/10000000)
                count_data_allsum = count_data.groupby(by=['chr'])['count'].sum()
                count_data_df = count_data_allsum.to_frame()
                sorted_index = sorted(count_data_df.index, key=natural_keys)
                count_data_df = count_data_df.loc[sorted_index]
                
                if autosomal:
                    count_data_df = count_data_df.iloc[0:22, :]
                
                count_sum = sum(count_data_df.loc[:, 'count'])
                norm_data = [
                    (count_data_df.iloc[i, count_data_df.columns.get_loc('count')] / float(count_sum)) * 100
                    for i in range(len(count_data_df))
                ]
                
                count_data_df.insert(1, 'ratio', norm_data)
                normalized_dict = {
                    j: [
                        count_data_df.iloc[i, count_data_df.columns.get_loc('ratio')] / 
                        count_data_df.iloc[j, count_data_df.columns.get_loc('ratio')]
                        for i in range(len(count_data_df))
                    ]
                    for j in range(len(count_data_df))
                }
                
                dfs[idx] = pd.DataFrame.from_dict(normalized_dict)
            except Exception as e:
                logger.warning(f"    Failed to process {filename}: {e}")
                continue
        
        if not dfs:
            logger.error(f"    No valid count files processed for {prefix}")
            return
        
        # Mean & SD 계산
        all_dfs = pd.concat(dfs, axis=0)
        mean_df = all_dfs.groupby(level=1).mean()
        sd_df = all_dfs.groupby(level=1).std()
        
        # 파일 저장
        mean_df.to_csv(mean_file, sep='\t', index=False, header=False)
        sd_df.to_csv(sd_file, sep='\t', index=False, header=False)
        
        logger.info(f"    Saved: {mean_file}")
        logger.info(f"    Saved: {sd_file}")
    
    def _calculate_10mb_mean_sd(self, count_files: List[str], output_dir: str, prefix: str, autosomal: bool):
        """10mb bin level Mean & SD 계산"""
        import re
        
        logger.info(f"  Calculating {prefix} 10mb-level mean & sd...")
        
        mean_10mb_file = os.path.join(output_dir, f"{prefix}_10mb_mean.csv")
        sd_10mb_file = os.path.join(output_dir, f"{prefix}_10mb_sd.csv")
        
        def natural_keys(text):
            def atoi(text):
                return int(text) if text.isdigit() else text
            return [atoi(c) for c in re.split(r'(\d+)', text)]
        
        dfs = {}
        for idx, filename in enumerate(count_files):
            try:
                count_data = pd.read_csv(filename, sep='\t')
                count_data = count_data[['chr', 'start', 'reads']].rename(
                    columns={'start': 'bin', 'reads': 'count'}
                )
                
                count_data['bin'] = count_data['bin'].apply(lambda x: x/10000000)
                count_data_allsum = count_data.groupby(by=['chr'])['count'].sum()
                count_data_df = count_data_allsum.to_frame()
                count_data_df.insert(0, 'bin', 0)
                sorted_index = sorted(count_data_df.index, key=natural_keys)
                
                # Exclude specific chromosomes from aggregation
                excluded_chromosomes = ['chr9', 'chr13', 'chr16', 'chr18', 'chr21', 'chr22']
                if not autosomal:
                    excluded_chromosomes.extend(['chrX', 'chrY'])
                
                count_data_df = count_data_df.drop(excluded_chromosomes, errors='ignore')
                count_data_10mb = count_data.set_index('chr').loc[
                    [c for c in excluded_chromosomes if c in count_data.set_index('chr').index]
                ]
                
                count_data_mix = pd.concat([count_data_df, count_data_10mb])
                count_data_mix = count_data_mix.loc[sorted_index[:-2] if autosomal else sorted_index]
                
                count_sum = count_data_mix.loc[:, 'count'].sum()
                norm_data = [
                    (count_data_mix.iloc[i, count_data_mix.columns.get_loc('count')] / float(count_sum)) * 100
                    for i in range(len(count_data_mix))
                ]
                count_data_mix.insert(2, 'ratio', norm_data)
                
                normalized_dict = {
                    j: [
                        count_data_mix.iloc[i, count_data_mix.columns.get_loc('ratio')] /
                        count_data_mix.iloc[j, count_data_mix.columns.get_loc('ratio')]
                        if count_data_mix.iloc[i, count_data_mix.columns.get_loc('ratio')] != 0.0
                        and count_data_mix.iloc[j, count_data_mix.columns.get_loc('ratio')] != 0.0
                        else 0
                        for i in range(len(count_data_mix))
                    ]
                    for j in range(len(count_data_mix))
                }
                
                dfs[idx] = pd.DataFrame.from_dict(normalized_dict)
            except Exception as e:
                logger.warning(f"    Failed to process {filename}: {e}")
                continue
        
        if not dfs:
            logger.error(f"    No valid count files processed for {prefix} 10mb")
            return
        
        all_10mb_dfs = pd.concat(dfs, axis=0)
        mean_10mb_df = all_10mb_dfs.groupby(level=1).mean()
        sd_10mb_df = all_10mb_dfs.groupby(level=1).std()
        
        mean_10mb_df.to_csv(mean_10mb_file, sep='\t', index=False, header=False)
        sd_10mb_df.to_csv(sd_10mb_file, sep='\t', index=False, header=False)
        
        logger.info(f"    Saved: {mean_10mb_file}")
        logger.info(f"    Saved: {sd_10mb_file}")
    
    def _calculate_10mb_all_mean_sd(self, count_files: List[str], output_dir: str, prefix: str, autosomal: bool):
        """10mb_all level Mean & SD 계산"""
        import re
        
        logger.info(f"  Calculating {prefix} 10mb_all-level mean & sd...")
        
        mean_10mb_all_file = os.path.join(output_dir, f"{prefix}_10mb_all_mean.csv")
        sd_10mb_all_file = os.path.join(output_dir, f"{prefix}_10mb_all_sd.csv")
        
        def natural_keys(text):
            def atoi(text):
                return int(text) if text.isdigit() else text
            return [atoi(c) for c in re.split(r'(\d+)', text)]
        
        def getNormalized10mbData(filename, autosomal, key_start=0):
            df = pd.read_csv(filename, sep='\t')
            df = df[['chr', 'start', 'reads']].rename(columns={'start': 'bin', 'reads': 'count'})
            df['bin'] = df['bin'].apply(lambda x: (x-1)/10000000)
            df.bin = df.bin.astype(int)
            
            # Sort by chromosome order
            df = df.set_index('chr').loc[
                ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
                 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
                 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
            ]
            df.reset_index(level=0, inplace=True)
            
            count_data_allsum = df.groupby(by=['chr'])['count'].sum()
            sum_df = count_data_allsum.to_frame()
            sorted_index = sorted(sum_df.index, key=natural_keys)
            sum_df = sum_df.loc[sorted_index]
            
            sum_df = sum_df.iloc[0:22, :]
            count_sum = sum_df['count'].sum()
            
            norm_data = [
                (sum_df.loc[sum_df.index[i], 'count']/float(count_sum))*100
                for i in range(len(sum_df))
            ]
            sum_df.insert(1, 'ratio', norm_data)
            
            if autosomal == True:
                df = df.loc[~df['chr'].isin(['chrX', 'chrY'])]
            else:
                df = df.loc[df['chr'].isin(['chrX', 'chrY'])]
            
            norm_10mb_data = [
                (df.loc[df.index[i], 'count']/float(count_sum))*100
                for i in range(len(df))
            ]
            df.insert(2, 'ratio', norm_10mb_data)
            
            normalized_dict = {
                j + key_start: [
                    df.loc[df.index[j], 'ratio'] / sum_df.loc[sum_df.index[i], 'ratio']
                    if sum_df.loc[sum_df.index[i], 'ratio'] != 0.0
                    and df.loc[df.index[j], 'ratio'] != 0.0
                    else 0
                    for i in range(len(sum_df))
                ]
                for j in range(len(df))
            }
            
            return normalized_dict
        
        dfs = {}
        for idx, filename in enumerate(count_files):
            try:
                normalized_dict = getNormalized10mbData(filename, autosomal)
                dfs[idx] = pd.DataFrame.from_dict(normalized_dict)
            except Exception as e:
                logger.warning(f"    Failed to process {filename}: {e}")
                continue
        
        if not dfs:
            logger.error(f"    No valid count files processed for {prefix} 10mb_all")
            return
        
        all_dfs = pd.concat(dfs, axis=0, keys=range(len(dfs)))
        mean_df = all_dfs.groupby(level=1).mean()
        sd_df = all_dfs.groupby(level=1).std()
        
        mean_df.to_csv(mean_10mb_all_file, sep='\t', index=False, header=False)
        sd_df.to_csv(sd_10mb_all_file, sep='\t', index=False, header=False)
        
        logger.info(f"    Saved: {mean_10mb_all_file}")
        logger.info(f"    Saved: {sd_10mb_all_file}")
    
    def _calculate_10mb_all_combined(self, all_count_files: List[str], male_count_files: List[str], 
                                     female_count_files: List[str], output_dir: str):
        """
        10mb_all level Mean & SD 계산 (autosomal + sex chromosomes combined)
        
        Generate 10mb_all references by combining:
        - Autosomal bins (from all files)
        - Sex chromosome bins (from male/female files separately)
        """
        import re
        
        def natural_keys(text):
            def atoi(text):
                return int(text) if text.isdigit() else text
            return [atoi(c) for c in re.split(r'(\d+)', text)]
        
        def getNormalized10mbData(filename, autosomal, key_start=0):
            df = pd.read_csv(filename, sep='\t')
            df = df[['chr', 'start', 'reads']].rename(columns={'start': 'bin', 'reads': 'count'})
            df['bin'] = df['bin'].apply(lambda x: (x-1)/10000000)
            df.bin = df.bin.astype(int)
            
            # Sort by chromosome order
            df = df.set_index('chr').loc[
                ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
                 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
                 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
            ]
            df.reset_index(level=0, inplace=True)
            
            count_data_allsum = df.groupby(by=['chr'])['count'].sum()
            sum_df = count_data_allsum.to_frame()
            sorted_index = sorted(sum_df.index, key=natural_keys)
            sum_df = sum_df.loc[sorted_index]
            
            sum_df = sum_df.iloc[0:22, :]
            count_sum = sum_df['count'].sum()
            
            norm_data = [
                (sum_df.loc[sum_df.index[i], 'count']/float(count_sum))*100
                for i in range(len(sum_df))
            ]
            sum_df.insert(1, 'ratio', norm_data)
            
            if autosomal == True:
                df = df.loc[~df['chr'].isin(['chrX', 'chrY'])]
            else:
                df = df.loc[df['chr'].isin(['chrX', 'chrY'])]
            
            norm_10mb_data = [
                (df.loc[df.index[i], 'count']/float(count_sum))*100
                for i in range(len(df))
            ]
            df.insert(2, 'ratio', norm_10mb_data)
            
            normalized_dict = {
                j + key_start: [
                    df.loc[df.index[j], 'ratio'] / sum_df.loc[sum_df.index[i], 'ratio']
                    if sum_df.loc[sum_df.index[i], 'ratio'] != 0.0
                    and df.loc[df.index[j], 'ratio'] != 0.0
                    else 0
                    for i in range(len(sum_df))
                ]
                for j in range(len(df))
            }
            
            return normalized_dict
        
        def makeReference_10mball(file_list, autosomal):
            dfs = {}
            for idx, filename in enumerate(file_list):
                try:
                    normalized_dict = getNormalized10mbData(filename, autosomal)
                    dfs[idx] = pd.DataFrame.from_dict(normalized_dict)
                except Exception as e:
                    logger.warning(f"    Failed to process {filename}: {e}")
                    continue
            
            if not dfs:
                return None, None
            
            all_dfs = pd.concat(dfs, axis=0, keys=range(len(dfs)))
            mean_df = all_dfs.groupby(level=1).mean()
            sd_df = all_dfs.groupby(level=1).std()
            
            return mean_df, sd_df
        
        # 1. Calculate autosomal bins (from all files)
        logger.info(f"  Calculating total 10mb_all (autosomal)...")
        total_mean_df, total_sd_df = makeReference_10mball(all_count_files, autosomal=True)
        
        if total_mean_df is None:
            logger.error("  Failed to calculate autosomal 10mb_all reference")
            return
        
        # 2. Calculate male sex chromosome bins
        logger.info(f"  Calculating male 10mb_all (sex chromosomes)...")
        male_sex_mean_df, male_sex_sd_df = makeReference_10mball(male_count_files, autosomal=False)
        
        # 3. Calculate female sex chromosome bins
        logger.info(f"  Calculating female 10mb_all (sex chromosomes)...")
        female_sex_mean_df, female_sex_sd_df = makeReference_10mball(female_count_files, autosomal=False)
        
        # 4. Combine autosomal + sex chromosomes for male
        if male_sex_mean_df is not None:
            male_mean_10mball_df = pd.concat([total_mean_df, male_sex_mean_df], axis=1)
            male_sd_10mball_df = pd.concat([total_sd_df, male_sex_sd_df], axis=1)
            
            male_mean_file = os.path.join(output_dir, 'male_10mb_all_mean.csv')
            male_sd_file = os.path.join(output_dir, 'male_10mb_all_sd.csv')
            
            male_mean_10mball_df.to_csv(male_mean_file, sep='\t', index=False, header=False)
            male_sd_10mball_df.to_csv(male_sd_file, sep='\t', index=False, header=False)
            
            logger.info(f"    Saved: {male_mean_file} (shape: {male_mean_10mball_df.shape})")
            logger.info(f"    Saved: {male_sd_file}")
        
        # 5. Combine autosomal + sex chromosomes for female
        if female_sex_mean_df is not None:
            female_mean_10mball_df = pd.concat([total_mean_df, female_sex_mean_df], axis=1)
            female_sd_10mball_df = pd.concat([total_sd_df, female_sex_sd_df], axis=1)
            
            female_mean_file = os.path.join(output_dir, 'female_10mb_all_mean.csv')
            female_sd_file = os.path.join(output_dir, 'female_10mb_all_sd.csv')
            
            female_mean_10mball_df.to_csv(female_mean_file, sep='\t', index=False, header=False)
            female_sd_10mball_df.to_csv(female_sd_file, sep='\t', index=False, header=False)
            
            logger.info(f"    Saved: {female_mean_file} (shape: {female_mean_10mball_df.shape})")
            logger.info(f"    Saved: {female_sd_file}")


class WCReferenceCreator:
    """Wisecondor Reference 생성 클래스 (Docker 기반)"""
    
    def __init__(self, config: ReferenceConfig, samples: Dict[str, List[str]], sample_df: pd.DataFrame = None):
        self.config = config
        self.samples = samples
        self.sample_df = sample_df
        self.docker_image = "nipt_docker_v1.3"
        self.wc_docker_script = os.path.join(SCRIPT_DIR, "create_wc_docker.py")
    
    def create_reference(self, group: str, binsize: int = 200000):
        """
        Wisecondor reference 생성 (Docker 기반)
        
        Args:
            group: 'orig', 'fetus', 'mom'
            binsize: bin 크기 (기본 200kb)
        """
        logger.info(f"\n{'='*60}")
        logger.info(f"Creating WC Reference for {group} (Docker)")
        logger.info(f"{'='*60}")
        
        # Sample list를 임시 필터링된 파일로 준비 (이미 필터링된 경우)
        # Docker 스크립트가 sample_df를 사용하므로 기존 로직 활용
        
        # Docker 명령 실행
        output_dir_docker = f"/refs/{os.path.basename(self.config.output_base_dir)}/WC"
        
        # 필터링된 샘플 리스트 경로 (create_reference.py에서 이미 생성됨)
        # sample_list_file을 찾거나 임시로 생성
        sample_list_docker = "/refs/ucl/reference_make/reference_sample_list_UCL_filtered.tsv"
        
        cmd = [
            "docker", "run", "--rm",
            "--entrypoint", "bash",
            "-v", f"{ANALYSIS_DIR}:/analysis:ro",
            "-v", f"{DATA_DIR}/refs:/refs",
            "-v", f"{SCRIPT_DIR}:/scripts:ro",
            self.docker_image,
            "-c",
            f"python3 /scripts/create_wc_docker.py {sample_list_docker} {group} {output_dir_docker}"
        ]
        
        logger.info(f"Running Wisecondor newref in Docker...")
        logger.info(f"Docker command: {' '.join(cmd[:10])}...")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            logger.info(f"WC reference created for {group}")
            logger.info(result.stdout)
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to create WC reference: {e}")
            logger.error(f"stdout: {e.stdout}")
            logger.error(f"stderr: {e.stderr}")
    
    def _collect_npz_files(self, sample_ids: List[str], group: str, tool: str) -> List[str]:
        """NPZ 파일 수집"""
        npz_files = []
        
        # sample_df가 있으면 효율적으로 파일 찾기
        if self.sample_df is not None:
            sample_subset = self.sample_df[self.sample_df['sample_id'].isin(sample_ids)]
            
            for _, row in sample_subset.iterrows():
                sample_id = row['sample_id']
                sample_dir = row['sample_dir']
                
                # tool과 group에 따른 파일 경로 구성
                if tool.upper() == 'WC':
                    if group == 'orig':
                        # orig는 proper_paired를 사용
                        filename = f"{sample_id}.wc.of_orig.npz"
                    else:
                        filename = f"{sample_id}.wc.of_{group}.npz"
                    output_subdir = "Output_WC"
                elif tool.upper() == 'WCX':
                    if group == 'orig':
                        filename = f"{sample_id}.wcx.of_orig.npz"
                    else:
                        filename = f"{sample_id}.wcx.of_{group}.npz"
                    output_subdir = "Output_WCX"
                else:
                    continue
                
                file_path = os.path.join(sample_dir, output_subdir, filename)
                
                if os.path.exists(file_path):
                    npz_files.append(file_path)
        else:
            # 기존 방식 (느림)
            for sample_id in sample_ids:
                npz_file = self._find_npz_file(sample_id, group, tool)
                if npz_file and os.path.exists(npz_file):
                    npz_files.append(npz_file)
        
        return npz_files
    
    def _find_npz_file(self, sample_id: str, group: str, tool: str) -> str:
        """샘플의 NPZ 파일 찾기 (느린 fallback 방법)"""
        for root, dirs, files in os.walk(ANALYSIS_DIR):
            if tool.lower() in root.lower() and sample_id in root:
                for file in files:
                    if file.endswith(".npz") and group in file:
                        return os.path.join(root, file)
        
        return None


class WCXReferenceCreator:
    """WisecondorX Reference 생성 클래스 (Docker 기반)"""
    
    def __init__(self, config: ReferenceConfig, samples: Dict[str, List[str]], sample_df: pd.DataFrame = None):
        self.config = config
        self.samples = samples
        self.sample_df = sample_df
        self.docker_image = "nipt_docker_v1.3"
        self.wcx_docker_script = os.path.join(SCRIPT_DIR, "create_wcx_docker.py")
    
    def create_reference(self, group: str, binsize: int = 200000):
        """
        WisecondorX reference 생성 (Docker 기반)
        
        Args:
            group: 'orig', 'fetus', 'mom'
            binsize: bin 크기 (기본 200kb)
            
        Note:
            Docker 스크립트가 다음을 자동 생성:
            - orig/fetus: M, F, combined (3개)
            - mom: combined만 (1개)
        """
        logger.info(f"\n{'='*60}")
        logger.info(f"Creating WCX Reference for {group} (Docker)")
        logger.info(f"{'='*60}")
        
        # Docker 명령 실행
        output_dir_docker = f"/refs/{os.path.basename(self.config.output_base_dir)}/WCX"
        
        # 필터링된 샘플 리스트 경로
        sample_list_docker = "/refs/ucl/reference_make/reference_sample_list_UCL_filtered.tsv"
        
        cmd = [
            "docker", "run", "--rm",
            "--entrypoint", "bash",
            "-v", f"{ANALYSIS_DIR}:/analysis:ro",
            "-v", f"{DATA_DIR}/refs:/refs",
            "-v", f"{SCRIPT_DIR}:/scripts:ro",
            self.docker_image,
            "-c",
            f"python3 /scripts/create_wcx_docker.py {sample_list_docker} {group} {output_dir_docker}"
        ]
        
        logger.info(f"Running WisecondorX newref in Docker...")
        logger.info(f"Docker command: {' '.join(cmd[:10])}...")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            logger.info(f"WCX reference created for {group}")
            if group == 'mom':
                logger.info(f"  - {group}_200k_of.npz (combined)")
            else:
                logger.info(f"  - {group}_M_200k_{'proper_paired' if group=='orig' else 'of'}.npz")
                logger.info(f"  - {group}_F_200k_{'proper_paired' if group=='orig' else 'of'}.npz")
                logger.info(f"  - {group}_200k_{'proper_paired' if group=='orig' else 'of'}.npz (combined)")
            logger.info(result.stdout)
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to create WCX reference: {e}")
            logger.error(f"stdout: {e.stdout}")
            logger.error(f"stderr: {e.stderr}")
    
    def _collect_npz_files(self, sample_ids: List[str], group: str, tool: str) -> List[str]:
        """NPZ 파일 수집"""
        npz_files = []
        
        # sample_df가 있으면 효율적으로 파일 찾기
        if self.sample_df is not None:
            sample_subset = self.sample_df[self.sample_df['sample_id'].isin(sample_ids)]
            
            for _, row in sample_subset.iterrows():
                sample_id = row['sample_id']
                sample_dir = row['sample_dir']
                
                # tool과 group에 따른 파일 경로 구성
                if tool.upper() == 'WC':
                    if group == 'orig':
                        filename = f"{sample_id}.wc.of_orig.npz"
                    else:
                        filename = f"{sample_id}.wc.of_{group}.npz"
                    output_subdir = "Output_WC"
                elif tool.upper() == 'WCX':
                    if group == 'orig':
                        filename = f"{sample_id}.wcx.of_orig.npz"
                    else:
                        filename = f"{sample_id}.wcx.of_{group}.npz"
                    output_subdir = "Output_WCX"
                else:
                    continue
                
                file_path = os.path.join(sample_dir, output_subdir, filename)
                
                if os.path.exists(file_path):
                    npz_files.append(file_path)
        else:
            # 기존 방식 (느림)
            for sample_id in sample_ids:
                npz_file = self._find_npz_file(sample_id, group, tool)
                if npz_file and os.path.exists(npz_file):
                    npz_files.append(npz_file)
        
        return npz_files
    
    def _find_npz_file(self, sample_id: str, group: str, tool: str) -> str:
        """샘플의 NPZ 파일 찾기 (느린 fallback 방법)"""
        for root, dirs, files in os.walk(ANALYSIS_DIR):
            if tool.lower() in root.lower() and sample_id in root:
                for file in files:
                    if file.endswith(".npz") and group in file:
                        return os.path.join(root, file)
        
        return None


def main():
    parser = argparse.ArgumentParser(
        description="NIPT Reference Creation Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate sample list and create all references (one step)
  python create_reference.py --analysis-dirs /path/2507 /path/2508 --prefix GNCI --labcode ucl --ref-type all
  
  # Generate sample list only
  python create_reference.py --generate-sample-list --analysis-dirs /path/2507 /path/2508 --prefix GNCI --sample-list output.tsv
  
  # Create references from existing sample list
  python create_reference.py --sample-list samples.tsv --labcode ucl --ref-type all
  
  # Create only WC and WCX references (multiple types at once)
  python create_reference.py --sample-list samples.tsv --labcode ucl --ref-type wc wcx
  
  # Create only EZD reference
  python create_reference.py --sample-list samples.tsv --labcode ucl --ref-type ezd
  
  # Create for specific groups
  python create_reference.py --sample-list samples.tsv --labcode ucl --ref-type prizm --groups orig fetus
  
  # Preview sample filtering
  python create_reference.py --sample-list samples.tsv --preview-only
        """
    )
    
    # Sample list generation options
    parser.add_argument(
        '--generate-sample-list',
        action='store_true',
        help='Generate sample list from analysis directories (requires --analysis-dirs)'
    )
    
    parser.add_argument(
        '--analysis-dirs',
        nargs='+',
        help='Analysis directories to scan for samples (e.g., /path/2507 /path/2508)'
    )
    
    parser.add_argument(
        '--prefix',
        default='GNMF',
        help='Sample ID prefix for directory scanning (default: GNMF)'
    )
    
    parser.add_argument(
        '--output-root',
        help='Root directory containing report JSONs (default: WORKSPACE_ROOT/output)'
    )
    
    # Sample list input
    parser.add_argument(
        '--sample-list',
        help='Sample list TSV file (generated or from Excel). If --analysis-dirs provided, will be auto-generated.'
    )
    
    parser.add_argument(
        '--labcode',
        default='ucl',
        help='Laboratory code (default: ucl)'
    )
    
    parser.add_argument(
        '--ref-type',
        nargs='+',
        choices=['all', 'ezd', 'prizm', 'wc', 'wcx'],
        default=['all'],
        help='Reference type(s) to create (default: all). Can specify multiple: --ref-type wc wcx'
    )
    
    parser.add_argument(
        '--groups',
        nargs='+',
        choices=['orig', 'fetus', 'mom'],
        default=['orig', 'fetus', 'mom'],
        help='Groups to process (default: all)'
    )
    
    parser.add_argument(
        '--output-dir',
        help='Output directory (default: DATA_DIR/refs/LABCODE)'
    )
    
    parser.add_argument(
        '--reference-source',
        help='Source directory for existing Positive samples (default: DATA_DIR/refs/LABCODE)'
    )
    
    parser.add_argument(
        '--min-seqff',
        type=float,
        default=4.0,
        help='Minimum SeqFF value (default: 4.0)'
    )
    
    parser.add_argument(
        '--max-seqff',
        type=float,
        default=30.0,
        help='Maximum SeqFF value (default: 30.0)'
    )
    
    parser.add_argument(
        '--min-mapping-rate',
        type=float,
        default=95.0,
        help='Minimum mapping rate %% (default: 95.0)'
    )
    
    parser.add_argument(
        '--max-duplication-rate',
        type=float,
        help='Maximum duplication rate %% (default: no limit)'
    )
    
    parser.add_argument(
        '--min-gc-content',
        type=float,
        help='Minimum GC content %% (default: no limit)'
    )
    
    parser.add_argument(
        '--max-gc-content',
        type=float,
        help='Maximum GC content %% (default: no limit)'
    )
    
    parser.add_argument(
        '--min-reads',
        type=int,
        help='Minimum number of reads (default: no limit)'
    )
    
    parser.add_argument(
        '--preview-only',
        action='store_true',
        help='Preview sample filtering only (no reference creation)'
    )
    
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Enable debug logging'
    )
    
    args = parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    
    logger.info("="*60)
    logger.info("NIPT Reference Creation Pipeline")
    logger.info("="*60)
    
    # 0. Generate sample list if analysis-dirs provided
    if args.analysis_dirs:
        logger.info("\n" + "="*60)
        logger.info("Step 1: Generating Sample List")
        logger.info("="*60)
        
        # Auto-generate sample list file name if not provided
        if not args.sample_list:
            args.sample_list = f"reference_sample_list_{args.prefix}.tsv"
        
        SampleListGenerator.generate_sample_list(
            analysis_dirs=args.analysis_dirs,
            output_file=args.sample_list,
            prefix=args.prefix,
            output_root=args.output_root
        )
        
        # If only generating sample list, exit
        if args.generate_sample_list:
            logger.info("\n" + "="*60)
            logger.info("Sample list generation complete (--generate-sample-list mode)")
            logger.info("="*60)
            sys.exit(0)
    
    # Validate sample list exists
    if not args.sample_list:
        logger.error("Error: --sample-list is required (or provide --analysis-dirs to generate)")
        sys.exit(1)
    
    if not os.path.exists(args.sample_list):
        logger.error(f"Sample list file not found: {args.sample_list}")
        sys.exit(1)
    
    # 1. Load and filter samples
    logger.info("\n" + "="*60)
    logger.info("Step 2: Loading and Filtering Samples")
    logger.info("="*60)
    
    selector = SampleSelector(args.sample_list)
    
    if not selector.load_samples():
        logger.error("Failed to load samples")
        sys.exit(1)
    
    # Check if already filtered
    is_already_filtered = '_filtered' in os.path.basename(args.sample_list)
    
    if is_already_filtered:
        logger.info("\n" + "="*60)
        logger.info("Using pre-filtered sample list (skipping filter step)")
        logger.info("="*60)
        # Still need to classify samples by gender
        selector._classify_by_gender()
    else:
        if not selector.filter_samples(
            min_seqff=args.min_seqff,
            max_seqff=args.max_seqff,
            min_mapping_rate=args.min_mapping_rate,
            max_duplication_rate=args.max_duplication_rate,
            min_gc_content=args.min_gc_content,
            max_gc_content=args.max_gc_content,
            min_reads=args.min_reads
        ):
            logger.error("Failed to filter samples")
            sys.exit(1)
        
        # Save filtered list
        filtered_output = args.sample_list.replace('.tsv', '_filtered.tsv')
        selector.save_filtered_list(filtered_output)
    
    if args.preview_only:
        logger.info("\n" + "="*60)
        logger.info("Preview mode - no reference creation")
        logger.info("="*60)
        sys.exit(0)
    
    # 2. Initialize configuration
    logger.info("\n" + "="*60)
    logger.info("Step 3: Reference Generation Configuration")
    logger.info("="*60)
    
    # 'all'이 있으면 모든 타입으로 확장
    if 'all' in args.ref_type:
        ref_types_to_create = ['ezd', 'prizm', 'wc', 'wcx']
    else:
        ref_types_to_create = args.ref_type
    
    logger.info(f"Reference types to create: {', '.join(ref_types_to_create)}")
    logger.info(f"Groups to process: {', '.join(args.groups)}")
    logger.info(f"Output directory: {args.output_dir if args.output_dir else f'data/refs/{args.labcode}'}")
    
    config = ReferenceConfig(args.labcode, args.output_dir, args.reference_source)
    
    # 3. Create references
    for group in args.groups:
        logger.info(f"\n{'#'*60}")
        logger.info(f"Processing group: {group}")
        logger.info(f"{'#'*60}")
        
        samples = selector.get_samples_by_group(group)
        
        logger.info(f"Samples for {group}: M={len(samples['M'])}, F={len(samples['F'])}")
        
        # ref_type은 이제 list이므로 'all' 포함 여부 또는 특정 타입 포함 여부 확인
        if 'all' in args.ref_type or 'ezd' in args.ref_type:
            creator = EZDReferenceCreator(config, samples, selector.df)
            creator.create_reference(group)
        
        if 'all' in args.ref_type or 'prizm' in args.ref_type:
            creator = PRIZMReferenceCreator(config, samples, selector.df)
            creator.create_reference(group)
        
        if 'all' in args.ref_type or 'wc' in args.ref_type:
            creator = WCReferenceCreator(config, samples, selector.df)
            creator.create_reference(group)
        
        if 'all' in args.ref_type or 'wcx' in args.ref_type:
            creator = WCXReferenceCreator(config, samples, selector.df)
            creator.create_reference(group)
    
    logger.info("\n" + "="*60)
    logger.info("Reference Creation Completed")
    logger.info("="*60)


if __name__ == "__main__":
    main()

