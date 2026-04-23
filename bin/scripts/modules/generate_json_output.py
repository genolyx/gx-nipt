#!/usr/bin/miniconda3/bin/python3
"""
---------------------------------------------
Generate the final json output file for NIPT Review Page
Reads actual data files and builds complete JSON structure

Author: Hyukjung Kwon
Contact: joykwon77@gmail.com
Updated: 2025-06-02
---
"""

import argparse
import glob
import json
import logging
import os
import sys

import pandas as pd

__author__ = "Hyukjung Kwon"
__email__ = "joykwon77@gmail.com"
__version__ = "1.0"

logger = logging.getLogger(__name__)
ch = logging.StreamHandler(sys.stdout)
ch.setFormatter(
    logging.Formatter("%(asctime)s -%(lineno)s- [%(levelname) 9s] - %(message)s")
)
logger.addHandler(ch)

# ==============================================================
# Configuration and Constants
# ==============================================================
CHROMOSOME_MAPPING = {
    "Chromosome 1": "chr1",
    "Chromosome 2": "chr2",
    "Chromosome 3": "chr3",
    "Chromosome 4": "chr4",
    "Chromosome 5": "chr5",
    "Chromosome 6": "chr6",
    "Chromosome 7": "chr7",
    "Chromosome 8": "chr8",
    "Chromosome 9": "chr9",
    "Chromosome 10": "chr10",
    "Chromosome 11": "chr11",
    "Chromosome 12": "chr12",
    "Chromosome 13": "chr13",
    "Chromosome 14": "chr14",
    "Chromosome 15": "chr15",
    "Chromosome 16": "chr16",
    "Chromosome 17": "chr17",
    "Chromosome 18": "chr18",
    "Chromosome 19": "chr19",
    "Chromosome 20": "chr20",
    "Chromosome 21": "chr21",
    "Chromosome 22": "chr22",
    "Chromosome X": "chrX",
    "Chromosome Y": "chrY",
}

MD8_DISEASE_MAPPING = {
    "1p36 deletion syndrome": {"item": "md1", "location": "1p36"},
    "2q33.1 deletion syndrome": {"item": "md2", "location": "2q33.1"},
    "Wolf-Hirschhorn syndrome": {"item": "md3", "location": "4p16.3"},
    "Cri Du Chat syndrome": {"item": "md4", "location": "5p-"},
    "Williams-Beuren syndrome": {"item": "md5", "location": "7q11.23"},
    "Jacobsen syndrome": {"item": "md6", "location": "11q23"},
    "Prader-willi/Angelman syndrome": {"item": "md7", "location": "15q11.2-q13"},
    "DiGeorge syndrome": {"item": "md8", "location": "22q11.2"},
}

APPID = "NIPT"


# ==============================================================
# Utility Functions
# ==============================================================
def file_check(parser, arg):
    if not os.path.exists(arg):
        parser.error(f"The file {arg} does not exist!")
    else:
        return str(arg)


def safe_read_csv(file_path, **kwargs):
    """Safely read CSV file with error handling"""
    try:
        if os.path.exists(file_path):
            return pd.read_csv(file_path, **kwargs)
        else:
            logger.warning(f"File {file_path} does not exist")
            return None
    except Exception as e:
        logger.error(f"Error reading file {file_path}: {e}")
        return None


def get_all_md_target_beds(config_file):
    with open(config_file) as f:
        cfg_json = json.load(f)
        bed_files = {
            key: value["bed"]
            for key, value in cfg_json.items()
            if key.startswith("MD_Target_") and "bed" in value
        }
    return bed_files


# ==============================================================
# Final Results Table Functions
# ==============================================================
def read_trisomy_detection_results(file_path):
    """Read trisomy detection results from TSV file"""
    results = []
    try:
        df = safe_read_csv(file_path, sep="\t")
        if df is not None and len(df.columns) >= 2:
            # Second column contains detection results
            detection_col = df.iloc[:, 1]  # Second column (index 1)
            disease_col = df.iloc[:, 0]  # First column (disease names)

            for i, detection in enumerate(detection_col):
                if detection != "Not Detected":
                    disease = (
                        disease_col.iloc[i] if i < len(disease_col) else f"Disease_{i}"
                    )
                    results.append(f"{disease} {detection}")

            if not results:
                results = ["Low Risk"]
    except Exception as e:
        logger.error(f"Error reading trisomy detection file {file_path}: {e}")
        results = ["Low Risk"]

    return results


def read_fetal_fraction_data(ff_path, gender_path):
    """
    Read fetal fraction data from TSV files

    Args:
        ff_path (str): Fetal fraction TSV 파일 경로 (sample.fetal_fraction.txt)
        gender_path (str): Gender detection TSV 파일 경로 (sample.gender.txt)

    Returns:
        dict: {
            'gender': str,     # gd_2의 gender 값
            'yff': float,      # YFF_2 값
            'seqff': float,    # SeqFF 값
            'ff_ratio': float  # YFF_2 / SeqFF 비율
        }
    """
    yff = 0.0
    seqff = 0.0
    gender = "Unknown"
    ff_ratio = 0.0

    try:
        # 1. Fetal fraction 파일 읽기
        try:
            ff_df = pd.read_csv(ff_path, sep="\t")

            # YFF_2와 SeqFF 값 추출
            for idx, row in ff_df.iterrows():
                ff_type = str(
                    row.iloc[0]
                )  # 첫 번째 컬럼 (Unnamed: 0 또는 첫 번째 컬럼)
                ff_value = float(row.iloc[1])  # 두 번째 컬럼 (value)

                if ff_type == "YFF_2":
                    yff = ff_value
                elif ff_type == "M-SeqFF":
                    seqff = ff_value

        except Exception as e:
            logger.error(f"Error reading fetal fraction file {ff_path}: {e}")
            yff = 0.0
            seqff = 0.0

        logger.info(f"yff : {yff}, seqff : {seqff}")
        # 2. Gender 파일 읽기
        try:
            gender_df = pd.read_csv(gender_path, sep="\t")

            # gd_2의 gender 값 추출
            for idx, row in gender_df.iterrows():
                gd_type = str(row.iloc[0])  # 첫 번째 컬럼
                logger.info(f"gd_type : {gd_type}")
                if (
                    gd_type == "gd_2" and len(row) >= 3
                ):  # gd_2 행이고 gender 컬럼이 있는 경우
                    gender = "Female" if str(row.iloc[2]) == "XX" else "Male"
                    break

        except Exception as e:
            logger.error(f"Error reading gender file {gender_path}: {e}")
            gender = "Unknown"

        logger.info(f"gender : {gender}")
        # 3. FF ratio 계산 (YFF_2 / SeqFF)
        try:
            if seqff != 0:
                ff_ratio = yff / seqff
            else:
                ff_ratio = 0.0
        except Exception as e:
            logger.error(f"Error calculating FF ratio: {e}")
            ff_ratio = 0.0

        logger.info(f"ff_ratio : {ff_ratio}")
        return {
            "gender": gender,
            "yff": yff,
            "seqff": seqff,
            "ff_ratio": round(ff_ratio, 2),
        }

    except Exception as e:
        logger.error(f"General error in read_fetal_fraction_data: {e}")
        return {"gender": "Unknown", "yff": 0.0, "seqff": 0.0, "ff_ratio": 0.0}


def read_sample_bias_qc(file_path):
    """Read sample bias QC from PRIZM QC file"""
    try:
        with open(file_path, "r") as f:
            content = f.read().strip()
            parts = content.split()
            if len(parts) >= 2:
                return parts[1]  # Return PASS/FAIL status
            return "Unknown"
    except Exception as e:
        logger.error(f"Error reading sample bias QC file {file_path}: {e}")
        return "Unknown"


def determine_final_trisomy_result(trisomy_details_chromosomes):
    """Determine final trisomy result based on chromosomes results"""
    final_results = []

    # Map chromosomes to trisomy names
    chr_to_trisomy = {
        "chr21": "T21",
        "chr18": "T18",
        "chr13": "T13",
        "chr9": "T9",
        "chr16": "T16",
        "chr22": "T22",
    }

    # Check each chromosome for High Risk status
    for chr_name, chr_data in trisomy_details_chromosomes.items():
        if chr_data.get("result") == "High Risk":
            if chr_name in chr_to_trisomy:
                final_results.append(f"{chr_to_trisomy[chr_name]} Detected")
            elif chr_name == "chrX":
                # For chrX, need to determine specific SCA condition based on data
                final_results.append("High Risk")
            elif chr_name == "chrY":
                final_results.append("High Risk")

    return final_results if final_results else ["Low Risk"]


def read_md_detection_results(wc_file, wcx_file, md_type):
    """Read microdeletion detection results from WC and WCX files"""
    results = []

    # Disease name mapping for md8
    md8_disease_names = {
        "1p36 deletion syndrome": "1p36 deletion syndrome",
        "2q33.1 deletion syndrome": "2q33.1 deletion syndrome",
        "Wolf-Hirschhorn syndrome": "Wolf-Hirschhorn syndrome",
        "Cri Du Chat syndrome": "Cri Du Chat syndrome",
        "Williams-Beuren syndrome": "Williams-Beuren syndrome",
        "Jacobsen syndrome": "Jacobsen syndrome",
        "Prader-willi/Angelman syndrome": "Prader-willi/Angelman syndrome",
        "DiGeorge syndrome": "DiGeorge syndrome",
    }

    detected_diseases = set()

    # Read WC results
    wc_df = safe_read_csv(wc_file, sep="\t")
    if wc_df is not None and "Disease-name" in wc_df.columns:
        for _, row in wc_df.iterrows():
            disease = row["Disease-name"]
            if md_type == "md8":
                if disease in md8_disease_names:
                    detected_diseases.add(f"{disease} Detected")
            else:
                detected_diseases.add(f"{md_type} Detected")
                break  # For md87, md116, md320, just indicate the category

    # Read WCX results
    wcx_df = safe_read_csv(wcx_file, sep="\t")
    if wcx_df is not None and "Disease-name" in wcx_df.columns:
        for _, row in wcx_df.iterrows():
            disease = row["Disease-name"]
            if md_type == "md8":
                if disease in md8_disease_names:
                    detected_diseases.add(f"{disease} Detected")
            else:
                detected_diseases.add(f"{md_type} Detected")
                break

    return list(detected_diseases) if detected_diseases else []


def build_final_results_table(analysis_dir, sample_name):
    """Build final results table from actual data files"""

    # 1. Read trisomy detection results
    orig_trisomy_file = f"{analysis_dir}/{sample_name}/Output_EZD/orig/Trisomy_detect_result_orig_with_SCA.tsv"
    fetus_trisomy_file = f"{analysis_dir}/{sample_name}/Output_EZD/fetus/Trisomy_detect_result_fetus_with_SCA.tsv"
    mom_trisomy_file = f"{analysis_dir}/{sample_name}/Output_EZD/mom/Trisomy_detect_result_mom_with_SCA.tsv"

    orig_results = read_trisomy_detection_results(orig_trisomy_file)
    fetus_results = read_trisomy_detection_results(fetus_trisomy_file)
    mom_results = read_trisomy_detection_results(mom_trisomy_file)

    # 2. Read fetal fraction data
    ff_file = f"{analysis_dir}/{sample_name}/Output_FF/{sample_name}.fetal_fraction.txt"
    gender_file = f"{analysis_dir}/{sample_name}/Output_FF/{sample_name}.gender.txt"
    ff_gender_data = read_fetal_fraction_data(ff_file, gender_file)

    # 3. Read sample bias QC
    sample_bias_file = f"{analysis_dir}/{sample_name}/Output_PRIZM/orig/{sample_name}.of_orig.prizm.qc.txt"
    sample_bias = read_sample_bias_qc(sample_bias_file)

    # 4. Read trisomy details chromosomes to determine final result (simplified)
    """
    final_trisomy_results = []
    try:
        # Just check if any detection files indicate positive results
        #final_trisomy_results = []
        
        for result_list in [orig_results, fetus_results]:
            for result in result_list:
                logger.info(f">>>>> trisomy_result : {result}")
                #if result != "Not Detected" and "Detected" in result:
                logger.info(f"result = {result}")
                #if result != "Not Detected" in result:
                if "Low Risk" not in result:
                    final_trisomy_results.append("High Risk")
                    logger.info(f"result : {result}, final_trisomy_results : {final_trisomy_results}")
        
        if not final_trisomy_results:
            logger.info(f"final_trisomy_results : {final_trisomy_results}")
            final_trisomy_results = ["Low Risk"]
            
    except Exception as e:
        logger.error(f"Error determining final trisomy result: {e}")
        final_trisomy_results = ["Low Risk"]

    # 5. Read MD detection results
    md_results = []
    
    # MD8 results
    wc_md8_file = f"{analysis_dir}/{sample_name}/Output_WC/orig/{sample_name}_WC_orig_md8.tsv"
    wcx_md8_file = f"{analysis_dir}/{sample_name}/Output_WCX/orig/{sample_name}_WCX_orig_md8.tsv"
    md8_results = read_md_detection_results(wc_md8_file, wcx_md8_file, 'md8')
    logger.info(md8_results)
    md_results.extend(md8_results)
    
    # MD87, MD116, MD320 results
    for md_type in ['md87', 'md116', 'md320']:
        wc_file = f"{analysis_dir}/{sample_name}/Output_WC/orig/{sample_name}_WC_orig_{md_type}.tsv"
        wcx_file = f"{analysis_dir}/{sample_name}/Output_WCX/orig/{sample_name}_WCX_orig_{md_type}.tsv"
        md_type_results = read_md_detection_results(wc_file, wcx_file, md_type)
        logger.info(md_type_results)
        md_results.extend(md_type_results)
   
    logger.info(f"md_results : {md_results}")
    """

    return {
        "original": orig_results,
        "fetus": fetus_results,
        "mom": mom_results,
        "fetal_gender": ff_gender_data["gender"],
        # 250610 : fetal_fraction_yff is set to "N/A" for gender Female
        "fetal_fraction_yff": "N/A"
        if ff_gender_data["gender"] == "Female"
        else ff_gender_data["yff"],
        "fetal_fraction_seqff": ff_gender_data["seqff"],
        "ff_ratio": ff_gender_data["ff_ratio"],
        "sample_bias_qc": sample_bias,
        #'final_trisomy_result': final_trisomy_results,
        #'md_results': md_results if md_results else ["Low Risk"]
    }


# ==============================================================
# Data Reading Functions
# ==============================================================
def read_risk_before_data(age, common_data_dir):
    """Read risk_before data based on age from CSV files"""

    # Ensure age is integer and clamp to valid range
    try:
        age = int(age)
        age = max(25, min(45, age))
    except (ValueError, TypeError):
        logger.warning(f"Invalid age value: {age}, using default age 35")
        age = 35

    try:
        # Read single pregnancy risk data
        single_file = os.path.join(common_data_dir, "Single_risk_before.csv")
        single_df = safe_read_csv(single_file, sep="\t")

        # Read twin pregnancy risk data
        twin_file = os.path.join(common_data_dir, "Twin_risk_before.csv")
        twin_df = safe_read_csv(twin_file, sep="\t")

        if single_df is None or twin_df is None:
            logger.warning(f"Could not read risk data files from {common_data_dir}")
            return {}

        # Ensure Age column is integer type
        if "Age" in single_df.columns:
            single_df["Age"] = pd.to_numeric(single_df["Age"], errors="coerce")
        if "Age" in twin_df.columns:
            twin_df["Age"] = pd.to_numeric(twin_df["Age"], errors="coerce")

        # Find row for the specified age
        single_row = single_df[single_df["Age"] == age]
        twin_row = twin_df[twin_df["Age"] == age]

        if single_row.empty or twin_row.empty:
            logger.warning(f"Could not find risk data for age {age}")
            return {}

        # Extract risk values for both single and twin pregnancies
        risk_data = {
            "T21": {
                "single": single_row.iloc[0]["Trisomy_21_risk"],
                "twin": twin_row.iloc[0]["T21_twin_odibo"],
            },
            "T18": {
                "single": single_row.iloc[0]["Trisomy_18_risk"],
                "twin": twin_row.iloc[0]["T18_twin_odibo"],
            },
            "T13": {
                "single": single_row.iloc[0]["Trisomy_13_risk"],
                "twin": twin_row.iloc[0]["T13_twin_odibo"],
            },
        }

        return risk_data

    except Exception as e:
        logger.error(f"Error reading risk_before data: {e}")
        return {}


def build_trisomy_results(
    analysis_dir, sample_name, fetus_gender, age, common_data_dir
):
    """Build trisomy results from EZD detection result files"""

    trisomy_order = [
        "Trisomy21",
        "Trisomy18",
        "Trisomy13",
        "Trisomy9",
        "Trisomy16",
        "Trisomy22",
        "XO",
        "XXX",
        "XXY",
        "XYY",
        "other",
    ]

    item_mapping = {
        "Trisomy21": "T21",
        "Trisomy18": "T18",
        "Trisomy13": "T13",
        "Trisomy9": "T9",
        "Trisomy16": "T16",
        "Trisomy22": "T22",
        "XO": "XO",
        "XXX": "XXX",
        "XXY": "XXY",
        "XYY": "XYY",
        "other": "other",
    }

    # Chromosome to trisomy mapping (기존 주요 trisomy들)
    chr_to_trisomy = {
        "chr21": "Trisomy21",
        "chr18": "Trisomy18",
        "chr13": "Trisomy13",
        "chr9": "Trisomy9",
        "chr16": "Trisomy16",
        "chr22": "Trisomy22",
        "chrX": ["XO", "XXX"],
        "chrY": ["XXY", "XYY"],  # chrX and chrY can have multiple conditions
    }

    # Other trisomies (주요 trisomy에 포함되지 않은 염색체들)
    other_trisomy_order = [
        "Trisomy1",
        "Trisomy2",
        "Trisomy3",
        "Trisomy4",
        "Trisomy5",
        "Trisomy6",
        "Trisomy7",
        "Trisomy8",
        "Trisomy10",
        "Trisomy11",
        "Trisomy12",
        "Trisomy14",
        "Trisomy15",
        "Trisomy17",
        "Trisomy19",
        "Trisomy20",
    ]

    # Other chromosome mapping
    other_chr_to_trisomy = {
        "chr1": "Trisomy1",
        "chr2": "Trisomy2",
        "chr3": "Trisomy3",
        "chr4": "Trisomy4",
        "chr5": "Trisomy5",
        "chr6": "Trisomy6",
        "chr7": "Trisomy7",
        "chr8": "Trisomy8",
        "chr10": "Trisomy10",
        "chr11": "Trisomy11",
        "chr12": "Trisomy12",
        "chr14": "Trisomy14",
        "chr15": "Trisomy15",
        "chr17": "Trisomy17",
        "chr19": "Trisomy19",
        "chr20": "Trisomy20",
    }

    results = []

    # Read risk_before data based on age
    risk_before_data = read_risk_before_data(age, common_data_dir)

    # Read both original and fetus results for comparison
    orig_file = f"{analysis_dir}/{sample_name}/Output_EZD/orig/Trisomy_detect_result_orig_with_SCA.tsv"
    fetus_file = f"{analysis_dir}/{sample_name}/Output_EZD/fetus/Trisomy_detect_result_fetus_with_SCA.tsv"

    orig_df = safe_read_csv(orig_file, sep="\t")
    fetus_df = safe_read_csv(fetus_file, sep="\t")

    if orig_df is None:
        logger.warning(f"Could not read original trisomy results from {orig_file}")
        orig_df = pd.DataFrame(columns=["chr", "result", "Z", "UAR"])

    if fetus_df is None:
        logger.warning(f"Could not read fetus trisomy results from {fetus_file}")
        fetus_df = pd.DataFrame(columns=["chr", "result", "Z", "UAR"])

    # Process each trisomy in order
    required_columns = ["chr", "result", "Z", "UAR"]
    if not all(col in orig_df.columns for col in required_columns):
        logger.error(
            f"Missing required columns in {orig_file}. Available: {list(orig_df.columns)}"
        )
        orig_df = pd.DataFrame(columns=required_columns)

    if not all(col in fetus_df.columns for col in required_columns):
        logger.error(
            f"Missing required columns in {fetus_file}. Available: {list(fetus_df.columns)}"
        )
        fetus_df = pd.DataFrame(columns=required_columns)

    # Process main trisomies (기존 로직)
    for trisomy in trisomy_order:
        # Skip "other" - will be processed separately
        if trisomy == "other":
            continue

        # Find corresponding chromosome(s)
        target_chrs = []
        for chr_name, trisomy_list in chr_to_trisomy.items():
            if isinstance(trisomy_list, list):
                if trisomy in trisomy_list:
                    target_chrs.append(chr_name)
            else:
                if trisomy == trisomy_list:
                    target_chrs.append(chr_name)

        # If no chromosome mapping found, skip
        if not target_chrs:
            continue

        # Get results from both dataframes
        orig_detected = False
        fetus_detected = False

        logger.info(f"target_chrs : {target_chrs}")
        for chr_name in target_chrs:
            logger.info(f"chr_name = {chr_name}")
            # Check original results
            if chr_name in orig_df["chr"].values:
                orig_row = orig_df[orig_df["chr"] == chr_name].iloc[0]

                # SCA 처리: chrX / chrY 분리
                if chr_name == "chrX":
                    # XO, XXX 은 여아에서만 계산
                    if (
                        fetus_gender in ("F", "Female")
                        and trisomy in orig_row["result"]
                    ):
                        orig_detected = True

                elif chr_name == "chrY":
                    # XXY, XYY 은 남아에서만 계산
                    if fetus_gender in ("M", "Male") and trisomy in orig_row["result"]:
                        orig_detected = True
                else:
                    # For autosomal trisomies, check if detected
                    # if "Not Detected" not in orig_row['result']:
                    if orig_row["result"] != "Not Detected":
                        orig_detected = True

            logger.info(f"chr_name : {chr_name}, orig_detected = {orig_detected}")

            # Check fetus results
            if chr_name in fetus_df["chr"].values:
                fetus_row = fetus_df[fetus_df["chr"] == chr_name].iloc[0]

                # Check for specific SCA conditions in result column
                # SCA 처리: chrX / chrY 분리
                if chr_name == "chrX":
                    # XO, XXX 은 여아에서만 계산
                    if (
                        fetus_gender in ("F", "Female")
                        and trisomy in fetus_row["result"]
                    ):
                        orig_detected = True

                elif chr_name == "chrY":
                    # XXY, XYY 은 남아에서만 계산
                    if fetus_gender in ("M", "Male") and trisomy in fetus_row["result"]:
                        orig_detected = True
                else:
                    # For autosomal trisomies, check if detected
                    if fetus_row["result"] != "Not Detected":
                        # if "Not Detected" not in fetus_row['result']:
                        fetus_detected = True

            logger.info(f"chr_name : {chr_name}, fetus_detected = {fetus_detected}")

        # Determine final result based on both orig and fetus results
        if orig_detected and fetus_detected:
            result_status = "High Risk"
        elif orig_detected or fetus_detected:
            result_status = "High Risk"  # Suspected
        else:
            result_status = "Low Risk"

        # Determine risk_after based on detection result (single value, not list)
        if result_status == "High Risk":
            risk_after = "90/100"
        else:
            risk_after = "<2/10,000"

        # Get risk_before from CSV data (both single and twin)
        trisomy_item = item_mapping[trisomy]
        risk_before_single = None
        risk_before_twin = None

        if trisomy_item in risk_before_data:
            risk_before_single = risk_before_data[trisomy_item]["single"]
            risk_before_twin = risk_before_data[trisomy_item]["twin"]

        # Set PPV/NPV for specific trisomies only
        ppv_npv_trisomies = [
            "Trisomy13",
            "Trisomy18",
            "Trisomy21",
            "XXX",
            "XO",
            "XYY",
            "XXY",
        ]
        if trisomy in ppv_npv_trisomies:
            ppv = ">99"
            npv = ">99"
        else:
            ppv = None
            npv = None

        # Create result entry (removed z_score and uar_percent)
        result = {
            "item": item_mapping[trisomy],
            "disease_name": trisomy,
            "result": result_status,
            "risk_before_single": risk_before_single,  # From Single CSV file based on age
            "risk_before_twin": risk_before_twin,  # From Twin CSV file based on age
            "risk_after": risk_after,  # Single value based on detection result
            "ppv": ppv,  # ">99" for T13/18/21 and SCA, null for others
            "npv": npv,  # ">99" for T13/18/21 and SCA, null for others
        }
        results.append(result)

    # ================================
    # Process "other" trisomies
    # ================================
    other_detected = False

    # Check for any detection in other chromosomes
    for chr_name, trisomy_name in other_chr_to_trisomy.items():
        # Check original results
        if chr_name in orig_df["chr"].values:
            orig_row = orig_df[orig_df["chr"] == chr_name].iloc[0]
            if "Not Detected" not in orig_row["result"]:
                other_detected = True
                logger.info(
                    f"Other trisomy detected in orig: {chr_name} - {orig_row['result']}"
                )
                break

        # Check fetus results
        if chr_name in fetus_df["chr"].values:
            fetus_row = fetus_df[fetus_df["chr"] == chr_name].iloc[0]
            if "Not Detected" not in fetus_row["result"]:
                other_detected = True
                logger.info(
                    f"Other trisomy detected in fetus: {chr_name} - {fetus_row['result']}"
                )
                break

    # Create "other" result entry
    if other_detected:
        other_result_status = "High Risk"
        other_risk_after = "90/100"
    else:
        other_result_status = "Low Risk"
        other_risk_after = "<2/10,000"

    # Get risk_before for "other" if available
    other_risk_before_single = None
    other_risk_before_twin = None
    if "other" in risk_before_data:
        other_risk_before_single = risk_before_data["other"]["single"]
        other_risk_before_twin = risk_before_data["other"]["twin"]

    # "other" does not have PPV/NPV
    other_result = {
        "item": "other",
        "disease_name": "other",
        "result": other_result_status,
        "risk_before_single": other_risk_before_single,
        "risk_before_twin": other_risk_before_twin,
        "risk_after": other_risk_after,
        "ppv": None,
        "npv": None,
    }
    results.append(other_result)

    logger.info("------------------------------------------------------")
    # logger.info(results)
    return results


def read_chromosome_analysis_from_ezd(file_path):
    """Read chromosome analysis results from EZD detection files - For compatibility with existing code"""
    df = safe_read_csv(file_path, sep="\t")
    if df is None:
        return {}

    chromosomes = {}

    for _, row in df.iterrows():
        chr_name = row["chr"]

        # Convert "Not Detected" to "Low Risk" and "Detected" to "High Risk"
        detection_result = row["result"]
        if detection_result == "Not Detected":
            risk_level = "Low Risk"
        elif detection_result == "Detected":
            risk_level = "High Risk"
        else:
            risk_level = "Low Risk"  # Default fallback

        chromosomes[chr_name] = {
            "result": risk_level,
            "z_score": row["Z"] if pd.notna(row["Z"]) else None,
            "uar_percent": row["UAR"] if pd.notna(row["UAR"]) else None,
            "z_threshold": "N<1.7~2.7<D",  # Standard threshold
            "uar_threshold": "N<8.48~8.53<D",  # Standard threshold
        }

    return chromosomes


def read_chromosome_analysis_from_ezd_prizm_detailed(
    ezd_file_path, prizm_file_path, fetus_gender
):
    """Read chromosome analysis results from EZD & PRIZM detection files with full chromosome table"""
    logger.info(f"Reading EZD detailed data from: {ezd_file_path}")
    logger.info(f"Reading PRIZM detailed data from: {prizm_file_path}")

    ezd_df = safe_read_csv(ezd_file_path, sep="\t")
    if ezd_df is None:
        logger.warning(f"Failed to read EZD file: {ezd_file_path}")
        return {}

    logger.info(
        f"EZD file read successfully. Shape: {ezd_df.shape}, Columns: {list(ezd_df.columns)}"
    )

    prizm_df = pd.read_csv(prizm_file_path, sep="\t", comment="#")
    if prizm_df is None:
        logger.warning(f"Failed to read PRIZM file: {prizm_file_path}")
        return {}

    logger.info(
        f"PRIZM file read successfully. Shape: {prizm_df.shape}, Columns: {list(prizm_df.columns)}"
    )

    chromosomes = {}
    merged_df = pd.merge(
        ezd_df,
        prizm_df[["Chromosome", "Decision"]],
        left_on="chr",
        right_on="Chromosome",
        how="inner",
    )

    for _, row in merged_df.iterrows():
        chr_name = row["chr"]

        # Convert chromosome name to full name for display
        chr_display_name = chr_name.replace("chr", "Chromosome ")
        if chr_display_name == "Chromosome X":
            chr_display_name = "Chromosome X"
        elif chr_display_name == "Chromosome Y":
            chr_display_name = "Chromosome Y"

        # EZD result
        if chr_name == "chrX" and fetus_gender in ("M", "Male"):
            ezd_detection = "Low Risk"
        elif chr_name == "chrY" and fetus_gender in ("F", "Female"):
            ezd_detection = "Low Risk"
        else:
            ezd_result = row["result"]
            if ezd_result == "Not Detected":
                ezd_detection = "Low Risk"  # or "Not Detected" based on preference
            elif ezd_result == "Detected" or "Suspected":
                ezd_detection = "High Risk"
            else:
                ezd_detection = "Low Risk"  # Default fallback

        # PRIZM result
        if chr_name == "chrX" and fetus_gender in ("M", "Male"):
            prizm_detection = "Low Risk"
        if chr_name == "chrY":
            prizm_detection = "Low Risk"
        else:
            prizm_decision = row["Decision"]
            prizm_detection = (
                "High Risk"
                if prizm_decision in ("Detected", "Suspected")
                else "Low Risk"
            )

        logger.debug(
            f"Processing {chr_name} -> {chr_display_name}: EZD={ezd_result} -> {ezd_detection}, PRIZM={prizm_decision} -> {prizm_detection}"
        )

        # For chrX and chrY, set thresholds to null initially
        if chr_name in ["chrX", "chrY"]:
            z_threshold = None
            uar_threshold = None
        else:
            # Will be updated later from threshold file
            z_threshold = None  # To be filled from threshold file
            uar_threshold = None  # To be filled from threshold file

        chromosomes[chr_display_name] = {
            "EZD Detection": ezd_detection,
            "PRIZM Detection": prizm_detection,
            "Z-score": row["Z"] if pd.notna(row["Z"]) else None,
            "UAR(%)": row["UAR"] if pd.notna(row["UAR"]) else None,
            "Z-score threshold": z_threshold,
            "UAR threshold": uar_threshold,
            "checked": ezd_detection == "High Risk" or prizm_detection == "High Risk",
        }

    logger.info(f"Processed {len(chromosomes)} chromosomes from EZD + PRIZM data")
    return chromosomes

def read_threshold_data(threshold_file_path):
    """Read threshold data from external file"""
    try:
        threshold_df = safe_read_csv(threshold_file_path, sep="\t")
        if threshold_df is None:
            return {}

        # Convert to dictionary for easy lookup
        threshold_dict = {}
        for _, row in threshold_df.iterrows():
            chr_name = row.get("chr", "")

            # Convert chromosome name to display format
            chr_display_name = chr_name.replace("chr", "Chromosome ")
            if chr_display_name == "Chromosome X":
                chr_display_name = "Chromosome X"
            elif chr_display_name == "Chromosome Y":
                chr_display_name = "Chromosome Y"

            # For chrX and chrY, set thresholds to null
            if chr_name in ["chrX", "chrY"]:
                threshold_dict[chr_display_name] = {
                    "z_threshold": None,
                    "uar_threshold": None,
                }
            else:
                # Format: L<UAR_min~UAR_max<H and L<Z_min~Z_max<H
                # Convert to string to handle both numeric and string values safely
                uar_min = (
                    str(row.get("UAR_min", ""))
                    if pd.notna(row.get("UAR_min", ""))
                    else ""
                )
                uar_max = (
                    str(row.get("UAR_max", ""))
                    if pd.notna(row.get("UAR_max", ""))
                    else ""
                )
                z_min = (
                    str(row.get("Z_min", "")) if pd.notna(row.get("Z_min", "")) else ""
                )
                z_max = (
                    str(row.get("Z_max", "")) if pd.notna(row.get("Z_max", "")) else ""
                )

                uar_threshold = (
                    f"L < {uar_min}~{uar_max} < H"
                    if uar_min and uar_max
                    else "L < 8.48~8.53 < H"
                )
                z_threshold = (
                    f"L < {z_min}~{z_max} < H" if z_min and z_max else "L < 1.7~2.7 < H"
                )

                threshold_dict[chr_display_name] = {
                    "z_threshold": z_threshold,
                    "uar_threshold": uar_threshold,
                }

        return threshold_dict

    except Exception as e:
        logger.warning(f"Could not read threshold data from {threshold_file_path}: {e}")
        return {}


def read_qc_data(file_path):
    """Read QC data with proper error handling"""
    try:
        df = safe_read_csv(file_path, sep="\t", header=None, index_col=0)
        if df is None or df.empty:
            logger.warning(f"QC file {file_path} is empty or could not be read")
            return pd.DataFrame()  # 빈 DataFrame 반환

        df.columns = ["value", "status"]
        return df
    except Exception as e:
        logger.error(f"Error reading QC data from {file_path}: {e}")
        return pd.DataFrame()  # 빈 DataFrame 반환


# ==============================================================
# Microdeletion Functions
# ==============================================================


def process_md_detection(wc_file, wcx_file, data_src, target_bed_file, md_type="md8"):
    """Process WC and WCX detection results with default values - Updated for md108, md87, md320"""

    results = {}
    wcx_chr_list = []

    # md8의 경우 고정된 8개 질병을 사용
    if md_type == "md8":
        # 고정된 md8 질병 목록
        fixed_md8_diseases = [
            {"name": "1p36 deletion syndrome", "location": "1p36"},
            {"name": "2q33.1 deletion syndrome", "location": "2q33.1"},
            {"name": "Wolf-Hirschhorn syndrome", "location": "4p16.3"},
            {"name": "Cri Du Chat syndrome", "location": "5p-"},
            {"name": "Williams-Beuren syndrome", "location": "7q11.23"},
            {"name": "Jacobsen syndrome", "location": "11q23"},
            {"name": "Prader-willi/Angelman syndrome", "location": "15q11.2-q13"},
            {"name": "DiGeorge syndrome", "location": "22q11.2"},
        ]

        # md1~md8까지 고정으로 생성
        for i, disease_info in enumerate(fixed_md8_diseases, 1):
            md_key = f"md{i}"
            results[md_key] = {
                "disease_name": disease_info["name"],
                "target_region": disease_info["location"],  # 간략한 위치 정보
                "detection": {"WC": None, "WCX": None},
                "detected_region": {"WC": None, "WCX": None},
                "length": {"WC": None, "WCX": None},
                "z_score": {"WC": None, "WCX": None},
                "detected_region_link": {"WC": "", "WCX": ""},
                "image": {"WC": "", "WCX": ""},
                "checked": False,
            }

        # 이제 실제 검출 결과로 업데이트 (WC 파일)
        wc_df = safe_read_csv(wc_file, sep="\t")
        if wc_df is not None and "Disease-name" in wc_df.columns:
            for _, row in wc_df.iterrows():
                disease = row["Disease-name"]

                # 해당 질병에 맞는 md key 찾기
                for md_key, md_data in results.items():
                    if md_data["disease_name"] == disease:
                        detected_region = f"{row['chr']}:{row['start']}-{row['end']}"
                        results[md_key]["detection"]["WC"] = "High Risk"
                        results[md_key]["detected_region"]["WC"] = detected_region
                        results[md_key]["length"]["WC"] = str(row["length"])
                        results[md_key]["z_score"]["WC"] = str(row["zscore"])
                        results[md_key]["detected_region_link"]["WC"] = (
                            f"https://deciphergenomics.org/browser#q/grch37:{detected_region}"
                        )
                        results[md_key]["checked"] = True
                        break

        # WCX 파일 처리
        wcx_df = safe_read_csv(wcx_file, sep="\t")
        if wcx_df is not None and "Disease-name" in wcx_df.columns:
            for _, row in wcx_df.iterrows():
                disease = row["Disease-name"]

                # 해당 질병에 맞는 md key 찾기
                for md_key, md_data in results.items():
                    if md_data["disease_name"] == disease:
                        detected_region = f"{row['chr']}:{row['start']}-{row['end']}"
                        results[md_key]["detection"]["WCX"] = "High Risk"
                        results[md_key]["detected_region"]["WCX"] = detected_region
                        results[md_key]["length"]["WCX"] = str(row["length"])
                        results[md_key]["z_score"]["WCX"] = str(row["zscore"])
                        results[md_key]["detected_region_link"]["WCX"] = (
                            f"https://deciphergenomics.org/browser#q/grch37:{detected_region}"
                        )
                        results[md_key]["image"]["WCX"] = (
                            (f"Output_WCX/chr_plots/{data_src}/chr{row['chr']}.png"),
                        )
                        results[md_key]["checked"] = True
                        wcx_chr_list.append(row["chr"])
                        break

    else:
        # md87, md108, md320의 경우: 검출된 결과를 순서대로 md1, md2... 로 매핑

        # 모든 검출 결과 수집 (WC와 WCX 모두)
        all_detections = {}  # Disease-name을 키로 하여 WC, WCX 정보를 모음

        # WC 파일에서 검출 결과 읽기
        wc_df = safe_read_csv(wc_file, sep="\t")
        if wc_df is not None and "Disease-name" in wc_df.columns:
            for _, row in wc_df.iterrows():
                disease = row["Disease-name"]
                if disease not in all_detections:
                    all_detections[disease] = {"WC": None, "WCX": None}

                detected_region = f"{row['chr']}:{row['start']}-{row['end']}"
                all_detections[disease]["WC"] = {
                    "detected_region": detected_region,
                    "length": str(row["length"]),
                    "z_score": str(row["zscore"]),
                    "chr": row["chr"],
                    "start": row["start"],
                    "end": row["end"],
                }

        # WCX 파일에서 검출 결과 읽기
        wcx_df = safe_read_csv(wcx_file, sep="\t")
        if wcx_df is not None and "Disease-name" in wcx_df.columns:
            for _, row in wcx_df.iterrows():
                disease = row["Disease-name"]
                if disease not in all_detections:
                    all_detections[disease] = {"WC": None, "WCX": None}

                detected_region = f"{row['chr']}:{row['start']}-{row['end']}"
                all_detections[disease]["WCX"] = {
                    "detected_region": detected_region,
                    "length": str(row["length"]),
                    "z_score": str(row["zscore"]),
                    "chr": row["chr"],
                    "start": row["start"],
                    "end": row["end"],
                }

        # 검출된 결과가 있으면 순서대로 md1, md2... 로 매핑
        if all_detections:
            # 검출된 질병들을 정렬 (일관된 순서를 위해)
            sorted_diseases = sorted(all_detections.keys())

            for i, disease in enumerate(sorted_diseases, 1):
                md_key = f"md{i}"
                detection_data = all_detections[disease]

                # 기본 구조 생성
                results[md_key] = {
                    "disease_name": disease,
                    "target_region": detection_data["WC"]["detected_region"]
                    if detection_data["WC"]
                    else detection_data["WCX"]["detected_region"],
                    "detection": {"WC": None, "WCX": None},
                    "detected_region": {"WC": None, "WCX": None},
                    "length": {"WC": None, "WCX": None},
                    "z_score": {"WC": None, "WCX": None},
                    "detected_region_link": {"WC": "", "WCX": ""},
                    "image": {"WC": "", "WCX": ""},
                    "checked": True,
                }

                # WC 데이터 채우기
                if detection_data["WC"]:
                    wc_data = detection_data["WC"]
                    results[md_key]["detection"]["WC"] = "High Risk"
                    results[md_key]["detected_region"]["WC"] = wc_data[
                        "detected_region"
                    ]
                    results[md_key]["length"]["WC"] = wc_data["length"]
                    results[md_key]["z_score"]["WC"] = wc_data["z_score"]
                    results[md_key]["detected_region_link"]["WC"] = (
                        f"https://deciphergenomics.org/browser#q/grch37:{wc_data['detected_region']}"
                    )

                # WCX 데이터 채우기
                if detection_data["WCX"]:
                    wcx_data = detection_data["WCX"]
                    results[md_key]["detection"]["WCX"] = "High Risk"
                    results[md_key]["detected_region"]["WCX"] = wcx_data[
                        "detected_region"
                    ]
                    results[md_key]["length"]["WCX"] = wcx_data["length"]
                    results[md_key]["z_score"]["WCX"] = wcx_data["z_score"]
                    results[md_key]["detected_region_link"]["WCX"] = (
                        f"https://deciphergenomics.org/browser#q/grch37:{wcx_data['detected_region']}"
                    )
                    results[md_key]["image"]["WCX"] = (
                        f"Output_WCX/chr_plots/{data_src}/chr{wcx_data['chr']}.png"
                    )
                    wcx_chr_list.append(wcx_data["chr"])

        else:
            # 검출된 결과가 없으면 기본 md1만 null 값으로 생성
            results["md1"] = {
                "disease_name": None,
                "target_region": None,
                "detection": {"WC": None, "WCX": None},
                "detected_region": {"WC": None, "WCX": None},
                "length": {"WC": None, "WCX": None},
                "z_score": {"WC": None, "WCX": None},
                "detected_region_link": {"WC": "", "WCX": ""},
                "image": {"WC": "", "WCX": ""},
                "checked": False,
            }

    return results, wcx_chr_list


# build_nipt_json 함수에서 md_details 섹션 수정
def build_md_details_section(analysis_dir, sample_name, target_bed_dir):
    """Build md_details section with proper default values"""

    md_sections = {
        "md8_results": ("TargetDB_md8.bed", "md8", "md8"),
        "md108_results": ("TargetDB_md108.bed", "md108", "other_md108"),
        "md320_results": ("TargetDB_md320.bed", "md320", "other_md320"),
        "md87_results": ("TargetDB_md87.bed", "md87", "other_md87"),
        "md141_results": ("TargetDB_md141.bed", "md141", "other_md141"),
    }

    data_sources = ["orig", "fetus", "mom"]

    wcx_chr_list_total = {g: [] for g in data_sources}
    _seen = {g: set() for g in data_sources}

    md_details = {}
    detected_md8_list = set()
    detected_others = set()

    detected_md_set = set()  # Set to track detected microdeletions

    for md_section, (bed_file, md_type, present_name) in md_sections.items():
        md_details[md_section] = {}

        # Process each data source
        for data_src in data_sources:
            # Image files
            md_details[md_section][data_src] = {
                "image": {
                    "WC": f"Output_WC/{sample_name}.wc.{data_src}_z.png",
                    "WCX": f"Output_WCX/{sample_name}.wcx.{data_src}.png",
                }
            }

            # MD detection files
            wc_file = f"{analysis_dir}/{sample_name}/Output_WC/{data_src}/{sample_name}_WC_{data_src}_{md_type}.tsv"
            wcx_file = f"{analysis_dir}/{sample_name}/Output_WCX/{data_src}/{sample_name}_WCX_{data_src}_{md_type}.tsv"
            target_bed_path = f"{target_bed_dir}/{bed_file}"

            detections, wcx_chr_list = process_md_detection(
                wc_file, wcx_file, data_src, target_bed_path, md_type
            )

            logger.info(f"wcx_chr_list for {data_src}: {wcx_chr_list}")
            # wcx_chr_list_total에 중복 없이 추가
            for chrom in wcx_chr_list:
                if chrom not in _seen[data_src]:
                    _seen[data_src].add(chrom)
                    wcx_chr_list_total[data_src].append(chrom)

            # Add detection results
            for md_key, md_data in detections.items():
                md_details[md_section][data_src][md_key] = md_data

                if (
                    isinstance(md_data, dict)
                    and "High Risk" in md_data.get("detection", {}).values()
                ):
                    if md_section == "md8_results":
                        disease_name = md_data.get("disease_name")
                        if disease_name:
                            detected_md8_list.add(disease_name)
                    else:
                        detected_others.add(present_name)

    logger.info(f"[WCX chromosome plots-orig] {wcx_chr_list_total['orig']}")
    logger.info(f"[WCX chromosome plots-fetus] {wcx_chr_list_total['fetus']}")
    logger.info(f"[WCX chromosome plots-mom] {wcx_chr_list_total['mom']}")
    logger.info(f"Detected microdeletions: {detected_md_set}")

    return (
        md_details,
        sorted(detected_md8_list),
        sorted(detected_others),
        wcx_chr_list_total,
    )


def find_fastqc_reports(qc_dir, sample_subdir=""):
    """
    qc_dir: Output_QC 가 있는 최상위 경로
    sample_subdir: sample 폴더 안에 있으면 그 이름(예: sample_name)
    """
    base = os.path.join(qc_dir, sample_subdir)
    # R1/R2 HTML 파일 전부 스캔
    all_htmls = glob.glob(os.path.join(base, "*_fastqc.html"))

    fastqc_r1 = next(
        (os.path.relpath(p, base) for p in all_htmls if "_R1_" in os.path.basename(p)),
        None,
    )
    fastqc_r2 = next(
        (os.path.relpath(p, base) for p in all_htmls if "_R2_" in os.path.basename(p)),
        None,
    )

    return fastqc_r1, fastqc_r2


# ==============================================================
# Main JSON Building Function
# ==============================================================
def build_nipt_json(
    analysis_dir,
    output_dir,
    ref_dir,
    sample_name,
    fetus_gender,
    age,
    version,
    target_bed_dir,
    config,
):
    """Build complete NIPT JSON structure matching output.json"""

    # Initialize output structure
    output = {
        "OrderInfo": {},
        APPID: {
            "algorithm_version": version,
            "Summary": {},
            "review": {
                "reviewer1": {
                    "Trisomy_result": "",
                    "Trisomy_comment": None,
                    "MD_result": None,
                    "MD_comment": None,
                    "saved": False,
                    "username": None,
                    "name": None
                },
                "reviewer2": {
                    "Trisomy_result": "",
                    "Trisomy_comment": None,
                    "MD_result": None,
                    "MD_comment": None,
                    "saved": False,
                    "username": None,
                    "name": None
                },
            },
            "lab_test": {
                "sample_suitability": "Pass",
                "dna_quality": "Pass",
                "library_quality": "Pass",
                "ngs_data_quality": "Pass",
                "reference_material_test": "Pass",
            },
        },
    }

    # Define directory paths
    dirs = {
        "result": f"{analysis_dir}/{sample_name}/Output_Result/",
        "ff": f"{analysis_dir}/{sample_name}/Output_FF/",
        "wc_data": f"{output_dir}/{sample_name}/Output_WC/",
        "wcx_data": f"{output_dir}/{sample_name}/Output_WCX/",
    }

    # 1. Build final results table from actual data files
    final_results = build_final_results_table(analysis_dir, sample_name)

    # 250610 : trisomy_result, md_result should have "Disease" list. So, I changed the return.
    if final_results:
        output[APPID]["final_results"] = {
            "order_id": sample_name,
            "fetal_fraction_yff": f"{final_results['fetal_fraction_yff']}",
            "fetal_fraction_seqff": f"{final_results['fetal_fraction_seqff']}",
            "ff_ratio": str(final_results["ff_ratio"]),
            #"sample_bias": final_results["sample_bias_qc"],
            "fetal_gender": final_results["fetal_gender"],
            # "trisomy_result": "Low Risk" if final_results['final_trisomy_result'] == ["Low Risk"] else "High Risk",
            # "md_result": "Low Risk" if final_results['md_results'] == ["Low Risk"] else "High Risk"
        }

    # 2. Read trisomy results
    trisomy_results = build_trisomy_results(
        analysis_dir, sample_name, fetus_gender, age, target_bed_dir
    )
    output[APPID]["trisomy_results"] = trisomy_results

    # ----------------------------------------------
    # 250610 : Trisomy detected list added.
    # put final_trisomy_result data
    # high_risk_results = [f"{entry['disease_name']} {entry['result']}" for entry in trisomy_results if entry["result"] == "High Risk"]

    # ---------------------------------------------------------------
    # 250716 : 아래 코드들을 참고해서 PRIZM detection에도 사용하자.
    # 일단은 EZD 위주로 trisomy_results, reviewer1,2 result를 채우고
    # PRIZM은 무시하자.
    # ---------------------------------------------------------------
    high_risk_results = [
        f"{entry['disease_name']}"
        for entry in trisomy_results
        if entry["result"] == "High Risk"
    ]
    final_result_trisomy_output = (
        #high_risk_results if high_risk_results else ["Low Risk"]
        high_risk_results if high_risk_results else []
    )
    output[APPID]["final_results"]["trisomy_result"] = final_result_trisomy_output

    output[APPID]["review"]["reviewer1"]["Trisomy_result"] = (
        "High Risk" if high_risk_results else "Low Risk"
    )

    Trisoy_LowRisk_comment = "The results are consistent with a low risk pregnancy."
    if not high_risk_results:
        output[APPID]["review"]["reviewer1"]["Trisomy_comment"] = Trisoy_LowRisk_comment
    
    output[APPID]["review"]["reviewer2"]["Trisomy_result"] = (
        "High Risk" if high_risk_results else "Low Risk"
    )
    if not high_risk_results:
        output[APPID]["review"]["reviewer2"]["Trisomy_comment"] = Trisoy_LowRisk_comment
    
    # ----------------------------------------------

    # 3. Build trisomy details
    groups = ["orig", "fetus", "mom"]

    output[APPID]["trisomy_details"] = {}

    for group in groups:
        # File paths for each analysis method
        files = {
            f"{group}_ezd_plot": f"Output_EZD/{group}_EZD_grid.png",
            f"{group}_prizm_chr_plot": f"Output_PRIZM/{sample_name}_{group}_chromosome_line.png",
            f"{group}_prizm_10mb_plot": f"Output_PRIZM/{sample_name}_{group}_10mb_line.png",
            f"{group}_wc_plot": f"Output_WC/{sample_name}.wc.{group}_z.png",
            f"{group}_wc_result": f"Output_WC/{sample_name}.wc.{group}.report.txt",
            f"{group}_wcx_plot": f"Output_WCX/{sample_name}.wcx.{group}.png",
            f"{group}_wcx_result": f"Output_WCX/{sample_name}.wcx.{group}_aberrations.bed",
        }

        # Read chromosome analysis
        ezd_chr_file = f"{analysis_dir}/{sample_name}/Output_EZD/{group}/Trisomy_detect_result_{group}_with_SCA.tsv"
        prizm_chr_file = f"{analysis_dir}/{sample_name}/Output_PRIZM/{group}/{sample_name}_{group}.trisomy_detection.tsv"
        result_table = read_chromosome_analysis_from_ezd_prizm_detailed(
            ezd_chr_file, prizm_chr_file, fetus_gender
        )

        threshold_file = f"{ref_dir}/EZD/{group}/{group}_thresholds_new.tsv"
        threshold_data = read_threshold_data(threshold_file)

        # None 값들을 실제 threshold 값으로 교체
        if threshold_data:
            for chr_name, chr_data in result_table.items():
                if chr_name in threshold_data:
                    chr_data["Z-score threshold"] = threshold_data[chr_name][
                        "z_threshold"
                    ]  # None → "L<2.18~2.7<H"
                    chr_data["UAR threshold"] = threshold_data[chr_name][
                        "uar_threshold"
                    ]  # None → "L<8.45~8.50<H"

        output[APPID]["trisomy_details"][group] = {
            **files,
            "result_table": result_table,
        }

        # -----------------------------------------------------------------------
        # 250716 : Add PRIZM result to "trisomy_results" and "final_results"
        '''
        if group in {"orig", "fetus"}:
            chromosome_to_item = {
                "Chromosome 21": "T21",
                "Chromosome 18": "T18",
                "Chromosome 13": "T13",
                "Chromosome 9":  "T9",
                "Chromosome 16": "T16",
                "Chromosome 22": "T22"
            }

            # PRIZM에서는 일단 chrX, Y Disease는 보고하지 말자.
            excluded_chromosomes = {"Chromosome X", "Chromosome Y"}

            for chrom, result_info in prizm_results.items():
                if chrom in excluded_chromosomes:
                    continue  # X, Y는 무시

                if result_info.get("PRIZM Detection") == "High Risk":
                    # 주요 염색체(T21~T22)인지 확인
                    item_code = chromosome_to_item.get(chrom, "other")

                    for result in output[APPID]["trisomy_results"]:
                        if result["item"] == item_code:
                            result["result"] = "High Risk"
                            break 
        '''
        # -----------------------------------------------------------------------


    # 4. Build md_results (summary table)
    md_results_table = []

    # Add 8 common syndromes
    for disease, info in MD8_DISEASE_MAPPING.items():
        md_results_table.append(
            {
                "item": info["item"],
                "location": info["location"],
                "disease_name": disease,
                "result": "Low Risk",  # Default, will be updated if detected
            }
        )

    # Add other syndromes
    md_results_table.extend(
        [
            {
                "item": "other_md108",
                "location": "other",
                "disease_name": "other_md108",
                "result": "Low Risk",
            },
            {
                "item": "other_md320",
                "location": "other",
                "disease_name": "other_md320",
                "result": "Low Risk",
            },
            {
                "item": "other_md87",
                "location": "other",
                "disease_name": "other_md87",
                "result": "Low Risk",
            },
            {
                "item": "other_md141",
                "location": "other",
                "disease_name": "other_md141",
                "result": "Low Risk",
            },
        ]
    )

    output[APPID]["md_results"] = {"result_table": md_results_table}

    # 5. Build md_details (detailed results)
    logger.info("Calling build_md_details_section ...")
    md_details, md8_detected, other_detected, wcx_chr_list_total = (
        build_md_details_section(analysis_dir, sample_name, target_bed_dir)
    )
    output[APPID]["md_details"] = md_details

    logger.info(f"Detected MD8 diseases: {md8_detected}")
    logger.info(f"Detected other MD diseases: {other_detected}")
    logger.info(f"WCX detectedchromosome list: {wcx_chr_list_total}")

    # final_md_result update with Disease names
    detected_md_output = [disease for disease in md8_detected + other_detected]
    output[APPID]["final_results"]["md_result"] = detected_md_output

    output[APPID]["review"]["reviewer1"]["MD_result"] = (
        "High Risk" if len(detected_md_output) > 0 else "Low Risk"
    )

    MD_LowRisk_comment = "The results are consistent with low risk of microdeletion/duplications in the regions of interest."
    if len(detected_md_output) == 0:
        output[APPID]["review"]["reviewer1"]["MD_comment"] = MD_LowRisk_comment

    output[APPID]["review"]["reviewer2"]["MD_result"] = (
        "High Risk" if len(detected_md_output) > 0 else "Low Risk"
    )
    if len(detected_md_output) == 0:
        output[APPID]["review"]["reviewer2"]["MD_comment"] = MD_LowRisk_comment

    # md_results_table update. It's initialized as "Low Risk" for all diseases.
    for row in md_results_table:
        if row["disease_name"] in md8_detected or row["item"] in other_detected:
            row["result"] = "High Risk"

    # 6. Build quality control section
    Final_QC_result = "PASS"
    No_Call_reason = [] 
    try:
        qc_file = f"{analysis_dir}/{sample_name}/Output_QC/{sample_name}.qc.filter.txt"
        logger.info(f"Looking for QC file: {qc_file}")

        if not os.path.exists(qc_file):
            Final_QC_result = "FAIL"
            No_Call_reason.append(f"QC file not found: {qc_file}")
            logger.warning(f"QC file not found: {qc_file}")
            logger.warning("Using default QC values")

            # 기본 QC 데이터 생성
            sequencing_metrics = {
                "total_reads": {
                    "value": 0,
                    "status": "UNKNOWN",
                    "unit": "reads",
                    "threshold": ">10M",
                },
                "mapped_reads": {
                    "value": 0,
                    "status": "UNKNOWN",
                    "unit": "reads",
                    "threshold": ">85%",
                },
                "mapping_rate": {
                    "value": 0,
                    "status": "UNKNOWN",
                    "unit": "%",
                    "threshold": ">85%",
                },
                "duplicated_reads": {
                    "value": 0,
                    "status": "UNKNOWN",
                    "unit": "reads",
                    "threshold": "<30%",
                },
                "duplication_rate": {
                    "value": 0,
                    "status": "UNKNOWN",
                    "unit": "%",
                    "threshold": "<40%",
                },
                "mean_mapping_quality": {
                    "value": 0,
                    "status": "UNKNOWN",
                    "unit": "score",
                    "threshold": ">20",
                },
                "mean_coverage": {
                    "value": 0,
                    "status": "UNKNOWN",
                    "unit": "X",
                    "threshold": ">0.1X",
                },
                "gc_content": {
                    "value": 0,
                    "status": "UNKNOWN",
                    "unit": "%",
                    "threshold": "33-55%",
                },
            }
        else:
            qc_df = read_qc_data(qc_file)

            # Sequencing metrics
            sequencing_metrics = {}

            if isinstance(qc_df, pd.DataFrame) and not qc_df.empty:
                qc_mapping = {
                    "number_of_reads": "total_reads",
                    "number_of_mapped_reads": "mapped_reads",
                    "mapping_rate": "mapping_rate",
                    "number_of_duplicated_reads": "duplicated_reads",  # 실제 파일명에 맞춤
                    "duplication_rate": "duplication_rate",
                    "mean_mapping_quality": "mean_mapping_quality",
                    "mean_coverageData": "mean_coverage",
                    "GC_content": "gc_content",
                }

                for qc_key, json_key in qc_mapping.items():
                    try:
                        if qc_key in qc_df.index:
                            row = qc_df.loc[qc_key]
                            # Parse value - remove units and convert to appropriate type
                            value_str = str(row.value).replace(" ", "").replace("%", "")
                            try:
                                if "." in value_str:
                                    value = float(value_str)
                                else:
                                    value = int(value_str)
                            except:
                                value = row.value

                            # Determine unit based on original value
                            unit = ""
                            if "%" in str(row.value):
                                unit = "%"
                            elif "reads" in json_key:
                                unit = "reads"
                            elif "quality" in json_key:
                                unit = "score"
                            elif "coverage" in json_key:
                                unit = "X"

                            # Determine threshold
                            threshold = ""
                            if json_key == "total_reads":
                                threshold = ">10M"
                            elif "mapping_rate" in json_key:
                                threshold = ">85%"
                            elif "duplication_rate" in json_key:
                                threshold = "<40%"
                            elif "quality" in json_key:
                                threshold = ">20"
                            elif "coverage" in json_key:
                                threshold = ">0.1X"
                            elif "gc_content" in json_key:
                                threshold = "33-55%"

                            sequencing_metrics[json_key] = {
                                "value": value,
                                "status": row.status,
                                "unit": unit,
                                "threshold": threshold,
                            }

                            if row.status != "PASS":
                                Final_QC_result = "FAIL"

                                fail_reason_map = {
                                    "total_reads": "Total reads are less than 10 Million",
                                    "mapped_reads": "Mapped reads are less than 9.5 Million",
                                    "mapping_rate": "Mapping rate is below 85%",
                                    "duplication_rate": "Duplication rate is over 40%",
                                    "mean_mapping_quality": "Quality is below 20",
                                    "mean_coverage": "Coverage is below 0.1X",
                                    "gc_content": "GC content % is out of range (33-55%)",
                                }

                                fail_reason = fail_reason_map.get(json_key)
                                if fail_reason and fail_reason not in No_Call_reason:
                                    No_Call_reason.append(fail_reason)


                    except Exception as e:
                        logger.error(f"Error processing QC metric {qc_key}: {e}")
                        continue

            else:
                logger.warning("QC DataFrame is empty, using default values")
                sequencing_metrics = {}

        # Analysis QC
        analysis_qc = {}

        qc_config = config.get("QC", {})
        orig_biqc = qc_config.get("orig_biqc", 4.0)
        yff_threshold = qc_config.get("YFF", 4.0)
        seqff_threshold = qc_config.get("seqFF", 4.0)

        try:
            prizm_qc_file = f"{analysis_dir}/{sample_name}/Output_PRIZM/orig/{sample_name}.of_orig.prizm.qc.txt"
            with open(prizm_qc_file, "r") as qc_f:
                parts = qc_f.read().strip().split()
                sample_bias_val = round(float(parts[0]), 3)
                sample_bias_status = parts[1]

            if sample_bias_status != "PASS":
                Final_QC_result = "FAIL"
                No_Call_reason.append(f"Sample bias is greater than {orig_biqc}")

            logger.info(f"sample_bias_status : {sample_bias_status}")

            # Get fetal gender to determine which FF to use for QC
            fetal_gender = final_results.get("fetal_gender", "Unknown")
            logger.info(f"fetal_gender : {fetal_gender}")

            yff_val = final_results.get("fetal_fraction_yff", "N/A")
            if yff_val != "N/A":
                yff_val_float = float(yff_val)
                yff_status = "PASS" if yff_val_float >= yff_threshold else "FAIL"
            else:
                yff_status = "PASS"  # 예외로 인정

            seqff_val = round(float(final_results["fetal_fraction_seqff"]), 3)
            seqff_status = "PASS" if seqff_val >= seqff_threshold else "FAIL"

            # Apply FF threshold based on gender
            # Male: Use YFF for QC (SeqFF is not reliable for males)
            # Female: Use SeqFF for QC (YFF is N/A for females)
            if fetal_gender == "Male":
                # For Male, only check YFF
                if yff_status != "PASS":
                    Final_QC_result = "FAIL"
                    No_Call_reason.append(f"Low Fetal Fraction (<{yff_threshold}%)")
                logger.info(f"Male sample - using YFF for QC: {yff_status}")
            else:
                # For Female (or Unknown), check SeqFF
                if seqff_status != "PASS":
                    Final_QC_result = "FAIL"
                    No_Call_reason.append(f"Low Fetal Fraction (<{seqff_threshold}%)")
                logger.info(f"Female sample - using SeqFF for QC: {seqff_status}")

            logger.info(f"yff_status : {yff_status}")
            logger.info(f"seqff_status : {seqff_status}")

            ff_ratio_threshold = qc_config.get("FF_Ratio", 2.5)
            ff_ratio_val = float(final_results["ff_ratio"])
            ff_ratio_status = "PASS" if ff_ratio_val < ff_ratio_threshold else "FAIL"

            if ff_ratio_status != "PASS":
                Final_QC_result = "FAIL"
                No_Call_reason.append(f"FF_ratio is greater than {ff_ratio_threshold}")

            logger.info(f"ff_ratio_status : {ff_ratio_status}")

            analysis_qc = {
                "fetal_fraction_yff": {
                    "value": yff_val,
                    "unit": "%",
                    "status": yff_status,
                    "threshold": f">{yff_threshold}%",
                },
                "fetal_fraction_seqff": {
                    "value": seqff_val,
                    "unit": "%",
                    "status": seqff_status,
                    "threshold": f">{seqff_threshold}%",
                },
                "ff_ratio": {
                    "value": ff_ratio_val,
                    "unit": "%",
                    "status": ff_ratio_status,
                    "threshold": f"<{ff_ratio_threshold}",
                },
                "sample_bias_qc": {
                    "value": sample_bias_val,
                    "status": sample_bias_status,
                    "threshold": f"<{orig_biqc}",
                },
            }

        except Exception as e:
            logger.error(f"Error building analysis QC: {e}")
            analysis_qc = {}

        qc_dir = os.path.join(analysis_dir, sample_name, "Output_QC")
        fastqc_r1_report, fastqc_r2_report = find_fastqc_reports(qc_dir)

        output[APPID]["quality_control"] = {
            "sequencing_metrics": sequencing_metrics,
            "analysis_qc": analysis_qc,
            "qc_files": {
                "Fastqc_R1_report": f"Output_QC/{fastqc_r1_report}",
                "Fastqc_R2_report": f"Output_QC/{fastqc_r2_report}",
                # "qualimap_report": f"Output_QC/{sample_name}.Qualimap.zip",
                "Qualimap_report": "Output_QC/qualimapReport.html",
            },
        }

    except Exception as e:
        logger.error(f"Error building quality control section: {e}")
        # 완전히 기본값으로 설정
        output[APPID]["quality_control"] = {
            "sequencing_metrics": {},
            "analysis_qc": {},
            "qc_files": {
                "qualimap_report": f"Output_QC/{sample_name}.Qualimap.zip",
                "qc_summary_report": "Output_QC/qualimapReport.html",
            },
        }

    # 250713 : added to show QC Results in Final Results Summary section
    # 250809 : MD_result, MD_comment are also set as No call
    output[APPID]["final_results"]["QC_result"] = Final_QC_result

    if Final_QC_result != "PASS":
        output[APPID]["review"]["reviewer1"]["Trisomy_result"] = "No call"
        output[APPID]["review"]["reviewer1"]["Trisomy_comment"] = ", ".join(No_Call_reason)
        output[APPID]["review"]["reviewer1"]["MD_result"] = "No call"
        output[APPID]["review"]["reviewer1"]["MD_comment"] = ", ".join(No_Call_reason)

        output[APPID]["review"]["reviewer2"]["Trisomy_result"] = "No call"
        output[APPID]["review"]["reviewer2"]["Trisomy_comment"] = ", ".join(No_Call_reason)
        output[APPID]["review"]["reviewer2"]["MD_result"] = "No call"
        output[APPID]["review"]["reviewer2"]["MD_comment"] = ", ".join(No_Call_reason)

    # 7. Add QC section outside NIPT (matching output.json structure)
    #output["quality_control"] = output[APPID]["quality_control"]

    # 8. Build S3 upload files
    # I don't need to put this data. Just make tar file
    """
    output["S3_upload_files"] = {
        #"qualimap_report_button": f"Output_QC/{sample_name}.Qualimap.zip",
        "fastqc_r1_report_button": f"Output_QC/{fastqc_r1_report}",
        "fastqc_r2_report_button": f"Output_QC/{fastqc_r2_report}",
        "qualimap_report_button": f"Output_QC/qualimapReport.html",

        "original_ezd_image": f"Output_EZD/orig_EZD_grid.png",
        "original_prizm_chromosome_image": f"Output_PRIZM/{sample_name}_orig_chromosome_line.png",
        "original_prizm_10mb_image": f"Output_PRIZM/{sample_name}_orig_10mb_line.png",
        "original_wisecondor_image": f"Output_WC/{sample_name}.wc.orig_z.png",
        "original_wisecondor_txt": f"Output_WC/{sample_name}.wc.orig.report.txt",
        "original_wisecondorx_image": f"Output_WCX/{sample_name}.wcx.orig.png",
        "original_wisecondorx_txt": f"Output_WCX/{sample_name}.wcx.orig_aberrations.bed",

        "fetus_ezd_image": f"Output_EZD/fetus_EZD_grid.png",
        "fetus_prizm_chromosome_image": f"Output_PRIZM/{sample_name}_fetus_chromosome_line.png",
        "fetus_prizm_10mb_image": f"Output_PRIZM/{sample_name}_fetus_10mb_line.png",
        "fetus_wisecondor_image": f"Output_WC/{sample_name}.wc.fetus_z.png",
        "fetus_wisecondor_txt": f"Output_WC/{sample_name}.wc.fetus.report.txt",
        "fetus_wisecondorx_image": f"Output_WCX/{sample_name}.wcx.fetus.png",
        "fetus_wisecondorx_txt": f"Output_WCX/{sample_name}.wcx.fetus_aberrations.bed",

        "mom_ezd_image": f"Output_EZD/mom_EZD_grid.png",
        "mom_prizm_chromosome_image": f"Output_PRIZM/{sample_name}_mom_chromosome_line.png",
        "mom_prizm_10mb_image": f"Output_PRIZM/{sample_name}_mom_10mb_line.png",
        "mom_wisecondor_image": f"Output_WC/{sample_name}.wc.mom_z.png",
        "mom_wisecondor_txt": f"Output_WC/{sample_name}.wc.mom.report.txt",
        "mom_wisecondorx_image": f"Output_WCX/{sample_name}.wcx.mom.png",
        "mom_wisecondorx_txt": f"Output_WCX/{sample_name}.wcx.mom_aberrations.bed",
    }
    """

    json_output_path = f"{output_dir}/{sample_name}/{sample_name}.json"

    try:
        # 출력 디렉토리 생성
        os.makedirs(os.path.dirname(json_output_path), exist_ok=True)

        # JSON 파일 저장
        with open(json_output_path, "w") as json_file:
            json.dump(output, json_file, indent=2)

        # 파일 생성 확인
        if os.path.exists(json_output_path):
            file_size = os.path.getsize(json_output_path)
            logger.info(
                f"NIPT JSON successfully saved: {json_output_path} ({file_size} bytes)"
            )

            # 파일 경로 반환 (JSON 데이터가 아닌!)
            return json_output_path, wcx_chr_list_total
        else:
            logger.error(f"JSON file was not created: {json_output_path}")
            return None, None

    except Exception as e:
        logger.error(f"Failed to save JSON file to {json_output_path}: {e}")
        return None, None


# ==============================================================
# Main Function
# ==============================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
    )

    pgroup = parser.add_argument_group("Input")
    pgroup.add_argument(
        "-analysis_dir",
        dest="analysis_dir",
        type=lambda x: file_check(parser, x),
        help="Analysis directory",
    )
    pgroup.add_argument("-sample_name", dest="sample_name", help="Sample name")
    # pgroup.add_argument('-config_file', dest='config_file', help='Configuration file')
    pgroup.add_argument(
        "-output_dir",
        dest="output_dir",
        type=lambda x: file_check(parser, x),
        help="Output directory",
    )
    pgroup.add_argument(
        "-target_bed_dir",
        dest="target_bed_dir",
        type=lambda x: file_check(parser, x),
        help="Target bed directory",
    )
    pgroup.add_argument("-version", dest="version", help="NIPT version (V1.0)")
    pgroup.add_argument(
        "-ref_dir",
        dest="ref_dir",
        help="Reference directory (for EZD thresholds etc.)",
        default="",
    )
    pgroup.add_argument(
        "-age",
        dest="age",
        type=int,
        default=35,
        help="Maternal age (default: 35)",
    )
    pgroup.add_argument(
        "-config_file",
        dest="config_file",
        help="Pipeline config JSON file",
        default=None,
    )

    ogroup = parser.add_argument_group("Options")
    ogroup.add_argument(
        "-h", "--help", action="help", help="show this help message and exit"
    )
    ogroup.add_argument(
        "--debug",
        dest="debug",
        action="store_true",
        default=False,
        help=argparse.SUPPRESS,
    )

    args = parser.parse_args()
    if args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    # Generate JSON
    # data = build_nipt_json(args.analysis_dir, args.output_dir, args.sample_name, args.version, args.config_file, args.target_bed_dir)

    # Derive fetus_gender from the published gender.txt produced by FF_GENDER_WORKFLOW
    _gender_txt = os.path.join(
        args.analysis_dir, args.sample_name, "Output_FF", f"{args.sample_name}.gender.txt"
    )
    fetus_gender = "F"  # safe default
    if os.path.isfile(_gender_txt):
        for _line in open(_gender_txt):
            _parts = _line.strip().split()
            if len(_parts) >= 2 and _parts[0].lower() == "final_gender":
                _g = _parts[1].upper()
                if _g in ("XY", "MALE", "M"):
                    fetus_gender = "M"
                else:
                    fetus_gender = "F"
                break

    # Load pipeline config dict (optional)
    _config = {}
    if args.config_file and os.path.isfile(args.config_file):
        import json as _json
        with open(args.config_file) as _f:
            _config = _json.load(_f)

    # build_nipt_json() saves the file internally and returns (json_path, wcx_chr_dict).
    # Do NOT write again — that would overwrite the real JSON with the path string.
    json_path, wcx_chr_dict = build_nipt_json(
        args.analysis_dir,
        args.output_dir,
        args.ref_dir,
        args.sample_name,
        fetus_gender,
        args.age,
        args.version,
        args.target_bed_dir,
        _config,
    )

    if not json_path:
        raise RuntimeError("build_nipt_json() failed to produce a JSON file")

    logger.info(f"NIPT JSON has been saved at {json_path}")
