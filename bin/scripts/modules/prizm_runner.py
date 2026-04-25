#!/usr/bin/env python

"""
---------------------------------------------
PRIZM (Prenatal Risk Z-score Matrix) - Complete Analysis Pipeline
Z-score calculation for trisomy detection with integrated visualization
Includes chromosome, 10mb, and 10mb_all Z-score calculations with decision logic
Author: Hyukjung Kwon
---------------------------------------------
"""

import argparse
import logging
import os
import re
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.patches as patches
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)

# Constants
DATA_DIR = os.environ.get("DATA_DIR", "/Work/NIPT/data")
ANALYSIS_DIR = os.environ.get("ANALYSIS_DIR", "/Work/NIPT/analysis")

# ============================================================================
# DATACLASSES
# ============================================================================


@dataclass
class PRIZMResult:
    """Container for PRIZM analysis results"""

    zscore_chr: pd.DataFrame
    zscore_10mb: pd.DataFrame
    zscore_10mb_all: pd.DataFrame
    chr_index_dict: Dict[int, str]
    mb10_index_dict: Dict[int, str]
    row_count: int


@dataclass
class PRIZMConfig:
    """Configuration for PRIZM analysis parameters"""

    qc_cutoff: float = 3.0
    single_output: bool = False


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================


def setup_logging(debug: bool = False) -> logging.Logger:
    """Setup logging configuration"""
    logger = logging.getLogger(__name__)
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter(
        "%(asctime)s -%(lineno)s- [%(levelname) 5s] - %(message)s"
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG if debug else logging.INFO)
    return logger


def natural_sort_key(text: str) -> List:
    """Natural sorting key for chromosome names"""
    return [int(c) if c.isdigit() else c for c in re.split(r"(\d+)", text)]


def natural_keys(text):
    """Natural sorting key for chromosome names (MultiZ style)"""

    def atoi(text):
        return int(text) if text.isdigit() else text

    return [atoi(c) for c in re.split(r"(\d+)", text)]


def file_exists(parser: argparse.ArgumentParser, filepath: str) -> str:
    """Validate file existence"""
    if not os.path.exists(filepath):
        parser.error(f"File {filepath} does not exist!")
    return str(filepath)


def create_output_dir(dir_name="plots"):
    """Create output directory for plots"""
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    return dir_name


def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z


# ============================================================================
# CORE Z-SCORE CALCULATION FUNCTIONS
# ============================================================================


def load_count_data(count_file: str) -> pd.DataFrame:
    """Load and preprocess 10mb count data"""
    logger.info(f"Loading count data: {count_file}")

    count_data = pd.read_csv(count_file, sep="\t")
    count_data = count_data[["chr", "start", "reads"]].rename(
        columns={"start": "bin", "reads": "count"}
    )
    count_data["bin"] = count_data["bin"] / 10000000

    return count_data


def prepare_chromosome_data(
    count_data: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Prepare chromosome-level and 10mb bin-level data"""

    # 1. Chromosome-level aggregation
    chr_summary = count_data.groupby("chr")["count"].sum().to_frame()
    sorted_index = sorted(chr_summary.index, key=natural_sort_key)
    chr_summary = chr_summary.loc[sorted_index]

    # Calculate ratios for chromosome data
    total_count = chr_summary["count"].sum()
    chr_summary["ratio"] = (chr_summary["count"] / total_count) * 100

    # 2. Prepare 10mb bin data
    # Separate special chromosomes and regular chromosomes
    special_chroms = [
        "chr9",
        "chr13",
        "chr16",
        "chr18",
        "chr21",
        "chr22",
        "chrX",
        "chrY",
    ]
    regular_chroms = [chr for chr in sorted_index if chr not in special_chroms]

    # Regular chromosomes (aggregated to chromosome level)
    regular_data = chr_summary.loc[regular_chroms].copy()
    regular_data.insert(0, "bin", 0)  # Add bin column

    # Special chromosomes (keep 10mb resolution)
    special_data = count_data.set_index("chr").loc[special_chroms].reset_index()

    # Combine data
    mb10_data = pd.concat([regular_data.reset_index(), special_data], ignore_index=True)

    # Sort by chromosome order
    mb10_data = mb10_data.set_index("chr").loc[sorted_index].reset_index()

    # Calculate ratios for 10mb data
    total_count_10mb = mb10_data["count"].sum()
    mb10_data["ratio"] = (mb10_data["count"] / total_count_10mb) * 100

    return chr_summary, mb10_data


def calculate_normalized_matrix(data_df: pd.DataFrame) -> np.ndarray:
    """Calculate normalized ratio matrix"""
    n_samples = len(data_df)
    normalized_matrix = np.zeros((n_samples, n_samples))

    ratios = data_df["ratio"].values

    for i in range(n_samples):
        for j in range(n_samples):
            if ratios[j] != 0:
                normalized_matrix[i, j] = ratios[i] / ratios[j]
            else:
                normalized_matrix[i, j] = 0

    return normalized_matrix


def calculate_zscore_from_normalized(
    normalized_matrix: np.ndarray, mean_file: str, sd_file: str
) -> np.ndarray:
    """Calculate Z-score matrix from normalized matrix"""

    mean_data = pd.read_csv(mean_file, sep="\t", header=None).values
    sd_data = pd.read_csv(sd_file, sep="\t", header=None).values

    n_samples = normalized_matrix.shape[0]
    zscore_matrix = np.zeros_like(normalized_matrix)

    for i in range(n_samples):
        for j in range(n_samples):
            if sd_data[i, j] != 0:
                zscore_matrix[i, j] = (
                    normalized_matrix[i, j] - mean_data[i, j]
                ) / sd_data[i, j]
            else:
                zscore_matrix[i, j] = 0.0

    return zscore_matrix


def create_zscore_dataframe(
    zscore_matrix: np.ndarray, index_names: List[str]
) -> pd.DataFrame:
    """Create labeled Z-score DataFrame"""

    return pd.DataFrame(zscore_matrix, index=index_names, columns=index_names)


# ============================================================================
# 10MB_ALL Z-SCORE CALCULATION (FROM MultiZ)
# ============================================================================


def getNormalized10mbData(filename, autosomal, key_start=0):
    """Get normalized 10mb data (from MultiZ multiZ_core_all_v4_10mball.py)"""

    df = pd.read_csv(filename, sep="\t")
    df = df[["chr", "start", "reads"]].rename(
        columns={"start": "bin", "reads": "count"}
    )
    df["bin"] = df["bin"].apply(lambda x: (x - 1) / 10000000)
    df.bin = df.bin.astype(int)

    # HMMcopy.R makes some different order, so I need to sort them again
    df = df.set_index("chr").loc[
        [
            "chr1",
            "chr2",
            "chr3",
            "chr4",
            "chr5",
            "chr6",
            "chr7",
            "chr8",
            "chr9",
            "chr10",
            "chr11",
            "chr12",
            "chr13",
            "chr14",
            "chr15",
            "chr16",
            "chr17",
            "chr18",
            "chr19",
            "chr20",
            "chr21",
            "chr22",
            "chrX",
            "chrY",
        ]
    ]
    df.reset_index(level=0, inplace=True)

    # sum function doesn't follow the original order, so I need to sort it too
    count_data_allsum = df.groupby(by=["chr"])["count"].sum()
    sum_df = count_data_allsum.to_frame()
    sorted_index = sorted(sum_df.index, key=natural_keys)
    sum_df = sum_df.loc[sorted_index]

    # df : 10mb dataframe
    # sum_df : chromosome dataframe
    sum_df = sum_df.iloc[0:22, :]
    count_sum = sum_df["count"].sum()

    # Calculate normalized ratios for chromosomes
    norm_data = [
        (sum_df.loc[sum_df.index[i], "count"] / float(count_sum)) * 100
        for i in range(len(sum_df))
    ]
    sum_df.insert(1, "ratio", norm_data)

    if autosomal == True:
        df = df.loc[~df["chr"].isin(["chrX", "chrY"])]
    else:
        df = df.loc[df["chr"].isin(["chrX", "chrY"])]

    norm_10mb_data = [
        (df.loc[df.index[i], "count"] / float(count_sum)) * 100 for i in range(len(df))
    ]
    df.insert(2, "ratio", norm_10mb_data)

    normalized_dict = {
        j + key_start: [
            df.loc[df.index[j], "ratio"] / sum_df.loc[sum_df.index[i], "ratio"]
            if sum_df.loc[sum_df.index[i], "ratio"] != 0.0
            and df.loc[df.index[j], "ratio"] != 0.0
            else 0
            for i in range(len(sum_df))
        ]
        for j in range(len(df))
    }

    return normalized_dict


def calc_zscore_10mb_all(count_file, mean_df, sd_df):
    """Calculate 10mb_all zscore (from MultiZ multiZ_core_all_v4_10mball.py)"""

    logger.debug("calc_zscore_10mb_all run")

    # Calculate normalized data for autosomal and sex chromosomes
    norm_dict = getNormalized10mbData(count_file, True)
    norm_XY_dict = getNormalized10mbData(
        count_file, False, key_start=len(norm_dict.keys())
    )

    logger.debug("norm_dict length %d" % len(norm_dict))
    logger.debug("norm_XY_dict length %d" % len(norm_XY_dict))

    merged_norm_dict = merge_two_dicts(norm_dict, norm_XY_dict)
    normalized_df = pd.DataFrame(merged_norm_dict)

    logger.debug("Merged Normalized DF shape : %s " % str(normalized_df.shape))

    # Calculate Z-scores
    z = {
        j: [
            (normalized_df.iloc[i, j] - mean_df.iloc[i, j]) / sd_df.iloc[i, j]
            if sd_df.iloc[i, j] != 0.0
            else 0.0
            for i in range(len(mean_df))
        ]
        for j in range(len(mean_df.iloc[0]))
    }

    zscore_10mb_all = pd.DataFrame(z)

    return zscore_10mb_all


# ============================================================================
# MAIN CALCULATION PIPELINE
# ============================================================================


def calculate_prizm_zscores(
    count_file: str,
    mean_file: str,
    sd_file: str,
    mean_10mb_file: str,
    sd_10mb_file: str,
    mean_10mb_all_file: str,
    sd_10mb_all_file: str,
) -> PRIZMResult:
    """Main PRIZM Z-score calculation pipeline"""

    # 1. Load count data
    count_data = load_count_data(count_file)

    # 2. Prepare chromosome and 10mb data
    chr_summary, mb10_data = prepare_chromosome_data(count_data)

    # 3. Calculate chromosome-level Z-scores
    logger.info("Calculating chromosome-level Z-scores")
    chr_normalized = calculate_normalized_matrix(chr_summary)
    chr_zscore_matrix = calculate_zscore_from_normalized(
        chr_normalized, mean_file, sd_file
    )

    # Create chromosome index mapping
    chr_index_dict = {i: f"chr{i + 1}" for i in range(22)}
    chr_index_dict.update({22: "chrX", 23: "chrY"})

    chr_zscore_df = create_zscore_dataframe(
        chr_zscore_matrix, list(chr_index_dict.values())
    )

    # 4. Calculate 10mb-level Z-scores
    logger.info("Calculating 10mb-level Z-scores")
    mb10_normalized = calculate_normalized_matrix(mb10_data)
    mb10_zscore_matrix = calculate_zscore_from_normalized(
        mb10_normalized, mean_10mb_file, sd_10mb_file
    )

    # Create 10mb index mapping
    mb10_index_dict = {i: mb10_data.iloc[i]["chr"] for i in range(len(mb10_data))}

    mb10_zscore_df = create_zscore_dataframe(
        mb10_zscore_matrix, list(mb10_index_dict.values())
    )

    # 5. Calculate 10mb_all Z-scores
    logger.info("Calculating 10mb_all Z-scores")
    mean_10mb_all_df = pd.read_csv(mean_10mb_all_file, sep="\t", header=None)
    sd_10mb_all_df = pd.read_csv(sd_10mb_all_file, sep="\t", header=None)

    zscore_10mb_all_df = calc_zscore_10mb_all(
        count_file, mean_10mb_all_df, sd_10mb_all_df
    )

    logger.info("PRIZM Z-score calculation completed successfully")

    return PRIZMResult(
        zscore_chr=chr_zscore_df,
        zscore_10mb=mb10_zscore_df,
        zscore_10mb_all=zscore_10mb_all_df,
        chr_index_dict=chr_index_dict,
        mb10_index_dict=mb10_index_dict,
        row_count=len(mb10_data),
    )


# ============================================================================
# STATISTICAL TRISOMY DETECTION (NEW)
# ============================================================================


def statistical_trisomy_detection(
    zscore_chr: pd.DataFrame,
    target_chromosomes: List[str] = None,
    confidence_level: float = 0.95,
    effect_size_threshold: float = 0.8,
    min_high_scores: int = 15,
) -> Dict[str, Dict]:
    """
    Statistical trisomy detection based on Z-score distribution analysis

    Args:
        zscore_chr: 24x24 chromosome Z-score matrix
        target_chromosomes: Chromosomes to test (default: all autosomal + sex chromosomes)
        confidence_level: Statistical confidence level for detection
        effect_size_threshold: Minimum Cohen's d effect size for detection
        min_high_scores: Minimum number of high Z-scores required

    Returns:
        Dict with detection results for each chromosome
    """

    if target_chromosomes is None:
        # Analyze all chromosomes by default
        target_chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

    results = {}

    # Get all chromosomes except sex chromosomes for baseline comparison
    autosomal_chroms = [f"chr{i}" for i in range(1, 23)]

    for target_chr in target_chromosomes:
        if target_chr not in zscore_chr.index:
            logger.warning(
                f"Target chromosome {target_chr} not found in Z-score matrix"
            )
            continue

        # logger.info(f"=== Analyzing {target_chr} ===")

        # Get Z-scores for target chromosome (excluding self-comparison)
        target_scores = []
        comparison_chroms = []

        for other_chr in zscore_chr.columns:
            if other_chr != target_chr:  # Exclude self-comparison (which should be ~0)
                score = zscore_chr.loc[target_chr, other_chr]
                target_scores.append(score)
                comparison_chroms.append(other_chr)

        target_scores = np.array(target_scores)

        # Get baseline Z-scores from other autosomal chromosomes (for comparison)
        baseline_scores = []

        for baseline_chr in autosomal_chroms:
            if baseline_chr != target_chr and baseline_chr in zscore_chr.index:
                for other_chr in zscore_chr.columns:
                    if other_chr != baseline_chr:  # Exclude self-comparison
                        score = zscore_chr.loc[baseline_chr, other_chr]
                        baseline_scores.append(score)

        baseline_scores = np.array(baseline_scores)

        # Statistical analysis
        analysis_result = analyze_chromosome_distribution(
            target_chr,
            target_scores,
            baseline_scores,
            comparison_chroms,
            confidence_level,
            effect_size_threshold,
            min_high_scores,
        )

        results[target_chr] = analysis_result

        # Log results for key chromosomes
        if target_chr in ["chr13", "chr18", "chr21"]:
            decision = analysis_result["decision"]
            confidence = analysis_result["confidence"]
            effect_size = analysis_result["effect_size"]

            logger.info(f"{target_chr}: {decision}")
            logger.info(f"  - Confidence: {confidence:.1%}")
            logger.info(f"  - Effect size (Cohen's d): {effect_size:.3f}")
            logger.info(
                f"  - High Z-scores: {analysis_result['high_count']}/{len(target_scores)}"
            )
            logger.info(f"  - Mean Z-score: {analysis_result['mean_zscore']:.3f}")

    return results


def analyze_chromosome_distribution(
    target_chr: str,
    target_scores: np.ndarray,
    baseline_scores: np.ndarray,
    comparison_chroms: List[str],
    confidence_level: float,
    effect_size_threshold: float,
    min_high_scores: int,
) -> Dict:
    """
    Analyze the statistical distribution of Z-scores for a specific chromosome
    """

    # Basic statistics
    target_mean = np.mean(target_scores)
    target_std = np.std(target_scores)
    baseline_mean = np.mean(baseline_scores)
    baseline_std = np.std(baseline_scores)

    # Statistical tests
    from scipy import stats

    # 1. Two-sample t-test: Are target scores significantly higher than baseline?
    t_stat, p_value = stats.ttest_ind(
        target_scores, baseline_scores, alternative="greater"
    )

    # 2. Effect size (Cohen's d): How large is the difference?
    pooled_std = np.sqrt(
        (
            (len(target_scores) - 1) * target_std**2
            + (len(baseline_scores) - 1) * baseline_std**2
        )
        / (len(target_scores) + len(baseline_scores) - 2)
    )

    cohens_d = (target_mean - baseline_mean) / pooled_std if pooled_std > 0 else 0

    # 3. Count high Z-scores (adaptive threshold based on baseline distribution)
    baseline_q95 = np.percentile(baseline_scores, 95)  # 95th percentile of baseline
    high_scores = target_scores[target_scores > baseline_q95]
    high_count = len(high_scores)

    # 4. Outlier analysis: How many target scores are outliers?
    baseline_q75 = np.percentile(baseline_scores, 75)
    baseline_iqr = np.percentile(baseline_scores, 75) - np.percentile(
        baseline_scores, 25
    )
    outlier_threshold = baseline_q75 + 1.5 * baseline_iqr
    outlier_count = np.sum(target_scores > outlier_threshold)

    # 5. Consistency check: Are the high scores consistent across comparisons?
    target_sorted_indices = np.argsort(target_scores)[::-1]  # Descending order
    top_chromosomes = [
        comparison_chroms[i]
        for i in target_sorted_indices[: min(10, len(target_sorted_indices))]
    ]

    # Decision logic
    decision, confidence = make_statistical_decision(
        p_value,
        cohens_d,
        high_count,
        outlier_count,
        len(target_scores),
        confidence_level,
        effect_size_threshold,
        min_high_scores,
    )

    # Detailed results
    result = {
        "decision": decision,
        "confidence": 1 - p_value
        if p_value < 0.5
        else confidence,  # Statistical confidence
        "p_value": p_value,
        "effect_size": cohens_d,
        "mean_zscore": target_mean,
        "std_zscore": target_std,
        "baseline_mean": baseline_mean,
        "baseline_std": baseline_std,
        "high_count": high_count,
        "outlier_count": outlier_count,
        "total_comparisons": len(target_scores),
        "high_score_ratio": high_count / len(target_scores),
        "adaptive_threshold": baseline_q95,
        "outlier_threshold": outlier_threshold,
        "top_chromosomes": top_chromosomes[
            :5
        ],  # Top 5 chromosomes with highest Z-scores
        "statistical_summary": {
            "t_statistic": t_stat,
            "cohens_d_interpretation": interpret_effect_size(cohens_d),
            "baseline_95th_percentile": baseline_q95,
        },
    }

    return result


def make_statistical_decision(
    p_value: float,
    cohens_d: float,
    high_count: int,
    outlier_count: int,
    total_comparisons: int,
    confidence_level: float,
    effect_size_threshold: float,
    min_high_scores: int,
) -> Tuple[str, float]:
    """
    Make final decision based on multiple statistical criteria
    """

    alpha = 1 - confidence_level
    high_ratio = high_count / total_comparisons
    outlier_ratio = outlier_count / total_comparisons

    # Multiple criteria for robust detection
    criteria_met = 0
    total_criteria = 5

    # Criterion 1: Statistical significance
    if p_value < alpha:
        criteria_met += 1

    # Criterion 2: Large effect size
    if cohens_d > effect_size_threshold:
        criteria_met += 1

    # Criterion 3: Sufficient high scores
    if high_count >= min_high_scores:
        criteria_met += 1

    # Criterion 4: High proportion of elevated scores
    if high_ratio > 0.4:  # At least 40% of comparisons are high
        criteria_met += 1

    # Criterion 5: Presence of outliers
    if outlier_count >= 3:  # At least 3 outlier comparisons
        criteria_met += 1

    # Decision based on criteria met
    if criteria_met >= 4:
        decision = "Detected"
        confidence = min(0.99, 1 - p_value)
    elif criteria_met >= 3:
        decision = "Suspected"
        confidence = min(0.90, 1 - p_value)
    elif criteria_met >= 2:
        decision = "Possible"
        confidence = min(0.80, 1 - p_value)
    else:
        decision = "Normal"
        confidence = min(0.95, 1 - (alpha - p_value)) if p_value > alpha else 0.5

    return decision, confidence


def interpret_effect_size(cohens_d: float) -> str:
    """
    Interpret Cohen's d effect size
    """
    if cohens_d < 0.2:
        return "negligible"
    elif cohens_d < 0.5:
        return "small"
    elif cohens_d < 0.8:
        return "medium"
    else:
        return "large"


def create_detection_report(
    results: Dict[str, Dict], output_file: str, sample_name: str
):
    """
    Create a detailed statistical detection report
    """

    with open(output_file, "w") as f:
        f.write("Statistical Trisomy Detection Report\n")
        f.write(f"Sample: {sample_name}\n")
        f.write(f"Analysis Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("=" * 60 + "\n\n")

        # Summary section for key chromosomes
        key_chroms = ["chr13", "chr18", "chr21"]
        f.write("KEY CHROMOSOMES SUMMARY\n")
        f.write("-" * 30 + "\n")
        for chrom in key_chroms:
            if chrom in results:
                result = results[chrom]
                f.write(
                    f"{chrom}: {result['decision']} (confidence: {result['confidence']:.1%})\n"
                )
        f.write("\n")

        # Detailed analysis for all chromosomes
        f.write("DETAILED ANALYSIS (ALL CHROMOSOMES)\n")
        f.write("-" * 40 + "\n")
        for chrom, result in results.items():
            f.write(f"{chrom.upper()} ANALYSIS\n")
            f.write("-" * 20 + "\n")
            f.write(f"Decision: {result['decision']}\n")
            f.write(f"Confidence: {result['confidence']:.1%}\n")
            f.write(f"P-value: {result['p_value']:.6f}\n")
            f.write(
                f"Effect Size (Cohen's d): {result['effect_size']:.3f} ({result['statistical_summary']['cohens_d_interpretation']})\n"
            )
            f.write(
                f"Mean Z-score: {result['mean_zscore']:.3f} ± {result['std_zscore']:.3f}\n"
            )
            f.write(
                f"Baseline Mean: {result['baseline_mean']:.3f} ± {result['baseline_std']:.3f}\n"
            )
            f.write(
                f"High Z-scores: {result['high_count']}/{result['total_comparisons']} ({result['high_score_ratio']:.1%})\n"
            )
            f.write(
                f"Outliers: {result['outlier_count']}/{result['total_comparisons']}\n"
            )
            f.write(f"Adaptive Threshold: {result['adaptive_threshold']:.3f}\n")
            f.write(f"Top Chromosomes: {', '.join(result['top_chromosomes'])}\n")
            f.write("\n")

        f.write("INTERPRETATION GUIDE\n")
        f.write("-" * 20 + "\n")
        f.write("Detected:  High confidence trisomy detection\n")
        f.write("Suspected: Moderate confidence, requires clinical correlation\n")
        f.write("Possible:  Low confidence, monitor or retest\n")
        f.write("Normal:    No evidence of trisomy\n")
        f.write("\n")
        f.write("Effect Size Interpretation:\n")
        f.write("- Negligible: < 0.2\n")
        f.write("- Small: 0.2 - 0.5\n")
        f.write("- Medium: 0.5 - 0.8\n")
        f.write("- Large: > 0.8\n")

    logger.info(f"Statistical detection report saved to: {output_file}")


def create_detection_table(
    results: Dict[str, Dict], output_file: str, sample_name: str
):
    """
    Create a comprehensive table with all chromosomes and statistical metrics
    """

    # Define chromosome order
    chr_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

    # Prepare data for table
    table_data = []

    for chrom in chr_order:
        if chrom in results:
            result = results[chrom]
            row = {
                "Chromosome": chrom,
                "Decision": result["decision"],
                "Confidence": f"{result['confidence']:.3f}",
                "P_value": f"{result['p_value']:.6f}",
                "Effect_Size_Cohens_d": f"{result['effect_size']:.3f}",
                "Effect_Size_Interpretation": result["statistical_summary"][
                    "cohens_d_interpretation"
                ],
                "Mean_Z_score": f"{result['mean_zscore']:.3f}",
                "Std_Z_score": f"{result['std_zscore']:.3f}",
                "Baseline_Mean": f"{result['baseline_mean']:.3f}",
                "Baseline_Std": f"{result['baseline_std']:.3f}",
                "High_Count": result["high_count"],
                "Total_Comparisons": result["total_comparisons"],
                "High_Score_Ratio": f"{result['high_score_ratio']:.3f}",
                "Outlier_Count": result["outlier_count"],
                "Adaptive_Threshold": f"{result['adaptive_threshold']:.3f}",
                "Outlier_Threshold": f"{result['outlier_threshold']:.3f}",
                "T_Statistic": f"{result['statistical_summary']['t_statistic']:.3f}",
                "Top_5_Chromosomes": "|".join(result["top_chromosomes"]),
            }
        else:
            # If chromosome not analyzed, fill with NA
            row = {
                "Chromosome": chrom,
                "Decision": "NA",
                "Confidence": "NA",
                "P_value": "NA",
                "Effect_Size_Cohens_d": "NA",
                "Effect_Size_Interpretation": "NA",
                "Mean_Z_score": "NA",
                "Std_Z_score": "NA",
                "Baseline_Mean": "NA",
                "Baseline_Std": "NA",
                "High_Count": "NA",
                "Total_Comparisons": "NA",
                "High_Score_Ratio": "NA",
                "Outlier_Count": "NA",
                "Adaptive_Threshold": "NA",
                "Outlier_Threshold": "NA",
                "T_Statistic": "NA",
                "Top_5_Chromosomes": "NA",
            }

        table_data.append(row)

    # Create DataFrame and save as TSV
    df = pd.DataFrame(table_data)

    # Add metadata header
    with open(output_file, "w") as f:
        f.write("# Statistical Trisomy Detection Results\n")
        f.write(f"# Sample: {sample_name}\n")
        f.write(f"# Analysis Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("# \n")
        f.write("# Column Descriptions:\n")
        f.write("# Decision: Detected/Suspected/Possible/Normal\n")
        f.write("# Confidence: Statistical confidence (0-1)\n")
        f.write("# P_value: Two-sample t-test p-value\n")
        f.write("# Effect_Size_Cohens_d: Cohen's d effect size\n")
        f.write("# High_Score_Ratio: Proportion of Z-scores above adaptive threshold\n")
        f.write(
            "# Top_5_Chromosomes: Chromosomes with highest Z-scores (separated by |)\n"
        )
        f.write("# \n")

    # Append the DataFrame
    df.to_csv(output_file, sep="\t", index=False, mode="a")

    logger.info(f"Statistical detection table saved to: {output_file}")

    return df


def run_statistical_trisomy_detection(
    zscore_chr: pd.DataFrame, output_dir: str, sample_name: str
) -> Dict[str, Dict]:
    """
    Run the complete statistical trisomy detection pipeline
    """

    logger.info("=== Running Statistical Trisomy Detection ===")

    # Run detection for all chromosomes
    results = statistical_trisomy_detection(
        zscore_chr=zscore_chr,
        target_chromosomes=None,  # Analyze all chromosomes
        confidence_level=0.95,
        effect_size_threshold=0.8,
        min_high_scores=15,
    )

    # Create detailed report
    report_file = f"{output_dir}/{sample_name}.trisomy_detection.txt"
    create_detection_report(results, report_file, sample_name)

    # Create comprehensive table
    table_file = f"{output_dir}/{sample_name}.trisomy_detection.tsv"
    detection_table = create_detection_table(results, table_file, sample_name)

    # Create summary for main pipeline (focus on key chromosomes)
    summary = {}
    key_chromosomes = ["chr13", "chr18", "chr21"]

    for chrom in key_chromosomes:
        if chrom in results:
            summary[chrom] = results[chrom]["decision"]
        else:
            summary[chrom] = "NA"

    logger.info("Statistical Trisomy Detection Summary (Key Chromosomes):")
    for chrom, decision in summary.items():
        if chrom in results:
            confidence = results[chrom]["confidence"]
            logger.info(f"  {chrom}: {decision} (confidence: {confidence:.1%})")
        else:
            logger.info(f"  {chrom}: {decision}")

    return results


# ============================================================================
# FILE SAVE FUNCTIONS
# ============================================================================


def save_zscore_results(zscore_df: pd.DataFrame, output_file: str):
    """Save Z-score results to file"""
    # Ensure directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    zscore_df.to_csv(output_file, sep="\t")
    logger.info(f"Z-score saved to: {output_file}")


def perform_qc_analysis_old(
    zscore_df: pd.DataFrame, cutoff: float, output_file: str
) -> bool:
    """Perform quality control analysis"""

    # Flatten Z-score matrix
    zscore_values = zscore_df.values.flatten()
    zscore_sorted = sorted(zscore_values)

    # Calculate 75th percentile
    qc_point = int(len(zscore_sorted) * 0.75)
    qc_value = zscore_sorted[qc_point - 1]

    passed = qc_value <= cutoff
    status = "PASS" if passed else "FAIL"

    # Ensure directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, "w") as f:
        f.write(f"{qc_value:.6f}\t{status}\n")

    logger.info(f"QC result: {status} (value: {qc_value:.6f})")
    return passed


# ============================================================================
# PLOTTING FUNCTIONS (UPDATED)
# ============================================================================


def plot_zscore_by_chromosome_filtered(
    df, output_dir, sample_name="sample", extremes_count=2, file_prefix=None, dpi=200
):
    # Data copy
    df_mod = df.copy()

    # Exclude '22' and '23' columns
    cols_to_use = [col for col in df_mod.columns if col not in ["22", "23"]]
    df_filtered = df_mod.loc[:, cols_to_use]

    # Chromosome name mapping
    chr_mapping = {}
    for idx in df_mod.index:
        if idx == "22":
            chr_mapping[idx] = "chrX"
        elif idx == "23":
            chr_mapping[idx] = "chrY"
        elif "chr" in str(idx).lower():
            chr_label = str(idx).lower()
            if chr_label == "chrx":
                chr_label = "chrX"
            elif chr_label == "chry":
                chr_label = "chrY"
            chr_mapping[idx] = chr_label
        else:
            if idx in ["X", "Y", "x", "y"]:
                chr_mapping[idx] = f"chr{idx.upper()}"
            else:
                chr_mapping[idx] = f"chr{idx}"

    # Analyze extreme values for each chromosome
    high_value_counts = Counter()
    low_value_counts = Counter()

    for i, row_idx in enumerate(df_filtered.index):
        row_data = df_filtered.loc[row_idx].values

        # Find top extreme values
        max_indices = np.argsort(row_data)[-extremes_count:]
        for idx in max_indices:
            col_name = df_filtered.columns[idx]
            if col_name in chr_mapping:
                mapped_name = chr_mapping[col_name]
            else:
                if "chr" not in str(col_name).lower():
                    if col_name in ["X", "Y", "x", "y"]:
                        mapped_name = f"chr{col_name.upper()}"
                    else:
                        mapped_name = f"chr{col_name}"
                else:
                    mapped_name = col_name
            high_value_counts[mapped_name] += 1

        # Find bottom extreme values
        min_indices = np.argsort(row_data)[:extremes_count]
        for idx in min_indices:
            col_name = df_filtered.columns[idx]
            if col_name in chr_mapping:
                mapped_name = chr_mapping[col_name]
            else:
                if "chr" not in str(col_name).lower():
                    if col_name in ["X", "Y", "x", "y"]:
                        mapped_name = f"chr{col_name.upper()}"
                    else:
                        mapped_name = f"chr{col_name}"
                else:
                    mapped_name = col_name
            low_value_counts[mapped_name] += 1

    # Filtering criteria
    total_chromosomes = len(df_filtered.index)
    threshold = total_chromosomes - 1

    logger.info(
        f"Most frequent high Z-score chromosomes: {high_value_counts.most_common(3)}"
    )
    logger.info(
        f"Most frequent low Z-score chromosomes: {low_value_counts.most_common(3)}"
    )

    # Identify chromosomes to filter
    chromosomes_to_filter = []
    for chrom, count in high_value_counts.items():
        if count >= threshold:
            chromosomes_to_filter.append(chrom)
            logger.info(
                f"Filtering: {chrom} - high values found in {count} chromosomes"
            )

    for chrom, count in low_value_counts.items():
        if count >= threshold:
            chromosomes_to_filter.append(chrom)
            logger.info(f"Filtering: {chrom} - low values found in {count} chromosomes")

    # Apply filtering (set extreme values to 0)
    filtered_data = df_filtered.copy()
    for col in filtered_data.columns:
        if col in chr_mapping:
            mapped_col = chr_mapping[col]
        elif "chr" not in str(col).lower():
            if col in ["X", "Y", "x", "y"]:
                mapped_col = f"chr{col.upper()}"
            else:
                mapped_col = f"chr{col}"
        else:
            mapped_col = col

        if mapped_col in chromosomes_to_filter:
            filtered_data[col] = 0.0

    # Create plot
    plt.figure(figsize=(20, 6))

    z_threshold = 3
    all_z_scores = []
    chr_boundaries = [0]
    chr_labels = []

    # Concatenate all chromosome Z-score values
    for i, row_idx in enumerate(filtered_data.index):
        z_scores = filtered_data.loc[row_idx].values

        clean_z_scores = []
        for z in z_scores:
            try:
                clean_z_scores.append(float(z))
            except (ValueError, TypeError):
                clean_z_scores.append(0.0)

        all_z_scores.extend(clean_z_scores)
        chr_boundaries.append(chr_boundaries[-1] + len(clean_z_scores))
        chr_labels.append((chr_boundaries[-2] + chr_boundaries[-1]) / 2)

    x_pos = np.arange(len(all_z_scores))
    ax = plt.subplot(111)

    # Plot Z-score line
    ax.plot(
        x_pos, all_z_scores, color="blue", linewidth=1, label="PRIZM Z-score (filtered)"
    )

    # Threshold lines
    ax.axhline(
        y=z_threshold,
        color="gray",
        linestyle="--",
        linewidth=0.8,
        alpha=0.7,
        label="Threshold",
    )
    ax.axhline(y=-z_threshold, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
    ax.axhline(y=0, color="gray", linestyle="-", linewidth=0.5, alpha=0.5)

    # Grid
    ax.grid(True, axis="y", linestyle=":", alpha=0.3)
    ax.grid(True, axis="x", linestyle="-", alpha=0.2)

    # Highlight threshold violations
    above_threshold_x = []
    above_threshold_y = []
    below_threshold_x = []
    below_threshold_y = []

    for i, z in enumerate(all_z_scores):
        if z > z_threshold:
            above_threshold_x.append(i)
            above_threshold_y.append(z)
        elif z < -z_threshold:
            below_threshold_x.append(i)
            below_threshold_y.append(z)

    if above_threshold_x:
        ax.scatter(
            above_threshold_x,
            above_threshold_y,
            color="red",
            s=15,
            marker="o",
            label="Values above threshold",
        )

    if below_threshold_x:
        ax.scatter(
            below_threshold_x,
            below_threshold_y,
            color="darkred",
            s=15,
            marker="s",
            label="Values below threshold",
        )

    # Chromosome boundaries
    for boundary in chr_boundaries[1:-1]:
        ax.axvline(x=boundary, color="black", linestyle="-", linewidth=0.5, alpha=0.5)

    # Chromosome labels
    row_indices = list(filtered_data.index)
    for i, pos in enumerate(chr_labels):
        if i < len(row_indices):
            row_idx = row_indices[i]
            mapped_name = chr_mapping.get(row_idx, row_idx)
            ax.text(pos, -5.8, mapped_name, ha="center", fontsize=10, fontweight="bold")

            if mapped_name in chromosomes_to_filter:
                ax.text(
                    pos,
                    -5.2,
                    "(filtered)",
                    ha="center",
                    fontsize=8,
                    color="red",
                    style="italic",
                )

    # Highlight special chromosomes
    special_chrs = ["chr13", "chr18", "chr21", "chrX", "chrY"]
    for i, row_idx in enumerate(row_indices):
        mapped_name = chr_mapping.get(row_idx, row_idx)
        if mapped_name in special_chrs and i < len(chr_boundaries) - 1:
            start = chr_boundaries[i]
            end = chr_boundaries[i + 1]
            highlight = patches.Rectangle(
                (start, -6), end - start, 12, facecolor="pink", alpha=0.2, zorder=0
            )
            ax.add_patch(highlight)

    # Axis settings
    ax.set_xlim(0, len(all_z_scores))
    ax.set_ylim(-6, 6)
    ax.set_ylabel("PRIZM Z-score", fontsize=12)

    filtered_info = ""
    if chromosomes_to_filter:
        filtered_info = f" (Zeroed: {', '.join(chromosomes_to_filter)})"

    ax.set_title(
        f"PRIZM Z-score versus chromosomal position - Sample {sample_name}{filtered_info}",
        fontsize=14,
    )
    ax.set_xticks([])
    ax.set_yticks([-3, 0, 3])
    ax.legend(loc="upper right", frameon=True, fontsize=9)

    plt.tight_layout()

    # Generate filename
    if file_prefix:
        filename = f"{output_dir}/{file_prefix}_chromosome_line.png"
    else:
        filename = f"{output_dir}/prizm_chromosome_zscore_plot.png"

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(filename, dpi=dpi, bbox_inches="tight")
    plt.close()
    logger.info(f"PRIZM chromosome line plot saved to: {filename}")


def plot_10mb_zscore_by_chromosome_filtered(
    df, output_dir, sample_name="10mb", extremes_count=2, file_prefix=None, dpi=200
):
    df_mod = df.copy()

    logger.info(f"Matrix: {df_mod.shape}")
    logger.info(f"Row count: {len(df_mod.index)}")

    # 염색체 경계 위치
    x_ticks_loc = [
        0,
        25,
        50,
        70,
        90,
        109,
        127,
        143,
        158,
        173,
        187,
        201,
        215,
        227,
        238,
        249,
        259,
        268,
        276,
        282,
        289,
        294,
        300,
        316,
        322,
    ]

    # 염색체별 bin 수 계산
    chr_bin_counts = []
    for i in range(len(x_ticks_loc) - 1):
        chr_bin_counts.append(x_ticks_loc[i + 1] - x_ticks_loc[i])

    chr_labels = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

    # -3 이하 값 추적을 위한 딕셔너리
    below_threshold_chroms = defaultdict(Counter)

    # 각 bin 위치(컬럼)에서 -3 이하 값 분석
    for col_idx in range(len(df_mod.columns)):
        for row_idx, row_name in enumerate(df_mod.index):
            if col_idx < len(df_mod.columns):
                try:
                    val = float(df_mod.loc[row_name, df_mod.columns[col_idx]])
                    if val < -3:
                        # 행 이름(chr1, chr2 등)을 키로 사용하고, 행 인덱스를 값으로 카운트
                        below_threshold_chroms[row_name][row_idx] += 1
                except (ValueError, TypeError):
                    continue

    # -3 이하 값 분석 결과 출력
    # logger.info("\n===== -3 이하 Z-score 값 분석 =====")
    # for target_chr, idx_counter in below_threshold_chroms.items():
    #    if idx_counter:
    #        logger.info(f"{target_chr}: {dict(idx_counter)}")

    # 가장 자주 -3 이하 값으로 나타나는 인덱스 찾기
    all_neg_indices = Counter()
    for target_chr, idx_counter in below_threshold_chroms.items():
        for idx, count in idx_counter.items():
            all_neg_indices[idx] += count

    # logger.info("\nMost frequently shown below -3 valued chromosome indices:", all_neg_indices.most_common(3))

    # 필터링 대상 인덱스 선정 (전체 염색체 중 절반 이상에서 -3 이하 값을 가진 인덱스)
    threshold_count = len(df_mod.index) // 2
    indices_to_filter = [
        idx for idx, count in all_neg_indices.items() if count >= threshold_count
    ]
    # logger.info(f"\n필터링할 염색체 인덱스: {indices_to_filter}")

    # 필터링 적용 (극단값을 0으로 치환)
    filtered_data = df_mod.copy()

    # 필터링 대상 인덱스에 해당하는 행의 모든 값을 0으로 치환
    rows_to_filter = []
    for idx in indices_to_filter:
        if idx < len(filtered_data.index):
            row_name = filtered_data.index[idx]
            rows_to_filter.append(row_name)
            logger.info(
                f"Filter : {row_name} (Index {idx}) - {all_neg_indices[idx]} Found below -3"
            )
            filtered_data.loc[row_name, :] = 0.0

    # 그림 설정
    plt.figure(figsize=(20, 6))

    # Z-score 임계값
    z_threshold = 3

    # 모든 데이터를 염색체별로 재구성
    all_z_scores = []
    chr_boundaries = [0]  # 각 염색체의 경계 위치
    chr_labels_pos = []  # 각 염색체의 레이블 위치

    # 각 염색체에 대해
    for chr_idx in range(len(chr_bin_counts)):
        if chr_idx < len(chr_labels):
            target_chr = chr_labels[chr_idx]
        else:
            target_chr = f"chr{chr_idx + 1}"

        chr_start = len(all_z_scores)

        # 현재 염색체의 bin 수
        bin_count = chr_bin_counts[chr_idx]
        col_offset = x_ticks_loc[chr_idx]

        # 각 bin에 대해
        for bin_idx in range(bin_count):
            # 컬럼 인덱스 계산
            col_idx = col_offset + bin_idx

            # 모든 행(염색체)에 대해 해당 bin의 Z-score 값 추가
            for row_idx in filtered_data.index:
                if col_idx < len(filtered_data.columns):
                    try:
                        val = float(
                            filtered_data.loc[row_idx, filtered_data.columns[col_idx]]
                        )
                        all_z_scores.append(val)
                    except (ValueError, TypeError):
                        all_z_scores.append(
                            0.0
                        )  # 숫자로 변환할 수 없는 값은 0으로 처리

        # 염색체 경계 위치 추가
        chr_boundaries.append(len(all_z_scores))

        # 염색체 레이블 위치 (중간 지점)
        chr_labels_pos.append((chr_boundaries[-2] + chr_boundaries[-1]) / 2)

    # x축 생성 (각 값의 위치)
    x_pos = np.arange(len(all_z_scores))

    # 메인 플롯 생성
    ax = plt.subplot(111)

    # Z-score 선 그리기
    ax.plot(x_pos, all_z_scores, color="blue", linewidth=1, label="Z-score")

    # Z-score 임계값 선 그리기
    ax.axhline(
        y=z_threshold,
        color="gray",
        linestyle="--",
        linewidth=0.8,
        alpha=0.7,
        label="Threshold",
    )
    ax.axhline(y=-z_threshold, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
    ax.axhline(y=0, color="gray", linestyle="-", linewidth=0.5, alpha=0.5)

    # 격자 추가
    ax.grid(True, axis="y", linestyle=":", alpha=0.3)

    # 임계값을 넘는 영역 강조
    above_threshold = [(i, z) for i, z in enumerate(all_z_scores) if z > z_threshold]
    below_threshold = [(i, z) for i, z in enumerate(all_z_scores) if z < -z_threshold]

    if above_threshold:
        above_x, above_y = zip(*above_threshold)
        ax.scatter(
            above_x,
            above_y,
            color="red",
            s=15,
            marker="o",
            label="Values above threshold",
        )

    if below_threshold:
        below_x, below_y = zip(*below_threshold)
        ax.scatter(
            below_x,
            below_y,
            color="darkred",
            s=15,
            marker="s",
            label="Values below threshold",
        )

    # 염색체 경계 표시 (수직선)
    for boundary in chr_boundaries[1:-1]:  # 첫 번째와 마지막은 제외
        ax.axvline(x=boundary, color="black", linestyle="-", linewidth=0.5, alpha=0.5)

    # 염색체 레이블 추가
    for i, pos in enumerate(chr_labels_pos):
        if i < len(chr_labels):
            # 염색체 이름 표시
            ax.text(
                pos, -5.8, chr_labels[i], ha="center", fontsize=8, fontweight="bold"
            )

    # 특정 염색체 하이라이트 (예: 13, 18, 21번 염색체)
    special_indices = [12, 17, 20]  # chr13, chr18, chr21의 인덱스 (0-기반)
    for idx in special_indices:
        if idx < len(chr_boundaries) - 1:
            start = chr_boundaries[idx]
            end = chr_boundaries[idx + 1]
            highlight = patches.Rectangle(
                (start, -6), end - start, 12, facecolor="pink", alpha=0.2, zorder=0
            )
            ax.add_patch(highlight)

    # 축 설정
    ax.set_xlim(0, len(all_z_scores))
    ax.set_ylim(-6, 6)
    ax.set_ylabel("Z-score", fontsize=12)

    # 필터링 정보 제목에 추가
    # filtered_info = ""
    # if rows_to_filter:
    #    filtered_info = f" (Filtered rows: {', '.join(str(r) for r in rows_to_filter)})"

    ax.set_title(f"Z-score by chromosome (10mb bins) - {sample_name}", fontsize=14)

    # x축 눈금 제거 (염색체 이름으로 대체)
    ax.set_xticks([])

    # y축 눈금
    ax.set_yticks([-3, 0, 3])

    # 범례 추가
    ax.legend(loc="upper right", frameon=True, fontsize=9)

    # 그림 저장
    plt.tight_layout()

    # Generate filename
    if file_prefix:
        filename = f"{output_dir}/{file_prefix}_10mb_line.png"
    else:
        filename = f"{output_dir}/prizm_chromosome_zscore_plot.png"

    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(filename, dpi=dpi, bbox_inches="tight")
    plt.close()
    logger.info(f"PRIZM 10mb chromosome line plot saved to: {filename}")

    return


def generate_prizm_plots(
    results: PRIZMResult, output_dir: str, sample_name: str, group: str = "", dpi=200
):
    """Generate PRIZM plots - both heatmaps and line plots"""

    logger.info("Generating PRIZM visualization plots...")

    try:
        os.makedirs(output_dir, exist_ok=True)

        # Generate file name prefix
        if group:
            file_prefix = f"{sample_name}_{group}"
        else:
            file_prefix = sample_name

        # ===== HEATMAP PLOTS =====

        # 1. Chromosome heatmap (24x24 matrix)
        plt.figure(figsize=(12, 10))
        plt.imshow(results.zscore_chr.values, cmap="RdBu_r", vmin=-7, vmax=7)
        plt.colorbar(label="Z-score")
        plt.title(f"PRIZM Chromosome Z-scores Heatmap - {sample_name}")
        plt.xticks(
            range(len(results.zscore_chr.columns)),
            results.zscore_chr.columns,
            rotation=45,
        )
        plt.yticks(range(len(results.zscore_chr.index)), results.zscore_chr.index)
        plt.tight_layout()
        plt.savefig(
            f"{output_dir}/{file_prefix}_chromosome_heatmap.png",
            dpi=dpi,
            bbox_inches="tight",
        )
        plt.close()

        # 2. 10mb heatmap (using 10mb_all if available, otherwise regular 10mb)
        if results.zscore_10mb_all is not None:
            # 10mb_all heatmap (22x322 matrix)
            plt.figure(figsize=(50, 5))  # Wide format for 322 columns
            plt.imshow(results.zscore_10mb_all.values, cmap="RdBu_r", vmin=-7, vmax=7)
            plt.colorbar(label="Z-score")
            plt.title(f"PRIZM 10mb_all Z-scores Heatmap - {sample_name}")

            # Set chromosome labels on Y-axis (22 chromosomes)
            chr_labels_y = [f"chr{i}" for i in range(1, 23)]
            plt.yticks(range(22), chr_labels_y)

            # Set 10mb bin positions on X-axis with correct positions
            x_ticks_loc = [
                0,
                25,
                50,
                70,
                90,
                109,
                127,
                143,
                158,
                173,
                187,
                201,
                215,
                227,
                238,
                249,
                259,
                268,
                276,
                282,
                289,
                294,
                300,
                316,
            ]
            x_ticks_name = [
                "chr1",
                "chr2",
                "chr3",
                "chr4",
                "chr5",
                "chr6",
                "chr7",
                "chr8",
                "chr9",
                "chr10",
                "chr11",
                "chr12",
                "chr13",
                "chr14",
                "chr15",
                "chr16",
                "chr17",
                "chr18",
                "chr19",
                "chr20",
                "chr21",
                "chr22",
                "chrX",
                "chrY",
            ]

            # Only use positions that are within the data range
            valid_ticks_loc = [
                pos for pos in x_ticks_loc if pos < results.zscore_10mb_all.shape[1]
            ]
            valid_ticks_name = x_ticks_name[: len(valid_ticks_loc)]

            plt.xticks(valid_ticks_loc, valid_ticks_name, rotation=60)

            # Add vertical lines at chromosome boundaries
            vline_loc = [
                x - 0.5 for x in valid_ticks_loc if x > 0
            ]  # Exclude the first position
            for vline in vline_loc:
                plt.axvline(
                    x=vline, color="black", linestyle="-", linewidth=1.0, alpha=0.7
                )

            # Minor ticks for better visualization
            plt.xticks(np.arange(-0.5, results.zscore_10mb_all.shape[1], 1), minor=True)
            plt.yticks(np.arange(-0.5, 22, 1), minor=True)
            plt.grid(which="minor", linestyle=":", linewidth=0.3)

            plt.xlabel("10mb bins (by chromosome)")
            plt.ylabel("Chromosomes")
            plt.tight_layout()
            plt.savefig(
                f"{output_dir}/{file_prefix}_10mb_heatmap.png",
                dpi=dpi,
                bbox_inches="tight",
            )
            plt.close()

            logger.info(
                f"10mb_all heatmap saved to: {output_dir}/{file_prefix}_10mb_heatmap.png"
            )
        else:
            # Regular 10mb heatmap (fallback)
            plt.figure(figsize=(15, 12))
            plt.imshow(results.zscore_10mb.values, cmap="RdBu_r", vmin=-7, vmax=7)
            plt.colorbar(label="Z-score")
            plt.title(f"PRIZM 10mb Z-scores Heatmap - {sample_name}")
            plt.xticks(
                range(len(results.zscore_10mb.columns)),
                results.zscore_10mb.columns,
                rotation=45,
            )
            plt.yticks(range(len(results.zscore_10mb.index)), results.zscore_10mb.index)
            plt.tight_layout()
            plt.savefig(
                f"{output_dir}/{file_prefix}_10mb_heatmap.png",
                dpi=dpi,
                bbox_inches="tight",
            )
            plt.close()

        # ===== LINE PLOTS =====

        # 3. Chromosome line plot (filtered extreme values)
        plot_zscore_by_chromosome_filtered(
            results.zscore_chr,
            output_dir,
            sample_name=sample_name,
            file_prefix=file_prefix,
        )

        # 4. 10mb line plot (filtered)
        if results.zscore_10mb_all is not None:
            # Use 10mb_all for line plot
            plot_10mb_zscore_by_chromosome_filtered(
                results.zscore_10mb_all,
                output_dir,
                sample_name=sample_name,
                file_prefix=file_prefix,
            )
        else:
            # Fallback to regular 10mb
            plot_10mb_zscore_by_chromosome_filtered(
                results.zscore_10mb,
                output_dir,
                sample_name=sample_name,
                file_prefix=file_prefix,
            )

        logger.info(f"All PRIZM plots generated successfully in: {output_dir}")

    except Exception as e:
        logger.error(f"Error generating PRIZM plots: {e}")
        raise


# ============================================================================
# PIPELINE INTEGRATION FUNCTIONS
# ============================================================================


def calculate_extreme_ratio_qc(zscore_matrix, threshold=3.0, trisomy_candidates=None):
    """
    절댓값 3 이상인 값들의 비율로 QC 판단
    trisomy 후보 염색체는 제외하고 계산
    """
    values = zscore_matrix.flatten()

    # Trisomy 후보가 있으면 해당 염색체 제외
    if trisomy_candidates:
        filtered_values = []
        n_chroms = zscore_matrix.shape[0]

        for i in range(len(values)):
            row_idx = i // n_chroms
            col_idx = i % n_chroms

            # 해당 염색체의 row나 column이면 제외
            if row_idx not in trisomy_candidates and col_idx not in trisomy_candidates:
                filtered_values.append(values[i])

        analysis_values = np.array(filtered_values) if filtered_values else values
    else:
        analysis_values = values

    extreme_count = np.sum(np.abs(analysis_values) > threshold)
    extreme_ratio = extreme_count / len(analysis_values)

    if extreme_ratio <= 0.05:
        status = "PASS"
    elif extreme_ratio <= 0.10:
        status = "WARNING"
    else:
        status = "FAIL"

    return {
        "extreme_ratio": extreme_ratio,
        "extreme_count": extreme_count,
        "total_count": len(analysis_values),
        "excluded_trisomy": trisomy_candidates is not None,
        "status": status,
    }


def calculate_variance_qc(zscore_matrix, trisomy_candidates=None):
    """
    Z-score 분포의 표준편차로 QC 판단
    trisomy 후보 염색체는 제외하고 계산
    """
    values = zscore_matrix.flatten()

    # Trisomy 후보가 있으면 해당 염색체 제외
    if trisomy_candidates:
        filtered_values = []
        n_chroms = zscore_matrix.shape[0]

        for i in range(len(values)):
            row_idx = i // n_chroms
            col_idx = i % n_chroms

            # 해당 염색체의 row나 column이면 제외
            if row_idx not in trisomy_candidates and col_idx not in trisomy_candidates:
                filtered_values.append(values[i])

        analysis_values = np.array(filtered_values) if filtered_values else values
    else:
        analysis_values = values

    std_dev = np.std(analysis_values)
    mean_val = np.mean(analysis_values)

    # 정규분포에서 벗어난 정도
    deviation_score = abs(std_dev - 1.0) + abs(mean_val)

    if deviation_score <= 1.0:
        status = "PASS"
    elif deviation_score <= 2.0:
        status = "WARNING"
    else:
        status = "FAIL"

    return {
        "std_deviation": std_dev,
        "mean_value": mean_val,
        "deviation_score": deviation_score,
        "excluded_trisomy": trisomy_candidates is not None,
        "status": status,
    }


def detect_potential_trisomy_for_qc(zscore_matrix, chromosome_names, threshold=2.5):
    """
    QC를 위한 trisomy 후보 염색체 검출
    """
    trisomy_candidates = []

    for i, chrom in enumerate(chromosome_names):
        # 해당 염색체 row의 평균 Z-score (자기 자신 제외)
        row_values = zscore_matrix[i, :]
        row_mean = np.mean([val for j, val in enumerate(row_values) if j != i])

        # 해당 염색체 column의 평균 Z-score (자기 자신 제외)
        col_values = zscore_matrix[:, i]
        col_mean = np.mean([val for j, val in enumerate(col_values) if j != i])

        # 양방향 모두 높으면 trisomy 후보
        if row_mean > threshold and col_mean > threshold:
            trisomy_candidates.append(i)
            logger.info(
                f"QC: Potential trisomy detected for {chrom} (row_mean: {row_mean:.3f}, col_mean: {col_mean:.3f})"
            )

    return trisomy_candidates


def calculate_diagonal_qc(zscore_matrix, trisomy_candidates=None):
    """
    대각선 값들(self-comparison)의 안정성으로 QC 판단
    trisomy 후보는 대각선에서 제외
    """
    diagonal_values = np.diag(zscore_matrix)

    # Trisomy 후보가 있으면 해당 대각선 값 제외
    if trisomy_candidates:
        filtered_diagonal = [
            val for i, val in enumerate(diagonal_values) if i not in trisomy_candidates
        ]
        analysis_diagonal = (
            np.array(filtered_diagonal) if filtered_diagonal else diagonal_values
        )
    else:
        analysis_diagonal = diagonal_values

    diagonal_std = np.std(analysis_diagonal)
    diagonal_mean = np.mean(np.abs(analysis_diagonal))

    # 대각선 값들이 얼마나 0에 가까운지
    if diagonal_std <= 0.3 and diagonal_mean <= 0.2:
        status = "PASS"
    elif diagonal_std <= 0.6 and diagonal_mean <= 0.5:
        status = "WARNING"
    else:
        status = "FAIL"

    return {
        "diagonal_std": diagonal_std,
        "diagonal_mean_abs": diagonal_mean,
        "diagonal_values": diagonal_values.tolist(),
        "excluded_trisomy": trisomy_candidates is not None,
        "excluded_count": len(trisomy_candidates) if trisomy_candidates else 0,
        "status": status,
    }


def perform_qc_analysis(
    zscore_df: pd.DataFrame, cutoff: float, output_file: str
) -> bool:
    """
    Enhanced QC analysis with trisomy-aware evaluation
    """
    zscore_matrix = zscore_df.values
    chromosome_names = list(zscore_df.index)

    # 0. Detect potential trisomy for QC adjustment
    trisomy_candidates = detect_potential_trisomy_for_qc(
        zscore_matrix, chromosome_names
    )

    # 1. Legacy QC (75th percentile) - 호환성을 위해 유지
    zscore_values = zscore_matrix.flatten()
    zscore_sorted = sorted(zscore_values)
    qc_point = int(len(zscore_sorted) * 0.75)
    legacy_qc_value = zscore_sorted[qc_point - 1]
    legacy_passed = legacy_qc_value <= cutoff

    # 2. Enhanced QC Methods (trisomy-aware)
    extreme_qc = calculate_extreme_ratio_qc(
        zscore_matrix, trisomy_candidates=trisomy_candidates
    )
    variance_qc = calculate_variance_qc(
        zscore_matrix, trisomy_candidates=trisomy_candidates
    )
    diagonal_qc = calculate_diagonal_qc(
        zscore_matrix, trisomy_candidates=trisomy_candidates
    )

    # 3. Additional metrics
    additional_metrics = {
        "max_absolute_zscore": np.max(np.abs(zscore_values)),
        "values_above_5": np.sum(np.abs(zscore_values) > 5),
        "negative_extreme_ratio": np.sum(zscore_values < -3) / len(zscore_values),
        "positive_extreme_ratio": np.sum(zscore_values > 3) / len(zscore_values),
        "potential_trisomy_count": len(trisomy_candidates),
    }

    # 4. Scoring system
    scores = []
    if extreme_qc["status"] == "PASS":
        scores.append(3)
    elif extreme_qc["status"] == "WARNING":
        scores.append(2)
    else:
        scores.append(0)

    if variance_qc["status"] == "PASS":
        scores.append(3)
    elif variance_qc["status"] == "WARNING":
        scores.append(2)
    else:
        scores.append(0)

    if diagonal_qc["status"] == "PASS":
        scores.append(3)
    elif diagonal_qc["status"] == "WARNING":
        scores.append(2)
    else:
        scores.append(0)

    total_score = sum(scores)

    # Trisomy가 검출된 경우 점수 조정 (더 관대한 기준)
    if trisomy_candidates:
        if total_score >= 6:  # Trisomy 있을 때는 더 관대하게
            final_status = "PASS"
            final_passed = True
        elif total_score >= 4:
            final_status = "WARNING"
            final_passed = True
        else:
            final_status = "FAIL"
            final_passed = False
    else:
        # 일반적인 기준
        if total_score >= 8:
            final_status = "PASS"
            final_passed = True
        elif total_score >= 5:
            final_status = "WARNING"
            final_passed = True
        else:
            final_status = "FAIL"
            final_passed = False

    # 5. Write detailed report
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Legacy format file (호환성)
    with open(output_file, "w") as f:
        f.write(f"{legacy_qc_value:.6f}\t{'PASS' if legacy_passed else 'FAIL'}\n")

    # Enhanced QC report
    enhanced_report_file = output_file.replace(".qc.txt", ".enhanced_qc.txt")
    with open(enhanced_report_file, "w") as f:
        f.write("PRIZM Enhanced Quality Control Report (Trisomy-Aware)\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Overall Status: {final_status}\n")
        f.write(f"Quality Score: {total_score}/9\n")
        f.write(
            f"Legacy QC (75th percentile): {legacy_qc_value:.6f} - {'PASS' if legacy_passed else 'FAIL'}\n"
        )

        # Trisomy detection results
        if trisomy_candidates:
            f.write("\nTrisomy Detection for QC:\n")
            f.write(
                f"Potential trisomy chromosomes detected: {len(trisomy_candidates)}\n"
            )
            for idx in trisomy_candidates:
                chrom_name = (
                    chromosome_names[idx]
                    if idx < len(chromosome_names)
                    else f"chr{idx + 1}"
                )
                f.write(f"  - {chrom_name} (index {idx})\n")
            f.write("Note: These chromosomes were excluded from QC calculations\n")
        else:
            f.write("\nTrisomy Detection for QC: No potential trisomy detected\n")

        f.write("\nDetailed Analysis:\n")
        f.write("-" * 30 + "\n")
        f.write(f"1. Extreme Values Check: {extreme_qc['status']}\n")
        f.write(
            f"   - Extreme ratio: {extreme_qc['extreme_ratio']:.4f} ({extreme_qc['extreme_count']}/{extreme_qc['total_count']})\n"
        )
        f.write(
            f"   - Trisomy exclusion: {'Yes' if extreme_qc['excluded_trisomy'] else 'No'}\n"
        )
        f.write("   - Threshold: |Z-score| > 3.0\n")
        f.write("   - Recommendation: Keep below 5%\n\n")

        f.write(f"2. Variance Stability: {variance_qc['status']}\n")
        f.write(f"   - Standard deviation: {variance_qc['std_deviation']:.4f}\n")
        f.write(f"   - Mean deviation: {variance_qc['mean_value']:.4f}\n")
        f.write(f"   - Deviation score: {variance_qc['deviation_score']:.4f}\n")
        f.write(
            f"   - Trisomy exclusion: {'Yes' if variance_qc['excluded_trisomy'] else 'No'}\n"
        )
        f.write("   - Recommendation: Std should be close to 1.0, mean close to 0\n\n")

        f.write(f"3. Diagonal Stability: {diagonal_qc['status']}\n")
        f.write(f"   - Diagonal std: {diagonal_qc['diagonal_std']:.4f}\n")
        f.write(f"   - Diagonal mean(abs): {diagonal_qc['diagonal_mean_abs']:.4f}\n")
        f.write(
            f"   - Trisomy exclusion: {'Yes' if diagonal_qc['excluded_trisomy'] else 'No'}\n"
        )
        if diagonal_qc["excluded_trisomy"]:
            f.write(f"   - Excluded diagonal count: {diagonal_qc['excluded_count']}\n")
        f.write("   - Recommendation: Self-comparisons should be close to 0\n\n")

        f.write("Additional Metrics:\n")
        f.write("-" * 20 + "\n")
        for key, value in additional_metrics.items():
            if isinstance(value, int):
                f.write(f"   {key}: {value}\n")
            else:
                f.write(f"   {key}: {value:.4f}\n")

        f.write("\nInterpretation Guide:\n")
        f.write("-" * 20 + "\n")
        f.write("PASS: High quality, reliable for analysis\n")
        f.write("WARNING: Acceptable but monitor closely\n")
        f.write("FAIL: Quality issues, review sample preparation\n\n")

        f.write("Scoring Details:\n")
        f.write("-" * 15 + "\n")
        f.write("Each criterion scores 0-3 points (FAIL=0, WARNING=2, PASS=3)\n")
        f.write("Total possible: 9 points\n")
        if trisomy_candidates:
            f.write(
                "Trisomy-adjusted thresholds: ≥6 points (PASS), ≥4 points (WARNING), <4 points (FAIL)\n"
            )
        else:
            f.write(
                "Standard thresholds: ≥8 points (PASS), ≥5 points (WARNING), <5 points (FAIL)\n"
            )

        f.write("\nTrisomy-Aware QC Notes:\n")
        f.write("-" * 25 + "\n")
        f.write(
            "- Chromosomes with consistently high Z-scores are automatically detected\n"
        )
        f.write("- These are excluded from QC calculations to avoid false negatives\n")
        f.write(
            "- Scoring thresholds are adjusted when potential trisomy is detected\n"
        )
        f.write("- This ensures quality samples with true trisomy are not rejected\n")

    # Log results
    logger.info("Trisomy-Aware QC Analysis:")
    logger.info(f"  - Overall Status: {final_status} (Score: {total_score}/9)")
    logger.info(
        f"  - Legacy QC: {'PASS' if legacy_passed else 'FAIL'} (75th percentile: {legacy_qc_value:.6f})"
    )
    if trisomy_candidates:
        trisomy_names = [
            chromosome_names[i] if i < len(chromosome_names) else f"chr{i + 1}"
            for i in trisomy_candidates
        ]
        logger.info(
            f"  - Potential Trisomy: {', '.join(trisomy_names)} (excluded from QC)"
        )
    logger.info(
        f"  - Extreme Values: {extreme_qc['status']} (ratio: {extreme_qc['extreme_ratio']:.4f})"
    )
    logger.info(
        f"  - Variance: {variance_qc['status']} (std: {variance_qc['std_deviation']:.4f})"
    )
    logger.info(
        f"  - Diagonal: {diagonal_qc['status']} (std: {diagonal_qc['diagonal_std']:.4f})"
    )
    logger.info(f"Enhanced QC report saved to: {enhanced_report_file}")

    return final_passed


def calculate_zscore_stats(zscore_df):
    """Calculate summary statistics for Z-score matrix"""
    values = zscore_df.values.flatten()

    return {
        "mean": np.mean(values),
        "std": np.std(values),
        "max": np.max(values),
        "min": np.min(values),
        "high_count": np.sum(values > 3),
        "low_count": np.sum(values < -3),
    }


def identify_potential_trisomies(zscore_df, threshold=3.0):
    """Identify potential trisomies based on Z-scores"""
    potential_trisomies = []

    # Check diagonal values (self-comparison should be close to 0)
    # High diagonal values might indicate issues
    # Check off-diagonal values for consistent high scores

    for i, chrom in enumerate(zscore_df.index):
        # Get Z-scores for this chromosome vs all others
        chrom_scores = zscore_df.iloc[i, :]
        high_scores = chrom_scores[chrom_scores > threshold]

        if (
            len(high_scores) > len(zscore_df.columns) * 0.7
        ):  # If >70% of comparisons are high
            max_score = high_scores.max()
            potential_trisomies.append((chrom, max_score))

    return sorted(potential_trisomies, key=lambda x: x[1], reverse=True)


def create_prizm_summary_report(sample_name, results, output_dir, trisomy_results=None):
    """Create PRIZM analysis summary report"""

    report_file = f"{output_dir}/{sample_name}.prizm_summary.txt"

    try:
        # Ensure directory exists
        os.makedirs(output_dir, exist_ok=True)

        with open(report_file, "w") as f:
            f.write("PRIZM Analysis Summary Report\n")
            f.write(f"Sample: {sample_name}\n")
            f.write(f"Analysis Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("=" * 50 + "\n\n")

            f.write(f"Chromosome Z-scores Matrix Shape: {results.zscore_chr.shape}\n")
            f.write(f"10mb Z-scores Matrix Shape: {results.zscore_10mb.shape}\n")
            f.write(
                f"10mb_all Z-scores Matrix Shape: {results.zscore_10mb_all.shape}\n"
            )
            f.write(f"Total Chromosomes Analyzed: {len(results.chr_index_dict)}\n")
            f.write(f"Total 10mb Bins Analyzed: {results.row_count}\n\n")

            # Calculate summary statistics
            chr_zscore_stats = calculate_zscore_stats(results.zscore_chr)
            mb10_zscore_stats = calculate_zscore_stats(results.zscore_10mb)
            mb10_all_zscore_stats = calculate_zscore_stats(results.zscore_10mb_all)

            f.write("Chromosome-level Z-score Statistics:\n")
            f.write(f"  Mean: {chr_zscore_stats['mean']:.4f}\n")
            f.write(f"  Std Dev: {chr_zscore_stats['std']:.4f}\n")
            f.write(f"  Max: {chr_zscore_stats['max']:.4f}\n")
            f.write(f"  Min: {chr_zscore_stats['min']:.4f}\n")
            f.write(f"  Values > 3: {chr_zscore_stats['high_count']}\n")
            f.write(f"  Values < -3: {chr_zscore_stats['low_count']}\n\n")

            f.write("10mb-level Z-score Statistics:\n")
            f.write(f"  Mean: {mb10_zscore_stats['mean']:.4f}\n")
            f.write(f"  Std Dev: {mb10_zscore_stats['std']:.4f}\n")
            f.write(f"  Max: {mb10_zscore_stats['max']:.4f}\n")
            f.write(f"  Min: {mb10_zscore_stats['min']:.4f}\n")
            f.write(f"  Values > 3: {mb10_zscore_stats['high_count']}\n")
            f.write(f"  Values < -3: {mb10_zscore_stats['low_count']}\n\n")

            f.write("10mb_all-level Z-score Statistics:\n")
            f.write(f"  Mean: {mb10_all_zscore_stats['mean']:.4f}\n")
            f.write(f"  Std Dev: {mb10_all_zscore_stats['std']:.4f}\n")
            f.write(f"  Max: {mb10_all_zscore_stats['max']:.4f}\n")
            f.write(f"  Min: {mb10_all_zscore_stats['min']:.4f}\n")
            f.write(f"  Values > 3: {mb10_all_zscore_stats['high_count']}\n")
            f.write(f"  Values < -3: {mb10_all_zscore_stats['low_count']}\n\n")

            # Trisomy detection results
            if trisomy_results:
                f.write("Statistical Trisomy Detection Results:\n")
                for chrom, result in trisomy_results.items():
                    decision = result["decision"]
                    confidence = result["confidence"]
                    f.write(f"  {chrom}: {decision} (confidence: {confidence:.1%})\n")
                f.write("\n")

            # Identify potential trisomies using statistical approach
            potential_trisomies = identify_potential_trisomies(results.zscore_chr)
            if potential_trisomies:
                f.write("Potential Trisomies (Statistical Z-score > 3):\n")
                for chrom, zscore in potential_trisomies:
                    f.write(f"  {chrom}: {zscore:.4f}\n")
            else:
                f.write("No potential trisomies detected (all Z-scores <= 3)\n")

        logger.info(f"PRIZM summary report created: {report_file}")

    except Exception as e:
        logger.error(f"Error creating PRIZM summary report: {e}")


def run_multiple_prizm_analysis(
    sample_name, gender, labcode, config, analysis_dir, data_dir
):
    """
    Run comprehensive PRIZM analysis for orig, fetus, and mom with gender-specific references

    This is the main function for pipeline integration.
    It handles all necessary PRIZM analyses including wig conversion, reference file validation,
    and multiple analysis types (orig, fetus, mom, plus 10mb_all variants).

    Args:
        sample_name: Sample identifier
        gender: Gender determined from fetal fraction ('MALE' or 'FEMALE')
        labcode: Laboratory code for reference files
        config: Configuration dictionary
        analysis_dir: Analysis directory path
        data_dir: Data directory path for reference files

    Returns:
        bool: True if all analyses completed successfully
    """

    logger.info(
        f"=== Starting comprehensive PRIZM analysis for {sample_name} (Gender: {gender}) ==="
    )

    # Determine gender prefix for reference files
    gender_prefix = "male" if gender == "MALE" else "female"

    # Create PRIZM output directory
    hmmcopy_output_dir = f"{analysis_dir}/{sample_name}/Output_hmmcopy"
    prizm_output_dir = f"{analysis_dir}/{sample_name}/Output_PRIZM"
    os.makedirs(prizm_output_dir, exist_ok=True)

    # EZD (HMMcopy normalization) input file paths
    ezd_input_paths = {
        "orig": f"{hmmcopy_output_dir}/{sample_name}.of_orig.10mb.wig.Normalization.txt",
        "fetus": f"{hmmcopy_output_dir}/{sample_name}.of_fetus.10mb.wig.Normalization.txt",
        "mom": f"{hmmcopy_output_dir}/{sample_name}.of_mom.10mb.wig.Normalization.txt",
    }

    # Define analysis types and their corresponding input files
    analyses = [
        {
            "type": "orig",
            "input_file": ezd_input_paths["orig"],
            "qc_cutoff": config.get("QC", {}).get("orig_biqc", 4.0),
            "description": "Original (no size filter)",
            "enable_plots": True,
            "enable_10mb_all": True,
            "enable_trisomy_detection": True,  # Enable for all analysis types now
        },
        {
            "type": "fetus",
            "input_file": ezd_input_paths["fetus"],
            "qc_cutoff": config.get("QC", {}).get("fetus_biqc", 4.5),
            "description": "Fetal (size filtered 100-160bp)",
            "enable_plots": True,
            "enable_10mb_all": True,
            "enable_trisomy_detection": True,  # Main analysis for trisomy detection
        },
        {
            "type": "mom",
            "input_file": ezd_input_paths["mom"],
            "qc_cutoff": config.get("QC", {}).get("mom_biqc", 4.5),
            "description": "Maternal (size filtered >160bp)",
            "enable_plots": True,
            "enable_10mb_all": True,
            "enable_trisomy_detection": True,  # For comparison/validation
        },
    ]

    all_success = True

    # Run main analyses (orig, fetus, mom)
    for analysis in analyses:
        analysis_type = analysis["type"]
        input_file = analysis["input_file"]
        qc_cutoff = analysis["qc_cutoff"]
        description = analysis["description"]
        enable_plots = analysis["enable_plots"]
        enable_10mb_all = analysis["enable_10mb_all"]
        enable_trisomy_detection = analysis["enable_trisomy_detection"]

        logger.info(f"Running PRIZM {analysis_type} analysis: {description}")

        # Check if input file exists
        if not os.path.exists(input_file):
            logger.warning(f"Input file not found: {input_file}")
            logger.warning(f"Skipping PRIZM {analysis_type} analysis")
            continue

        # Create analysis-specific output directory
        analysis_output_dir = f"{prizm_output_dir}/{analysis_type}"
        os.makedirs(analysis_output_dir, exist_ok=True)

        # Define reference file paths based on gender and analysis type
        # gx-nipt reference layout: <ref_dir>/labs/<labcode>/PRIZM/<analysis_type>/
        ref_base_path = f"{data_dir}/labs/{labcode}/PRIZM/{analysis_type}"

        reference_files = {
            "mean_file": f"{ref_base_path}/{gender_prefix}_mean.csv",
            "sd_file": f"{ref_base_path}/{gender_prefix}_sd.csv",
            "mean_10mb_file": f"{ref_base_path}/{gender_prefix}_10mb_mean.csv",
            "sd_10mb_file": f"{ref_base_path}/{gender_prefix}_10mb_sd.csv",
        }

        # Add 10mb_all reference files if needed
        if enable_10mb_all:
            reference_files.update(
                {
                    "mean_10mb_all_file": f"{ref_base_path}/{gender_prefix}_10mb_all_mean.csv",
                    "sd_10mb_all_file": f"{ref_base_path}/{gender_prefix}_10mb_all_sd.csv",
                }
            )

        # Check if all reference files exist
        missing_refs = []
        for ref_name, ref_path in reference_files.items():
            if not os.path.exists(ref_path):
                missing_refs.append(ref_path)

        if missing_refs:
            logger.warning(
                f"Missing reference files for PRIZM {analysis_type}: {missing_refs}"
            )
            logger.warning(f"Skipping PRIZM {analysis_type} analysis")
            all_success = False
            continue

        try:
            # Run PRIZM analysis
            logger.info(f"Executing PRIZM analysis for {sample_name}_{analysis_type}")

            results = run_prizm_analysis(
                count_file_10mb=input_file,
                mean_file=reference_files["mean_file"],
                sd_file=reference_files["sd_file"],
                mean_10mb_file=reference_files["mean_10mb_file"],
                sd_10mb_file=reference_files["sd_10mb_file"],
                mean_10mb_all_file=reference_files.get("mean_10mb_all_file"),
                sd_10mb_all_file=reference_files.get("sd_10mb_all_file"),
                sample_name=sample_name,
                qc_cutoff=qc_cutoff,
                skip_plots=not enable_plots,
                output_dir=analysis_output_dir,
                enable_10mb_all=enable_10mb_all,
                group=analysis_type,  # Pass the group for proper filename generation
                dpi=config.get("PRIZM", {}).get("resolution_dpi", 200),
            )

            if results:
                # Save Z-score results to analysis-specific directory
                save_zscore_results(
                    results.zscore_chr,
                    f"{analysis_output_dir}/{sample_name}_{analysis_type}.prizm.zscore.txt",
                )
                save_zscore_results(
                    results.zscore_10mb,
                    f"{analysis_output_dir}/{sample_name}_{analysis_type}.prizm.zscore.10mb.txt",
                )

                if enable_10mb_all and results.zscore_10mb_all is not None:
                    save_zscore_results(
                        results.zscore_10mb_all,
                        f"{analysis_output_dir}/{sample_name}_{analysis_type}.prizm.zscore.10mb_all.txt",
                    )

                logger.info(f"PRIZM {analysis_type} analysis completed successfully")
                logger.info(
                    f"  - Chromosome Z-scores matrix: {results.zscore_chr.shape}"
                )
                logger.info(f"  - 10mb Z-scores matrix: {results.zscore_10mb.shape}")
                if enable_10mb_all and results.zscore_10mb_all is not None:
                    logger.info(
                        f"  - 10mb_all Z-scores matrix: {results.zscore_10mb_all.shape}"
                    )
                logger.info(f"  - Results saved to: {analysis_output_dir}")

                # Perform statistical trisomy detection if enabled
                trisomy_results = None
                if enable_trisomy_detection:
                    try:
                        trisomy_results = run_statistical_trisomy_detection(
                            results.zscore_chr,
                            analysis_output_dir,
                            f"{sample_name}_{analysis_type}",
                        )
                        logger.info(
                            f"Statistical trisomy detection completed for {analysis_type}"
                        )
                    except Exception as e:
                        logger.error(
                            f"Error in statistical trisomy detection for {analysis_type}: {e}"
                        )

                # Create summary report for this analysis
                create_prizm_summary_report(
                    f"{sample_name}_{analysis_type}",
                    results,
                    analysis_output_dir,
                    trisomy_results,
                )
            else:
                logger.error(f"PRIZM {analysis_type} analysis returned no results")
                all_success = False

        except Exception as e:
            logger.error(f"Error in PRIZM {analysis_type} analysis: {str(e)}")
            all_success = False

    # Final summary
    if all_success:
        logger.info("=== All PRIZM analyses completed successfully ===")
    else:
        logger.warning("=== Some PRIZM analyses failed or were skipped ===")

    return all_success


# ============================================================================
# MAIN FUNCTION
# ============================================================================


def main():
    """Main function"""

    # Setup argument parser
    parser = argparse.ArgumentParser(
        description="PRIZM (Prenatal Risk Z-score Matrix) Calculator with Plotting",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Input arguments
    input_group = parser.add_argument_group("Input")
    input_group.add_argument(
        "-c10",
        dest="count_file_10mb",
        required=True,
        type=lambda x: file_exists(parser, x),
        help="10mb count file",
    )
    input_group.add_argument(
        "-m",
        dest="mean",
        required=True,
        type=lambda x: file_exists(parser, x),
        help="Reference mean file",
    )
    input_group.add_argument(
        "-s",
        dest="sd",
        required=True,
        type=lambda x: file_exists(parser, x),
        help="Reference standard deviation file",
    )
    input_group.add_argument(
        "-m10",
        dest="mean_10mb",
        required=True,
        type=lambda x: file_exists(parser, x),
        help="10mb reference mean file",
    )
    input_group.add_argument(
        "-s10",
        dest="sd_10mb",
        required=True,
        type=lambda x: file_exists(parser, x),
        help="10mb reference standard deviation file",
    )
    input_group.add_argument(
        "-m10all",
        dest="mean_10mb_all",
        type=lambda x: file_exists(parser, x),
        help="10mb_all reference mean file",
    )
    input_group.add_argument(
        "-s10all",
        dest="sd_10mb_all",
        type=lambda x: file_exists(parser, x),
        help="10mb_all reference standard deviation file",
    )
    input_group.add_argument(
        "-q", dest="qc_cutoff", type=float, default=3.0, help="QC cutoff threshold"
    )
    input_group.add_argument(
        "-single",
        dest="single_output",
        default="N",
        choices=["Y", "N"],
        help="Single output mode",
    )
    input_group.add_argument(
        "--sample_name", dest="sample_name", help="Sample name for plot titles"
    )
    input_group.add_argument(
        "--skip_plots",
        dest="skip_plots",
        action="store_true",
        help="Skip plot generation",
    )
    input_group.add_argument(
        "--enable_10mb_all",
        dest="enable_10mb_all",
        action="store_true",
        help="Enable 10mb_all Z-score calculation",
    )
    input_group.add_argument(
        "--enable_trisomy_detection",
        dest="enable_trisomy_detection",
        action="store_true",
        help="Enable trisomy detection",
    )

    # Options
    options_group = parser.add_argument_group("Options")
    options_group.add_argument(
        "--debug", action="store_true", default=False, help="Enable debug mode"
    )
    options_group.add_argument(
        "-v", "--version", action="version", version="%(prog)s 1.0.0"
    )

    args = parser.parse_args()

    # Setup logging and configuration
    global logger
    logger = setup_logging(args.debug)

    config = PRIZMConfig(
        qc_cutoff=args.qc_cutoff, single_output=(args.single_output == "Y")
    )

    # Extract sample name from count file if not provided
    sample_name = args.sample_name
    if not sample_name:
        sample_name = os.path.basename(args.count_file_10mb).split(".")[0]

    try:
        # Main PRIZM Z-score calculation
        if args.enable_10mb_all and args.mean_10mb_all and args.sd_10mb_all:
            results = calculate_prizm_zscores(
                args.count_file_10mb,
                args.mean,
                args.sd,
                args.mean_10mb,
                args.sd_10mb,
                args.mean_10mb_all,
                args.sd_10mb_all,
            )
        else:
            # Basic calculation without 10mb_all
            count_data = load_count_data(args.count_file_10mb)
            chr_summary, mb10_data = prepare_chromosome_data(count_data)

            # Calculate chromosome-level Z-scores
            chr_normalized = calculate_normalized_matrix(chr_summary)
            chr_zscore_matrix = calculate_zscore_from_normalized(
                chr_normalized, args.mean, args.sd
            )

            chr_index_dict = {i: f"chr{i + 1}" for i in range(22)}
            chr_index_dict.update({22: "chrX", 23: "chrY"})

            chr_zscore_df = create_zscore_dataframe(
                chr_zscore_matrix, list(chr_index_dict.values())
            )

            # Calculate 10mb-level Z-scores
            mb10_normalized = calculate_normalized_matrix(mb10_data)
            mb10_zscore_matrix = calculate_zscore_from_normalized(
                mb10_normalized, args.mean_10mb, args.sd_10mb
            )

            mb10_index_dict = {
                i: mb10_data.iloc[i]["chr"] for i in range(len(mb10_data))
            }
            mb10_zscore_df = create_zscore_dataframe(
                mb10_zscore_matrix, list(mb10_index_dict.values())
            )

            results = PRIZMResult(
                zscore_chr=chr_zscore_df,
                zscore_10mb=mb10_zscore_df,
                zscore_10mb_all=None,
                chr_index_dict=chr_index_dict,
                mb10_index_dict=mb10_index_dict,
                row_count=len(mb10_data),
            )

        # QC analysis
        filename = os.path.basename(args.count_file_10mb).split(".10mb")[0]
        output_dir = (
            ""
            if config.single_output
            else "/".join(args.count_file_10mb.split("/")[:-1])
        )
        qc_file = (
            os.path.join(output_dir, f"{filename}.prizm.qc.txt")
            if output_dir
            else f"{filename}.prizm.qc.txt"
        )

        qc_passed = perform_qc_analysis(results.zscore_chr, config.qc_cutoff, qc_file)

        # Generate plots if not skipped
        if not args.skip_plots:
            plots_output_dir = output_dir if output_dir else "."
            generate_prizm_plots(results, plots_output_dir, sample_name)

        # Perform statistical trisomy detection if enabled
        if args.enable_trisomy_detection:
            trisomy_results = run_statistical_trisomy_detection(
                results.zscore_chr, output_dir if output_dir else ".", sample_name
            )

        logger.info("PRIZM analysis completed successfully!")
        logger.info(f"Chromosome Z-scores shape: {results.zscore_chr.shape}")
        logger.info(f"10mb Z-scores shape: {results.zscore_10mb.shape}")
        if results.zscore_10mb_all is not None:
            logger.info(f"10mb_all Z-scores shape: {results.zscore_10mb_all.shape}")

    except Exception as e:
        logger.error(f"PRIZM analysis failed: {str(e)}")
        sys.exit(1)


# ============================================================================
# CALLABLE FUNCTIONS FOR PIPELINE INTEGRATION
# ============================================================================


def run_prizm_analysis(
    count_file_10mb: str,
    mean_file: str,
    sd_file: str,
    mean_10mb_file: str,
    sd_10mb_file: str,
    mean_10mb_all_file: str = None,
    sd_10mb_all_file: str = None,
    sample_name: str = None,
    qc_cutoff: float = 3.0,
    skip_plots: bool = False,
    output_dir: str = None,
    enable_10mb_all: bool = True,
    group: str = "",
    dpi=200,
) -> PRIZMResult:
    """
    Run PRIZM analysis programmatically from other scripts

    This is the core analysis function that performs the actual Z-score calculations.
    It's called by run_multiple_prizm_analysis for each analysis type.

    Args:
        count_file_10mb: Path to 10mb count file
        mean_file: Path to reference mean file
        sd_file: Path to reference standard deviation file
        mean_10mb_file: Path to 10mb reference mean file
        sd_10mb_file: Path to 10mb reference standard deviation file
        mean_10mb_all_file: Path to 10mb_all reference mean file
        sd_10mb_all_file: Path to 10mb_all reference standard deviation file
        sample_name: Sample name for plot titles
        qc_cutoff: QC cutoff threshold
        skip_plots: Whether to skip plot generation
        output_dir: Output directory for QC files and plots
        enable_10mb_all: Whether to calculate 10mb_all Z-scores
        group: Analysis group (orig, fetus, mom) for filename generation

    Returns:
        PRIZMResult object containing analysis results
    """

    # Extract sample name if not provided
    if not sample_name:
        sample_name = os.path.basename(count_file_10mb).split(".")[0]

    logger.info(f"Starting PRIZM analysis for sample: {sample_name}")

    # Main PRIZM Z-score calculation
    if enable_10mb_all and mean_10mb_all_file and sd_10mb_all_file:
        results = calculate_prizm_zscores(
            count_file_10mb,
            mean_file,
            sd_file,
            mean_10mb_file,
            sd_10mb_file,
            mean_10mb_all_file,
            sd_10mb_all_file,
        )
    else:
        # Basic calculation without 10mb_all
        count_data = load_count_data(count_file_10mb)
        chr_summary, mb10_data = prepare_chromosome_data(count_data)

        # Calculate chromosome-level Z-scores
        chr_normalized = calculate_normalized_matrix(chr_summary)
        chr_zscore_matrix = calculate_zscore_from_normalized(
            chr_normalized, mean_file, sd_file
        )

        chr_index_dict = {i: f"chr{i + 1}" for i in range(22)}
        chr_index_dict.update({22: "chrX", 23: "chrY"})

        chr_zscore_df = create_zscore_dataframe(
            chr_zscore_matrix, list(chr_index_dict.values())
        )

        # Calculate 10mb-level Z-scores
        mb10_normalized = calculate_normalized_matrix(mb10_data)
        mb10_zscore_matrix = calculate_zscore_from_normalized(
            mb10_normalized, mean_10mb_file, sd_10mb_file
        )

        mb10_index_dict = {i: mb10_data.iloc[i]["chr"] for i in range(len(mb10_data))}
        mb10_zscore_df = create_zscore_dataframe(
            mb10_zscore_matrix, list(mb10_index_dict.values())
        )

        results = PRIZMResult(
            zscore_chr=chr_zscore_df,
            zscore_10mb=mb10_zscore_df,
            zscore_10mb_all=None,
            chr_index_dict=chr_index_dict,
            mb10_index_dict=mb10_index_dict,
            row_count=len(mb10_data),
        )

    # Determine output directory for QC and plots
    if output_dir is None:
        output_dir = "/".join(count_file_10mb.split("/")[:-1])

    # QC analysis
    filename = os.path.basename(count_file_10mb).split(".10mb")[0]
    qc_file = os.path.join(output_dir, f"{filename}.prizm.qc.txt")

    qc_passed = perform_qc_analysis(results.zscore_chr, qc_cutoff, qc_file)

    # Generate plots if not skipped
    if not skip_plots:
        generate_prizm_plots(results, output_dir, sample_name, group, dpi)

    logger.info(f"PRIZM analysis completed for sample: {sample_name}")

    return results


# =========================================================
#  Adapter CLI: Nextflow-compatible entry point
#
#  Used by modules/prizm.nf. Coexists with the legacy `main()`
#  which consumes `-c10/-m/-s/-m10/-s10/...`. The adapter selects
#  reference files from a flat ref-dir based on the gender parsed
#  from the FF/gender file and then delegates to run_prizm_analysis().
# =========================================================
def _parse_gender_from_file(path: str, default: str = "female") -> str:
    """Return 'male' or 'female' prefix from a gender.txt / fetal_fraction.txt file.

    Expected formats:
      - gender.txt  (preferred): contains a line `final_gender\t<MALE|FEMALE>`
      - fetal_fraction.txt      : no gender — raises ValueError

    Also accepts lowercase/mixed case and karyotype notation (XY/XX).
    When final_gender is UNKNOWN (low FF sample), falls back to *default* ('female').
    Returns 'male' or 'female'.
    """
    if not path or not os.path.exists(path):
        raise ValueError(f"Gender file not found: {path}")

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            # final_gender<TAB>MALE|FEMALE  (or karyotype XY|XX written by ff_gender_improved.py)
            if line.lower().startswith("final_gender"):
                parts = re.split(r"[\t\s]+", line, maxsplit=1)
                if len(parts) == 2:
                    val = parts[1].strip().upper()
                    if val in ("MALE", "M", "XY"):
                        return "male"
                    if val in ("FEMALE", "F", "XX"):
                        return "female"
                    if val in ("UNKNOWN", "NA", "N/A", ""):
                        logger.warning(
                            f"final_gender is '{val}' in {path} (low FF sample?). "
                            f"Defaulting to '{default}' for PRIZM reference selection."
                        )
                        return default
    raise ValueError(
        f"Could not parse final_gender from {path}. "
        "Expected a line like 'final_gender\\tMALE' or 'final_gender\\tFEMALE' (or XY/XX)."
    )


def _nf_main():
    import argparse
    import json

    parser = argparse.ArgumentParser(
        description=(
            "PRIZM adapter CLI (Nextflow). Runs PRIZM for a single "
            "(sample, group) and writes results flat into --outdir."
        )
    )
    parser.add_argument("--sample", required=True)
    parser.add_argument(
        "--group", required=True, choices=["orig", "fetus", "mom"]
    )
    parser.add_argument(
        "--count-file",
        required=True,
        help="HMMcopy 10mb normalization TXT "
        "(e.g. <sample>.of_<group>.10mb.wig.Normalization.txt)",
    )
    parser.add_argument(
        "--ref-dir",
        required=True,
        help=(
            "PRIZM reference directory for this group, expected to contain "
            "{male|female}_{mean,sd,10mb_mean,10mb_sd,10mb_all_mean,10mb_all_sd}.csv. "
            "Typically <top_ref>/labs/<labcode>/PRIZM/<group>."
        ),
    )
    parser.add_argument(
        "--gender-file",
        default=None,
        help="gender.txt from GENDER_DECISION (contains 'final_gender\\tMALE|FEMALE')",
    )
    parser.add_argument(
        "--ff-file",
        default=None,
        help=(
            "Legacy alias for --gender-file. Must contain a 'final_gender' line "
            "for gender resolution."
        ),
    )
    parser.add_argument("--labcode", default=None, help="Lab code (informational)")
    parser.add_argument(
        "--outdir", required=True, help="Flat output dir for QC/plots"
    )
    parser.add_argument(
        "--config",
        default=None,
        help="pipeline_config.json (optional; used for QC cutoff / plot DPI)",
    )
    args = parser.parse_args()

    global logger
    logger = setup_logging(False)

    gender_src = args.gender_file or args.ff_file
    if not gender_src:
        logger.error("--gender-file (or --ff-file with gender info) is required.")
        sys.exit(2)

    try:
        gender_prefix = _parse_gender_from_file(gender_src)
    except ValueError as e:
        logger.error(str(e))
        sys.exit(2)

    logger.info(
        f"[PRIZM adapter] {args.sample}/{args.group} gender={gender_prefix}"
    )

    cfg = {}
    if args.config and os.path.exists(args.config):
        try:
            with open(args.config) as fh:
                cfg = json.load(fh)
        except Exception as e:
            logger.warning(f"Failed to parse config {args.config}: {e}")

    qc_cutoff = (
        cfg.get("QC", {}).get(f"{args.group}_biqc")
        or cfg.get("QC", {}).get("biqc")
        or 3.0
    )
    dpi = cfg.get("PRIZM", {}).get("resolution_dpi", 200)

    ref_dir = args.ref_dir.rstrip("/")
    mean_file         = f"{ref_dir}/{gender_prefix}_mean.csv"
    sd_file           = f"{ref_dir}/{gender_prefix}_sd.csv"
    mean_10mb_file    = f"{ref_dir}/{gender_prefix}_10mb_mean.csv"
    sd_10mb_file      = f"{ref_dir}/{gender_prefix}_10mb_sd.csv"
    mean_10mb_all_fp  = f"{ref_dir}/{gender_prefix}_10mb_all_mean.csv"
    sd_10mb_all_fp    = f"{ref_dir}/{gender_prefix}_10mb_all_sd.csv"

    missing = [p for p in (mean_file, sd_file, mean_10mb_file, sd_10mb_file)
               if not os.path.exists(p)]
    if missing:
        logger.error(f"[PRIZM adapter] Missing required reference files: {missing}")
        sys.exit(3)

    enable_10mb_all = (
        os.path.exists(mean_10mb_all_fp) and os.path.exists(sd_10mb_all_fp)
    )
    if not enable_10mb_all:
        logger.warning(
            f"[PRIZM adapter] 10mb_all references not found under {ref_dir}; "
            "proceeding without 10mb_all computation."
        )

    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    results = run_prizm_analysis(
        count_file_10mb=args.count_file,
        mean_file=mean_file,
        sd_file=sd_file,
        mean_10mb_file=mean_10mb_file,
        sd_10mb_file=sd_10mb_file,
        mean_10mb_all_file=mean_10mb_all_fp if enable_10mb_all else None,
        sd_10mb_all_file=sd_10mb_all_fp if enable_10mb_all else None,
        sample_name=args.sample,
        qc_cutoff=qc_cutoff,
        skip_plots=False,
        output_dir=outdir,
        enable_10mb_all=enable_10mb_all,
        group=args.group,
        dpi=dpi,
    )

    if results is None:
        logger.error("[PRIZM adapter] run_prizm_analysis returned None")
        sys.exit(1)

    try:
        run_statistical_trisomy_detection(
            results.zscore_chr, outdir, f"{args.sample}_{args.group}"
        )
    except Exception as e:
        logger.warning(
            f"[PRIZM adapter] Statistical trisomy detection skipped: {e}"
        )

    try:
        create_prizm_summary_report(
            f"{args.sample}_{args.group}", results, outdir
        )
    except Exception as e:
        logger.warning(f"[PRIZM adapter] Summary report skipped: {e}")

    logger.info(
        f"[PRIZM adapter] {args.sample}/{args.group} -> {outdir} done."
    )


if __name__ == "__main__":
    # Detect Nextflow-style invocation by the presence of long-form
    # arguments that the legacy CLI never used.
    argv = sys.argv[1:]
    _nf_flags = ("--sample", "--group", "--count-file", "--ref-dir", "--outdir")
    if any(flag in argv for flag in _nf_flags):
        _nf_main()
    else:
        main()
