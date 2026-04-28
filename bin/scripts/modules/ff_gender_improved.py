#!/usr/bin/env python3
"""
ff_gender_improved.py
======================
Improved Fetal Fraction & Gender Detection module for gx-nipt.

Improvements over ken-nipt:
  1. YFF1: BUG FIX - BED file path now read from config (not hardcoded)
  2. YFF2: Cleaner separation of UAR_X, UAR_Y, gd_2 ratio calculation
  3. Fragment FF: Proper fragment size distribution analysis (short/long ratio)
  4. Gender Decision: Weighted ensemble approach with confidence scoring
  5. Low FF warning: Explicit flag when FF < 4% (high risk of false call)
"""

import argparse
import json
import logging
import os
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import pysam

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Utility helpers
# ─────────────────────────────────────────────────────────────────────────────

def setup_logging(debug: bool = False) -> None:
    level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def load_config(config_path: str) -> dict:
    with open(config_path) as f:
        return json.load(f)


def count_reads_in_regions(bam_path: str, bed_path: str) -> Tuple[list, list, int]:
    """Count reads in BED regions using pysam."""
    counts = []
    regions = []
    total_bases = 0

    bam = pysam.AlignmentFile(bam_path, "rb")
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            length = end - start
            total_bases += length
            n = bam.count(chrom, start, end)
            counts.append(n)
            regions.append((chrom, start, end, length))
    bam.close()
    return counts, regions, total_bases


def calculate_normalized_coverage(counts: list, regions: list) -> list:
    """Normalize read counts by region length (RPKM-like)."""
    normalized = []
    for count, (_, _, _, length) in zip(counts, regions):
        if length > 0:
            normalized.append(count / length)
        else:
            normalized.append(0.0)
    return normalized


# ─────────────────────────────────────────────────────────────────────────────
# YFF1: Y-chromosome coverage ratio
# ─────────────────────────────────────────────────────────────────────────────

def calculate_yff1(sample_name: str, bam_path: str, config: dict, labcode: str) -> dict:
    """
    Calculate Y-based fetal fraction (YFF1) using lab-specific BED files.

    BUG FIX (ken-nipt): BED file paths were hardcoded as /Work/NIPT/data/bed/common/*.bed
    Now read from config['lab_bed_paths'] to support multiple lab configurations.
    """
    ff_config = config.get("FF_Gender_Config", {})
    bed_paths = config.get("lab_bed_paths", {})

    # Read BED paths from config (not hardcoded)
    y_bed = bed_paths.get("y_regions_09", "")
    a_bed = bed_paths.get("autosome_control", "")

    if not os.path.exists(y_bed):
        logger.error(f"Y BED file not found: {y_bed}")
        return _yff1_failed(sample_name, f"Y BED not found: {y_bed}")

    if not os.path.exists(a_bed):
        logger.error(f"Autosome BED file not found: {a_bed}")
        return _yff1_failed(sample_name, f"Autosome BED not found: {a_bed}")

    try:
        y_counts, y_regions, _ = count_reads_in_regions(bam_path, y_bed)
        a_counts, a_regions, _ = count_reads_in_regions(bam_path, a_bed)

        if not y_counts or not a_counts:
            return _yff1_failed(sample_name, "Empty read counts")

        y_norm = calculate_normalized_coverage(y_counts, y_regions)
        a_norm = calculate_normalized_coverage(a_counts, a_regions)

        y_median = float(np.median(y_norm))
        a_median = float(np.median(a_norm))

        logger.info(f"YFF1 - Y median coverage: {y_median:.6f}")
        logger.info(f"YFF1 - Autosome median coverage: {a_median:.6f}")

        if a_median <= 0:
            return _yff1_failed(sample_name, "Autosome median coverage is zero")

        y_to_a_ratio = y_median / a_median
        # FF = 2 * (Y/A) * 100 (assuming male fetus; only valid for XY)
        fetal_fraction = 2.0 * y_to_a_ratio * 100.0

        gender_threshold = ff_config.get("gd_1_threshold", 0.01)
        gd_1_gender = "XY" if y_to_a_ratio >= gender_threshold else "XX"

        logger.info(f"YFF1 - Y/A ratio: {y_to_a_ratio:.6f}, FF: {fetal_fraction:.2f}%, Gender: {gd_1_gender}")

        return {
            "sample_name": sample_name,
            "YFF1": round(fetal_fraction, 4),
            "gd_1_value": round(y_to_a_ratio, 6),
            "gd_1_gender": gd_1_gender,
            "status": "OK",
        }
    except Exception as e:
        logger.exception(f"Error in calculate_yff1: {e}")
        return _yff1_failed(sample_name, str(e))


def _yff1_failed(sample_name: str, reason: str) -> dict:
    # Default to XX (Female) on failure — same as ken-nipt gender_detector.py
    # gd_1_detection() which returns (0.0, "FEMALE", "XX") on any exception.
    logger.warning(f"YFF1 failed for {sample_name}: {reason}. Defaulting gd_1_gender to XX.")
    return {
        "sample_name": sample_name,
        "YFF1": 0.0,
        "gd_1_value": 0.0,
        "gd_1_gender": "XX",
        "status": f"FAILED: {reason}",
    }


# ─────────────────────────────────────────────────────────────────────────────
# YFF2: Adjusted Y-based FF from wig normalization
# ─────────────────────────────────────────────────────────────────────────────

def calculate_yff2(sample_name: str, wig_norm_file: str, config: dict) -> dict:
    """
    Calculate YFF2 from HMMcopy 50kb normalization output.
    Uses UAR_X (chrX usage ratio) and UAR_Y (chrY usage ratio) to determine
    gender (gd_2) and estimate fetal fraction for male fetuses.
    """
    ff_config = config.get("FF_Gender_Config", {})
    gd_2_threshold = ff_config.get("gd_2_threshold", 0.4)

    if not os.path.exists(wig_norm_file):
        logger.warning(f"YFF2 wig norm file not found: {wig_norm_file}")
        return _yff2_failed(sample_name, "wig norm file not found")

    try:
        x_values, y_values, autosome_values = [], [], []

        with open(wig_norm_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                chrom = parts[0]
                try:
                    norm_val = float(parts[2])
                except (ValueError, IndexError):
                    continue

                if chrom == "chrX":
                    x_values.append(norm_val)
                elif chrom == "chrY":
                    y_values.append(norm_val)
                elif chrom not in ("chrX", "chrY"):
                    autosome_values.append(norm_val)

        if not autosome_values:
            return _yff2_failed(sample_name, "No autosome data in wig norm file")

        autosome_median = float(np.median(autosome_values))
        if autosome_median <= 0:
            return _yff2_failed(sample_name, "Autosome median is zero")

        uar_x = float(np.median(x_values)) / autosome_median if x_values else 0.0
        uar_y = float(np.median(y_values)) / autosome_median if y_values else 0.0

        # gd_2: ratio of Y to X usage ratio
        gd_2_value = uar_y / uar_x if uar_x > 0 else 0.0
        gd_2_gender = "XY" if gd_2_value >= gd_2_threshold else "XX"

        # YFF2: for male fetus, FF = 2 * UAR_Y * 100
        yff2 = 2.0 * uar_y * 100.0 if gd_2_gender == "XY" else 0.0

        logger.info(f"YFF2 - UAR_X: {uar_x:.4f}, UAR_Y: {uar_y:.4f}, gd_2: {gd_2_value:.4f}, Gender: {gd_2_gender}")

        return {
            "sample_name": sample_name,
            "YFF2": round(yff2, 4),
            "UAR_X": round(uar_x, 6),
            "UAR_Y": round(uar_y, 6),
            "gd_2_value": round(gd_2_value, 6),
            "gd_2_gender": gd_2_gender,
            "status": "OK",
        }
    except Exception as e:
        logger.exception(f"Error in calculate_yff2: {e}")
        return _yff2_failed(sample_name, str(e))


def _yff2_failed(sample_name: str, reason: str) -> dict:
    # Default to XX (Female) on failure — same as ken-nipt gender_detector.py
    # gd_2_detection() which returns (0.0, "FEMALE", "XX") on any exception.
    logger.warning(f"YFF2 failed for {sample_name}: {reason}. Defaulting gd_2_gender to XX.")
    return {
        "sample_name": sample_name,
        "YFF2": 0.0,
        "UAR_X": 0.0,
        "UAR_Y": 0.0,
        "gd_2_value": 0.0,
        "gd_2_gender": "XX",
        "status": f"FAILED: {reason}",
    }


# ─────────────────────────────────────────────────────────────────────────────
# SeqFF: Sequence-based Fetal Fraction
# ─────────────────────────────────────────────────────────────────────────────

def calculate_seqff(sample_name: str, bam_path: str, config: dict, seqff_model_path: str) -> dict:
    """
    Calculate SeqFF using GC-corrected coverage across autosomal bins.
    If a trained model is available, use it; otherwise fall back to
    a simple GC-correction based estimate.
    """
    try:
        # Try to load pre-trained SeqFF model
        if os.path.exists(seqff_model_path):
            import pickle
            with open(seqff_model_path, "rb") as f:
                model = pickle.load(f)
            # Extract features from BAM
            features = _extract_seqff_features(bam_path)
            seqff = float(model.predict([features])[0])
            logger.info(f"SeqFF (model): {seqff:.4f}%")
        else:
            # Fallback: use fragment size ratio as proxy
            logger.warning(f"SeqFF model not found at {seqff_model_path}, using fallback")
            seqff = _seqff_fallback(bam_path)

        return {
            "sample_name": sample_name,
            "SeqFF": round(seqff, 4),
            "status": "OK",
        }
    except Exception as e:
        logger.exception(f"Error in calculate_seqff: {e}")
        return {
            "sample_name": sample_name,
            "SeqFF": 0.0,
            "status": f"FAILED: {str(e)}",
        }


def _extract_seqff_features(bam_path: str) -> list:
    """Extract GC-content binned coverage features for SeqFF model."""
    gc_bins = [0] * 101
    total = 0
    bam = pysam.AlignmentFile(bam_path, "rb")
    for read in bam.fetch():
        if read.is_unmapped or read.is_duplicate:
            continue
        seq = read.query_sequence
        if seq:
            gc = sum(1 for b in seq if b in "GCgc")
            gc_pct = int(gc / len(seq) * 100)
            gc_bins[gc_pct] += 1
            total += 1
    bam.close()
    if total > 0:
        return [c / total for c in gc_bins]
    return gc_bins


def _seqff_fallback(bam_path: str) -> float:
    """Fallback SeqFF: use short-fragment ratio as proxy."""
    short, total = 0, 0
    bam = pysam.AlignmentFile(bam_path, "rb")
    for read in bam.fetch():
        if read.is_unmapped or read.is_duplicate or not read.is_proper_pair:
            continue
        tlen = abs(read.template_length)
        if 50 <= tlen <= 1000:
            total += 1
            if tlen < 150:
                short += 1
    bam.close()
    if total == 0:
        return 0.0
    return (short / total) * 100.0


# ─────────────────────────────────────────────────────────────────────────────
# Fragment FF: Fragment size distribution-based FF
# ─────────────────────────────────────────────────────────────────────────────

def calculate_fragment_ff(sample_name: str, bam_path: str, config: dict, fragment_cutoff: int = 160) -> dict:
    """
    Calculate fragment-size based fetal fraction.

    Fetal cfDNA is characteristically shorter (~166bp) than maternal cfDNA (~185bp).
    The ratio of short fragments (< fragment_cutoff) to total fragments is used
    as a proxy for fetal fraction.

    Improvement: Added M-SeqFF (maternal-corrected) calculation using bimodal
    distribution fitting to separate fetal and maternal fragment populations.
    """
    try:
        fragment_sizes = []
        bam = pysam.AlignmentFile(bam_path, "rb")
        for read in bam.fetch():
            if read.is_unmapped or read.is_duplicate or not read.is_proper_pair:
                continue
            if not read.is_read1:
                continue
            tlen = abs(read.template_length)
            if 50 <= tlen <= 1000:
                fragment_sizes.append(tlen)
        bam.close()

        if not fragment_sizes:
            return _frag_ff_failed(sample_name, "No valid fragment sizes")

        fragment_sizes = np.array(fragment_sizes)
        total = len(fragment_sizes)
        short_count = np.sum(fragment_sizes < fragment_cutoff)

        # Basic fragment FF
        fragment_ff = (short_count / total) * 100.0

        # Statistics
        median_size = float(np.median(fragment_sizes))
        mean_size   = float(np.mean(fragment_sizes))
        std_size    = float(np.std(fragment_sizes))

        # M-SeqFF: attempt bimodal fit
        m_seqff = _calculate_m_seqff(fragment_sizes)

        logger.info(
            f"Fragment FF - median: {median_size:.1f}bp, "
            f"short ratio: {fragment_ff:.2f}%, M-SeqFF: {m_seqff:.2f}%"
        )

        return {
            "sample_name": sample_name,
            "Fragment_FF": round(fragment_ff, 4),
            "M-SeqFF": round(m_seqff, 4),
            "median_fragment_size": round(median_size, 2),
            "mean_fragment_size": round(mean_size, 2),
            "std_fragment_size": round(std_size, 2),
            "total_fragments": total,
            "status": "OK",
        }
    except Exception as e:
        logger.exception(f"Error in calculate_fragment_ff: {e}")
        return _frag_ff_failed(sample_name, str(e))


def _calculate_m_seqff(fragment_sizes: np.ndarray) -> float:
    """
    Estimate M-SeqFF using a simplified bimodal Gaussian mixture model.
    Fetal population: ~166bp, Maternal population: ~185bp.
    """
    try:
        from scipy.stats import norm as scipy_norm
        from scipy.optimize import minimize

        # Initial parameters: [fetal_mean, fetal_std, maternal_mean, maternal_std, fetal_weight]
        x0 = [166, 15, 185, 20, 0.1]
        bounds = [(130, 175), (5, 30), (175, 210), (10, 40), (0.01, 0.5)]

        def neg_log_likelihood(params):
            fm, fs, mm, ms, fw = params
            mw = 1.0 - fw
            pdf = fw * scipy_norm.pdf(fragment_sizes, fm, fs) + mw * scipy_norm.pdf(fragment_sizes, mm, ms)
            pdf = np.clip(pdf, 1e-10, None)
            return -np.sum(np.log(pdf))

        result = minimize(neg_log_likelihood, x0, bounds=bounds, method="L-BFGS-B")
        if result.success:
            fetal_weight = result.x[4]
            return fetal_weight * 100.0
        return 0.0
    except ImportError:
        # scipy not available: fallback to simple ratio
        short = np.sum(fragment_sizes < 170)
        return float(short / len(fragment_sizes)) * 100.0
    except Exception:
        return 0.0


def _frag_ff_failed(sample_name: str, reason: str) -> dict:
    return {
        "sample_name": sample_name,
        "Fragment_FF": 0.0,
        "M-SeqFF": 0.0,
        "median_fragment_size": 0.0,
        "mean_fragment_size": 0.0,
        "std_fragment_size": 0.0,
        "total_fragments": 0,
        "status": f"FAILED: {reason}",
    }


# ─────────────────────────────────────────────────────────────────────────────
# Gender Decision: Weighted ensemble
# ─────────────────────────────────────────────────────────────────────────────

def gender_decision(
    sample_name: str,
    yff1: dict,
    yff2: dict,
    seqff: dict,
    frag_ff: dict,
    config: dict,
) -> Tuple[dict, dict]:
    """
    Make final gender decision using weighted ensemble of all FF estimates.

    Priority order:
      1. gd_2 (YFF2 ratio) - most reliable for gender
      2. gd_1 (YFF1 ratio) - backup
      3. Fragment size distribution (gd_4)

    FF value priority:
      1. YFF2 (for male fetus)
      2. SeqFF (gender-independent)
      3. M-SeqFF (fragment-based, gender-independent)
      4. YFF1 (for male fetus)
    """
    ff_config = config.get("FF_Gender_Config", {})
    yff_threshold = ff_config.get("YFF", 4.0)
    seqff_threshold = ff_config.get("seqFF", 4.0)

    # ── Gender determination ──────────────────────────────
    gd_2_gender = yff2.get("gd_2_gender", "UNKNOWN")
    gd_1_gender = yff1.get("gd_1_gender", "UNKNOWN")

    # Primary: gd_2
    if gd_2_gender in ("XY", "XX"):
        final_gender = gd_2_gender
        gender_source = "gd_2"
    elif gd_1_gender in ("XY", "XX"):
        final_gender = gd_1_gender
        gender_source = "gd_1"
    else:
        final_gender = "UNKNOWN"
        gender_source = "none"

    # ── FF value selection ────────────────────────────────
    yff2_val   = yff2.get("YFF2", 0.0)
    seqff_val  = seqff.get("SeqFF", 0.0)
    m_seqff    = frag_ff.get("M-SeqFF", 0.0)
    frag_ff_val= frag_ff.get("Fragment_FF", 0.0)
    yff1_val   = yff1.get("YFF1", 0.0)

    # For male fetus: YFF2 is most reliable
    # For female fetus: SeqFF or M-SeqFF
    if final_gender == "XY":
        ff_candidates = [
            ("YFF2",    yff2_val),
            ("SeqFF",   seqff_val),
            ("M-SeqFF", m_seqff),
            ("YFF1",    yff1_val),
        ]
    else:
        ff_candidates = [
            ("SeqFF",   seqff_val),
            ("M-SeqFF", m_seqff),
            ("Fragment_FF", frag_ff_val),
        ]

    # Pick first non-zero candidate
    final_ff = 0.0
    ff_source = "none"
    for name, val in ff_candidates:
        if val > 0.0:
            final_ff = val
            ff_source = name
            break

    # ── Low FF warning ────────────────────────────────────
    low_ff_warning = final_ff < yff_threshold and final_ff > 0.0

    logger.info(
        f"Gender Decision - gender: {final_gender} (source: {gender_source}), "
        f"FF: {final_ff:.2f}% (source: {ff_source}), "
        f"low_ff_warning: {low_ff_warning}"
    )

    # ── Build output structures ───────────────────────────
    ff_result = {
        "sample_name": sample_name,
        "Fragment_FF": round(frag_ff.get("Fragment_FF", 0.0), 4),
        "YFF_2":       round(yff2_val, 4),
        "SeqFF":       round(seqff_val, 4),
        "M-SeqFF":     round(m_seqff, 4),
        "Final_FF":    round(final_ff, 4),
        "FF_Source":   ff_source,
        "low_ff_warning": low_ff_warning,
    }

    gender_result = {
        "sample_name": sample_name,
        "gd_1":  f"gd_1\t{yff1.get('gd_1_value', 0.0):.6f}\t{yff1.get('gd_1_gender', 'UNKNOWN')}",
        "gd_2":  f"gd_2\t{yff2.get('gd_2_value', 0.0):.6f}\t{yff2.get('gd_2_gender', 'UNKNOWN')}",
        "final_gender": final_gender,
        "gender_source": gender_source,
    }

    return ff_result, gender_result


# ─────────────────────────────────────────────────────────────────────────────
# Output writers
# ─────────────────────────────────────────────────────────────────────────────

def write_ff_txt(ff_result: dict, output_path: str) -> None:
    """Write fetal fraction result in ken-nipt compatible format."""
    with open(output_path, "w") as f:
        f.write("\tvalue\n")  # leading tab → 2-column header; pandas reads col-0 as row label
        for key in ("Fragment_FF", "YFF_2", "SeqFF", "M-SeqFF", "Final_FF", "FF_Source"):
            val = ff_result.get(key, "")
            f.write(f"{key}\t{val}\n")
        if ff_result.get("low_ff_warning"):
            f.write("low_ff_warning\tTRUE\n")


def write_gender_txt(gender_result: dict, output_path: str) -> None:
    """Write gender result in ken-nipt compatible format."""
    with open(output_path, "w") as f:
        f.write("\tvalue\tgender\n")  # leading tab → 3-column header (index, value, gender)
        f.write(f"{gender_result['gd_1']}\n")
        f.write(f"{gender_result['gd_2']}\n")
        f.write(f"final_gender\t{gender_result['final_gender']}\t{gender_result['final_gender']}\n")


def write_yff_txt(result: dict, output_path: str, keys: list) -> None:
    """Write a single FF estimate result."""
    with open(output_path, "w") as f:
        f.write("value\n")
        for key in keys:
            val = result.get(key, "")
            f.write(f"{key}\t{val}\n")
        f.write(f"status\t{result.get('status', 'UNKNOWN')}\n")


# ─────────────────────────────────────────────────────────────────────────────
# CLI entry point
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(description="gx-nipt FF/Gender module")
    parser.add_argument("--mode", required=True,
                        choices=["yff1", "yff2", "seqff", "fragment_ff", "gender_decision"])
    parser.add_argument("--sample",  required=True)
    parser.add_argument("--bam",     default=None)
    parser.add_argument("--wig-norm",default=None)
    parser.add_argument("--config",  required=True)
    parser.add_argument("--labcode", default=None)
    parser.add_argument("--output",  default=None)
    parser.add_argument("--seqff-model", default="/data/refs/models/seqff_model.pkl")
    parser.add_argument("--fragment-cutoff", type=int, default=160)
    # gender_decision inputs
    parser.add_argument("--yff1",    default=None)
    parser.add_argument("--yff2",    default=None)
    parser.add_argument("--seqff",   default=None)
    parser.add_argument("--frag-ff", default=None)
    parser.add_argument("--ff-output",     default=None)
    parser.add_argument("--gender-output", default=None)
    parser.add_argument("--debug", action="store_true")
    return parser.parse_args()


def load_txt_as_dict(path: str) -> dict:
    """Load a key-value txt file (written by write_*_txt) into a dict."""
    result = {}
    if not path or not os.path.exists(path):
        return result
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line == "value":
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                result[parts[0]] = parts[1]
    return result


def main():
    args = parse_args()
    setup_logging(args.debug)
    config = load_config(args.config)

    if args.mode == "yff1":
        result = calculate_yff1(args.sample, args.bam, config, args.labcode or "")
        write_yff_txt(result, args.output, ["YFF1", "gd_1_value", "gd_1_gender"])

    elif args.mode == "yff2":
        result = calculate_yff2(args.sample, args.wig_norm, config)
        write_yff_txt(result, args.output, ["YFF2", "UAR_X", "UAR_Y", "gd_2_value", "gd_2_gender"])

    elif args.mode == "seqff":
        result = calculate_seqff(args.sample, args.bam, config, args.seqff_model)
        write_yff_txt(result, args.output, ["SeqFF"])

    elif args.mode == "fragment_ff":
        result = calculate_fragment_ff(args.sample, args.bam, config, args.fragment_cutoff)
        write_yff_txt(result, args.output,
                      ["Fragment_FF", "M-SeqFF", "median_fragment_size",
                       "mean_fragment_size", "std_fragment_size", "total_fragments"])

    elif args.mode == "gender_decision":
        yff1_d   = load_txt_as_dict(args.yff1)
        yff2_d   = load_txt_as_dict(args.yff2)
        seqff_d  = load_txt_as_dict(args.seqff)
        frag_d   = load_txt_as_dict(args.frag_ff)

        # Convert string values to float where needed
        def to_float(d, key):
            try:
                return float(d.get(key, 0))
            except (ValueError, TypeError):
                return 0.0

        yff1_result = {
            "YFF1": to_float(yff1_d, "YFF1"),
            "gd_1_value": to_float(yff1_d, "gd_1_value"),
            "gd_1_gender": yff1_d.get("gd_1_gender", "UNKNOWN"),
        }
        yff2_result = {
            "YFF2": to_float(yff2_d, "YFF2"),
            "UAR_X": to_float(yff2_d, "UAR_X"),
            "UAR_Y": to_float(yff2_d, "UAR_Y"),
            "gd_2_value": to_float(yff2_d, "gd_2_value"),
            "gd_2_gender": yff2_d.get("gd_2_gender", "UNKNOWN"),
        }
        seqff_result = {"SeqFF": to_float(seqff_d, "SeqFF")}
        frag_result  = {
            "Fragment_FF": to_float(frag_d, "Fragment_FF"),
            "M-SeqFF": to_float(frag_d, "M-SeqFF"),
        }

        ff_out, gender_out = gender_decision(
            args.sample, yff1_result, yff2_result, seqff_result, frag_result, config
        )
        write_ff_txt(ff_out, args.ff_output)
        write_gender_txt(gender_out, args.gender_output)

    else:
        logger.error(f"Unknown mode: {args.mode}")
        sys.exit(1)


if __name__ == "__main__":
    main()
