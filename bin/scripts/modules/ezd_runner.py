import logging
import os
from os.path import join

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from sca_detector import SCADetector

matplotlib.use("Agg")

sca_detector = None
logger = logging.getLogger(__name__)

# Constants
DATA_DIR = os.environ.get("DATA_DIR", "/Work/NIPT/data")
ANALYSIS_DIR = os.environ.get("ANALYSIS_DIR", "/Work/NIPT/analysis")

# -----------------------------
# Core EZD Functions
# -----------------------------


def analyze_chromosomes(normalization_file, sample_name, output_dir):
    logger.info("1. Loading and preprocessing data...")

    normal = pd.read_csv(normalization_file, sep="\t")

    # cor.gc NA row deletion
    filtered_normal = normal.dropna(subset=["cor.gc"])
    all_filter = filtered_normal[["chr", "start", "end", "cor.gc"]].copy()
    all_data = all_filter[all_filter["cor.gc"] != 0].copy()

    chr_data = {}
    for i in range(1, 23):
        chr_name = f"chr{i}"
        chr_data[chr_name] = all_data[all_data["chr"] == chr_name].copy()

    chrX = all_data[all_data["chr"] == "chrX"].copy()
    chrY = all_data[all_data["chr"] == "chrY"].copy()

    # centromere exclusion
    chrX = chrX[(chrX["start"] < 58100001) | (chrX["start"] > 63000000)]
    chrY = chrY[(chrY["start"] < 11600001) | (chrY["start"] > 14000000)]

    logger.info("2. Calculating median values for each chromosome...")

    all_med = []
    for i in range(1, 23):
        chr_name = f"chr{i}"
        if len(chr_data[chr_name]) > 0:
            med_val = chr_data[chr_name]["cor.gc"].median()
        else:
            med_val = 0
        all_med.append(med_val)

    # 3. First Z-score calculation
    logger.info("3. Calculating first Z-scores...")

    median_chrX = chrX["cor.gc"].median() if len(chrX) > 0 else 0
    all_med_array = np.array(all_med)
    median_all_med_sorted = np.sort(all_med_array)[10]
    sd_all_med = np.std(all_med_array, ddof=1)

    first_zscore = (all_med_array - median_all_med_sorted) / sd_all_med

    # 4. Reference set for 2nd zscore (Excluding chr13, chr17-22)
    logger.info("4. Creating reference set...")

    # chr13(index 12), chr17-22(index 16-21) are excluded
    # chr1 ~ 12 + chr 14, chr15, chr16
    ref_indices = list(range(12)) + list(range(13, 16))
    # chr1 ~ 12 + chr 14, chr15
    # ref_indices = list(range(12)) + list(range(13, 15))
    ref_set = first_zscore[ref_indices]

    # 5. Second Z-score 계산
    logger.info("5. Calculating second Z-scores...")

    median_ref_set = np.sort(ref_set)[7]
    sd_ref_set = np.std(ref_set, ddof=1)

    second_zscore = (first_zscore - median_ref_set) / sd_ref_set

    first_zscore_X = (median_chrX - median_all_med_sorted) / sd_all_med
    second_zscore_X = (first_zscore_X - median_ref_set) / sd_ref_set

    zscore_output = normalization_file.replace(".txt", ".zscore.txt")
    np.savetxt(zscore_output, second_zscore, fmt="%.6f", delimiter="\t")

    # 6. UAR% (Unique Alignee Read Percent) calculation
    logger.info("5. Calculating UAR percentages...")

    each_chr_sum = []
    for i in range(1, 23):
        chr_name = f"chr{i}"
        if len(chr_data[chr_name]) > 0:
            chr_sum = chr_data[chr_name]["cor.gc"].sum()
        else:
            chr_sum = 0
        each_chr_sum.append(chr_sum)

    autosome_sum = sum(each_chr_sum)

    UAR = np.array(each_chr_sum) / autosome_sum * 100
    UAR_X = chrX["cor.gc"].sum() / autosome_sum * 100 if len(chrX) > 0 else 0
    UAR_Y = chrY["cor.gc"].sum() / autosome_sum * 100 if len(chrY) > 0 else 0

    results = {
        "sample_name": sample_name,
        "median_values": all_med,
        "first_zscore": first_zscore.tolist(),
        "second_zscore": second_zscore.tolist(),
        "first_zscore_X": first_zscore_X,
        "second_zscore_X": second_zscore_X,
        "UAR": UAR.tolist(),
        "UAR_X": UAR_X,
        "UAR_Y": UAR_Y,
        "autosome_sum": autosome_sum,
    }

    # 결과를 TSV 파일로 저장
    save_results_to_tsv(results, output_dir)

    # UAR, Z-score만 포함하는 간단한 테이블 저장
    save_simple_ur_zscore_table(results, output_dir)

    # 결과 출력
    logger.info(f"\nAnalysis Results for {sample_name}:")
    logger.info(f"UAR_X: {UAR_X:.2f}%")
    logger.info(f"UAR_Y: {UAR_Y:.2f}%")

    return results


def run_ezd_pipeline(sample_name, wig_path, labcode, output_dir, group):
    """
    Parameters:
    - sample_name: 샘플 이름
    - wig_path: wig 파일 경로 (normalization.txt로 변환 필요)
    - labcode: 실험실 코드
    - output_dir: 출력 디렉토리
    - group: 그룹명 (orig, fetus, mom)

    Returns:
    - pd.DataFrame: UAR과 Z-score를 포함한 데이터프레임
    """

    # analyze_chromosomes 함수 호출
    results = analyze_chromosomes(wig_path, sample_name, output_dir)

    data_rows = []

    # Autosomal (chr1-chr22)
    for i in range(22):
        chr_name = f"chr{i + 1}"
        ur_val = results["UAR"][i]
        z_val = results["second_zscore"][i]
        data_rows.append({"Chromosome": chr_name, "UAR": ur_val, "Z": z_val})

    # X, Y added
    data_rows.append(
        {"Chromosome": "chrX", "UAR": results["UAR_X"], "Z": results["second_zscore_X"]}
    )

    data_rows.append(
        {
            "Chromosome": "chrY",
            "UAR": results["UAR_Y"],
            "Z": 0.0,  # Y 염색체는 Z-score 없음
        }
    )

    ezd_df = pd.DataFrame(data_rows)

    ezd_output_file = os.path.join(output_dir, f"{group}_ezd_results.tsv")
    ezd_df.to_csv(ezd_output_file, sep="\t", index=False)

    logger.info(f"EZD pipeline completed for {sample_name} - {group}")

    return ezd_df


def run_decision(ezd_df, threshold_path, config_type="orig"):
    """
    임계값을 기반으로 각 염색체의 EZD 결과를 판정하는 함수

    Parameters:
    - ezd_df: run_ezd_pipeline의 결과 DataFrame
    - threshold_path: 임계값 파일 경로
    - config_type: SCA 설정 타입 ('orig', 'fetus', 'mom')

    Returns:
    - pd.DataFrame: 각 염색체별 판정 결과 (chr, result, UAR, Z, config)
    """
    threshold_df = pd.read_csv(threshold_path, sep="\t")
    logger.info(threshold_df)
    decisions = []

    for i in range(22):
        chr_name = f"chr{i + 1}"
        row = threshold_df[threshold_df["chr"] == chr_name]

        # ezd_df에서 해당 염색체 데이터 찾기
        chr_row = ezd_df[ezd_df["Chromosome"] == chr_name]

        if chr_row.empty or row.empty:
            decision = "Not Detected"
            uar_value = np.nan
            z_value = np.nan
        else:
            uar_value = chr_row["UAR"].values[0]
            z_value = chr_row["Z"].values[0]

            ur_min = row["UAR_min"].values[0]
            ur_max = row["UAR_max"].values[0]
            z_min = row["Z_min"].values[0]
            z_max = row["Z_max"].values[0]

            logger.info(f"--------- {chr_name} ----------")
            logger.info(f"ur_min:{ur_min}, z_min:{z_min}")
            logger.info(f"ur_max:{ur_max}, z_max:{z_max}")
            logger.info(f"uar_value:{uar_value}, z_value:{z_value}")
            if z_value >= z_max and uar_value >= ur_max:
                decision = "Detected"
            # 250615 : AND --> OR
            elif z_min <= z_value or ur_min <= uar_value:
                decision = "Suspected"
            else:
                decision = "Not Detected"

        # 5개 컬럼 (config 추가)
        decisions.append(
            {
                "chr": chr_name,
                "result": decision,
                "UAR": uar_value,
                "Z": z_value,
                # ,
                #'config': config_type
            }
        )

    # DataFrame 생성
    decision_df = pd.DataFrame(decisions)

    logger.info(f"Trisomy decision 완료: {decision_df.shape}")

    # 결과 요약
    detected_count = (decision_df["result"] == "Detected").sum()
    suspected_count = (decision_df["result"] == "Suspected").sum()

    logger.info(
        f"Trisomy 결과 ({config_type}): Detected={detected_count}, Suspected={suspected_count}"
    )

    return decision_df


def add_sca_to_decision_df(decision_df, sca_results, config_type="orig"):
    """
    간단한 구조의 decision_df에 SCA 결과 추가

    Args:
        decision_df (DataFrame): run_decision에서 생성된 DataFrame (chr, result, UAR, Z, config)
        sca_results (dict): SCA detection 결과
        config_type (str): 사용한 설정 타입

    Returns:
        DataFrame: SCA 결과가 추가된 decision_df
    """
    try:
        # 기존 decision_df 복사
        result_df = decision_df.copy()

        # SCA 결과 행들 생성
        sca_rows = []

        # Female SCA (chrX) 결과 추가
        if "female" in sca_results:
            female_data = sca_results["female"]

            sca_rows.append(
                {
                    "chr": "chrX",
                    "result": female_data["detection"],
                    "UAR": female_data["ur_x"],
                    "Z": female_data["z_score"],
                    # ,
                    #'config': config_type
                }
            )
            logger.info(
                f"Female SCA 추가: chrX, {female_data['detection']} ({config_type})"
            )

        # Male SCA (chrY) 결과 추가
        if "male" in sca_results:
            male_data = sca_results["male"]

            sca_rows.append(
                {
                    "chr": "chrY",
                    "result": male_data["detection"],
                    "UAR": male_data["ur_y"],
                    "Z": 0.0,
                    #'Z': np.nan
                    # ,  # chrY는 Z-score 없음
                    #'config': config_type
                }
            )
            logger.info(
                f"Male SCA 추가: chrY, {male_data['detection']} ({config_type})"
            )

        if sca_rows:
            # SCA 결과를 DataFrame으로 변환
            sca_df = pd.DataFrame(sca_rows)

            # 기존 DataFrame과 합치기
            final_df = pd.concat([result_df, sca_df], ignore_index=True)

            logger.info(f"SCA 결과 추가 완료: {len(sca_rows)}개 행")
            logger.info(f"최종 DataFrame 크기: {final_df.shape}")

            return final_df
        else:
            logger.warning("추가할 SCA 결과가 없습니다.")
            return result_df

    except Exception as e:
        logger.error(f"SCA 결과 추가 오류: {e}")
        return decision_df


def init_sca_detector(data_dir, labcode):
    """SCA Detector 초기화"""
    global sca_detector
    if sca_detector is None:
        try:
            sca_detector = SCADetector(data_dir, labcode)
            sca_detector.load_all_configs()
            logger.info(f"SCA Detector 초기화 완료: {labcode}")
        except Exception as e:
            logger.error(f"SCA Detector 초기화 실패: {e}")
            sca_detector = None
    return sca_detector


def run_sca_detection(ezd_df, config, config_type="orig"):
    if sca_detector is None:
        logger.error("SCA Detector has not been initialized")
        return {}

    sca_results = {}

    try:
        # chrX, chrY 데이터 추출
        chrx_row = ezd_df[ezd_df["Chromosome"] == "chrX"]
        chry_row = ezd_df[ezd_df["Chromosome"] == "chrY"]

        if len(chrx_row) == 0:
            logger.warning("chrX data is not found")
            return {}

        # chrX 데이터 추출
        test_ur_x = chrx_row["UAR"].iloc[0]
        test_z_score = chrx_row["Z"].iloc[0]

        # chrY 데이터 추출 (있는 경우에만)
        test_ur_y = None
        if len(chry_row) > 0:
            test_ur_y = chry_row["UAR"].iloc[0]
            # chrY의 Z 값은 보통 NaN이므로 사용하지 않음

        logger.info("SCA Detection data extraction:")
        logger.info(f"  chrX: UAR={test_ur_x:.6f}, Z={test_z_score:.6f}")
        if test_ur_y is not None:
            logger.info(f"  chrY: UAR={test_ur_y:.6f}")

        # Female SCA Detection (chrX 기반)
        if test_ur_x is not None and test_z_score is not None:
            female_results = sca_detector.compare_configs(
                test_ur_x, z_score=test_z_score, gender="female"
            )
            female_detection = female_results.get(config_type, {}).get(
                "result", "Not Detected"
            )

            sca_results["female"] = {
                "results": female_results,
                "ur_x": test_ur_x,
                "z_score": test_z_score,
                "detection": female_detection,
            }
            logger.info(f"Female SCA ({config_type}): {female_detection}")

        # Male SCA Detection (chrX + chrY 기반)
        if test_ur_x is not None and test_ur_y is not None:
            male_results = sca_detector.compare_configs(
                test_ur_x, ur_y=test_ur_y, gender="male"
            )
            male_detection = male_results.get(config_type, {}).get(
                "result", "Not Detected"
            )

            sca_results["male"] = {
                "results": male_results,
                "ur_x": test_ur_x,
                "ur_y": test_ur_y,
                "detection": male_detection,
            }
            logger.info(f"Male SCA ({config_type}): {male_detection}")

            # Male 결과 상세 출력
            logger.info("Male SCA Detection Result:")
            sca_detector.print_male_results(male_results, test_ur_x, test_ur_y)
        else:
            logger.warning("No chrY data. Male SCA detection cannot be performed")

        return sca_results

    except Exception as e:
        logger.error(f"SCA detection 오류: {e}")
        logger.error(f"ezd_df 구조: {ezd_df.columns.tolist()}")
        logger.error(f"ezd_df 크기: {ezd_df.shape}")
        return {}


def run_ezd_group(
    sample_name: str,
    group: str,
    wig_path: str,
    labcode: str,
    analysis_dir: str,
    data_dir: str,
    config,
):
    """
    Parameters:
    - sample_name: 샘플 이름
    - group: 그룹명 (orig, fetus, mom)
    - wig_path: wig 파일 경로
    - labcode: 실험실 코드
    - analysis_dir: 분석 디렉토리
    - data_dir: 데이터 디렉토리
    """

    init_sca_detector(data_dir, labcode)

    logger.info("\n=== Male configuration ===")
    orig_male = sca_detector.get_male_params("orig")
    fetus_male = sca_detector.get_male_params("fetus")
    mom_male = sca_detector.get_male_params("mom")  # None이어야 함

    logger.info(f"Orig Male: {orig_male}")
    logger.info(f"Fetus Male: {fetus_male}")
    logger.info(f"Mom Male: {mom_male}")

    logger.info("\n=== Female configuration ===")
    orig_female = sca_detector.get_female_params("orig")
    fetus_female = sca_detector.get_female_params("fetus")
    mom_female = sca_detector.get_female_params("mom")

    logger.info(f"Orig Female XO threshold: {orig_female['xo_z_threshold']}")
    logger.info(f"Fetus Female XO threshold: {fetus_female['xo_z_threshold']}")
    logger.info(f"Mom Female XO threshold: {mom_female['xo_z_threshold']}")

    # Output directory
    output_dir = join(analysis_dir, sample_name, "Output_EZD", group)
    os.makedirs(output_dir, exist_ok=True)

    threshold_path = join(
        data_dir, "refs", labcode, "EZD", group, f"{group}_thresholds_new.tsv"
    )

    # For plotting
    chr_table_dir = join(data_dir, "refs", labcode, "EZD", f"{group}")

    logger.info(f"Starting EZD pipeline for {sample_name} - {group}")

    try:
        # 1. EZD main pipeline (analyze_chromosomes)
        ezd_df = run_ezd_pipeline(sample_name, wig_path, labcode, output_dir, group)

        # 2. Decision
        decision_df = run_decision(ezd_df, threshold_path)

        # 2.1 SCA Decision
        sca_results = run_sca_detection(ezd_df, config, group)
        decision_df_extended = add_sca_to_decision_df(decision_df, sca_results, group)

        # 3. Save results
        decision_output = join(
            output_dir, f"Trisomy_detect_result_{group}_with_SCA.tsv"
        )

        decision_df_rounded = decision_df_extended.round({"UAR": 3, "Z": 3})
        decision_df_rounded.to_csv(decision_output, sep="\t", index=False)

        # UAR, Z를 float으로 변환 후 반올림
        # df_round['UAR'] = pd.to_numeric(decision_df_extended['UAR'], errors='coerce').round(3)
        # df_round['Z']   = pd.to_numeric(decision_df_extended['Z'],   errors='coerce').round(3)

        output = join(output_dir, f"Trisomy_detect_result_{group}_with_SCA.tsv")
        decision_df_rounded.to_csv(decision_output, sep="\t", index=False)

        dpi = config.get("EZD", {}).get("resolution_dpi", 200)

        # 4. Plotting
        try:
            plot_chr_scatter_grid(
                chr_table_dir,
                sample_name,
                threshold_path,
                join(output_dir, f"{group}_EZD_grid.png"),
                group,
                output_dir,
                dpi,
            )
            # plot_ezd_interactive(ezd_df, join(output_dir, f'{group}_EZD_plot.html'))
        except Exception as plot_error:
            logger.warning(f"Plotting failed: {plot_error}")

        logger.info(f"EZD group analysis completed for {sample_name} - {group}")

        return ezd_df, decision_df

    except Exception as e:
        logger.error(f"EZD group analysis failed for {sample_name} - {group}: {e}")
        raise


def save_simple_ur_zscore_table(results, output_dir):
    """
    UAR과 Z-score만 포함하는 간단한 테이블을 TSV로 저장

    Parameters:
    - results: analyze_chromosomes 함수의 결과 딕셔너리
    - output_dir: 출력 디렉토리
    """

    # 출력 파일 경로 설정
    output_file = os.path.join(output_dir, "ur_zscore_table.tsv")
    os.makedirs(output_dir, exist_ok=True)

    # TSV 파일 작성
    with open(output_file, "w") as f:
        # 헤더 작성
        f.write("Chromosome\tUAR\tZ\n")

        # 상염색체 데이터 (chr1-chr22)
        for i in range(22):
            chr_name = f"chr{i + 1}"
            ur_val = results["UAR"][i]
            z_val = results["second_zscore"][i]

            f.write(f"{chr_name}\t{ur_val:.4f}\t{z_val:.6f}\n")

        # X 염색체 데이터
        f.write(f"chrX\t{results['UAR_X']:.4f}\t{results['second_zscore_X']:.6f}\n")

        # Y 염색체 데이터 (Z-score는 0.0)
        f.write(f"chrY\t{results['UAR_Y']:.4f}\t0.0\n")

    logger.info(f"Simple UAR/Z-score table saved to: {output_file}")


def save_results_to_tsv(results, output_dir):
    # 출력 파일 경로 설정
    output_file = os.path.join(output_dir, "ezd_result.tsv")
    os.makedirs(output_dir, exist_ok=True)

    # TSV 파일 작성
    with open(output_file, "w") as f:
        # 헤더 작성
        f.write(
            "Category\tChromosome\tMedian_Value\tFirst_Zscore\tSecond_Zscore\tUAR_Percent\n"
        )

        # 상염색체 데이터 (chr1-chr22)
        for i in range(22):
            chr_name = f"chr{i + 1}"
            median_val = results["median_values"][i]
            first_z = results["first_zscore"][i]
            second_z = results["second_zscore"][i]
            ur_pct = results["UAR"][i]

            f.write(
                f"Autosome\t{chr_name}\t{median_val:.6f}\t{first_z:.6f}\t{second_z:.6f}\t{ur_pct:.4f}\n"
            )

        # X 염색체 데이터
        f.write(
            f"Sex_Chromosome\tchrX\t-\t{results['first_zscore_X']:.6f}\t{results['second_zscore_X']:.6f}\t{results['UAR_X']:.4f}\n"
        )

        # Y 염색체 데이터
        f.write(f"Sex_Chromosome\tchrY\t-\t-\t-\t{results['UAR_Y']:.4f}\n")

        # 요약 정보
        f.write("\n# Summary Information\n")
        f.write(f"Sample_Name\t{results['sample_name']}\n")

    logger.info(f"Results saved to: {output_file}")


# -----------------------------
# Plotting
# -----------------------------


def add_sca_lines_to_chrx_plot(ax, config_type="orig"):
    """
    chrX (Female) subplot에 SCA detection 라인 추가

    Args:
        ax: matplotlib axes 객체
        config_type (str): 사용할 설정 타입 ('orig', 'fetus', 'mom')
    """
    if sca_detector is None:
        # SCA Manager가 없으면 기본값 사용
        ax.axhline(
            y=-6, color="gray", linestyle="--", alpha=0.8, linewidth=1
        )  # XO threshold
        ax.axhline(
            y=4.5, color="gray", linestyle="--", alpha=0.8, linewidth=1
        )  # XXX threshold
        ax.axvline(x=5.2, color="gray", linestyle="--", alpha=0.8, linewidth=1)
        ax.axvline(x=5.6, color="gray", linestyle="--", alpha=0.8, linewidth=1)
        return

    try:
        # Female 파라미터 가져오기
        female_params = sca_detector.get_female_params(config_type)
        if female_params is None:
            return

        logger.info("---------- Threshold ------------------")
        logger.info(f"xo_z_threshold : {female_params['xo_z_threshold']}")
        logger.info(f"xxx_z_threshold : {female_params['xxx_z_threshold']}")
        logger.info(f"ur_x_low : {female_params['ur_x_low']}")
        logger.info(f"ur_x_high : {female_params['ur_x_high']}")
        logger.info(f"z_normal_low : {female_params['z_normal_low']}")
        logger.info(f"z_normal_high : {female_params['z_normal_high']}")
        logger.info(f"xo_ur_x_min : {female_params['xo_ur_x_min']}")
        logger.info(f"xo_ur_x_max : {female_params['xo_ur_x_max']}")
        logger.info(f"xxx_ur_x_min : {female_params['xxx_ur_x_min']}")
        logger.info(f"xxx_ur_x_max : {female_params['xxx_ur_x_max']}")
        logger.info("---------------------------------------")

        # XO threshold
        ax.axhline(
            y=female_params["xo_z_threshold"],
            color="gray",
            linestyle="--",
            alpha=0.8,
            linewidth=1,
        )
        ax.axvline(
            x=female_params["xo_ur_x_max"],
            color="gray",
            linestyle="--",
            alpha=0.8,
            linewidth=1,
        )

        # XXX threshold
        ax.axhline(
            y=female_params["xxx_z_threshold"],
            color="gray",
            linestyle="--",
            alpha=0.8,
            linewidth=1,
        )
        ax.axvline(
            x=female_params["xxx_ur_x_min"],
            color="gray",
            linestyle="--",
            alpha=0.8,
            linewidth=1,
        )

        # Normal range
        ax.axhline(
            y=female_params["z_normal_high"],
            color="gray",
            linestyle="-",
            alpha=0.8,
            linewidth=1,
        )
        ax.axhline(
            y=female_params["z_normal_low"],
            color="gray",
            linestyle="-",
            alpha=0.8,
            linewidth=1,
        )
        ax.axvline(
            x=female_params["ur_x_low"],
            color="gray",
            linestyle="-",
            alpha=0.8,
            linewidth=1,
        )
        ax.axvline(
            x=female_params["ur_x_high"],
            color="gray",
            linestyle="-",
            alpha=0.8,
            linewidth=1,
        )

        # plot 설정에서 axis limit 가져오기
        try:
            plot_settings = sca_manager.get_plot_settings(config_type, "female")
            if plot_settings and "gender_settings" in plot_settings:
                female_plot_settings = plot_settings["gender_settings"]
                if "x_axis_range" in female_plot_settings:
                    ax.set_xlim(female_plot_settings["x_axis_range"])
                if "y_axis_range" in female_plot_settings:
                    ax.set_ylim(female_plot_settings["y_axis_range"])
        except:
            # plot 설정이 없으면 기본값
            ax.set_xlim(4.5, 6)
            ax.set_ylim(-30, 30)

    except Exception as e:
        logger.error(f"chrX SCA 라인 추가 오류: {e}")


def add_sca_lines_to_chry_plot(ax, config_type="orig"):
    """
    chrY (Male) subplot에 SCA detection 라인 추가

    Args:
        ax: matplotlib axes 객체
        config_type (str): 사용할 설정 타입 ('orig', 'fetus', 'mom')
    """
    if sca_detector is None:
        # SCA Manager가 없으면 기본값 사용
        ax.axvline(x=4.7, color="gray", linestyle="--", alpha=0.8, linewidth=1)
        ax.axvline(x=5.32, color="gray", linestyle="--", alpha=0.8, linewidth=1)
        ax.axvline(x=5.48, color="gray", linestyle="--", alpha=0.8, linewidth=1)
        ax.axhline(y=0.035, color="gray", linestyle="--", alpha=0.8, linewidth=1)
        # 기본 대각선
        slope = -0.08367685038621694
        intercept = 0.4876543737788628
        x_line = np.linspace(4.0, 6.0, 100)
        y_boundary = slope * x_line + intercept
        ax.plot(x_line, y_boundary, color="purple", linestyle="-", linewidth=2)
        return

    try:
        # Male 파라미터 가져오기
        male_params = sca_detector.get_male_params(config_type)
        if male_params is None:
            return

        # 기존 임계값 라인들 (회색, 얇게)
        # ax.axvline(x=4.7, color='gray', linestyle='--', alpha=0.5, linewidth=1)
        # ax.axvline(x=5.32, color='gray', linestyle='--', alpha=0.5, linewidth=1)
        # ax.axvline(x=5.48, color='gray', linestyle='--', alpha=0.5, linewidth=1)
        # ax.axhline(y=0.035, color='gray', linestyle='--', alpha=0.5, linewidth=1)

        # XYY/XXY 구분 세로선 (config 기반)
        ax.axvline(
            x=male_params["ur_x_threshold"],
            color="gray",
            linestyle="-",
            alpha=0.8,
            linewidth=1,
        )

        # SCA 경계선 (대각선, config 기반)
        x_line = np.linspace(4.0, 6.0, 100)
        y_boundary = male_params["slope"] * x_line + male_params["intercept"]
        ax.plot(x_line, y_boundary, color="red", linestyle="--", linewidth=1)

        try:
            plot_settings = sca_manager.get_plot_settings(config_type, "male")
            if plot_settings and "gender_settings" in plot_settings:
                male_plot_settings = plot_settings["gender_settings"]
                if "x_axis_range" in male_plot_settings:
                    ax.set_xlim(male_plot_settings["x_axis_range"])
                if "y_axis_range" in male_plot_settings:
                    ax.set_ylim(male_plot_settings["y_axis_range"])
        except:
            # plot 설정이 없으면 기본값
            ax.set_xlim(4, 6)

    except Exception as e:
        logger.error(f"chrY SCA 라인 추가 오류: {e}")


def plot_chr_scatter_grid(
    chr_table_dir,
    sample_name,
    threshold_path,
    output_path,
    group,
    test_data_dir=None,
    dpi=200,
):
    """
    염색체별 UAR vs Z-score scatter plot grid 생성

    Parameters:
    - chr_table_dir: chr1.txt~chr22.txt 파일들이 있는 디렉토리
    - sample_name: 테스트 샘플 이름
    - threshold_path: 임계값 파일 경로 (orig_thresholds.tsv 등)
    - output_path: 출력 이미지 경로
    - test_data_dir: 테스트 샘플 데이터 디렉토리 (예: Output_EZD/orig)
    """

    try:
        logger.info("Creating scatter grid plot...")
        logger.info(f"  Chr data dir: {chr_table_dir}")
        logger.info(f"  Sample name: {sample_name}")
        logger.info(f"  Threshold file: {threshold_path}")
        logger.info(f"  Test data dir: {test_data_dir}")
        logger.info(f"  Output: {output_path}")

        # 테스트 샘플 데이터 로드
        test_sample_data = {}
        if test_data_dir and os.path.exists(test_data_dir):
            # 그룹명 추출 (orig, fetus, mom)
            group_name = os.path.basename(threshold_path).replace(
                "_thresholds_new.tsv", ""
            )
            test_file = os.path.join(test_data_dir, f"{group_name}_ezd_results.tsv")

            logger.info(f"  Loading test data from: {test_file}")

            if os.path.exists(test_file):
                try:
                    test_df = pd.read_csv(test_file, sep="\t")
                    # logger.info(f"    Test data loaded: {len(test_df)} chromosomes")
                    # logger.info(f"    Test data columns: {list(test_df.columns)}")

                    # 염색체별로 UAR, Z 값 저장
                    for _, row in test_df.iterrows():
                        chr_name = row["Chromosome"]
                        # if chr_name.startswith('chr') and chr_name[3:].isdigit():  # chr1, chr2, ... chr22만
                        test_sample_data[chr_name] = {"UAR": row["UAR"], "Z": row["Z"]}

                    # logger.info(f"    Test data for chromosomes: {list(test_sample_data.keys())}")

                except Exception as e:
                    logger.info(f"    Error loading test data: {e}")
            else:
                logger.info(f"    Test data file not found: {test_file}")
        else:
            logger.info("    No test data directory provided or not found")

        # 디렉토리 존재 확인
        if not os.path.exists(chr_table_dir):
            logger.info(f"  ERROR: Directory does not exist: {chr_table_dir}")
            return

        # 디렉토리 내 파일 목록 확인
        files_in_dir = os.listdir(chr_table_dir)
        chr_files = [
            f for f in files_in_dir if f.startswith("chr") and f.endswith(".txt")
        ]
        # logger.info(f"  Found chr files: {chr_files}")

        if not chr_files:
            logger.info(f"  ERROR: No chr*.txt files found in {chr_table_dir}")
            logger.info(f"  Files in directory: {files_in_dir}")
            return

        # 임계값 데이터 로드
        if os.path.exists(threshold_path):
            threshold_df = pd.read_csv(threshold_path, sep="\t")
            # logger.info(f"  Loaded thresholds for {len(threshold_df)} chromosomes")
            # logger.info(f"  Threshold columns: {list(threshold_df.columns)}")
        else:
            logger.info(f"  Warning: Threshold file not found: {threshold_path}")
            threshold_df = pd.DataFrame()

        # 4x6 그리드로 subplot 생성 (22개 염색체)
        fig, axes = plt.subplots(4, 6, figsize=(24, 16))
        plt.subplots_adjust(wspace=0.3, hspace=0.4)

        # 범례 표시 여부 체크
        legend_added = False
        total_plots_with_data = 0

        # ----------------------------------------------------
        # 각 염색체에 대해 플롯 생성
        # ----------------------------------------------------
        for i in range(22):
            row, col = divmod(i, 6)
            ax = axes[row][col]
            chr_name = f"chr{i + 1}"

            # 염색체 파일 경로
            chr_file = os.path.join(chr_table_dir, f"{chr_name}.txt")

            # logger.info(f"    Checking {chr_file}")

            # 파일이 없는 경우
            if not os.path.exists(chr_file):
                logger.info(f"      File not found: {chr_file}")
                ax.text(
                    0.5,
                    0.5,
                    f"{chr_name}\nNo Data",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                    fontsize=10,
                    color="gray",
                )
                ax.set_title(chr_name, fontsize=12, fontweight="bold")
                ax.set_xlabel("UAR (%)", fontsize=10)
                ax.set_ylabel("Z-score", fontsize=10)
                ax.grid(True, alpha=0.3)
                continue

            try:
                # 염색체 데이터 로드
                # logger.info(f"      Loading {chr_file}")
                df = pd.read_csv(chr_file, sep="\t")
                # logger.info(f"      Loaded {len(df)} rows")
                # logger.info(f"      Columns: {list(df.columns)}")

                # 처음 몇 행 출력
                # if len(df) > 0:
                #    logger.info(f"      First few rows:")
                #    logger.info(df.head(3).to_string())

                # 데이터가 비어있는지 확인
                if df.empty:
                    #    logger.info(f"      Empty dataframe")
                    ax.text(
                        0.5,
                        0.5,
                        f"{chr_name}\nEmpty",
                        ha="center",
                        va="center",
                        transform=ax.transAxes,
                        fontsize=10,
                        color="gray",
                    )
                    ax.set_title(chr_name, fontsize=12, fontweight="bold")
                    ax.set_xlabel("UAR (%)", fontsize=10)
                    ax.set_ylabel("Z-score", fontsize=10)
                    ax.grid(True, alpha=0.3)
                    continue

                # 필수 컬럼 확인
                required_columns = ["sample", "type", "UAR", "Z"]
                missing_columns = [
                    col for col in required_columns if col not in df.columns
                ]
                if missing_columns:
                    logger.info(f"      Missing columns: {missing_columns}")
                    ax.text(
                        0.5,
                        0.5,
                        f"{chr_name}\nMissing columns:\n{missing_columns}",
                        ha="center",
                        va="center",
                        transform=ax.transAxes,
                        fontsize=8,
                        color="red",
                    )
                    ax.set_title(chr_name, fontsize=12, fontweight="bold")
                    ax.set_xlabel("UAR (%)", fontsize=10)
                    ax.set_ylabel("Z-score", fontsize=10)
                    ax.grid(True, alpha=0.3)
                    continue

                total_plots_with_data += 1
                # logger.info(f"      Processing {chr_name}: {len(df)} samples")

                # 각 타입별 샘플 수 확인
                type_counts = df["type"].value_counts()
                # logger.info(f"      Type counts: {dict(type_counts)}")

                # N 샘플 (Normal - 검정색)
                df_N = df[df["type"] == "N"]
                if not df_N.empty:
                    ax.scatter(
                        df_N["UAR"],
                        df_N["Z"],
                        color="black",
                        s=8,
                        alpha=0.6,
                        facecolors="none",
                        edgecolors="black",
                        linewidths=0.5,
                        label="Not detected" if not legend_added else "",
                    )
                #    logger.info(f"        Plotted {len(df_N)} Normal samples")

                # P 샘플 (Positive - 녹색)
                df_P = df[df["type"] == "P"]
                if not df_P.empty:
                    ax.scatter(
                        df_P["UAR"],
                        df_P["Z"],
                        color="green",
                        s=8,
                        alpha=0.6,
                        facecolors="none",
                        edgecolors="green",
                        linewidths=0.5,
                        label="Positive" if not legend_added else "",
                    )
                    # logger.info(f"        Plotted {len(df_P)} Positive samples")

                # 테스트 샘플 그리기 (빨간색) - 별도 파일에서 로드
                test_ur, test_z = None, None
                if chr_name in test_sample_data:
                    test_ur = test_sample_data[chr_name]["UAR"]
                    test_z = test_sample_data[chr_name]["Z"]

                    # 유효한 숫자 값인지 확인
                    try:
                        test_ur_num = float(test_ur)
                        test_z_num = float(test_z)

                        ax.scatter(
                            test_ur_num,
                            test_z_num,
                            color="red",
                            edgecolor="white",
                            s=50,
                            linewidths=2,
                            label="Test Sample" if not legend_added else "",
                            zorder=10,
                        )

                        # logger.info(f"        Test sample plotted: UAR={test_ur_num:.3f}, Z={test_z_num:.3f}")

                    except (ValueError, TypeError):
                        logger.info(
                            f"        Test sample data invalid: UAR={test_ur}, Z={test_z}"
                        )
                else:
                    logger.info(f"        No test sample data for {chr_name}")

                # 임계값 라인 그리기
                if not threshold_df.empty:
                    th_row = threshold_df[threshold_df["chr"] == chr_name]

                    if not th_row.empty:
                        th_data = th_row.iloc[0]
                        # logger.info(f"        Drawing thresholds: UAR({th_data.get('UAR_min', 'N/A')}-{th_data.get('UAR_max', 'N/A')}), Z({th_data.get('Z_min', 'N/A')}-{th_data.get('Z_max', 'N/A')})")

                        # UAR 임계값 (수직선)
                        if "UAR_min" in th_data and pd.notna(th_data["UAR_min"]):
                            ax.axvline(
                                th_data["UAR_min"],
                                color="gray",
                                linestyle="--",
                                alpha=0.8,
                                linewidth=1,
                            )
                        if "UAR_max" in th_data and pd.notna(th_data["UAR_max"]):
                            ax.axvline(
                                th_data["UAR_max"],
                                color="gray",
                                linestyle="--",
                                alpha=0.8,
                                linewidth=1,
                            )

                        # Z 임계값 (수평선)
                        if "Z_min" in th_data and pd.notna(th_data["Z_min"]):
                            ax.axhline(
                                th_data["Z_min"],
                                color="gray",
                                linestyle="--",
                                alpha=0.8,
                                linewidth=1,
                            )
                        if "Z_max" in th_data and pd.notna(th_data["Z_max"]):
                            ax.axhline(
                                th_data["Z_max"],
                                color="gray",
                                linestyle="--",
                                alpha=0.8,
                                linewidth=1,
                            )

                # 축 레이블 및 제목 설정
                ax.set_xlabel("UAR (%)", fontsize=10)
                ax.set_ylabel("Z-score", fontsize=10)
                ax.set_title(chr_name, fontsize=12, fontweight="bold")

                # 범례 추가 (첫 번째 subplot에만)
                if not legend_added and (
                    not df_N.empty or not df_P.empty or test_ur is not None
                ):
                    ax.legend(loc="lower right", fontsize=8)
                    legend_added = True

                # 그리드 추가
                ax.grid(True, alpha=0.3)

                # 축 범위 자동 조정
                if not df.empty:
                    # UAR과 Z 값의 범위를 기반으로 축 범위 설정
                    ur_values = df["UAR"].dropna()
                    z_values = df["Z"].dropna()

                    if len(ur_values) > 0 and len(z_values) > 0:
                        ur_margin = (ur_values.max() - ur_values.min()) * 0.1
                        z_margin = (z_values.max() - z_values.min()) * 0.1

                        ax.set_xlim(
                            ur_values.min() - ur_margin, ur_values.max() + ur_margin
                        )
                        ax.set_ylim(
                            z_values.min() - z_margin, z_values.max() + z_margin
                        )

            except Exception as e:
                logger.info(f"      Error processing {chr_name}: {e}")
                ax.text(
                    0.5,
                    0.5,
                    f"{chr_name}\nError:\n{str(e)[:30]}...",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                    fontsize=8,
                    color="red",
                )
                ax.set_title(chr_name, fontsize=12, fontweight="bold")
                ax.set_xlabel("UAR (%)", fontsize=10)
                ax.set_ylabel("Z-score", fontsize=10)
                ax.grid(True, alpha=0.3)

        # ----------------------------------------------------
        # SCA Female 플롯 (chr23 위치)
        # ----------------------------------------------------
        row, col = divmod(22, 6)  # chr23 위치
        ax = axes[row][col]

        try:
            # female.txt 파일 로드
            female_file = os.path.join(chr_table_dir, "female.txt")

            logger.info(f"female_file : {female_file}")
            if os.path.exists(female_file):
                female_df = pd.read_csv(female_file, sep="\t")

                # logger.info(female_df.head())
                # XX Normal 샘플
                xx_samples = female_df[female_df["type"] == "XX"]
                if not xx_samples.empty:
                    ax.scatter(
                        xx_samples["UAR"],
                        xx_samples["Z"],
                        color="black",
                        s=8,
                        alpha=0.5,
                        facecolors="none",
                        edgecolors="black",
                        linewidths=0.5,
                        label="XX",
                    )

                # XXX 샘플들
                xxx = female_df[female_df["type"] == "XXX"]
                if not xxx.empty:
                    ax.scatter(
                        xxx["UAR"],
                        xxx["Z"],
                        color="orange",
                        s=8,
                        alpha=0.5,
                        facecolors="none",
                        edgecolors="orange",
                        linewidths=0.5,
                        label="XXX",
                    )
                    # marker='^', label='XXX')

                # XO 샘플들
                xo = female_df[female_df["type"] == "XO"]
                if not xo.empty:
                    ax.scatter(
                        xo["UAR"],
                        xo["Z"],
                        color="coral",
                        s=8,
                        alpha=0.5,
                        facecolors="none",
                        edgecolors="coral",
                        linewidths=0.5,
                        label="XO",
                    )
                    # marker='^', label='XO')

                # 테스트 샘플 (X 염색체 데이터)
                if "chrX" in test_sample_data:
                    test_ur_x = test_sample_data["chrX"]["UAR"]
                    test_z_x = test_sample_data["chrX"]["Z"]
                    try:
                        test_ur_x_num = float(test_ur_x)
                        test_z_x_num = float(test_z_x)
                        ax.scatter(
                            test_ur_x_num,
                            test_z_x_num,
                            color="red",
                            edgecolor="white",
                            s=50,
                            linewidths=2,
                            zorder=10,
                        )
                    except (ValueError, TypeError):
                        pass

                # 임계값 라인 (female 기준)
                add_sca_lines_to_chrx_plot(ax, group)

                ax.set_xlabel("UAR[X] (%)", fontsize=10)
                ax.set_ylabel("Z-score[X]", fontsize=10)
                # ax.set_title('SCA (Female)', fontsize=12, fontweight='bold')
                ax.set_title("chrX", fontsize=12, fontweight="bold")
                ax.legend(loc="upper right", fontsize=6)
                ax.grid(True, alpha=0.3)

            else:
                ax.text(
                    0.5,
                    0.5,
                    "chrX\nNo Data",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                    fontsize=10,
                    color="gray",
                )
                ax.set_title("chrX", fontsize=12, fontweight="bold")

        except Exception as e:
            print(f"        Error plotting SCA Female: {e}")
            ax.text(
                0.5,
                0.5,
                "SCA (Female)\nError",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=10,
                color="red",
            )
            ax.set_title("SCA (Female)", fontsize=12, fontweight="bold")

        # ----------------------------------------------------
        # SCA Male 플롯 (chr24 위치)
        # ----------------------------------------------------
        row, col = divmod(23, 6)  # chr24 위치
        ax = axes[row][col]

        try:
            # male.txt 파일 로드
            male_file = os.path.join(chr_table_dir, "male.txt")
            # male_file = os.path.join(os.path.dirname(chr_table_dir), '..', 'male.txt')

            if os.path.exists(male_file):
                male_df = pd.read_csv(male_file, sep="\t")

                # XY Normal 샘플
                xy_samples = male_df[male_df["type"] == "XY"]
                if not xy_samples.empty:
                    ax.scatter(
                        xy_samples["UAR.X"],
                        xy_samples["UAR.Y"],
                        color="black",
                        s=8,
                        alpha=0.5,
                        facecolors="none",
                        edgecolors="black",
                        linewidths=0.5,
                        label="XY",
                    )

                # XXY 샘플들
                xxy = male_df[male_df["type"] == "XXY"]
                if not xxy.empty:
                    ax.scatter(
                        xxy["UAR.X"],
                        xxy["UAR.Y"],
                        color="orange",
                        s=8,
                        alpha=0.5,
                        facecolors="none",
                        edgecolors="orange",
                        linewidths=0.5,
                        label="XXY",
                    )
                    # marker='^', label='XXY')

                # XYY 샘플들
                xyy = male_df[male_df["type"] == "XYY"]
                if not xyy.empty:
                    ax.scatter(
                        xyy["UAR.X"],
                        xyy["UAR.Y"],
                        color="green",
                        s=8,
                        alpha=0.5,
                        facecolors="none",
                        edgecolors="green",
                        linewidths=0.5,
                        label="XYY",
                    )
                    # marker='^', label='XYY')

                # XXYY 샘플들
                xxyy = male_df[male_df["type"] == "XXYY"]
                if not xxyy.empty:
                    ax.scatter(
                        xxyy["UAR.X"],
                        xxyy["UAR.Y"],
                        color="purple",
                        s=8,
                        alpha=0.5,
                        facecolors="none",
                        edgecolors="purple",
                        linewidths=0.5,
                        label="XXYY",
                    )
                    # marker='P', label='XXYY')

                # 테스트 샘플 (X, Y 염색체 데이터)
                test_ur_x = None
                test_ur_y = None
                if "chrX" in test_sample_data:
                    test_ur_x = test_sample_data["chrX"]["UAR"]
                if "chrY" in test_sample_data:
                    test_ur_y = test_sample_data["chrY"]["UAR"]

                logger.info(f"test_ur_x: {test_ur_x}, test_ur_y: {test_ur_y}")
                if test_ur_x is not None and test_ur_y is not None:
                    try:
                        test_ur_x_num = float(test_ur_x)
                        test_ur_y_num = float(test_ur_y)
                        ax.scatter(
                            test_ur_x_num,
                            test_ur_y_num,
                            color="red",
                            edgecolor="darkred",
                            s=50,
                            linewidths=2,
                            zorder=10,
                        )
                    except (ValueError, TypeError):
                        pass

                # 임계값 라인 (male 기준)
                add_sca_lines_to_chry_plot(ax, group)

                ax.set_xlabel("UAR[X] (%)", fontsize=10)
                ax.set_ylabel("UAR[Y] (%)", fontsize=10)
                # ax.set_title('SCA (Male)', fontsize=12, fontweight='bold')
                ax.set_title("chrY", fontsize=12, fontweight="bold")
                ax.legend(loc="upper right", fontsize=6)
                ax.grid(True, alpha=0.3)

            else:
                ax.text(
                    0.5,
                    0.5,
                    "chrY\nNo Data",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                    fontsize=10,
                    color="gray",
                )
                ax.set_title("chrY", fontsize=12, fontweight="bold")

        except Exception as e:
            print(f"        Error plotting SCA Male: {e}")
            ax.text(
                0.5,
                0.5,
                "SCA (Male)\nError",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=10,
                color="red",
            )
            ax.set_title("SCA (Male)", fontsize=12, fontweight="bold")

        # 전체 제목 설정
        group_name = (
            os.path.basename(threshold_path).replace("_thresholds.tsv", "").upper()
        )
        if not group_name:
            group_name = "EZD"

        plt.suptitle(
            f"Chromosome UAR vs Z-score Analysis\n{sample_name} - {group_name})",
            fontsize=16,
            fontweight="bold",
        )

        # 출력 디렉토리 생성
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        # 레이아웃 조정 및 저장
        plt.tight_layout(rect=[0, 0, 1, 0.92])
        plt.savefig(output_path, dpi=dpi, bbox_inches="tight", facecolor="white")
        plt.close()

        logger.info("  Scatter grid plot saved successfully!")

    except Exception as e:
        logger.info(f"Error creating scatter grid plot: {e}")
        import traceback

        traceback.logger.info_exc()

        # 에러 발생 시에도 빈 이미지 생성
        try:
            fig, ax = plt.subplots(1, 1, figsize=(10, 6))
            ax.text(
                0.5,
                0.5,
                f"Error creating plot:\n{str(e)}",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=12,
                color="red",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray"),
            )
            ax.set_title("Plot Generation Error", fontsize=14, fontweight="bold")
            ax.axis("off")

            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            plt.savefig(output_path, dpi=dpi, bbox_inches="tight", facecolor="white")
            plt.close()
            logger.info(f"  Error plot saved to: {output_path}")
        except:
            logger.info("  Failed to save error plot")


# Plotly. Not necessary now
def plot_ezd_interactive(ezd_df, out_path):
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=ezd_df["UAR"],
            y=ezd_df["Z"],
            mode="markers+text",
            text=ezd_df["chr"],
            textposition="top center",
            marker=dict(size=12, color="blue"),
        )
    )
    fig.update_layout(
        title="EZD Interactive Plot",
        xaxis_title="UAR",
        yaxis_title="Z-score",
        width=800,
        height=600,
    )
    fig.write_html(out_path)
