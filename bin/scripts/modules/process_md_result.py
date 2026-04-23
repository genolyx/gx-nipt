import pandas as pd
import sys
import json
import logging
import os
from pathlib import Path

logger = logging.getLogger(__name__)

# Function to check overlap between two regions
def check_overlap(start1, end1, start2, end2):
    #print(start1, end1, start2, end2)
    return max(start1, start2) <= min(end1, end2)

def check_full_inclusion(start1, end1, start2, end2):
    return start1 <= start2 and end1 >= end2

def process_wc(input_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if "# Test results" in line:
                data_start = i + 1
                break

    df = pd.read_csv(input_file, sep=r"\s+", skiprows=data_start)
    if len(df) > 0:
        df[['chr', 'start_end']] = df['location'].str.split(':', expand=True)
        df[['start', 'end']] = df['start_end'].str.split('-', expand=True)
        df['start'] = pd.to_numeric(df['start'])
        df['end'] = pd.to_numeric(df['end'])
        df['length'] = df['end'] - df['start']
        df = df[['chr', 'start', 'end', 'effect', 'z-score', 'length']]
        df.columns = ['chr', 'start', 'end', 'effect', 'zscore', 'length']

    return df

def process_wcx(input_file):
    df = pd.read_csv(input_file, sep=r"\s+")
    df['chr'] = df['chr'].astype('str')
    df['length'] = df['end'] - df['start'] + 1
    df['effect'] = df['ratio']
    df = df[['chr', 'start', 'end', 'effect', 'zscore', 'length']]

    return df

def process_wcff(input_file):
    df = pd.read_csv(input_file, sep=r"\s+")
    df['chr'] = df['chr'].astype('str')
    df['length'] = df['end'] - df['start'] + 1
    df['effect'] = df['ratio']
    df = df[['chr', 'start', 'end', 'effect', 'zscore', 'length']]

    return df

# Kenneth 250311 : Apply "Minimum length"
def compare_with_db_old(df, target_file, threshold):
    # Load the DB.bed file
    # Kenneth 250407 : Apply critical genes (Currently, it's been applied to md108 only)
    db_df = pd.read_csv(target_file, sep="\t", header=None, names=['chr', 'start', 'end', 'Disease-name', 'gain_loss', 'overlap_inside', 'minimum_length', 'cri_genes', 'cri_start', 'cri_end'])

    # Create empty lists to store results
    matches_first_8 = []
    matches_beyond_8 = []

    # Iterate through the output dataframe and check for overlaps with the DB.bed file
    for index, row in df.iterrows():
        match_found = False
        # Check first 8 rows of DB.bed
        for _, db_row in db_df.iloc[:8].iterrows():
            if db_row['overlap_inside'] == 'overlap':
                if row['chr'] == db_row['chr'] and check_overlap(row['start'], row['end'], db_row['start'], db_row['end']):
                    if (row['zscore'] < -float(threshold) and db_row['gain_loss'] == 'loss') or (row['zscore'] > float(threshold) and db_row['gain_loss'] == 'gain'):
                        row['Disease-name'] = db_row['Disease-name']
                        matches_first_8.append(row)
                        match_found = True
                        break
            elif db_row['overlap_inside'] == 'inside':
                if row['chr'] == db_row['chr'] and check_full_inclusion(row['start'], row['end'], db_row['start'], db_row['end']):
                    if (row['zscore'] < -float(threshold) and db_row['gain_loss'] == 'loss') or (row['zscore'] > float(threshold) and db_row['gain_loss'] == 'gain'):
                        row['Disease-name'] = db_row['Disease-name']
                        matches_first_8.append(row)
                        match_found = True
                        break
            else:
                continue

        
        # If no match found in first 8 rows, check the rest of DB.bed
        if not match_found:
            for _, db_row in db_df.iloc[8:].iterrows():
                if db_row['overlap_inside'] == 'overlap':
                    if row['chr'] == db_row['chr'] and check_overlap(row['start'], row['end'], db_row['start'], db_row['end']):
                        if (row['zscore'] < -float(threshold) and db_row['gain_loss'] == 'loss') or (row['zscore'] > float(threshold) and db_row['gain_loss'] == 'gain'):
                            row['Disease-name'] = db_row['Disease-name']
                            matches_beyond_8.append(row)
                            break

                elif db_row['overlap_inside'] == 'inside':
                    if row['chr'] == db_row['chr'] and check_full_inclusion(row['start'], row['end'], db_row['start'], db_row['end']):
                        if (row['zscore'] < -float(threshold) and db_row['gain_loss'] == 'loss') or (row['zscore'] > float(threshold) and db_row['gain_loss'] == 'gain'):
                            row['Disease-name'] = db_row['Disease-name']
                            matches_beyond_8.append(row)
                            break
                else:
                    continue

    # Convert result lists to DataFrames
    first_8_df = pd.DataFrame(matches_first_8)
    beyond_8_df = pd.DataFrame(matches_beyond_8)

    return first_8_df, beyond_8_df

# Kenneth 250311 : Apply "Minimum length"

def get_min_length(length, md_cfg):
    if length <= md_cfg["region1"]["less_than"]:
        return md_cfg["region1"]["min_length"]
    elif length <= md_cfg["region2"]["less_than"]:
        return md_cfg["region2"]["min_length"]
    elif length <= md_cfg["region3"]["less_than"]:
        return md_cfg["region3"]["min_length"]
    elif length > md_cfg["region4"]["greater_than"]:
        return md_cfg["region4"]["min_length"]
    else:
        return md_cfg["region3"]["min_length"]


def compare_with_target(result_df, target_file, md_cfg, threshold):

    # Load the DB.bed file
    # Kenneth 250407 : Apply critical genes (Currently, it's been applied to md108 only)
    db_df = pd.read_csv(target_file, sep="\t", header=None, names=['chr', 'start', 'end', 'Disease-name', 'gain_loss', 'overlap_inside', 'minimum_length', 'cri_genes', 'cri_start', 'cri_end'])

    md_output_df = pd.DataFrame()

    for _, result_row in result_df.iterrows():
        #print(result_row)
        chr1 = result_row["chr"]
        start1 = int(result_row["start"])
        end1 = int(result_row["end"])
        z = float(result_row["zscore"])
        detected_length = end1 - start1 + 1

        match_found = False

        for _, db_row in db_df.iterrows():
            #print(db_row)
            chr2 = str(db_row["chr"])
            start2 = int(db_row["start"])
            end2 = int(db_row["end"])
            min_len_field = db_row["minimum_length"]
            mode = db_row["overlap_inside"]
            gain_loss = db_row["gain_loss"]
            target_length = end2 - start2 + 1

            overlap = False
            if mode == "overlap":
                detected = check_overlap(start1, end1, start2, end2)
            elif mode == "inside":
                detected = check_full_inclusion(start1, end1, start2, end2)
            else:
                continue  # skip invalid rows

            if chr1 == chr2 and detected:
                #print(db_row, detected)
                #print(chr1, start1, end1, chr2, start2, end2)
                #print(f"z:{z}, threshold:{threshold}")
                # Check zscore match with gain/loss type
                if ((gain_loss == 'loss' and z < -float(threshold)) or
                    (gain_loss == 'gain' and z > float(threshold)) or
                    (gain_loss == 'both' and (z < -float(threshold) or z > float(threshold)))):

                    # Check cri_genes condition
                    cri_genes = db_row.get("cri_genes", "-")
                    if cri_genes != "-":
                        try:
                            cri_start = int(db_row["cri_start"])
                            cri_end = int(db_row["cri_end"])
                        except ValueError:
                            continue  # skip if invalid cri range

                        # Make sure the entire region is within critical range
                        detected = check_overlap(start1, end1, cri_start, cri_end)
                        if detected == True:
                            result_row["Disease-name"] = db_row["Disease-name"]
                            md_output_df = pd.concat([md_output_df, result_row.to_frame().T], ignore_index=True)
                            match_found = True
                            break  # match found
                        else:
                            continue  # doesn't fit cri_gene bounds

                    else:

                        # Determine minimum length requirement
                        if min_len_field != "-":
                            try:
                                min_len = int(min_len_field)
                            except ValueError:
                                continue
                        else:
                            min_len = get_min_length(target_length, md_cfg)

                        if detected_length >= min_len:
                            logger.info(f"[MATCHED] detected_length:{detected_length}, min_len:{min_len}")
                            result_row["Disease-name"] = db_row["Disease-name"]
                            md_output_df = pd.concat([md_output_df, result_row.to_frame().T], ignore_index=True)
                            match_found = True
                            break  # One good match is enough
                        else:
                            logger.info(f"[SHORT LENGTH] detected_length:{detected_length}, min_len:{min_len}")

    return md_output_df


def process_microdeletion_result(sample_name, tool_name, input_file, tag, target_file, config_file, md_select, output_dir, threshold):
    if tool_name == "WC":
        result_df = process_wc(input_file)
    elif tool_name == "WCX":
        result_df = process_wcx(input_file)
    elif tool_name == "WCFF":
        result_df = process_wcff(input_file)
    else:
        logger.info(f"[ERROR] Invalid tool name: {tool_name}. Use 'WC' or 'WCX' or 'WCFF'.")
        return False

    md_file_dict = {
        "MD_Target_8": "md8",
        "MD_Target_108": "md108",
        "MD_Target_320": "md320",
        "MD_Target_87": "md87",
        "MD_Target_141": "md141"
    }

    if md_select not in md_file_dict:
        logger.info(f"[ERROR] Invalid md_select: {md_select}.")
        return False

    output_filename_md = f"{sample_name}_{tool_name}_{tag}_{md_file_dict[md_select]}.tsv"

    try:
        with open(config_file) as f:
            cfg_json = json.load(f)
    except Exception as e:
        logger.info(f"[ERROR] Failed to load config file: {e}")
        return False

    if md_select not in cfg_json:
        logger.info(f"[ERROR] {md_select} not found in config file.")
        return False

    md_cfg = cfg_json[md_select]

    if result_df.empty:
        logger.info(f"[INFO] No result found for {tool_name} in {md_select} : {sample_name}.")
        return False

    output_df = compare_with_target(result_df, target_file, md_cfg, threshold)

    if output_df.empty:
        logger.info(f"[INFO] No matching microdeletion targets found for {tool_name} in {md_select}:{sample_name}.")
        return False

    output_path = os.path.join(output_dir, tag, output_filename_md)
    output_df.to_csv(output_path, sep='\t', index=False)
    logger.info(f"[SUCCESS] {output_path} was generated successfully!")

    return True

def run_microdeletion_decision_pipeline(
    sample_name, labcode, config, analysis_dir, output_dir, bed_dir,
    config_file=None,
):
    """
    Run microdeletion decision for all methods (WC, WCX) and types (orig, fetus, mom)
    across all MD targets defined in the config.

    Args:
        config_file: Optional explicit path to pipeline_config.json. If not
                     provided, falls back to the ken-nipt legacy location.
    """

    methods = ["WC", "WCX"]
    types = ["orig", "fetus", "mom"]
    md_targets = ["MD_Target_8", "MD_Target_87", "MD_Target_108", "MD_Target_320", "MD_Target_141"]

    wc_input_files = {
        "orig": f"{analysis_dir}/{sample_name}/Output_WC/orig/{sample_name}.wc.orig.report.txt",
        "fetus": f"{analysis_dir}/{sample_name}/Output_WC/fetus/{sample_name}.wc.fetus.report.txt",
        "mom": f"{analysis_dir}/{sample_name}/Output_WC/mom/{sample_name}.wc.mom.report.txt",
    }

    wcx_input_files = {
        "orig": f"{analysis_dir}/{sample_name}/Output_WCX/orig/{sample_name}.wcx.orig_aberrations.bed",
        "fetus": f"{analysis_dir}/{sample_name}/Output_WCX/fetus/{sample_name}.wcx.fetus_aberrations.bed",
        "mom": f"{analysis_dir}/{sample_name}/Output_WCX/mom/{sample_name}.wcx.mom_aberrations.bed",
    }

    if config_file is None:
        # Legacy ken-nipt fallback
        config_file = f"/Work/NIPT/config/{labcode}/pipeline_config.json"

    for md_key in md_targets:
        bed_file = f"{bed_dir}/common/{config[md_key]['bed']}"

        for method in methods:
            for tag in types:
                # Determine input file and output dir
                if method == "WC":
                    input_file = wc_input_files[tag]
                    output_subdir = "Output_WC"
                elif method == "WCX":
                    input_file = wcx_input_files[tag]
                    output_subdir = "Output_WCX"
                else:
                    continue

                output_path = f"{analysis_dir}/{sample_name}/{output_subdir}"
                threshold = config[method].get(f"{tag}_call_threshold")

                if not Path(input_file).exists():
                    logging.warning(f"[SKIP] Input file missing: {input_file}")
                    continue

                logging.info(f"[MD] {method}-{tag} on {md_key} (Threshold: {threshold})")
                success = process_microdeletion_result(
                    sample_name=sample_name,
                    tool_name=method,
                    input_file=input_file,
                    tag=tag,
                    target_file=bed_file,
                    config_file=config_file,
                    md_select=md_key,
                    output_dir=output_path,
                    threshold=threshold
                )

                if success:
                    logging.info(f"[PASS] MD call: {method}-{tag} {md_key}")
                #else:
                #    logging.info(f"[NOT MATCHED] MD call: {method}-{tag} {md_key}")
    return True


def run_microdeletion_test_decision_pipeline(
    sample_name, labcode, config, analysis_dir, output_dir, bed_dir,
    types=None, md_targets=None
):
    """
    Run microdeletion decision for specified methods, types, and MD targets.
    This is a flexible version for testing/artificial samples.
    
    Args:
        sample_name: Sample ID
        labcode: Lab code
        config: Configuration dictionary
        analysis_dir: Analysis directory path
        output_dir: Output directory path
        bed_dir: BED files directory path
        types: List of BAM types to process (default: ["orig", "fetus"])
        md_targets: List of MD targets to process (default: ["MD_Target_8"])
    
    Returns:
        bool: True if successful
    """
    
    # Default values for testing
    if types is None:
        types = ["orig", "fetus"]
    
    if md_targets is None:
        md_targets = ["MD_Target_8"]
    
    methods = ["WC", "WCX"]

    wc_input_files = {
        "orig": f"{analysis_dir}/{sample_name}/Output_WC/orig/{sample_name}.wc.orig.report.txt",
        "fetus": f"{analysis_dir}/{sample_name}/Output_WC/fetus/{sample_name}.wc.fetus.report.txt",
        "mom": f"{analysis_dir}/{sample_name}/Output_WC/mom/{sample_name}.wc.mom.report.txt",
    }

    wcx_input_files = {
        "orig": f"{analysis_dir}/{sample_name}/Output_WCX/orig/{sample_name}.wcx.orig_aberrations.bed",
        "fetus": f"{analysis_dir}/{sample_name}/Output_WCX/fetus/{sample_name}.wcx.fetus_aberrations.bed",
        "mom": f"{analysis_dir}/{sample_name}/Output_WCX/mom/{sample_name}.wcx.mom_aberrations.bed",
    }

    config_file = f"/Work/NIPT/config/{labcode}/pipeline_config.json"
    
    # Load config file to get MD target and threshold information
    try:
        with open(config_file) as f:
            config_data = json.load(f)
    except Exception as e:
        logging.error(f"[ERROR] Failed to load config file {config_file}: {e}")
        return False

    for md_key in md_targets:
        if md_key not in config_data:
            logging.warning(f"[SKIP] MD target {md_key} not found in config")
            continue
        
        bed_file = f"{bed_dir}/common/{config_data[md_key]['bed']}"

        for method in methods:
            for tag in types:
                # Determine input file and output dir
                if method == "WC":
                    input_file = wc_input_files.get(tag)
                    output_subdir = "Output_WC"
                elif method == "WCX":
                    input_file = wcx_input_files.get(tag)
                    output_subdir = "Output_WCX"
                else:
                    continue

                if input_file is None:
                    logging.warning(f"[SKIP] Unknown type: {tag}")
                    continue

                output_path = f"{analysis_dir}/{sample_name}/{output_subdir}"
                threshold = config_data.get(method, {}).get(f"{tag}_call_threshold")

                if not Path(input_file).exists():
                    logging.warning(f"[SKIP] Input file missing: {input_file}")
                    continue

                logging.info(f"[MD] {method}-{tag} on {md_key} (Threshold: {threshold})")
                success = process_microdeletion_result(
                    sample_name=sample_name,
                    tool_name=method,
                    input_file=input_file,
                    tag=tag,
                    target_file=bed_file,
                    config_file=config_file,
                    md_select=md_key,
                    output_dir=output_path,
                    threshold=threshold
                )

                if success:
                    logging.info(f"[PASS] MD call: {method}-{tag} {md_key}")
                #else:
                #    logging.info(f"[NOT MATCHED] MD call: {method}-{tag} {md_key}")
    
    return True


def copy_md_output_files(sample_id, analysis_dir, output_dir):
    """Copy MD output images and reports to unified output directory."""
    groups = ['orig', 'fetus', 'mom']
    wc_dir = os.path.join(analysis_dir, sample_id, 'Output_WC')
    wcx_dir = os.path.join(analysis_dir, sample_id, 'Output_WCX')
    for group in groups:
        wc_img = os.path.join(wc_dir, group, f"{sample_id}.{group}.png")
        wc_rep = os.path.join(wc_dir, group, f"{sample_id}.{group}.report.txt")
        wcx_bed = os.path.join(wcx_dir, group, f"{sample_id}.wcx.{group}_aberrations.bed")
        wcx_img = os.path.join(wcx_dir, group, f"{sample_id}.wcx.{group}.plots", "genome_wide.png")

        dest_dir = os.path.join(output_dir, sample_id, "Output_MD", group)
        os.makedirs(dest_dir, exist_ok=True)
        for src in [wc_img, wc_rep, wcx_bed, wcx_img]:
            if os.path.exists(src):
                shutil.copy2(src, os.path.join(dest_dir, os.path.basename(src)))
   
    return True

# =========================================================
#  Adapter CLI: Nextflow-compatible entry point
#
#  Used by modules/microdeletion.nf. Orchestrates
#  run_microdeletion_decision_pipeline() across all WC/WCX
#  results and concatenates the per-(method,group,md_key) TSVs
#  into a single <sample>.md_result.tsv for downstream reporting.
# =========================================================
def _legacy_positional_main():
    """Legacy ken-nipt positional CLI kept for backwards compatibility."""
    if len(sys.argv) != 10:
        logger.error(
            "Usage: python process_md_result.py <sample_name> <tool_name> "
            "<input_file> <tag> <target_file> <config_file> <md_key> "
            "<output_dir> <threshold>"
        )
        sys.exit(1)

    sample_name = sys.argv[1]
    tool_name   = sys.argv[2]
    input_file  = sys.argv[3]
    tag         = sys.argv[4]
    target_file = sys.argv[5]
    config_file = sys.argv[6]
    md_select   = sys.argv[7]
    output_dir  = sys.argv[8]
    threshold   = sys.argv[9]

    process_microdeletion_result(
        sample_name, tool_name, input_file, tag, target_file,
        config_file, md_select, output_dir, threshold,
    )


def _nf_main():
    import argparse
    import glob

    parser = argparse.ArgumentParser(
        description=(
            "Microdeletion adapter CLI (Nextflow). Runs per-(method,group,"
            "MD_Target) decisions and emits a single summary TSV."
        )
    )
    parser.add_argument("--sample",       required=True)
    parser.add_argument("--labcode",      required=True)
    parser.add_argument("--config",       required=True, help="pipeline_config.json")
    parser.add_argument("--analysis-dir", required=True,
                        help="Root analysis dir where Output_WC / Output_WCX live")
    parser.add_argument("--bed-dir",      required=True,
                        help="Directory containing MD target BED files "
                             "(expects a 'common/' subdir)")
    parser.add_argument("--output",       required=True,
                        help="Path to write unified md_result.tsv")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    if not os.path.exists(args.config):
        logger.error(f"Config file not found: {args.config}")
        sys.exit(2)
    with open(args.config) as fh:
        cfg = json.load(fh)

    # The per-md/method/group TSVs are written under
    # <analysis_dir>/<sample>/Output_WC(X)/<group>/*.tsv by
    # process_microdeletion_result().
    run_microdeletion_decision_pipeline(
        sample_name=args.sample,
        labcode=args.labcode,
        config=cfg,
        analysis_dir=args.analysis_dir,
        output_dir=args.analysis_dir,
        bed_dir=args.bed_dir,
        config_file=args.config,
    )

    # Aggregate all md*.tsv files that the per-call routine produced
    # into a single summary TSV (Nextflow's declared single output).
    summary_rows = []
    sample_out_root = os.path.join(args.analysis_dir, args.sample)
    patterns = [
        os.path.join(sample_out_root, "Output_WC",  "*", f"{args.sample}_WC_*_md*.tsv"),
        os.path.join(sample_out_root, "Output_WCX", "*", f"{args.sample}_WCX_*_md*.tsv"),
    ]
    for pat in patterns:
        for fp in sorted(glob.glob(pat)):
            try:
                sub = pd.read_csv(fp, sep="\t")
            except Exception as e:
                logger.warning(f"Skipping unreadable {fp}: {e}")
                continue
            if sub.empty:
                continue
            # Tag each row with its provenance (method, group, md_key)
            basename = os.path.basename(fp).replace(".tsv", "")
            parts = basename.split("_")
            # <sample>_<method>_<tag>_<md_key>
            if len(parts) >= 4:
                method = parts[-3]
                tag    = parts[-2]
                mdkey  = parts[-1]
            else:
                method, tag, mdkey = "NA", "NA", "NA"
            sub.insert(0, "md_key", mdkey)
            sub.insert(0, "group",  tag)
            sub.insert(0, "method", method)
            sub.insert(0, "sample", args.sample)
            summary_rows.append(sub)

    if summary_rows:
        out_df = pd.concat(summary_rows, ignore_index=True)
    else:
        out_df = pd.DataFrame(
            columns=["sample", "method", "group", "md_key",
                     "chr", "start", "end", "effect", "zscore", "length"]
        )

    os.makedirs(os.path.dirname(os.path.abspath(args.output)) or ".", exist_ok=True)
    out_df.to_csv(args.output, sep="\t", index=False)
    logger.info(
        f"[MD adapter] wrote {len(out_df)} rows to {args.output} "
        f"from {len(summary_rows)} per-MD files."
    )


if __name__ == "__main__":
    argv = sys.argv[1:]
    _nf_flags = ("--sample", "--labcode", "--analysis-dir", "--bed-dir", "--output")
    if any(flag in argv for flag in _nf_flags):
        _nf_main()
        sys.exit(0)
    # Legacy positional fallback (kept for backwards compatibility)
    _legacy_positional_main()
    sys.exit(0)
