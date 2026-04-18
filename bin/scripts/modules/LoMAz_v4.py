import os
import glob
import sys
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import logging

class HeaderColorFormatter(logging.Formatter):
    COLORS = {
        'DEBUG': '\033[36m',    # Cyan
        'INFO': '\033[32m',     # Green
        'WARNING': '\033[33m',  # Yellow
        'ERROR': '\033[31m',    # Red
        'CRITICAL': '\033[41m', # Red background
    }
    RESET = '\033[0m'

    def format(self, record):
        color = self.COLORS.get(record.levelname, self.RESET)
        header = f"{color}[{record.levelname} {record.filename}:{record.lineno}]{self.RESET}"
        return f"{header} {record.getMessage()}"

# Define standard chromosome order
CHROMS = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']

# ---------------------------------------------------
# chromosome order
# ---------------------------------------------------
def chrom_order_key(chrom):
    try:
        return int(chrom.replace('chr', '').replace('X', '23').replace('Y', '24'))
    except:
        return 99

# ---------------------------------------------------
# Heatmap plot for MAD matrix visualization (chroms)
# ---------------------------------------------------
def plot_chrom_mad_heatmap(mad_mat, chroms, out_path):
    """
    Plot heatmap of MAD values between chromosome ratios.
    Helps diagnose low-MAD or zero-MAD issues.
    """
    plt.figure(figsize=(9, 8))
    # Plot heatmap with chromosomes as labels on both axes
    sns.heatmap(mad_mat, xticklabels=chroms, yticklabels=chroms,
                cmap='YlGnBu', vmin=0, vmax=0.5, square=True, annot=False,
                linewidths=0.3, linecolor='gray')
    plt.title("MAD Matrix Heatmap", loc='center')
    plt.tight_layout()
    plt.savefig(out_path)
    logging.info(f"Saved MAD heatmap to {out_path}")

# ------------------------------------------------
# Heatmap plot for MAD matrix visualization (bin)
# ------------------------------------------------
def plot_bin_mad_vector(mad_vector, chr_ids, bin_ids, output_path):
    """
    Plot a bin-level MAD heatmap with sorted chromosome order and clean x-ticks.
    """
    df = pd.DataFrame({'chr': chr_ids, 'bin': bin_ids, 'mad': mad_vector})

    # Sort chromosome and bin order
    df['chr'] = pd.Categorical(df['chr'], ordered=True,
                               categories=sorted(df['chr'].unique(), key=chrom_order_key))
    df['bin'] = pd.Categorical(df['bin'], ordered=True,
                               categories=sorted(df['bin'].unique(), key=lambda b: (
                                   chrom_order_key(b.split(':')[0]),
                                   int(b.split(':')[1].replace('Mb', '')))
                               ))

    # Pivot for heatmap
    pivot = df.pivot(index='chr', columns='bin', values='mad')

    # Plot
    plt.figure(figsize=(24, 6))
    ax = sns.heatmap(pivot, cmap='YlGnBu', linewidths=0.1, linecolor='gray')

    # X-ticks: one per chromosome
    #xtick_positions = [i for i, b in enumerate(pivot.columns) if b.endswith('0Mb')]
    #xtick_labels = [b.split(':')[0] for b in pivot.columns if b.endswith('0Mb')]

    # Set xticks: one per chromosome (first appearance)
    xtick_positions = []
    xtick_labels = []
    seen_chr = set()

    for i, bin_label in enumerate(pivot.columns):
        chrom = bin_label.split(':')[0]
        if chrom not in seen_chr:
            xtick_positions.append(i)
            xtick_labels.append(chrom)
            seen_chr.add(chrom)

    ax.set_xticks(xtick_positions)
    ax.set_xticklabels(xtick_labels, rotation=90)

    plt.title("Bin-level MAD Heatmap", loc='center')
    plt.xlabel("Chromosome")
    plt.ylabel("Chromosome")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


# ---------------------------------------------
# Barplot of mean Z-scores per chromosome
# ---------------------------------------------
def plot_mean_z_per_chr(zmat, chroms, output):
    """
    Compute and plot the mean Z-score for each chromosome.
    Excludes diagonal self-comparisons (i == j).
    """
    # Remove diagonal elements by deleting columns where i == j and compute mean across rows
    #mean_z = np.mean(np.delete(zmat, np.arange(len(chroms)), axis=1), axis=1)
    #mean_z = np.nanmean(np.delete(zmat, np.arange(len(chroms)), axis=1), axis=1)
    mean_z = np.array([
        np.nanmean([zmat[i, j] for j in range(len(chroms)) if i != j])
        for i in range(len(chroms))
    ])
    plt.figure(figsize=(12, 4))
    # Barplot of mean Z-scores per chromosome
    sns.barplot(x=chroms, y=mean_z)
    plt.xticks(rotation=90)
    plt.ylabel("Mean Z-score vs other chromosomes")
    plt.title("LoMA-Z Mean Z-score per Chromosome", loc='center')
    plt.tight_layout()
    plt.savefig(output)
    logging.info(f"Saved mean Z-score barplot to {output}")

# ---------------------------------------------
# Load .Normalization.txt file and filter chromosomes
# ---------------------------------------------
def parse_normalization(file_path, value_column):
    """
    Load Normalization.txt data and return only CHROMS chromosomes.
    """
    df = pd.read_csv(file_path, sep='\t')
    # Filter to only standard chromosomes defined in CHROMS
    df = df[df['chr'].isin(CHROMS)]
    return df

# ---------------------------------------------
# Compute log2 ratio of chromosome/bin coverage
# ---------------------------------------------
def log_ratios(counts):
    """
    Compute log2(x / (total - x)) ratio for normalization.
    """
    total = counts.sum()
    # Add small constant to avoid division by zero or log of zero
    return np.log2(counts / (total - counts) + 1e-6)  # Small epsilon to prevent division by zero)

# ---------------------------------------------
# Collapse counts to chromosome level
# ---------------------------------------------
def summarize_chrom_level(df, value_column):
    """
    Sum per-chromosome counts for chromosome-level analysis.
    """
    # Group by chromosome and sum the counts, reindex to ensure all CHROMS present
    return df.groupby('chr')[value_column].sum().reindex(CHROMS).fillna(0)

# ---------------------------------------------
# Build Reference data
# ---------------------------------------------
def build_reference(ref_dir, value_column, output_dir):
    """
    Unified entry point to build both chromosome-level and bin-level references.
    """
    logging.info("Building chromosome-level reference...")
    build_reference_matrix_chr(ref_dir, value_column, output_dir=output_dir)

    #logging.info("Building bin-level reference...")
    #build_reference_matrix_bin(ref_dir, value_column, output_dir=output_dir)

# ---------------------------------------------
# Build Reference data (chromosome-wise)
# ---------------------------------------------
def build_reference_matrix_chr(ref_dir, value_column, output_dir, z_mode='loma'):
    groups = ['orig', 'fetus', 'mom']
    sexes = ['M', 'F']

    for group in groups:
        all_files = []
        for sex_dir in sexes:
            group_path = os.path.join(ref_dir, group, sex_dir)
            if os.path.exists(group_path):
                all_files += glob.glob(os.path.join(group_path, '*.10mb.txt'))
        if not all_files:
            logging.warning(f"No shared files found for group {group}")
            continue

        # Insert z_mode conditional here
        if z_mode == 'covar':
            # Compute pairwise ratio matrix: ratio[i]/ratio[j]
            shared_ratios = []
            for f in all_files:
                df = parse_normalization(f, value_column)
                counts = summarize_chrom_level(df, value_column)
                ratios = counts.values / (counts.sum() - counts.values + 1e-6)  # Small epsilon to prevent division by zero)
                shared_ratios.append(ratios)
            shared_ratios = np.stack(shared_ratios)
            mean_mat = np.zeros((24, 24))
            sd_mat = np.zeros((24, 24))
            for i in range(24):
                for j in range(24):
                    if i == j:
                        continue
                    pairwise = shared_ratios[:, i] / shared_ratios[:, j]
                    mean = np.mean(pairwise)
                    std = np.std(pairwise)
                    mean_mat[i, j] = mean
                    sd_mat[i, j] = std + 1e-6
            for sex in sexes:
                col_tag = value_column.replace('.', '_')
                out_file = os.path.join(output_dir, f'reference_covar_{group}_{sex}_chrom_{col_tag}.npz')
                np.savez(out_file, mean=mean_mat, sd=sd_mat, chroms=CHROMS)
                logging.info(f"Saved covar-level reference to: {out_file}")
            continue  # skip LoMAZ generation for covar mode

        shared_data = []
        for f in all_files:
            #logging.info(f)
            df = parse_normalization(f, value_column)
            counts = summarize_chrom_level(df, value_column)
            ratios = log_ratios(counts)
            shared_data.append(ratios.values)
        shared_data = np.stack(shared_data)

        auto_med = np.zeros((24, 24))
        auto_mad = np.zeros((24, 24))
        for i in range(22):
            for j in range(22):
                if i == j:
                    continue
                diffs = shared_data[:, i] - shared_data[:, j]
                med = np.median(diffs)
                mad = np.median(np.abs(diffs - med)) + 0.05 * abs(med)  # Added to stabilize MAD against low variability bins
                auto_med[i, j] = med
                auto_mad[i, j] = mad

        for sex in sexes:
            sex_dir = os.path.join(ref_dir, group, sex)
            sex_files = glob.glob(os.path.join(sex_dir, '*.10mb.txt')) if os.path.exists(sex_dir) else []
            if not sex_files:
                logging.warning(f"No sex-specific files for {group}/{sex}")
                continue

            sex_data = []
            for f in sex_files:
                df = parse_normalization(f, value_column)
                counts = summarize_chrom_level(df, value_column)
                ratios = log_ratios(counts)
                sex_data.append(ratios.values)
            sex_data = np.stack(sex_data)

            med_mat = np.copy(auto_med)
            mad_mat = np.copy(auto_mad)

            for i in range(24):
                for j in range(24):
                    if i < 22 and j < 22:
                        continue
                    if i == j:
                        continue

                    diffs = sex_data[:, i] - sex_data[:, j]
                    med = np.median(diffs)

                    if sex == 'F' and (i == 23 or j == 23):
                        mad = 1.0
                        med = 0.0
                    else:
                        mad = np.median(np.abs(diffs - med))
                        if sex == 'M' and (i == 23 or j == 23):
                            mad = np.clip(mad, 0.05, None)
                        else:
                            mad = mad + 0.05 * abs(med)  # Added to stabilize MAD against low variability bins

                    med_mat[i, j] = med
                    mad_mat[i, j] = mad

            col_tag = value_column.replace('.', '_')
            out_file = os.path.join(output_dir, f'reference_lomaz_{group}_{sex}_chrom_{col_tag}.npz')
            np.savez(out_file, median=med_mat, mad=mad_mat, chroms=CHROMS)
            plot_chrom_mad_heatmap(mad_mat, CHROMS, out_file.replace('.npz', '_mad.png'))
            plot_chrom_mad_heatmap(med_mat, CHROMS, out_file.replace('.npz', '_med.png'))
            logging.info(f"Saved chrom-level reference to: {out_file}")

# ---------------------------------------------
# Build Reference data (10mb-bin)
# ---------------------------------------------
def build_reference_matrix_bin(ref_dir, value_column, output_dir):
    groups = ['orig', 'fetus', 'mom']
    sexes = ['M', 'F']

    for group in groups:
        shared_files = []
        for sex_dir in sexes:
            group_path = os.path.join(ref_dir, group, sex_dir)
            if os.path.exists(group_path):
                shared_files += glob.glob(os.path.join(group_path, '*.10mb.txt'))
        if not shared_files:
            logging.warning(f"No shared bin files for group {group}")
            continue

        all_ratios = []
        bin_ids = []
        chr_ids = []

        for f in shared_files:
            df = parse_normalization(f, value_column)
            df = df.sort_values(['chr', 'start'])
            df['log_ratio'] = log_ratios(df[value_column])
            all_ratios.append(df['log_ratio'].values)
            if not bin_ids:
                bin_ids = [f"{row['chr']}:{row['start']//1000000}Mb" for _, row in df.iterrows()]
                chr_ids = df['chr'].tolist()

        all_ratios = np.stack(all_ratios)

        shared_med = np.median(all_ratios, axis=0)
        shared_mad = np.median(np.abs(all_ratios - shared_med), axis=0)
        shared_mad += 0.05 * np.abs(shared_med)

        for sex in sexes:
            sex_dir = os.path.join(ref_dir, group, sex)
            sex_files = glob.glob(os.path.join(sex_dir, '*.10mb.txt')) if os.path.exists(sex_dir) else []
            if not sex_files:
                logging.warning(f"No sex-specific bin files for {group}/{sex}")
                continue

            sex_data = []
            for f in sex_files:
                df = parse_normalization(f, value_column)
                df = df.sort_values(['chr', 'start'])
                df['log_ratio'] = log_ratios(df[value_column])
                sex_data.append(df['log_ratio'].values)
            sex_data = np.stack(sex_data)

            #logging.info(df.head())

            final_med = np.copy(shared_med)
            final_mad = np.copy(shared_mad)

            for i, chrom in enumerate(chr_ids):
                if chrom == 'chrY':
                    diffs = sex_data[:, i]
                    med = np.median(diffs)
                    if sex == 'F':
                        final_med[i] = 0.0
                        final_mad[i] = 1.0
                    else:
                        mad = np.median(np.abs(diffs - med))
                        final_med[i] = med
                        final_mad[i] = np.clip(mad, 0.05, None)

            col_tag = value_column.replace('.', '_')
            out_file = os.path.join(output_dir, f'reference_lomaz_{group}_{sex}_bin_{col_tag}.npz')
            np.savez(out_file, median=final_med, mad=final_mad,
                     bin_ids=np.array(bin_ids), chr_ids=np.array(chr_ids))
            plot_bin_mad_vector(final_mad, chr_ids, bin_ids, out_file.replace('.npz', '_mad.png'))
            plot_bin_mad_vector(final_med, chr_ids, bin_ids, out_file.replace('.npz', '_med.png'))
            logging.info(f"Saved bin-level reference to: {out_file}")

# ---------------------------------------------
# Compute chromosome-level Z matrix
# ---------------------------------------------
def compute_chrom_z_matrix(test_file, value_column, ref_file):
    """
    Compute Z-score matrix between chromosomes using reference medians and MADs or covariance-type reference.
    """
    ref = np.load(ref_file)
    chroms = ref['chroms']
    df = parse_normalization(test_file, value_column)
    counts = summarize_chrom_level(df, value_column)

    # Detect file type from keys
    if 'median' in ref and 'mad' in ref:
        ratios = log_ratios(counts)
        med_mat = ref['median']
        mad_mat = ref['mad']
        zmat = np.zeros((24, 24))
        for i in range(24):
            for j in range(24):
                if i == j:
                    continue
                diff = ratios.iloc[i] - ratios.iloc[j]
                mad = mad_mat[i, j]
                z = (diff - med_mat[i, j]) / mad if mad > 0 else 0
                zmat[i, j] = z
    elif 'mean' in ref and 'sd' in ref:
        mean_mat = ref['mean']
        sd_mat = ref['sd']
        values = counts.values
        zmat = np.zeros((24, 24))
        for i in range(24):
            for j in range(24):
                if i == j:
                    continue
                ratio = values[i] / (values[j] + 1e-6)  # Small epsilon to prevent division by zero)
                z = (ratio - mean_mat[i, j]) / sd_mat[i, j] if sd_mat[i, j] > 0 else 0
                zmat[i, j] = z
    else:
        raise ValueError("Unknown reference file format")
    return zmat, chroms

# ---------------------------------------------
# Plot heatmap of chromosome-level Z-scores
# ---------------------------------------------
def plot_chrom_zscore_heatmap(zmat, chroms, args):
    """
    Plot chromosome-level Z-score heatmap with colorbar.
    """

    # Extract sample name from test_file
    sample_basename = os.path.basename(args.test_file).split('.')[0]

    fig, ax = plt.subplots(figsize=(9, 9))
    sns.heatmap(zmat,
                xticklabels=chroms,
                yticklabels=chroms,
                cmap='RdBu_r',
                vmin=-3,
                vmax=3,
                center=0,
                square=True,
                linewidths=0.1,
                linecolor='lightgray',
                cbar=False,
                ax=ax)

    # Custom colorbar on the right
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.2)
    im = ax.collections[0]
    plt.colorbar(im, cax=cax, label="Z-score")

    ax.set_title(f"{sample_basename} [{args.group}, {args.value_column}]")
    #plt.title("LoMA-Z Chromosome Heatmap", loc='center')
    #plt.xlabel("Chromosome")
    #plt.ylabel("Chromosome")

    plt.tight_layout()

    # Create descriptive output filename
    output_name = f"{sample_basename}.chr_{args.group}_{args.value_column.replace('.', '_')}.png"
    output_path = os.path.join(args.output_dir, output_name)

    plt.savefig(output_path)
    plt.close()
    logging.info(f"Saved chromosome-level zscore heatmap: {output_path}")

# ---------------------------------------------
# Compute bin-level Z-scores
# ---------------------------------------------
def compute_bin_z_matrix(test_file, value_column, ref_file):
    """
    Compute Z-scores for each bin using reference median and MAD.
    chrY values are masked for female samples.
    """
    ref = np.load(ref_file)
    med = ref['median']
    mad = ref['mad']
    bin_ids = ref['bin_ids']
    chr_ids = ref['chr_ids']

    df = parse_normalization(test_file, value_column)
    # Sort bins by chromosome and start position to align with reference
    df = df.sort_values(['chr', 'start'])
    # Compute log ratios for test sample bins
    log_rat = log_ratios(df[value_column])
    # Compute Z-scores per bin
    z = (log_rat - med) / mad
    # Replace divisions by zero MAD with zero Z-score
    z = np.where(mad == 0, 0, z)
    # Replace any NaNs (e.g., from invalid operations) with zero
    z[np.isnan(z)] = 0
    # Mask chrY bins in female samples by setting Z to zero (chrY not present in females)
    if 'F' in ref_file:
        z = np.where(np.array(chr_ids) == 'chrY', 0, z)
    # Add Z-scores and bin labels to dataframe for downstream analysis or plotting
    df['z'] = z
    df['bin_label'] = [f"{row['chr']}:{row['start']//1000000}Mb" for _, row in df.iterrows()]
    return df[['chr', 'bin_label', 'z']]

# ------------------------------------------------
# Heatmap plot for Z-score of 10mb bin
# ------------------------------------------------
def plot_bin_zscore_heatmap(z_df, args):
    """
    Plot bin-level Z-score heatmap with chromosome rows and bins grouped on x-axis.
    Only one x-tick per chromosome is shown.
    """
    logging.info("plot_bin_zscore_heatmap")

    sample_basename = os.path.basename(args.test_file).split('.')[0]

    # Pivot dataframe: rows = chromosomes, columns = bin labels
    pivot = z_df.pivot(index='chr', columns='bin_label', values='z')

    # Sort chromosome order naturally
    def chrom_order_key(chrom):
        base = chrom.replace("chr", "")
        if base == "X":
            return 23
        elif base == "Y":
            return 24
        try:
            return int(base)
        except ValueError:
            return 99

    pivot = pivot.reindex(sorted(pivot.index, key=chrom_order_key))

    # Sort columns (bins) by chromosome and numeric position
    pivot.columns = pd.Categorical(
        pivot.columns,
        categories=sorted(pivot.columns, key=lambda b: (
            chrom_order_key(b.split(':')[0]),
            int(b.split(':')[1].replace('Mb', '')))
        ),
        ordered=True
    )
    pivot = pivot[pivot.columns]  # Apply new column order

    # Create the heatmap
    plt.figure(figsize=(24, 6))
    ax = sns.heatmap(
        pivot,
        cmap='RdBu_r',
        center=0,
        vmin=-1,
        vmax=1,
        linewidths=0.1,
        linecolor='gray',
        cbar_kws={"label": "Z-score"}
    )

    # Show only one x-tick per chromosome
    xtick_positions = []
    xtick_labels = []
    seen_chr = set()
    for i, b in enumerate(pivot.columns):
        chrom = b.split(':')[0]
        if chrom not in seen_chr:
            xtick_positions.append(i)
            xtick_labels.append(chrom)
            seen_chr.add(chrom)

    ax.set_xticks(xtick_positions)
    ax.set_xticklabels(xtick_labels, rotation=90)

    ax.set_title(f"{sample_basename} [{args.group}, {args.value_column}]")
    #plt.title("LoMA-Z Bin Heatmap", loc='center')
    #plt.xlabel("Chromosome")
    #plt.ylabel("Chromosome")
    plt.tight_layout()

     # Auto-generate output path
    output_name = f"{sample_basename}.bin_{args.group}_{args.value_column.replace('.', '_')}.png"
    output_path = os.path.join(args.output_dir, output_name)

    plt.savefig(output_path)
    plt.close()
    logging.info(f"Saved bin-level zscore heatmap: {output_path}")

# ---------------------------------------------
# CLI entry point and parser
# ---------------------------------------------
def main():
    if len(sys.argv) == 1:
        logging.info("Use --help for usage")
        sys.exit(1)

    parser = argparse.ArgumentParser(description='LoMA-Z: Log & MAD-adjusted Z-score based CNV/NIPT detection tool')
    parser.add_argument('--mode', choices=['reference', 'analyze_chrom', 'analyze_bin'], required=True)
    parser.add_argument('--ref_dir', help='Directory with reference sample folders')
    parser.add_argument('--test_file', help='Sample Normalization.txt file')
    parser.add_argument('--group', choices=['orig', 'fetus', 'mom'])
    parser.add_argument('--gender', choices=['M', 'F'])
    parser.add_argument('--ref_file', help='Reference .npz file')
    parser.add_argument('--value_column', default='cor.gc')
    parser.add_argument('--ref_type', choices=['chrom', 'bin'])
    parser.add_argument('--output_dir', default='.')
    parser.add_argument('--output', default='lomaz_output.png')
    parser.add_argument('--debug', action='store_true')

    args = parser.parse_args()
    #logging.basicConfig(level=logging.DEBUG if args.debug else logging.INFO, format='[%(levelname)s] %(message)s')
    # Set up color handler
    handler = logging.StreamHandler()
    handler.setFormatter(HeaderColorFormatter())
    logging.root.handlers = [handler]
    logging.root.setLevel(logging.DEBUG if args.debug else logging.INFO)

    if args.mode == 'reference':
        # Build reference matrices from provided directory and type
        assert args.ref_dir and args.ref_type, 'Error: --ref_dir and --ref_type are required for reference mode.'
        build_reference(args.ref_dir, args.value_column, args.output_dir)

    elif args.mode == 'analyze_chrom':
        if not args.ref_file:
            # Auto-select reference file if not provided, based on group and gender
            assert args.ref_dir and args.group and args.gender, 'Error: --ref_dir, --group, and --gender are required to auto-select reference.'
            col_tag = args.value_column.replace('.', '_')
            ref_file = os.path.join(args.ref_dir, f"reference_lomaz_{args.group}_{args.gender}_chrom_{col_tag}.npz")
            logging.info(f"Auto-selected reference: {ref_file}")
        else:
            ref_file = args.ref_file

        # Compute chrom-level Z-scores dataframe
        zmat, chroms = compute_chrom_z_matrix(args.test_file, args.value_column, ref_file)
        # Plot chrom-level Z-score heatmap
        plot_chrom_zscore_heatmap(zmat, chroms, args)

    elif args.mode == 'analyze_chrom_old':
        if not args.ref_file:
            # Auto-select reference file if not provided, based on group and gender
            assert args.ref_dir and args.group and args.gender, 'Error: --ref_dir, --group, and --gender are required to auto-select reference.'
            col_tag = args.value_column.replace('.', '_')
            ref_file = os.path.join(args.ref_dir, f"reference_lomaz_{args.group}_{args.gender}_chrom_{col_tag}.npz")
            logging.info(f"Auto-selected reference: {ref_file}")
        else:
            ref_file = args.ref_file

        sample_basename = os.path.basename(args.test_file).split('.')[0]

        # Compute chromosome-level Z-score matrix
        zmat, chroms = compute_chrom_z_matrix(args.test_file, args.value_column, ref_file)
        #plt.figure(figsize=(9, 9))
        fig, ax = plt.subplots(figsize=(9, 9))
        # Plot heatmap of chromosome Z-scores
        sns.heatmap(zmat,
                    xticklabels=chroms,
                    yticklabels=chroms,
                    cmap='RdBu_r',
                    vmin=-5,
                    vmax=5,
                    center=0,
                    square=True,
                    linewidths=0.1,
                    linecolor='lightgray',
                    cbar=False,
                    #cbar_kws={'label': 'Z-score'},
                    ax=ax)
        # Add colorbar manually using a new axes
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.2)

        # Create the colorbar using the heatmap collection
        im = ax.collections[0]
        plt.colorbar(im, cax=cax, label="Z-score")

        #plt.title("LoMA-Z Chromosome Heatmap")
        ax.set_title(f"{sample_basename} [{args.group}, {args.value_column}]")

        plt.tight_layout()

        # Extract sample name from test_file
        sample_basename = os.path.basename(args.test_file).split('.')[0]

        # Create descriptive output filename
        output_name = f"{sample_basename}.chr_{args.group}_{args.value_column.replace('.', '_')}.png"
        output_path = os.path.join(args.output_dir, output_name)

        # Save plot
        plt.savefig(output_path)
        logging.info(f"Saved chromosome heatmap: {output_path}")

        # Plot mean Z-score barplot per chromosome
        plot_mean_z_per_chr(zmat, chroms, args.output.replace('.png', '_meanZ.png'))

    elif args.mode == 'analyze_bin':
        if not args.ref_file:
            # Auto-select reference file if not provided, based on group and gender
            assert args.ref_dir and args.group and args.gender, 'Error: --ref_dir, --group, and --gender are required to auto-select reference.'
            col_tag = args.value_column.replace('.', '_')
            ref_file = os.path.join(args.ref_dir, f"reference_lomaz_{args.group}_{args.gender}_bin_{col_tag}.npz")
            logging.info(f"Auto-selected reference: {ref_file}")
        else:
            ref_file = args.ref_file

        # Compute bin-level Z-scores dataframe
        z_df = compute_bin_z_matrix(args.test_file, args.value_column, ref_file)
        # Plot bin-level Z-score heatmap
        plot_bin_zscore_heatmap(z_df, args)

if __name__ == '__main__':
    main()
