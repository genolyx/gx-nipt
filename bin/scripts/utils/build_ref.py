# ---------------------------------------------
# Build reference medians and MADs for chrom/bin
# ---------------------------------------------
def build_reference_matrix(ref_dir, value_column, mode, output_dir):
    """
    Construct reference matrices (median and MAD) from multiple reference samples.
    Includes chrY handling for female samples and MAD stabilization.
    """
    groups = ['orig', 'fetus', 'mom']
    sexes = ['M', 'F']

    for group in groups:
        for sex in sexes:
            group_dir = os.path.join(ref_dir, group, sex)
            logging.debug(f"Searching in: {group_dir}")
            if not os.path.exists(group_dir):
                logging.warning(f"Skipping missing group: {group_dir}")
                continue

            files = glob.glob(os.path.join(group_dir, '*.10mb.txt'))
            if not files:
                logging.warning(f"No files in group: {group_dir}")
                continue

            logging.debug(f"Found files: {files}")


            # ---------------------------------------------
            # Chromosome-wise
            # ---------------------------------------------
            if mode == 'chrom':
                for group in groups:
                    # Gather all M and F files for chr1–chr22
                    all_files = []
                    for sex_dir in ['M', 'F']:
                        group_dir = os.path.join(ref_dir, group, sex_dir)
                        if os.path.exists(group_dir):
                            all_files += glob.glob(os.path.join(group_dir, '*.10mb.txt'))
                    if not all_files:
                        logging.warning(f"No files found for group: {group}")
                        continue

                    # Compute shared autosome reference med/mad
                    shared_data = []
                    for f in all_files:
                        df = parse_normalization(f, value_column)
                        counts = summarize_chrom_level(df, value_column)
                        ratios = log_ratios(counts)
                        shared_data.append(ratios.values)
                    shared_data = np.stack(shared_data)

                    # Now compute per-sex med/mad with sex-specific chrX/chrY
                    for sex in sexes:
                        # Get sex-specific files (for chrX/chrY only)
                        sex_dir = os.path.join(ref_dir, group, sex)
                        sex_files = glob.glob(os.path.join(sex_dir, '*.10mb.txt')) if os.path.exists(sex_dir) else []
                        if not sex_files:
                            logging.warning(f"No sex-specific files for {group}/{sex}, skipping.")
                            continue

                        sex_data = []
                        for f in sex_files:
                            df = parse_normalization(f, value_column)
                            counts = summarize_chrom_level(df, value_column)
                            ratios = log_ratios(counts)
                            sex_data.append(ratios.values)
                        sex_data = np.stack(sex_data)

                        # Initialize full matrix from shared autosome med/mad
                        med_mat = np.zeros((24, 24))
                        mad_mat = np.zeros((24, 24))
                        for i in range(24):
                            for j in range(24):
                                if i == j:
                                    continue
                                # chr1–22: use shared data
                                if i < 22 and j < 22:
                                    diffs = shared_data[:, i] - shared_data[:, j]
                                else:
                                    diffs = sex_data[:, i] - sex_data[:, j]

                                med = np.median(diffs)
                                if sex == 'F' and (i == 23 or j == 23):
                                    # Female: chrY is irrelevant
                                    med = 0.0
                                    mad = 1.0
                                else:
                                    mad = np.median(np.abs(diffs - med))
                                    if sex == 'M' and (i == 23 or j == 23):
                                        mad = np.clip(mad, 0.05, None)
                                    else:
                                        mad = mad + 0.05 * abs(med)

                                med_mat[i, j] = med
                                mad_mat[i, j] = mad

                        col_tag = value_column.replace('.', '_')
                        out_file = os.path.join(output_dir, f'reference_lomaz_{group}_{sex}_chrom_{col_tag}.npz')
                        np.savez(out_file, median=med_mat, mad=mad_mat, chroms=CHROMS)
                        plot_mad_heatmap(mad_mat, CHROMS, out_file.replace('.npz', '_mad.png'))
                        logging.info(f"Saved full reference for {group}/{sex}: {out_file}")

            # ---------------------------------------------
            # 10mb bin-wise
            # ---------------------------------------------
            elif mode == 'bin':
                for group in groups:
                    # Step 1: Collect all shared files for chr1–22
                    shared_files = []
                    for sex_dir in ['M', 'F']:
                        group_dir = os.path.join(ref_dir, group, sex_dir)
                        if os.path.exists(group_dir):
                            shared_files += glob.glob(os.path.join(group_dir, '*.10mb.txt'))
                    if not shared_files:
                        logging.warning(f"No shared bin-level files for group: {group}")
                        continue

                    # Load all shared files for chr1–chr22 bins
                    all_ratios = []
                    for f in shared_files:
                        df = parse_normalization(f, value_column)
                        df = df.sort_values(['chr', 'start'])
                        df['log_ratio'] = log_ratios(df[value_column])
                        all_ratios.append(df['log_ratio'].values)
                        if 'bin_ids' not in locals():
                            chr_ids = df['chr'].tolist()
                            bin_ids = [f"{row['chr']}:{row['start']//1000000}Mb" for _, row in df.iterrows()]
                    all_ratios = np.stack(all_ratios)

                    # Step 2: Process per-sex med/MAD with chrX/Y from sex-specific files
                    for sex in sexes:
                        sex_dir = os.path.join(ref_dir, group, sex)
                        sex_files = glob.glob(os.path.join(sex_dir, '*.10mb.txt')) if os.path.exists(sex_dir) else []
                        if not sex_files:
                            logging.warning(f"No sex-specific bin files for {group}/{sex}")
                            continue

                        sex_ratios = []
                        for f in sex_files:
                            df = parse_normalization(f, value_column)
                            df = df.sort_values(['chr', 'start'])
                            df['log_ratio'] = log_ratios(df[value_column])
                            sex_ratios.append(df['log_ratio'].values)
                        sex_ratios = np.stack(sex_ratios)

                        # Now compute med/mad vector for each bin
                        med = np.zeros(all_ratios.shape[1])
                        mad = np.zeros(all_ratios.shape[1])
                        for idx in range(len(bin_ids)):
                            chrom = chr_ids[idx]
                            if chrom in CHROMS[:22]:  # chr1~22 → use all_ratios (M+F)
                                values = all_ratios[:, idx]
                            else:  # chrX or chrY → use sex-specific values
                                values = sex_ratios[:, idx]

                            med[idx] = np.median(values)
                            if sex == 'F' and chrom == 'chrY':
                                mad[idx] = 1.0
                                med[idx] = 0.0
                            else:
                                mad_val = np.median(np.abs(values - med[idx]))
                                if sex == 'M' and chrom == 'chrY':
                                    mad[idx] = np.clip(mad_val, 0.05, None)
                                else:
                                    mad[idx] = mad_val + 0.05 * abs(med[idx])

                        # Save reference
                        col_tag = value_column.replace('.', '_')
                        out_file = os.path.join(output_dir, f'reference_lomaz_{group}_{sex}_bin_{col_tag}.npz')
                        np.savez(out_file, median=med, mad=mad, bin_ids=np.array(bin_ids), chr_ids=np.array(chr_ids))
                        plot_mad_heatmap(np.reshape(mad[:24], (1, -1)), CHROMS[:24], out_file.replace('.npz', '_mad.png'))
                        logging.info(f"Saved bin-level reference for {group}/{sex}: {out_file}")

