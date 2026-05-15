#!/bin/bash
# run_v10.sh
# 843 samples (male+female), coverage+fragment, 32 workers
# WIG path에서 BAM_PATH 자동 추가
set -euo pipefail

REFS=/home/ken/gx-nipt/refs/gxff
PATCH=/home/ken/gx-nipt/bin/gxff_patch/pipeline.py
PATCH_DEST=/opt/conda/envs/nipt/lib/python3.11/site-packages/gxff/core/pipeline.py
GXFF=/opt/conda/envs/nipt/bin/python
THREADS=32

echo "=== Step 1: model_v10 config 생성 (v6 config + BAM_PATH 추가) ==="
mkdir -p "$REFS/model_v10"

python3 - <<'PYEOF'
import os
import pandas as pd
from pathlib import Path

src = "/home/ken/gx-nipt/refs/gxff/training_config_v6.tsv"
dst = "/home/ken/gx-nipt/refs/gxff/model_v10/training_config.tsv"

df = pd.read_csv(src, sep="\t")
print(f"Loaded {len(df)} samples. FF range: {df.FF_REFERENCE.min():.4f}~{df.FF_REFERENCE.max():.4f}")

# WIG path 예시: /home/ken/ken-nipt/analysis/2603/GNCI26030020/Output_hmmcopy/GNCI26030020.of_fetus.50kb.wig.Normalization.txt
# BAM path:     /home/ken/ken-nipt/analysis/2603/GNCI26030020/GNCI26030020.proper_paired.bam
def wig_to_bam(wig_path, sample_id):
    # wig_path의 부모 부모 디렉토리 = sample 디렉토리
    sample_dir = str(Path(wig_path).parent.parent)
    return os.path.join(sample_dir, f"{sample_id}.proper_paired.bam")

bam_paths = []
missing_bam = 0
for _, row in df.iterrows():
    bam = wig_to_bam(row["FILEPATH"], row["SAMPLE_ID"])
    if os.path.exists(bam):
        bam_paths.append(bam)
    else:
        bam_paths.append("")  # missing: fragment features will be skipped
        missing_bam += 1

df["BAM_PATH"] = bam_paths
print(f"BAM found: {len(df) - missing_bam}/{len(df)}  missing: {missing_bam}")

df.to_csv(dst, sep="\t", index=False)
print(f"Saved: {dst}")
PYEOF

echo ""
echo "=== Step 2: gxff train 시작 (843 samples, coverage+fragment, ${THREADS} workers) ==="
nohup docker run --rm \
  -v /home/ken:/home/ken \
  -v "$PATCH:$PATCH_DEST:ro" \
  -e HOME=/root \
  gx-nipt:latest \
  "$GXFF" -m gxff train \
    --config  "$REFS/model_v10/training_config.tsv" \
    --out     "$REFS/model_v10" \
    --genome  hg19 \
    --features coverage fragment \
    --threads "$THREADS" \
  > "$REFS/model_v10.train.log" 2>&1 &

echo "Started PID: $!"
echo ""
echo "=== 로그 확인: tail -f $REFS/model_v10.train.log ==="
sleep 5
tail -5 "$REFS/model_v10.train.log"
