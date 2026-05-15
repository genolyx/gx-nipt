#!/bin/bash
# run_v9_experiments.sh
# v9a: 170 samples, coverage-only  (fragment 영향 분리)
# v9b: 442 samples, coverage+frag  (샘플 수 영향 검증)
set -euo pipefail

REFS=/home/ken/gx-nipt/refs/gxff
PATCH=/home/ken/gx-nipt/bin/gxff_patch/pipeline.py
PATCH_DEST=/opt/conda/envs/nipt/lib/python3.11/site-packages/gxff/core/pipeline.py
GXFF=/opt/conda/envs/nipt/bin/python

echo "=== Preparing directories ==="
mkdir -p "$REFS/model_v9a" "$REFS/model_v9b"

# v9a: model_v8 config 재사용 (FF 이미 분율 스케일, 170 samples)
cp "$REFS/model_v8/training_config.tsv" "$REFS/model_v9a/training_config.tsv"
echo "[v9a] config copied (170 samples, coverage-only)"

# v9b: model_v7 config 복사 후 FF /100 변환 (442 samples)
sudo chmod 644 "$REFS/model_v7/training_config.tsv"
cp "$REFS/model_v7/training_config.tsv" "$REFS/model_v9b/training_config.tsv"

python3 - <<'PYEOF'
import pandas as pd
path = "/home/ken/gx-nipt/refs/gxff/model_v9b/training_config.tsv"
df = pd.read_csv(path, sep="\t")
if df["FF_REFERENCE"].max() > 1.0:
    df["FF_REFERENCE"] = (df["FF_REFERENCE"] / 100.0).round(6)
    df.to_csv(path, sep="\t", index=False)
    print(f"[v9b] FF fixed → range: {df.FF_REFERENCE.min():.4f}~{df.FF_REFERENCE.max():.4f}  n={len(df)}")
else:
    print(f"[v9b] already fraction scale  n={len(df)}")
PYEOF

echo ""
echo "=== Starting model_v9a (170 samples, coverage only) ==="
nohup docker run --rm \
  -v /home/ken:/home/ken \
  -v "$PATCH:$PATCH_DEST:ro" \
  -e HOME=/root \
  gx-nipt:latest \
  "$GXFF" -m gxff train \
    --config "$REFS/model_v9a/training_config.tsv" \
    --out    "$REFS/model_v9a" \
    --genome hg19 \
    --features coverage \
  > "$REFS/model_v9a.train.log" 2>&1 &
echo "v9a PID: $!"

echo ""
echo "=== Starting model_v9b (442 samples, coverage+fragment) ==="
nohup docker run --rm \
  -v /home/ken:/home/ken \
  -v "$PATCH:$PATCH_DEST:ro" \
  -e HOME=/root \
  gx-nipt:latest \
  "$GXFF" -m gxff train \
    --config "$REFS/model_v9b/training_config.tsv" \
    --out    "$REFS/model_v9b" \
    --genome hg19 \
    --features coverage fragment \
  > "$REFS/model_v9b.train.log" 2>&1 &
echo "v9b PID: $!"

echo ""
echo "=== Both started. Monitor with: ==="
echo "  tail -f $REFS/model_v9a.train.log"
echo "  tail -f $REFS/model_v9b.train.log"
