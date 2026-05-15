#!/bin/bash
# restart_model_v8.sh
# 1. 현재 model_v8 훈련 컨테이너 중단
# 2. training_config.tsv FF_REFERENCE를 분율 스케일로 수정 (% → /100)
# 3. 수정된 config로 gxff train 재시작
set -euo pipefail

LOG=/home/ken/gx-nipt/refs/gxff/model_v8.train.log
CONFIG=/home/ken/gx-nipt/refs/gxff/model_v8/training_config.tsv
OUTDIR=/home/ken/gx-nipt/refs/gxff/model_v8

echo "=== Step 1: 실행 중인 training 컨테이너 중단 ==="
# model_v8 관련 컨테이너 찾아서 중단
CONTAINER=$(docker ps --format "{{.Names}}" | grep -v daemon | head -1)
if [ -n "$CONTAINER" ]; then
    echo "Stopping container: $CONTAINER"
    docker stop "$CONTAINER"
else
    echo "No training container found (already stopped)"
fi

echo ""
echo "=== Step 2: FF_REFERENCE 분율 스케일 변환 (% → /100) ==="
python3 - <<'PYEOF'
import pandas as pd

config_path = "/home/ken/gx-nipt/refs/gxff/model_v8/training_config.tsv"
df = pd.read_csv(config_path, sep="\t")

ff_before = df["FF_REFERENCE"].describe()
print(f"[Before] FF_REFERENCE range: {df['FF_REFERENCE'].min():.4f} ~ {df['FF_REFERENCE'].max():.4f}")

# 1보다 크면 퍼센트 스케일이므로 /100 변환
if df["FF_REFERENCE"].max() > 1.0:
    df["FF_REFERENCE"] = (df["FF_REFERENCE"] / 100.0).round(6)
    print(f"[After]  FF_REFERENCE range: {df['FF_REFERENCE'].min():.6f} ~ {df['FF_REFERENCE'].max():.6f}")
    df.to_csv(config_path, sep="\t", index=False)
    print(f"Saved: {config_path}")
else:
    print("FF_REFERENCE is already in fraction scale, no change needed.")
PYEOF

echo ""
echo "=== Step 3: gxff train 재시작 (model_v8) ==="
nohup docker run --rm \
  -v /home/ken:/home/ken \
  -v /home/ken/gx-nipt/bin/gxff_patch/pipeline.py:/opt/conda/envs/nipt/lib/python3.11/site-packages/gxff/core/pipeline.py:ro \
  -e HOME=/root \
  gx-nipt:latest \
  /opt/conda/envs/nipt/bin/python -m gxff train \
    --config "$CONFIG" \
    --out    "$OUTDIR" \
    --genome hg19 \
    --features coverage fragment \
  >> "$LOG" 2>&1 &

echo "Started PID: $!"
echo ""
echo "=== 로그 확인: tail -f $LOG ==="
sleep 3
tail -5 "$LOG"
