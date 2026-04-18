# NIPT Reference Creation Pipeline

통합 reference 생성 파이프라인 - 단일 스크립트로 모든 reference 생성 가능

## Overview

이 스크립트는 NIPT 분석을 위한 4가지 reference를 생성합니다:

1. **EZD** - UAR & Z-score 기반 threshold
2. **PRIZM** - Mean & SD 기반 Z-score  
3. **WC** - Wisecondor npz reference
4. **WCX** - WisecondorX npz reference

각 reference는 **orig**, **fetus**, **mom** 세 그룹으로 분류되어 생성됩니다.

## Prerequisites

### Required Tools
- Python 3.6+
- WisecondorX (`/usr/bin/miniconda3/bin/WisecondorX`)
- Wisecondor (`/Work/NIPT/bin/wisecondor/wisecondor.py`)
- Python 2.7 (Wisecondor 실행용)

### Required Python Packages
```bash
pip install pandas numpy
```

### Input Data Requirements

1. **샘플 리스트 파일** (TSV 형식)
   - 엑셀에서 추출한 샘플 정보
   - 필수 컬럼:
     - `sample_id`: 샘플 ID
     - `fetal_gender(gd_2)`: 태아 성별 (XY 또는 XX)
     - `Result`: EZD 결과 (Low Risk, High Risk, No Call 등)
     - `MDResult`: MD 결과  
     - `SeqFF`: Fetal fraction
     - `mapping_rate(%)`: Mapping rate
   
2. **분석 완료된 샘플들**
   - EZD 결과 파일: `analysis/YYMM/SAMPLE_ID/Output_EZD/GROUP/GROUP_ezd_results.tsv`
   - HMMcopy 결과: `analysis/YYMM/SAMPLE_ID/Output_hmmcopy/SAMPLE_ID.of_GROUP.10mb.wig.Normalization.txt`
   - NPZ 파일: `analysis/YYMM/SAMPLE_ID/.../SAMPLE_ID*.npz`

## Usage

### 1. 샘플 필터링 미리보기

실제 reference를 생성하기 전에 어떤 샘플들이 선택되는지 확인:

```bash
python create_reference.py \
  --sample-list reference_sample_list.tsv \
  --preview-only
```

### 2. 모든 Reference 생성

```bash
python create_reference.py \
  --sample-list reference_sample_list.tsv \
  --labcode ucl \
  --ref-type all
```

### 3. 특정 Reference만 생성

#### EZD만 생성
```bash
python create_reference.py \
  --sample-list reference_sample_list.tsv \
  --labcode ucl \
  --ref-type ezd
```

#### PRIZM만 생성
```bash
python create_reference.py \
  --sample-list reference_sample_list.tsv \
  --labcode ucl \
  --ref-type prizm
```

#### WC만 생성
```bash
python create_reference.py \
  --sample-list reference_sample_list.tsv \
  --labcode ucl \
  --ref-type wc
```

#### WCX만 생성
```bash
python create_reference.py \
  --sample-list reference_sample_list.tsv \
  --labcode ucl \
  --ref-type wcx
```

### 4. 특정 Group만 생성

```bash
python create_reference.py \
  --sample-list reference_sample_list.tsv \
  --labcode ucl \
  --ref-type all \
  --groups orig fetus
```

### 5. SeqFF 범위 지정

```bash
python create_reference.py \
  --sample-list reference_sample_list.tsv \
  --labcode ucl \
  --ref-type all \
  --min-seqff 5.0 \
  --max-seqff 25.0
```

### 6. 출력 디렉토리 지정

```bash
python create_reference.py \
  --sample-list reference_sample_list.tsv \
  --labcode cordlife \
  --ref-type all \
  --output-dir /Work/NIPT/data/refs/cordlife_new
```

## Sample Filtering Logic

스크립트는 다음 기준으로 샘플을 자동 필터링합니다:

### 제외 기준
1. `Result` 컬럼이 "High Risk" 또는 "No Call"인 샘플 **(대소문자 무관)**
2. `MDResult` 컬럼이 "High Risk" 또는 "No Call"인 샘플 **(대소문자 무관)**
3. `SeqFF` 값이 범위 밖인 샘플 (기본: 4.0-30.0%)
4. `mapping_rate` < 95%인 샘플

**참고:** "High Risk", "high risk", "High risk" 등 대소문자가 다르게 표기된 경우도 모두 감지하여 제외합니다.

### 분류 기준
- **Male (M)**: `fetal_gender(gd_2)` == "XY"
- **Female (F)**: `fetal_gender(gd_2)` == "XX"

## Output Structure

```
/Work/NIPT/data/refs/LABCODE/
├── EZD/
│   ├── orig/
│   │   ├── chr1.txt ~ chr22.txt
│   │   ├── male.txt, female.txt
│   │   └── orig_thresholds_new.tsv
│   ├── fetus/
│   │   └── ...
│   └── mom/
│       └── ...
├── PRIZM/
│   ├── orig/
│   │   ├── total_mean.csv, total_sd.csv
│   │   ├── male_mean.csv, male_sd.csv
│   │   ├── female_mean.csv, female_sd.csv
│   │   ├── total_10mb_mean.csv, total_10mb_sd.csv
│   │   └── ...
│   ├── fetus/
│   └── mom/
├── WC/
│   ├── orig_200k_of.npz
│   ├── fetus_200k_of.npz
│   └── mom_200k_of.npz
└── WCX/
    ├── orig_M_200k_of.npz
    ├── orig_F_200k_of.npz
    ├── fetus_M_200k_of.npz
    ├── fetus_F_200k_of.npz
    └── mom_200k_of.npz
```

## Command Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--sample-list` | 샘플 리스트 TSV 파일 (필수) | - |
| `--labcode` | Laboratory code | ucl |
| `--ref-type` | Reference 타입 (all/ezd/prizm/wc/wcx) | all |
| `--groups` | 처리할 그룹 (orig/fetus/mom) | all |
| `--output-dir` | 출력 디렉토리 | DATA_DIR/refs/LABCODE |
| `--min-seqff` | 최소 SeqFF 값 | 4.0 |
| `--max-seqff` | 최대 SeqFF 값 | 30.0 |
| `--preview-only` | 필터링 결과만 미리보기 | False |
| `--debug` | 디버그 로깅 활성화 | False |

## Reference Types

### 1. EZD Reference

EZD (Enhanced Z-score Detection) reference는 각 염색체별로:
- UAR (Unique Alignment Read) 값
- Z-score 값

을 기반으로 threshold를 설정합니다.

**생성 파일:**
- `chr1.txt ~ chr22.txt`: 각 염색체별 샘플 데이터
- `male.txt, female.txt`: 성염색체 데이터
- `orig_thresholds_new.tsv`: 계산된 threshold

### 2. PRIZM Reference

PRIZM (Prenatal Risk Z-score Matrix) reference는:
- Chromosome level (24x24 matrix)
- 10mb bin level
- 10mb_all level (22x322 matrix)

각 레벨에서 mean과 standard deviation을 계산합니다.

**생성 파일:**
- `{gender}_mean.csv`, `{gender}_sd.csv`
- `{gender}_10mb_mean.csv`, `{gender}_10mb_sd.csv`
- `{gender}_10mb_all_mean.csv`, `{gender}_10mb_all_sd.csv`

### 3. WC (Wisecondor) Reference

Original Wisecondor를 위한 npz reference 파일을 생성합니다.

**생성 파일:**
- `orig_200k_of.npz`
- `fetus_200k_of.npz`
- `mom_200k_of.npz`

### 4. WCX (WisecondorX) Reference

WisecondorX를 위한 npz reference 파일을 생성합니다.

**생성 파일:**
- `orig_M_200k_of.npz`, `orig_F_200k_of.npz`
- `fetus_M_200k_of.npz`, `fetus_F_200k_of.npz`
- `mom_200k_of.npz` (female only)

## Example Workflow

### Step 1: 샘플 리스트 준비

엑셀에서 샘플 정보를 추출하여 TSV로 저장:

```bash
# Excel → TSV 변환 (Excel에서 "다른 이름으로 저장" → "텍스트(탭으로 분리)")
```

### Step 2: 필터링 미리보기

```bash
python create_reference.py \
  --sample-list reference_sample_list_GNCI.tsv \
  --preview-only
```

출력 예시:
```
========== Sample Filtering ==========
Initial sample count: 500
Exclude Result='High Risk': 500 -> 485
Exclude Result='No Call': 485 -> 480
Exclude MDResult='High Risk': 480 -> 475
SeqFF range (4.0-30.0%): 475 -> 450
Mapping rate >= 95%: 450 -> 445

Gender distribution:
  Male (XY): 225 samples
  Female (XX): 220 samples

Filtering completed: 500 -> 445 samples
```

### Step 3: Reference 생성

```bash
python create_reference.py \
  --sample-list reference_sample_list_GNCI.tsv \
  --labcode ucl \
  --ref-type all \
  --groups orig fetus mom
```

### Step 4: 생성 확인

```bash
# EZD 확인
ls -lh /Work/NIPT/data/refs/ucl/EZD/orig/

# PRIZM 확인
ls -lh /Work/NIPT/data/refs/ucl/PRIZM/orig/

# WC/WCX 확인
ls -lh /Work/NIPT/data/refs/ucl/WC/
ls -lh /Work/NIPT/data/refs/ucl/WCX/
```

## Troubleshooting

### 문제: NPZ 파일을 찾을 수 없음

**원인:** 분석이 완료되지 않았거나 파일 경로가 다름

**해결:**
1. 해당 샘플의 분석이 완료되었는지 확인
2. NPZ 파일이 실제로 존재하는지 확인:
   ```bash
   find /Work/NIPT/analysis -name "SAMPLE_ID*.npz"
   ```

### 문제: Count 파일을 찾을 수 없음

**원인:** HMMcopy 결과가 없음

**해결:**
1. HMMcopy 분석이 완료되었는지 확인
2. Normalization.txt 파일 확인:
   ```bash
   find /Work/NIPT/analysis -name "*Normalization.txt"
   ```

### 문제: Memory Error

**원인:** PRIZM reference 생성 시 메모리 부족

**해결:**
1. 샘플 수를 줄임 (필터링 기준 강화)
2. 그룹별로 나누어 실행

### 문제: Permission Denied

**원인:** 출력 디렉토리 쓰기 권한 없음

**해결:**
```bash
sudo chown -R $USER:$USER /Work/NIPT/data/refs/
```

## Notes

### 성능 고려사항

- **EZD**: 빠름 (수 분)
- **PRIZM**: 느림 (수십 분 ~ 수 시간, 샘플 수에 따라)
- **WC/WCX**: 중간 (수 분 ~ 수십 분)

### 권장 샘플 수

- **최소**: 각 그룹당 50개 이상 (M/F 각각)
- **권장**: 각 그룹당 100-200개 (M/F 각각)
- **이상적**: 각 그룹당 300개 이상 (M/F 각각)

### SeqFF 범위 권장

- **일반**: 4.0% - 30.0%
- **엄격**: 5.0% - 25.0%
- **느슨**: 3.0% - 35.0%

## Advanced Usage

### 커스텀 필터링

스크립트를 수정하여 추가 필터링 기준을 적용할 수 있습니다:

```python
# create_reference.py의 SampleSelector.filter_samples() 메서드 수정

# 예: GC content 필터링 추가
if 'GC_content(%)' in self.df.columns:
    before = len(self.df)
    self.df = self.df[(self.df['GC_content(%)'] >= 40) & (self.df['GC_content(%)'] <= 45)]
    logger.info(f"GC content range (40-45%): {before} -> {len(self.df)}")
```

### 배치 실행

여러 labcode에 대해 한 번에 실행:

```bash
#!/bin/bash
for labcode in ucl cordlife vn testclient; do
    echo "Creating reference for $labcode..."
    python create_reference.py \
      --sample-list reference_sample_list_${labcode}.tsv \
      --labcode $labcode \
      --ref-type all
done
```

## Support

문제가 발생하면:
1. `--debug` 옵션으로 자세한 로그 확인
2. 필터링된 샘플 리스트 (`*_filtered.tsv`) 확인
3. 분석 결과 파일들의 존재 여부 확인

## Changelog

### Version 1.0 (2025-01-06)
- 초기 버전 릴리즈
- EZD, PRIZM, WC, WCX 지원
- orig, fetus, mom 그룹 지원
- 자동 샘플 필터링
- Gender 기반 분류

