"""
GxFFPipeline: Main orchestrator for feature extraction, model inference,
and model training.
"""

from __future__ import annotations

import logging
import os
import pickle
import time
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from gxff.core.config import GxFFConfig
from gxff.features.coverage import CoverageFeatureExtractor
from gxff.features.fragment import FragmentFeatureExtractor
from gxff.features.nucleosome import NucleosomeFeatureExtractor
from gxff.models.ensemble import GxFFEnsemble
from gxff.utils.io import load_bincount_file, load_training_config, save_tsv
from gxff.utils.augmentation import ArtificialSampleAugmentor
from gxff.validation.qc import QCChecker

logger = logging.getLogger(__name__)


class PredictionResult:
    """Container for a single-sample FF prediction result."""

    def __init__(
        self,
        sample_id: str,
        ff_gxff: float,
        ff_lgbm: float,
        ff_dnn: float,
        qc_flags: List[str],
        metadata: Optional[Dict] = None,
    ) -> None:
        self.sample_id = sample_id
        self.ff_gxff = ff_gxff
        self.ff_lgbm = ff_lgbm
        self.ff_dnn = ff_dnn
        self.qc_flags = qc_flags
        self.metadata = metadata or {}

    def to_dict(self) -> Dict:
        return {
            "SAMPLE_ID": self.sample_id,
            "FF_GXFF": round(self.ff_gxff, 4),
            "FF_LGBM": round(self.ff_lgbm, 4),
            "FF_DNN": round(self.ff_dnn, 4),
            "QC_FLAGS": ";".join(self.qc_flags) if self.qc_flags else "PASS",
            **self.metadata,
        }

    def __repr__(self) -> str:
        return (
            f"PredictionResult(sample_id={self.sample_id!r}, "
            f"ff_gxff={self.ff_gxff:.4f}, qc={self.qc_flags})"
        )


class GxFFPipeline:
    """
    End-to-end pipeline for gx-FF fetal fraction estimation.

    Supports:
    - Single-sample prediction from BAM/CRAM or pre-computed bin count file
    - Batch prediction from a sample list
    - Model training from a labeled training config

    Parameters
    ----------
    config : GxFFConfig
        Pipeline configuration.
    model_path : str or Path, optional
        Path to a pre-trained model file (.pkl). Required for prediction.
    """

    def __init__(
        self,
        config: Optional[GxFFConfig] = None,
        model_path: Optional[str | Path] = None,
    ) -> None:
        self.config = config or GxFFConfig()
        self.model: Optional[GxFFEnsemble] = None

        if model_path is not None:
            self._load_model(model_path)

        # Feature extractors (lazily initialized)
        self._cov_extractor: Optional[CoverageFeatureExtractor] = None
        self._frag_extractor: Optional[FragmentFeatureExtractor] = None
        self._nuc_extractor: Optional[NucleosomeFeatureExtractor] = None

    # ── Model I/O ─────────────────────────────────────────────────────────────

    def _load_model(self, path: str | Path) -> None:
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Model file not found: {path}")
        logger.info("Loading model from %s", path)
        with open(path, "rb") as f:
            self.model = pickle.load(f)
        logger.info("Model loaded successfully.")

    def _save_model(self, path: str | Path) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "wb") as f:
            pickle.dump(self.model, f)
        logger.info("Model saved to %s", path)

    # ── Feature extractors ────────────────────────────────────────────────────

    @property
    def cov_extractor(self) -> CoverageFeatureExtractor:
        if self._cov_extractor is None:
            self._cov_extractor = CoverageFeatureExtractor(self.config)
        return self._cov_extractor

    @property
    def frag_extractor(self) -> FragmentFeatureExtractor:
        if self._frag_extractor is None:
            self._frag_extractor = FragmentFeatureExtractor(self.config)
        return self._frag_extractor

    @property
    def nuc_extractor(self) -> NucleosomeFeatureExtractor:
        if self._nuc_extractor is None:
            self._nuc_extractor = NucleosomeFeatureExtractor(self.config)
        return self._nuc_extractor

    # ── Feature assembly ──────────────────────────────────────────────────────

    def _extract_features(
        self,
        bam_path: Optional[str] = None,
        bincount_path: Optional[str] = None,
        ref_path: Optional[str] = None,
        bmi: Optional[float] = None,
        ga_weeks: Optional[int] = None,
    ) -> np.ndarray:
        """
        Extract and concatenate all configured feature types into a 1-D
        feature vector for a single sample.
        """
        feature_parts = []

        if self.config.use_coverage:
            logger.debug("Extracting coverage features...")
            if bincount_path:
                bincount_df = load_bincount_file(bincount_path)
                cov_feat = self.cov_extractor.from_dataframe(bincount_df)
            elif bam_path:
                cov_feat = self.cov_extractor.from_bam(bam_path, ref_path=ref_path)
            else:
                raise ValueError("Either bam_path or bincount_path must be provided.")
            feature_parts.append(cov_feat)

        if self.config.use_fragment:
            logger.debug("Extracting fragment length features...")
            if bam_path:
                frag_feat = self.frag_extractor.from_bam(bam_path, ref_path=ref_path)
            else:
                logger.warning(
                    "Fragment features require BAM input; skipping for bincount input."
                )
                frag_feat = self.frag_extractor.empty_features()
            feature_parts.append(frag_feat)

        if self.config.use_nucleosome:
            logger.debug("Extracting nucleosome profiling features...")
            if bam_path:
                nuc_feat = self.nuc_extractor.from_bam(bam_path, ref_path=ref_path)
            else:
                logger.warning(
                    "Nucleosome features require BAM input; skipping for bincount input."
                )
                nuc_feat = self.nuc_extractor.empty_features()
            feature_parts.append(nuc_feat)

        # Append clinical covariates if available
        cov_vec = []
        cov_vec.append(float(bmi) if bmi is not None else np.nan)
        cov_vec.append(float(ga_weeks) if ga_weeks is not None else np.nan)
        feature_parts.append(np.array(cov_vec, dtype=np.float64))

        return np.concatenate(feature_parts)

    # ── Prediction ────────────────────────────────────────────────────────────

    def predict_single(
        self,
        bam_path: str,
        ref_path: Optional[str] = None,
        bmi: Optional[float] = None,
        ga_weeks: Optional[int] = None,
    ) -> PredictionResult:
        """Predict FF for a single BAM/CRAM file."""
        if self.model is None:
            raise RuntimeError("No model loaded. Call _load_model() first.")

        sample_id = Path(bam_path).stem
        logger.info("Predicting FF for sample: %s", sample_id)
        t0 = time.time()

        features = self._extract_features(
            bam_path=bam_path,
            ref_path=ref_path,
            bmi=bmi,
            ga_weeks=ga_weeks,
        )

        ff_gxff_arr, ff_lgbm_arr, ff_dnn_arr = self.model.predict(features.reshape(1, -1))
        ff_gxff = float(np.squeeze(ff_gxff_arr))
        ff_lgbm = float(np.squeeze(ff_lgbm_arr))
        ff_dnn = float(np.squeeze(ff_dnn_arr))
        qc_flags = QCChecker.check(ff_gxff, features)

        elapsed = time.time() - t0
        logger.info(
            "Sample %s: FF_GXFF=%.4f (lgbm=%.4f, dnn=%.4f) | QC=%s | %.1fs",
            sample_id, ff_gxff, ff_lgbm, ff_dnn,
            qc_flags if qc_flags else "PASS", elapsed,
        )

        return PredictionResult(
            sample_id=sample_id,
            ff_gxff=ff_gxff,
            ff_lgbm=ff_lgbm,
            ff_dnn=ff_dnn,
            qc_flags=qc_flags,
            metadata={"elapsed_sec": round(elapsed, 2)},
        )

    def predict_from_bincount(
        self,
        bincount_path: str,
        bmi: Optional[float] = None,
        ga_weeks: Optional[int] = None,
    ) -> PredictionResult:
        """Predict FF from a pre-computed bin count file."""
        if self.model is None:
            raise RuntimeError("No model loaded. Call _load_model() first.")

        sample_id = Path(bincount_path).stem
        logger.info("Predicting FF from bincount: %s", sample_id)
        t0 = time.time()

        features = self._extract_features(
            bincount_path=bincount_path,
            bmi=bmi,
            ga_weeks=ga_weeks,
        )

        ff_gxff_arr, ff_lgbm_arr, ff_dnn_arr = self.model.predict(features.reshape(1, -1))
        ff_gxff = float(np.squeeze(ff_gxff_arr))
        ff_lgbm = float(np.squeeze(ff_lgbm_arr))
        ff_dnn = float(np.squeeze(ff_dnn_arr))
        qc_flags = QCChecker.check(ff_gxff, features)

        elapsed = time.time() - t0
        logger.info(
            "Sample %s: FF_GXFF=%.4f | QC=%s | %.1fs",
            sample_id, ff_gxff,
            qc_flags if qc_flags else "PASS", elapsed,
        )

        return PredictionResult(
            sample_id=sample_id,
            ff_gxff=ff_gxff,
            ff_lgbm=ff_lgbm,
            ff_dnn=ff_dnn,
            qc_flags=qc_flags,
            metadata={"elapsed_sec": round(elapsed, 2)},
        )

    def predict_batch(
        self,
        sample_list_path: str,
        ref_path: Optional[str] = None,
    ) -> List[PredictionResult]:
        """Predict FF for all samples listed in a text file."""
        sample_list = Path(sample_list_path).read_text().strip().splitlines()
        sample_list = [s.strip() for s in sample_list if s.strip()]
        logger.info("Batch prediction: %d samples", len(sample_list))

        results = []
        for i, bam_path in enumerate(sample_list, 1):
            logger.info("Processing sample %d/%d: %s", i, len(sample_list), bam_path)
            try:
                result = self.predict_single(bam_path=bam_path, ref_path=ref_path)
                results.append(result)
            except Exception as exc:
                logger.error("Failed to process %s: %s", bam_path, exc)

        return results

    # ── Training ──────────────────────────────────────────────────────────────

    def train(self, config_path: str, output_dir: str) -> None:
        """
        Train the gx-FF ensemble model from a labeled training configuration.

        Parameters
        ----------
        config_path : str
            Path to TSV training config (SAMPLE_ID, FILEPATH, SEX_FETUS,
            FF_REFERENCE, [BMI, GA_WEEKS]).
        output_dir : str
            Directory to write the trained model and evaluation outputs.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        logger.info("Loading training configuration from %s", config_path)
        train_df = load_training_config(config_path)
        n_samples = len(train_df)
        logger.info("Training samples: %d", n_samples)

        # ── Feature extraction ────────────────────────────────────────────────
        logger.info("Extracting features for all training samples...")
        X_list, y_list, weights_list = [], [], []

        for _, row in train_df.iterrows():
            sample_id = row["SAMPLE_ID"]
            filepath = row["FILEPATH"]
            ff_ref = row.get("FF_REFERENCE", np.nan)
            bmi = row.get("BMI", None)
            ga_weeks = row.get("GA_WEEKS", None)
            # Optional BAM_PATH column: when FILEPATH is a bincount/wig file,
            # BAM_PATH provides the corresponding BAM for fragment feature extraction.
            bam_for_frag = row.get("BAM_PATH", None)
            if bam_for_frag is not None and (pd.isna(bam_for_frag) or str(bam_for_frag).strip() == ""):
                bam_for_frag = None

            # Skip samples without reference FF (e.g., female fetuses with no SNP-FF)
            if pd.isna(ff_ref):
                logger.debug("Skipping %s: no reference FF", sample_id)
                continue

            try:
                is_bam = filepath.endswith((".bam", ".cram"))
                feat = self._extract_features(
                    bam_path=bam_for_frag if (bam_for_frag and not is_bam) else (filepath if is_bam else None),
                    bincount_path=filepath if not is_bam else None,
                    bmi=float(bmi) if bmi and not pd.isna(bmi) else None,
                    ga_weeks=int(ga_weeks) if ga_weeks and not pd.isna(ga_weeks) else None,
                )
                X_list.append(feat)
                y_list.append(float(ff_ref))

                # Apply higher weight to low-FF samples
                w = self.config.low_ff_weight if float(ff_ref) < 0.05 else 1.0
                weights_list.append(w)

            except Exception as exc:
                logger.error("Feature extraction failed for %s: %s", sample_id, exc)

        if len(X_list) == 0:
            raise RuntimeError("No valid training samples with reference FF found.")

        X = np.vstack(X_list)
        y = np.array(y_list, dtype=np.float64)
        weights = np.array(weights_list, dtype=np.float64)
        logger.info("Feature matrix shape: %s | FF range: [%.3f, %.3f]",
                    X.shape, y.min(), y.max())

        # ── Artificial sample augmentation ────────────────────────────────────
        if self.config.augment:
            logger.info("Applying artificial sample augmentation for low-FF range...")
            augmentor = ArtificialSampleAugmentor(
                low_ff_threshold=0.05,
                target_low_ff_ratio=0.25,
            )
            X, y, weights = augmentor.augment(X, y, weights)
            logger.info("After augmentation: %d samples", len(y))

        # ── Model training ────────────────────────────────────────────────────
        logger.info("Training ensemble model...")
        self.model = GxFFEnsemble(
            n_pca_components=self.config.n_pca_components,
            cv_folds=self.config.cv_folds,
            n_jobs=self.config.threads,
        )
        cv_metrics = self.model.fit(X, y, sample_weight=weights)

        logger.info(
            "Cross-validation results: Pearson r=%.4f, MAE=%.4f, RMSE=%.4f",
            cv_metrics["pearson_r"], cv_metrics["mae"], cv_metrics["rmse"],
        )

        # ── Save outputs ──────────────────────────────────────────────────────
        model_path = output_dir / "gxff_model.pkl"
        self._save_model(model_path)

        metrics_path = output_dir / "cv_metrics.tsv"
        save_tsv(pd.DataFrame([cv_metrics]), metrics_path)
        logger.info("Training complete. Model saved to %s", model_path)

    # ── Output ────────────────────────────────────────────────────────────────

    def save_results(
        self, results: List[PredictionResult], output_dir: str
    ) -> None:
        """Write prediction results to a TSV file."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        records = [r.to_dict() for r in results]
        df = pd.DataFrame(records)
        out_path = output_dir / "ff_predictions.tsv"
        save_tsv(df, out_path)
        logger.info("Predictions written to %s", out_path)

        # Print summary to stdout
        print(df.to_string(index=False))
