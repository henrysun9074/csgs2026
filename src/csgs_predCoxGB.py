#!/usr/bin/env python3
"""
Runs global hyperparameter tuning once using Pearson correlation rather than MSE
"""

import os
import json
import logging
import argparse
import numpy as np
import pandas as pd
from collections import defaultdict
from scipy.stats import pearsonr

from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    make_scorer,
    mean_squared_error,
    mean_absolute_error,
    r2_score,
)
from sklearn.linear_model import Ridge, ElasticNet
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor

from skopt import BayesSearchCV
from skopt.space import Real, Integer, Categorical


parser = argparse.ArgumentParser(
    description="Regression GS Pipeline for phenotype_gwas"
)
parser.add_argument("--genofile", "-g", type=str, required=True)
parser.add_argument("--phenofile", "-p", type=str, required=True)
parser.add_argument("--outdir", "-o", type=str, required=True)
parser.add_argument("--hpt_iter", type=int, default=50)
parser.add_argument("--verbose", "-v", action="store_true")
args = parser.parse_args()

logging.basicConfig(
    level=logging.DEBUG if args.verbose else logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
)
logger = logging.getLogger(__name__)

BASE_OUT_DIR = os.path.join("/work/hs325/csgs2026", args.outdir)
os.makedirs(BASE_OUT_DIR, exist_ok=True)

def get_search_spaces():
    return {
        "GB": (
            XGBRegressor(
                tree_method="hist",
                device="cpu",
                objective="reg:squarederror",
                eval_metric="rmse",
                random_state=123            
                ),
            {
                "n_estimators": Integer(100, 2000),
                "max_depth": Integer(3, 15),
                "learning_rate": Real(1e-3, 0.3, prior="log-uniform"),
                "subsample": Real(0.5, 1.0),
                "colsample_bytree": Real(0.5, 1.0),
                "min_child_weight": Integer(1, 10),
                "gamma": Real(0, 5),
            },
        ),
    }

# Pearson correlation scoring function
def pearson_corr_func(y_true, y_pred):
    if y_pred.ndim == 2:  
        y_pred = y_pred[:, 1] 
    if np.allclose(y_pred, y_pred[0]) or len(np.unique(y_true)) < 2:
        return 0.0
    corr, _ = pearsonr(y_pred, y_true)
    return corr if not np.isnan(corr) else 0.0

# Define Pearson scorer for optimization
pearson_scorer = make_scorer(pearson_corr_func, greater_is_better=True)

def load_and_align_data(geno_path, pheno_path):
    logger.info("Loading genotype and phenotype files...")

    geno_df = pd.read_csv(geno_path)
    pheno_df = pd.read_csv(pheno_path)

    if "ID" not in geno_df.columns:
        raise ValueError("Genotype file must contain an 'ID' column.")
    if "SampleID" not in pheno_df.columns:
        raise ValueError("Phenotype file must contain a 'SampleID' column.")
    if "phenotype_gwas" not in pheno_df.columns:
        raise ValueError("Phenotype file must contain a 'phenotype_gwas' column.")

    geno_temp = geno_df.set_index("ID")
    pheno_temp = pheno_df.set_index("SampleID")

    common_ids = geno_temp.index.intersection(pheno_temp.index)
    logger.info(f"Intersecting samples discovered: {len(common_ids)}")

    if len(common_ids) == 0:
        raise ValueError(
            "No matching entries found between Genotype 'ID' and Phenotype 'SampleID'."
        )

    common_ids_sorted = sorted(list(common_ids))
    geno_temp = geno_temp.loc[common_ids_sorted]
    pheno_temp = pheno_temp.loc[common_ids_sorted]

    ax_cols = [col for col in geno_temp.columns if col.startswith("AX")]
    if len(ax_cols) == 0:
        raise ValueError("No genotype marker columns found starting with 'AX'.")

    y = pd.to_numeric(pheno_temp["phenotype_gwas"], errors="coerce")
    valid_mask = y.notna().to_numpy()

    X = geno_temp[ax_cols].to_numpy()[valid_mask]
    y = y.to_numpy()[valid_mask].astype(float)
    ids = np.array(common_ids_sorted)[valid_mask]

    logger.info(f"Final aligned shape -> Features X: {X.shape}, Target y: {y.shape}")
    return X, y, ids

def run_global_tuning(X, y, model_name, n_iter):
    base_model, search_space = get_search_spaces()[model_name]
    logger.info(
        f"--- Launching Global Tuning Sequence for {model_name} "
        f"({n_iter} iterations) ---"
    )

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    cv = KFold(n_splits=5, shuffle=True, random_state=123)

    opt = BayesSearchCV(
        estimator=base_model,
        search_spaces=search_space,
        n_iter=n_iter,
        cv=cv,
        scoring=pearson_scorer,
        n_jobs=-1,
        verbose=0,
        random_state=123,
    )

    opt.fit(X_scaled, y)

    logger.info(
        f"[GLOBAL_OPTIMAL] {model_name} -> Best CV Pearson R: {opt.best_score_:.6f}"
    )

    return opt.best_params_


def build_model_with_params(name, params):
    if name == "GB":
        return XGBRegressor(
            tree_method="hist",
            device="cpu",
            objective="reg:squarederror",
            eval_metric="rmse",
            random_state=123,
            **params,
        )
    raise ValueError(f"Unknown model profile: {name}")

def regression_metrics(y_true, y_pred):
    mse = mean_squared_error(y_true, y_pred)
    rmse = np.sqrt(mse)
    mae = mean_absolute_error(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)
    pearson_r = pearson_corr_func(y_true, y_pred)
    return mse, rmse, mae, r2, pearson_r

def main():
    X, y, ids = load_and_align_data(args.genofile, args.phenofile)
    model_names = list(get_search_spaces().keys())

    global_optimal_parameters = {}

    for name in model_names:
        best_params = run_global_tuning(X, y, name, args.hpt_iter)
        global_optimal_parameters[name] = best_params
        print(f"[TUNED_PARAMS] Global Optimization Complete | {name} -> {best_params}")

    with open(os.path.join(BASE_OUT_DIR, "global_best_hyperparams.json"), "w") as f:
        json.dump(global_optimal_parameters, f, indent=4)

    logger.info("10x Repeated 5-Fold Cross Validation")

    test_predictions = defaultdict(list)
    all_metrics = []

    for repeat in range(10):
        logger.info(f"\n{'=' * 20} STARTING REPETITION {repeat + 1}/10 {'=' * 20}")

        kf_outer = KFold(n_splits=5, shuffle=True, random_state=100 + repeat)

        for fold, (train_idx, test_idx) in enumerate(kf_outer.split(X)):
            logger.info(f"--- Repetition {repeat + 1} | Outer Fold {fold + 1}/5 ---")

            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]
            ids_test = ids[test_idx]

            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)
            X_test = scaler.transform(X_test)

            for name in model_names:
                params = global_optimal_parameters[name]
                model = build_model_with_params(name, params)
                model.fit(X_train, y_train)

                preds = model.predict(X_test)

                for idx, sample_id in enumerate(ids_test):
                    test_predictions[(sample_id, name)].append(preds[idx])

                mse, rmse, mae, r2, pearson_r = regression_metrics(y_test, preds)

                all_metrics.append(
                    {
                        "Repeat": repeat + 1,
                        "Fold": fold + 1,
                        "Model": name,
                        "MSE": mse,
                        "RMSE": rmse,
                        "MAE": mae,
                        "R2": r2,
                        "PearsonR": pearson_r,
                    }
                )

    logger.info("Processing complete. Summarizing predictions and empirical metrics...")

    final_records = []

    for sample_id in np.unique(ids):
        phenotype_val = y[ids == sample_id][0]

        record = {
            "ID": sample_id,
            "phenotype_gwas": phenotype_val,
        }

        for name in model_names:
            preds = test_predictions[(sample_id, name)]

            if preds:
                record[f"{name}_GEBV_Mean"] = np.mean(preds)
                record[f"{name}_GEBV_SD"] = (
                    np.std(preds, ddof=1) if len(preds) > 1 else 0.0
                )
            else:
                record[f"{name}_GEBV_Mean"] = np.nan
                record[f"{name}_GEBV_SD"] = np.nan

        final_records.append(record)

    gebv_df = pd.DataFrame(final_records).sort_values("ID")
    gebv_out_path = os.path.join(BASE_OUT_DIR, "GB_Final_GEBVs_Summary.csv")
    gebv_df.to_csv(gebv_out_path, index=False)
    logger.info(f"Saved aggregated breeding entries matrix to: {gebv_out_path}")

    metrics_df = pd.DataFrame(all_metrics)
    metrics_out_path = os.path.join(BASE_OUT_DIR, "GB_CrossValidation_Metrics.csv")
    metrics_df.to_csv(metrics_out_path, index=False)
    logger.info(f"Saved metric tracking parameters sheet to: {metrics_out_path}")

if __name__ == "__main__":
    main()