#!/usr/bin/env python3
"""
Updated script for genomic selection status prediction Cape Shore GS 2026.
Optimized to run global hyperparameter tuning ONCE, then apply optimal parameters 
across the 10x repeated 5-fold cross-validation loops.
"""

import os
import json
import logging
import argparse
import numpy as np
import pandas as pd
from collections import defaultdict
from scipy.stats import pearsonr

from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import make_scorer, roc_auc_score, log_loss, brier_score_loss
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier

from skopt import BayesSearchCV
from skopt.space import Real, Integer, Categorical

# --- CLI Configuration ---
parser = argparse.ArgumentParser(description="Combined HPT and Cross-Validation Pipeline")
parser.add_argument("--genofile", "-g", type=str, required=True,
                    help="Path to the Genotype CSV file")
parser.add_argument("--phenofile", "-p", type=str, required=True,
                    help="Path to the Phenotype CSV file")
parser.add_argument("--outdir", "-o", type=str, required=True,
                    help="Subfolder name inside /work/hs325/csgs2026/ to save results")
parser.add_argument("--hpt_iter", type=int, default=50,
                    help="Number of Global Bayesian Optimization iterations, defaults 50 for speed")
parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
args = parser.parse_args()

# Setup logging
logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                    format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Base output directory matching requirements
BASE_OUT_DIR = os.path.join("/work/hs325/csgs2026", args.outdir)
os.makedirs(BASE_OUT_DIR, exist_ok=True)


# --- Search Spaces Definition ---
def get_search_spaces():
    return {
        "LR": (
            LogisticRegression(max_iter=1000, solver="saga"),
            {
                "C": Real(1e-5, 10, prior="log-uniform"),
                "penalty": Categorical(["l1", "l2"]),
            },
        ),
        "RF": (
            RandomForestClassifier(n_jobs=-1),
            {
                "n_estimators": Integer(100, 2000),
                "max_depth": Integer(3, 50),
                "max_features": Categorical(["sqrt", "log2"]),
                "min_samples_split": Integer(2, 20),
            },
        ),
        "GB": (
            XGBClassifier(tree_method="hist", device="cuda", eval_metric="logloss"),
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

pearson_scorer = make_scorer(pearson_corr_func, response_method="predict_proba")


# Align genetic data to phenotypic sample ID using explicitly saved columns
def load_and_align_data(geno_path, pheno_path):
    logger.info("Loading genotype and phenotype files...")
    
    # Load files directly (retaining their native numeric range indices)
    geno_df = pd.read_csv(geno_path)
    pheno_df = pd.read_csv(pheno_path)
    
    # Set ID columns as temporary lookup indices to isolate clean overlaps
    geno_temp = geno_df.set_index("ID")
    pheno_temp = pheno_df.set_index("SampleID")
    
    # Find intersecting Sample IDs
    common_ids = geno_temp.index.intersection(pheno_temp.index)
    logger.info(f"Intersecting samples discovered: {len(common_ids)}")
    
    if len(common_ids) == 0:
        raise ValueError("No matching entries found between Genotype 'ID' and Phenotype 'SampleID' columns.")
        
    common_ids_sorted = sorted(list(common_ids))
    geno_temp = geno_temp.loc[common_ids_sorted]
    pheno_temp = pheno_temp.loc[common_ids_sorted]
    
    # Isolate status targets and marker features
    y = pheno_temp["status_01"].to_numpy()
    ax_cols = [col for col in geno_temp.columns if col.startswith("AX")]
    X = geno_temp[ax_cols].to_numpy()
    ids = np.array(common_ids_sorted)
    
    logger.info(f"Final aligned shape -> Features X: {X.shape}, Target y: {y.shape}")
    return X, y, ids


# Global Hyperparameter Tuning Execution Block
def run_global_tuning(X, y, model_name, n_iter):
    base_model, search_space = get_search_spaces()[model_name]
    logger.info(f"--- Launching Global Tuning Sequence for {model_name} ({n_iter} iterations) ---")
    
    # Run dynamic scaling strictly on the tuning data representation
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    opt = BayesSearchCV(
        estimator=base_model,
        search_spaces=search_space,
        n_iter=n_iter,
        cv=5,
        scoring=pearson_scorer,
        n_jobs=-1,
        verbose=0,
        random_state=123
    )
    opt.fit(X_scaled, y)
    logger.info(f"[GLOBAL_OPTIMAL] {model_name} -> Best Score: {opt.best_score_:.4f}")
    return opt.best_params_


# Helper to re-instantiate optimized model frameworks inside folds
def build_model_with_params(name, params):
    if name == "LR":
        return LogisticRegression(max_iter=1000, solver="saga", **params)
    elif name == "RF":
        return RandomForestClassifier(n_jobs=-1, **params)
    elif name == "GB":
        return XGBClassifier(tree_method="hist", device="cuda", eval_metric="logloss", **params)
    raise ValueError(f"Unknown framework model profile: {name}")


# --- Main Pipeline ---
def main():
    X, y, ids = load_and_align_data(args.genofile, args.phenofile)
    model_names = list(get_search_spaces().keys())
    
    # Run Global Hyperparameter tuning ONCE per architecture framework
    global_optimal_parameters = {}
    for name in model_names:
        best_params = run_global_tuning(X, y, name, args.hpt_iter)
        global_optimal_parameters[name] = best_params
        print(f"[TUNED_PARAMS] Global Optimization Complete | {name} -> {best_params}")
        
    # Save the global parameters out for user audit records
    with open(os.path.join(BASE_OUT_DIR, "global_best_hyperparams.json"), "w") as f:
        json.dump(global_optimal_parameters, f, indent=4)
        
    # Production CV 
    logger.info("10x Repeated 5-Fold Cross Validation")
    
    test_predictions = defaultdict(list)
    all_metrics = []
    
    # 10 Repetitions of Outer 5-Fold Cross Validation
    for repeat in range(10):
        logger.info(f"\n{'='*20} STARTING REPETITION {repeat + 1}/10 {'='*20}")
        skf_outer = StratifiedKFold(n_splits=5, shuffle=True, random_state=100 + repeat)
        
        for fold, (train_idx, test_idx) in enumerate(skf_outer.split(X, y)):
            logger.info(f"--- Repetition {repeat + 1} | Outer Fold {fold + 1}/5 ---")
            
            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]
            ids_test = ids[test_idx]
            
            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)
            X_test = scaler.transform(X_test)
            
            for name in model_names:
                # Build model directly using static parameter configs discovered in Phase 1
                params = global_optimal_parameters[name]
                model = build_model_with_params(name, params)
                model.fit(X_train, y_train)
                
                # Generate out-of-fold predictions on isolated test subset
                probs = model.predict_proba(X_test)[:, 1]
                
                for idx, sample_id in enumerate(ids_test):
                    test_predictions[(sample_id, name)].append(probs[idx])
                
                # Performance reporting per model framework block
                auc = roc_auc_score(y_test, probs)
                logloss_val = log_loss(y_test, probs, labels=[0, 1])
                brier_val = brier_score_loss(y_test, probs)
                r_val = pearson_corr_func(y_test, probs)
                
                all_metrics.append({
                    "Repeat": repeat + 1,
                    "Fold": fold + 1,
                    "Model": name,
                    "AUC": auc,
                    "LogLoss": logloss_val,
                    "Brier": brier_val,
                    "PearsonR": r_val
                })

    logger.info("Processing complete. Summarizing predictions and empirical metrics...")
    
    # Generate final GEBVs 
    final_records = []
    for sample_id in np.unique(ids):
        status_val = y[ids == sample_id][0]
        record = {"ID": sample_id, "Status": status_val}
        
        for name in model_names:
            preds = test_predictions[(sample_id, name)]
            if preds:
                record[f"{name}_GEBV_Mean"] = np.mean(preds)
                record[f"{name}_GEBV_SD"] = np.std(preds, ddof=1) if len(preds) > 1 else 0.0
            else:
                record[f"{name}_GEBV_Mean"] = np.nan
                record[f"{name}_GEBV_SD"] = np.nan
                
        final_records.append(record)

    # Export genomic estimated breeding dataframe
    gebv_df = pd.DataFrame(final_records).sort_values("ID")
    gebv_out_path = os.path.join(BASE_OUT_DIR, "Final_GEBVs_Summary.csv")
    gebv_df.to_csv(gebv_out_path, index=False)
    logger.info(f"Saved aggregated breeding entries matrix to: {gebv_out_path}")

    # Export metrics dataframe
    metrics_df = pd.DataFrame(all_metrics)
    metrics_out_path = os.path.join(BASE_OUT_DIR, "CrossValidation_Metrics.csv")
    metrics_df.to_csv(metrics_out_path, index=False)
    logger.info(f"Saved metric tracking parameters sheet to: {metrics_out_path}")

if __name__ == "__main__":
    main()