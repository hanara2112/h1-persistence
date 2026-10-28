"""Calibrate a longest-H1 detector on disks, then test independent clouds.

Run from the repository root: python src/hole_detection.py
Outputs are separate from the original study, in results/hole_detection/.
"""

import json
import platform
from importlib.metadata import version
from math import ceil
from pathlib import Path
from time import perf_counter

from analysis import h1_summary, make_samples, save_figure
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ripser import ripser
from scipy.stats import binomtest
from sklearn.metrics import roc_auc_score

SAMPLE_SIZES = [30, 100, 300]
NOISE_LEVELS = [0.0, 0.1, 0.3, 0.5, 1.0]
N_CALIBRATION = 199
N_TEST = 200  # Per shape and condition.
ALPHA = 0.05
BASE_SEED = 424242
OUTPUT_DIR = Path(__file__).resolve().parents[1] / "results" / "hole_detection"


def simulate_scores(shape, n_points, noise_level, seeds):
    records = []
    for seed in seeds:
        samples, noise = make_samples(n_points, seed)
        points = samples[shape] + noise_level * noise
        h1 = ripser(points, maxdim=1)["dgms"][1]
        records.append({
            "shape": shape, "n_points": n_points, "noise": noise_level,
            "seed": seed, "score": h1_summary(h1)["max_H1_persistence"],
        })
    return pd.DataFrame(records)


def null_threshold(scores, alpha=ALPHA):
    """Upper order statistic; a new score must strictly exceed it to detect."""
    rank = ceil((len(scores) + 1) * (1 - alpha))
    if rank > len(scores):
        raise ValueError("More calibration simulations are needed for this alpha.")
    return float(np.sort(scores)[rank - 1])


def rate_interval(successes, total):
    interval = binomtest(successes, total).proportion_ci(method="exact")
    return successes / total, interval.low, interval.high


def evaluate_condition(calibration, disks, circles):
    threshold = null_threshold(calibration["score"].to_numpy())
    evaluation = pd.concat([disks, circles], ignore_index=True)
    evaluation["threshold"] = threshold
    evaluation["detected"] = evaluation["score"] > threshold
    false_positives = int((disks["score"] > threshold).sum())
    true_positives = int((circles["score"] > threshold).sum())
    fpr, fpr_low, fpr_high = rate_interval(false_positives, len(disks))
    tpr, tpr_low, tpr_high = rate_interval(true_positives, len(circles))
    summary = {
        "n_points": int(disks["n_points"].iloc[0]),
        "noise": float(disks["noise"].iloc[0]),
        "n_calibration": len(calibration), "n_test_per_shape": len(disks),
        "alpha": ALPHA, "threshold": threshold,
        "false_positives": false_positives, "true_positives": true_positives,
        "false_positive_rate": fpr, "fpr_ci_low": fpr_low, "fpr_ci_high": fpr_high,
        "detection_rate": tpr, "detection_ci_low": tpr_low, "detection_ci_high": tpr_high,
        "roc_auc": roc_auc_score(evaluation["shape"] == "circle", evaluation["score"]),
        "disk_mean": disks["score"].mean(), "disk_sd": disks["score"].std(),
        "circle_mean": circles["score"].mean(), "circle_sd": circles["score"].std(),
    }
    return evaluation, summary


def plot_detection(summary):
    figure, axes = plt.subplots(1, 3, figsize=(13, 4), layout="constrained")
    for n_points in SAMPLE_SIZES:
        values = summary[summary["n_points"] == n_points]
        for ax, metric, prefix in zip(
            axes[:2], ["false_positive_rate", "detection_rate"], ["fpr", "detection"]
        ):
            error = [values[metric] - values[f"{prefix}_ci_low"],
                     values[f"{prefix}_ci_high"] - values[metric]]
            ax.errorbar(values["noise"], values[metric], yerr=error,
                        marker="o", capsize=3, label=f"n = {n_points}")
        axes[2].plot(values["noise"], values["roc_auc"], marker="o", label=f"n = {n_points}")
    axes[0].axhline(ALPHA, color="black", linestyle=":", label="Nominal 5%")
    axes[2].axhline(0.5, color="black", linestyle=":", label="Chance ranking")
    for ax, title in zip(axes, ["False-positive rate", "Detection rate", "ROC-AUC"]):
        ax.set(xlabel="Gaussian noise SD", ylabel=title, ylim=(-0.03, 1.03), xticks=NOISE_LEVELS)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    axes[0].set(ylim=(0, 0.16), yticks=[0, 0.05, 0.10, 0.15])
    figure.suptitle("Independent test clouds; rate bars are pointwise 95% exact intervals")
    save_figure(figure, OUTPUT_DIR / "detection_metrics.png")


def plot_scores(summary):
    figure, axes = plt.subplots(1, len(SAMPLE_SIZES), figsize=(13, 4), layout="constrained")
    for ax, n_points in zip(axes, SAMPLE_SIZES):
        values = summary[summary["n_points"] == n_points]
        for shape in ["circle", "disk"]:
            ax.errorbar(values["noise"], values[f"{shape}_mean"],
                        yerr=values[f"{shape}_sd"], marker="o", capsize=3, label=shape)
        ax.plot(values["noise"], values["threshold"], "k--", label="Null threshold")
        ax.set(xlabel="Gaussian noise SD", ylabel="Longest H1 persistence",
               title=f"n = {n_points}", ylim=(0, 1.8), xticks=NOISE_LEVELS)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    figure.suptitle("Test-score means ± SD and independently calibrated thresholds")
    save_figure(figure, OUTPUT_DIR / "scores_and_thresholds.png")


def main():
    start = perf_counter()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    calibration_tables, evaluation_tables, summaries = [], [], []
    next_seed = BASE_SEED
    for n_points in SAMPLE_SIZES:
        for noise_level in NOISE_LEVELS:
            calibration = simulate_scores("disk", n_points, noise_level,
                                          range(next_seed, next_seed + N_CALIBRATION))
            next_seed += N_CALIBRATION
            disks = simulate_scores("disk", n_points, noise_level,
                                    range(next_seed, next_seed + N_TEST))
            next_seed += N_TEST
            circles = simulate_scores("circle", n_points, noise_level,
                                      range(next_seed, next_seed + N_TEST))
            next_seed += N_TEST
            evaluation, summary = evaluate_condition(calibration, disks, circles)
            calibration_tables.append(calibration)
            evaluation_tables.append(evaluation)
            summaries.append(summary)
            print(f"n={n_points:3d}, noise={noise_level:.1f}: "
                  f"FPR={summary['false_positive_rate']:.3f}, "
                  f"detection={summary['detection_rate']:.3f}, "
                  f"AUC={summary['roc_auc']:.3f}", flush=True)

    calibration = pd.concat(calibration_tables, ignore_index=True)
    evaluation = pd.concat(evaluation_tables, ignore_index=True)
    summary = pd.DataFrame(summaries)
    calibration.to_csv(OUTPUT_DIR / "null_calibration.csv", index=False)
    evaluation.to_csv(OUTPUT_DIR / "test_scores.csv", index=False)
    summary.to_csv(OUTPUT_DIR / "summary.csv", index=False)
    plot_detection(summary)
    plot_scores(summary)
    environment = {
        "python": platform.python_version(), "platform": platform.platform(),
        "packages": {name: version(name) for name in
                     ["numpy", "pandas", "matplotlib", "scipy", "scikit-learn", "ripser", "persim"]},
        "settings": {
            "sample_sizes": SAMPLE_SIZES, "noise_levels": NOISE_LEVELS,
            "n_calibration": N_CALIBRATION, "n_test_per_shape": N_TEST, "alpha": ALPHA,
            "first_seed": BASE_SEED, "last_seed": next_seed - 1,
            "seed_assignment": "Unique consecutive seeds: size, noise, calibration disk, test disk, test circle",
            "score": "Longest finite H1 persistence, zero if empty",
            "rule": "score > sorted_null[ceil((m+1)*(1-alpha))-1]",
            "maxdim": 1, "metric": "euclidean", "coeff": 2, "distance_cutoff": None,
        },
        "elapsed_seconds": perf_counter() - start,
    }
    (OUTPUT_DIR / "environment.json").write_text(json.dumps(environment, indent=2) + "\n")
    print(f"Saved {len(calibration) + len(evaluation)} simulations to {OUTPUT_DIR}")
    print(f"Elapsed: {environment['elapsed_seconds']:.1f} seconds")


if __name__ == "__main__":
    main()
