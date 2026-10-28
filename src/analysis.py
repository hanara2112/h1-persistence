"""Evaluate circle/disk persistence, noise sensitivity, and Ripser runtime."""

import json
import platform
from importlib.metadata import version
from pathlib import Path
from time import perf_counter

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from persim import bottleneck, plot_diagrams
from ripser import ripser

NOISE_LEVELS = [0.0, 0.1, 0.3, 0.5, 1.0]
N_POINTS = 100
SEEDS = list(range(42, 52))
EXAMPLE_SEED = 42
SAMPLE_SIZES = [30, 100, 300]
BENCHMARK_NOISE = 0.1
RUNTIME_REPEATS = 3
SHAPES = ["circle", "disk"]
COLORS = {"circle": "tab:blue", "disk": "tab:orange"}
RESULTS_DIR = Path(__file__).resolve().parents[1] / "results"


def make_samples(n_points, seed):
    """Return clean circle/disk samples and shared standard Gaussian noise."""
    rng = np.random.default_rng(seed)
    angles = rng.uniform(0, 2 * np.pi, n_points)
    circle = np.column_stack([np.cos(angles), np.sin(angles)])
    noise = rng.standard_normal(circle.shape)
    radius = np.sqrt(rng.uniform(0, 1, n_points))
    disk = radius[:, None] * circle
    return {"circle": circle, "disk": disk}, noise


def h1_summary(h1):
    """Summarize finite H1 intervals; use zero when an interval is absent."""
    finite = h1[np.isfinite(h1).all(axis=1)]
    lifetimes = np.sort(finite[:, 1] - finite[:, 0])[::-1]
    longest = lifetimes[0] if len(lifetimes) else 0.0
    second = lifetimes[1] if len(lifetimes) > 1 else 0.0
    return {
        "num_H1_intervals": len(finite),
        "max_H1_persistence": longest,
        "second_H1_persistence": second,
        "dominance_gap": longest - second,
    }


def run_noise_experiment():
    records = []
    examples = {}
    for seed in SEEDS:
        samples, noise = make_samples(N_POINTS, seed)
        for shape in SHAPES:
            clean_diagrams = ripser(samples[shape], maxdim=1)["dgms"]
            for level in NOISE_LEVELS:
                points = samples[shape] + level * noise
                diagrams = (clean_diagrams if level == 0 else
                            ripser(points, maxdim=1)["dgms"])
                metrics = h1_summary(diagrams[1])
                centered = points - points.mean(axis=0)
                rms_radius = np.sqrt(np.mean(np.sum(centered ** 2, axis=1)))
                records.append({
                    "shape": shape, "seed": seed, "noise": level, **metrics,
                    "rms_radius": rms_radius,
                    "relative_persistence": metrics["max_H1_persistence"] / rms_radius,
                    "bottleneck_to_clean": bottleneck(clean_diagrams[1], diagrams[1]),
                })
                if seed == EXAMPLE_SEED:
                    examples[shape, level] = (points, diagrams)
    return pd.DataFrame(records), examples


def summarize_noise(trials):
    return trials.groupby(["shape", "noise"]).agg(
        max_mean=("max_H1_persistence", "mean"),
        max_sd=("max_H1_persistence", "std"),
        second_mean=("second_H1_persistence", "mean"),
        gap_mean=("dominance_gap", "mean"),
        gap_sd=("dominance_gap", "std"),
        intervals_mean=("num_H1_intervals", "mean"),
        bottleneck_mean=("bottleneck_to_clean", "mean"),
        bottleneck_sd=("bottleneck_to_clean", "std"),
    ).reset_index()


def compare_trials(trials):
    circle = trials[trials["shape"] == "circle"]
    raw = circle.pivot(index="seed", columns="noise", values="max_H1_persistence")
    relative = circle.pivot(index="seed", columns="noise", values="relative_persistence")
    rebound = pd.DataFrame({
        "raw_at_0.5": raw[0.5],
        "raw_at_1.0": raw[1.0],
        "raw_change": raw[1.0] - raw[0.5],
        "relative_change": relative[1.0] - relative[0.5],
    }).reset_index()
    paired = trials.pivot(index=["noise", "seed"], columns="shape",
                          values="max_H1_persistence")
    circle_larger = (paired["circle"] > paired["disk"]).groupby("noise").sum()
    comparison = pd.DataFrame({
        "circle_larger_count": circle_larger, "total_pairs": len(SEEDS),
    }).reset_index()
    return rebound, comparison


def run_benchmark():
    clean, _ = make_samples(N_POINTS, EXAMPLE_SEED)
    ripser(clean["circle"], maxdim=1)  # Warm-up outside the timer.
    records = []
    for n_points in SAMPLE_SIZES:
        for seed in SEEDS:
            samples, noise = make_samples(n_points, seed)
            for shape in SHAPES:
                points = samples[shape] + BENCHMARK_NOISE * noise
                timings = []
                for _ in range(RUNTIME_REPEATS):
                    start = perf_counter()
                    diagrams = ripser(points, maxdim=1)["dgms"]
                    timings.append(1000 * (perf_counter() - start))
                records.append({
                    "shape": shape, "n_points": n_points, "seed": seed,
                    "median_ms": np.median(timings),
                    "max_H1_persistence": h1_summary(diagrams[1])["max_H1_persistence"],
                })
    trials = pd.DataFrame(records)
    summary = trials.groupby(["shape", "n_points"]).agg(
        runtime_mean_ms=("median_ms", "mean"),
        runtime_sd_ms=("median_ms", "std"),
        persistence_mean=("max_H1_persistence", "mean"),
        persistence_sd=("max_H1_persistence", "std"),
    ).reset_index()
    return trials, summary


def save_figure(figure, path):
    figure.savefig(path, dpi=160)
    plt.close(figure)


def plot_clean_shapes(examples, output_dir):
    figure, axes = plt.subplots(2, 2, figsize=(9, 7), layout="constrained")
    for column, shape in enumerate(SHAPES):
        points, diagrams = examples[shape, 0.0]
        ax = axes[0, column]
        ax.scatter(*points.T, s=16, color=COLORS[shape])
        ax.set(xlim=(-1.2, 1.2), ylim=(-1.2, 1.2),
               title=f"Clean {shape}", xlabel="x", ylabel="y")
        ax.set_aspect("equal")
        plot_diagrams(diagrams, ax=axes[1, column], show=False,
                      xy_range=[-0.1, 2.1, -0.1, 2.1],
                      title=f"{shape.capitalize()}: H0 and H1")
    save_figure(figure, output_dir / "01_clean_shapes.png")


def plot_noise_examples(examples, output_dir):
    circle_examples = [examples["circle", level] for level in NOISE_LEVELS]
    coordinate_limit = 1.1 * max(np.abs(points).max() for points, _ in circle_examples)
    diagram_limit = 1.1 * max(
        diagram[np.isfinite(diagram)].max()
        for _, diagrams in circle_examples for diagram in diagrams if diagram.size
    )
    figure, axes = plt.subplots(2, len(NOISE_LEVELS), figsize=(16, 6),
                               layout="constrained")
    for column, level in enumerate(NOISE_LEVELS):
        points, diagrams = examples["circle", level]
        ax = axes[0, column]
        ax.scatter(*points.T, s=12)
        ax.set(xlim=(-coordinate_limit, coordinate_limit),
               ylim=(-coordinate_limit, coordinate_limit),
               title=f"Noise = {level}", xlabel="x", ylabel="y")
        ax.set_aspect("equal")
        plot_diagrams(diagrams, ax=axes[1, column], show=False,
                      xy_range=[-0.1, diagram_limit, -0.1, diagram_limit])
    figure.suptitle(f"One circle trial (seed {EXAMPLE_SEED}); shared axes within each row")
    save_figure(figure, output_dir / "02_circle_noise_examples.png")


def plot_errorbars(ax, summary, x_column, mean_column, sd_column):
    for shape in SHAPES:
        values = summary[summary["shape"] == shape]
        ax.errorbar(values[x_column], values[mean_column], yerr=values[sd_column],
                    marker="o", capsize=3, color=COLORS[shape], label=shape)
    ax.grid(alpha=0.25)
    ax.legend()


def plot_noise_summary(summary, output_dir):
    figure, axes = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
    for ax, metric, title in zip(axes, ["max", "gap"],
                                ["Longest H1 interval", "Dominance gap: longest − second"]):
        plot_errorbars(ax, summary, "noise", f"{metric}_mean", f"{metric}_sd")
        ax.set(xlabel="Noise standard deviation", ylabel="Persistence (distance units)",
               title=title)
    save_figure(figure, output_dir / "03_noise_summary.png")


def plot_bottleneck(summary, output_dir):
    figure, ax = plt.subplots(figsize=(7, 4), layout="constrained")
    plot_errorbars(ax, summary, "noise", "bottleneck_mean", "bottleneck_sd")
    ax.set(xlabel="Noise standard deviation", ylabel="H1 bottleneck distance",
           title="Change from each trial's clean reference (mean ± SD)")
    save_figure(figure, output_dir / "04_bottleneck.png")


def plot_benchmark(summary, output_dir):
    figure, axes = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
    plot_errorbars(axes[0], summary, "n_points", "runtime_mean_ms", "runtime_sd_ms")
    plot_errorbars(axes[1], summary, "n_points", "persistence_mean", "persistence_sd")
    axes[0].set(ylabel="Ripser time (ms)", title="Computation time (mean ± SD)")
    axes[1].set(ylabel="Longest H1 persistence", title="Sampling density (mean ± SD)")
    for ax in axes:
        ax.set(xlabel="Number of points", xticks=SAMPLE_SIZES)
    save_figure(figure, output_dir / "05_density_runtime.png")


def save_environment(output_dir):
    packages = ["numpy", "pandas", "matplotlib", "scipy", "scikit-learn", "ripser", "persim"]
    environment = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "packages": {name: version(name) for name in packages},
        "settings": {
            "noise_levels": NOISE_LEVELS, "n_points": N_POINTS, "seeds": SEEDS,
            "example_seed": EXAMPLE_SEED, "sample_sizes": SAMPLE_SIZES,
            "benchmark_noise": BENCHMARK_NOISE, "runtime_repeats": RUNTIME_REPEATS,
            "maxdim": 1, "metric": "euclidean", "coeff": 2,
        },
    }
    (output_dir / "environment.json").write_text(json.dumps(environment, indent=2) + "\n")


def main():
    RESULTS_DIR.mkdir(exist_ok=True)
    plt.rcParams.update({"figure.dpi": 110, "font.size": 10})
    noise_trials, examples = run_noise_experiment()
    noise_summary = summarize_noise(noise_trials)
    rebound, comparison = compare_trials(noise_trials)
    runtime_trials, runtime_summary = run_benchmark()

    tables = {
        "noise_trials": noise_trials, "noise_summary": noise_summary,
        "high_noise_changes": rebound, "circle_disk_comparison": comparison,
        "runtime_trials": runtime_trials, "runtime_summary": runtime_summary,
    }
    for name, table in tables.items():
        table.to_csv(RESULTS_DIR / f"{name}.csv", index=False)
    plot_clean_shapes(examples, RESULTS_DIR)
    plot_noise_examples(examples, RESULTS_DIR)
    plot_noise_summary(noise_summary, RESULTS_DIR)
    plot_bottleneck(noise_summary, RESULTS_DIR)
    plot_benchmark(runtime_summary, RESULTS_DIR)
    save_environment(RESULTS_DIR)

    print(noise_summary.round(3).to_string(index=False))
    print(runtime_summary.round(3).to_string(index=False))
    print(f"High-noise increases: raw {(rebound['raw_change'] > 0).sum()}/{len(SEEDS)}, "
          f"relative {(rebound['relative_change'] > 0).sum()}/{len(SEEDS)}")
    print(f"Saved {len(tables)} CSV tables and 5 figures to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
