"""Plot production energy diagnostics from GROMACS XVG files.

This script reads the production ``energy.xvg`` files for both Task 3 solvent
systems with pandas.  Each XVG file is expected to be generated from
``prod.edr`` with ``gmx energy`` and to contain the columns:

* Temperature
* Pressure
* Potential

The only output is the report-ready figure:

* ``report/static/prod_energy_comparison.png``

There are intentionally no command-line arguments, following the project
workflow.  Run it from the repository root with:

    python scripts/task3/plot_prod_energy.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
PNG_PATH = ROOT / "report" / "static" / "prod_energy_comparison.png"
SYSTEMS = {
    "water": {
        "label": "Pure water",
        "xvg_path": ROOT / "systems" / "water" / "PROD" / "energy.xvg",
        "color": "#1f77b4",
    },
    "water_ethanol": {
        "label": "Water/ethanol",
        "xvg_path": ROOT / "systems" / "water_ethanol" / "PROD" / "energy.xvg",
        "color": "#ff7f0e",
    },
}


def read_prod_energy_xvg(xvg_path: Path) -> pd.DataFrame:
    """Read the production energy XVG file into a DataFrame.

    GROMACS XVG files contain metadata lines beginning with ``#`` or ``@``.
    The numeric table in this workflow has four whitespace-separated columns:
    time in ps, potential energy in kJ/mol, temperature in K, and pressure in
    bar.  The order follows the legends written by ``gmx energy``.

    Returns:
        A DataFrame with ``time_ns``, ``temperature_K``, ``pressure_bar``, and
        ``potential_kj_mol`` columns.
    """

    skiprows = 0
    for line in xvg_path.read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.lstrip()
        if stripped.startswith(("#", "@")) or not stripped:
            skiprows += 1
            continue
        break

    data = pd.read_csv(
        xvg_path,
        sep=r"\s+",
        skiprows=skiprows,
        names=["time_ps", "potential_kj_mol", "temperature_K", "pressure_bar"],
        usecols=[0, 1, 2, 3],
        engine="python",
    )
    if data.empty:
        raise ValueError(f"No production energy samples were parsed from {xvg_path}")

    data["time_ns"] = data["time_ps"] / 1000.0
    return data


def plot_prod_energy_comparison(datasets: dict[str, pd.DataFrame], png_path: Path) -> None:
    """Create a three-panel production diagnostic comparison plot.

    Args:
        datasets: Mapping from system key to a DataFrame with ``time_ns``,
            ``temperature_K``, ``pressure_bar``, and ``potential_kj_mol``
            columns.  The keys must match ``SYSTEMS`` for consistent labels and
            colors.
        png_path: Destination path for the report-ready comparison figure.
    """
    import matplotlib.pyplot as plt

    png_path.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(3, 1, figsize=(7.8, 7.4), dpi=180, sharex=True)

    axes[0].axhline(300.0, color="#d62728", linestyle="--", linewidth=1.0, label="Target 300 K")
    axes[0].set_ylabel("Temperature (K)")
    axes[0].set_title("Production Diagnostics")

    axes[1].set_ylabel("Potential (kJ/mol)")

    axes[2].set_xlabel("Time (ns)")
    axes[2].set_ylabel("Pressure (bar)")

    min_time = None
    max_time = None
    for system_name, data in datasets.items():
        config = SYSTEMS[system_name]
        label = config["label"]
        color = config["color"]

        # Keep one color per solvent system across all observables.
        axes[0].plot(data["time_ns"], data["temperature_K"], color=color, linewidth=1.0, alpha=0.88, label=label)
        axes[1].plot(data["time_ns"], data["potential_kj_mol"], color=color, linewidth=1.0, alpha=0.88, label=label)
        axes[2].plot(data["time_ns"], data["pressure_bar"], color=color, linewidth=0.9, alpha=0.82, label=label)

        current_min_time = data["time_ns"].min()
        current_max_time = data["time_ns"].max()
        min_time = current_min_time if min_time is None else min(min_time, current_min_time)
        max_time = current_max_time if max_time is None else max(max_time, current_max_time)

    for axis in axes:
        axis.grid(True, linewidth=0.5, alpha=0.35)

    axes[0].legend(frameon=False, loc="best", ncol=3)
    axes[-1].set_xlim(min_time, max_time)
    fig.tight_layout()
    fig.savefig(png_path)
    plt.show()


def main() -> None:
    """Read available production XVG files and write one comparison figure."""

    datasets: dict[str, pd.DataFrame] = {}
    for system_name, config in SYSTEMS.items():
        xvg_path = config["xvg_path"]
        if not xvg_path.exists():
            print(f"Skipped {system_name}: {xvg_path} does not exist")
            continue

        datasets[system_name] = read_prod_energy_xvg(xvg_path)

    if not datasets:
        raise FileNotFoundError("No production energy XVG files were found")

    plot_prod_energy_comparison(datasets, PNG_PATH)

    for system_name, data in datasets.items():
        config = SYSTEMS[system_name]
        xvg_path = config["xvg_path"]
        print(f"[{system_name}] Parsed {len(data)} production energy samples from {xvg_path}")
        print(f"[{system_name}] Simulation length: {data['time_ns'].max():.2f} ns")
        print(f"[{system_name}] Mean temperature: {data['temperature_K'].mean():.2f} K")
        print(f"[{system_name}] Mean pressure: {data['pressure_bar'].mean():.2f} bar")
        print(f"[{system_name}] Mean potential: {data['potential_kj_mol'].mean():.2f} kJ/mol")
    print(f"Wrote {PNG_PATH}")


if __name__ == "__main__":
    main()
