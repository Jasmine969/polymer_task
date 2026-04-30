"""Plot NPT temperature, pressure, and volume from Task 3 GROMACS XVG files.

This script reads the NPT ``energy.xvg`` files for both solvent systems with
pandas.  Each XVG file is expected to be generated from the corresponding
``npt.edr`` with ``gmx energy`` and to contain the columns:

* Temperature
* Pressure
* Volume

The only output is the report-ready figure:

* ``report/static/npt_energy_comparison.png``

There are intentionally no command-line arguments, following the project
workflow.  Run it from the repository root with:

    python scripts/task3/plot_npt_energy.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
PNG_PATH = ROOT / "report" / "static" / "npt_energy_comparison.png"
STABLE_VOLUME_WINDOW_PS = 20.0
SYSTEMS = {
    "water": {
        "label": "Pure water",
        "xvg_path": ROOT / "systems" / "water" / "NPT" / "energy.xvg",
        "color": "#1f77b4",
    },
    "water_ethanol": {
        "label": "Water/ethanol",
        "xvg_path": ROOT / "systems" / "water_ethanol" / "NPT" / "energy.xvg",
        "color": "#ff7f0e",
    },
}


def read_npt_energy_xvg(xvg_path: Path) -> pd.DataFrame:
    """Read the NPT energy XVG file into a DataFrame.

    GROMACS XVG files contain metadata lines beginning with ``#`` or ``@``.
    The numeric table in this workflow has four whitespace-separated columns:
    time in ps, temperature in K, pressure in bar, and volume in nm^3.

    Returns:
        A DataFrame with ``time_ps``, ``temperature_K``, ``pressure_bar``, and
        ``volume_nm3`` columns.
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
        names=["time_ps", "temperature_K", "pressure_bar", "volume_nm3"],
        usecols=[0, 1, 2, 3],
        engine="python",
    )
    if data.empty:
        raise ValueError(f"No NPT energy samples were parsed from {xvg_path}")

    return data


def plot_npt_energy_comparison(datasets: dict[str, pd.DataFrame], png_path: Path) -> None:
    """Create a three-panel NPT diagnostic comparison plot for the report.

    Args:
        datasets: Mapping from system key to a DataFrame with ``time_ps``,
            ``temperature_K``, ``pressure_bar``, and ``volume_nm3`` columns.
            The keys must match ``SYSTEMS`` so labels and colors stay in one
            place.
        png_path: Destination path for the report-ready comparison figure.
    """

    import matplotlib.pyplot as plt

    png_path.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(3, 1, figsize=(7.8, 7.4), dpi=180, sharex=True)

    axes[0].axhline(300.0, color="#d62728", linestyle="--", linewidth=1.0, label="Target 300 K")
    axes[0].set_ylabel("Temperature (K)")
    axes[0].set_title("NPT Equilibration")

    axes[1].axhline(1.0, color="#d62728", linestyle="--", linewidth=1.0, label="Target 1 bar")
    axes[1].set_ylabel("Pressure (bar)")

    axes[2].set_xlabel("Time (ps)")
    axes[2].set_ylabel("Volume (nm$^3$)")

    min_time = None
    max_time = None
    for system_name, data in datasets.items():
        config = SYSTEMS[system_name]
        label = config["label"]
        color = config["color"]

        # Use the same color for all three observables so each solvent system
        # can be tracked across panels at a glance.
        axes[0].plot(data["time_ps"], data["temperature_K"], color=color, linewidth=1.15, alpha=0.88, label=label)
        axes[1].plot(data["time_ps"], data["pressure_bar"], color=color, linewidth=1.0, alpha=0.82, label=label)
        axes[2].plot(data["time_ps"], data["volume_nm3"], color=color, linewidth=1.15, alpha=0.88, label=label)

        current_min_time = data["time_ps"].min()
        current_max_time = data["time_ps"].max()
        min_time = current_min_time if min_time is None else min(min_time, current_min_time)
        max_time = current_max_time if max_time is None else max(max_time, current_max_time)

        stable_start_time = current_max_time - STABLE_VOLUME_WINDOW_PS
        stable_window = data[data["time_ps"] >= stable_start_time]
        stable_volume = stable_window["volume_nm3"].mean()

        axes[2].axhline(stable_volume, color=color, linestyle=":", linewidth=1.0, alpha=0.9)
        axes[2].annotate(
            f"{label}: {stable_volume:.1f} nm$^3$",
            xy=(current_max_time, stable_volume),
            xytext=(-4, 5),
            textcoords="offset points",
            color=color,
            fontsize=8.4,
            ha="right",
            va="bottom",
        )

    for axis in axes:
        axis.grid(True, linewidth=0.5, alpha=0.35)

    axes[0].legend(frameon=False, loc="lower center", ncol=3)
    axes[1].legend(frameon=False, loc="lower center", ncol=3)
    axes[-1].set_xlim(min_time, max_time)
    fig.tight_layout()
    fig.savefig(png_path)
    plt.show()


def main() -> None:
    """Read both NPT XVG files and write one comparison figure."""

    datasets = {
        system_name: read_npt_energy_xvg(config["xvg_path"])
        for system_name, config in SYSTEMS.items()
    }
    plot_npt_energy_comparison(datasets, PNG_PATH)

    for system_name, data in datasets.items():
        config = SYSTEMS[system_name]
        print(f"[{config['label']}] Parsed {len(data)} NPT energy samples from {config['xvg_path']}")
        print(f"[{config['label']}] Mean temperature: {data['temperature_K'].mean():.2f} K")
        print(f"[{config['label']}] Mean pressure: {data['pressure_bar'].mean():.2f} bar")
        print(
            f"[{config['label']}] Volume range: "
            f"{data['volume_nm3'].min():.2f}-{data['volume_nm3'].max():.2f} nm^3"
        )
        stable_start_time = data["time_ps"].max() - STABLE_VOLUME_WINDOW_PS
        stable_volume = data[data["time_ps"] >= stable_start_time]["volume_nm3"].mean()
        print(
            f"[{config['label']}] Final {STABLE_VOLUME_WINDOW_PS:.0f} ps mean volume: "
            f"{stable_volume:.2f} nm^3"
        )
    print(f"Wrote {PNG_PATH}")


if __name__ == "__main__":
    main()
