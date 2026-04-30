"""Plot NVT temperature profiles from the Task 3 GROMACS XVG files.

This script reads the ``temperature_edr.xvg`` files for both solvent systems
with pandas.  Each XVG file should be generated from the corresponding
``nvt.edr`` with ``gmx energy``.  Using the ``.edr``-derived XVG files is more
standard than parsing the human-readable ``nvt.log`` files.

The only output is the report-ready figure:

* ``report/static/nvt_temperature_comparison.png``

There are intentionally no command-line arguments, following the project
workflow.  Run it from the repository root with:

    python scripts/task3/plot_nvt_temperature.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
PNG_PATH = ROOT / "report" / "static" / "nvt_temperature_comparison.png"
SYSTEMS = {
    "water": {
        "label": "Pure water",
        "xvg_path": ROOT / "systems" / "water" / "NVT" / "temperature_edr.xvg",
        "color": "#1f77b4",
    },
    "water_ethanol": {
        "label": "Water/ethanol",
        "xvg_path": ROOT / "systems" / "water_ethanol" / "NVT" / "temperature_edr.xvg",
        "color": "#ff7f0e",
    },
}


def read_temperature_xvg(xvg_path: Path) -> pd.DataFrame:
    """Read a two-column GROMACS XVG temperature file into a DataFrame.

    GROMACS XVG files contain metadata lines beginning with ``#`` or ``@``.
    In the files produced by ``gmx energy``, those metadata lines appear before
    the numeric table.  This helper counts that header block and then lets
    pandas parse the whitespace-separated numeric columns directly.

    Returns:
        A DataFrame with ``time_ps`` and ``temperature_K`` columns.
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
        names=["time_ps", "temperature_K"],
        usecols=[0, 1],
        engine="python",
    )
    if data.empty:
        raise ValueError(f"No temperature samples were parsed from {xvg_path}")

    return data


def plot_temperature_comparison(datasets: dict[str, pd.DataFrame], png_path: Path) -> None:
    """Create a combined PNG plot for comparing NVT temperature equilibration.

    Args:
        datasets: Mapping from system key to a DataFrame with ``time_ps`` and
            ``temperature_K`` columns.  The keys must match ``SYSTEMS`` so that
            labels and colors remain centralized.
        png_path: Destination path for the report-ready comparison figure.
    """

    import matplotlib.pyplot as plt

    png_path.parent.mkdir(parents=True, exist_ok=True)

    fig, axis = plt.subplots(figsize=(7.6, 4.4), dpi=180)
    axis.axhline(300.0, color="#d62728", linestyle="--", linewidth=1.0, label="Target 300 K")

    min_time = None
    max_time = None
    for system_name, data in datasets.items():
        config = SYSTEMS[system_name]
        average_temperature = data["temperature_K"].mean()
        label = f"{config['label']} (mean {average_temperature:.1f} K)"

        # Plot the raw temperature trace so short-time thermostat fluctuations
        # remain visible in the report figure.
        axis.plot(
            data["time_ps"],
            data["temperature_K"],
            color=config["color"],
            linewidth=1.2,
            alpha=0.86,
            label=label,
        )

        current_min_time = data["time_ps"].min()
        current_max_time = data["time_ps"].max()
        min_time = current_min_time if min_time is None else min(min_time, current_min_time)
        max_time = current_max_time if max_time is None else max(max_time, current_max_time)

    axis.set_title("NVT Temperature Equilibration")
    axis.set_xlabel("Time (ps)")
    axis.set_ylabel("Temperature (K)")
    axis.set_xlim(min_time, max_time)
    axis.grid(True, linewidth=0.5, alpha=0.35)
    axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(png_path)
    plt.show()


def main() -> None:
    """Read both NVT XVG files and write one comparison figure."""

    datasets = {
        system_name: read_temperature_xvg(config["xvg_path"])
        for system_name, config in SYSTEMS.items()
    }
    plot_temperature_comparison(datasets, PNG_PATH)

    for system_name, data in datasets.items():
        config = SYSTEMS[system_name]
        print(f"[{config['label']}] Parsed {len(data)} temperature samples from {config['xvg_path']}")
        print(
            f"[{config['label']}] Temperature range: "
            f"{data['temperature_K'].min():.2f}-{data['temperature_K'].max():.2f} K"
        )
        print(f"[{config['label']}] Mean temperature: {data['temperature_K'].mean():.2f} K")
    print(f"Wrote {PNG_PATH}")


if __name__ == "__main__":
    main()
