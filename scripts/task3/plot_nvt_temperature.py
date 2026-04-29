"""Plot the temperature profile from the pure-water NVT GROMACS XVG file.

This script reads ``systems/water/NVT/temperature_edr.xvg`` with pandas.  The
XVG file should be generated from ``systems/water/NVT/nvt.edr`` with
``gmx energy``.  Using the ``.edr``-derived XVG file is more standard than
parsing the human-readable ``nvt.log`` file.

The only output is the report-ready figure:

* ``report/static/water_nvt_temperature.png``

There are intentionally no command-line arguments, following the project
workflow.  Run it from the repository root with:

    python scripts/task3/plot_nvt_temperature.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
XVG_PATH = ROOT / "systems" / "water" / "NVT" / "temperature_edr.xvg"
PNG_PATH = ROOT / "report" / "static" / "water_nvt_temperature.png"


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


def plot_temperature(data: pd.DataFrame, png_path: Path) -> None:
    """Create a PNG plot for checking NVT temperature equilibration."""
    import matplotlib.pyplot as plt

    average_temperature = data["temperature_K"].mean()

    png_path.parent.mkdir(parents=True, exist_ok=True)

    fig, axis = plt.subplots(figsize=(7.0, 4.2), dpi=180)
    axis.plot(
        data["time_ps"],
        data["temperature_K"],
        color="#1f77b4",
        linewidth=1.4,
        label="NVT temperature",
    )
    axis.axhline(300.0, color="#d62728", linestyle="--", linewidth=1.0, label="Target 300 K")
    axis.axhline(
        average_temperature,
        color="#2ca02c",
        linestyle=":",
        linewidth=1.2,
        label=f"Mean {average_temperature:.1f} K",
    )

    axis.set_title("Pure Water NVT Temperature")
    axis.set_xlabel("Time (ps)")
    axis.set_ylabel("Temperature (K)")
    axis.set_xlim(data["time_ps"].min(), data["time_ps"].max())
    axis.grid(True, linewidth=0.5, alpha=0.35)
    axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(png_path)
    plt.show()


def main() -> None:
    """Read the pure-water NVT XVG file and write the report figure."""

    data = read_temperature_xvg(XVG_PATH)
    plot_temperature(data, PNG_PATH)

    print(f"Parsed {len(data)} temperature samples from {XVG_PATH}")
    print(
        "Temperature range: "
        f"{data['temperature_K'].min():.2f}-{data['temperature_K'].max():.2f} K"
    )
    print(f"Mean temperature: {data['temperature_K'].mean():.2f} K")
    print(f"Wrote {PNG_PATH}")


if __name__ == "__main__":
    main()
