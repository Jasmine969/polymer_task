"""Plot pure-water NPT temperature, pressure, and volume from GROMACS XVG.

This script reads ``systems/water/NPT/energy.xvg`` with pandas.  The XVG file
is expected to be generated from ``systems/water/NPT/npt.edr`` with
``gmx energy`` and to contain the columns:

* Temperature
* Pressure
* Volume

The only output is the report-ready figure:

* ``report/static/water_npt_energy.png``

There are intentionally no command-line arguments, following the project
workflow.  Run it from the repository root with:

    python scripts/task3/plot_npt_energy.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
XVG_PATH = ROOT / "systems" / "water" / "NPT" / "energy.xvg"
PNG_PATH = ROOT / "report" / "static" / "water_npt_energy.png"


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


def plot_npt_energy(data: pd.DataFrame, png_path: Path) -> None:
    """Create a three-panel NPT diagnostic plot for the report."""

    import matplotlib.pyplot as plt

    png_path.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(3, 1, figsize=(7.2, 7.2), dpi=180, sharex=True)

    axes[0].plot(data["time_ps"], data["temperature_K"], color="#1f77b4", linewidth=1.2)
    axes[0].axhline(300.0, color="#d62728", linestyle="--", linewidth=1.0)
    axes[0].set_ylabel("Temperature (K)")
    axes[0].set_title("Pure Water NPT Equilibration")

    axes[1].plot(data["time_ps"], data["pressure_bar"], color="#9467bd", linewidth=1.0)
    axes[1].axhline(1.0, color="#d62728", linestyle="--", linewidth=1.0)
    axes[1].set_ylabel("Pressure (bar)")

    axes[2].plot(data["time_ps"], data["volume_nm3"], color="#2ca02c", linewidth=1.2)
    axes[2].set_xlabel("Time (ps)")
    axes[2].set_ylabel("Volume (nm$^3$)")

    for axis in axes:
        axis.grid(True, linewidth=0.5, alpha=0.35)

    axes[-1].set_xlim(data["time_ps"].min(), data["time_ps"].max())
    fig.tight_layout()
    fig.savefig(png_path)
    plt.show()


def main() -> None:
    """Read the pure-water NPT XVG file and write the report figure."""

    data = read_npt_energy_xvg(XVG_PATH)
    plot_npt_energy(data, PNG_PATH)

    print(f"Parsed {len(data)} NPT energy samples from {XVG_PATH}")
    print(f"Mean temperature: {data['temperature_K'].mean():.2f} K")
    print(f"Mean pressure: {data['pressure_bar'].mean():.2f} bar")
    print(
        "Volume range: "
        f"{data['volume_nm3'].min():.2f}-{data['volume_nm3'].max():.2f} nm^3"
    )
    print(f"Wrote {PNG_PATH}")


if __name__ == "__main__":
    main()
