"""使用 pandas 绘制任务四的 Rg(t) 与 Ree(t) 组合对比图。

本脚本读取两个体系在 ``ANALYSIS`` 目录中已经由 GROMACS 生成的 XVG 文件：

* ``systems/water/ANALYSIS/rg.xvg``
* ``systems/water/ANALYSIS/ree.xvg``
* ``systems/water_ethanol/ANALYSIS/rg.xvg``
* ``systems/water_ethanol/ANALYSIS/ree.xvg``

其中 ``rg.xvg`` 来自 ``gmx gyrate``，第二列为总回转半径 Rg；
``ree.xvg`` 来自 ``gmx distance``，第二列为末端距 Ree。脚本使用
``pandas`` 完成 XVG 解析、统计汇总和滑动平均，并生成：

* ``report/static/task4_rg_ree_comparison.png``
* ``report/static/task4_rg_ree_summary.csv``

项目脚本约定是不使用命令行参数，因此路径均在脚本中固定。请从项目根目录运行：

    python scripts/task4/plot_rg_ree.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
STATIC_DIR = ROOT / "report" / "static"
COMBINED_PNG_PATH = STATIC_DIR / "task4_rg_ree_comparison.png"
SUMMARY_PATH = STATIC_DIR / "task4_rg_ree_summary.csv"
ROLLING_WINDOW = 21

SYSTEMS = {
    "water": {
        "label": "Pure water",
        "analysis_dir": ROOT / "systems" / "water" / "ANALYSIS",
        "color": "#2563eb",
    },
    "water_ethanol": {
        "label": "Water/ethanol",
        "analysis_dir": ROOT / "systems" / "water_ethanol" / "ANALYSIS",
        "color": "#d97706",
    },
}

OBSERVABLES = {
    "rg": {
        "filename": "rg.xvg",
        "column": "rg_nm",
        "ylabel": r"$R_g$ (nm)",
        "title": "Radius of Gyration",
    },
    "ree": {
        "filename": "ree.xvg",
        "column": "ree_nm",
        "ylabel": r"$R_{ee}$ (nm)",
        "title": "End-to-End Distance",
    },
}


def count_xvg_header_lines(xvg_path: Path) -> int:
    """统计 XVG 文件开头需要跳过的说明行数。

    GROMACS 的 XVG 文件通常先写入以 ``#`` 和 ``@`` 开头的元信息，
    然后才是纯数值表。这里逐行扫描，遇到第一行数值数据时停止，
    后续交给 ``pandas.read_csv`` 读取。

    Args:
        xvg_path: 待读取的 XVG 文件路径。

    Returns:
        文件开头元信息行和空行的总数。

    Raises:
        FileNotFoundError: 当 XVG 文件不存在时抛出。
    """

    if not xvg_path.exists():
        raise FileNotFoundError(f"Missing XVG file: {xvg_path}")

    skiprows = 0
    for line in xvg_path.read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.lstrip()
        if stripped.startswith(("#", "@")) or not stripped:
            skiprows += 1
            continue
        break

    return skiprows


def read_observable_xvg(xvg_path: Path, value_column: str) -> pd.DataFrame:
    """用 pandas 读取 XVG 文件的时间列和第二列目标数值。

    对 ``rg.xvg`` 而言，文件中还包含三个方向分量，但任务四需要比较的是
    总 Rg，因此只读取第一列时间和第二列总 Rg。对 ``ree.xvg`` 而言，
    第二列即为 ``gmx distance`` 计算出的末端距。

    Args:
        xvg_path: 待读取的 XVG 文件路径。
        value_column: 第二列数值在返回表中的列名，例如 ``rg_nm``。

    Returns:
        包含 ``time_ps``、``time_ns`` 和目标数值列的 DataFrame。

    Raises:
        ValueError: 当文件中没有读到数值数据时抛出。
    """

    skiprows = count_xvg_header_lines(xvg_path)
    data = pd.read_csv(
        xvg_path,
        sep=r"\s+",
        skiprows=skiprows,
        names=["time_ps", value_column],
        usecols=[0, 1],
        engine="python",
    )
    if data.empty:
        raise ValueError(f"No numeric samples were parsed from {xvg_path}")

    data["time_ns"] = data["time_ps"] / 1000.0
    return data


def load_task4_datasets() -> dict[str, dict[str, pd.DataFrame]]:
    """读取两个体系的 Rg 与 Ree 数据。

    Returns:
        两级字典。第一级键为体系名，第二级键为 ``rg`` 或 ``ree``，
        值为包含时间和目标数值的 DataFrame。
    """

    datasets: dict[str, dict[str, pd.DataFrame]] = {}
    for system_name, system_config in SYSTEMS.items():
        analysis_dir = system_config["analysis_dir"]
        datasets[system_name] = {}

        for observable_name, observable_config in OBSERVABLES.items():
            xvg_path = analysis_dir / observable_config["filename"]
            value_column = observable_config["column"]
            datasets[system_name][observable_name] = read_observable_xvg(
                xvg_path,
                value_column,
            )

    return datasets


def summarize_observable(data: pd.DataFrame, value_column: str) -> dict[str, float]:
    """计算一个指标的全程、前半程和后半程统计量。

    任务书要求对 ``Rg`` 做简单收敛性检查。这里对 ``Rg`` 和 ``Ree`` 都采用
    同一套透明指标：全程均值、标准差、最小值、最大值，以及后半程均值
    与前半程均值的差值。后半程和前半程越接近，说明该指标的整体漂移越小。

    Args:
        data: 包含 ``time_ns`` 和目标数值列的 DataFrame。
        value_column: 需要统计的目标数值列名。

    Returns:
        可直接写入 CSV 的统计量字典。
    """

    midpoint_ns = data["time_ns"].max() / 2.0
    first_half = data[data["time_ns"] <= midpoint_ns]
    second_half = data[data["time_ns"] > midpoint_ns]

    return {
        "n_samples": float(len(data)),
        "time_start_ns": data["time_ns"].min(),
        "time_end_ns": data["time_ns"].max(),
        "mean_nm": data[value_column].mean(),
        "std_nm": data[value_column].std(),
        "min_nm": data[value_column].min(),
        "max_nm": data[value_column].max(),
        "first_half_mean_nm": first_half[value_column].mean(),
        "second_half_mean_nm": second_half[value_column].mean(),
        "second_minus_first_nm": second_half[value_column].mean() - first_half[value_column].mean(),
    }


def build_summary_table(datasets: dict[str, dict[str, pd.DataFrame]]) -> pd.DataFrame:
    """把两个体系、两个指标的统计量整理为一张 pandas 表。

    Args:
        datasets: ``load_task4_datasets`` 返回的数据字典。

    Returns:
        每一行对应一个体系的一个指标，列为统计量。
    """

    rows: list[dict[str, object]] = []
    for system_name, observable_map in datasets.items():
        for observable_name, data in observable_map.items():
            value_column = OBSERVABLES[observable_name]["column"]
            row: dict[str, object] = {
                "system": system_name,
                "label": SYSTEMS[system_name]["label"],
                "observable": observable_name,
            }
            row.update(summarize_observable(data, value_column))
            rows.append(row)

    return pd.DataFrame(rows)


def plot_task4_combined_figure(
    datasets: dict[str, dict[str, pd.DataFrame]],
    summary: pd.DataFrame,
    png_path: Path,
) -> None:
    """把 Rg(t) 和 Ree(t) 画在同一张图的两个子图中。

    上子图为回转半径 Rg，下子图为末端距 Ree。每个子图中浅色线表示原始
    数据，深色线表示 pandas 滑动平均；滑动平均只用于辅助观察趋势，
    统计表中的均值和标准差仍然来自原始数据。

    Args:
        datasets: 两个体系的 Rg/Ree 数据。
        summary: 统计汇总表，用于在图例中标出全程均值。
        png_path: 输出 PNG 文件路径。
    """

    import matplotlib.pyplot as plt

    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 1, figsize=(7.8, 7.0), dpi=180, sharex=True)

    min_time = None
    max_time = None
    for axis, observable_name in zip(axes, OBSERVABLES):
        observable_config = OBSERVABLES[observable_name]
        value_column = observable_config["column"]

        for system_name, observable_map in datasets.items():
            data = observable_map[observable_name]
            system_config = SYSTEMS[system_name]
            summary_row = summary[
                (summary["system"] == system_name)
                & (summary["observable"] == observable_name)
            ].iloc[0]

            rolling = data[value_column].rolling(
                window=ROLLING_WINDOW,
                center=True,
                min_periods=1,
            ).mean()

            # 原始曲线显示短时波动，滑动平均显示更容易比较的构象变化趋势。
            axis.plot(
                data["time_ns"],
                data[value_column],
                color=system_config["color"],
                linewidth=0.6,
                alpha=0.22,
            )
            axis.plot(
                data["time_ns"],
                rolling,
                color=system_config["color"],
                linewidth=1.7,
                label=f"{system_config['label']}",
            )

            current_min_time = data["time_ns"].min()
            current_max_time = data["time_ns"].max()
            min_time = current_min_time if min_time is None else min(min_time, current_min_time)
            max_time = current_max_time if max_time is None else max(max_time, current_max_time)

        axis.set_title(observable_config["title"])
        axis.set_ylabel(observable_config["ylabel"])
        axis.grid(True, linewidth=0.5, alpha=0.35)
        axis.legend(frameon=False, loc="best")

    axes[-1].set_xlabel("Time (ns)")
    axes[-1].set_xlim(min_time, max_time)
    fig.tight_layout()
    fig.savefig(png_path)
    plt.show()


def main() -> None:
    """读取 XVG、绘制组合图、写出统计 CSV，并在终端打印摘要。"""

    datasets = load_task4_datasets()
    summary = build_summary_table(datasets)

    STATIC_DIR.mkdir(parents=True, exist_ok=True)
    summary.to_csv(SUMMARY_PATH, index=False, float_format="%.6f")
    plot_task4_combined_figure(datasets, summary, COMBINED_PNG_PATH)

    print(f"Wrote {COMBINED_PNG_PATH}")
    print(f"Wrote {SUMMARY_PATH}")

    for _, row in summary.iterrows():
        print(
            f"[{row['label']}] {row['observable']}: "
            f"mean={row['mean_nm']:.3f} nm, std={row['std_nm']:.3f} nm, "
            f"second-half minus first-half={row['second_minus_first_nm']:.3f} nm"
        )


if __name__ == "__main__":
    main()
