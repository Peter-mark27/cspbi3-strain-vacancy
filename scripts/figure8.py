"""Create a 3-by-2 paper figure for the three parameter studies.

Each row is one parameter group: (a,b) strain, (c,d) temperature, (e,f) initial
concentration. Left: full COMSOL concentration profiles. Right: COMSOL data
markers versus independently solved nonlinear PB theory, with matching
coordinates, limits and case colors. No interpolation or fitting is used.

Suggested caption: Steady concentration profiles and comparison with the
full nonlinear PB equilibrium model. (a,b) Strain scan at T=600 K and
c0/c_ref=1; (c,d) temperature scan at |epsilon0|=8% and c0/c_ref=1; (e,f) initial
concentration scan at T=600 K and |epsilon0|=8%. Left panels show complete
COMSOL profiles; right panels show selected COMSOL markers and full nonlinear
PB numerical solutions. The colors identify the same cases across each row.

Redraws the packaged frozen data snapshot without modifying it.
Exports a 600-dpi PNG and vector PDF.
"""

from __future__ import annotations

import csv
import argparse
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FormatStrFormatter
from matplotlib.transforms import Bbox


ROOT = Path(__file__).resolve().parents[1]
SNAPSHOT = ROOT / "data" / "figures" / "figure8_data.csv"
DEFAULT_OUTPUT = ROOT / "reproduced_figures"
MASTER_SCALE = 2.0
# Keep both legends' handles equally short, fitting three complete theory dashes.
LEGEND_HANDLE_LENGTH = 2.3
THEORY_LEGEND_DASH_SCALE = LEGEND_HANDLE_LENGTH*10 / ((3*3.2 + 2*2.2)*1.9)
# Fixed export bounds for Figure 8.
EXPORT_BOX = Bbox.from_bounds(.468, .8576344401041663, 16.264, 21.82527261594722)

COLORS = ("#7E57C2", "#3C5488", "#00A087", "#E64B35")
TEMP_COLORS = ("#3C5488", "#00A087", "#E69F00", "#E64B35")
PANEL_ORDER = ("strain", "temperature", "concentration")
CASE_ORDER = {
    "strain": ("600K_eps05", "600K_eps06", "600K_eps07", "600K_eps08"),
    "temperature": ("450K_eps08", "500K_eps08", "550K_eps08", "600K_eps08"),
    "concentration": ("c0r0.04", "c0r0.16", "c0r0.40", "c0r1.00"),
}
LABELS = {
    ("strain", "600K_eps05"): r"$|\varepsilon_0|=5\%$",
    ("strain", "600K_eps06"): r"$|\varepsilon_0|=6\%$",
    ("strain", "600K_eps07"): r"$|\varepsilon_0|=7\%$",
    ("strain", "600K_eps08"): r"$|\varepsilon_0|=8\%$",
    ("temperature", "450K_eps08"): r"$T=450\,\mathrm{K}$",
    ("temperature", "500K_eps08"): r"$T=500\,\mathrm{K}$",
    ("temperature", "550K_eps08"): r"$T=550\,\mathrm{K}$",
    ("temperature", "600K_eps08"): r"$T=600\,\mathrm{K}$",
    ("concentration", "c0r0.04"): r"$c_0/c_{\mathrm{ref}}=0.04$",
    ("concentration", "c0r0.16"): r"$c_0/c_{\mathrm{ref}}=0.16$",
    ("concentration", "c0r0.40"): r"$c_0/c_{\mathrm{ref}}=0.40$",
    ("concentration", "c0r1.00"): r"$c_0/c_{\mathrm{ref}}=1.00$",
}
PANEL_COLORS = {
    "strain": dict(zip(CASE_ORDER["strain"], COLORS)),
    "temperature": dict(zip(CASE_ORDER["temperature"], TEMP_COLORS)),
    "concentration": dict(zip(CASE_ORDER["concentration"], COLORS)),
}


def load_table(path: Path) -> dict[tuple[str, str], dict[str, np.ndarray | str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    rows: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            rows[(row["panel"], row["case_id"])].append(row)

    profiles: dict[tuple[str, str], dict[str, np.ndarray | str]] = {}
    for key, case_rows in rows.items():
        case_rows.sort(key=lambda item: float(item["x_nm"]))
        x_nm = np.array([float(item["x_nm"]) for item in case_rows])
        comsol = np.array([float(item["comsol_c_over_c0"]) for item in case_rows])
        theory = np.array([float(item["theory_c_over_c0"]) for item in case_rows])
        if x_nm.size != 532 or np.any(np.diff(x_nm) <= 0):
            raise ValueError(f"Unexpected profile grid for {key}")
        source = case_rows[0]["comsol_source_file"]
        # Historical source labels are metadata; the frozen values are self-contained.
        profiles[key] = {
            "x_nm": x_nm,
            "comsol": comsol,
            "theory": theory,
            "source": source,
        }
    return profiles


def marker_indices(x_nm: np.ndarray) -> np.ndarray:
    requested = np.array(
        [0, 3, 8, 16, 30, 45, 60, 80, 100, 125, 150, 175, 200,
         220, 240, 255, 270, 284, 292, 297, 300],
        dtype=float,
    )
    return np.unique([np.abs(x_nm - value).argmin() for value in requested])


def style_axes(ax: plt.Axes) -> None:
    ax.tick_params(
        "both", which="major", direction="in", length=4.5*MASTER_SCALE, width=1.10*MASTER_SCALE,
        labelsize=13.5*MASTER_SCALE, pad=3*MASTER_SCALE, top=False, right=False,
    )
    ax.tick_params(
        "both", which="minor", direction="in", length=2.3*MASTER_SCALE, width=.80*MASTER_SCALE,
        top=False, right=False,
    )
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    for spine in ax.spines.values():
        spine.set_linewidth(1.30*MASTER_SCALE)


def limits_for_panel(
    profiles: dict[tuple[str, str], dict[str, np.ndarray | str]],
    panel: str,
) -> tuple[float, float]:
    arrays = []
    for case_id in CASE_ORDER[panel]:
        profile = profiles[(panel, case_id)]
        arrays.extend((profile["comsol"], profile["theory"]))
    values = np.concatenate(arrays)
    span = float(values.max() - values.min())
    return float(values.min() - 0.055 * span), float(values.max() + 0.055 * span)


def main(output_dir: Path = DEFAULT_OUTPUT) -> None:
    profiles = load_table(SNAPSHOT)
    expected = {(panel, case_id) for panel in PANEL_ORDER for case_id in CASE_ORDER[panel]}
    if set(profiles) != expected:
        raise ValueError(f"Unexpected case set: {set(profiles) ^ expected}")
    data_path = SNAPSHOT

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Liberation Serif", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.unicode_minus": False,
            "axes.linewidth": 1.30*MASTER_SCALE,
            "font.size": 13.5*MASTER_SCALE,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "savefig.dpi": 600,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, axes = plt.subplots(3, 2, figsize=(8.8*MASTER_SCALE,11.5*MASTER_SCALE))

    for row, panel in enumerate(PANEL_ORDER):
        numerical = axes[row, 0]
        comparison = axes[row, 1]
        cases = CASE_ORDER[panel]

        for case_id in cases:
            profile = profiles[(panel, case_id)]
            x_nm = profile["x_nm"]
            comsol = profile["comsol"]
            theory = profile["theory"]
            color = PANEL_COLORS[panel][case_id]
            numerical.plot(
                x_nm, comsol, color=color, linewidth=2.0*MASTER_SCALE,
                solid_capstyle="butt", label=LABELS[(panel, case_id)],
            )
            pick = marker_indices(x_nm)
            comparison.plot(
                x_nm, theory, color=color, linewidth=1.9*MASTER_SCALE,
                linestyle=(0, (3.2, 2.2)), dash_capstyle="butt", zorder=1,
            )
            comparison.plot(
                x_nm[pick], comsol[pick], linestyle="None", marker="o",
                markersize=5.6*MASTER_SCALE, markerfacecolor=color, markeredgecolor="black",
                markeredgewidth=0.55*MASTER_SCALE, zorder=2,
            )

        numerical.legend(
            loc="upper right", frameon=False, fontsize=10*MASTER_SCALE,
            handlelength=LEGEND_HANDLE_LENGTH, handletextpad=0.48, labelspacing=0.28,
            borderaxespad=0.42,
        )
        handles = [
            Line2D([], [], color=PANEL_COLORS[panel][case_id], linewidth=2.0*MASTER_SCALE,
                   label=LABELS[(panel,case_id)])
            for case_id in cases
        ]
        handles.extend([
            Line2D([], [], linestyle="None", marker="o", markersize=5.6*MASTER_SCALE,
                   markerfacecolor="0.18", markeredgecolor="black",
                   markeredgewidth=0.55*MASTER_SCALE, label="Data"),
            Line2D([], [], color="0.18", linewidth=1.9*MASTER_SCALE,
                   linestyle=(0, (3.2*THEORY_LEGEND_DASH_SCALE,
                                  2.2*THEORY_LEGEND_DASH_SCALE)),
                   dash_capstyle="butt", label="Theory"),
        ])
        comparison.legend(
            handles=handles, loc="upper right", frameon=False, fontsize=10*MASTER_SCALE,
            handlelength=LEGEND_HANDLE_LENGTH, handletextpad=0.48, labelspacing=0.23,
            borderaxespad=0.4,
        )

        y_limits = limits_for_panel(profiles, panel)
        for ax in (numerical, comparison):
            # Five-percent padding keeps endpoint curves and markers clear of the frame.
            ax.set_xlim(-15.0, 315.0)
            ax.set_ylim(*y_limits)
            ax.set_xticks(np.arange(0,301,50))
            ax.get_xticklabels()[-1].set_horizontalalignment('center')
            ax.set_xlabel(r"$x$ (nm)", fontsize=18*MASTER_SCALE, labelpad=3*MASTER_SCALE)
            ax.set_ylabel(r"$c/c_0$", fontsize=18*MASTER_SCALE, labelpad=4*MASTER_SCALE)
            style_axes(ax)
            ax.yaxis.set_major_locator(MultipleLocator(.2 if panel=='concentration' else .05))
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f' if panel=='concentration' else '%.2f'))
        for col, ax in enumerate((numerical, comparison)):
            ax.text(
                -.20, 1.035, f'({"abcdef"[2 * row + col]})', transform=ax.transAxes,
                ha="left", va="bottom", fontsize=18*MASTER_SCALE, clip_on=False,
            )

    fig.subplots_adjust(
        left=.105, right=.945, bottom=.135, top=.925,
        wspace=.38, hspace=.46,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = output_dir / "figure8"
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches=EXPORT_BOX, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches=EXPORT_BOX, facecolor="white")
    plt.close(fig)

    print(f"WROTE {stem.with_suffix('.png')}")
    print(f"WROTE {stem.with_suffix('.pdf')}")
    print(f"READ {data_path}")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUTPUT)
    main(output_dir=parser.parse_args().output_dir)
