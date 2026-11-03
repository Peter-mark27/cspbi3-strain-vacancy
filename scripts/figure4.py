"""Create the final mesh/time-step independence figure from source CSV data.

Panels:
  (a) Final concentration profiles for three mapped-mesh spacings.
  (b) Final concentration profiles for three maximum time steps.

No difference curves, interpolation, smoothing, or synthetic samples are used.
Reads only the published frozen CSV in data/figures/.
Outputs go to reproduced_figures/ unless --output-dir is specified.

Suggested caption:
Mesh and time-step independence of the steady concentration distribution.
(a) Concentration profiles obtained with three mesh spacings. (b)
Concentration profiles obtained with three maximum time steps. The overlap of
the profiles demonstrates that the adopted spatial and temporal resolutions
are sufficient for the subsequent simulations.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path
from types import MethodType

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FormatStrFormatter
from matplotlib.lines import Line2D
from matplotlib.offsetbox import TextArea
from matplotlib.transforms import Bbox


ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "data" / "figures" / "figure4_data.csv"
DEFAULT_OUTPUT = ROOT / "reproduced_figures"

MASTER_SCALE = 2.0
# Match the approved Figure 4 export widths; its PNG/PDF tight bounds differ
# slightly because their font renderers measure text at different resolutions.
PDF_EXPORT_WIDTH = 20.717215711805554
PNG_EXPORT_WIDTH = (12429 + 0.25) / 600
LEGEND_FONT_SIZE = 11.0
LEGEND_HANDLE_LENGTH = 2.3
PROFILE_LINE_WIDTH = 2.35
COLORS = ("#3C5488", "#E64B35", "#00A087")
LINESTYLES = ("--", "-", ":")

CURVES = {
    "mesh": (
        ("h_coarse_nm", r"$h=1.304\ \mathrm{nm}$"),
        ("h_base_nm", r"$h=0.654\ \mathrm{nm}$"),
        ("h_fine_nm", r"$h=0.327\ \mathrm{nm}$"),
    ),
    "timestep": (
        ("dtmax_base", r"$\Delta t_{\max}=5.206\times10^{-5}\ \mathrm{s}$"),
        ("half_dtmax", r"$\Delta t_{\max}=2.603\times10^{-5}\ \mathrm{s}$"),
        ("quarter_dtmax", r"$\Delta t_{\max}=1.302\times10^{-5}\ \mathrm{s}$"),
    ),
}

def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_snapshot() -> list[dict[str, str]]:
    if not DATA_PATH.exists():
        raise FileNotFoundError(
            f"Missing published figure data: {DATA_PATH}"
        )
    with DATA_PATH.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError("Empty figure4 data snapshot")
    return rows


def extract_curve(
    rows: list[dict[str, str]], panel: str, curve_id: str
) -> tuple[np.ndarray, np.ndarray]:
    selected = [
        row
        for row in rows
        if row["panel"] == panel and row["curve_id"] == curve_id
    ]
    if not selected:
        raise ValueError(f"Missing curve {panel}/{curve_id}")
    x = np.asarray([float(row["x_nm"]) for row in selected])
    y = np.asarray([float(row["c_over_c0"]) for row in selected])
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError(f"Non-finite snapshot values: {panel}/{curve_id}")
    order = np.argsort(x, kind="stable")
    x, y = x[order], y[order]
    if np.any(np.diff(x) <= 0.0):
        raise ValueError(f"Non-increasing x coordinates: {panel}/{curve_id}")
    return x, y


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Liberation Serif", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.unicode_minus": False,
            "axes.linewidth": 1.30*MASTER_SCALE,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "savefig.dpi": 600,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def style_axis(ax: plt.Axes) -> None:
    ax.tick_params(
        "both",
        which="major",
        direction="in",
        length=4.0*MASTER_SCALE,
        width=1.1*MASTER_SCALE,
        labelsize=14*MASTER_SCALE,
        pad=5*MASTER_SCALE,
        top=False,
        right=False,
    )
    ax.tick_params(
        "both",
        which="minor",
        direction="in",
        length=2.5*MASTER_SCALE,
        width=0.85*MASTER_SCALE,
        top=False,
        right=False,
    )
    ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.yaxis.set_major_locator(MultipleLocator(0.05))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    for spine in ax.spines.values():
        spine.set_linewidth(1.30*MASTER_SCALE)


def legend_handle(color: str, linestyle: str, label: str) -> Line2D:
    """Keep full dash/dot samples inside the original short legend handles."""
    handle = Line2D([], [], color=color, linewidth=PROFILE_LINE_WIDTH*MASTER_SCALE,
                    linestyle=linestyle, label=label, dash_capstyle="butt")
    if linestyle != "-":
        on, off, count = (3.7, 1.6, 3) if linestyle == "--" else (1.0, 1.65, 6)
        length = LEGEND_HANDLE_LENGTH*LEGEND_FONT_SIZE / PROFILE_LINE_WIDTH
        factor = length / (count*on + (count-1)*off)
        handle.set_linestyle((0, (on*factor, off*factor)))
    return handle


def align_legend_rows(axes: np.ndarray) -> None:
    """Use common text ascent/descent so superscripts do not shift legend rows.

    Metrics are evaluated for each renderer, preserving alignment in PNG and PDF.
    Labels, font sizes and dash patterns remain native legend artists.
    """
    areas = [row.get_children()[1]
             for ax in axes
             for column in ax.get_legend()._legend_handle_box.get_children()
             for row in column.get_children()]

    def common_bbox(area: TextArea, renderer) -> Bbox:
        bounds = [TextArea.get_bbox(item, renderer) for item in areas]
        own = bounds[areas.index(area)]
        return Bbox.from_extents(own.x0, min(b.y0 for b in bounds),
                                 own.x1, max(b.y1 for b in bounds))

    for area in areas:
        area.get_bbox = MethodType(common_bbox, area)


def export_boxes(fig: plt.Figure) -> tuple[Bbox, Bbox]:
    """Fit the layout to the approved export width without stretching artwork."""
    pad = 0.14*MASTER_SCALE
    for _ in range(6):
        fig.canvas.draw()
        tight = fig.get_tightbbox(fig.canvas.get_renderer())
        correction = PDF_EXPORT_WIDTH - (tight.width + 2*pad)
        if abs(correction) < 1e-7:
            break
        fig.set_size_inches(fig.get_figwidth() + correction/0.9,
                            fig.get_figheight(), forward=False)
    fig.canvas.draw()
    tight = fig.get_tightbbox(fig.canvas.get_renderer())
    origin = (tight.x0-pad, tight.y0-pad)
    height = tight.height + 2*pad
    return (Bbox.from_bounds(*origin, PNG_EXPORT_WIDTH, height),
            Bbox.from_bounds(*origin, PDF_EXPORT_WIDTH, height))


def draw(rows: list[dict[str, str]], output_dir: Path) -> None:
    configure_style()
    # Increase panel height while matching the approved Figure 4 overall width.
    fig, axes = plt.subplots(1, 2, figsize=(10.4*MASTER_SCALE, 4.20*MASTER_SCALE))

    for ax, panel, panel_label in zip(
        axes, ("mesh", "timestep"), ("(a)", "(b)")
    ):
        values: list[np.ndarray] = []
        for (curve_id, label), color, linestyle in zip(
            CURVES[panel], COLORS, LINESTYLES
        ):
            x_nm, concentration = extract_curve(rows, panel, curve_id)
            values.append(concentration)
            ax.plot(
                x_nm,
                concentration,
                color=color,
                linestyle=linestyle,
                linewidth=PROFILE_LINE_WIDTH*MASTER_SCALE,
                solid_capstyle="butt",
                dash_capstyle="butt",
                label=label,
                zorder={"-": 1, "--": 2, ":": 3}[linestyle],
            )

        all_values = np.concatenate(values)
        span = float(all_values.max() - all_values.min())
        ax.set_xlim(-6.0, 306.0)
        ax.set_ylim(
            float(all_values.min() - max(0.08 * span, 0.02)),
            float(all_values.max() + max(0.08 * span, 0.02)),
        )
        ax.set_xlabel(r"$x$ (nm)", fontsize=19*MASTER_SCALE, labelpad=3*MASTER_SCALE)
        ax.set_ylabel(r"$c/c_0$", fontsize=19*MASTER_SCALE, labelpad=4*MASTER_SCALE)
        ax.legend(
            handles=[legend_handle(color, linestyle, label)
                     for (_, label), color, linestyle in zip(CURVES[panel], COLORS, LINESTYLES)],
            loc="upper right",
            frameon=False,
            fontsize=LEGEND_FONT_SIZE*MASTER_SCALE,
            handlelength=LEGEND_HANDLE_LENGTH,
            handletextpad=0.5,
            labelspacing=0.32,
            borderaxespad=0.5,
        )
        ax.text(
            -0.115,
            1.015,
            panel_label,
            transform=ax.transAxes,
            fontsize=15.5*MASTER_SCALE,
            ha="left",
            va="bottom",
            clip_on=False,
        )
        style_axis(ax)

    fig.subplots_adjust(left=0.085, right=0.98, bottom=0.185, top=0.92, wspace=0.29)
    align_legend_rows(axes)
    png_box, pdf_box = export_boxes(fig)
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = output_dir / "figure4"
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches=png_box, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches=pdf_box, facecolor="white")
    plt.close(fig)

    for panel in ("mesh", "timestep"):
        reference_id = CURVES[panel][-1][0]
        _, reference = extract_curve(rows, panel, reference_id)
        for curve_id, _ in CURVES[panel][:-1]:
            _, values = extract_curve(rows, panel, curve_id)
            maximum = np.max(np.abs((values - reference) / reference)) * 100.0
            print(f"{panel}/{curve_id} vs {reference_id}: max difference={maximum:.6g}%")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    original_hash = sha256(DATA_PATH)
    draw(read_snapshot(), args.output_dir)
    assert sha256(DATA_PATH) == original_hash
    print(f"Unchanged CSV SHA256: {original_hash}")


if __name__ == "__main__":
    main()
