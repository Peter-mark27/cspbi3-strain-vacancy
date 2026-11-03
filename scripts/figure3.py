"""Three stacked atomistic-calibration panels from the frozen original data.

The 6 formation-energy, 5 uniform-strain NEB and 8 MD samples are unchanged.
(a) retains the adopted Ef0=2.51 eV, Af=4.500 eV relation for the original
    gradient specimen epsilon(x)=-0.08(1-x/L'), with L_atom=94.374 angstrom.
    This line is the adopted model relation, not a new free-intercept fit.
(b) uses fractional strain in Em=0.543-1.047*epsilon; the axis is percent.
(c) retains the original 0.001 ps/step, 2 ns-to-end MSD calibration.
    The display and CSV use SI diffusivities in m^2/s.
    D0 is the vacancy diffusion prefactor obtained from the vacancy MSD.

Panel (a) displays the original sampling positions in angstrom and the original
formation energies in the gradient specimen. Panel (b) displays the uniform
NEB strain in percent; formula coefficients multiply fractional strain.
The model links the two forms through epsilon(x)=-0.08(1-x/L_atom).
An equivalent position-dependent migration barrier is derived in the model,
not presented as additional atomistic data or as an independent validation.
Times New Roman/STIX, grey hollow points, red dashed curves and equal axes
sizes are retained. The enlarged vertical layout is intended for approximately
170 mm full-width placement.
Panel (c) combines Data, the fit, D0 and Ea in one aligned legend column.
The fit coefficient of determination R2=0.9923 is retained in the caption
information rather than printed inside the axes.
No synthetic points, smoothing, or parameter refitting are introduced.
Reads only the published frozen CSV in data/figures/.
Outputs go to reproduced_figures/ unless --output-dir is specified.
PNG is 600 dpi; PDF uses vector lines, markers and embedded TrueType fonts.

Suggested caption:
Atomistic calibration of the continuum model. (a) Vacancy formation energies
at six original positions in the 94.374 angstrom gradient specimen under
8% maximum compression, and the adopted formation-energy relation.
(b) The five uniform-strain NEB migration barriers and the fixed-intercept fit.
(c) Arrhenius fit to the eight vacancy diffusivities obtained from the
three-dimensional vacancy MSD. The log axis uses D/(1 m^2/s).
The apparent MD activation energy is distinct from the NEB migration barrier.
The Arrhenius coefficient of determination is R2=0.9923.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnchoredOffsetbox, DrawingArea
from matplotlib.text import Text
from matplotlib.transforms import Bbox
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "data" / "figures" / "figure3_data.csv"
DEFAULT_OUTPUT = ROOT / "reproduced_figures"

L_ATOM_ANGSTROM = 94.374
EPSILON0 = -0.08
EF0_EV = 2.51
AF_ADOPTED_EV = 4.500
EM0_EV = 0.543
KB_EV_PER_K = 8.617333262e-5

COLOR_DATA = "#424949"
COLOR_FIT = "#C53A32"

def source_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_snapshot() -> list[dict[str, str]]:
    with DATA_PATH.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    return rows


def panel_values(rows: list[dict], panel: str, x_column: str, y_column: str):
    selected = [row for row in rows if row["panel"] == panel]
    x = np.array([float(row[x_column]) for row in selected])
    y = np.array([float(row[y_column]) for row in selected])
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError(f"Non-finite snapshot values: {panel}")
    order = np.argsort(x)
    return x[order], y[order]


def calibration(rows: list[dict]):
    xa, ef = panel_values(rows, "formation", "x_over_L", "energy_eV")
    epsilon, em = panel_values(rows, "migration", "strain_fraction", "energy_eV")
    inverse_t, log_d = panel_values(
        rows, "arrhenius", "inverse_temperature_1000_per_K", "ln_D_V_SI"
    )
    if (len(xa), len(epsilon), len(inverse_t)) != (6, 5, 8):
        raise ValueError("Expected exactly 6 formation, 5 NEB, and 8 MD samples")
    af_fit = np.sum(EPSILON0 * (1-xa) * (ef-EF0_EV)) / np.sum((EPSILON0*(1-xa))**2)
    am = float(np.sum(epsilon*(em-EM0_EV)) / np.sum(epsilon**2))
    slope, intercept = np.polyfit(inverse_t, log_d, 1)
    activation = float(-slope * 1000 * KB_EV_PER_K)
    prefactor = float(np.exp(intercept))
    # This operation redraws the existing calibration, not a new parameter set.
    if not np.isclose(am, -1.047, rtol=0, atol=1e-9):
        raise ValueError("NEB data no longer reproduce the existing calibration")
    if not np.isclose(prefactor, 5.997569823971534e-9, rtol=1e-7, atol=0.0):
        raise ValueError("MSD data no longer reproduce the existing prefactor")
    if not np.isclose(activation, 0.5319959819896776, rtol=0, atol=1e-7):
        raise ValueError("MSD data no longer reproduce the existing activation energy")
    return xa, ef, epsilon, em, inverse_t, log_d, af_fit, am, slope, intercept, activation, prefactor




MASTER_SCALE = 2.0
LEGEND_FONT_SIZE = 11.5
FORMULA_FONT_SIZE = 15.0  # Enlarge only the equations in panels (a) and (b).
PDF_EXPORT_WIDTH = 12.80  # Balanced portrait layout; avoid stretching stacked panels.
PNG_EXPORT_WIDTH = (7680 + 0.25) / 600.0


def check_archived_points(rows):
    xa, ef = panel_values(rows, "formation", "position_angstrom", "energy_eV")
    epsilon, em = panel_values(rows, "migration", "strain_fraction", "energy_eV")
    expected = (
        (xa, [14.1561, 26.7392, 39.3224, 51.9055, 64.4887, 77.0719]),
        (ef, [2.169621, 2.241870, 2.306140, 2.372112, 2.436047, 2.492500]),
        (epsilon, [-0.02, -0.01, 0, 0.01, 0.02]),
        (em, [0.5601, 0.5530, 0.5431, 0.5317, 0.5184]),
    )
    for values, reference in expected:
        if not np.allclose(values, reference, rtol=0.0, atol=5e-13):
            raise ValueError("Snapshot no longer matches the original samples")


def paper_style():
    q = MASTER_SCALE
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Liberation Serif", "DejaVu Serif"],
        "mathtext.fontset": "stix", "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.labelsize": 18*q, "xtick.labelsize": 14*q,
        "ytick.labelsize": 14*q, "legend.fontsize": LEGEND_FONT_SIZE*q,
        "axes.linewidth": 1.8*q, "axes.grid": False,
        "axes.unicode_minus": True, "savefig.dpi": 600,
    })


def draw_pair(ax, x, y, xfit, yfit, fit_label, loc, extra_labels=()):
    q = MASTER_SCALE
    fit, = ax.plot(xfit, yfit, color=COLOR_FIT, lw=2.0*q,
                   linestyle=(0, (3.7, 1.6)), dash_capstyle="butt", zorder=2)
    data, = ax.plot(x, y, linestyle="none", marker="o", ms=7.3*q,
                    color=COLOR_DATA, mfc="white", mew=1.0*q, zorder=3)
    # Three complete dashes in the legend; the plotted curve is unchanged.
    fit_handle = Line2D([], [], color=COLOR_FIT, lw=2.0*q, dash_capstyle="butt")
    factor = (2.3*LEGEND_FONT_SIZE/2.0) / (3*3.7 + 2*1.6)
    fit_handle.set_linestyle((0, (3.7*factor, 1.6*factor)))
    handles = [data, fit_handle] + [Line2D([], [], linestyle="none") for _ in extra_labels]
    labels = ["Data", fit_label, *extra_labels]
    if extra_labels:
        legend = aligned_parameter_legend(ax, data, fit_handle, labels, loc)
        return data, fit, legend
    legend = ax.legend(handles, labels, loc=loc,
                      frameon=False, handlelength=2.3, handletextpad=0.55,
                      labelspacing=0.38, borderaxespad=0.8)
    legend.get_texts()[1].set_fontsize(FORMULA_FONT_SIZE*q)
    return data, fit, legend



def aligned_parameter_legend(ax, data, fit_handle, labels, loc):
    """Align complete symbol-and-text rows with the lower parameter block."""
    fontsize = LEGEND_FONT_SIZE*MASTER_SCALE
    renderer = ax.figure.canvas.get_renderer()
    widths = []
    for label in labels:
        probe = ax.text(0, 0, label, fontsize=fontsize)
        widths.append(probe.get_window_extent(renderer).width*72/ax.figure.dpi)
        probe.remove()
    handle_width = 2.3*fontsize
    text_x = handle_width + 0.55*fontsize
    pitch = 1.5*fontsize
    baseline0 = 0.35*fontsize
    area_width = max(2*text_x+max(widths[:2]), text_x+max(widths[2:]))
    area = DrawingArea(area_width, (len(labels)-1)*pitch+1.35*fontsize)
    texts = []
    for index, label in enumerate(labels):
        baseline = baseline0+(len(labels)-1-index)*pitch
        row_text_x = 2*text_x if index < 2 else text_x
        text = Text(row_text_x, baseline, label, fontsize=fontsize,
                    ha="left", va="baseline")
        area.add_artist(text)
        texts.append(text)
        center_y = baseline+0.3*fontsize
        if index == 0:
            area.add_artist(Line2D([text_x+handle_width/2], [center_y], linestyle="none",
                marker="o", ms=data.get_markersize(), color=COLOR_DATA,
                mfc="white", mew=data.get_markeredgewidth()))
        elif index == 1:
            line = Line2D([text_x, text_x+handle_width], [center_y, center_y],
                color=COLOR_FIT, lw=2.0*MASTER_SCALE, dash_capstyle="butt")
            factor = (2.3*LEGEND_FONT_SIZE/2.0)/(3*3.7+2*1.6)
            line.set_linestyle((0, (3.7*factor, 1.6*factor)))
            area.add_artist(line)
    legend = AnchoredOffsetbox(loc=loc, child=area, frameon=False, pad=0,
                               borderpad=0.8, prop={"size": fontsize})
    legend.figure3_texts = texts
    ax.add_artist(legend)
    return legend


def build_figure(rows):
    paper_style()
    xa, ef, epsilon, em, inverse_t, log_d, af_fit, am, slope, intercept, ea, d0 = calibration(rows)
    q = MASTER_SCALE
    fig = plt.figure(figsize=(12.9, 26.3))
    axes = [fig.add_axes([0.145, bottom, 0.835, 0.254])
            for bottom in (0.715, 0.392, 0.069)]
    for ax, label in zip(axes, ("(a)", "(b)", "(c)")):
        ax.tick_params(which="major", direction="in", length=4*q, width=1.35*q, pad=5*q)
        ax.tick_params(which="minor", direction="in", length=2.5*q, width=1.0*q)
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.text(-0.065, 1.018, label, transform=ax.transAxes,
                fontsize=17*q, fontweight="normal", ha="left", va="bottom")

    ax = axes[0]
    # Retain the measured sampling coordinates of the gradient specimen.
    x, _ = panel_values(rows, "formation", "position_angstrom", "energy_eV")
    xf = np.linspace(x.min(), x.max(), 400)
    yf = EF0_EV + AF_ADOPTED_EV*EPSILON0*(1-xf/L_ATOM_ANGSTROM)
    draw_pair(ax, x, ef, xf, yf,
              r"$E_\mathrm{f}=2.51-0.360\left(1-\frac{x}{L'}\right)$", "upper left")
    ax.set(xlim=(9, 82), ylim=(2.14, 2.535))
    ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.set_xlabel(r"$x$ ($\mathrm{\AA}$)", labelpad=1*q)
    ax.set_ylabel(r"$E_\mathrm{f}$ (eV)", labelpad=5*q)

    ax = axes[1]
    xf = np.linspace(epsilon.min(), epsilon.max(), 400)
    draw_pair(ax, epsilon*100, em, xf*100, EM0_EV+am*xf,
              r"$E_\mathrm{m}=0.543-1.047\varepsilon$", "upper right")
    ax.set(xlim=(-2.3, 2.3), ylim=(0.514, 0.5665))
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(0.01))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.set_xlabel(r"$\varepsilon$ (%)", labelpad=5*q)
    ax.set_ylabel(r"$E_\mathrm{m}$ (eV)", labelpad=5*q)

    ax = axes[2]
    exponent = int(np.floor(np.log10(d0)))
    mantissa = d0/10**exponent
    parameter_labels = (
        rf"$D_0={mantissa:.3f}\times10^{{{exponent}}}\ "
        rf"\mathrm{{m}}^2\,\mathrm{{s}}^{{-1}}$",
        rf"$E_\mathrm{{a}}={ea:.3f}\ \mathrm{{eV}}$",
    )
    xf = np.linspace(inverse_t.min(), inverse_t.max(), 400)
    draw_pair(ax, inverse_t, log_d, xf, slope*xf+intercept,
              "Arrhenius fit", "upper right", extra_labels=parameter_labels)
    ax.set(xlim=(1.425, 2.045),
           ylim=(-22.5-np.log(1e4), -18.5-np.log(1e4)))
    ax.xaxis.set_major_locator(MultipleLocator(0.1))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.set_xlabel(r"$1000/T$ ($\mathrm{K}^{-1}$)", labelpad=5*q)
    ax.set_ylabel(r"$\ln D$ ($\mathrm{m}^2\,\mathrm{s}^{-1}$)", labelpad=5*q)
    return fig, axes


def export_boxes(fig):
    padding = 0.12*MASTER_SCALE
    for _ in range(8):
        fig.canvas.draw()
        tight = fig.get_tightbbox(fig.canvas.get_renderer())
        correction = PDF_EXPORT_WIDTH-(tight.width+2*padding)
        if abs(correction) < 1e-7:
            break
        fig.set_size_inches(fig.get_figwidth()+correction/0.92,
                            fig.get_figheight(), forward=True)
    fig.canvas.draw()
    tight = fig.get_tightbbox(fig.canvas.get_renderer())
    origin = tight.x0-padding, tight.y0-padding
    height = tight.height+2*padding
    return (Bbox.from_bounds(*origin, PNG_EXPORT_WIDTH, height),
            Bbox.from_bounds(*origin, PDF_EXPORT_WIDTH, height))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    rows = read_snapshot()
    check_archived_points(rows)
    original_hash = source_hash(DATA_PATH)
    fig, axes = build_figure(rows)
    if [len(ax.lines[1].get_xdata()) for ax in axes] != [6, 5, 8]:
        raise AssertionError("Every original point must be displayed")
    png_box, pdf_box = export_boxes(fig)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pdf = args.output_dir / "figure3.pdf"
    png = args.output_dir / "figure3.png"
    fig.savefig(pdf, bbox_inches=pdf_box, facecolor="white",
                metadata={"Title": "Atomistic calibration of the continuum model"})
    fig.savefig(png, bbox_inches=png_box, dpi=600, facecolor="white")
    plt.close(fig)
    assert source_hash(DATA_PATH) == original_hash
    print(f"All 19 original points retained; unchanged CSV SHA256: {original_hash}")
    print(f"PDF width: {PDF_EXPORT_WIDTH:.12f} in; PNG: 600 dpi")
    print(f"WROTE {pdf}\nWROTE {png}")


if __name__ == "__main__":
    main()
