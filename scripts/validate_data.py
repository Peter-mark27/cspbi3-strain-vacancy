"""Validate internal scientific consistency of the published frozen data.

Uses only the Python standard library, resolves paths from this script, and
writes nothing. This checks the released tables and their recorded provenance;
it does not rerun COMSOL or reconstruct unbundled raw simulation exports.
If SHA256SUMS is present, its non-generated file entries are verified as well.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_ROWS = {3: 19, 4: 2760, 5: 286000, 6: 5520,
                 7: 8236, 8: 6384, 9: 96, 10: 529}
STAGES = ("initial", "early", "late", "steady")
GENERATED_DIRS = {"reproduced_figures", "generated", "__pycache__", ".git"}


class ValidationError(Exception):
    """A missing file or inconsistent scientific record."""


def require(condition, message):
    if not condition:
        raise ValidationError(message)


def near(actual, expected, message, *, rel=1e-11, absolute=1e-14):
    require(math.isfinite(actual) and math.isfinite(expected),
            f"{message}: non-finite value")
    require(math.isclose(actual, expected, rel_tol=rel, abs_tol=absolute),
            f"{message}: {actual:.17g} != {expected:.17g}")


def number(row, column):
    try:
        value = float(row[column])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValidationError(f"Invalid numeric column {column!r}") from exc
    require(math.isfinite(value), f"Non-finite value in {column}")
    return value


def rows_at(relative_path):
    path = ROOT / relative_path
    require(path.is_file(), f"Missing released file: {relative_path}")
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    require(bool(rows), f"Empty released table: {relative_path}")
    return rows


def group_by(rows, *columns):
    result = defaultdict(list)
    for row in rows:
        result[tuple(row[column] for column in columns)].append(row)
    return result


def check_counts(figures):
    for figure, expected in EXPECTED_ROWS.items():
        require(len(figures[figure]) == expected,
                f"Figure {figure}: expected {expected} rows, found {len(figures[figure])}")
    require(Counter(r["panel"] for r in figures[3]) ==
            {"formation": 6, "migration": 5, "arrhenius": 8},
            "Figure 3: expected 6 formation, 5 migration and 8 Arrhenius samples")
    groups = group_by(figures[4], "panel", "curve_id")
    require(len(groups) == 6 and all(len(g) == 460 for g in groups.values()),
            "Figure 4: expected six 460-point convergence curves")
    print("PASS: all eight figure CSV row counts; atomistic and convergence sample counts")


def check_benchmark(figures):
    snapshots = rows_at("data/benchmark/paper_flux_snapshot_times.csv")
    require([r["stage"] for r in snapshots] == list(STAGES),
            "Benchmark snapshot stages or their order changed")
    indices = [int(r["index_zero_based"]) for r in snapshots]
    require(indices == [0, 8, 63, 291], "Benchmark stored snapshot indices changed")
    for row in snapshots:
        require(int(row["solnum"]) == int(row["index_zero_based"]) + 1,
                "Benchmark solnum/index mismatch")
    times = [number(r, "time_s") for r in snapshots]
    provenance = json.loads((ROOT / "data/benchmark/paper_flux_plot_provenance.json")
                            .read_text(encoding="utf-8-sig"))
    steady_time = float(provenance["steady_time_s"])
    near(times[3], steady_time, "Snapshot/provenance steady time")
    near(steady_time, 0.009202752029406014, "Published steady time")
    require(provenance["spatial_nodes"] == 460, "Provenance spatial node count")
    require(len(provenance["snapshots"]) == 4, "Provenance snapshot count")
    for index, item in enumerate(provenance["snapshots"]):
        require(item["stage"] == STAGES[index] and item["stored_index"] == indices[index],
                "Snapshot/provenance stage or index mismatch")
        near(float(item["time_s"]), times[index], "Snapshot/provenance time")

    series = group_by(figures[5], "series")
    expected = {str(i): 460 for i in range(4)}
    expected.update({"4": 800, **{str(i): 70840 for i in range(5, 9)}})
    require({key[0]: len(value) for key, value in series.items()} == expected,
            "Figure 5: expected 4 x 460 profiles, 800 history states, 4 x 70840 domain nodes")
    profiles = []
    for index in range(4):
        profile = sorted(series[(str(index),)], key=lambda r: number(r, "x_m"))
        require(len({r["node"] for r in profile}) == 460, "Figure 5: duplicate profile nodes")
        for row in profile:
            require(row["stage"] == STAGES[index] and
                    int(row["index_zero_based"]) == indices[index],
                    "Figure 5: profile stage/index mismatch")
            near(number(row, "time_s"), times[index], "Figure 5 profile time")
            near(number(row, "tss_s"), steady_time, "Figure 5 declared steady time")
            require(number(row, "c_over_c0") > 0, "Figure 5: non-positive concentration")
        profiles.append(profile)
        cloud = series[(str(index + 5),)]
        require(len({r["node"] for r in cloud}) == 70840, "Figure 5: duplicate cloud nodes")
        metadata = [r for r in cloud if r["time_s"]]
        require(bool(metadata), "Figure 5: missing cloud time metadata")
        for row in metadata:
            near(number(row, "time_s"), times[index], "Figure 5 cloud time")
            near(number(row, "tss_s"), steady_time, "Figure 5 cloud steady time")
            require(int(row["index_zero_based"]) == indices[index], "Figure 5 cloud index")
        require(all(r["stage"] == "domain_" + STAGES[index] for r in cloud),
                "Figure 5 cloud stage")
        require(all(number(r, "c_over_c0") > 0 for r in cloud),
                "Figure 5: non-positive domain concentration")

    history = sorted(series[("4",)], key=lambda r: int(r["index_zero_based"]))
    require([int(r["index_zero_based"]) for r in history] == list(range(800)),
            "Figure 5: missing/duplicate history indices")
    history_times = [number(r, "time_s") for r in history]
    require(all(a < b for a, b in zip(history_times, history_times[1:])),
            "Figure 5: non-increasing history times")
    # Check the stored persistent trace against the backward maximum of the
    # stored profile residual, without claiming to recompute absent raw profiles.
    running_max = 0.0
    for row in reversed(history):
        running_max = max(running_max, number(row, "profile_relative_residual"))
        near(number(row, "persistent_profile_residual"), running_max,
             "Figure 5 persistent residual")
        near(number(row, "tss_s"), steady_time, "Figure 5 history steady time")
    accepted = [i for i, r in enumerate(history)
                if number(r, "persistent_profile_residual") <= 1e-4]
    require(bool(accepted) and accepted[0] == indices[3],
            "Figure 5: first persistent residual <= 1e-4 is not the declared steady state")
    near(history_times[accepted[0]], steady_time, "Figure 5 residual-based steady time")
    for index, profile in enumerate(profiles):
        record = history[indices[index]]
        near(number(record, "c_over_c0"), number(profile[0], "c_over_c0"),
             "Figure 5 snapshot/history left endpoint")
        near(number(record, "profile_relative_residual"),
             float(provenance["snapshots"][index]["residual"]),
             "Figure 5 residual/provenance snapshot")
    ratio = number(profiles[3][0], "c_over_c0") / number(profiles[3][-1], "c_over_c0")
    require(f"{ratio:.3f}" == "1.355", "Figure 5: steady endpoint ratio does not round to 1.355")

    groups = group_by(figures[6], "panel", "series")
    require(set(groups) == {(panel, str(i)) for panel in "abc" for i in range(4)},
            "Figure 6 panel/series set changed")
    for (panel, index_text), group in groups.items():
        index = int(index_text)
        require(len(group) == 460, "Figure 6: expected 460 nodes per field/time")
        for row in group:
            require(row["stage"] == STAGES[index] and
                    int(row["index_zero_based"]) == indices[index],
                    "Figure 6 snapshot stage/index mismatch")
            near(number(row, "time_s"), times[index], "Figure 6/5 snapshot time")
    print(f"PASS: Figures 5-6 snapshots, domain nodes and stored steady criterion; "
          f"t_ss={steady_time:.12g} s, endpoint ratio={ratio:.12g}")
    return times, steady_time


def check_flux(rows, times, steady_time):
    snapshots = [r for r in rows if r["panel"] in "abcd"]
    groups = group_by(snapshots, "panel", "sample")
    require(len(groups) == 4 * 460, "Figure 7: expected 4 x 460 spatial flux records")
    require(Counter(r["panel"] for r in snapshots) == dict.fromkeys("abcd", 1840),
            "Figure 7 spatial panel counts")
    maximum_error = 0.0
    for (panel, sample), group in groups.items():
        components = {r["component"]: r for r in group}
        require(len(group) == 4 and set(components) == {"Jdiff", "Jstrain", "Jelec", "Jx"},
                "Figure 7: missing/duplicate signed flux component")
        coordinate = number(group[0], "x_SI")
        for row in group:
            near(number(row, "time_s"), times["abcd".index(panel)], "Figure 7/5 snapshot time")
            near(number(row, "x_SI"), coordinate, "Figure 7 component coordinate", absolute=1e-20)
            near(number(row, "tss_s"), steady_time, "Figure 7 steady time")
        parts = [number(components[c], "value_SI") for c in ("Jdiff", "Jstrain", "Jelec")]
        total = number(components["Jx"], "value_SI")
        scale = sum(abs(value) for value in parts)
        error = abs(total - math.fsum(parts))
        require(error <= 1e-10 * scale + 1e-18,
                f"Figure 7: Jx != Jdiff + Jstrain + Jelec at {panel}/{sample}")
        maximum_error = max(maximum_error, error / max(scale, 1e-300))

    history_times = rows_at("data/benchmark/paper_flux_history_times.csv")
    require(len(history_times) == 292, "Benchmark flux-history time count")
    expected_times = {int(r["index_zero_based"]): number(r, "time_s") for r in history_times}
    shares = group_by([r for r in rows if r["panel"] == "e"], "sample")
    require({int(k[0]) for k in shares} == set(range(292)), "Figure 7 share-history indices")
    require(set(expected_times) == set(range(292)), "Benchmark flux-history indices")
    for (sample,), group in shares.items():
        require(len(group) == 3 and {r["component"] for r in group} ==
                {"Jdiff", "Jstrain", "Jelec"}, "Figure 7: expected three flux shares")
        values = [number(r, "value_SI") for r in group]
        require(all(0 <= value <= 1 for value in values), "Figure 7: share outside [0, 1]")
        near(math.fsum(values), 1.0, "Figure 7 flux-share sum", rel=0, absolute=1e-12)
        for row in group:
            near(number(row, "x_SI"), expected_times[int(sample)], "Figure 7 share time")
            near(number(row, "tss_s"), steady_time, "Figure 7 share steady time")
    print(f"PASS: Figure 7 signed flux closure (max relative residual {maximum_error:.3g}); "
          "292 three-component share sums")


def check_profiles(rows):
    groups = group_by(rows, "panel", "case_id")
    require(len(groups) == 12 and all(len(g) == 532 for g in groups.values()),
            "Figure 8: expected 12 profiles with 532 points each")
    require(Counter(key[0] for key in groups) ==
            {"strain": 4, "temperature": 4, "concentration": 4}, "Figure 8 parameter groups")
    for key, group in groups.items():
        ordered = sorted(group, key=lambda r: number(r, "x_nm"))
        x = [number(r, "x_nm") for r in ordered]
        require(all(a < b for a, b in zip(x, x[1:])), f"Figure 8 duplicate grid point: {key}")
        near(x[0], 0.0, "Figure 8 left boundary")
        near(x[-1], 300.0, "Figure 8 right boundary")
        for field, ratio_column in (("comsol_c_over_c0", "comsol_endpoint_ratio"),
                                    ("theory_c_over_c0", "theory_endpoint_ratio")):
            require(all(number(r, field) > 0 for r in ordered), "Figure 8 non-positive concentration")
            ratio = number(ordered[0], field) / number(ordered[-1], field)
            for row in ordered:
                near(number(row, ratio_column), ratio, f"Figure 8 endpoint ratio: {key}")
    print("PASS: Figure 8 twelve 532-point profiles and endpoint ratios")


def check_parameter_matrix(figures):
    records = rows_at("data/parameter_matrix/matrix_summary.csv")
    matrix = {r["case"]: r for r in records}
    require(len(records) == len(matrix) == 64, "Parameter matrix: expected 64 unique cases")
    grids = group_by(figures[9], "panel")
    require(set(grids) == {(p,) for p in "abcdef"}, "Figure 9 heatmap panels")
    for group in grids.values():
        cells = {(int(r["cell_row"]), int(r["cell_column"])) for r in group}
        require(len(group) == 16 and cells == {(i, j) for i in range(4) for j in range(4)},
                "Figure 9: expected 16 unique cells per heatmap")
    points = [r for r in figures[10] if r["type"] == "COMSOL"]
    for panel in "ab":
        subset = [r for r in points if r["panel"] == panel]
        require(len(subset) == 64 and {r["case"] for r in subset} == set(matrix),
                "Figure 10: each panel must contain the same 64 matrix cases")
    require(len(points) == 128, "Figure 10: expected two representations of 64 COMSOL cases")
    for row in [*figures[9], *points]:
        require(row["case"] in matrix, "Figure 9/10 case absent from parameter matrix")
        source = matrix[row["case"]]
        for column in ("endpoint_ratio", "temperature_K", "abs_strain", "c0_ratio", "Pi", "delta", "alpha"):
            near(number(row, column), number(source, column), f"Figure 9/10 matrix {column}")
        near(number(row, "alpha"), number(row, "Pi") * number(row, "delta"),
             "Figure 9/10 alpha=Pi*delta")
    curve = [r for r in figures[10] if r["type"] == "short_screening"]
    require(len(curve) == 401 and all(r["panel"] == "b" for r in curve),
            "Figure 10: expected 401 short-screening curve points in panel b")
    for row in figures[10]:
        near(number(row, "short_screening_ratio"), math.exp(2 * number(row, "alpha")),
             "Figure 10 short-screening exp(2*alpha)")
    for row in points:
        prediction = math.exp(2 * number(row, "alpha") *
                              math.tanh(1 / (2 * number(row, "delta"))))
        near(number(row, "linear_pb_ratio"), prediction, "Figure 10 finite-domain linear PB ratio")
    print("PASS: Figure 9 96 cells; Figure 10 64 cases per panel match the matrix; theory formulas")


def check_convergence(figure4):
    rows = rows_at("data/convergence/convergence_metrics.csv")
    metrics = {r["comparison"]: r for r in rows}
    mesh = metrics["mesh_base_vs_fine"]
    time = metrics["dtmax_vs_quarter_dtmax"]
    mesh_percent = number(mesh, "profile_max_difference_tmax_percent")
    time_percent = number(time, "profile_max_difference_all_times_percent")
    require(f"{mesh_percent:.4f}" == "0.0054", "Mesh final-profile error must round to 0.0054%")
    require(f"{time_percent:.4f}" == "0.0063", "Time-step full-history error must round to 0.0063%")
    require(mesh["required_acceptance"].lower() == time["required_acceptance"].lower() == "true",
            "Convergence required-acceptance flags changed")
    require(mesh["passes_1percent"].lower() == time["passes_1percent"].lower() == "true",
            "Convergence acceptance did not pass")
    curves = group_by(figure4, "panel", "curve_id")
    for panel, baseline, reference, record in (
        ("mesh", "h_base_nm", "h_fine_nm", mesh),
        ("timestep", "dtmax_base", "quarter_dtmax", time),
    ):
        first = sorted(curves[(panel, baseline)], key=lambda r: number(r, "x_nm"))
        second = sorted(curves[(panel, reference)], key=lambda r: number(r, "x_nm"))
        differences = []
        for a, b in zip(first, second):
            near(number(a, "x_nm"), number(b, "x_nm"), "Convergence comparison grid")
            differences.append(abs(number(a, "c_over_c0") / number(b, "c_over_c0") - 1) * 100)
        near(max(differences), number(record, "profile_max_difference_tmax_percent"),
             f"Figure 4/{panel} final-profile convergence metric", rel=1e-7, absolute=2e-11)
    # The 0.0063% value is a recorded all-times maximum, not the near-zero
    # final-time difference of the displayed time-step curves.
    require(number(time, "profile_max_difference_tmax_percent") < time_percent,
            "Time-step final/all-times convergence metrics are inconsistent")
    print(f"PASS: mesh final-profile error {mesh_percent:.10g}% -> 0.0054%; "
          f"recorded time-step full-history error {time_percent:.10g}% -> 0.0063%")


def check_checksums():
    manifest = ROOT / "SHA256SUMS"
    if not manifest.exists():
        print("SKIP: SHA256SUMS is absent (optional package-integrity check)")
        return
    checked = 0
    skipped = 0
    seen = set()
    for line_number, line in enumerate(manifest.read_text(encoding="utf-8-sig").splitlines(), 1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        match = re.fullmatch(r"([0-9a-fA-F]{64})\s+\*?(.+)", line)
        require(match is not None, f"Malformed SHA256SUMS line {line_number}")
        expected, relative = match.groups()
        relative = relative.replace("\\", "/")
        require(not re.match(r"^[A-Za-z]:", relative), "Absolute checksum path is not portable")
        path = (ROOT / relative).resolve()
        require(path.is_relative_to(ROOT) and path != ROOT, "Checksum path escapes release root")
        parts = path.relative_to(ROOT).parts
        if any(part in GENERATED_DIRS or part.startswith("generated_") for part in parts):
            skipped += 1
            continue
        require(path not in seen, f"Duplicate checksum entry: {relative}")
        seen.add(path)
        require(path.is_file(), f"Checksum file missing: {relative}")
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
        require(digest.hexdigest() == expected.lower(), f"Checksum mismatch: {relative}")
        checked += 1
    require(checked > 0, "SHA256SUMS contains no verifiable released files")
    print(f"PASS: SHA256SUMS {checked} files verified; {skipped} generated entries skipped")


def main():
    figures = {i: rows_at(f"data/figures/figure{i}_data.csv") for i in EXPECTED_ROWS}
    check_counts(figures)
    times, steady_time = check_benchmark(figures)
    check_flux(figures[7], times, steady_time)
    check_profiles(figures[8])
    check_parameter_matrix(figures)
    check_convergence(figures[4])
    check_checksums()
    print("All released-data checks passed. No files were modified.")


if __name__ == "__main__":
    try:
        main()
    except (ValidationError, OSError, KeyError, ValueError, TypeError, OverflowError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        sys.exit(1)
