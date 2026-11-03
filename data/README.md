# Data dictionary

These tables contain the processed data used by the manuscript figures and selected supporting summaries. CSV files have a header row and use UTF-8 encoding; some have a UTF-8 byte-order mark. Row counts below exclude the header. Empty fields normally mean that a quantity does not apply to that row type. The exception is the shared metadata in Figure 5's two-dimensional series 5-8: those values are stored only in each series' first row and inherited by the remaining nodes, as specified below. The literal value `nan` means that a diagnostic was unavailable or its criterion was not reached; it is not zero.

## Common conventions

- `C = c/c0` is normalized vacancy concentration. `c0` is the initial concentration for the particular case; `c_ref = 2.002 mol m^-3` is the baseline reference concentration.
- `temperature_K` is in kelvin. `abs_strain` and `strain_abs` are positive magnitudes of the applied compressive strain as fractions, so `0.08` means 8% compression. The signed applied strain is negative.
- `c0_ratio` and `c0_over_c_ref` are `c0/c_ref`; `c0_mol_m3` is in mol m^-3.
- `x_m` and `y_m` are reference coordinates in metres, `x_nm` is in nanometres, and `x_over_L` is the dimensionless longitudinal coordinate. Positive longitudinal flux points toward increasing x. The compressed end is x = 0.
- `Pi` is the formation-energy change across the sample divided by thermal energy; `delta = lambda_D/L`; `alpha = Pi*delta`. These quantities are dimensionless.
- `endpoint_ratio = C(0)/C(L)`. All logarithms are natural logarithms.
- `source_file`, `comsol_source_file`, and file-valued `theory_parameter_source` entries are portable provenance labels. Multiple source labels and their hashes are separated by `|` in corresponding order. `source_data/` labels do not imply that those upstream files are included. See `../PROVENANCE.csv`.
- `source_sha256` hashes identify upstream source bytes, not necessarily the edited release copies. Distributed-file checksums are in `../SHA256SUMS`.
- In `../PROVENANCE.csv`, `availability = included` identifies a supplied counterpart, `not_included` identifies an upstream file outside this package, and `unavailable_at_preparation` identifies a source that could not be located when the package was assembled. `original_sha256` records the upstream checksum. `hash_basis = computed_from_available_source` means that this checksum was computed from the available upstream file during preparation. `hash_basis = recorded_in_snapshot_not_reverified` means that it was retained from an existing snapshot because the upstream file was unavailable; it is a historical record, not a newly verified source checksum.
- `processing_notes` and `provenance` describe extraction, fitting, or numerical conventions. They are retained with their numerical records.

## Figure data

| File in `figures/` | Rows | Record structure |
|---|---:|---|
| `figure3_data.csv` | 19 | 6 formation-energy points, 5 migration-barrier points, 8 diffusion points |
| `figure4_data.csv` | 2,760 | 6 profiles, each with 460 positions |
| `figure5_data.csv` | 286,000 | 4 line profiles, 800 history states, 4 two-dimensional fields |
| `figure6_data.csv` | 5,520 | 3 quantities x 4 times x 460 positions |
| `figure7_data.csv` | 8,236 | 4 fluxes x 4 times x 460 positions, plus 3 shares x 292 times |
| `figure8_data.csv` | 6,384 | 3 parameter scans x 4 cases x 532 positions |
| `figure9_data.csv` | 96 | 6 heatmaps x 16 cells; some cases appear in multiple slices |
| `figure10_data.csv` | 529 | 64 cases in each of 2 panels, plus 401 analytical-curve points |

### Figure 3: atomistic calibration

The `panel` field selects `formation`, `migration`, or `arrhenius`. `point_index` is a one-based index within that group.

| Column | Meaning and units |
|---|---|
| `position_angstrom` | Longitudinal vacancy position in the atomistic strain-gradient specimen, angstrom |
| `x_over_L` | Position divided by atomistic specimen length 94.374 angstrom |
| `strain_fraction` | Signed local strain; negative values denote compression |
| `temperature_K` | MD temperature, K |
| `energy_eV` | Formation energy or migration barrier, selected by `panel`, eV |
| `D_V_m2_s` | Vacancy diffusivity from the three-dimensional vacancy MSD slope divided by 6, m^2 s^-1 |
| `inverse_temperature_1000_per_K` | 1000 divided by temperature in K, as used on the Arrhenius horizontal axis |
| `ln_D_V_SI` | Logarithm of the numerical diffusivity expressed in m^2 s^-1 |
| `msd_timestep_ps` | MD time per recorded step, ps |
| `msd_fit_start_ns`, `msd_fit_end_ns` | MSD fitting interval, ns |
| `msd_fit_npoints` | Number of points in the MSD fit |
| `msd_fit_slope_A2_per_ns` | Fitted MSD slope, angstrom^2 ns^-1 |
| `msd_fit_r2` | Coefficient of determination of that MSD fit |
| `source_file`, `source_sha256`, `processing_notes` | Source identification and fitting/extraction conventions |

The formation points come from a strain-gradient specimen, not a newly generated uniform-strain series. The migration points are from uniform-strain NEB calculations. The `D_V_m2_s` values are vacancy diffusivities evaluated from the vacancy MSD. Raw atomistic trajectories and input decks are outside this package.

### Figure 4: convergence profiles

| Column | Meaning |
|---|---|
| `panel` | `mesh` or `timestep` |
| `curve_id` | Source-table column identifying the mesh or maximum-step case |
| `legend_label` | Figure legend text, including the displayed physical spacing or time step |
| `x_nm` | Position in nm on the common 460-point grid |
| `c_over_c0` | Normalized vacancy concentration |
| `source_file`, `source_sha256` | Supporting profile table and its original checksum |

The grid profiles have already been evaluated on the common comparison coordinates. The columns named `h_coarse_nm`, `h_base_nm`, and `h_fine_nm` in the supporting wide-format table contain concentrations for those cases, not numerical mesh sizes. The mesh spacings are approximately 1.304, 0.654, and 0.327 nm. The maximum internal steps are approximately 5.206e-5, 2.603e-5, and 1.302e-5 s.

### Figure 5: baseline concentration

| `series` | `stage` | Rows | Interpretation |
|---|---|---:|---|
| 0-3 | `initial`, `early`, `late`, `steady` | 460 each | Native-node longitudinal profiles at y = 0 |
| 4 | `left_end_history` | 800 | Concentration at x = 0 and profile residuals over saved times |
| 5-8 | `domain_initial`, `domain_early`, `domain_late`, `domain_steady` | 70,840 each | Native two-dimensional concentration fields |

**Metadata inheritance for series 5-8:** each group contains 70,840 node rows. Only its first row (`node = 0`) stores `index_zero_based`, `time_s`, `temperature_K`, `abs_strain`, `c0_ratio`, `c0_mol_m3`, `tss_s`, `source_file`, `source_sha256`, and `provenance`. Blank entries in these ten fields on subsequent rows inherit that first-row value within the same series. The four stored-solution indices are 0, 8, 63, and 291 for series 5, 6, 7, and 8, respectively. Do not carry these values across series boundaries. Every node row explicitly contains `series`, `stage`, `node`, `x_m`, `y_m`, and `c_over_c0`. The two residual columns remain blank throughout these four groups because residual histories apply only to series 4; they must not be filled from another series.

| Column | Meaning and units |
|---|---|
| `index_zero_based` | Zero-based stored-solution index |
| `node` | Zero-based spatial-node index within the series |
| `x_m`, `y_m` | Reference coordinates in m; blank y for line-profile records is not a missing concentration |
| `c_over_c0` | Normalized concentration at the node or at the left end |
| `time_s` | Actual saved time, s |
| `temperature_K`, `abs_strain`, `c0_ratio`, `c0_mol_m3` | Case metadata, as defined above |
| `tss_s` | Baseline settling time, 0.009202752029406014 s |
| `profile_relative_residual` | `norm(C(t)-C_final)/norm(C_initial-C_final)` over the 460 line nodes; populated for series 4 |
| `persistent_profile_residual` | Maximum of that residual from the current saved time through the final saved time |
| `source_file`, `source_sha256`, `provenance` | Upstream identifiers, original hashes, and extraction notes |

The first history row with `persistent_profile_residual <= 1e-4` is index 291, giving the manuscript value 9.203 ms. The displayed early and late snapshots are actual stored times, approximately 0.253 and 1.992 ms. The plotting script spatially interpolates the native two-dimensional nodes only to render inset maps. It does not replace the archived node values. The upstream standalone two-dimensional CSV is unavailable; all four native-node fields used for the figure are included here. Its recorded `source_sha256` is retained in this table and in `../PROVENANCE.csv` with `hash_basis = recorded_in_snapshot_not_reverified`, rather than being presented as a checksum reverified against the missing file.

### Figure 6: electrostatic fields

Each `panel` contains four `series` (0-3) and `stage` values (`initial`, `early`, `late`, `steady`). `index_zero_based` is the stored-solution index and `node` is the spatial index.

| `panel` | `quantity` | `value_SI` |
|---|---|---|
| a | `rho_C_m3` | Space-charge density, C m^-3 |
| b | `phi_V` | Electrostatic potential, V |
| c | `E_x_V_m` | Longitudinal electric field, V m^-1 |

`x_SI` is position in m. `time_s` and the case/provenance columns follow the common conventions. The potential uses the model's gauge; physical interpretation should use potential differences.

### Figure 7: flux components

| Columns or values | Meaning |
|---|---|
| `panel` a-d | Initial, early, late, and steady spatial profiles, respectively |
| `component` | `Jdiff`, `Jstrain`, `Jelec`, or total `Jx` |
| `quantity` ending in `_mol_m2_s` | `value_SI` is flux in mol m^-2 s^-1; `x_SI` is position in m |
| `panel` e | Time-dependent shares of the three component fluxes; total `Jx` has no share series |
| `quantity = spatial_mean_absolute_flux_share` | `value_SI` is the dimensionless share `integral(abs(J_i),dx)/sum_j integral(abs(J_j),dx)`; `x_SI` is time in s |
| `sample` | Zero-based sample index within a panel/component |
| `time_s` | Snapshot time for panels a-d; blank for panel e because its times are stored in `x_SI` |
| `tss_s` | Same baseline settling time as Figure 5 |

Case and provenance columns follow the common definitions. The shares use absolute flux magnitudes; they are not signed fractions of the net flux. The total flux is the signed sum `Jdiff + Jstrain + Jelec`.

### Figure 8: scans and nonlinear equilibrium profiles

`panel` is `strain`, `temperature`, or `concentration`; `case_id` identifies one of four cases in that scan. Each case has 532 positions. The repeated baseline across scans is intentional.

| Column | Meaning and units |
|---|---|
| `temperature_K`, `strain_abs`, `c0_over_c_ref`, `Pi`, `delta` | Case parameters and dimensionless groups |
| `x_m`, `x_nm`, `x_over_L` | The same position expressed in m, nm, and as a fraction of L |
| `comsol_c_over_c0` | Normalized concentration from COMSOL |
| `theory_c_over_c0` | Normalized concentration from the independently calculated nonlinear Poisson-Boltzmann equilibrium |
| `pb_iterations` | Number of nonlinear solver iterations used upstream |
| `theory_mass_relative_error` | Signed `integral(C_theory,d(x/L)) - 1` |
| `L2_relative_error` | Discrete Euclidean norm of theory-minus-COMSOL divided by the COMSOL norm |
| `max_abs_error` | Maximum absolute difference in normalized concentration |
| `comsol_endpoint_ratio`, `theory_endpoint_ratio` | Endpoint ratios of the respective profiles |
| `comsol_source_file` | Upstream profile identifier |
| `theory_parameter_source` | Parameter-table identifier or literal description of the parameters used |

Case-level diagnostics repeat on every position in that case; do not count them as independent samples. The supplied script reads the stored theoretical profiles; it does not rerun the Poisson-Boltzmann solver.

### Figure 9: heatmap slices

| Column | Meaning |
|---|---|
| `panel` | Individual heatmap identifier a-f |
| `panel_group`, `slice_position` | Manuscript group a-c and left/right member of that pair |
| `case` | Parameter-matrix case identifier |
| `temperature_K`, `abs_strain`, `c0_ratio`, `Pi`, `delta`, `alpha`, `endpoint_ratio` | Physical case and response |
| `fixed_parameter`, `fixed_value` | Variable held fixed for that slice, and its value in the named variable's units |
| `x_parameter`, `y_parameter` | Variables placed on the heatmap axes |
| `cell_column`, `cell_row` | Zero-based column and row in the 4 x 4 grid |
| `color_min`, `color_max`, `colormap` | Shared color normalization and Matplotlib colormap |
| `source_file`, `source_sha256` | Source matrix identifier and original checksum |

All cells are actual cases. Unequally spaced concentration ratios are shown as equally spaced categorical cells; there is no interpolation between cases.

### Figure 10: dimensionless response

`type = COMSOL` supplies 64 cases in each panel. `type = short_screening` supplies 401 curve points in panel b; case-specific columns are blank for those points.

| Column | Meaning |
|---|---|
| `panel`, `type`, `case` | Plot membership, row type, and case identifier |
| `temperature_K`, `abs_strain`, `c0_ratio`, `Pi`, `delta`, `alpha`, `endpoint_ratio` | Physical case and numerical response |
| `linear_pb_ratio` | Finite-length linearized prediction `exp(2*alpha*tanh(1/(2*delta)))` |
| `linear_ratio_error_percent` | Signed relative error, `100*(linear_pb_ratio/endpoint_ratio - 1)` |
| `short_screening_ratio` | Short-screening prediction `exp(2*alpha)` |
| `short_screening_error_percent` | Absolute relative error, `100*abs(short_screening_ratio/endpoint_ratio - 1)`; blank on curve-only rows |
| `color_min`, `color_max` | Shared numerical response color range |
| `source_file`, `source_sha256` | Source matrix identifier and original checksum |

## Baseline supporting files

`benchmark/paper_flux_model_parameters.csv` has eight rows, with columns `parameter` and `value_SI`:

| Parameter | Meaning | Unit |
|---|---|---|
| T | Temperature | K |
| eps0 | Signed applied strain | 1 |
| c0 | Initial vacancy concentration | mol m^-3 |
| L, H | Domain length and height | m |
| D0 | Diffusion prefactor used by the continuum model | m^2 s^-1 |
| epsr | Relative permittivity | 1 |
| zV | Vacancy charge number | 1 |

`benchmark/paper_flux_snapshot_times.csv` has four rows. `stage` identifies the displayed time; `index_zero_based` is the zero-based stored state; `solnum` is COMSOL's corresponding one-based solution number; `time_s` is the actual time in seconds.

`benchmark/paper_flux_history_times.csv` has 292 rows, with `index_zero_based` and `time_s`. It lists the saved states through the selected settling state used for the flux-share history. The complete 800-state left-end concentration history is in Figure 5's table.

`benchmark/paper_flux_plot_provenance.json` records the baseline parameters, y = 0 cut line, 460 spatial nodes, source hashes, the persistent Euclidean criterion, four snapshot records, and historical plotting metadata. Each snapshot records its saved index/time, profile residual, mean absolute net flux in mol m^-2 s^-1, relative component-sum error, and net flux relative to its initial magnitude. The `old_absolute_max_profile_criterion_time_s` field is retained as historical metadata; it is not the settling time reported by the current manuscript. Historical `output` labels and plotting-style text are provenance, not required filenames or plot configuration for the supplied scripts.

## Convergence supporting files

`convergence/mesh_profiles_tmax.csv` and `convergence/timestep_profiles_tmax.csv` each have 460 rows. `x_nm` is the common position coordinate. Every other column contains normalized concentration: `h_coarse_nm`, `h_base_nm`, and `h_fine_nm` select mesh cases; `dtmax_base`, `half_dtmax`, and `quarter_dtmax` select time-step cases. Profiles are evaluated at the final time of the convergence runs, 0.1 s.

`convergence/case_summary.csv` has seven rows:

| Column | Meaning and units |
|---|---|
| `case_key`, `case_name` | Short and full identifiers; coarse/base/fine mesh, half/quarter maximum step, or loose/tight solver tolerance |
| `node_count`, `saved_time_count` | Number of line nodes and saved time states |
| `strict_steady_time_s` | First time after which the maximum pointwise relative difference from the final profile remains <= 1e-4, s |
| `final_mass_relative_error` | Signed relative departure of total vacancy content from its initial value |
| `final_charge_relative_error` | Net charge divided by the initial vacancy-charge scale |
| `final_min_cVI_mol_m3` | Minimum final vacancy concentration, mol m^-3 |
| `final_Cratio` | Final endpoint concentration ratio |
| `final_DeltaPhi_V` | Final potential difference phi(L)-phi(0), V |
| `final_etaE` | Screening fraction, `-DeltaPhiND/Pi` |
| `final_PiEff` | Effective driving parameter, `Pi + DeltaPhiND` |
| `profile_residual_at_comparison_time` | Maximum pointwise relative profile residual at 0.007 s |

These are 721-state convergence runs. Their baseline `strict_steady_time_s = 0.00675` uses a different residual and saved-time protocol from the 800-state Euclidean settling time in Figures 5-7.

`convergence/convergence_metrics.csv` has six comparison rows:

| Column | Meaning |
|---|---|
| `comparison` | Tested case followed by its reference case |
| `required_acceptance` | Whether the source convergence analysis designated this comparison as required |
| `profile_max_difference_tmax_percent` | Maximum pointwise relative concentration difference at final time, percent |
| `profile_max_difference_at_comparison_time_percent` | Same at the common comparison time 0.007 s |
| `profile_max_difference_all_times_percent` | Largest value over all common saved times, percent |
| `time_of_max_difference_s` | Saved time at which that largest profile difference occurs, s |
| `Cratio_difference_tmax_percent`, `DeltaPhi_difference_tmax_percent`, `etaE_difference_tmax_percent`, `PiEff_difference_tmax_percent` | Absolute relative differences in the indicated final diagnostic, percent |
| `strict_steady_time_difference_percent` | Absolute relative difference in the convergence-study steady time, percent |
| `passes_1percent` | Source acceptance flag requiring final/comparison-time profile errors, final diagnostic differences, and steady-time difference below 1%; it is not a flag for every available diagnostic |

The mesh profiles are compared on the tested case's coordinates with interpolation of the reference where needed. Absolute relative errors use the reference quantity in the denominator. The manuscript's 0.0054% value is row `mesh_base_vs_fine`, column `profile_max_difference_tmax_percent` (0.005410631980461907). Its 0.0063% value is row `dtmax_vs_quarter_dtmax`, column `profile_max_difference_all_times_percent` (0.006258973277387292).

## Full parameter matrix

`parameter_matrix/matrix_summary.csv` has 64 rows, one per combination of temperatures 450/500/550/600 K, compressive-strain magnitudes 0.05/0.06/0.07/0.08, and initial-concentration ratios 0.04/0.16/0.40/1.00. This table retains additional diagnostics beyond those plotted in Figures 9-10.

### Physical parameters and steady response

| Columns | Meaning and units |
|---|---|
| `case` | Source run identifier; its embedded mesh label describes the parameter-scan run, not the separate baseline benchmark |
| `temperature_K`, `abs_strain`, `c0_ratio`, `c0_mol_m3` | Temperature, compressive magnitude, normalized and absolute initial concentration |
| `Pi`, `delta`, `alpha` | Dimensionless driving, screening, and their product |
| `lambdaD_nm` | Reference Debye length, nm |
| `beta_diffusivity_gradient` | `ln(Dmax/Dmin)`, dimensionless |
| `Dref_geometric_m2_s` | Geometric-mean reference diffusivity `sqrt(Dmin*Dmax)`, m^2 s^-1 |
| `Dharmonic_m2_s` | Harmonic spatial mean diffusivity, m^2 s^-1 |
| `Dmin_m2_s`, `Dmax_m2_s` | Minimum and maximum local diffusivity, m^2 s^-1 |
| `tauD_s` | `lambda_D^2/Dref_geometric`, s |
| `endpoint_ratio`, `ln_endpoint_ratio` | Final COMSOL endpoint ratio and its logarithm |
| `global_max_C`, `global_min_C` | Maximum and minimum normalized concentration in the final domain |
| `DeltaPhi_V` | Final phi(L)-phi(0), V |
| `DeltaPhiND` | Dimensionless potential difference `zV*F*DeltaPhi/(R*T)` |
| `eta_screening` | `-DeltaPhiND/Pi` |
| `PiEff` | `Pi + DeltaPhiND` |
| `electrochemical_log_closure_abs` | `abs(ln(endpoint_ratio)-PiEff)` |

### Equilibrium comparisons

Unless stated otherwise, an error in this section is signed and equals `100*(prediction/reference - 1)`.

| Columns | Meaning |
|---|---|
| `linear_pb_ratio`, `linear_pb_ln_ratio` | `exp(2*alpha*tanh(1/(2*delta)))` and its logarithm |
| `linear_ratio_error_percent`, `linear_ln_error_percent` | Errors against the COMSOL ratio or its logarithm, respectively |
| `pb_endpoint_ratio`, `pb_ln_endpoint_ratio` | Refined finite-domain nonlinear Poisson-Boltzmann result (2,125 nodes) and its logarithm |
| `pb_ratio_error_percent` | Nonlinear PB ratio error against COMSOL |
| `nonlinear_thin_layer_ratio` | Non-overlapping nonlinear PB layer prediction: `exp(y_plus-y_minus)`, where the positive and negative roots satisfy `exp(y)-y-1=alpha^2/2` |
| `thin_layer_ratio_error_percent` | Nonlinear thin-layer ratio error against COMSOL |
| `thin_layer_vs_finite_pb_percent` | Nonlinear thin-layer ratio error against the refined finite-domain PB ratio |
| `cubic_ln_ratio` | Cubic approximation `2*alpha + alpha^3/18` to the logarithmic ratio |
| `cubic_ratio_error_percent` | Error of `exp(cubic_ln_ratio)` against COMSOL |
| `pb_grid_refinement_ratio_percent` | Difference of the 1,063-node PB ratio relative to the 2,125-node ratio |
| `pb_mass_abs_error` | `abs(integral(C_pb,d(x/L))-1)` on the refined grid |
| `pb_iterations` | Iterations of the 1,063-node PB solve |
| `neutral_ratio` | Unscreened prediction `exp(Pi)` |
| `neutral_overcharged_ratio` | `neutral_ratio/endpoint_ratio` |
| `pb_Qplus_C_m2` | PB integrated positive space charge per transverse area, C m^-2 |
| `pb_separated_percent` | `100*pb_Qplus/(F*c0*L)` |
| `pb_voltage_abs_V` | Magnitude of PB endpoint potential difference, V |
| `pb_Qplus_error_percent` | PB positive-charge error against the spatial COMSOL estimate, when available |

The nonlinear thin-layer and cubic predictions are retained supporting diagnostics. They are not the same as Figure 10's short-screening linearized prediction `exp(2*alpha)`.

### Time and numerical diagnostics

| Columns | Meaning and units |
|---|---|
| `t95_ratio_s`, `t95_voltage_s` | First 95% completion of the initial-to-final change in endpoint ratio or absolute voltage, with linear interpolation between saved times, s |
| `t95_ratio_over_tauD`, `t95_voltage_over_tauD` | Corresponding times divided by `tauD_s` |
| `voltage_to_ratio_t95` | Voltage 95% time divided by ratio 95% time |
| `tss_flux_s` | First saved time after which `max(abs(Jx))/max(abs(Jstrain)) <= 1e-4` persists, s; `nan` if not reached |
| `final_flux_residual_relative` | Final value of that flux residual |
| `tmax_s` | Last saved simulation time, s |
| `max_massRel` | Maximum absolute domain-integrated relative vacancy-content error over time |
| `endpoint_final20_relative_drift`, `voltage_final20_relative_drift` | Largest relative deviation from the final endpoint ratio or voltage over the final 20 saved states |
| `max_chargeRel` | Maximum absolute net charge normalized by the initial vacancy-charge scale |
| `max_strain_error` | Maximum absolute deviation of solved strain from prescribed strain |
| `energy_drop_J_per_m` | Initial minus final model free energy per out-of-plane length, J m^-1 |
| `energy_rise_max_fraction` | Largest positive single-step energy increase divided by the absolute total energy drop |
| `energy_increases_above_1e-8_drop` | Count of single-step increases exceeding 1e-8 of the absolute total drop |

### Spatial diagnostics with incomplete upstream coverage

`spatial_export_complete` is 1 where the source analysis had the required spatial exports and 0 otherwise. In the distributed `matrix_summary.csv`, all 64 rows have `spatial_export_complete = 0`, and the COMSOL spatial diagnostics listed below are `nan`. This describes the supporting matrix table; the separate figure-specific concentration and field tables do contain their documented profiles. Availability of an entry does not mean its full upstream spatial files are included in this package.

| Columns | Meaning and units |
|---|---|
| `strict_profile_tss_s` | First saved time after which the normalized Euclidean profile residual remains <= 1e-4, s |
| `Qplus_C_m2` | `F*c0*integral(max(C-1,0),dx)` at final time, C m^-2 |
| `separated_percent` | `100*Qplus/(F*c0*L)` |
| `t95_Q_s` | 95% completion time for positive-charge separation, s |
| `w90_left_over_lambdaD`, `w90_right_over_lambdaD` | Distance from each end enclosing 90% of the corresponding excess/depleted concentration, divided by Debye length |
| `pb_profile_L2_relative` | Spatially integrated L2 error of the interpolated refined PB concentration relative to COMSOL |
| `mu_span_over_RT` | Range of the relative electrochemical potential divided by RT across the final profile |
| `direct_flux_closure_relative` | Maximum absolute difference between stored net flux and sum of components, normalized by the maximum sum of absolute component magnitudes |
| `line_mass_max_error` | Largest absolute departure of `integral(C,dx)/L` from 1 over saved times |

Missing diagnostics must be excluded explicitly from any secondary statistical analysis. Do not replace `nan` values with zeros or interpret them as successful convergence checks.
