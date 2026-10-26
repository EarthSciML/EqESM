# EarthSciModels / ESD / EqIce / EarthSciAST survey — reuse opportunities for MOSART

Research deliverable (read-only) for the EqESM manager, to be folded into
PLAN.md §3/§5 and stage1/boundaries.md §2.6. Surveys:

1. EarthSciModels corpus (full) — components, lib, couplings, registered_functions
2. EarthSciDiscretizations (grids, rules, stencils, problems)
3. EqIce `stage2/` — the stage-2 factoring/test pattern
4. EarthSciAST — esm 1.3.0 spec features (`function_tables`, `data_sources`,
   `solver`, `reaction_systems`, `tests`)

Reference model context: MOSART mosart1.1.13 (`stage1/cesm/components/mosart/src/riverroute/mosart_physics.F90`):
hillslope overland flow via Manning-type law `CREHT_nosqrt`, tributary/subnetwork
routing and main-channel kinematic wave `Routing_KW` via `CRVRMAN_nosqrt`, the
`GR*` cross-section geometry functions, per-cell CFL substepping
(`numDT_r`/`numDT_t`, `DLevelH2R`), irrigation/flood/direct-transfer terms, and an
ID/dnID network topology (`frivinp` netCDF) that is *not* a lat-lon PDE grid.

---

## 1. Executive summary

**What can be reused:**

- **Nothing computes river routing, kinematic wave, or network routing anywhere
  in the surveyed corpora.** `environmental_transport/surface_runoff/` is a full
  1-D Saint-Venant **PDE** (Wang et al. 2020, mass + momentum) with a
  non-dimensionalized Manning friction slope — related physics (Manning
  friction, overland flow) but a different model class than MOSART's
  storage-based kinematic approximation. Verdict for the PLAN row: **confirm
  "related but NOT river routing; treat as new"** — with the caveat that the
  Saint-Venant files *do* provide the reusable structural patterns below.
- **Structural patterns to reuse (not equations):**
  - EqIce stage-2 factoring — pointwise laws as `expression_templates`
    (params/match/body) in a `*_templates.esm` library imported by reference
    from a physics component; continuous-ODE form with the Fortran per-substep
    Euler documented as a "scheme gap" in test tolerances; per-cell/pointwise
    assertions via `coords` at index-space positions. This is the closest
    existing analog to MOSART's structure, especially
    `eqice/stage2/hydrology.esm` + `hydrology_templates.esm` (Routing/Shreve
    flux law, staggered operators, upwind flux — all as templates).
  - The `faq`/`index` machinery (esm-spec §4.3) expresses the two SMM route
    handles (`rh_eroutUp`, `rh_direct`) and basin accumulation as gather
    expressions — `faq[i] = index(erout, index(dnID, i))` — with ragged index
    sets for basin membership. No new grid or operator is needed for this.
  - `data_sources` (esm-spec §8) with `kind: "static"` + per-variable
    `file_variable` parameter bindings — the `frivinp` routing-network netCDF
    maps directly onto this (era5/ncep_ncar/landfire precedents; netCDF is
    read by the runtime reader — era5 uses `.nc`).
  - The `couplings/` param_to_var pattern for `mosart_coupling.esm`.
- **Must be authored new:** the MOSART physics itself — `mosart_routing.esm`
  (hillslope + subnetwork + main-channel storage ODEs), a
  `mosart_physics_templates.esm`-style template library for the pointwise laws
  (`CRVRMAN_nosqrt`, `CREHT_nosqrt`, `GRMR/GRHR/GRPR/GRRR/GRHT/GRPT`),
  irrigation/flood/direct-transfer rules, and a `data_sources` entry for the
  network file. There is no kinematic-wave template, no network/graph routing
  grid in EarthSciDiscretizations, and no CFL-substepping template anywhere.

---

## 2. Candidate-by-candidate findings

All paths under `/projects/illinois/eng/cee/ctessum/ctessum/code/`.

### 2.1 EarthSciModels (esm corpus: 276 files at 1.0.0, 4 at 1.1.0)

| Path | esm ver | What it computes | Imports / exports | Verdict | Notes |
|---|---|---|---|---|---|
| `EarthSciModels/components/environmental_transport/surface_runoff/saint_venant.esm` | 1.0.0 | Full 1-D Saint-Venant PDE: ∂h̃/∂t = −∂q/∂l + (P−I); ∂q/∂t = −∂F/∂l + g·h̃·(S₀−Sᶠ); observed `aggregate` fields F_mom, S_f (Manning friction slope, non-dimensionalized via n_ref/q_ref/h_ref). Discretization-agnostic leaf (spatial `D` ops, no grid). | Imports: none (tests inject `saint_venant_upwind1_inflow.esm` per §9.7.10 form C). Exports: h_tilde(l), q_flux(l), F_mom, S_f; params g, n_mann, S_0_slope, P_rate, I_rate, h_min, h_bc, q_bc, dx, refs. | **Related-but-distinct** | The PLAN's "Saint-Venant overland flow" row: confirmed correct. Manning friction concept transfers; equations do NOT (MOSART hillslope is storage-kinematic: `ehout = −CREHT_nosqrt(...)`, no momentum equation; MOSART Manning is dimensional, `v = (rr²)^(1/3)·√s/n`). Reuse the *structure*: observed algebraic closure field, discretization-agnostic leaf + per-test injection. Note: uses `aggregate` (deprecated alias for `faq` since 1.1.0, esm-spec §4.2) — new code must use `faq`. |
| `.../surface_runoff/surface_runoff_component.esm` | 1.0.0 | Single-node ODE reduction of the same Saint-Venant system: dqdl/dFdl exposed as **input parameters** (spatial divergences supplied by a discretization or coupled model); Manning closure as observed S_f. | Params in: P, I_infil, S_0, n_mann, dqdl, dFdl, ...; states h_tilde, q. | **Related-but-distinct** | Structurally the closest analog to MOSART's per-cell storage ODEs (pointwise ODE + algebraic closure, divergence terms as inputs). The *shape* of this component (states + observed closure, params for lateral terms) is the right template for `mosart_routing.esm`'s cells; the physics differs (kinematic wave vs momentum equation). |
| `.../surface_runoff/saint_venant_upwind1_inflow.esm` | 1.0.0 | Discretization **adapter**: first-order upwind spatial rules for the Saint-Venant leaf, per-variable inflow BCs (q_bc, F_bc), interior via ESD `upwind1_D_interior` stencil imported with `${ESD_ROOT}` env-var ref and `rename {x: l}`. | Imports `esd://...upwind1_D_interior.esm` (via `${ESD_ROOT}`); exports ground rules D(q_flux,l), D(F_mom,l). | **No overlap for MOSART physics**; pattern reuse | MOSART is not a regular-grid PDE, so the adapter does not transfer. But the adapter pattern (per-variable ground rules, ESD import with rename, test-time injection) is the model for any future ESD PR (e.g. a kinematic-wave 1-D rule). |
| `.../surface_runoff/heaviside_bc.esm` | 1.0.0 | BC template for the Saint-Venant inflow. | — | **No overlap** | Saint-Venant-specific BC. |
| `EarthSciModels/components/urban_canopy/hydro/impervious_runoff.esm` (and `surface_runoff_infiltration.esm`, `pervious_road_water_balance.esm`, `snow_capping_runoff.esm`, ...) | 1.0.0 | CLMU (CTSM) runoff *generation*: ponding-limit excess runoff, infiltration, pervious/impervious balances (Oleson et al. 2010 CLM tech note). Observed-only algebraic components. | Inputs w_liq_1, q_liq_0, q_seva, has_snow → q_over, w_liq_1_new. | **No overlap with routing**; relevant to the CTSM milestone | These produce the `Flrl_rofsur`-family inputs MOSART consumes in coupled runs; they are NOT routing. Keep in mind for CTSM work, not MOSART. |
| `EarthSciModels/components/earthsci_data/era5.esm` | 1.0.0 | Reanalysis data loader: `data_sources` (kind grid, `.nc`, url_template with `{date:...}` substitutions, temporal block, time_variable) + `function_tables` (1-axis linear pressure-level table) + loader fields consumed as `ERA5.pl.<var>`. | Exports loader fields via `update: {kind: data, source, from: {file_variable}}` bindings; regridding at coupling edge. | **Reuse as-is as the data_source pattern** | This is the canonical `data_sources` + `function_tables` exemplar in the corpus. Directly applicable: the MOSART `frivinp` netCDF becomes a `data_sources` entry with per-variable `file_variable` bindings. |
| `EarthSciModels/components/earthsci_data/ncep_ncar.esm`, `landfire.esm`, `usgs3dep.esm`, `mtbs_severity_loader.esm` | 1.0.0 | Data loaders, several with `kind: "static"` (time-invariant) sources — landfire (GeoTIFF), usgs3dep, ncep_ncar. | Same update-binding pattern. | **Pattern reuse** | Confirms `kind: "static"` is an established corpus pattern for time-invariant datasets — the right kind for the routing network file. |
| `EarthSciModels/components/wildland_fire/fuel_model_lookup.esm` | 1.0.0 | Anderson 13 fuel-model parameter lookup: `ifelse` + `interp.linear` over integer codes. | code → σ, w₀, δ, Mₓ, h, T_f. | **No overlap** (idiom note) | Shows the corpus idiom for small integer-keyed tables without `function_tables`. MOSART mosart1.1.13 has no such tables (no Manning-n map, no reservoir tables). |
| `EarthSciModels/couplings/earthscidata_spatial_params.esm` | 1.0.0 | Coupling-interface stubs + `variable_map` param_to_var edges (GEOSFP/ERA5 → CEDS/EDGAR/NEI). | variable_map lat/lon/lev. | **Pattern reuse** | The stub + variable_map idiom for `mosart_coupling.esm` (map `Flrl_rofsur/rofsub/rofgwl/rofi/irrig` onto qsur/qsub/qgwl/qirrig parameters; export `Forr_rofl/rofi`, `Flrr_flood/volr/volrmch`). |
| `EarthSciModels/components/environmental_transport/puff.esm`, `gaussian/*`, `boundary_layer_mixing_kc.esm`, `plume_rise/*` | 1.0.0 | Plume/dispersion advection. | — | **No overlap** | Grepped; no river/overland content. |
| `EarthSciModels/lib/solar.esm`, `components/geodynamics/*`, `components/vegetation/*`, `components/wildland_fire/terrain_slope.esm` | 1.0.0 | Various. | — | **No overlap** | terrain_slope computes slope from elevation — MOSART slopes come from the network file instead. |

Full-corpus grep for `manning|kinematic|saint-venant|river|routing|overland|runoff|watershed|channel|hydrology`
confirms the only hydrology hits are the surface_runoff trio, the CLMU hydro
set, and false positives (chemistry "runoff"-adjacent names, deposition).
**There is no river-routing component in EarthSciModels.**

### 2.2 EarthSciDiscretizations (all files esm 1.2.0)

| Item | What exists | Relevance |
|---|---|---|
| `grids/cartesian_uniform_1d/` (+rules/stencils: upwind1, central, PPM, WENO, flux limiters, godunov, varcoeff laplacian) | Regular-grid 1-D rules. `upwind1_D_interior` is the stencil the Saint-Venant adapter imports. | A kinematic-wave *PDE* along a river (if ever needed, e.g. mizuRoute later) would build on the upwind1 family. MOSART itself is storage-based per cell — **no PDE grid needed**. |
| `grids/cartesian_uniform_2d/`, `cartesian_nonuniform_1d/`, `column_nonuniform_1d/`, `latlon/`, `latlon3d/`, `duo/`, `mpas/` | Grid families with rules/stencils; MPAS has ragged `edgesOnCell`-style member gathers (`index(dvEdge, index(edgesOnCell, i, k))` idiom in problems/). | The MPAS ragged-gather idiom is the **exact mechanism** MOSART's dnID gather and basin membership need (esm-spec §4.3.1 ragged ranges). |
| `problems/*` | MMS drivers: advection, heat, diffusion, SWE on duo, MPAS operators, eikonal. | `duo_swe` is the shallow-water system — not river routing. |
| grep `cfl|substep` | Only false positives (PPM docs). | **No CFL-substepping template exists** — the Fortran's `numDT_r/numDT_t` subcycling has no ESD analog. |
| grep `river|routing|kinematic|manning|channel` in docs/problems | False positives only (eikonal/godunov "channel"-adjacent text). | **No network/drainage/river grid family exists.** |

Gap: a `grids/network` (or `grids/drainage`) family — ID/dnID topology, cell→downstream gather, basin/outlet ragged sets — does not exist and would be a candidate ESD PR *if* we want grid-level reuse; today the `faq` gather + ragged index sets express everything without a new grid.

### 2.3 EqIce `stage2/` (the stage-2 pattern reference)

- **Factoring:** every physics component (`hydrology.esm`, `sia.esm`, `ssa.esm`,
  `ice_energy.esm`, ...) is one component importing a sibling
  `*_templates.esm` library by reference (`expression_template_imports` with
  `ref: "./hydrology_templates.esm"`). Templates carry `params`/`match`/`body`
  per pointwise law (e.g. `routing_conductivity`: `K = k W^(α−1) (Π+ε²)^((β−2)/2)`,
  `routing_upwind_flux`, `routing_flow`, mask-aware edge averaging). Free names
  in template bodies resolve against the importing model's parameters
  (consumer-supplied-geometry contract). Template libs are 1.1.0.
- **Component form:** continuous ODE/PDE (e.g. `dW/dt = surface + basal −
  dWtill/dt − div(V W) + div(D grad W)`), *not* the Fortran's per-substep Euler
  — the C++'s 19 internal CFL-capped Euler substeps are documented as an
  "Euler-vs-ODE scheme gap" measured at up to ~4.4% and folded into the test
  tolerances. This is the precedent for MOSART's `Euler`/`DLevelH2R`/`numDT_*`
  subcycling.
- **Tests:** `id`/`description`/`time_span`/`assertions` (plus
  `parameter_overrides`, `initial_conditions`); assertions are
  `{variable, time, coords: {i, j, o}, expected, tolerance}` — per-cell
  **1-based index-space** coordinates (e.g. dome center `(16,16)`), `rel 1e-12`
  for pure functions, looser for scheme gaps. Instantaneous/derivative tests
  use `time_span: {start: 0, end: 0}` with `time: 0` (e.g.
  `staggered_operators_initial`). One-step evolution tests assert selected
  cells against dumped columns with the scheme gap documented per cell.
- **No `function_tables`/`data_sources`/`reaction_systems`** anywhere in
  EqIce stage2 — all reference data enters as test tuples extracted by
  `stage2/tools/*.py` scripts from `stage1/dumps/`. (Scripts used only for test
  extraction — matches AGENTS.md.)
- **Closest conceptual analog to MOSART:** `hydrology.esm`'s Routing (Shreve)
  block — a flux law + staggered discrete operators + upwind flux + mask
  branches — but on a regular staggered grid, whereas MOSART is ID/dnID
  network based. The **template-factoring pattern transfers; the operators do
  not.**

### 2.4 EarthSciAST (esm 1.3.0 spec)

Version state: corpus is at 1.0.0/1.1.0 (EarthSciModels) and 1.2.0 (ESD);
PLAN.md §2 requires 1.3.0 for new documents. `aggregate` is deprecated in
favor of `faq` (esm-spec §4.2) — new MOSART files must use `faq`.

---

## 3. esm 1.3.0 feature applicability for MOSART

### 3.1 `function_tables` (esm-spec §9.5) — **NOT needed; use expression templates**

Schema requirements (quoted): each entry carries `axes`, `data`, and optionally
`description`, `interpolation`, `out_of_bounds`, `outputs`, `shape`,
`schema_version`; "`axes` is an ordered list of named `FunctionTableAxis`
entries... v0.4.0 caps axis count at 2"; each axis declares `name`, `values`
("strictly-increasing finite floats, ≥ 2 entries"); "`interpolation` is
`linear` (1 axis), `bilinear` (2 axes), or `nearest` (1 axis...)"; data shape
`[len(outputs), len(axes[0].values), len(axes[1].values), ...]`; lookup via the
`table_lookup` op with `axes` map and optional named `output`, lowering to
`interp.linear`/`interp.bilinear`/`index` (§9.5.3).

**Recommendation:** the GR* geometry functions and the Manning laws are
**closed-form analytic functions of storage** — they are expression templates
in a `mosart_physics_templates.esm` library (EqIce style), not tables. Nothing
in mosart1.1.13 routing is tabulated: Manning's n is a *per-cell parameter from
the network file* (a data source, §3.2), not a lookup. (If a later MOSART
version or a reference compset adds reservoir/release lookup tables, §9.5 is
ready for them — but the 2-axis cap and monotone-axis requirement fit such
tables fine.)

### 3.2 `data_sources` (esm-spec §8) — **YES: the `frivinp` network file**

- Kind (quoted): "`kind` | ✓ | Structural kind: `"grid"` (gridded array
  source), `"points"` (scattered point/station source), or `"static"`
  (time-invariant source). Any grid geometry the source reads — coordinates,
  connectivity, metric arrays — arrives as ordinary parameters bound to its
  file variables...". The routing network is time-invariant → **`kind: "static"`**
  (corpus precedent: landfire/usgs3dep/ncep_ncar). `temporal` is optional and
  should be omitted.
- File discovery (§8.2): `source.url_template` (with `{var}`-style and custom
  substitutions; relative paths resolve against the document directory).
- **netCDF support: yes.** The schema is format-agnostic; format decoding is
  the runtime reader's job (`reader_options` is "a free-form object handed to
  the format reader verbatim", §8.9.1). The corpus reads `.nc` directly (era5:
  `url_template: "...era5_pl_{date:%Y}_{date:%m}.nc"` with `time_variable`), so
  a static netCDF routing file is within the established reader capability.
- Parameter binding (§8.5, quoted): "`<param_name>: {type: parameter,
  units: string # required; units as seen by the model, shape: [index_set,
  ...] # required for a data update (§5.4), update: {kind: data, source:
  string # required; a data_sources key, from: {file_variable: string #
  required; name in the source file, unit_conversion: number | Expression #
  optional, ...}}}`".
- **Recommendation:** one `data_sources` entry `mosart_routing_network`
  (kind static) with each frivinp variable bound to a shaped parameter:
  `frac, fdir (→ mask), ID, dnID, area, areaTotal, hslp, nh, tslp, twid, tlen,
  nt, rslp, rwid, rwid0, rdep, rlen, nr, gxr` — all `shape: ["cell"]`
  (`unit_conversion` where the file stores e.g. degrees vs radians — none
  expected here; slopes/lengths are SI already). The `cell` index set
  (`index_sets`, §11.1) is sized by the file via `extent:
  {metaparameter: "N_CELLS"}` (§8.9.4) or a fixed size + `select` range.
  The two route handles (`rh_eroutUp`, `rh_direct`) and `areaTotal2` upstream
  accumulation are then **`faq` gathers over the shaped parameters**, not
  separate data: `eroutUp[i] = index(erout, index(dnID, i))`; basin membership
  for `rh_direct` is a **ragged index set** (`kind: ragged`, §4.3.1) whose
  members/offsets can be built at build time by value invention
  (`skolem`/`distinct`/`rank` over `dnID`/`fdir` — the MPAS `edgesOnCell`
  idiom: "index(dvEdge, index(edgesOnCell, i, k))"). Note the normative
  `ragged_values_not_gathered` rule — a ragged-range body must read the
  set's `values` array explicitly.

### 3.3 `solver` (esm-spec §2.2, `external`/`time_span` at 1.3.0) — **declarative hints; step control stays out**

- Fields (quoted table): `stiffness` `"low"|"moderate"|"high"`; `abstol`,
  `reltol` (integration tolerances, distinct from assertion `tolerance` —
  §2.2.2); `splitting` `"none"|"lie"|"strang"` ("Carries **no** prescribed
  substep structure"); `external` (1.3.0): "A reference to a run-configuration
  document... Resolved by §2.2.6. Advisory for every EarthSciAST binding,
  which never opens it"; `time_span` (1.3.0): "{start, end}... Each bound is a
  number in the document's time units or an RFC 3339 UTC timestamp".
- §2.2.3 explicitly excludes "A schedule, inline" — the algorithm per operator,
  splitting schedule and per-operator tolerances live in the referenced
  run-configuration document. **The MOSART CFL subcycling (`numDT_r/numDT_t`,
  `DLevelH2R`, `delt_mosart`) therefore cannot be declared as solver step
  control in the .esm.**
- **Recommendation:** `mosart_routing.esm` declares `solver:
  {stiffness: "moderate", abstol: ..., reltol: ..., time_span: {start: 0,
  end: <coupling period s>}}` (numbers in seconds avoid the time-origin
  machinery of §2.2.7). Two implementation routes for the subcycling, both
  with precedent: (a) **continuous-ODE form** (EqIce pattern) — storage
  equations `dwh/dt = qsur + ehout`, `dwt/dt = etin + etout`,
  `dwr/dt = erlateral + erin + erout + qgwl·area·frac` as ODEs, integrator
  step control handles stability, and the Euler-vs-ODE scheme gap is measured
  into the test tolerances; or (b) at stage 3, an EarthSciIntegrator run
  config with the per-cell subcycling expressed as operator partitioning,
  referenced via `solver.external`. Route (a) is the proven EqIce route for
  stage 2.

### 3.4 `reaction_systems` (esm-spec §7) — **confirmed NOT relevant**

Quoted: "Reaction systems provide a declarative representation of chemical or
biological reaction networks... the ODE form is derived automatically from the
reaction stoichiometry and rate laws." MOSART routing has no reactions;
**exclude** reaction_systems from the MOSART files (keep in mind only for CAM
chemistry, as planned).

### 3.5 `tests` (esm-spec §6.6) — **per-cell time-series assertions**

- Test fields (quoted table): `id` (required, unique), `description`,
  `initial_conditions`, `parameter_overrides` (shaped values may be inline
  nested arrays; scalar broadcasts), `time_span` (required `{start, end}`),
  `tolerance` (per-test default), `expression_template_imports` (per-test
  injection, §9.7.10 — survives round-trip), `assertions` (required, ≥1).
- Assertion fields: `variable` (local or scoped), `time` (must lie in
  `[start, end]`), `expected`, `tolerance`, `coords` (PDE only; "Map from
  spatial index-set (dimension) name to a **1-based, fractional index-space
  position**"; `idx = ceil(c − 1/2)` nearest-index sampling), `reduce`
  (`integral|mean|max|min|L2_error|Linf_error`; `integral` = uniform-cell
  Riemann sum under unit total measure = mean over interval index sets).
- Pass predicate: `|actual − expected| ≤ max(abs, rel·max(|actual|,
  |expected|))`, non-finite actual fails (§6.6.3).
- **Per-cell time-series assertion on wt/wr/wh — the spelling:** declare
  `wt, wr, wh` as `unknown` with `shape: ["cell"]`; assert
  `{"variable": "wr", "coords": {"cell": <k>}, "time": <T>, "expected": ...,
  "tolerance": {...}}` for selected cells k (1-based; cell order must be the
  dumped index order — pin it in the dumps' README), or element-name form
  `"wr[<k>]"` for a pointwise assertion. Instantaneous pointwise-law tests
  (CRVRMAN/CREHT/GR*) use `time_span: {start: 0, end: 0}`, `time: 0`, and
  `parameter_overrides`/`initial_conditions` carrying the input tuples — the
  EqIce `staggered_operators_initial` shape. Note §6.6: **tests are dropped at
  a mount edge** — the subassembly `mosart_coupling.esm` must declare its own
  assertions; leaf tests run only when the leaf file itself is the target.
  §6.6.5 also allows `reference` from_file JSON snapshots (relative to the
  .esm) for field-scale checks — an option for full-network one-step dumps.

### 3.6 Recommended MOSART .esm structure (stage 2 draft)

```
stage2/
  mosart_physics_templates.esm     # expression_template library: CRVRMAN_nosqrt,
                                    #   CREHT_nosqrt, GRMR/GRHR/GRPR/GRRR/GRHT/GRPT,
                                    #   dry-cell clamps (ifelse vs TINYVALUE),
                                    #   negative-storage guards, irrigation/flood
                                    #   rate rules (EqIce *_templates.esm style)
  mosart_routing.esm               # esm 1.3.0: index_sets {cell}; unknowns
                                    #   wh, wt, wr (shape [cell]); observed
                                    #   geometry/velocity fields; data-fed shaped
                                    #   parameters (network); faq gathers for
                                    #   eroutUp/direct transfer; solver block;
                                    #   tests (pointwise D tests + per-cell
                                    #   one-step storage assertions + budget)
  mosart_coupling.esm              # coupling stubs: import Flrl_* → qsur/qsub/
                                    #   qgwl/qirrig; export Forr_*/Flrr_*;
                                    #   direct-transfer/bypass option logic
```

---

## 4. Gaps and risks

1. **No kinematic-wave or network-routing physics exists anywhere** — the MOSART
   routing component and its template library must be authored from scratch
   (against instrumented dumps per boundaries.md §2.6). Nothing in
   EarthSciModels computes a routing step; the PLAN's "related but NOT river
   routing" verdict for `surface_runoff` is confirmed.
2. **No network/graph grid in EarthSciDiscretizations.** The ID/dnID topology
   is expressible today with `faq` gathers + ragged index sets (no grid
   needed for the storage-based scheme), but if the team wants reusable,
   reviewable network primitives (cell→downstream gather, basin ragged sets,
   upstream accumulation), submit a `grids/network` (drainage) family to ESD
   as a PR — the MPAS ragged-gather idioms are the precedent to lift.
3. **No CFL-substepping template or step-control mechanism.** The Fortran's
   per-cell `numDT_*`/`DLevelH2R` subcycling cannot be declared in a `.esm`
   (solver block is advisory-only, §2.2.3). Use the continuous-ODE route
   (EqIce's documented Euler-vs-ODE scheme gap, measured tolerances) at stage
   2, and EarthSciIntegrator partitioning via `solver.external` at stage 3.
   The instrumentation must dump per-level values to disentangle the
   DLevelH2R averaging (boundaries.md §2.2 already flags this).
4. **Version drift:** corpus is esm 1.0.0/1.1.0 (EarthSciModels), 1.2.0 (ESD);
   new files are 1.3.0 per PLAN.md §2. `aggregate` is deprecated (use `faq`).
   Confirm the EarthSciModels CI (tools/run_esm_inline_tests.py / .jl) accepts
   1.3.0 documents before the MOSART PR (M2/M3).
5. **Branch-heavy Fortran code:** the TINYVALUE dry-cell clamps and the
   negative-storage guards (`wt/wr < TINYVALUE` re-clamps) must be reproduced
   as `ifelse` templates and tested on both branches; the Fortran's
   "negative channel storage" *warning* (non-aborting) suggests potential
   non-conservation under some `bypass_routing_option` settings — candidate
   for `bugs.md` at instrumentation time.
6. **`data_sources` runtime dependency:** the frivinp source makes local
   `./esm test` runs depend on the network file being present at the resolved
   URL/path (§8.2: a source that cannot be located "MUST NOT be skipped").
   Tests should either ship a tiny toy network fixture next to the .esm or
   override the file via the ingest-time URL override, so the leaf tests stay
   self-contained; full-network assertions use from_file snapshots (§6.6.5).
7. **`solver.time_span` timestamps need a time origin** (§2.2.7: domain
   temporal reference_time/start, else POSIX epoch). Use plain numbers
   (seconds) in MOSART tests to sidestep this.
8. **No lakes/reservoirs/groundwater module in mosart1.1.13** — confirmed
   against the physics file and boundaries.md §2.5; nothing to reuse, and any
   future lake/reservoir component (later MOSART versions) would be new.

---

*Prepared for the EqESM manager (read-only survey; no commits).*
