# MOM6 KPP stage 2 — column prefix sums as ESD cumulative integrals (`mom6_kpp.esm`)

Refactor slice of the MOM6 KPP stage-2 component: the column **prefix sums**
over the level index — the interface heights `z_inter`, the cell heights
`z_cell`, and the hydrostatic reference pressure `pRef` — were written as
`faq` prefix reductions; this slice rewrites them as **ESD `integral`
operators** lowered by the EarthSciDiscretizations `column_nonuniform_1d`
cumulative-integral rules, per the house rule that *integrals are
`integral`, not aggregate sums or `faq` index gathers*.  All 158 existing
assertions pass **unchanged** (same expected values, same rel 1e-12
tolerances); the whole `./esm test stage2/` suite is unchanged at
**3073P / 0F / 0E** (baseline and after).

The parent component (mom6_kpp.esm, KPPBLD + KPPCoefs, 8 tests / 158
assertions) and its pointwise law library (mom6_kpp_templates.esm, all
pointwise, no sums) are described in `README-mom6.md`; this README covers
only the integral refactor.  F90 line numbers cite the instrumented fork
`stage1/EqESM-cesm-mom6` (branch `instrument/mom6-dumps`,
MOM_CVMix_KPP.F90, run configuration C_JRA_t62_mom6_inst).

## File map

| File | Change |
|---|---|
| `mom6_kpp_discretization.esm` | **new**: the layer/interface grid of the MOM6 column (index sets `layer`/`iface`, metaparameters `N_LAYER`/`N_IFACE`) and the two EarthSciDiscretizations rules imported by reference at a 9.7.7 renaming edge |
| `mom6_kpp.esm` | the three prefix sums rewritten as `integral`; `term` dropped; new named intermediate `depth`; `unit_layer` integrand; `z_sfc`/`z_top` bounds declared; file-level `index_sets`/`metaparameters` moved into the discretization file |
| `mom6_kpp_templates.esm` | unchanged — audited: every template is pointwise (no `faq`, no prefix sums) |
| `README-mom6-kpp-integrals.md` | this file |

The rules themselves are **not** copied: `mom6_kpp_discretization.esm`
imports, by reference, `integral_lev_nodes_cumulative_from_bottom` and
`integral_lev_cumulative_from_bottom` from
`$ESD_ROOT/grids/column_nonuniform_1d/rules/` (the same library
`cism_velo_discretization.esm` imports its DIVA integrals from; here the
surface-to-interface / surface-to-layer-centre forms instead of the
interface-to-top form).

## What was rewritten, and how

The Fortran (KPP_compute_BLD, MOM_CVMix_KPP.F90) computes the geometry and
pressure as forward running loops over `k`:

```
iFaceHeight(1) = 0
cellHeight(k)    = iFaceHeight(k) - 0.5*dh(k)      (F90:1237)
iFaceHeight(k+1) = iFaceHeight(k) - dh(k)          (F90:1238-1243, dh = dz here: min_thickness = 0)
pRef = p_surf                                      (F90:1282-1283)
... per k: use pRef, then pRef = pRef + PRESS_STEP*h(k)   (F90:1452-1456)
```

The `.esm` equations now read as the integrals these loops discretize:

```
depth(i)  = ∫_{z_sfc}^{iface i} dz          = sum_{j<i} dz(j)        (integral_lev_nodes_cumulative_from_bottom)
z_inter(i)= -depth(i)
z_cell(k) = -∫_{z_sfc}^{zc(k)} dz           = -(sum_{j<k} dz(j) + dz(k)/2)   (integral_lev_cumulative_from_bottom)
pRef(i)   = p_surf + PRESS_STEP·depth(i)    (Boussinesq hydrostatic pressure at the interface depth)
z_top     = depth(N_IFACE)                  (the bottom bound name, measure-free in the lowering)
```

Mechanics, all following the `cism_velo_discretization.esm` precedent:

- **Grid and metaparameters.** The file-level `index_sets` (`layer`,
  `iface`) and `metaparameters` (`N_LAYER`, `N_IFACE`) moved from
  `mom6_kpp.esm` into `mom6_kpp_discretization.esm`, which is imported by
  both models (KPPBLD and KPPCoefs — the latter only needs the `iface`
  axis).  The `iface` axis is sized `N_LAYER + 1` symbolically so it is
  deep-equal to the renamed grid axis of the imported rules.
- **Renaming edge (esm-spec 9.7.7).** The rules are imported with
  `rename: {lev: layer, lev_nodes: iface}` and `bindings: {NLEV: N_LAYER}`,
  so the integral operators are spelled with the MOM6 axis names
  (`var: layer, lower: z_sfc, upper: iface`); esm-spec 9.7.7 rewrites the
  rule match patterns (`var`/`lower`/`upper`), bodies and where-shapes at
  the import edge (verified with a scratch probe before authoring).
- **Measure.** The rules' free-name measure `dz` is the model's existing
  `dz` parameter (the layer heights, thickness_to_dz), which is exactly the
  Fortran's `dh` in the instrumented run (`min_thickness = 0`, no thickness
  inflation).  The Fortran's `pRef` loop actually accumulates `h`; `h = dz`
  in the run, so the measure carries the same values (noted in the
  variable descriptions).
- **Unit integrand.** The rules' `where` clauses require the integrand to
  be a declared layer-shaped variable (a scalar `1` matches nothing and the
  `integral` survives unlowered — probed).  The depth *is* the integral of
  the measure, so the integrand is the new `unit_layer` field (identically
  1 on every layer): `depth = integral(unit_layer, ...)` is the ESD idiom
  for `∫ dz`.  Because the integrand is exactly 1.0, the lowered products
  `f[j]·dz[j] = dz[j]` are exact and the prefix scan is bit-identical to
  the Fortran loops (see next section).
- **Named intermediate.** `depth` is the shared integrand result of
  `z_inter` and `pRef` (the house rule: never repeat a subexpression
  inline).  The old `term` array (a `makearray` pinning `p_surf` at index
  1) is dropped — `pRef` is now the Boussinesq relation
  `p_surf + PRESS_STEP·depth`.

## Summation order vs the Fortran — measured residuals

The task rule: *if the rewritten form differs from the Fortran's summation
order beyond rel 1e-12 on any pin, report the measured residual and the
reason*.  Measured with `scratch/verify_sumorder.py` (not committed — the
numbers below; it re-derives both accumulation orders in float64 from the
five KPPBLD tests' inline inputs, 75 layers × 5 columns):

| quantity | worst |relative| residual vs the Fortran order |
|---|---|---|
| `z_inter` | **0.0** — bit-identical (negation and the 1.0 product are exact; the forward exclusive scan is the Fortran's loop order) |
| `z_cell` | **0.0** — bit-identical (`cellHeight(k) = iFaceHeight(k) − 0.5·dh(k)` is exactly `−(depth(k) + 0.5·dz(k))`) |
| `pRef` | **6.5e-16** relative — the only deviation |

The `pRef` deviation has a single cause: the Fortran accumulates the
per-layer products `PRESS_STEP·h(k)` onto the running `p_surf`, while the
integral form accumulates the layer heights first and applies `PRESS_STEP`
and `p_surf` after the prefix — both factored out of the per-layer
additions (the products `PRESS_STEP·dz(j)` are then not formed per layer,
and `p_surf` is added last instead of first).  The residual is ~7 orders
of magnitude inside the rel 1e-12 pin tolerance, and propagates through the
EOS (wright_rho at `pRef`) and into `BulkRi`/`N2`/`deltaU2` at the same
order — all 158 pins pass unchanged.  Nothing was loosened.

## What could NOT be rewritten (dump-coverage / rule-coverage gaps)

- **The surface-layer means `surfHtemp`/`surfHsalt`/`surfHu`/`surfHv` stay
  `faq` sums.**  Mathematically they are exact integrals of T/S/u/v over
  the depth interval `[0, SLdepth(i)]`, with `delH(i,j) =
  min(max(SLdepth(i) + z_inter(j), 0), dz(j))` the exact layer-overlap
  measure — but the upper bound `SLdepth(i)` is an *i-dependent interior
  depth* (the assumed OBL depth of each layer), and the five
  `column_nonuniform_1d` rules match only the ground bounds `z_sfc`,
  `z_top`, `lev`, `lev_nodes` (an integral over any other bound survives
  lowering and trips the `unlowered_operator` gate — esm-spec 9.6.1/9.6.8).
  They are also whole-column sums in the aggregation sense, so
  `integral_lev_whole` was considered: it cannot express the `SLdepth(i)`
  mask either (and the where-clause needs a layer-shaped integrand, while
  the masked integrand is `layer × layer`).  The `delH` closed form is
  documented in the file as bit-identical to the Fortran's running `hTot`
  accumulation (Sterbenz).  An ESD PR that generalizes the integral rules
  to a parameter-valued interior bound would be the future home.
- **Whole-column sums (`integral_lev_whole`)**: audited — the only
  whole-column sums in the file are the four `surfH*` masked sums above;
  there are no other level-index sums (`kfirst`/`kwfirst` are `min`
  reductions, not sums).
- **`z_cell`** was rewritten even though it is not itself a prefix sum: its
  closed form *is* the cumulative midpoint integral, and the rewritten form
  is bit-identical to the Fortran `cellHeight` (the old
  `0.5·(z_inter(i)+z_inter(i+1))` average differed by up to 1 ulp).

## frazil `p_top` — how the two hydrostatic prefixes should share one law

`mom6_frazil.esm` (merged by the mom6-frazil-bflux lane) computes `p_top =
p_surf + Σ_{j<i} dP(j)` with `dP = g_Earth·RHO_0·h` — the *same* hydrostatic
prefix as `pRef` (KPP's `PRESS_STEP = g_Earth·H_to_RZ = 9.8·1035` equals
frazil's `g_Earth·RHO_0` under Boussinesq).  When both components are
merged:

1. `mom6_frazil.esm` should import `mom6_kpp_discretization.esm` (or the
   discretization file should be renamed/relocated to a shared
   `mom6_discretization.esm`) and write `p_top = p_surf + integral(...)` in
   the same `∫_{z_sfc}^{iface} dP` form — the rules lower it to the exact
   forward scan, no new machinery.
2. The per-layer pressure increment should be factored once — the natural
   home is a single template in `mom6_templates.esm` (already imported by
   both components), e.g. a `hydrostatic_layer_increment(h)` template, or
   the two components should at least share the *same* named integral
   equation shape so a reviewer can diff them.  KPP's `pRef` currently
   absorbs the increment into `PRESS_STEP·depth`; frazil's `dP` keeps it as
   a layer field — the merged form should pick one (the layer-field form
   preserves the Fortran's per-layer products exactly and is the
   recommended one for a future shared law).

## Pins and coverage

Unchanged from the parent component: the 5 KPPBLD tests pin `mom6_kpp_obl.csv`
rows (cells 782/883/784/1184/684 at step 3600/10800) for `BulkRi`/`deltaU2`/
`N2` at layers 2..kOBL+3, spot `z_cell` at layers 1..6 and 75, `z_inter` at
ifaces 2..6 and 75, and `OBLdepth`/`kOBL`/`kfirst`; the 3 KPPCoefs tests pin
`mom6_kpp.csv` rows for `Kd_heat`/`Kv` at interfaces 2..kwup plus one
pass-through `Kv_out` at kwup+1.  The integral forms are exercised by every
`z_cell`/`z_inter` pin (and indirectly by every EOS-fed pin via `pRef`).
Unpinned branches (unchanged from before, see `README-mom6.md`): the
`OBL_limit` branch, the kfirst = 2/3 interpolated OBL depths (floored to
`MINIMUM_OBL_DEPTH` in every trace), the unstable `surfBuoyFlux < 0`
branches, the Langmuir paths, and the `MINIMUM_VT2` floor.  The rules'
`z_top` bound is never exercised numerically (measure-free), and the
`N_LAYER`/`N_IFACE` defaults (75/76) are used by every test.

## CLI notes

- A template-library file must declare a top-level `expression_templates`
  key even when it carries no templates of its own (the schema refuses a
  library without it); `mom6_kpp_discretization.esm` declares an empty
  block.  It has no tests; `./esm validate` + the consumers' tests cover it.
- `./esm units stage2/mom6_kpp.esm` reports no unit mismatches (`unit_layer`
  is dimensionless, the depth integral is `m`, `pRef` is `Pa`).
- Commit: `23d230f` on branch `wt/mom6-kpp-integrals` (whole suite
  3073P/0F/0E before and after).
