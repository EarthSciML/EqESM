# CICE stage-2 PR preparation — session notes (2026-11-01/02)

All work done in sibling worktrees under `/projects/illinois/eng/cee/ctessum/ctessum/code/`.
Nothing pushed; no PRs opened; no `main` edited. EqESM repo otherwise untouched
(no PROGRESS.md / jobs.json / bugs.md edits).

## Branches

| Repo | Branch | HEAD | Contents |
|---|---|---|---|
| EarthSciDiscretizations | `cice-bgrid-rules` | `a0c0421` | `755a289` (origin/main) + cherry-pick of `11fe246` (`feat(rules): CICE B-grid flux-divergence and stress-divergence lowering rules with conformance goldens`, author Christopher Tessum). `PR_DESCRIPTION.md` untracked. |
| EarthSciModels | `cice` | `f7e05c5` | `80f303c` (origin/main) + `components/cryosphere/cice/` (6 .esm + README.md, 20080 insertions). `PR_DESCRIPTION.md` untracked. |
| EqESM | `wt/cice-upstream-prep` | `754af41` (unchanged) | this notes file (scratch/) + nothing else |

## Placement and ESD-reference decisions

- **Placement**: `components/cryosphere/cice/` — new domain directory following the
  corpus `<domain>/<component>` pattern (`environmental_transport/mosart`); no
  pre-existing cryosphere content in the corpus.
- **ESD refs**: `${ESD_ROOT}` refs kept as-is (corpus convention). The rule files
  declare no metaparameters/index_sets, so bindings-free imports cannot collide
  with the components' local `N_I/N_J/NCAT` declarations; grid metrics are local
  parameters carrying the dump values (grid.esm is deliberately NOT imported).
- **Ordering dependency**: EarthSciModels CI clones ESD `main` for ESD_ROOT —
  which lacks `cice_bgrid` until the ESD PR merges — so the CICE PR's CI is red
  on the two lowered files + the assembly until then. Documented in both PR
  descriptions; merge the ESD PR first.
- **Metadata**: MOSART pattern (`tags`/`created`/`authors`/`references`; created
  `2026-11-01T00:00:00Z`; authors `["EqESM project contributors"]`; Hunke et al.
  2023 LA-UR-23-31851 + CESM URL with per-file F90 lists + EarthSciModels URL;
  transport adds Dukowicz & Baumgardner 2000 and Lipscomb & Hunke 2004). The two
  files that already carried stage-2-era license/tags/references (dynamics,
  transport) got their block replaced, `license: MIT` kept.
- **Description fixes**: `stage2/README-cice.md` → `README.md`;
  `README-cice.md` → `README.md`; `RangeSpec round-trip parse bug, see
  PROGRESS.md)` → `RangeSpec round-trip parse bug of the esm CLI)`. `bugs.md`
  citations kept verbatim; README gained the EqESM-repo clarification line.
- **Identity**: conform script's STRIP-COMPARE prints IDENTICAL for all 6 files
  (equations, tests, every number); raw diffs vs stage2 confined to description
  strings and the metadata block. `README.md` = stage2 README with the Testing
  section re-pointed at `esm test components/cryosphere/cice` (134 assertions,
  CICE-only) + the clarification line.

## Runner results (EarthSciAST `main`, 0.8.0-era bindings)

### Rust CLI (`esm test components/cryosphere/cice`, main-built 0.8.0 binary)
- ESD_ROOT = ESD `main` (CI default): **99P / 0F / 3E** — coupling 69P, thermo
  30P; dynamics, transport, and the assembly error with
  `template_import_unresolved` (ESD main has no cice_bgrid).
- ESD_ROOT = `cice-bgrid-rules` worktree: **134P / 0F / 0E** (69+30+19+16).
- `esm validate`: all 6 files pass on 0.8.0 with the ESD branch.

### Python gate (`tools/run_esm_inline_tests.py --files ...`, earthsci_ast 0.8.0)
- ESD_ROOT = ESD `main`: **24P / 0F / 79E** (assembly 0P/1E — load;
  coupling 0P/69E; dynamics 0P/1E — load; templates 1P;
  thermodynamics 23P/7E; transport 0P/1E — load).
- ESD_ROOT = `cice-bgrid-rules` branch: **42P / 0F / 94E** (assembly 1P;
  coupling 0P/69E; dynamics 2P/17E; templates 1P; thermodynamics 23P/7E;
  transport 15P/1E). The lowered files load and most transport assertions
  pass; every error is the runner bug below (the failing dynamics/transport
  tests override shaped parameters that the solve then cannot resolve).
- **Runner bug (minimal repro)**: a test `parameter_overrides` with an ARRAY
  value for a shaped parameter (esm-spec §6.6.2) is routed by
  `inline_tests.simulate_states` through `esm_problem(p=...)`, whose `p`
  channel carries only scalars — the array is silently dropped (the array
  channel is `const_arrays`, never fed by the runner). Consequences:
  * every coupling test errors at solve — `Unresolved symbol: 'IceCap.aicen'`
    or `index applied to scalar value` (69E);
  * every §6.6.5 array-state/observed assertion errors — `array state 'zTin'
    has no cells in var_map` (thermodynamics 7E); the observed-replay context
    additionally omits `input_arrays`, so observeds reading shaped parameters
    are dropped by `skip_unresolved`.
  Repro `p_probe.esm` (1 model, 1 shaped param, 1 faq-index equation, 1 test
  with array override): rust `esm test` 1/1 PASS; python `solve failed:
  Unresolved symbol: 'M.aicen'`. The form is blessed by the python package's
  own conformance tests (`shaped_parameter_broadcast`,
  `shaped_parameter_const_arrays`) and used corpus-wide (`slab.esm` same faq
  form).

### Julia tree-walk (`run_esm_tests`, EarthSciAST.jl 0.8.0 dev'd, ODE solver loaded)
Per-file walk of the 4 ESD-independent-of-coupling documents
(dynamics/thermodynamics/transport/templates in a scratch dir), plus the
assembly and coupling from the full-walk runs:

| document | ESD `main` | `cice-bgrid-rules` branch |
|---|---|---|
| cice.esm (assembly) | FAIL (1 row, load: missing cice_bgrid) | not walked (see coupling note); all mounts resolve — python/rust confirm |
| cice_coupling.esm | **NO COMPLETION** | same (ESD-independent) |
| cice_dynamics.esm | 0P/0F/1E (load) | 0P/0F/19E — `E_TREEWALK_UNSUPPORTED_SHAPE: Dynamics.sigma12` |
| cice_thermodynamics.esm | 0P/0F/30E — `E_TREEWALK_UNSUPPORTED_OP: index reached _compile unresolved` | 0P/0F/30E (same; ESD-independent) |
| cice_transport.esm | 0P/0F/1E (load) | **16P/0F/0E** (12.9 s) |
| cice_templates.esm | load OK (0 rows) | load OK (0 rows) |

Coupling measurement: coupling + templates in a scratch dir, 40-min hard cap —
`timeout` fired with **no per-file row emitted**: 16.2e9 allocations, 7048 GC
cycles, 5.6 GB RSS at the cap (the earlier full-walk attempt reached 14.6 GB
RSS / 63 CPU-min / 297M minor faults before being killed). The document is
ESD-independent, so this is a Julia-walker scalability finding for the file
itself: one document past the CI shard's 60-min job cap takes the shard down
— exactly the failure class the CI comments warn about. The Julia binding
also refuses the dynamics' shaped B-grid stress state and `index` ops in the
thermodynamics compile (`E_TREEWALK_UNSUPPORTED_SHAPE` / `_UNSUPPORTED_OP`),
while the Rust binding runs all of it (134P/0F/0E).

### ESD conformance (all five bindings, re-verified after the v0.8.0 flatten change)
- julia / python / rust / typescript / go, `ast` category, at EarthSciAST main
  (v0.8.0): **169 passed / 0 failed / 0 skipped** each, including the six
  cice_bgrid fixtures (golden-digest sha256 matches). Final cross-binding
  compare (all five in one invocation) re-run under 0.8.0.
- v0.8.0 `feat(flatten)!` (pointwise-lifted states name their index sets,
  EarthSciAST#623/#624) does NOT change the cice_bgrid goldens.
- validate: 645 L008 findings on the library are pre-existing schema-`$id`
  drift (1.2.0-declaring corpus vs 1.3.0 schema; identical at ESD main HEAD);
  zero findings on the cice_bgrid files. Lint fixtures 8/8, convergence 112/112.

## RangeSpec round-trip bug — re-check (task requirement)

**Still present on both CLIs.** `esm validate` and `esm test` fail with
`JSON parse error: data did not match any variant of untagged enum RangeSpec`
for a faq range bound as a string metaparameter reference (`[2, "N_I"]` or
`{"from": "N_I"}`) when the metaparameter is not declared in the importing
document itself. Verified on the local 0.7.0 binary and the main-built 0.8.0
binary (identical behavior).

Variant matrix (repro files in `/scratch.local/ctessum/cice-upstream-prep/rangespec/`):

| Variant | faq range | N_I | import | result |
|---|---|---|---|---|
| main.esm | `[2, "N_I"]` | local (default 3) | tpl (no metas) | ✓ |
| main2.esm | `[2, "N_I"]` | via import, bound 3 | lib with bindings | ✗ |
| main3.esm | `[2, 3]` literal | via import, bound 3 | lib with bindings | ✓ |
| main4.esm | `{"from": "N_I"}` | via import, bound 3 | lib with bindings | ✗ |
| main5.esm | `[2, "N_I"]` | via import, open | lib no bindings | ✗ |
| main6.esm | `[2, "N_I"]` | local (default 3) | lib with bindings | ✓ |
| main7.esm | `[2, "N_I"]` | local (default 3) | lib no bindings | ✗ |

Reading: the string ref round-trips only when the named metaparameter is
closed *in the importing document's own scope* AND not re-exported open by an
unbound import (main7: the unbound lib import re-exports its open N_I,
colliding with the local one). Import-edge bindings alone (main2) are not
enough. The CICE workaround — local declarations + imports of rule files that
declare no metaparameters — is therefore still REQUIRED on current CLIs, and
the current CICE files pass validate + test on 0.7.0 and 0.8.0.

## Upstream finding (no change made)

EarthSciModels `Project.toml` at HEAD: `[compat] EarthSciAST = "0.3"` while the
Julia package on `main` declares 0.8.0 (and now v0.8.0 is tagged) — Julia CI
resolve fails at HEAD for every PR, the same failure mode the file's own
comment documents. The CICE PR deliberately carries NO Project.toml change
(avoids conflict with MOSART PR #33's bump); verified locally via a temporary
uncommitted compat bump to `"0.8"` + `EARTHSCI_SERIALIZATION_PATH` dev from
the local EarthSciAST main checkout (restored afterwards).

## Scratch artifacts (all under /scratch.local/ctessum/cice-upstream-prep/)

- `conform_cice.py` — the stage2 → corpus conform script (STRIP-COMPARE
  IDENTICAL; handles both metadata layouts: insert-after-description and
  replace-existing-license/tags/references).
- `p_probe.esm` — python-runner array-override minimal repro (rust 1/1,
  python E).
- `rangespec/` — lib.esm + main.esm..main7.esm RangeSpec variant matrix.
- `julia-coupling-only/` — coupling + templates for the timed walk.
- Logs: julia-resolve-08.log, julia-walk-esdmain.log, julia-coupling-only.log,
  esd-conformance-all-08.log, rust-rebuild-08.log, ts-rebuild-08.log.

## All checks complete (no open items)

- ESD all-bindings conformance at 0.8.0: **1521 passed / 0 failed / 0 skipped**
  (julia, python, rust, typescript, go; golden-digest + cross-binding-bytes;
  identical to the 0.7.0-era result — the v0.8.0 flatten change did not move
  the cice_bgrid goldens).
- Julia table above complete; scratch notes committed to `wt/cice-upstream-prep`
  under `scratch/`.
