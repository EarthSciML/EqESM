# MOM6 implicit vertical tracer diffusion (mom6_vertical_diffusion.esm)

Stage-2 lane for the MOM6 vertical mixing of tracers by the implicit
tridiagonal solver of `tracer_vertdiff_Eulerian` (the no-sink branch,
MOM_tracer_diabatic.F90:384-413), driven by `MOM_diabatic_driver.F90:1728-1749`
(entrainments + the two solve calls for T and S), pinned to
`stage1/dumps/mom6/mom6_diabatic_tri.csv` (156 rows; columns
`cpl_step,ymd,tod,gi,gj,i,j,lon,lat,dt, ent_t_01..76, ent_s_01..76,
Kd_heat_01..76, Kd_salt_01..76, h_01..75, T_before_01..75, T_after_01..75,
S_before_01..75, S_after_01..75`).

## Component

`TridiagonalMixing` in `stage2/mom6_vertical_diffusion.esm` writes the
MOM6 operator as the PDE

    dT/dt = d/dz ( Kd_heat dT/dz ),   dS/dt = d/dz ( Kd_salt dS/dz )

discretized by the EarthSciDiscretizations rule
`varcoeff_face_laplacian_lev_flux_bc` (priority 10) with zero prescribed
face fluxes (`kdudz_bot = kdudz_top = 0`, units-free per the cam idiom
because the T and S fluxes carry different units).  On the esd column
(lev 1 = bottom .. lev NLEV = surface, `ze` increasing upward from
`ze[1] = -H` to `ze[NLEV+1] = 0`) the rule's lowering is

    G[n]   = Kd[n] (u[n] - u[n-1]) / (zc[n] - zc[n-1])     (interface flux)
    L[ell] = (G[ell+1] - G[ell]) / dz[ell]                  (cell tendency)

with the face fluxes zeroed, which is EXACTLY the MOM6 tridiagonal law

    h_k (u'_k - u_k) = ent_k (u'_{k-1} - u'_k) + ent_{k+1} (u'_{k+1} - u'_k)
    ent_k = dt Kd_k / (dz_neglect + 0.5 (h_{k-1} + h_k))

(the PDE's h divides out through the cell measure).  The esd column
reverses the MOM6 top-down order: esd lev ell <-> MOM6 layer NLEV+1-ell,
esd interface n <-> MOM6 interface NLEV+2-n.  `dt*Kd/dz_int` and the
tridiagonal residual are factored as `interface_entrainment` and
`tridiagonal_residual` in `stage2/mom6_templates.esm` (applied once per
field: T with Kd_heat/ent, S with Kd_salt/ent_s), so the T and S branches
differ only in the bound diffusivity field.  T and S also carry `D(u,t)`
ODEs (ptend) for stage-3 time stepping; the inline tests use
`time_span {0,0}` so the ODEs are not integrated.

## Pinned rows

| test id | row (tod, gi, gj) | dt | character |
|---|---|---|---|
| rowA_convective | (32400, 3, 79) | 3600 s | deep convection, ent ~ 3.6e6 m on interfaces 61-75, mixed layers 70-71 |
| rowB_moderate | (10800, 5, 79) | 3600 s | ent spike 831 at interface 73, mixed layers 72-73 |
| rowC_shallow | (3600, 3, 78) | 3600 s | shallow mixed layer at interface 24 (ent 719), deep column quiet |

Each row pins, per the wave-2 proposal:

- **ptend_T/ptend_S at the boundaries** (esd levs 74-75, surface; and the
  bottom, esd lev 1), **at a mid-column quiescent layer** (esd 40), and
  **at the mixed region** (esd 5 / 4 / 53).  The expected values are
  `(T_after - T_before)/dt` and `(S_after - S_before)/dt` from the dump
  columns — derived, not dump columns.
- **ent/ent_s at three esd interfaces** per row: the max-ent interface,
  a **T/S-differing-diffusivity interface** (Kd_heat != Kd_salt: esd 28
  MOM6 49, esd 40 MOM6 37 with ent 36% apart, esd 18 MOM6 59 with ent 46%
  apart), and a near-surface interface.
- **the tridiagonal identity** (`tri_res_T/S == 0`) at the mixed layer
  and the bottom layer, with abs tolerance `1e-15 * scale` where
  `scale = |lhs|+|e_in|+|e_out|+|ent_k u'_k|+|ent_{k+1} u'_k|` is the
  operand-magnitude scale of `validate_mom6_w2.py` (the actual scaled
  residuals at the pinned layers are <= 1.4e-16, below the dump-wide
  maxima 1.2e-12 T / 2.0e-12 S).

## Tolerances (all stated in the test descriptions)

Rel 1e-12 is used wherever the dump supports it (ent pins, surface and
mid-column tendencies, all of row C's shallow-region pins).  At the deep
convective layers the SOLVED doubles T'/S' satisfy the discrete identity
only to the solve's round-off, which convective entrainments
(ent ~ 3.6e6 m) amplify into a ~1e-4-relative mismatch of h (T'-T)
(residual ~ 1.6e-10 degC*m vs terms ~ 2.1e-6; perturbing T'_75 by one
ULP moves the ent term by ~ 3.6e6 * 2e-16 ~ 7e-10).  The scaled residual
stays <= 1.4e-16 — this is the double-precision conditioning limit of
the solve, not a defect, so the tendencies at those layers (esd 1-5
of row A) are not pinned (manager review: rel 1e-3 / abs 2e-8 K/s would
be uninformative); the scaled tri_res identity covers them.  The
deep-interface ent pins use rel 1e-9 because the esd column reconstructs
dz_int from ze, and subtracting ~1e-3 m from a ~743 m depth coordinate
rounds at the depth magnitude (~1e-13 m absolute -> ~1e-10 relative on
the deep spacings); `interface_entrainment` therefore forms the spacing
as `0.5 (ze[nlev+1] - ze[nlev-1])` (an exact subtraction) rather than
`zc[nlev] - zc[nlev-1]` (which would round at the cell-center magnitude),
and the residual pins are insensitive to the leftover ~1e-12 relative
spacing error.

## Coverage notes

- **Kd_heat != Kd_salt on most interfaces** (10-36 of 76 per column, up
  to ~50%): the T/S-differing-diffusivity branch IS exercised, both by
  the separate ent/ent_s pins at the differing interfaces and by the
  separate ptend_T/ptend_S pins.  This corrects the claim in
  `README-mom6-w2.md` that Kd_heat == Kd_salt on every dumped row — the
  dump shows otherwise (a documentation error there, not a CESM bug; no
  `bugs.md` entry).
- The surface/bottom **face fluxes are not pinned as separate values**;
  the zero-flux boundary conditions are exercised through the operator
  at the boundary layers (esd 74-75 and 1), which are exactly the
  `kdudz_top`/`kdudz_bot = 0` stencil points.
- `h_neglect`/`dz_neglect` (1e-30 m) are far below any tolerance here
  and are not carried.
- Missing inputs: none — the solve consumes only h, Kd_heat, Kd_salt,
  dt and the entering tracers, all in this dump; no other lane's output
  is required.

## CLI limits worked around

An implicit (backward-Euler) solve is not expressible in the esd inline
tests: algebraic self-referential equations are rejected (unresolved
symbol / observed_cycle), so the component carries the operator in ODE
form (`D(u,t) = ptend`), the dumped after-state as the initial
condition, and the tests assert the operator evaluated on that state
(ptend = (T_after - T_before)/dt, the ent identity, and the residual
identity) instead of asserting T' directly (which would only copy an
input back out).

## Files

- `stage2/mom6_vertical_diffusion.esm` — the component (NLEV = 75).
- `stage2/mom6_templates.esm` — appended `interface_entrainment` and
  `tridiagonal_residual` (append-only file; the Wright fragments are
  untouched).
- `stage1/analysis/extract_mom6_vertdiff_tuples.py` — extracts the test
  tuples from the dump (row selection, esd reversal, sanity checks:
  stencil-vs-dump at every pinned layer, ent identity, scaled residuals
  below the dump-wide maxima) and prints the parameter/IC arrays and
  expected values at full precision; run as
  `python3 extract_mom6_vertdiff_tuples.py [DUMP_DIR]`.
