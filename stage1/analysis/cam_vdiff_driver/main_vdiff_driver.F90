! Standalone reference driver for the CAM6 vertical-diffusion operator.
!
! Builds Coords1D from the step-1 dumped interface pressures pint (CAM
! order, lev 1 = model top), applies the finite-volume diffusion operator
! of linear_1d_operators.F90 with interface diffusivity K = dpidz_sq
! (kv = 1 m2/s; dpidz_sq is (dp/dz)^2 at the interfaces of the same dump)
! and the default BoundaryFixedFlux conditions at both ends (zero
! prescribed flux, matching the esm rule's kdudz_bot = kdudz_top = 0),
! to the step-1 entering wind u0 = u1 - ptend_u*ztodt (ztodt = 1800 s),
! and prints A*u0 in CAM order.
!
! The pinned esm rows reverse the output (esd lev i = CAM lev 33-i):
!   ptend_u(1)  = (A u0)[32] = 2.0633949733127771e-05
!   ptend_u(16) = (A u0)[17] = 1.7093444110831203e-07
!   ptend_u(32) = (A u0)[1]  = 8.6780995093086661e-67
!
! The DATA statements are transcribed by extract_step1_data.py from the
! step-1 rows of stage1/dumps/cam/cam_vdiff_prof.csv (17-digit dump
! precision, repr round-trip -- the same precision the .esm files use).
program vdiff_driver

  use shr_kind_mod, only: r8 => shr_kind_r8
  use coords_1d, only: Coords1D
  use linear_1d_operators, only: diffusion_operator, TriDiagOp

  implicit none

  integer, parameter :: pver = 32
  integer :: k
  real(r8) :: pint(pver + 1), dpidz_sq(pver + 1), u0(pver)
  real(r8) :: d_coef(1, pver + 1), u(1, pver), out(1, pver)
  type(Coords1D) :: coords
  type(TriDiagOp) :: op

! pint (33), dpidz_sq (33), u0 (32) -- step 1, CAM order:
  data pint / 225.5239523947239_r8, 503.16918920859297_r8, 1015.7947507913752_r8, 1855.53172375681_r8, 2973.467587052297_r8, 3927.3001624446247_r8, 4711.449936258011_r8, 5624.0499560208045_r8, 6680.049811850153_r8, 8070.141900625354_r8, 9494.104347324384_r8, 11169.321404306324_r8, 13140.144979858878_r8, 15458.726012713647_r8, 18186.422609487316_r8, 21102.669844630946_r8, 24533.534314275326_r8, 28569.892529385947_r8, 33318.839428376654_r8, 38906.68735723806_r8, 45482.53370868157_r8, 53223.66275062209_r8, 60975.2962748278_r8, 68106.55401087951_r8, 74141.05342121309_r8, 78630.18032864774_r8, 81230.48827047336_r8, 83671.13748645567_r8, 85932.03544898969_r8, 87994.35722354907_r8, 89840.70160218416_r8, 91455.24102284685_r8, 92848.96251956573_r8 /
  data dpidz_sq / 0.0009336295172983327_r8, 0.004947241250016872_r8, 0.022099894253797766_r8, 0.07724499263247359_r8, 0.20509861646498687_r8, 0.36712348714249965_r8, 0.5402358536147887_r8, 0.7862683097802157_r8, 1.1270109196163576_r8, 1.650387753496299_r8, 2.2657328412483975_r8, 3.2179501696022648_r8, 4.592483421476296_r8, 6.3105758356134185_r8, 8.540292915166466_r8, 10.95811449358474_r8, 13.81296648547239_r8, 17.336766583820914_r8, 21.817810549433215_r8, 27.573046994736767_r8, 35.128271573798365_r8, 45.29586785449068_r8, 56.76529478263576_r8, 68.21272417397428_r8, 78.05102141983093_r8, 85.42956372396445_r8, 89.5179241603476_r8, 93.59357425687253_r8, 97.39580873262373_r8, 100.90559719024924_r8, 104.13855413415965_r8, 107.28745785031799_r8, 110.41760768707236_r8 /
  data u0 / 3.676249032166472e-67_r8, 1.9244104698866517e-59_r8, 2.0577267759704973e-52_r8, 1.0198359578029224e-45_r8, 3.250702812670074e-39_r8, 4.401531211660759e-33_r8, 3.0750274997756505e-27_r8, 1.7712963947651545e-21_r8, 8.11760687085515e-16_r8, 2.8742026224385195e-10_r8, 5.884546319839207e-05_r8, 6.38887267816177_r8, 6.3443642690380315_r8, 4.000750061915931_r8, 1.1986831699513256_r8, 0.596798703981576_r8, 2.386571507304945_r8, 4.2390870207605005_r8, 4.3729102638376665_r8, 4.021306923395206_r8, 3.5339665817020194_r8, 3.3988040130833124_r8, 5.2962672024737465_r8, 4.946001742972865_r8, 4.674512907937394_r8, 4.034058594187169_r8, 3.5591567560448185_r8, 3.2117619336921175_r8, 2.8488663218424306_r8, 2.3450054899294184_r8, 1.821989658621705_r8, 1.418813430289268_r8 /

  ! 1-column wrappers (the operator's innermost index is the column).
  d_coef(1, :) = dpidz_sq
  u(1, :) = u0

  coords = Coords1D(reshape(pint, [1, pver + 1]))
  op = diffusion_operator(coords, d_coef)

  out = op%apply(u)

  write (*, '(a)') '# CAM-order A*u0 rows (lev, m/s^2):'
  do k = 1, pver
    write (*, '(i3, 1x, es26.18)') k, out(1, k)
  end do

  call op%finalize()
  call coords%finalize()

end program vdiff_driver
