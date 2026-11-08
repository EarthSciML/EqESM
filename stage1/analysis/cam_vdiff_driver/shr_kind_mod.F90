! Minimal shr_kind_mod stub for the standalone vdiff driver build.
! Faithful to share/src/shr_kind_mod.F90 of the fork (the r8 kind CAM uses).
module shr_kind_mod

  integer, parameter :: shr_kind_r8 = selected_real_kind(12)  ! 8 byte real

end module shr_kind_mod
