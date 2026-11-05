! Minimal shr_log_mod stub for the standalone vdiff driver build.
! Only errMsg is used, in "..."//errMsg(__FILE__, __LINE__) string
! concatenations inside linear_1d_operators.F90.
module shr_log_mod

  character(len=512), parameter :: shr_kind_cx = "512"

contains

  function shr_log_errMsg(file, line) result(msg)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    character(len=512) :: msg
    write (msg, '(a,a,i0)') 'ERROR in ', trim(file), ' at line ', line
  end function shr_log_errMsg

end module shr_log_mod
