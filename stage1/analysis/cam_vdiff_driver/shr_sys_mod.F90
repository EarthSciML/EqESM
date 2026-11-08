! Minimal shr_sys_mod stub for the standalone vdiff driver build.
module shr_sys_mod

contains

  subroutine shr_sys_abort(msg)
    character(len=*), optional, intent(in) :: msg
    if (present(msg)) then
      write (*, '(a)') 'ABORT: ', trim(msg)
    end if
    error stop
  end subroutine shr_sys_abort

end module shr_sys_mod
