module interface
  !use, intrinsic :: iso_c_binding
  use ISO_C_BINDING
  !implicit none

  interface
     !subroutine arrays(infile, outfile, omegar) bind(C, name='arrays_')
     subroutine arrays(infile, omegar) bind(C, name='arrays_')
       import C_CHAR, C_DOUBLE ! Make iso c binding visible here
!!$       character(kind=C_CHAR), dimension(*) :: infile, outfile
       character(kind=C_CHAR), dimension(*) :: infile
       real(kind=C_DOUBLE) :: omegar
    end subroutine arrays

    type(C_PTR) function handxofa(point) bind(C, name='handxofa_')
      import C_PTR, C_DOUBLE
      real(kind=C_DOUBLE) :: point
      end function handxofa
  end interface

end module interface

!!$PROGRAM appel_C
!!$  use background_C
!!$  IMPLICIT NONE
!!$
!!$  CHARACTER(len=40) :: infile = 'galileon.ini' //C_NULL_CHAR, outfile = 'background.txt' // C_NULL_CHAR
!!$  type(C_PTR) :: cptr_to_handx
!!$  real(kind=C_DOUBLE), pointer :: handx(:)
!!$  logical :: ex
!!$
!!$  inquire(file=infile, exist=ex)
!!$  IF (ex) THEN
!!$     PRINT *, "File exists"
!!$  END IF
!!$
!!$  call arrays(infile, outfile)
!!$  call C_F_POINTER(cptr_to_array, h_array, [27639])
!!$
!!$  cptr_to_handx = handxofa(1.0007491242d-6)
!!$  call C_F_POINTER(cptr_to_handx, handx, [2])
!!$
!!$  PRINT *, handx
!!$
!!$END PROGRAM appel_C