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

    function grhogal(point) bind(C, name='grhogal_')
      import C_DOUBLE
      real(kind=C_DOUBLE) :: point, grhogal
    end function grhogal

    function gpresgal(point) bind(C, name='gpresgal_')
      import C_DOUBLE
      real(kind=C_DOUBLE) :: point, gpresgal
      end function gpresgal

    function Chigal(dgrho, eta, dphi, dphiprime, point, k) bind(C, name='Chigal_')
      import C_DOUBLE
      real(kind=C_DOUBLE) :: dgrho, eta, dphi, dphiprime, point, k, Chigal
    end function Chigal

    function qgal(dgq, eta, dphi, dphiprime, point, k) bind(C, name='qgal_')
      import C_DOUBLE
      real(kind=C_DOUBLE) :: dgq, eta, dphi, dphiprime, point, k, qgal
    end function qgal

    function Pigal(dgrho, dgq, dgpi, eta, dphi, point, k) bind(C, name='Pigal_')
      import C_DOUBLE
      real(kind=C_DOUBLE) :: dgrho, dgq, dgpi, eta, dphi, point, k, Pigal
    end function Pigal

    function dphisecond(dgrho, dgq, z, eta, dphi, dphiprime, point, k, grho, gpres, grhob, clxb, clxbdot, grhoc, clxc, clxcdot, grhor, clxr, clxrdot, grhog, clxg, clxgdot) bind(C, name='dphisecond_')
      import C_DOUBLE
      real(kind=C_DOUBLE) :: dgrho, dgq, z, eta, dphi, dphiprime, point, k, dphisecond, grho, gpres, grhob, clxb, clxbdot, grhoc, clxc, clxcdot, grhor, clxr, clxrdot, grhog, clxg, clxgdot
    end function dphisecond
!!$
!!$    function dphisecond(grho, gpres, grhob, clxb, clxbdot, grhoc, clxc, clxcdot, &
!!$         grhor, clxr, clxrdot, grhog, clxg, clxgdot, dgrho, dgq, z, &
!!$         eta, dphi, dphiprime, point, k) bind(C, name='dphisecond_')
!!$      import C_DOUBLE
!!$      real(kind=C_DOUBLE) :: grho, gpres
!!$      real(kind=C_DOUBLE) :: grhob, clxb, clxbdot
!!$      real(kind=C_DOUBLE) :: grhoc, clxc, clxcdot
!!$      real(kind=C_DOUBLE) :: grhor, clxr, clxrdot
!!$      real(kind=C_DOUBLE) :: grhog, clxg, clxgdot
!!$      real(kind=C_DOUBLE) :: dgrho, dgq, z, eta, dphi, dphiprime, point, k
!!$    end function dphisecond

    subroutine zprime(grho, dgrho, dgq, grhob, clxb, clxbdot, grhoc, clxc, clxcdot, grhor, clxr, clxrdot, grhog, clxg, clxgdot, eta, dphi, dphiprime, dphisecond, point, k) bind(C, name='zprime_')
      import C_DOUBLE
      real(kind=C_DOUBLE) :: grho, dgrho, dgq, grhob, clxb, clxbdot, grhoc, clxc, clxcdot, grhor, clxr, clxrdot, grhog, clxg, clxgdot, eta, dphi, dphiprime, dphisecond, point, k
    end subroutine zprime

    type(C_PTR) function conservation(grho, gpres, dgrho, grhob, clxb, clxbdot, grhoc, clxc, clxcdot, grhor, clxr, clxrdot, grhog, clxg, clxgdot, dgq, qr, qrdot, qg, qgdot, dgpi, eta, dphi, dphiprime, dphisecond, point, k) bind(C, name='conservation_')
      import C_PTR, C_DOUBLE
      real(kind=C_DOUBLE) :: grho, gpres, dgrho, grhob, clxb, clxbdot, grhoc, clxc, clxcdot, grhor, clxr, clxrdot, grhog, clxg, clxgdot, dgq, qr, qrdot, qg, qgdot, dgpi, eta, dphi, dphiprime, dphisecond, point, k
    end function conservation

  end interface

end module interface
