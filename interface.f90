module interface
  !use, intrinsic :: iso_c_binding
  use ISO_C_BINDING
  implicit none

  interface
!!$     subroutine arrays(infile, outfile, omegar) bind(C, name='arrays_')
     function arrays(infile, omegar, omegam, H0in, c2in, c3in, c4in, c5in, cGin, grhormass, nu_masses, nu_mass_eigenstates) bind(C, name='arrays_')
       import C_CHAR, C_DOUBLE, C_INT, C_PTR ! Make iso c binding visible here
!!$       character(kind=C_CHAR), dimension(*) :: infile, outfile
       character(kind=C_CHAR), dimension(*) :: infile
       real(kind=C_DOUBLE) :: omegar, omegam, H0in, c2in, c3in, c4in, c5in, cGin
       integer(C_INT) :: nu_mass_eigenstates
       type(C_PTR), value :: grhormass, nu_masses
!!$       real(C_DOUBLE), intent(in), dimension(*) :: grhormass, nu_masses
       integer :: arrays
    end function arrays

    function handxofa(point) bind(C, name='handxofa_')
      import C_PTR, C_DOUBLE
      real(kind=C_DOUBLE) :: point
      type(C_PTR) :: handxofa
    end function handxofa

    function grhogal(point) bind(C, name='grhogal_')
      import C_DOUBLE
      real(kind=C_DOUBLE) :: point, grhogal
    end function grhogal

    function gpresgal() bind(C, name='gpresgal_')
      import C_DOUBLE
      real(kind=C_DOUBLE) :: gpresgal
      end function gpresgal

    function Chigal(dgrho, eta, dphi, dphiprime, k) bind(C, name='Chigal_')
      import C_DOUBLE
      real(kind=C_DOUBLE) :: dgrho, eta, dphi, dphiprime, k, Chigal
    end function Chigal

    function qgal(dgq, eta, dphi, dphiprime, k) bind(C, name='qgal_')
      import C_DOUBLE
      real(kind=C_DOUBLE) :: dgq, eta, dphi, dphiprime, k, qgal
    end function qgal

    function Pigal(dgrho, dgq, dgpi, eta, dphi, k) bind(C, name='Pigal_')
      import C_DOUBLE
      real(kind=C_DOUBLE) :: dgrho, dgq, dgpi, eta, dphi, k, Pigal
    end function Pigal

    function dphisecond(dgrho, dgq, eta, dphi, dphiprime, k, deltafprime) bind(C, name='dphisecond_')
      import C_DOUBLE
      real(kind=C_DOUBLE) :: dgrho, dgq, eta, dphi, dphiprime, k, dphisecond, deltafprime
    end function dphisecond

    function pigalprime(dgrho, dgq, dgpi, pidot, eta, dphi, dphiprime, k, grho, gpres) bind(C, name='pigalprime_')
      import C_DOUBLE
      real(kind=C_DOUBLE) :: dgrho, dgq, dgpi, pidot, eta, dphi, dphiprime, k, grho, gpres, pigalprime
    end function pigalprime

    function crosschecks(dgrho, dgq, dgpi, eta, dphi, dphiprime, dphiprimeprime, k, grho, gpres, deltafprime) bind(C, name='crosschecks_')
      import C_PTR, C_DOUBLE
      real(kind=C_DOUBLE) :: dgrho, dgq, dgpi, eta, dphi, dphiprime, dphiprimeprime, k, grho, gpres, deltafprime
      type(C_PTR) :: crosschecks
    end function crosschecks

    subroutine freegal() bind(C, name='freegal_')
      !subroutine to free non-necessary memory after perturbations calculation
    end subroutine freegal
  end interface
end module interface
