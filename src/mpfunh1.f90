!*****************************************************************************

!  MPFUN20-Fort: A thread-safe arbitrary precision package with special functions
!  Medium precision language interface module (module MPFUNH)
!  Variant 1: Precision level arguments are NOT required.

!  Note: !> and !>> comments delimit variant differences, and are used to generate
!  the variant files of this module using the gencodes.f90 program. Do not
!  change these comments.

!  Revision date:  13 May 2023

!  AUTHOR:
!    David H. Bailey
!    Lawrence Berkeley National Lab (retired) and University of California, Davis
!    Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2023 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!  PURPOSE OF PACKAGE:
!    This package permits one to perform floating-point computations (real and
!    complex) to arbitrarily high numeric precision, by making only relatively
!    minor changes to existing Fortran-90 programs. All basic arithmetic
!    operations and transcendental functions are supported, together with numerous
!    special functions.

!    In addition to fast execution times, one key feature of this package is a
!    100% THREAD-SAFE design, which means that user-level applications can be
!    easily converted for parallel execution, say using a threaded parallel
!    environment such as OpenMP.

!  DOCUMENTATION:
!    A detailed description of this package, and instructions for compiling
!    and testing this program on various specific systems are included in the
!    README file accompanying this package, and, in more detail, in the
!    following technical paper:

!    David H. Bailey, "MPFUN2020: A new thread-safe arbitrary precision package,"
!    available at http://www.davidhbailey.com/dhbpapers/mpfun2020.pdf.

!  DESCRIPTION OF THIS MODULE (MPFUNH):
!    This module contains all high-level Fortran-90 language interfaces for
!    operations using the medium precision datatype. This is not needed for
!    most programs, but provides significantly better performance (and
!    less memory) for certain large, demanding applications that can benefit
!    from two high-precision datatypes (e.g., tpslqm3.f90 and tpphix3.f90).

!    There are two variants of this module:
!      mpfunh1.f90	Mixed-mode operations ARE allowed; precision level
!                     arguments are NOT required in certain functions;
!                     quad precision (128-bit) is NOT supported.
!      mpfunh2.f90	Mixed-mode operations are NOT allowed; precision level
!                     arguments ARE required in certain functions.
!                     quad precision (128-bit) is NOT supported.
!    Compile and link whichever one of these is most appropriate for the
!    given applications. See documentation for details.

module mpfunh
use mpfuna
use mpfunb
use mpfunc
use mpfund
use mpfune
use mpfunf
use mpfung
implicit none

!  Assignments and the five arithmetic operators:

private &
  mp_eqrr, mp_eqdr, mp_eqrd, mp_eqir, mp_eqri, mp_eqra, mp_eqrz, &
  mp_eqzr, mp_eqzz, mp_eqdz, mp_eqzd, mp_eqdcz, mp_eqzdc, &
  mp_addrr, mp_adddr, mp_addrd, mp_addir, mp_addri, mp_addzz, &
  mp_adddz, mp_addzd, mp_adddcz, mp_addzdc, mp_addrz, mp_addzr, &
  mp_subrr, mp_subdr, mp_subrd, mp_subir, mp_subri, mp_subzz, &
  mp_subdz, mp_subzd, mp_subdcz, mp_subzdc, mp_subrz, mp_subzr, &
  mp_negr, mp_negz, &
  mp_mulrr, mp_muldr, mp_mulrd, mp_mulir, mp_mulri, mp_mulzz, &
  mp_muldz, mp_mulzd, mp_muldcz, mp_mulzdc, mp_mulrz, mp_mulzr, &
  mp_divrr, mp_divdr, mp_divrd, mp_divir, mp_divri, mp_divzz, &
  mp_divdz, mp_divzd, mp_divdcz, mp_divzdc, mp_divrz, mp_divzr, &
  mp_expri, mp_exprr, mp_expzi, mp_expzz, mp_exprz, mp_expzr

!  The six comparison tests:

private &
  mp_eqtrr, mp_eqtdr, mp_eqtrd, mp_eqtir, mp_eqtri, mp_eqtzz, &
  mp_eqtdz, mp_eqtzd, mp_eqtdcz, mp_eqtzdc, mp_eqtrz, mp_eqtzr, &
  mp_netrr, mp_netdr, mp_netrd, mp_netir, mp_netri, mp_netzz, &
  mp_netdz, mp_netzd, mp_netdcz, mp_netzdc, mp_netrz, mp_netzr, &
  mp_letrr, mp_letdr, mp_letrd, mp_letir, mp_letri, &
  mp_getrr, mp_getdr, mp_getrd, mp_getir, mp_getri, &
  mp_lttrr, mp_lttdr, mp_lttrd, mp_lttir, mp_lttri, &
  mp_gttrr, mp_gttdr, mp_gttrd, mp_gttir, mp_gttri

!  Algebraic, transcendental and type conversion functions:

private &
  mp_absr, mp_absz, mp_acos, mp_acosh, mp_agm, mp_aimag, mp_aint, &
  mp_anint, mp_asin, mp_asinh, mp_atan, mp_atan2, mp_atanh, mp_ator1, &
  mp_atorn, mp_berne, mp_bessel_i, mp_bessel_in, mp_bessel_j, &
  mp_bessel_j0, mp_bessel_j1, mp_bessel_jn, mp_bessel_k, mp_bessel_kn, &
  mp_bessel_y, mp_bessel_yn, mp_bessel_y0, mp_bessel_y1, mp_binmd, &
  mp_ccos, mp_cexp, mp_clog, mp_conjg, mp_cos, mp_cosh, mp_csin, &
  mp_csqrt, mp_cssh, mp_cssn, mp_dctoz, mp_dctoz2, mp_decmd, &
  mp_digamma_be, mp_dtor, mp_dtor2, mp_eform, mp_egamma, mp_erf, &
  mp_erfc, mp_exp, mp_expint, mp_fform, mp_gamma, mp_hurwitz_zetan, &
  mp_hurwitz_zetan_be, mp_hypergeom_pfq, mp_hypot, mp_incgamma, &
  mp_log, mp_log10, mp_log2, mp_max, mp_min, mp_mod, mp_nrt, mp_pi, &
  mp_polygamma, mp_polygamma_be, mp_polylog_inix, mp_polylog_negx, &
  mp_polylog_pos, mp_prodd, mp_prodq, mp_qtor, mp_qtor2, mp_quotd, &
  mp_quotq, mp_rand, mp_readr1, mp_readr2, mp_readr3, mp_readr4, &
  mp_readr5, mp_readz1, mp_readz2, mp_readz3, mp_readz4, mp_readz5, &
  mp_rtod, mp_rtom, mp_rtoq, mp_rtor, mp_rtoz, mp_setwp, mp_sign, &
  mp_sin, mp_sinh, mp_sqrt, mp_struve_hn, mp_tan, mp_tanh, mp_wprec, &
  mp_wprecz, mp_writer, mp_writez, mp_zeta, mp_zeta_be, mp_zeta_int, &
  mp_ztodc, mp_ztor, mp_ztoz, mp_ztozm

!  Compared with the corresponding statement in mpfung2.f90, the above list
!  omits mp_init, mp_mtor and mp_zmtoz, but includes mp_rtom and mp_ztozm.

!  Operator extension interface blocks:

interface assignment (=)
  module procedure mp_eqrr
  module procedure mp_eqdr
  module procedure mp_eqir
  module procedure mp_eqrz
  module procedure mp_eqzr
  module procedure mp_eqzz
  module procedure mp_eqdz
  module procedure mp_eqdcz

!>  In variant #1, uncomment these lines:
  module procedure mp_eqrd
  module procedure mp_eqri
  module procedure mp_eqra
  module procedure mp_eqzd
  module procedure mp_eqzdc
!>>
end interface

interface operator (+)
  module procedure mp_addrr
  module procedure mp_adddr
  module procedure mp_addrd
  module procedure mp_addir
  module procedure mp_addri
  module procedure mp_addzz
  module procedure mp_adddz
  module procedure mp_addzd
  module procedure mp_adddcz
  module procedure mp_addzdc
  module procedure mp_addrz
  module procedure mp_addzr
end interface

interface operator (-)
  module procedure mp_subrr
  module procedure mp_subdr
  module procedure mp_subrd
  module procedure mp_subir
  module procedure mp_subri
  module procedure mp_subzz
  module procedure mp_subdz
  module procedure mp_subzd
  module procedure mp_subdcz
  module procedure mp_subzdc
  module procedure mp_subrz
  module procedure mp_subzr
  module procedure mp_negr
  module procedure mp_negz
end interface

interface operator (*)
  module procedure mp_mulrr
  module procedure mp_muldr
  module procedure mp_mulrd
  module procedure mp_mulir
  module procedure mp_mulri
  module procedure mp_mulzz
  module procedure mp_muldz
  module procedure mp_mulzd
  module procedure mp_muldcz
  module procedure mp_mulzdc
  module procedure mp_mulrz
  module procedure mp_mulzr
end interface

interface operator (/)
  module procedure mp_divrr
  module procedure mp_divdr
  module procedure mp_divrd
  module procedure mp_divir
  module procedure mp_divri
  module procedure mp_divzz
  module procedure mp_divdz
  module procedure mp_divzd
  module procedure mp_divdcz
  module procedure mp_divzdc
  module procedure mp_divrz
  module procedure mp_divzr
end interface

interface operator (**)
   module procedure mp_expri
   module procedure mp_exprr
   module procedure mp_expzi
   module procedure mp_expzz
   module procedure mp_exprz
   module procedure mp_expzr
end interface

interface operator (==)
  module procedure mp_eqtrr
  module procedure mp_eqtdr
  module procedure mp_eqtrd
  module procedure mp_eqtir
  module procedure mp_eqtri
  module procedure mp_eqtzz
  module procedure mp_eqtdz
  module procedure mp_eqtzd
  module procedure mp_eqtdcz
  module procedure mp_eqtzdc
  module procedure mp_eqtrz
  module procedure mp_eqtzr
end interface

interface operator (/=)
  module procedure mp_netrr
  module procedure mp_netdr
  module procedure mp_netrd
  module procedure mp_netir
  module procedure mp_netri
  module procedure mp_netzz
  module procedure mp_netdz
  module procedure mp_netzd
  module procedure mp_netdcz
  module procedure mp_netzdc
  module procedure mp_netrz
  module procedure mp_netzr
end interface

interface operator (<=)
  module procedure mp_letrr
  module procedure mp_letdr
  module procedure mp_letrd
  module procedure mp_letir
  module procedure mp_letri
end interface

interface operator (>=)
  module procedure mp_getrr
  module procedure mp_getdr
  module procedure mp_getrd
  module procedure mp_getir
  module procedure mp_getri
end interface

interface operator (<)
  module procedure mp_lttrr
  module procedure mp_lttdr
  module procedure mp_lttrd
  module procedure mp_lttir
  module procedure mp_lttri
end interface

interface operator (>)
  module procedure mp_gttrr
  module procedure mp_gttdr
  module procedure mp_gttrd
  module procedure mp_gttir
  module procedure mp_gttri
end interface

!  MP generic function interface blocks, listed alphabetically by interface name:

interface abs
  module procedure mp_absr
  module procedure mp_absz
end interface

interface acos
  module procedure mp_acos
end interface

interface acosh
  module procedure mp_acosh
end interface

interface agm
  module procedure mp_agm
end interface

interface aimag
  module procedure mp_aimag
end interface

interface aint
  module procedure mp_aint
end interface

interface anint
  module procedure mp_anint
end interface

interface asin
  module procedure mp_asin
end interface

interface asinh
  module procedure mp_asinh
end interface

interface atan
  module procedure mp_atan
end interface

interface atan2
  module procedure mp_atan2
end interface

interface atanh
  module procedure mp_atanh
end interface

interface mpberne
  module procedure mp_berne
end interface

interface bessel_i
  module procedure mp_bessel_i
end interface

interface bessel_in
  module procedure mp_bessel_in
end interface

interface bessel_j
  module procedure mp_bessel_j
end interface

interface bessel_jn
  module procedure mp_bessel_jn
end interface

interface bessel_j0
  module procedure mp_bessel_j0
end interface

interface bessel_j1
  module procedure mp_bessel_j1
end interface

interface bessel_k
  module procedure mp_bessel_k
end interface

interface bessel_kn
  module procedure mp_bessel_kn
end interface

interface bessel_y
  module procedure mp_bessel_y
end interface

interface bessel_yn
  module procedure mp_bessel_yn
end interface

interface bessel_y0
  module procedure mp_bessel_y0
end interface

interface bessel_y1
  module procedure mp_bessel_y1
end interface

interface conjg
  module procedure mp_conjg
end interface

interface cos
  module procedure mp_cos
  module procedure mp_ccos
end interface

interface cosh
  module procedure mp_cosh
end interface

interface dble
  module procedure mp_rtod
end interface

interface dcmplx
  module procedure mp_ztodc
end interface

interface digamma_be
  module procedure mp_digamma_be
end interface

interface erf
  module procedure mp_erf
end interface

interface erfc
  module procedure mp_erfc
end interface

interface exp
  module procedure mp_exp
  module procedure mp_cexp
end interface

interface expint
  module procedure mp_expint
end interface

interface gamma
  module procedure mp_gamma
end interface

interface hurwitz_zetan
  module procedure mp_hurwitz_zetan
end interface

interface hurwitz_zetan_be
  module procedure mp_hurwitz_zetan_be
end interface

interface hypergeom_pfq
  module procedure mp_hypergeom_pfq
end interface

interface hypot
  module procedure mp_hypot
end interface

interface log
  module procedure mp_log
  module procedure mp_clog
end interface

interface log10
  module procedure mp_log10
end interface

interface max
  module procedure mp_max
end interface

interface min
  module procedure mp_min
end interface

interface mod
  module procedure mp_mod
end interface

interface mpbinmd
  module procedure mp_binmd
end interface

interface mpcmplxm
  module procedure mp_dctoz
  module procedure mp_rtoz
  module procedure mp_ztoz
  module procedure mp_ztozm
end interface

interface mpcmplxdcm
  module procedure mp_dctoz2
end interface

interface mpcssh
  module procedure mp_cssh
end interface

interface mpcssn
  module procedure mp_cssn
end interface

interface mpdecmd
  module procedure mp_decmd
end interface

interface mpeform
  module procedure mp_eform
end interface

interface mpegammam
  module procedure mp_egamma
end interface

interface mpfform
  module procedure mp_fform
end interface

interface incgamma
  module procedure mp_incgamma
end interface

interface mplog2m
  module procedure mp_log2
end interface

interface mpnrt
  module procedure mp_nrt
end interface

interface mppim
  module procedure mp_pi
end interface

interface mpprod
  module procedure mp_prodd
  module procedure mp_prodq
end interface

interface mpquot
  module procedure mp_quotd
  module procedure mp_quotq
end interface

interface mprand
  module procedure mp_rand
end interface

interface mpread
  module procedure mp_readr1
  module procedure mp_readr2
  module procedure mp_readr3
  module procedure mp_readr4
  module procedure mp_readr5
  module procedure mp_readz1
  module procedure mp_readz2
  module procedure mp_readz3
  module procedure mp_readz4
  module procedure mp_readz5
end interface

interface mprealm
  module procedure mp_ator1
  module procedure mp_atorn
  module procedure mp_dtor
  module procedure mp_rtor
  module procedure mp_ztor
  module procedure mp_rtom
  module procedure mp_qtor
end interface

interface mprealdm
  module procedure mp_dtor2
end interface

interface mprealqm
  module procedure mp_qtor2
end interface

interface mpwprec
  module procedure mp_wprec
  module procedure mp_wprecz
end interface

interface mpwrite
  module procedure mp_writer
  module procedure mp_writez
end interface

interface polygamma
  module procedure mp_polygamma
end interface

interface polygamma_be
  module procedure mp_polygamma_be
end interface

interface polylog_ini
  module procedure mp_polylog_inix
end interface

interface polylog_neg
  module procedure mp_polylog_negx
end interface

interface polylog_pos
  module procedure mp_polylog_pos
end interface

interface qreal
  module procedure mp_rtoq
end interface

interface sign
  module procedure mp_sign
end interface

interface sin
  module procedure mp_sin
  module procedure mp_csin
end interface

interface sinh
  module procedure mp_sinh
end interface

interface sqrt
  module procedure mp_sqrt
  module procedure mp_csqrt
end interface

interface struve_hn
  module procedure mp_struve_hn
end interface

interface tan
  module procedure mp_tan
end interface

interface tanh
  module procedure mp_tanh
end interface

interface zeta
  module procedure mp_zeta
end interface

interface zeta_be
  module procedure mp_zeta_be
end interface

interface zeta_intm
  module procedure mp_zeta_int
end interface

contains

!  This routine outputs an error message if iprec exceeds mpwdsm.

  function mp_setwp (iprec)
    integer mp_setwp
    integer, intent (in):: iprec
    if (iprec > mpwdsm) then
      write (mpldb, 1)
1       format ( &
        '*** MP_SETWP: requested precision level exceeds default medium precision.'/ &
        'Increase default medium precision in module MPFUNF.')
      call mpabrt ( 801)
    endif
    mp_setwp = iprec
  end function

!  Assignment routines:

  subroutine mp_eqrr (ra, rb)
    implicit none
    type (mp_realm), intent (out):: ra
    type (mp_realm), intent (in):: rb
    integer mpnw
    mpnw = min (int (rb%mpr(1)), mpwdsm)
    ra%mpr(0) = mpwdsm6
    call mpeq (rb%mpr, ra%mpr, mpnw)
    return
  end subroutine

  subroutine mp_eqdr (da, rb)
    implicit none
    real (mprknd), intent (out):: da
    type (mp_realm), intent (in):: rb
    integer ib, mpnw
    mpnw = min (int (rb%mpr(1)), mpwdsm)
    call mpmdc (rb%mpr, da, ib, mpnw)
    da = da * 2.d0 ** ib
    return
  end subroutine

  subroutine mp_eqrd (ra, db)
    implicit none
    type (mp_realm), intent (out):: ra
    real (mprknd), intent (in):: db
    integer i1, mpnw
    mpnw = mpwdsm
    ra%mpr(0) = mpwdsm6
    i1 = 0
    call mpdmc40 (db, i1, ra%mpr, mpnw)
    return
  end subroutine

  subroutine mp_eqir (ia, rb)
    implicit none
    integer, intent (out):: ia
    type (mp_realm), intent (in):: rb
    real (mprknd) da
    integer ib, mpnw
    mpnw = min (int (rb%mpr(1)), mpwdsm)
    call mpmdc (rb%mpr, da, ib, mpnw)
    ia = da * 2.d0 ** ib
    return
  end subroutine

  subroutine mp_eqri (ra, ib)
    implicit none
    type (mp_realm), intent (out):: ra
    integer, intent (in):: ib
    real (mprknd) db
    integer i1, mpnw
    mpnw = mpwdsm
    ra%mpr(0) = mpwdsm6
    db = ib
    i1 = 0
    call mpdmc40 (db, i1, ra%mpr, mpnw)
    return
  end subroutine

  subroutine mp_eqra (ra, ab)
    implicit none
    type (mp_realm), intent (out):: ra
    character(*), intent (in):: ab
    character(1) :: chr1(len(ab))
    integer i, l1, mpnw
    mpnw = mpwdsm
    l1 = len (ab)
    do i = 1, l1
      chr1(i) = ab(i:i)
    enddo
    ra%mpr(0) = mpwdsm6
    call mpctomp (chr1, l1, ra%mpr, mpnw)
    return
  end subroutine

  subroutine mp_eqzz (za, zb)
    implicit none
    type (mp_complexm), intent (out):: za
    type (mp_complexm), intent (in):: zb
    integer l1, l2, mpnw
    l1 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    za%mpc(0) = mpwdsm6
    za%mpc(l2) = mpwdsm6
    call mpceq (zb%mpc, za%mpc, mpnw)
    return
  end subroutine

  subroutine mp_eqdz (da, zb)
    implicit none
    real (mprknd), intent (out):: da
    type (mp_complexm), intent (in):: zb
    integer l1, mpnw, n1
    real (mprknd) d1
    l1 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    call mpmdc (zb%mpc, d1, n1, mpnw)
    da = d1 * 2.d0 ** n1
    return
  end subroutine

  subroutine mp_eqzd (za, db)
    implicit none
    type (mp_complexm), intent (out):: za
    real (mprknd), intent (in):: db
    real (mprknd) d1
    integer l1, mpnw
    mpnw = mpwdsm
    l1 = mpwdsm6
    za%mpc(0) = mpwdsm6
    za%mpc(l1) = mpwdsm6
    d1 = 0.d0
    call mpdmc40 (db, 0, za%mpc, mpnw)
    call mpdmc40 (d1, 0, za%mpc(l1:), mpnw)
    return
  end subroutine

  subroutine mp_eqdcz (dca, zb)
    implicit none
    complex (kind (0.d0)), intent (out):: dca
    type (mp_complexm), intent (in):: zb
    integer l1, mpnw, n1, n2
    real (mprknd) d1, d2
    l1 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    call mpmdc (zb%mpc, d1, n1, mpnw)
    d1 = d1 * 2.d0 ** n1
    call mpmdc (zb%mpc(l1:), d2, n2, mpnw)
    d2 = d2 * 2.d0 ** n2
    dca = cmplx (d1, d2, mprknd)
    return
  end subroutine

  subroutine mp_eqzdc (za, dcb)
    implicit none
    type (mp_complexm), intent (out):: za
    complex (kind (0.d0)), intent (in):: dcb
    integer l1, mpnw
    mpnw = mpwdsm
    l1 = mpwdsm6
    za%mpc(0) = mpwdsm6
    za%mpc(l1) = mpwdsm6
    call mpdmc40 (dble (dcb), 0, za%mpc, mpnw)
    call mpdmc40 (aimag (dcb), 0, za%mpc(l1:), mpnw)
    return
  end subroutine

  subroutine mp_eqrz (ra, zb)
    implicit none
    type (mp_realm), intent (out):: ra
    type (mp_complexm), intent (in):: zb
    integer l1, mpnw
    l1 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    ra%mpr(0) = mpwdsm6
    call mpeq (zb%mpc, ra%mpr, mpnw)
    return
  end subroutine

  subroutine mp_eqzr (za, rb)
    implicit none
    type (mp_complexm), intent (out):: za
    type (mp_realm), intent (in):: rb
    type (mp_realm) r1
    integer l1, mpnw
    mpnw = min (int (rb%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    call mpdmc (0.d0, 0, r1%mpr, mpnw)
    l1 = mpwdsm6
    za%mpc(0) = mpwdsm6
    za%mpc(l1) = mpwdsm6
    call mpeq (rb%mpr, za%mpc, mpnw)
    call mpeq (r1%mpr, za%mpc(l1:), mpnw)
    return
  end subroutine

!  Addition routines:

  function mp_addrr (ra, rb)
    implicit none
    type (mp_realm):: mp_addrr
    type (mp_realm), intent (in):: ra, rb
    integer mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwdsm)
    mp_addrr%mpr(0) = mpwdsm6
    call mpadd (ra%mpr, rb%mpr, mp_addrr%mpr, mpnw)
    return
  end function

  function mp_adddr (da, rb)
    implicit none
    type (mp_realm):: mp_adddr
    real (mprknd), intent (in):: da
    type (mp_realm), intent (in):: rb
    type (mp_realm) r1
    integer i1, mpnw
    mpnw = min (int (rb%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    mp_adddr%mpr(0) = mpwdsm6
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpadd (r1%mpr, rb%mpr, mp_adddr%mpr, mpnw)
    return
  end function

  function mp_addrd (ra, db)
    implicit none
    type (mp_realm):: mp_addrd
    type (mp_realm), intent (in):: ra
    real (mprknd), intent (in):: db
    type (mp_realm) r1
    integer i1, mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    mp_addrd%mpr(0) = mpwdsm6
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpadd (ra%mpr, r1%mpr, mp_addrd%mpr, mpnw)
    return
  end function

  function mp_addir (ia, rb)
    implicit none
    type (mp_realm):: mp_addir
    integer, intent (in):: ia
    type (mp_realm), intent (in):: rb
    type (mp_realm) r1
    real (mprknd) da
    integer i1, mpnw
    mpnw = min (int (rb%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    mp_addir%mpr(0) = mpwdsm6
    da = ia
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpadd (r1%mpr, rb%mpr, mp_addir%mpr, mpnw)
    return
  end function

  function mp_addri (ra, ib)
    implicit none
    type (mp_realm):: mp_addri
    type (mp_realm), intent (in):: ra
    integer, intent (in):: ib
    real (mprknd) db
    type (mp_realm) r1
    integer i1, mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    mp_addri%mpr(0) = mpwdsm6
    db = ib
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpadd (ra%mpr, r1%mpr, mp_addri%mpr, mpnw)
    return
  end function

  function mp_addzz (za, zb)
    implicit none
    type (mp_complexm):: mp_addzz
    type (mp_complexm), intent (in):: za, zb
    integer l1, l2, l3, mpnw
    l1 = za%mpc(0)
    l2 = zb%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (zb%mpc(1)), &
      int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwdsm)
    l3 = mpwdsm6
    mp_addzz%mpc(0) = mpwdsm6
    mp_addzz%mpc(l3) = mpwdsm6
    call mpcadd (za%mpc, zb%mpc, mp_addzz%mpc, mpnw)
    return
  end function

  function mp_adddz (da, zb)
    implicit none
    type (mp_complexm):: mp_adddz
    real (mprknd) , intent (in):: da
    type (mp_complexm), intent (in):: zb
    integer l1, l2, mpnw
    real (mprknd) d1
    type (mp_complexm) z1
    l1 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    z1%mpc(0) = mpwdsm6
    z1%mpc(l2) = mpwdsm6
    mp_adddz%mpc(0) = mpwdsm6
    mp_adddz%mpc(l2) = mpwdsm6
    d1 = 0.d0
    call mpdmc40 (da, 0, z1%mpc, mpnw)
    call mpdmc40 (d1, 0, z1%mpc(l2:), mpnw)
    call mpcadd (z1%mpc, zb%mpc, mp_adddz%mpc, mpnw)
    return
  end function

  function mp_addzd (za, db)
    implicit none
    type (mp_complexm):: mp_addzd
    type (mp_complexm), intent (in):: za
    real (mprknd), intent (in):: db
    integer l1, l2, mpnw
    real (mprknd) d1
    type (mp_complexm) z1
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    z1%mpc(0) = mpwdsm6
    z1%mpc(l2) = mpwdsm6
    mp_addzd%mpc(0) = mpwdsm6
    mp_addzd%mpc(l2) = mpwdsm6
    d1 = 0.d0
    call mpdmc40 (db, 0, z1%mpc, mpnw)
    call mpdmc40 (d1, 0, z1%mpc(l2:), mpnw)
    call mpcadd (za%mpc, z1%mpc, mp_addzd%mpc, mpnw)
    return
  end function

  function mp_adddcz (dca, zb)
    implicit none
    type (mp_complexm):: mp_adddcz
    complex (kind (0.d0)), intent (in):: dca
    type (mp_complexm), intent (in):: zb
    integer l1, l2, mpnw
    type (mp_complexm) z1
    l1 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    z1%mpc(0) = mpwdsm6
    z1%mpc(l2) = mpwdsm6
    mp_adddcz%mpc(0) = mpwdsm6
    mp_adddcz%mpc(l2) = mpwdsm6
    call mpdmc40 (dble (dca), 0, z1%mpc, mpnw)
    call mpdmc40 (aimag (dca), 0, z1%mpc(l2:), mpnw)
    call mpcadd (z1%mpc, zb%mpc, mp_adddcz%mpc, mpnw)
    return
  end function

  function mp_addzdc (za, dcb)
    implicit none
    type (mp_complexm):: mp_addzdc
    type (mp_complexm), intent (in):: za
    complex (kind (0.d0)), intent (in):: dcb
    integer l1, l2, mpnw
    type (mp_complexm) z1
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    z1%mpc(0) = mpwdsm6
    z1%mpc(l2) = mpwdsm6
    mp_addzdc%mpc(0) = mpwdsm6
    mp_addzdc%mpc(l2) = mpwdsm6
    call mpdmc40 (dble (dcb), 0, z1%mpc, mpnw)
    call mpdmc40 (aimag (dcb), 0, z1%mpc(l2:), mpnw)
    call mpcadd (za%mpc, z1%mpc, mp_addzdc%mpc, mpnw)
    return
  end function

  function mp_addrz (ra, zb)
    implicit none
    type (mp_complexm):: mp_addrz
    type (mp_realm), intent (in):: ra
    type (mp_complexm), intent (in):: zb
    integer l2, l3, mpnw
    l2 = zb%mpc(0)
    mpnw = max (int (ra%mpr(1)), int (zb%mpc(1)), int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwdsm)
    l3 = mpwdsm6
    mp_addrz%mpc(0) = mpwdsm6
    mp_addrz%mpc(l3) = mpwdsm6
    call mpadd (ra%mpr, zb%mpc, mp_addrz%mpc, mpnw)
    call mpeq (zb%mpc(l2:), mp_addrz%mpc(l3:), mpnw)
    return
  end function

  function mp_addzr (za, rb)
    implicit none
    type (mp_complexm):: mp_addzr
    type (mp_complexm), intent (in):: za
    type (mp_realm), intent (in):: rb
    integer l1, l3, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwdsm)
    l3 = mpwdsm6
    mp_addzr%mpc(0) = mpwdsm6
    mp_addzr%mpc(l3) = mpwdsm6
    call mpadd (za%mpc, rb%mpr, mp_addzr%mpc, mpnw)
    call mpeq (za%mpc(l1:), mp_addzr%mpc(l3:), mpnw)
    return
  end function

!  Subtraction routines:

  function mp_subrr (ra, rb)
    implicit none
    type (mp_realm):: mp_subrr
    type (mp_realm), intent (in):: ra, rb
    integer mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwdsm)
    mp_subrr%mpr(0) = mpwdsm6
    call mpsub (ra%mpr, rb%mpr, mp_subrr%mpr, mpnw)
    return
  end function

  function mp_subdr (da, rb)
    implicit none
    type (mp_realm):: mp_subdr
    real (mprknd), intent (in):: da
    type (mp_realm), intent (in):: rb
    type (mp_realm) r1
    integer i1, mpnw
    mpnw = min (int (rb%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    mp_subdr%mpr(0) = mpwdsm6
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpsub (r1%mpr, rb%mpr, mp_subdr%mpr, mpnw)
    return
  end function

  function mp_subrd (ra, db)
    implicit none
    type (mp_realm):: mp_subrd
    type (mp_realm), intent (in):: ra
    real (mprknd), intent (in):: db
    type (mp_realm) r1
    integer i1, mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    mp_subrd%mpr(0) = mpwdsm6
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpsub (ra%mpr, r1%mpr, mp_subrd%mpr, mpnw)
    return
  end function

  function mp_subir (ia, rb)
    implicit none
    type (mp_realm):: mp_subir
    integer, intent (in):: ia
    type (mp_realm), intent (in):: rb
    type (mp_realm) r1
    real (mprknd) da
    integer i1, mpnw
    mpnw = min (int (rb%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    mp_subir%mpr(0) = mpwdsm6
    da = ia
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpsub (r1%mpr, rb%mpr, mp_subir%mpr, mpnw)
    return
  end function

  function mp_subri (ra, ib)
    implicit none
    type (mp_realm):: mp_subri
    type (mp_realm), intent (in):: ra
    integer, intent (in):: ib
    type (mp_realm) r1
    real (mprknd) db
    integer i1, mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    mp_subri%mpr(0) = mpwdsm6
    db = ib
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpsub (ra%mpr, r1%mpr, mp_subri%mpr, mpnw)
    return
  end function

  function mp_subzz (za, zb)
    implicit none
    type (mp_complexm):: mp_subzz
    type (mp_complexm), intent (in):: za, zb
    integer l1, l2, l3, mpnw
    l1 = za%mpc(0)
    l2 = zb%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (zb%mpc(1)), &
      int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwdsm)
    l3 = mpwdsm6
    mp_subzz%mpc(0) = mpwdsm6
    mp_subzz%mpc(l3) = mpwdsm6
    call mpcsub (za%mpc, zb%mpc, mp_subzz%mpc, mpnw)
    return
  end function

  function mp_subdz (da, zb)
    implicit none
    type (mp_complexm):: mp_subdz
    real (mprknd) , intent (in):: da
    type (mp_complexm), intent (in):: zb
    integer l1, l2, mpnw
    real (mprknd) d1
    type (mp_complexm) z1
    l1 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    z1%mpc(0) = mpwdsm6
    z1%mpc(l2) = mpwdsm6
    mp_subdz%mpc(0) = mpwdsm6
    mp_subdz%mpc(l2) = mpwdsm6
    d1 = 0.d0
    call mpdmc40 (da, 0, z1%mpc, mpnw)
    call mpdmc40 (d1, 0, z1%mpc(l2:), mpnw)
    call mpcsub (z1%mpc, zb%mpc, mp_subdz%mpc, mpnw)
    return
  end function

  function mp_subzd (za, db)
    implicit none
    type (mp_complexm):: mp_subzd
    type (mp_complexm), intent (in):: za
    real (mprknd), intent (in):: db
    integer l1, l2, mpnw
    real (mprknd) d1
    type (mp_complexm) z1
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    z1%mpc(0) = mpwdsm6
    z1%mpc(l2) = mpwdsm6
    mp_subzd%mpc(0) = mpwdsm6
    mp_subzd%mpc(l2) = mpwdsm6
    d1 = 0.d0
    call mpdmc40 (db, 0, z1%mpc, mpnw)
    call mpdmc40 (d1, 0, z1%mpc(l2:), mpnw)
    call mpcsub (za%mpc, z1%mpc, mp_subzd%mpc, mpnw)
    return
  end function

  function mp_subdcz (dca, zb)
    implicit none
    type (mp_complexm):: mp_subdcz
    complex (kind (0.d0)), intent (in):: dca
    type (mp_complexm), intent (in):: zb
    integer l1, l2, mpnw
    type (mp_complexm) z1
    l1 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    z1%mpc(0) = mpwdsm6
    z1%mpc(l2) = mpwdsm6
    mp_subdcz%mpc(0) = mpwdsm6
    mp_subdcz%mpc(l2) = mpwdsm6
    call mpdmc40 (dble (dca), 0, z1%mpc, mpnw)
    call mpdmc40 (aimag (dca), 0, z1%mpc(l2:), mpnw)
    call mpcsub (z1%mpc, zb%mpc, mp_subdcz%mpc, mpnw)
    return
  end function

  function mp_subzdc (za, dcb)
    implicit none
    type (mp_complexm):: mp_subzdc
    type (mp_complexm), intent (in):: za
    complex (kind (0.d0)), intent (in):: dcb
    integer l1, l2, mpnw
    type (mp_complexm) z1
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    z1%mpc(0) = mpwdsm6
    z1%mpc(l2) = mpwdsm6
    mp_subzdc%mpc(0) = mpwdsm6
    mp_subzdc%mpc(l2) = mpwdsm6
    call mpdmc40 (dble (dcb), 0, z1%mpc, mpnw)
    call mpdmc40 (aimag (dcb), 0, z1%mpc(l2:), mpnw)
    call mpcsub (za%mpc, z1%mpc, mp_subzdc%mpc, mpnw)
    return
  end function

  function mp_subrz (ra, zb)
    implicit none
    type (mp_complexm):: mp_subrz
    type (mp_realm), intent (in):: ra
    type (mp_complexm), intent (in):: zb
    integer l2, l3, mpnw
    l2 = zb%mpc(0)
    mpnw = max (int (ra%mpr(1)), int (zb%mpc(1)), int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwdsm)
    l3 = mpwdsm6
    mp_subrz%mpc(0) = mpwdsm6
    mp_subrz%mpc(l3) = mpwdsm6
    call mpsub (ra%mpr, zb%mpc, mp_subrz%mpc, mpnw)
    call mpeq (zb%mpc(l2:), mp_subrz%mpc(l3:), mpnw)
    mp_subrz%mpc(l3+2) = - mp_subrz%mpc(l3+2)
    return
  end function

  function mp_subzr (za, rb)
    implicit none
    type (mp_complexm):: mp_subzr
    type (mp_complexm), intent (in):: za
    type (mp_realm), intent (in):: rb
    integer l1, l3, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwdsm)
    l3 = mpwdsm6
    mp_subzr%mpc(0) = mpwdsm6
    mp_subzr%mpc(l3) = mpwdsm6
    call mpsub (za%mpc, rb%mpr, mp_subzr%mpc, mpnw)
    call mpeq (za%mpc(l1:), mp_subzr%mpc(l3:), mpnw)
    return
  end function

!  Negation routines:

  function mp_negr (ra)
    implicit none
    type (mp_realm):: mp_negr
    type (mp_realm), intent (in):: ra
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_negr%mpr(0) = mpwdsm6
    call mpneg (ra%mpr, mp_negr%mpr, mpnw)
    return
  end function

  function mp_negz (za)
    implicit none
    type (mp_complexm):: mp_negz
    type (mp_complexm), intent (in):: za
    integer l1, l2, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    mp_negz%mpc(0) = mpwdsm6
    mp_negz%mpc(l2) = mpwdsm6
    call mpceq (za%mpc, mp_negz%mpc, mpnw)
    mp_negz%mpc(2) = - mp_negz%mpc(2)
    mp_negz%mpc(l2+2) = - mp_negz%mpc(l2+2)
    return
  end function

!  Multiplication routines:

  function mp_mulrr (ra, rb)
    implicit none
    type (mp_realm):: mp_mulrr
    type (mp_realm), intent (in):: ra, rb
    integer mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwdsm)
    mp_mulrr%mpr(0) = mpwdsm6
    call mpmul (ra%mpr, rb%mpr, mp_mulrr%mpr, mpnw)
    return
  end function

  function mp_muldr (da, rb)
    implicit none
    type (mp_realm):: mp_muldr
    real (mprknd), intent (in):: da
    type (mp_realm), intent (in):: rb
    integer mpnw
    mpnw = min (int (rb%mpr(1)), mpwdsm)
    mp_muldr%mpr(0) = mpwdsm6
    call mpmuld40 (rb%mpr, da, mp_muldr%mpr, mpnw)
    return
  end function

  function mp_mulrd (ra, db)
    implicit none
    type (mp_realm):: mp_mulrd
    type (mp_realm), intent (in):: ra
    real (mprknd), intent (in):: db
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_mulrd%mpr(0) = mpwdsm6
    call mpmuld40 (ra%mpr, db, mp_mulrd%mpr, mpnw)
    return
  end function

  function mp_mulir (ia, rb)
    implicit none
    type (mp_realm):: mp_mulir
    integer, intent (in):: ia
    type (mp_realm), intent (in):: rb
    real (mprknd) da
    integer mpnw
    mpnw = min (int (rb%mpr(1)), mpwdsm)
    mp_mulir%mpr(0) = mpwdsm6
    da = ia
    call mpmuld40 (rb%mpr, da, mp_mulir%mpr, mpnw)
    return
  end function

  function mp_mulri (ra, ib)
    implicit none
    type (mp_realm):: mp_mulri
    type (mp_realm), intent (in):: ra
    integer, intent (in):: ib
    real (mprknd) db
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_mulri%mpr(0) = mpwdsm6
    db = ib
    call mpmuld40 (ra%mpr, db, mp_mulri%mpr, mpnw)
    return
  end function

  function mp_mulzz (za, zb)
    implicit none
    type (mp_complexm):: mp_mulzz
    type (mp_complexm), intent (in):: za, zb
    integer l1, l2, l3, mpnw
    l1 = za%mpc(0)
    l2 = zb%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (zb%mpc(1)), &
      int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwdsm)
    l3 = mpwdsm6
    mp_mulzz%mpc(0) = mpwdsm6
    mp_mulzz%mpc(l3) = mpwdsm6
    call mpcmul (za%mpc, zb%mpc, mp_mulzz%mpc, mpnw)
    return
  end function

  function mp_muldz (da, zb)
    implicit none
    type (mp_complexm):: mp_muldz
    real (mprknd) , intent (in):: da
    type (mp_complexm), intent (in):: zb
    integer l1, l2, mpnw
    l1 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    mp_muldz%mpc(0) = mpwdsm6
    mp_muldz%mpc(l2) = mpwdsm6
    call mpmuld40 (zb%mpc, da, mp_muldz%mpc, mpnw)
    call mpmuld40 (zb%mpc(l1:), da, mp_muldz%mpc(l2:), mpnw)
    return
  end function

  function mp_mulzd (za, db)
    implicit none
    type (mp_complexm):: mp_mulzd
    type (mp_complexm), intent (in):: za
    real (mprknd), intent (in):: db
    integer l1, l2, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    mp_mulzd%mpc(0) = mpwdsm6
    mp_mulzd%mpc(l2) = mpwdsm6
    call mpmuld40 (za%mpc, db, mp_mulzd%mpc, mpnw)
    call mpmuld40 (za%mpc(l1:), db, mp_mulzd%mpc(l2:), mpnw)
    return
  end function

  function mp_muldcz (dca, zb)
    implicit none
    type (mp_complexm):: mp_muldcz
    complex (kind (0.d0)), intent (in):: dca
    type (mp_complexm), intent (in):: zb
    integer l1, l2, mpnw
    type (mp_complexm) z1
    l1 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    z1%mpc(0) = mpwdsm6
    z1%mpc(l2) = mpwdsm6
    mp_muldcz%mpc(0) = mpwdsm6
    mp_muldcz%mpc(l2) = mpwdsm6
    call mpdmc40 (dble (dca), 0, z1%mpc, mpnw)
    call mpdmc40 (aimag (dca), 0, z1%mpc(l2:), mpnw)
    call mpcmul (z1%mpc, zb%mpc, mp_muldcz%mpc, mpnw)
    return
  end function

  function mp_mulzdc (za, dcb)
    implicit none
    type (mp_complexm):: mp_mulzdc
    type (mp_complexm), intent (in):: za
    complex (kind (0.d0)), intent (in):: dcb
    integer l1, l2, mpnw
    type (mp_complexm) z1
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    z1%mpc(0) = mpwdsm6
    z1%mpc(l2) = mpwdsm6
    mp_mulzdc%mpc(0) = mpwdsm6
    mp_mulzdc%mpc(l2) = mpwdsm6
    call mpdmc40 (dble (dcb), 0, z1%mpc, mpnw)
    call mpdmc40 (aimag (dcb), 0, z1%mpc(l2:), mpnw)
    call mpcmul (za%mpc, z1%mpc, mp_mulzdc%mpc, mpnw)
    return
  end function

  function mp_mulrz (ra, zb)
    implicit none
    type (mp_complexm):: mp_mulrz
    type (mp_realm), intent (in):: ra
    type (mp_complexm), intent (in):: zb
    integer l2, l3, mpnw
    l2 = zb%mpc(0)
    mpnw = max (int (ra%mpr(1)), int (zb%mpc(1)), int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwdsm)
    l3 = mpwdsm6
    mp_mulrz%mpc(0) = mpwdsm6
    mp_mulrz%mpc(l3) = mpwdsm6
    call mpmul (ra%mpr, zb%mpc, mp_mulrz%mpc, mpnw)
    call mpmul (ra%mpr, zb%mpc(l2:), mp_mulrz%mpc(l3:), mpnw)
    return
  end function

  function mp_mulzr (za, rb)
    implicit none
    type (mp_complexm):: mp_mulzr
    type (mp_complexm), intent (in):: za
    type (mp_realm), intent (in):: rb
    integer l1, l3, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwdsm)
    l3 = mpwdsm6
    mp_mulzr%mpc(0) = mpwdsm6
    mp_mulzr%mpc(l3) = mpwdsm6
    call mpmul (za%mpc, rb%mpr, mp_mulzr%mpc, mpnw)
    call mpmul (za%mpc(l1:), rb%mpr, mp_mulzr%mpc(l3:), mpnw)
    return
  end function

!  Division routines:

  function mp_divrr (ra, rb)
    implicit none
    type (mp_realm):: mp_divrr
    type (mp_realm), intent (in):: ra, rb
    integer mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwdsm)
    mp_divrr%mpr(0) = mpwdsm6
    call mpdiv (ra%mpr, rb%mpr, mp_divrr%mpr, mpnw)
    return
  end function

  function mp_divdr (da, rb)
    implicit none
    type (mp_realm):: mp_divdr
    real (mprknd), intent (in):: da
    type (mp_realm), intent (in):: rb
    type (mp_realm) r1
    integer i1, mpnw
    mpnw = min (int (rb%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    mp_divdr%mpr(0) = mpwdsm6
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpdiv (r1%mpr, rb%mpr, mp_divdr%mpr, mpnw)
    return
  end function

  function mp_divrd (ra, db)
    implicit none
    type (mp_realm):: mp_divrd
    type (mp_realm), intent (in):: ra
    real (mprknd), intent (in):: db
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_divrd%mpr(0) = mpwdsm6
    call mpdivd40 (ra%mpr, db, mp_divrd%mpr, mpnw)
    return
  end function

  function mp_divir (ia, rb)
    implicit none
    type (mp_realm):: mp_divir
    integer, intent (in):: ia
    type (mp_realm), intent (in):: rb
    type (mp_realm) r1
    real (mprknd) da
    integer i1, mpnw
    mpnw = min (int (rb%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    mp_divir%mpr(0) = mpwdsm6
    da = ia
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpdiv (r1%mpr, rb%mpr, mp_divir%mpr, mpnw)
    return
  end function

  function mp_divri (ra, ib)
    implicit none
    type (mp_realm):: mp_divri
    type (mp_realm), intent (in):: ra
    integer, intent (in):: ib
    real (mprknd) db
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_divri%mpr(0) = mpwdsm6
    db = ib
    call mpdivd40 (ra%mpr, db, mp_divri%mpr, mpnw)
    return
  end function

  function mp_divzz (za, zb)
    implicit none
    type (mp_complexm):: mp_divzz
    type (mp_complexm), intent (in):: za, zb
    integer l1, l2, l3, mpnw
    l1 = za%mpc(0)
    l2 = zb%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (zb%mpc(1)), &
      int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwdsm)
    l3 = mpwdsm6
    mp_divzz%mpc(0) = mpwdsm6
    mp_divzz%mpc(l3) = mpwdsm6
    call mpcdiv (za%mpc, zb%mpc, mp_divzz%mpc, mpnw)
    return
  end function

  function mp_divdz (da, zb)
    implicit none
    type (mp_complexm):: mp_divdz
    real (mprknd), intent (in):: da
    type (mp_complexm), intent (in):: zb
    real (mprknd) d1
    integer l1, l2, mpnw
    type (mp_complexm) z1
    l1 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    z1%mpc(0) = mpwdsm6
    z1%mpc(l2) = mpwdsm6
    mp_divdz%mpc(0) = mpwdsm6
    mp_divdz%mpc(l2) = mpwdsm6
    d1 = 0.d0
    call mpdmc40 (da, 0, z1%mpc, mpnw)
    call mpdmc40 (d1, 0, z1%mpc(l2:), mpnw)
    call mpcdiv (z1%mpc, zb%mpc, mp_divdz%mpc, mpnw)
    return
  end function

  function mp_divzd (za, db)
    implicit none
    type (mp_complexm):: mp_divzd
    type (mp_complexm), intent (in):: za
    real (mprknd), intent (in):: db
    real (mprknd) d1
    integer l1, l2, mpnw
    type (mp_complexm) z1
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    z1%mpc(0) = mpwdsm6
    z1%mpc(l2) = mpwdsm6
    mp_divzd%mpc(0) = mpwdsm6
    mp_divzd%mpc(l2) = mpwdsm6
    d1 = 0.d0
    call mpdmc40 (db, 0, z1%mpc, mpnw)
    call mpdmc40 (d1, 0, z1%mpc(l2:), mpnw)
    call mpcdiv (za%mpc, z1%mpc, mp_divzd%mpc, mpnw)
    return
  end function

  function mp_divdcz (dca, zb)
    implicit none
    type (mp_complexm):: mp_divdcz
    complex (kind (0.d0)), intent (in):: dca
    type (mp_complexm), intent (in):: zb
    integer l1, l2, mpnw
    type (mp_complexm) z1
    l1 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    z1%mpc(0) = mpwdsm6
    z1%mpc(l2) = mpwdsm6
    mp_divdcz%mpc(0) = mpwdsm6
    mp_divdcz%mpc(l2) = mpwdsm6
    call mpdmc40 (dble (dca), 0, z1%mpc, mpnw)
    call mpdmc40 (aimag (dca), 0, z1%mpc(l2:), mpnw)
    call mpcdiv (z1%mpc, zb%mpc, mp_divdcz%mpc, mpnw)
    return
  end function

  function mp_divzdc (za, dcb)
    implicit none
    type (mp_complexm):: mp_divzdc
    type (mp_complexm), intent (in):: za
    complex (kind (0.d0)), intent (in):: dcb
    integer l1, l2, mpnw
    type (mp_complexm) z1
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    z1%mpc(0) = mpwdsm6
    z1%mpc(l2) = mpwdsm6
    mp_divzdc%mpc(0) = mpwdsm6
    mp_divzdc%mpc(l2) = mpwdsm6
    call mpdmc40 (dble (dcb), 0, z1%mpc, mpnw)
    call mpdmc40 (aimag (dcb), 0, z1%mpc(l2:), mpnw)
    call mpcdiv (za%mpc, z1%mpc, mp_divzdc%mpc, mpnw)
    return
  end function

  function mp_divrz (ra, zb)
    implicit none
    type (mp_complexm):: mp_divrz
    type (mp_realm), intent (in):: ra
    type (mp_complexm), intent (in):: zb
    type (mp_realm):: r1, r2, r3
    integer l2, l3, mpnw
    l2 = zb%mpc(0)
    mpnw = max (int (ra%mpr(1)), int (zb%mpc(1)), int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwdsm)
    r1%mpr(0) = mpwdsm6
    r2%mpr(0) = mpwdsm6
    r3%mpr(0) = mpwdsm6
    l3 = mpwdsm6
    mp_divrz%mpc(0) = mpwdsm6
    mp_divrz%mpc(l3) = mpwdsm6
    call mpmul (zb%mpc, zb%mpc, r1%mpr, mpnw)
    call mpmul (zb%mpc(l2:), zb%mpc(l2:), r2%mpr, mpnw)
    call mpadd (r1%mpr, r2%mpr, r3%mpr, mpnw)
    call mpmul (ra%mpr, zb%mpc, r1%mpr, mpnw)
    call mpdiv (r1%mpr, r3%mpr, mp_divrz%mpc, mpnw)
    call mpmul (ra%mpr, zb%mpc(l2:), r1%mpr, mpnw)
    call mpdiv (r1%mpr, r3%mpr, mp_divrz%mpc(l3:), mpnw)
    mp_divrz%mpc(l3+2) = - mp_divrz%mpc(l3+2)
    return
  end function

  function mp_divzr (za, rb)
    implicit none
    type (mp_complexm):: mp_divzr
    type (mp_complexm), intent (in):: za
    type (mp_realm), intent (in):: rb
    integer l1, l3, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwdsm)
    l3 = mpwdsm6
    mp_divzr%mpc(0) = mpwdsm6
    mp_divzr%mpc(l3) = mpwdsm6
    call mpdiv (za%mpc, rb%mpr, mp_divzr%mpc, mpnw)
    call mpdiv (za%mpc(l1:), rb%mpr, mp_divzr%mpc(l3:), mpnw)
    return
  end function

!  Exponentiation routines:

  function mp_expri (ra, ib)
    implicit none
    type (mp_realm):: mp_expri
    type (mp_realm), intent (in):: ra
    integer, intent (in):: ib
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_expri%mpr(0) = mpwdsm6
    call mpnpwr (ra%mpr, ib, mp_expri%mpr, mpnw)
    return
  end function

  function mp_exprr (ra, rb)
    implicit none
    type (mp_realm):: mp_exprr
    type (mp_realm), intent (in):: ra, rb
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_exprr%mpr(0) = mpwdsm6
    call mppower (ra%mpr, rb%mpr, mp_exprr%mpr, mpnw)
    return
  end function

  function mp_expzi (za, ib)
    implicit none
    type (mp_complexm):: mp_expzi
    type (mp_complexm), intent (in):: za
    integer, intent (in):: ib
    integer l1, l2, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    mp_expzi%mpc(0) = mpwdsm6
    mp_expzi%mpc(l2) = mpwdsm6
    call mpcnpwr (za%mpc, ib, mp_expzi%mpc, mpnw)
    return
  end function

  function mp_expzz (za, zb)
    implicit none
    type (mp_complexm):: mp_expzz
    type (mp_complexm), intent (in):: za, zb
    integer l1, l2, l3, mpnw
    l1 = za%mpc(0)
    l2 = zb%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (zb%mpc(1)), &
      int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwdsm)
    l3 = mpwdsm6
    mp_expzz%mpc(0) = mpwdsm6
    mp_expzz%mpc(l3) = mpwdsm6
    call mpcpowcc (za%mpc, zb%mpc, mp_expzz%mpc, mpnw)
    return
  end function

  function mp_exprz (ra, zb)
    implicit none
    type (mp_complexm):: mp_exprz
    type (mp_realm), intent (in):: ra
    type (mp_complexm), intent (in):: zb
    integer l2, l3, mpnw
    l2 = zb%mpc(0)
    mpnw = max (int (ra%mpr(1)), int (zb%mpc(1)), int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwdsm)
    l3 = mpwdsm6
    mp_exprz%mpc(0) = mpwdsm6
    mp_exprz%mpc(l3) = mpwdsm6
    call mpcpowrc (ra%mpr, zb%mpc, mp_exprz%mpc, mpnw)
    return
  end function

  function mp_expzr (za, rb)
    implicit none
    type (mp_complexm):: mp_expzr
    type (mp_complexm), intent (in):: za
    type (mp_realm), intent (in):: rb
    integer l1, l3, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwdsm)
    l3 = mpwdsm6
    mp_expzr%mpc(0) = mpwdsm6
    mp_expzr%mpc(l3) = mpwdsm6
    call mpcpowcr (za%mpc, rb%mpr, mp_expzr%mpc, mpnw)
    return
  end function

!  Equality test routines:

  function mp_eqtrr (ra, rb)
    implicit none
    logical mp_eqtrr
    type (mp_realm), intent (in):: ra, rb
    integer ic, mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwdsm)
    call mpcpr (ra%mpr, rb%mpr, ic, mpnw)
    if (ic == 0) then
      mp_eqtrr = .true.
    else
      mp_eqtrr = .false.
    endif
    return
  end function

  function mp_eqtdr (da, rb)
    implicit none
    logical mp_eqtdr
    real (mprknd), intent (in):: da
    type (mp_realm), intent (in):: rb
    type (mp_realm) r1
    integer i1, ic, mpnw
    mpnw = min (int (rb%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpcpr (r1%mpr, rb%mpr, ic, mpnw)
    if (ic == 0) then
      mp_eqtdr = .true.
    else
      mp_eqtdr = .false.
    endif
    return
  end function

  function mp_eqtrd (ra, db)
    implicit none
    logical mp_eqtrd
    type (mp_realm), intent (in):: ra
    real (mprknd), intent (in):: db
    type (mp_realm) r1
    integer i1, ic, mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpcpr (ra%mpr, r1%mpr, ic, mpnw)
    if (ic == 0) then
      mp_eqtrd = .true.
    else
      mp_eqtrd = .false.
    endif
    return
  end function

  function mp_eqtir (ia, rb)
    implicit none
    logical mp_eqtir
    integer, intent (in):: ia
    type (mp_realm), intent (in):: rb
    type (mp_realm) r1
    real (mprknd) da
    integer i1, ic, mpnw
    mpnw = min (int (rb%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    da = ia
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpcpr (r1%mpr, rb%mpr, ic, mpnw)
    if (ic == 0) then
      mp_eqtir = .true.
    else
      mp_eqtir = .false.
    endif
    return
  end function

  function mp_eqtri (ra, ib)
    implicit none
    logical mp_eqtri
    type (mp_realm), intent (in):: ra
    integer, intent (in):: ib
    type (mp_realm) r1
    real (mprknd) db
    integer i1, ic, mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    db = ib
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpcpr (ra%mpr, r1%mpr, ic, mpnw)
    if (ic == 0) then
      mp_eqtri = .true.
    else
      mp_eqtri = .false.
    endif
    return
  end function

  function mp_eqtzz (za, zb)
    implicit none
    logical mp_eqtzz
    type (mp_complexm), intent (in):: za, zb
    integer ic1, ic2, l1, l2, mpnw
    l1 = za%mpc(0)
    l2 = zb%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (zb%mpc(1)), &
      int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwdsm)
    call mpcpr (za%mpc, zb%mpc, ic1, mpnw)
    call mpcpr (za%mpc(l1:), zb%mpc(l2:), ic2, mpnw)
    if (ic1 == 0 .and. ic2 == 0) then
      mp_eqtzz = .true.
    else
      mp_eqtzz = .false.
    endif
    return
  end function

  function mp_eqtdz (da, zb)
    implicit none
    logical mp_eqtdz
    real (mprknd), intent (in):: da
    type (mp_complexm), intent (in):: zb
    integer ic1, ic2, l1, l2, mpnw
    type (mp_complexm) z1
    real (mprknd) d1
    l1 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    z1%mpc(0) = mpwdsm6
    z1%mpc(l2) = mpwdsm6
    d1 = 0.d0
    call mpdmc40 (da, 0, z1%mpc, mpnw)
    call mpdmc40 (d1, 0, z1%mpc(l2:), mpnw)
    call mpcpr (z1%mpc, zb%mpc, ic1, mpnw)
    call mpcpr (z1%mpc(l2:), zb%mpc(l1:), ic2, mpnw)
    if (ic1 == 0 .and. ic2 == 0) then
      mp_eqtdz = .true.
    else
      mp_eqtdz = .false.
    endif
    return
  end function

  function mp_eqtzd (za, db)
    implicit none
    logical mp_eqtzd
    type (mp_complexm), intent (in):: za
    real (mprknd), intent (in):: db
    integer ic1, ic2, l1, l2, mpnw
    type (mp_complexm) z1
    real (mprknd) d1
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    z1%mpc(0) = mpwdsm6
    z1%mpc(l2) = mpwdsm6
    d1 = 0.d0
    call mpdmc40 (db, 0, z1%mpc, mpnw)
    call mpdmc40 (d1, 0, z1%mpc(l2:), mpnw)
    call mpcpr (za%mpc, z1%mpc, ic1, mpnw)
    call mpcpr (za%mpc(l1:), z1%mpc(l2:), ic2, mpnw)
    if (ic1 == 0 .and. ic2 == 0) then
      mp_eqtzd = .true.
    else
      mp_eqtzd = .false.
    endif
    return
  end function

  function mp_eqtdcz (dca, zb)
    implicit none
    logical mp_eqtdcz
    complex (kind (0.d0)), intent (in):: dca
    type (mp_complexm), intent (in):: zb
    integer ic1, ic2, l1, l2, mpnw
    type (mp_complexm) z1
    l1 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    z1%mpc(0) = mpwdsm6
    z1%mpc(l2) = mpwdsm6
    call mpdmc40 (dble (dca), 0, z1%mpc, mpnw)
    call mpdmc40 (aimag (dca), 0, z1%mpc(l2:), mpnw)
    call mpcpr (z1%mpc, zb%mpc, ic1, mpnw)
    call mpcpr (z1%mpc(l2:), zb%mpc(l1:), ic2, mpnw)
    if (ic1 == 0 .and. ic2 == 0) then
      mp_eqtdcz = .true.
    else
      mp_eqtdcz = .false.
    endif
    return
  end function

  function mp_eqtzdc (za, dcb)
    implicit none
    logical mp_eqtzdc
    type (mp_complexm), intent (in):: za
    complex (kind (0.d0)), intent (in):: dcb
    integer ic1, ic2, l1, l2, mpnw
    type (mp_complexm) z1
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    z1%mpc(0) = mpwdsm6
    z1%mpc(l2) = mpwdsm6
    call mpdmc40 (dble (dcb), 0, z1%mpc, mpnw)
    call mpdmc40 (aimag (dcb), 0, z1%mpc(l2:), mpnw)
    call mpcpr (za%mpc, z1%mpc, ic1, mpnw)
    call mpcpr (za%mpc(l1:), z1%mpc(l2:), ic2, mpnw)
    if (ic1 == 0 .and. ic2 == 0) then
      mp_eqtzdc = .true.
    else
      mp_eqtzdc = .false.
    endif
    return
  end function

  function mp_eqtrz (ra, zb)
    implicit none
    logical mp_eqtrz
    type (mp_realm), intent (in):: ra
    type (mp_complexm), intent (in):: zb
    integer ic1, ic2, l2, mpnw
    l2 = zb%mpc(0)
    mpnw = max (int (ra%mpr(1)), int (zb%mpc(1)), int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwdsm)
    call mpcpr (ra%mpr, zb%mpc, ic1, mpnw)
    ic2 = int (zb%mpc(l2+2))
    if (ic1 == 0 .and. ic2 == 0) then
      mp_eqtrz = .true.
    else
      mp_eqtrz = .false.
    endif
    return
  end function

  function mp_eqtzr (za, rb)
    implicit none
    logical mp_eqtzr
    type (mp_complexm), intent (in):: za
    type (mp_realm), intent (in):: rb
    integer ic1, ic2, l1, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwdsm)
    call mpcpr (za%mpc, rb%mpr, ic1, mpnw)
    ic2 = int (za%mpc(l1+2))
    if (ic1 == 0 .and. ic2 == 0) then
      mp_eqtzr = .true.
    else
      mp_eqtzr = .false.
    endif
    return
  end function

!  Non-equality test routines:

  function mp_netrr (ra, rb)
    implicit none
    logical mp_netrr
    type (mp_realm), intent (in):: ra, rb
    integer ic, mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwdsm)
    call mpcpr (ra%mpr, rb%mpr, ic, mpnw)
    if (ic /= 0) then
      mp_netrr = .true.
    else
      mp_netrr = .false.
    endif
    return
  end function

  function mp_netdr (da, rb)
    implicit none
    logical mp_netdr
    real (mprknd), intent (in):: da
    type (mp_realm), intent (in):: rb
    type (mp_realm) r1
    integer i1, ic, mpnw
    mpnw = min (int (rb%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpcpr (r1%mpr, rb%mpr, ic, mpnw)
    if (ic /= 0) then
      mp_netdr = .true.
    else
      mp_netdr = .false.
    endif
    return
  end function

  function mp_netrd (ra, db)
    implicit none
    logical mp_netrd
    type (mp_realm), intent (in):: ra
    real (mprknd), intent (in):: db
    type (mp_realm) r1
    integer i1, ic, mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpcpr (ra%mpr, r1%mpr, ic, mpnw)
    if (ic /= 0) then
      mp_netrd = .true.
    else
      mp_netrd = .false.
    endif
    return
  end function

  function mp_netir (ia, rb)
    implicit none
    logical mp_netir
    integer, intent (in):: ia
    type (mp_realm), intent (in):: rb
    type (mp_realm) r1
    real (mprknd) da
    integer i1, ic, mpnw
    mpnw = min (int (rb%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    da = ia
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpcpr (r1%mpr, rb%mpr, ic, mpnw)
    if (ic /= 0) then
      mp_netir = .true.
    else
      mp_netir = .false.
    endif
    return
  end function

  function mp_netri (ra, ib)
    implicit none
    logical mp_netri
    type (mp_realm), intent (in):: ra
    integer, intent (in):: ib
    type (mp_realm) r1
    real (mprknd) db
    integer i1, ic, mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    db = ib
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpcpr (ra%mpr, r1%mpr, ic, mpnw)
    if (ic /= 0) then
      mp_netri = .true.
    else
      mp_netri = .false.
    endif
    return
  end function

  function mp_netzz (za, zb)
    implicit none
    logical mp_netzz
    type (mp_complexm), intent (in):: za, zb
    integer ic1, ic2, l1, l2, mpnw
    l1 = za%mpc(0)
    l2 = zb%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (zb%mpc(1)), &
      int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwdsm)
    call mpcpr (za%mpc, zb%mpc, ic1, mpnw)
    call mpcpr (za%mpc(l1:), zb%mpc(l2:), ic2, mpnw)
    if (ic1 /= 0 .or. ic2 /= 0) then
      mp_netzz = .true.
    else
      mp_netzz = .false.
    endif
    return
  end function

  function mp_netdz (da, zb)
    implicit none
    logical mp_netdz
    real (mprknd), intent (in):: da
    type (mp_complexm), intent (in):: zb
    integer ic1, ic2, l1, l2, mpnw
    type (mp_complexm) z1
    real (mprknd) d1
    l1 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    z1%mpc(0) = mpwdsm6
    z1%mpc(l2) = mpwdsm6
    d1 = 0.d0
    call mpdmc40 (da, 0, z1%mpc, mpnw)
    call mpdmc40 (d1, 0, z1%mpc(l2:), mpnw)
    call mpcpr (z1%mpc, zb%mpc, ic1, mpnw)
    call mpcpr (z1%mpc(l2:), zb%mpc(l1:), ic2, mpnw)
    if (ic1 /= 0 .or. ic2 /= 0) then
      mp_netdz = .true.
    else
      mp_netdz = .false.
    endif
    return
  end function

  function mp_netzd (za, db)
    implicit none
    logical mp_netzd
    type (mp_complexm), intent (in):: za
    real (mprknd), intent (in):: db
    integer ic1, ic2, l1, l2, mpnw
    type (mp_complexm) z1
    real (mprknd) d1
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    z1%mpc(0) = mpwdsm6
    z1%mpc(l2) = mpwdsm6
    d1 = 0.d0
    call mpdmc40 (db, 0, z1%mpc, mpnw)
    call mpdmc40 (d1, 0, z1%mpc(l2:), mpnw)
    call mpcpr (za%mpc, z1%mpc, ic1, mpnw)
    call mpcpr (za%mpc(l1:), z1%mpc(l2:), ic2, mpnw)
    if (ic1 /= 0 .or. ic2 /= 0) then
      mp_netzd = .true.
    else
      mp_netzd = .false.
    endif
    return
  end function

  function mp_netdcz (dca, zb)
    implicit none
    logical mp_netdcz
    complex (kind (0.d0)), intent (in):: dca
    type (mp_complexm), intent (in):: zb
    integer ic1, ic2, l1, l2, mpnw
    type (mp_complexm) z1
    l1 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    z1%mpc(0) = mpwdsm6
    z1%mpc(l2) = mpwdsm6
    call mpdmc40 (dble (dca), 0, z1%mpc, mpnw)
    call mpdmc40 (aimag (dca), 0, z1%mpc(l2:), mpnw)
    call mpcpr (z1%mpc, zb%mpc, ic1, mpnw)
    call mpcpr (z1%mpc(l2:), zb%mpc(l1:), ic2, mpnw)
    if (ic1 /= 0 .or. ic2 /= 0) then
      mp_netdcz = .true.
    else
      mp_netdcz = .false.
    endif
    return
  end function

  function mp_netzdc (za, dcb)
    implicit none
    logical mp_netzdc
    type (mp_complexm), intent (in):: za
    complex (kind (0.d0)), intent (in):: dcb
    integer ic1, ic2, l1, l2, mpnw
    type (mp_complexm) z1
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    z1%mpc(0) = mpwdsm6
    z1%mpc(l2) = mpwdsm6
    call mpdmc40 (dble (dcb), 0, z1%mpc, mpnw)
    call mpdmc40 (aimag (dcb), 0, z1%mpc(l2:), mpnw)
    call mpcpr (za%mpc, z1%mpc, ic1, mpnw)
    call mpcpr (za%mpc(l1:), z1%mpc(l2:), ic2, mpnw)
    if (ic1 /= 0 .or. ic2 /= 0) then
      mp_netzdc = .true.
    else
      mp_netzdc = .false.
    endif
    return
  end function

  function mp_netrz (ra, zb)
    implicit none
    logical mp_netrz
    type (mp_realm), intent (in):: ra
    type (mp_complexm), intent (in):: zb
    integer ic1, ic2, l2, mpnw
    l2 = zb%mpc(0)
    mpnw = max (int (ra%mpr(1)), int (zb%mpc(1)), int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwdsm)
    call mpcpr (ra%mpr, zb%mpc, ic1, mpnw)
    ic2 = int (zb%mpc(l2+2))
    if (ic1 /= 0 .or. ic2 /= 0) then
      mp_netrz = .true.
    else
      mp_netrz = .false.
    endif
    return
  end function

  function mp_netzr (za, rb)
    implicit none
    logical mp_netzr
    type (mp_complexm), intent (in):: za
    type (mp_realm), intent (in):: rb
    integer ic1, ic2, l1, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwdsm)
    call mpcpr (za%mpc, rb%mpr, ic1, mpnw)
    ic2 = int (za%mpc(l1+2))
    if (ic1 /= 0 .or. ic2 /= 0) then
      mp_netzr = .true.
    else
      mp_netzr = .false.
    endif
    return
  end function

!  Less-than-or-equal test routines:

  function mp_letrr (ra, rb)
    implicit none
    logical mp_letrr
    type (mp_realm), intent (in):: ra, rb
    integer ic, mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwdsm)
    call mpcpr (ra%mpr, rb%mpr, ic, mpnw)
    if (ic <= 0) then
      mp_letrr = .true.
    else
      mp_letrr = .false.
    endif
    return
  end function

  function mp_letdr (da, rb)
    implicit none
    logical mp_letdr
    real (mprknd), intent (in):: da
    type (mp_realm), intent (in):: rb
    type (mp_realm) r1
    integer i1, ic, mpnw
    mpnw = min (int (rb%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpcpr (r1%mpr, rb%mpr, ic, mpnw)
    if (ic <= 0) then
      mp_letdr = .true.
    else
      mp_letdr = .false.
    endif
    return
  end function

  function mp_letrd (ra, db)
    implicit none
    logical mp_letrd
    type (mp_realm), intent (in):: ra
    real (mprknd), intent (in):: db
    type (mp_realm) r1
    integer i1, ic, mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpcpr (ra%mpr, r1%mpr, ic, mpnw)
    if (ic <= 0) then
      mp_letrd = .true.
    else
      mp_letrd = .false.
    endif
    return
  end function

  function mp_letir (ia, rb)
    implicit none
    logical mp_letir
    integer, intent (in):: ia
    type (mp_realm), intent (in):: rb
    type (mp_realm) r1
    real (mprknd) da
    integer i1, ic, mpnw
    mpnw = min (int (rb%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    da = ia
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpcpr (r1%mpr, rb%mpr, ic, mpnw)
    if (ic <= 0) then
      mp_letir = .true.
    else
      mp_letir = .false.
    endif
    return
  end function

  function mp_letri (ra, ib)
    implicit none
    logical mp_letri
    type (mp_realm), intent (in):: ra
    integer, intent (in):: ib
    type (mp_realm) r1
    real (mprknd) db
    integer i1, ic, mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    db = ib
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpcpr (ra%mpr, r1%mpr, ic, mpnw)
    if (ic <= 0) then
      mp_letri = .true.
    else
      mp_letri = .false.
    endif
    return
  end function

!  Greater-than-or-equal test routines:

  function mp_getrr (ra, rb)
    implicit none
    logical mp_getrr
    type (mp_realm), intent (in):: ra, rb
    integer ic, mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwdsm)
    call mpcpr (ra%mpr, rb%mpr, ic, mpnw)
    if (ic >= 0) then
      mp_getrr = .true.
    else
      mp_getrr = .false.
    endif
    return
  end function

  function mp_getdr (da, rb)
    implicit none
    logical mp_getdr
    real (mprknd), intent (in):: da
    type (mp_realm), intent (in):: rb
    type (mp_realm) r1
    integer i1, ic, mpnw
    mpnw = min (int (rb%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpcpr (r1%mpr, rb%mpr, ic, mpnw)
    if (ic >= 0) then
      mp_getdr = .true.
    else
      mp_getdr = .false.
    endif
    return
  end function

  function mp_getrd (ra, db)
    implicit none
    logical mp_getrd
    type (mp_realm), intent (in):: ra
    real (mprknd), intent (in):: db
    type (mp_realm) r1
    integer i1, ic, mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpcpr (ra%mpr, r1%mpr, ic, mpnw)
    if (ic >= 0) then
      mp_getrd = .true.
    else
      mp_getrd = .false.
    endif
    return
  end function

  function mp_getir (ia, rb)
    implicit none
    logical mp_getir
    integer, intent (in):: ia
    type (mp_realm), intent (in):: rb
    type (mp_realm) r1
    real (mprknd) da
    integer i1, ic, mpnw
    mpnw = min (int (rb%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    da = ia
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpcpr (r1%mpr, rb%mpr, ic, mpnw)
    if (ic >= 0) then
      mp_getir = .true.
    else
      mp_getir = .false.
    endif
    return
  end function

  function mp_getri (ra, ib)
    implicit none
    logical mp_getri
    type (mp_realm), intent (in):: ra
    integer, intent (in):: ib
    type (mp_realm) r1
    real (mprknd) db
    integer i1, ic, mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    db = ib
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpcpr (ra%mpr, r1%mpr, ic, mpnw)
    if (ic >= 0) then
      mp_getri = .true.
    else
      mp_getri = .false.
    endif
    return
  end function

!  Less-than test routines:

  function mp_lttrr (ra, rb)
    implicit none
    logical mp_lttrr
    type (mp_realm), intent (in):: ra, rb
    integer ic, mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwdsm)
    call mpcpr (ra%mpr, rb%mpr, ic, mpnw)
    if (ic < 0) then
      mp_lttrr = .true.
    else
      mp_lttrr = .false.
    endif
    return
  end function

  function mp_lttdr (da, rb)
    implicit none
    logical mp_lttdr
    real (mprknd), intent (in):: da
    type (mp_realm), intent (in):: rb
    type (mp_realm) r1
    integer i1, ic, mpnw
    mpnw = min (int (rb%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpcpr (r1%mpr, rb%mpr, ic, mpnw)
    if (ic < 0) then
      mp_lttdr = .true.
    else
      mp_lttdr = .false.
    endif
    return
  end function

  function mp_lttrd (ra, db)
    implicit none
    logical mp_lttrd
    type (mp_realm), intent (in):: ra
    real (mprknd), intent (in):: db
    type (mp_realm) r1
    integer i1, ic, mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpcpr (ra%mpr, r1%mpr, ic, mpnw)
    if (ic < 0) then
      mp_lttrd = .true.
    else
      mp_lttrd = .false.
    endif
    return
  end function

  function mp_lttir (ia, rb)
    implicit none
    logical mp_lttir
    integer, intent (in):: ia
    type (mp_realm), intent (in):: rb
    type (mp_realm) r1
    real (mprknd) da
    integer i1, ic, mpnw
    mpnw = min (int (rb%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    da = ia
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpcpr (r1%mpr, rb%mpr, ic, mpnw)
    if (ic < 0) then
      mp_lttir = .true.
    else
      mp_lttir = .false.
    endif
    return
  end function

  function mp_lttri (ra, ib)
    implicit none
    logical mp_lttri
    type (mp_realm), intent (in):: ra
    integer, intent (in):: ib
    type (mp_realm) r1
    real (mprknd) db
    integer i1, ic, mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    db = ib
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpcpr (ra%mpr, r1%mpr, ic, mpnw)
    if (ic < 0) then
      mp_lttri = .true.
    else
      mp_lttri = .false.
    endif
    return
  end function

!  Greater-than test routines:

  function mp_gttrr (ra, rb)
    implicit none
    logical mp_gttrr
    type (mp_realm), intent (in):: ra, rb
    integer ic, mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwdsm)
    call mpcpr (ra%mpr, rb%mpr, ic, mpnw)
    if (ic > 0) then
      mp_gttrr = .true.
    else
      mp_gttrr = .false.
    endif
    return
  end function

  function mp_gttdr (da, rb)
    implicit none
    logical mp_gttdr
    real (mprknd), intent (in):: da
    type (mp_realm), intent (in):: rb
    type (mp_realm) r1
    integer i1, ic, mpnw
    mpnw = min (int (rb%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpcpr (r1%mpr, rb%mpr, ic, mpnw)
    if (ic > 0) then
      mp_gttdr = .true.
    else
      mp_gttdr = .false.
    endif
    return
  end function

  function mp_gttrd (ra, db)
    implicit none
    logical mp_gttrd
    type (mp_realm), intent (in):: ra
    real (mprknd), intent (in):: db
    type (mp_realm) r1
    integer i1, ic, mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpcpr (ra%mpr, r1%mpr, ic, mpnw)
    if (ic > 0) then
      mp_gttrd = .true.
    else
      mp_gttrd = .false.
    endif
    return
  end function

  function mp_gttir (ia, rb)
    implicit none
    logical mp_gttir
    integer, intent (in):: ia
    type (mp_realm), intent (in):: rb
    type (mp_realm) r1
    real (mprknd) da
    integer i1, ic, mpnw
    mpnw = min (int (rb%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    da = ia
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpcpr (r1%mpr, rb%mpr, ic, mpnw)
    if (ic > 0) then
      mp_gttir = .true.
    else
      mp_gttir = .false.
    endif
    return
  end function

  function mp_gttri (ra, ib)
    implicit none
    logical mp_gttri
    type (mp_realm), intent (in):: ra
    integer, intent (in):: ib
    type (mp_realm) r1
    real (mprknd) db
    integer i1, ic, mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    db = ib
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpcpr (ra%mpr, r1%mpr, ic, mpnw)
    if (ic > 0) then
      mp_gttri = .true.
    else
      mp_gttri = .false.
    endif
    return
  end function

!   Algebraic and transcendental function definitions, listed alphabetically:

  function mp_absr (ra)
    implicit none
    type (mp_realm):: mp_absr
    type (mp_realm), intent (in):: ra
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_absr%mpr(0) = mpwdsm6
    call mpabs (ra%mpr, mp_absr%mpr, mpnw)
    return
  end function

  function mp_absz (za)
    implicit none
    type (mp_realm):: mp_absz
    type (mp_complexm), intent (in):: za
    integer l1, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    mp_absz%mpr(0) = mpwdsm6
    call mpcabs (za%mpc, mp_absz%mpr, mpnw)
    return
  end function

  function mp_acos (ra)
    implicit none
    type (mp_realm):: mp_acos
    type (mp_realm), intent (in):: ra
    type (mp_realm) r1, r2, r3
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    r2%mpr(0) = mpwdsm6
    r3%mpr(0) = mpwdsm6
    call mpmul (ra%mpr, ra%mpr, r1%mpr, mpnw)
    call mpdmc (1.d0, 0, r2%mpr, mpnw)
    call mpsub (r2%mpr, r1%mpr, r3%mpr, mpnw)
    if (r3%mpr(2) < 0) then
      write (mpldb, 1)
1     format ('*** MP_ACOS: argument is not in (-1, 1).')
      call mpabrt ( 802)
    endif
    call mpsqrt (r3%mpr, r1%mpr, mpnw)
    mp_acos%mpr(0) = mpwdsm6
    call mpang (ra%mpr, r1%mpr, mp_acos%mpr, mpnw)
    return
  end function

  function mp_acosh (ra)
    implicit none
    type (mp_realm):: mp_acosh
    type (mp_realm), intent (in):: ra
    type (mp_realm) r1, r2, r3
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    r2%mpr(0) = mpwdsm6
    r3%mpr(0) = mpwdsm6
    mp_acosh%mpr(0) = mpwdsm6
    call mpmul (ra%mpr, ra%mpr, r1%mpr, mpnw)
    call mpdmc (1.d0, 0, r2%mpr, mpnw)
    call mpsub (r1%mpr, r2%mpr, r3%mpr, mpnw)
    if (r3%mpr(2) < 0) then
      write (mpldb, 1)
1     format ('*** MP_ACOSH: argument is not >= 1.')
      call mpabrt ( 803)
    endif
    call mpsqrt (r3%mpr, r1%mpr, mpnw)
    call mpadd (ra%mpr, r1%mpr, r2%mpr, mpnw)
    call mplog (r2%mpr, mp_acosh%mpr, mpnw)
    return
  end function

   function mp_agm (ra, rb)
    implicit none
    type (mp_realm):: mp_agm
    type (mp_realm), intent (in):: ra, rb
    integer mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwdsm)
    mp_agm%mpr(0) = mpwdsm6
    call mpagmr (rb%mpr, ra%mpr, mp_agm%mpr, mpnw)
    return
  end function

  function mp_aimag (za)
    implicit none
    type (mp_realm):: mp_aimag
    type (mp_complexm), intent (in):: za
    integer l1, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    mp_aimag%mpr(0) = mpwdsm6
    call mpeq (za%mpc(l1:), mp_aimag%mpr, mpnw)
    return
  end function

  function mp_aint (ra)
    implicit none
    type (mp_realm):: mp_aint
    type (mp_realm), intent (in):: ra
    type (mp_realm) r1
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    mp_aint%mpr(0) = mpwdsm6
    call mpinfr (ra%mpr, mp_aint%mpr, r1%mpr, mpnw)
    return
  end function

   function mp_anint (ra)
    implicit none
    type (mp_realm):: mp_anint
    type (mp_realm), intent (in):: ra
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_anint%mpr(0) = mpwdsm6
    call mpnint (ra%mpr, mp_anint%mpr, mpnw)
    return
  end function

   function mp_asin (ra)
    implicit none
    type (mp_realm):: mp_asin
    type (mp_realm), intent (in):: ra
    type (mp_realm) r1, r2, r3
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    r2%mpr(0) = mpwdsm6
    r3%mpr(0) = mpwdsm6
    call mpmul (ra%mpr, ra%mpr, r1%mpr, mpnw)
    call mpdmc (1.d0, 0, r2%mpr, mpnw)
    call mpsub (r2%mpr, r1%mpr, r3%mpr, mpnw)
    if (r3%mpr(2) < 0) then
      write (mpldb, 1)
1     format ('*** MP_ASIN: argument is not in (-1, 1).')
      call mpabrt ( 804)
    endif
    call mpsqrt (r3%mpr, r1%mpr, mpnw)
    mp_asin%mpr(0) = mpwdsm6
    call mpang (r1%mpr, ra%mpr, mp_asin%mpr, mpnw)
    return
  end function

  function mp_asinh (ra)
    implicit none
    type (mp_realm):: mp_asinh
    type (mp_realm), intent (in):: ra
    type (mp_realm) r1, r2, r3
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    r2%mpr(0) = mpwdsm6
    r3%mpr(0) = mpwdsm6
    mp_asinh%mpr(0) = mpwdsm6
    call mpmul (ra%mpr, ra%mpr, r1%mpr, mpnw)
    call mpdmc (1.d0, 0, r2%mpr, mpnw)
    call mpadd (r1%mpr, r2%mpr, r3%mpr, mpnw)
    call mpsqrt (r3%mpr, r1%mpr, mpnw)
    call mpadd (ra%mpr, r1%mpr, r2%mpr, mpnw)
    call mplog (r2%mpr, mp_asinh%mpr, mpnw)
    return
  end function

   function mp_atan (ra)
    implicit none
    type (mp_realm):: mp_atan
    type (mp_realm), intent (in):: ra
    type (mp_realm) r1
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    call mpdmc (1.d0, 0, r1%mpr, mpnw)
    mp_atan%mpr(0) = mpwdsm6
    call mpang (r1%mpr, ra%mpr, mp_atan%mpr, mpnw)
    return
  end function

   function mp_atan2 (ra, rb)
    implicit none
    type (mp_realm):: mp_atan2
    type (mp_realm), intent (in):: ra, rb
    integer mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwdsm)
    mp_atan2%mpr(0) = mpwdsm6
    call mpang (rb%mpr, ra%mpr, mp_atan2%mpr, mpnw)
    return
  end function

  function mp_atanh (ra)
    implicit none
    type (mp_realm):: mp_atanh
    type (mp_realm), intent (in):: ra
    type (mp_realm) r1, r2, r3
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    r2%mpr(0) = mpwdsm6
    r3%mpr(0) = mpwdsm6
    mp_atanh%mpr(0) = mpwdsm6
    call mpdmc (1.d0, 0, r1%mpr, mpnw)
    call mpadd (r1%mpr, ra%mpr, r2%mpr, mpnw)
    call mpsub (r1%mpr, ra%mpr, r3%mpr, mpnw)
    call mpdiv (r2%mpr, r3%mpr, r1%mpr, mpnw)
    if (r1%mpr(2) < 0) then
      write (mpldb, 1)
1     format ('*** MP_ATANH: argument is not in (-1, 1).')
      call mpabrt ( 805)
    endif
    call mplog (r1%mpr, r2%mpr, mpnw)
    call mpmuld (r2%mpr, 0.5d0, mp_atanh%mpr, mpnw)
    return
  end function

  function mp_ator1 (a, ib, iprec)
    implicit none
    type (mp_realm) mp_ator1
    integer, intent (in):: ib
    character(1), intent (in):: a(ib)
    integer mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    mp_ator1%mpr(0) = mpwdsm6
    call mpctomp (a, ib, mp_ator1%mpr, mpnw)
    return
  end function

  function mp_atorn (aa, iprec)
    implicit none
    character(*), intent (in):: aa
    type (mp_realm):: mp_atorn
    character(1) :: chr1(len(aa))
    integer i, l1, mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    l1 = len (aa)
    do i = 1, l1
      chr1(i) = aa(i:i)
    enddo
    mp_atorn%mpr(0) = mpwdsm6
    call mpctomp (chr1, l1, mp_atorn%mpr, mpnw)
    return
  end function

  subroutine mp_berne (nb, rb, iprec)
    implicit none
    integer, intent (in):: nb
    type (mp_realm), intent (out):: rb(nb)
    integer i, n1, mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    do i = 1, nb
      rb(i)%mpr(0) = mpwdsm6
    enddo
    n1 = mpwdsm
    call mpberner (n1, nb, rb(1)%mpr, mpnw)
    return
  end subroutine

  function mp_bessel_i (qa, ra)
    implicit none
    type (mp_realm):: mp_bessel_i
    type (mp_realm), intent (in):: qa, ra
    integer mpnw
    mpnw = min (int (qa%mpr(1)), int (ra%mpr(1)), mpwdsm)
    mp_bessel_i%mpr(0) = mpwdsm6
    call mpbesselir (qa%mpr, ra%mpr, mp_bessel_i%mpr, mpnw)
    return
  end function

  function mp_bessel_in (nu, ra)
    implicit none
    type (mp_realm):: mp_bessel_in
    integer, intent (in):: nu
    type (mp_realm), intent (in):: ra
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_bessel_in%mpr(0) = mpwdsm6
    call mpbesselinr (nu, ra%mpr, mp_bessel_in%mpr, mpnw)
    return
  end function

  function mp_bessel_j (qa, ra)
    implicit none
    type (mp_realm):: mp_bessel_j
    type (mp_realm), intent (in):: qa, ra
    integer mpnw
    mpnw = min (int (qa%mpr(1)), int (ra%mpr(1)), mpwdsm)
    mp_bessel_j%mpr(0) = mpwdsm6
    call mpbesseljr (qa%mpr, ra%mpr, mp_bessel_j%mpr, mpnw)
    return
  end function

  function mp_bessel_jn (nu, ra)
    implicit none
    type (mp_realm):: mp_bessel_jn
    integer, intent (in):: nu
    type (mp_realm), intent (in):: ra
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_bessel_jn%mpr(0) = mpwdsm6
    call mpbesseljnr (nu, ra%mpr, mp_bessel_jn%mpr, mpnw)
    return
  end function

  function mp_bessel_j0 (ra)
    implicit none
    type (mp_realm):: mp_bessel_j0
    type (mp_realm), intent (in):: ra
    integer mpnw, nu
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    nu = 0
    mp_bessel_j0%mpr(0) = mpwdsm6
    call mpbesseljnr (nu, ra%mpr, mp_bessel_j0%mpr, mpnw)
    return
  end function

  function mp_bessel_j1 (ra)
    implicit none
    type (mp_realm):: mp_bessel_j1
    type (mp_realm), intent (in):: ra
    integer mpnw, nu
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    nu = 1
    mp_bessel_j1%mpr(0) = mpwdsm6
    call mpbesseljnr (nu, ra%mpr, mp_bessel_j1%mpr, mpnw)
    return
  end function

  function mp_bessel_k (qa, ra)
    implicit none
    type (mp_realm):: mp_bessel_k
    type (mp_realm), intent (in):: qa, ra
    integer mpnw
    mpnw = min (int (qa%mpr(1)), int (ra%mpr(1)), mpwdsm)
    mp_bessel_k%mpr(0) = mpwdsm6
    call mpbesselkr (qa%mpr, ra%mpr, mp_bessel_k%mpr, mpnw)
    return
  end function

  function mp_bessel_kn (nu, ra)
    implicit none
    type (mp_realm):: mp_bessel_kn
    integer, intent (in):: nu
    type (mp_realm), intent (in):: ra
    integer mpnw
    mpnw = int (ra%mpr(1))
    mpnw = min (mpnw, mpwdsm)
    mp_bessel_kn%mpr(0) = mpwdsm6
    call mpbesselknr (nu, ra%mpr, mp_bessel_kn%mpr, mpnw)
    return
  end function

  function mp_bessel_y (qa, ra)
    implicit none
    type (mp_realm):: mp_bessel_y
    type (mp_realm), intent (in):: qa, ra
    integer mpnw
    mpnw = min (int (qa%mpr(1)), int (ra%mpr(1)), mpwdsm)
    mp_bessel_y%mpr(0) = mpwdsm6
    call mpbesselyr (qa%mpr, ra%mpr, mp_bessel_y%mpr, mpnw)
    return
  end function

  function mp_bessel_yn (nu, ra)
    implicit none
    type (mp_realm):: mp_bessel_yn
    integer, intent (in):: nu
    type (mp_realm), intent (in):: ra
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_bessel_yn%mpr(0) = mpwdsm6
    call mpbesselynr (nu, ra%mpr, mp_bessel_yn%mpr, mpnw)
    return
  end function

  function mp_bessel_y0 (ra)
    implicit none
    type (mp_realm):: mp_bessel_y0
    type (mp_realm), intent (in):: ra
    integer mpnw, nu
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    nu = 0
    mp_bessel_y0%mpr(0) = mpwdsm6
    call mpbesselynr (nu, ra%mpr, mp_bessel_y0%mpr, mpnw)
    return
  end function

  function mp_bessel_y1 (nu, ra)
    implicit none
    type (mp_realm):: mp_bessel_y1
    type (mp_realm), intent (in):: ra
    integer mpnw, nu
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    nu = 1
    mp_bessel_y1%mpr(0) = mpwdsm6
    call mpbesselynr (nu, ra%mpr, mp_bessel_y1%mpr, mpnw)
    return
  end function

  subroutine mp_binmd (ra, db, ic)
    implicit none
    type (mp_realm), intent (in):: ra
    real (mprknd), intent (out):: db
    integer, intent (out):: ic
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    call mpmdc (ra%mpr, db, ic, mpnw)
    return
  end subroutine

  function mp_ccos (za)
    implicit none
    type (mp_complexm):: mp_ccos
    type (mp_complexm), intent (in):: za
    integer l1, l2, l3, mpnw
    type (mp_complexm) z1, z2, z3
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    z1%mpc(0) = mpwdsm6
    z1%mpc(l2) = mpwdsm6
    z2%mpc(0) = mpwdsm6
    z2%mpc(l2) = mpwdsm6
    z3%mpc(0) = mpwdsm6
    z3%mpc(l2) = mpwdsm6
    l3 = mpwdsm6
    mp_ccos%mpc(0) = mpwdsm6
    mp_ccos%mpc(l3) = mpwdsm6
    call mpdmc (1.d0, 0, z1%mpc, mpnw)
    call mpdmc (0.d0, 0, z1%mpc(l2:), mpnw)
    call mpeq (za%mpc, z3%mpc(l2:), mpnw)
    call mpeq (za%mpc(l1:), z3%mpc, mpnw)
    z3%mpc(2) = - z3%mpc(2)
    call mpcexp (z3%mpc, z2%mpc, mpnw)
    call mpcdiv (z1%mpc, z2%mpc, z3%mpc, mpnw)
    call mpcadd (z2%mpc, z3%mpc, z1%mpc, mpnw)
    call mpmuld (z1%mpc, 0.5d0, mp_ccos%mpc, mpnw)
    call mpmuld (z1%mpc(l2:), 0.5d0, mp_ccos%mpc(l3:), mpnw)
    return
  end function

  function mp_cexp (za)
    implicit none
    type (mp_complexm):: mp_cexp
    type (mp_complexm), intent (in):: za
    integer l1, l2, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    mp_cexp%mpc(0) = mpwdsm6
    mp_cexp%mpc(l2) = mpwdsm6
    call mpcexp (za%mpc, mp_cexp%mpc, mpnw)
    return
  end function

  function mp_clog (za)
    implicit none
    type (mp_complexm):: mp_clog
    type (mp_complexm), intent (in):: za
    integer l1, l2, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    mp_clog%mpc(0) = mpwdsm6
    mp_clog%mpc(l2) = mpwdsm6
    call mpclog (za%mpc, mp_clog%mpc, mpnw)
    return
  end function

  function mp_conjg (za)
    implicit none
    type (mp_complexm):: mp_conjg
    type (mp_complexm), intent (in):: za
    integer l1, l2, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    mp_conjg%mpc(0) = mpwdsm6
    mp_conjg%mpc(l2) = mpwdsm6
    call mpconjg (za%mpc, mp_conjg%mpc, mpnw)
    return
  end function

  function mp_cos (ra)
    implicit none
    type (mp_realm):: mp_cos
    type (mp_realm), intent (in):: ra
    integer mpnw
    type (mp_realm) r1
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_cos%mpr(0) = mpwdsm6
    r1%mpr(0) = mpwdsm6
    call mpcssnr (ra%mpr, mp_cos%mpr, r1%mpr, mpnw)
    return
  end function

  function mp_cosh (ra)
    implicit none
    type (mp_realm):: mp_cosh
    type (mp_realm), intent (in):: ra
    integer mpnw
    type (mp_realm) r1
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_cosh%mpr(0) = mpwdsm6
    r1%mpr(0) = mpwdsm6
    call mpcsshr (ra%mpr, mp_cosh%mpr, r1%mpr, mpnw)
    return
  end function

  function mp_csin (za)
    implicit none
    type (mp_complexm):: mp_csin
    type (mp_complexm), intent (in):: za
    integer l1, l2, l3, mpnw
    type (mp_complexm) z1, z2, z3
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    z1%mpc(0) = mpwdsm6
    z1%mpc(l2) = mpwdsm6
    z2%mpc(0) = mpwdsm6
    z2%mpc(l2) = mpwdsm6
    z3%mpc(0) = mpwdsm6
    z3%mpc(l2) = mpwdsm6
    l3 = mpwdsm6
    mp_csin%mpc(0) = mpwdsm6
    mp_csin%mpc(l3) = mpwdsm6
    call mpdmc (1.d0, 0, z1%mpc, mpnw)
    call mpdmc (0.d0, 0, z1%mpc(l2:), mpnw)
    call mpeq (za%mpc, z3%mpc(l2:), mpnw)
    call mpeq (za%mpc(l1:), z3%mpc, mpnw)
    z3%mpc(2) = - z3%mpc(2)
    call mpcexp (z3%mpc, z2%mpc, mpnw)
    call mpcdiv (z1%mpc, z2%mpc, z3%mpc, mpnw)
    call mpcsub (z2%mpc, z3%mpc, z1%mpc, mpnw)
    call mpmuld (z1%mpc, 0.5d0, mp_csin%mpc(l3:), mpnw)
    call mpmuld (z1%mpc(l2:), 0.5d0, mp_csin%mpc, mpnw)
    mp_csin%mpc(l3+2) = - mp_csin%mpc(l3+2)
    return
  end function

  function mp_csqrt (za)
    implicit none
    type (mp_complexm):: mp_csqrt
    type (mp_complexm), intent (in):: za
    integer l1, l2, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    l2 = mpwdsm6
    mp_csqrt%mpc(0) = mpwdsm6
    mp_csqrt%mpc(l2) = mpwdsm6
    call mpcsqrt (za%mpc, mp_csqrt%mpc, mpnw)
    return
  end function

  subroutine mp_cssh (ra, rb, rc)
    implicit none
    type (mp_realm), intent (in):: ra
    type (mp_realm), intent (out):: rb, rc
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    rb%mpr(0) = mpwdsm6
    rc%mpr(0) = mpwdsm6
    call mpcsshr (ra%mpr, rb%mpr, rc%mpr, mpnw)
    return
  end subroutine

  subroutine mp_cssn (ra, rb, rc)
    implicit none
    type (mp_realm), intent (in):: ra
    type (mp_realm), intent (out):: rb, rc
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    rb%mpr(0) = mpwdsm6
    rc%mpr(0) = mpwdsm6
    call mpcssnr (ra%mpr, rb%mpr, rc%mpr, mpnw)
    return
  end subroutine

  function mp_dctoz (dca, iprec)
    implicit none
    type (mp_complexm):: mp_dctoz
    complex (kind(0.d0)), intent (in):: dca
    integer l1, mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    l1 = mpwdsm6
    mp_dctoz%mpc(0) = mpwdsm6
    mp_dctoz%mpc(l1) = mpwdsm6
    call mpdmc40 (dble (dca), 0, mp_dctoz%mpc, mpnw)
    call mpdmc40 (aimag (dca), 0, mp_dctoz%mpc(l1:), mpnw)
    return
  end function

  function mp_dctoz2 (dca, iprec)
    implicit none
    type (mp_complexm):: mp_dctoz2
    complex (kind(0.d0)), intent (in):: dca
    integer l1, mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    l1 = mpwdsm6
    mp_dctoz2%mpc(0) = mpwdsm6
    mp_dctoz2%mpc(l1) = mpwdsm6
    call mpdmc (dble (dca), 0, mp_dctoz2%mpc, mpnw)
    call mpdmc (aimag (dca), 0, mp_dctoz2%mpc(l1:), mpnw)
    return
  end function

  subroutine mp_decmd (ra, db, ib)
    implicit none
    type (mp_realm), intent (in):: ra
    real (mprknd), intent (out):: db
    integer, intent (out):: ib
    real, parameter:: alg102 = 0.301029995663981195d0
    real (mprknd) dt1, dt2
    integer i1, mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    call mpmdc (ra%mpr, dt1, i1, mpnw)
    if (dt1 /= 0.d0) then
      dt2 = alg102 * i1 + log10 (abs (dt1))
      ib = dt2
      if (dt2 < 0.d0) ib = ib - 1
      db = sign (10.d0 ** (dt2 - ib), dt1)
    else
      db = 0.d0
      ib = 0
    endif
  end subroutine

  function mp_digamma_be (nb, rb, rc)
    implicit none
    integer, intent (in):: nb
    type (mp_realm):: mp_digamma_be
    type (mp_realm), intent (in):: rb(nb), rc
    integer n1, mpnw
    mpnw = max (int (rb(1)%mpr(1)), int (rc%mpr(1)))
    mpnw = min (mpnw, mpwdsm)
    mp_digamma_be%mpr(0) = mpwdsm6
    n1 = mpwdsm
    call mpdigammabe (n1, nb, rb(1)%mpr, rc%mpr, mp_digamma_be%mpr, mpnw)
    return
  end function

  function mp_dtor (da, iprec)
    implicit none
    type (mp_realm):: mp_dtor
    real (mprknd), intent (in):: da
    integer i1, mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    mp_dtor%mpr(0) = mpwdsm6
    i1 = 0
    call mpdmc40 (da, i1, mp_dtor%mpr, mpnw)
    return
  end function

  function mp_dtor2 (da, iprec)
    implicit none
    type (mp_realm):: mp_dtor2
    real (mprknd), intent (in):: da
    integer i1, mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    mp_dtor2%mpr(0) = mpwdsm6
    i1 = 0
    call mpdmc (da, i1, mp_dtor2%mpr, mpnw)
    return
  end function

  subroutine mp_eform (ra, nb, nd, b)
    implicit none
    type (mp_realm), intent (in):: ra
    integer, intent (in):: nb, nd
    character(1), intent (out):: b(nb)
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    call mpeformat (ra%mpr, nb, nd, b, mpnw)
    return
  end subroutine

  function mp_egamma (iprec)
    implicit none
    type (mp_realm):: mp_egamma
    integer mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    mp_egamma%mpr(0) = mpwdsm6
    if (mpwprecr (mpegammacon) < mpnw) then
      write (mpldb, 1) mpnw
1     format ('*** MP_EGAMMA: Egamma must be precomputed to precision',i9,' words'/ &
      'by calling mpinit. See documentation for details.')
      call mpabrt ( 806)
    endif
    call mpeq (mpegammacon, mp_egamma%mpr, mpnw)

  end function

  function mp_erf (ra)
    implicit none
    type (mp_realm):: mp_erf
    type (mp_realm), intent (in):: ra
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_erf%mpr(0) = mpwdsm6
    call mperfr (ra%mpr, mp_erf%mpr, mpnw)
    return
  end function

  function mp_erfc (ra)
    implicit none
    type (mp_realm):: mp_erfc
    type (mp_realm), intent (in):: ra
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_erfc%mpr(0) = mpwdsm6
    call mperfcr (ra%mpr, mp_erfc%mpr, mpnw)
    return
  end function

  function mp_exp (ra)
    implicit none
    type (mp_realm):: mp_exp
    type (mp_realm), intent (in):: ra
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_exp%mpr(0) = mpwdsm6
    call mpexp (ra%mpr, mp_exp%mpr, mpnw)
    return
  end function

  function mp_expint (ra)
    implicit none
    type (mp_realm):: mp_expint
    type (mp_realm), intent (in):: ra
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_expint%mpr(0) = mpwdsm6
    call mpexpint (ra%mpr, mp_expint%mpr, mpnw)
    return
  end function

  subroutine mp_fform (ra, nb, nd, b)
    implicit none
    type (mp_realm), intent (in):: ra
    integer, intent (in):: nb, nd
    character(1), intent (out):: b(nb)
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    call mpfformat (ra%mpr, nb, nd, b, mpnw)
    return
  end subroutine

  function mp_gamma (ra)
    implicit none
    type (mp_realm):: mp_gamma
    type (mp_realm), intent (in):: ra
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_gamma%mpr(0) = mpwdsm6
    call mpgammar (ra%mpr, mp_gamma%mpr, mpnw)
    return
  end function

  function mp_hurwitz_zetan (ia, rb)
    implicit none
    type (mp_realm):: mp_hurwitz_zetan
    integer, intent (in):: ia
    type (mp_realm), intent (in):: rb
    integer mpnw
    mpnw = min (int (rb%mpr(1)), mpwdsm)
    mp_hurwitz_zetan%mpr(0) = mpwdsm6
    call mphurwitzzetan (ia, rb%mpr, mp_hurwitz_zetan%mpr, mpnw)
    return
  end function

  function mp_hurwitz_zetan_be (nb, rb, is, aa)
    implicit none
    type (mp_realm):: mp_hurwitz_zetan_be
    integer, intent (in):: nb, is
    type (mp_realm), intent (in):: rb(nb), aa
    integer mpnw, n1
    mpnw = min (int (aa%mpr(1)), mpwdsm)
    mp_hurwitz_zetan_be%mpr(0) = mpwdsm6
    n1 = mpwdsm
    call mphurwitzzetanbe (n1, nb, rb(1)%mpr, is, aa%mpr, &
      mp_hurwitz_zetan_be%mpr, mpnw)
    return
  end function

  function mp_hypergeom_pfq (np, nq, aa, bb, xx)
    implicit none
    type (mp_realm):: mp_hypergeom_pfq
    integer, intent (in):: np, nq
    type (mp_realm), intent (in):: aa(np), bb(nq), xx
    integer mpnw
    mpnw = min (int (xx%mpr(1)), mpwdsm)
    mp_hypergeom_pfq%mpr(0) = mpwdsm6
    call mphypergeompfq (mpwdsm, np, nq, aa(1)%mpr(0), bb(1)%mpr(0), &
      xx%mpr, mp_hypergeom_pfq%mpr, mpnw)
    return
  end function

  function mp_hypot (ra, rb)
    implicit none
    type (mp_realm):: mp_hypot
    type (mp_realm), intent (in):: ra, rb
    type (mp_realm) r1, r2, r3
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_hypot%mpr(0) = mpwdsm6
    r1%mpr(0) = mpwdsm6
    r2%mpr(0) = mpwdsm6
    r3%mpr(0) = mpwdsm6
    call mpmul (ra%mpr, ra%mpr, r1%mpr, mpnw)
    call mpmul (rb%mpr, rb%mpr, r2%mpr, mpnw)
    call mpadd (r1%mpr, r2%mpr, r3%mpr, mpnw)
    call mpsqrt (r3%mpr, mp_hypot%mpr, mpnw)
    return
  end function

  function mp_incgamma (ra, rb)
    implicit none
    type (mp_realm):: mp_incgamma
    type (mp_realm), intent (in):: ra, rb
    integer mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwdsm)
    mp_incgamma%mpr(0) = mpwdsm6
    call mpincgammar (ra%mpr, rb%mpr, mp_incgamma%mpr, mpnw)
    return
  end function

  function mp_log (ra)
    implicit none
    type (mp_realm):: mp_log
    type (mp_realm), intent (in):: ra
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_log%mpr(0) = mpwdsm6
    call mplog (ra%mpr, mp_log%mpr, mpnw)
    return
  end function

  function mp_log10 (ra)
    implicit none
    type (mp_realm):: mp_log10
    type (mp_realm), intent (in):: ra
    type (mp_realm) t1, t2
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    t1%mpr(0) = mpwdsm6
    t2%mpr(0) = mpwdsm6
    mp_log10%mpr(0) = mpwdsm6
    call mpdmc (10.d0, 0, t1%mpr, mpnw)
    call mplog (t1%mpr, t2%mpr, mpnw)
    call mplog (ra%mpr, t1%mpr, mpnw)
    call mpdiv (t1%mpr, t2%mpr, mp_log10%mpr, mpnw)
    return
  end function

  function mp_log2 (iprec)
    implicit none
    type (mp_realm):: mp_log2
    integer mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    mp_log2%mpr(0) = mpwdsm6
    if (mpwprecr (mplog2con) < mpnw) then
      write (mpldb, 1) mpnw
1     format ('*** MP_LOG2: Log(2) must be precomputed to precision',i9,' words'/ &
      'by calling mpinit. See documentation for details.')
      call mpabrt ( 807)
    endif
    call mpeq (mplog2con, mp_log2%mpr, mpnw)
  end function

  function mp_max (ra, rb, rc)
    implicit none
    type (mp_realm):: mp_max
    type (mp_realm), intent (in):: ra, rb
    type (mp_realm), intent (in), optional:: rc
    integer ic, mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    if (present (rc)) mpnw = max (mpnw, int (rc%mpr(1)))
    mpnw = min (mpnw, mpwdsm)
    mp_max%mpr(0) = mpwdsm6
    call mpcpr (ra%mpr, rb%mpr, ic, mpnw)
    if (ic >= 0) then
      call mpeq (ra%mpr, mp_max%mpr, mpnw)
    else
      call mpeq (rb%mpr, mp_max%mpr, mpnw)
    endif
    if (present (rc)) then
      call mpcpr (rc%mpr, mp_max%mpr, ic, mpnw)
      if (ic >= 0) call mpeq (rc%mpr, mp_max%mpr, mpnw)
    endif
    return
  end function

  function mp_min (ra, rb, rc)
    implicit none
    type (mp_realm):: mp_min
    type (mp_realm), intent (in):: ra, rb
    type (mp_realm), intent (in), optional:: rc
    integer ic, mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    if (present (rc)) mpnw = max (mpnw, int (rc%mpr(1)))
    mpnw = min (mpnw, mpwdsm)
    mp_min%mpr(0) = mpwdsm6
    call mpcpr (ra%mpr, rb%mpr, ic, mpnw)
    if (ic <= 0) then
      call mpeq (ra%mpr, mp_min%mpr, mpnw)
    else
      call mpeq (rb%mpr, mp_min%mpr, mpnw)
    endif
    if (present (rc)) then
      call mpcpr (rc%mpr, mp_min%mpr, ic, mpnw)
      if (ic <= 0) call mpeq (rc%mpr, mp_min%mpr, mpnw)
    endif
    return
  end function

  function mp_mod (ra, rb)
    implicit none
    type (mp_realm):: mp_mod
    type (mp_realm), intent (in):: ra, rb
    type (mp_realm) r1, r2, r3
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_mod%mpr(0) = mpwdsm6
    r1%mpr(0) = mpwdsm6
    r2%mpr(0) = mpwdsm6
    r3%mpr(0) = mpwdsm6
    call mpdiv (ra%mpr, rb%mpr, r1%mpr, mpnw)
    call mpinfr (r1%mpr, r2%mpr, r3%mpr, mpnw)
    call mpmul (rb%mpr, r2%mpr, r3%mpr, mpnw)
    call mpsub (ra%mpr, r3%mpr, mp_mod%mpr, mpnw)
    return
  end function

  function mp_nrt (ra, ib)
    implicit none
    type (mp_realm):: mp_nrt
    type (mp_realm), intent (in):: ra
    integer, intent (in):: ib
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_nrt%mpr(0) = mpwdsm6
    call mpnrtr (ra%mpr, ib, mp_nrt%mpr, mpnw)
    return
  end function

  function mp_pi (iprec)
    implicit none
    type (mp_realm):: mp_pi
    integer mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    mp_pi%mpr(0) = mpwdsm6
    if (mpwprecr (mppicon) < mpnw) then
      write (mpldb, 1) mpnw
1     format ('*** MP_PI: Pi must be precomputed to precision',i9,' words'/ &
      'by calling mpinit. See documentation for details.')
      call mpabrt ( 808)
    endif
    call mpeq (mppicon, mp_pi%mpr, mpnw)
  end function

  function mp_polygamma (nn, ra)
    implicit none
    integer, intent (in):: nn
    type (mp_realm), intent (in):: ra
    type (mp_realm) mp_polygamma
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_polygamma%mpr(0) = mpwdsm6
    call mppolygamma (nn, ra%mpr, mp_polygamma%mpr, mpnw)
    return
  end function

  function mp_polygamma_be (nb, rb, nn, ra)
    implicit none
    integer, intent (in):: nb, nn
    type (mp_realm), intent (in):: ra, rb(nb)
    type (mp_realm) mp_polygamma_be
    integer n1, mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_polygamma_be%mpr(0) = mpwdsm6
    n1 = mpwdsm
    call mppolygammabe (n1, nb, rb(1)%mpr, nn, ra%mpr, mp_polygamma_be%mpr, mpnw)
    return
  end function

  subroutine mp_polylog_inix (nn, arr, iprec)
    implicit none
    integer, intent (in):: nn
    type (mp_realm), intent (out):: arr(abs(nn))
    integer i, mpnw, n1

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>

    do i = 1, abs (nn)
      arr(i)%mpr(0) = mpwdsm6
    enddo

    n1 = mpwdsm
    call mppolylogini (n1, nn, arr(1)%mpr, mpnw)
    return
  end subroutine

  function mp_polylog_negx (nn, arr, ra)
    implicit none
    integer, intent (in):: nn
    type (mp_realm), intent (in):: arr(abs(nn))
    type (mp_realm), intent (in):: ra
    type (mp_realm) mp_polylog_negx
    integer mpnw, n1
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_polylog_negx%mpr(0) = mpwdsm6
    n1 = mpwdsm
    call mppolylogneg (n1, nn, arr(1)%mpr, ra%mpr, mp_polylog_negx%mpr, mpnw)
    return
  end function

  function mp_polylog_pos (nn, ra)
    implicit none
    integer, intent (in):: nn
    type (mp_realm), intent (in):: ra
    type (mp_realm) mp_polylog_pos
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_polylog_pos%mpr(0) = mpwdsm6
    call mppolylogpos (nn, ra%mpr, mp_polylog_pos%mpr, mpnw)
    return
  end function

  function mp_prodd (ra, db)
    implicit none
    type (mp_realm):: mp_prodd
    type (mp_realm), intent (in):: ra
    real (mprknd), intent (in):: db
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_prodd%mpr(0) = mpwdsm6
    call mpmuld (ra%mpr, db, mp_prodd%mpr, mpnw)
  end function

  function mp_prodq (ra, qb)
    implicit none
    type (mp_realm):: mp_prodq
    type (mp_realm), intent (in):: ra
    real (max (mprknd2, kind (1.0))), intent (in):: qb
    type (mp_realm) r1
    integer mpnw
    r1%mpr(0) = mpwdsm6
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_prodq%mpr(0) = mpwdsm6
    call mpqmc (qb, 0, r1%mpr, mpnw)
    call mpmul (ra%mpr, r1%mpr, mp_prodq%mpr, mpnw)
  end function

  function mp_qtor (qa, iprec)
    implicit none
    type (mp_realm):: mp_qtor
    integer mpnw
    real (max (mprknd2, kind (1.))), intent (in):: qa

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    mp_qtor%mpr(0) = mpwdsm6
    call mpqmc90 (qa, 0, mp_qtor%mpr, mpnw)
    return
  end function

  function mp_qtor2 (qa, iprec)
    implicit none
    type (mp_realm):: mp_qtor2
    integer mpnw
    real (max (mprknd2, kind (1.))), intent (in):: qa

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    mp_qtor2%mpr(0) = mpwdsm6
    call mpqmc (qa, 0, mp_qtor2%mpr, mpnw)
    return
  end function

  function mp_quotd (ra, db)
    implicit none
    type (mp_realm):: mp_quotd
    type (mp_realm), intent (in):: ra
    real (mprknd), intent (in):: db
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_quotd%mpr(0) = mpwdsm6
    call mpdivd (ra%mpr, db, mp_quotd%mpr, mpnw)
  end function

  function mp_quotq (ra, qb)
    implicit none
    type (mp_realm):: mp_quotq
    type (mp_realm), intent (in):: ra
    real (max (mprknd2, kind (1.0))), intent (in):: qb
    type (mp_realm) r1
    integer mpnw
    r1%mpr(0) = mpwdsm6
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_quotq%mpr(0) = mpwdsm6
    call mpqmc (qb, 0, r1%mpr, mpnw)
    call mpdiv (ra%mpr, r1%mpr, mp_quotq%mpr, mpnw)
  end function

  function mp_rand (ra)
    implicit none
    type (mp_realm):: mp_rand
    type (mp_realm), intent (in):: ra
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_rand%mpr(0) = mpwdsm6
    call mprandr (ra%mpr, mp_rand%mpr, mpnw)
    return
  end function

!   Five variations are necessary here due to Fortran rules about optional arguments.

  subroutine mp_readr1 (iu, r1, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_realm), intent (out):: r1
    integer mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    r1%mpr(0) = mpwdsm6
    call mpinp (iu, r1%mpr, mpnw)
    return
  end subroutine

  subroutine mp_readr2 (iu, r1, r2, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_realm), intent (out):: r1, r2
    integer mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    r1%mpr(0) = mpwdsm6
    r2%mpr(0) = mpwdsm6
    call mpinp (iu, r1%mpr, mpnw)
    call mpinp (iu, r2%mpr, mpnw)
    return
  end subroutine

  subroutine mp_readr3 (iu, r1, r2, r3, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_realm), intent (out):: r1, r2, r3
    integer mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    r1%mpr(0) = mpwdsm6
    r2%mpr(0) = mpwdsm6
    r3%mpr(0) = mpwdsm6
    call mpinp (iu, r1%mpr, mpnw)
    call mpinp (iu, r2%mpr, mpnw)
    call mpinp (iu, r3%mpr, mpnw)
    return
  end subroutine

  subroutine mp_readr4 (iu, r1, r2, r3, r4, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_realm), intent (out):: r1, r2, r3, r4
    integer mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    r1%mpr(0) = mpwdsm6
    r2%mpr(0) = mpwdsm6
    r3%mpr(0) = mpwdsm6
    r4%mpr(0) = mpwdsm6
    call mpinp (iu, r1%mpr, mpnw)
    call mpinp (iu, r2%mpr, mpnw)
    call mpinp (iu, r3%mpr, mpnw)
    call mpinp (iu, r4%mpr, mpnw)
    return
  end subroutine

  subroutine mp_readr5 (iu, r1, r2, r3, r4, r5, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_realm), intent (out):: r1, r2, r3, r4, r5
    integer mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    r1%mpr(0) = mpwdsm6
    r2%mpr(0) = mpwdsm6
    r3%mpr(0) = mpwdsm6
    r4%mpr(0) = mpwdsm6
    r5%mpr(0) = mpwdsm6
    call mpinp (iu, r1%mpr, mpnw)
    call mpinp (iu, r2%mpr, mpnw)
    call mpinp (iu, r3%mpr, mpnw)
    call mpinp (iu, r4%mpr, mpnw)
    call mpinp (iu, r5%mpr, mpnw)
    return
  end subroutine

  subroutine mp_readz1 (iu, z1, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_complexm), intent (out):: z1
    integer l1, mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    z1%mpc(0) = mpwdsm6
    l1 = mpwdsm6
    z1%mpc(l1) = mpwdsm6
    call mpinp (iu, z1%mpc, mpnw)
    call mpinp (iu, z1%mpc(l1:), mpnw)
    return
  end subroutine

  subroutine mp_readz2 (iu, z1, z2, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_complexm), intent (out):: z1, z2
    integer l1, l2, mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    z1%mpc(0) = mpwdsm6
    l1 = mpwdsm6
    z1%mpc(l1) = mpwdsm6
    z2%mpc(0) = mpwdsm6
    l2 = mpwdsm6
    z2%mpc(l2) = mpwdsm6
    call mpinp (iu, z1%mpc, mpnw)
    call mpinp (iu, z1%mpc(l1:), mpnw)
    call mpinp (iu, z2%mpc, mpnw)
    call mpinp (iu, z2%mpc(l2:), mpnw)
    return
  end subroutine

  subroutine mp_readz3 (iu, z1, z2, z3, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_complexm), intent (out):: z1, z2, z3
    integer l1, l2, l3, mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    z1%mpc(0) = mpwdsm6
    l1 = mpwdsm6
    z1%mpc(l1) = mpwdsm6
    z2%mpc(0) = mpwdsm6
    l2 = mpwdsm6
    z2%mpc(l2) = mpwdsm6
    z3%mpc(0) = mpwdsm6
    l3 = mpwdsm6
    z3%mpc(l3) = mpwdsm6
    call mpinp (iu, z1%mpc, mpnw)
    call mpinp (iu, z1%mpc(l1:), mpnw)
    call mpinp (iu, z2%mpc, mpnw)
    call mpinp (iu, z2%mpc(l2:), mpnw)
    call mpinp (iu, z3%mpc, mpnw)
    call mpinp (iu, z3%mpc(l3:), mpnw)
    return
  end subroutine

  subroutine mp_readz4 (iu, z1, z2, z3, z4, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_complexm), intent (out):: z1, z2, z3, z4
    integer l1, l2, l3, l4, mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    z1%mpc(0) = mpwdsm6
    l1 = mpwdsm6
    z1%mpc(l1) = mpwdsm6
    z2%mpc(0) = mpwdsm6
    l2 = mpwdsm6
    z2%mpc(l2) = mpwdsm6
    z3%mpc(0) = mpwdsm6
    l3 = mpwdsm6
    z3%mpc(l3) = mpwdsm6
    z4%mpc(0) = mpwdsm6
    l4 = mpwdsm6
    z4%mpc(l4) = mpwdsm6
    call mpinp (iu, z1%mpc, mpnw)
    call mpinp (iu, z1%mpc(l1:), mpnw)
    call mpinp (iu, z2%mpc, mpnw)
    call mpinp (iu, z2%mpc(l2:), mpnw)
    call mpinp (iu, z3%mpc, mpnw)
    call mpinp (iu, z3%mpc(l3:), mpnw)
    call mpinp (iu, z4%mpc, mpnw)
    call mpinp (iu, z4%mpc(l4:), mpnw)
    return
  end subroutine

  subroutine mp_readz5 (iu, z1, z2, z3, z4, z5, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_complexm), intent (out):: z1, z2, z3, z4, z5
    integer l1, l2, l3, l4, l5, mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    z1%mpc(0) = mpwdsm6
    l1 = mpwdsm6
    z1%mpc(l1) = mpwdsm6
    z2%mpc(0) = mpwdsm6
    l2 = mpwdsm6
    z2%mpc(l2) = mpwdsm6
    z3%mpc(0) = mpwdsm6
    l3 = mpwdsm6
    z3%mpc(l3) = mpwdsm6
    z4%mpc(0) = mpwdsm6
    l4 = mpwdsm6
    z4%mpc(l4) = mpwdsm6
    z5%mpc(0) = mpwdsm6
    l5 = mpwdsm6
    z5%mpc(l5) = mpwdsm6
    call mpinp (iu, z1%mpc, mpnw)
    call mpinp (iu, z1%mpc(l1:), mpnw)
    call mpinp (iu, z2%mpc, mpnw)
    call mpinp (iu, z2%mpc(l2:), mpnw)
    call mpinp (iu, z3%mpc, mpnw)
    call mpinp (iu, z3%mpc(l3:), mpnw)
    call mpinp (iu, z4%mpc, mpnw)
    call mpinp (iu, z4%mpc(l4:), mpnw)
    call mpinp (iu, z5%mpc, mpnw)
    call mpinp (iu, z5%mpc(l5:), mpnw)
    return
  end subroutine

  function mp_rtod (ra)
    implicit none
    real (mprknd):: mp_rtod
    type (mp_realm), intent (in):: ra
    integer n1, mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    call mpmdc (ra%mpr, mp_rtod, n1, mpnw)
    mp_rtod = mp_rtod * 2.d0**n1
    return
  end function

  function mp_rtom (ra, iprec)
    implicit none
    type (mp_realm):: mp_rtom
    type (mp_real), intent (in):: ra
    integer mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    mp_rtom%mpr(0) = mpwdsm6
    call mpeq (ra%mpr, mp_rtom%mpr, mpnw)
    return
  end function

  function mp_rtoq (ra)
    implicit none
    real (max (mprknd2, kind (1.))):: mp_rtoq, t1
    type (mp_realm), intent (in):: ra
    integer n1, mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    call mpmqc (ra%mpr, mp_rtoq, n1, mpnw)
    t1 = 2.d0
    mp_rtoq = mp_rtoq * t1**n1
    return
  end function

  function mp_rtor (ra, iprec)
    implicit none
    type (mp_realm):: mp_rtor
    type (mp_realm), intent (in):: ra
    integer mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    mp_rtor%mpr(0) = mpwdsm6
    call mpeq (ra%mpr, mp_rtor%mpr, mpnw)
    return
  end function

  function mp_rtoz (ra, rb, iprec)
    implicit none
    type (mp_complexm):: mp_rtoz
    type (mp_realm), intent (in):: ra, rb
    integer l1, mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    l1 = mpwdsm6
    mp_rtoz%mpc(0) = mpwdsm6
    mp_rtoz%mpc(l1) = mpwdsm6
    call mpeq (ra%mpr, mp_rtoz%mpc, mpnw)
    call mpeq (rb%mpr, mp_rtoz%mpc(l1:), mpnw)
    return
  end function

  function mp_sign (ra, rb)
    implicit none
    type (mp_realm):: mp_sign
    type (mp_realm), intent (in):: ra, rb
    integer mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwdsm)
    mp_sign%mpr(0) = mpwdsm6
    call mpeq (ra%mpr, mp_sign%mpr, mpnw)
    mp_sign%mpr(2) = sign (mp_sign%mpr(2), rb%mpr(2))
    return
  end function

  function mp_sin (ra)
    implicit none
    type (mp_realm):: mp_sin
    type (mp_realm), intent (in):: ra
    integer mpnw
    type (mp_realm) r1
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_sin%mpr(0) = mpwdsm6
    r1%mpr(0) = mpwdsm6
    call mpcssnr (ra%mpr, r1%mpr, mp_sin%mpr, mpnw)
    return
  end function

  function mp_sinh (ra)
    implicit none
    type (mp_realm):: mp_sinh
    type (mp_realm), intent (in):: ra
    integer mpnw
    type (mp_realm) r1
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_sinh%mpr(0) = mpwdsm6
    r1%mpr(0) = mpwdsm6
    call mpcsshr (ra%mpr, r1%mpr, mp_sinh%mpr, mpnw)
    return
  end function

  function mp_sqrt (ra)
    implicit none
    type (mp_realm):: mp_sqrt
    type (mp_realm), intent (in):: ra
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_sqrt%mpr(0) = mpwdsm6
    call mpsqrt (ra%mpr, mp_sqrt%mpr, mpnw)
    return
  end function

  function mp_struve_hn (nu, ra)
    implicit none
    integer, intent (in):: nu
    type (mp_realm):: mp_struve_hn
    type (mp_realm), intent (in):: ra
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_struve_hn%mpr(0) = mpwdsm6
    call mpstruvehn (nu, ra%mpr, mp_struve_hn%mpr, mpnw)
    return
  end function

  function mp_tan (ra)
    implicit none
    type (mp_realm):: mp_tan
    type (mp_realm), intent (in):: ra
    type (mp_realm) r1, r2
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    r2%mpr(0) = mpwdsm6
    mp_tan%mpr(0) = mpwdsm6
    call mpcssnr (ra%mpr, r1%mpr, r2%mpr, mpnw)
    if (r1%mpr(2) == 0) then
      write (mpldb, 1)
1     format ('*** MP_TAN: Cos of argument is zero.')
      call mpabrt ( 809)
    endif
    call mpdiv (r2%mpr, r1%mpr, mp_tan%mpr, mpnw)
    return
  end function

  function mp_tanh (ra)
    implicit none
    type (mp_realm):: mp_tanh
    type (mp_realm), intent (in):: ra
    type (mp_realm) r1, r2
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    r1%mpr(0) = mpwdsm6
    r2%mpr(0) = mpwdsm6
    mp_tanh%mpr(0) = mpwdsm6
    call mpcsshr (ra%mpr, r1%mpr, r2%mpr, mpnw)
    call mpdiv (r2%mpr, r1%mpr, mp_tanh%mpr, mpnw)
    return
  end function

  function mp_wprec (ra)
    implicit none
    integer mp_wprec
    type (mp_realm), intent (in):: ra
    mp_wprec = min (int (ra%mpr(1)), mpwdsm)
    return
  end function

  function mp_wprecz (za)
    implicit none
    integer mp_wprecz
    type (mp_complexm), intent (in):: za
    integer l1
    l1 = za%mpc(0)
    mp_wprecz = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mp_wprecz = min (mp_wprecz, mpwdsm)
    return
  end function

  subroutine mp_writer (iu, ln, ld, r1, r2, r3, r4, r5)
    implicit none
    integer, intent (in):: iu, ln, ld
    type (mp_realm), intent (in):: r1, r2, r3, r4, r5
    optional:: r2, r3, r4, r5
    integer mpnw

    mpnw = min (int (r1%mpr(1)), mpwdsm)
    call mpout (iu, ln, ld, r1%mpr, mpnw)
    if (present (r2)) then
      mpnw = min (int (r2%mpr(1)), mpwdsm)
      call mpout (iu, ln, ld, r2%mpr, mpnw)
    endif
    if (present (r3)) then
      mpnw = min (int (r3%mpr(1)), mpwdsm)
      call mpout (iu, ln, ld, r3%mpr, mpnw)
    endif
    if (present (r4)) then
      mpnw = min (int (r4%mpr(1)), mpwdsm)
      call mpout (iu, ln, ld, r4%mpr, mpnw)
    endif
    if (present (r5)) then
      mpnw = min (int (r5%mpr(1)), mpwdsm)
      call mpout (iu, ln, ld, r5%mpr, mpnw)
    endif

    return
  end subroutine

  subroutine mp_writez (iu, ln, ld, z1, z2, z3, z4, z5)
    implicit none
    integer, intent (in):: iu, ln, ld
    type (mp_complexm), intent (in):: z1, z2, z3, z4, z5
    optional:: z2, z3, z4, z5
    integer l1, l2, l3, l4, l5, mpnw

    l1 = z1%mpc(0)
    mpnw = min (int (z1%mpc(1)), mpwdsm)
    call mpout (iu, ln, ld, z1%mpc, mpnw)
    call mpout (iu, ln, ld, z1%mpc(l1:), mpnw)
    if (present (z2)) then
      l2 = z2%mpc(0)
      mpnw = min (int (z2%mpc(1)), mpwdsm)
      call mpout (iu, ln, ld, z2%mpc, mpnw)
      call mpout (iu, ln, ld, z2%mpc(l2:), mpnw)
    endif
    if (present (z3)) then
      l3 = z3%mpc(0)
      mpnw = min (int (z3%mpc(1)), mpwdsm)
      call mpout (iu, ln, ld, z3%mpc, mpnw)
      call mpout (iu, ln, ld, z3%mpc(l3:), mpnw)
    endif
    if (present (z4)) then
      l4 = z4%mpc(0)
      mpnw = min (int (z4%mpc(1)), mpwdsm)
      call mpout (iu, ln, ld, z4%mpc, mpnw)
      call mpout (iu, ln, ld, z4%mpc(l4:), mpnw)
    endif
    if (present (z5)) then
      l5 = z5%mpc(0)
      mpnw = min (int (z5%mpc(1)), mpwdsm)
      call mpout (iu, ln, ld, z5%mpc, mpnw)
      call mpout (iu, ln, ld, z5%mpc(l5:), mpnw)
    endif

    return
  end subroutine

  function mp_zeta (ra)
    implicit none
    type (mp_realm):: mp_zeta
    type (mp_realm), intent (in):: ra
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwdsm)
    mp_zeta%mpr(0) = mpwdsm6
    call mpzetar (ra%mpr, mp_zeta%mpr, mpnw)
    return
  end function

  function mp_zeta_be (nb, rb, rc)
    implicit none
    integer, intent (in):: nb
    type (mp_realm):: mp_zeta_be
    type (mp_realm), intent (in):: rb(nb), rc
    integer n1, mpnw
    mpnw = max (int (rb(1)%mpr(1)), int (rc%mpr(1)))
    mpnw = min (mpnw, mpwdsm)
    mp_zeta_be%mpr(0) = mpwdsm6
    n1 = mpwdsm
    call mpzetabe (n1, nb, rb(1)%mpr, rc%mpr, mp_zeta_be%mpr, mpnw)
    return
  end function

  function mp_zeta_int (ia, iprec)
    implicit none
    type (mp_realm):: mp_zeta_int
    integer, intent (in):: ia
    integer mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    mp_zeta_int%mpr(0) = mpwdsm6
    call mpzetaintr (ia, mp_zeta_int%mpr, mpnw)
    return
  end function

  function mp_ztodc (za)
    implicit none
    complex (kind(0.d0)):: mp_ztodc
    type (mp_complexm), intent (in):: za
    integer l1, mpnw, n1, n2
    real (mprknd) d1, d2
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    call mpmdc (za%mpc, d1, n1, mpnw)
    d1 = d1 * 2.d0 ** n1
    call mpmdc (za%mpc(l1:), d2, n2, mpnw)
    d2 = d2 * 2.d0 ** n2
    mp_ztodc = cmplx (d1, d2, mprknd)
    return
  end function

  function mp_ztor (za, iprec)
    implicit none
    type (mp_realm):: mp_ztor
    type (mp_complexm), intent (in):: za
    integer l1, mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwdsm)
    mp_ztor%mpr(0) = mpwdsm6
    call mpeq (za%mpc, mp_ztor%mpr, mpnw)
    return
  end function

  function mp_ztoz (za, iprec)
    implicit none
    type (mp_complexm):: mp_ztoz
    type (mp_complexm), intent (in):: za
    integer l2, mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    l2 = mpwdsm6
    mp_ztoz%mpc(0) = mpwdsm6
    mp_ztoz%mpc(l2) = mpwdsm6
    call mpceq (za%mpc, mp_ztoz%mpc, mpnw)
    return
  end function

  function mp_ztozm (za, iprec)
    implicit none
    type (mp_complexm):: mp_ztozm
    type (mp_complex), intent (in):: za
    integer l2, mpnw

!>  In variant #1, uncomment these lines:
    integer, optional, intent (in):: iprec
    if (present (iprec)) then
      mpnw = mp_setwp (iprec)
    else
      mpnw = mpwdsm
    endif
!  Otherwise in variant #2, uncomment these lines:
!    integer, intent (in):: iprec
!    mpnw = mp_setwp (iprec)
!>>
    l2 = mpwdsm6
    mp_ztozm%mpc(0) = mpwdsm6
    mp_ztozm%mpc(l2) = mpwdsm6
    call mpceq (za%mpc, mp_ztozm%mpc, mpnw)
    return
  end function
end module mpfunh

