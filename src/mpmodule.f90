!*****************************************************************************

!  MPFUN20-Fort: A thread-safe arbitrary precision package with special functions
!  Main module (module MPMODULE) -- references all other modules for user.

!  Revision date:  16 Mar 2023

!  AUTHOR:
!     David H. Bailey
!     Lawrence Berkeley National Lab (retired) and University of California, Davis
!     Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2022 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!  PURPOSE OF PACKAGE:
!    This package permits one to perform floating-point computations (real and
!    complex) to arbitrarily high numeric precision, by making only relatively
!    minor changes to existing Fortran-90 programs.  All basic arithmetic
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

!  DESCRIPTION OF THIS MODULE (MPMODULE):
!    This module links all lower-level modules and is the connection between
!    user codes and the lower modules. It also declares as private routines in
!    lower-level modules that are not intended to be called directly by the user.
!    See documentation for details. 

module mpmodule

use mpfuna
use mpfunb
use mpfunc
use mpfund
use mpfune
use mpfunf
use mpfung
use mpfunh

!   Private subroutine names in module MPFUNB:

private &
  mpabs, mpadd, mpcabs, mpcadd, mpcdiv, mpceq, mpcmul, mpcnpwr, mpconjg, &
  mpcsqrt, mpcsub, mpcpr, mpdiv, mpdivd, mpdivd40, mpdmc, mpdmc40, mpeq, &
  mpinfr, mpmdc, mpmul, mpmuld, mpmuld40, mpneg, mpnint, mpnorm, mpnpwr, &
  mpnrtr, mprandr, mproun, mpsgn, mpsqrt, mpsub, mpmqc,  mpqmc, mpqmc90, &
  mpfftcr, mpfftrc, mpfft1, mpfft2, mpfft3, mpinifft, mplconv, mpmulx

!   Private subroutine names in module MPFUNC:

private &
  mpctomp, mpeformat, mpfformat, mpinp, mpout

!   Private subroutine names in module MPFUND:

private &
  mpagmr, mpang, mpcagm, mpcexp,  mpclog, mpcpowcc, mpcpowcr, mpcpowrc, &
  mpcsshr, mpcssnr, mpegammaq, mpexp, mpinitran, mplog, mplog2q, mppiq, mppower

!   Private subroutine names in module MPFUNE:

private &
  mpinitwds, mpwprecr, mpspacer, mpberner, mppolyadd, mppolysub, mppolymul, &
  mpbesselinr, mpbesselir, mpbesseljnr, mpbesseljr, mpbesselknr, mpbesselkr, &
  mpbesselynr, mpbesselyr, mpdigammabe, mperfr, mperfcr, mpexpint, mpgammar, &
  mphurwitzzetan, mphurwitzzetanbe, mphypergeompfq, mpincgammar, mppolygamma, &
  mppolygammabe, mppolylogini, mppolylogneg, mppolylogpos, mpstruvehn, mpzetar, &
  mpzetaintr, mpzetabe

end module mpmodule

