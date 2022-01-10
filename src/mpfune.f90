!*****************************************************************************

!  MPFUN20-Fort: A thread-safe arbitrary precision computation package
!  Special functions module (module MPFUNE)

!  Revision date:  25 Dec 2021

!  AUTHOR:
!    David H. Bailey
!    Lawrence Berkeley National Lab (retired) and University of California, Davis
!    Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2021 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!  PURPOSE OF PACKAGE:
!    This package permits one to perform floating-point computations (real and
!    complex) to arbitrarily high numeric precision, by making only relatively
!    minor changes to existing Fortran-90 programs.  All basic arithmetic
!    operations and transcendental functions are supported, together with several
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
 
!  DESCRIPTION OF THIS MODULE (MPFUNE):
!    This module contains subroutines to perform special functions. Additional
!    functions will be added as they are completed.

!  NOTE ON PROGRAMMING CONVENTION FOR THIS MODULE:
!    This module is designed to facilitate easy translation (using the program
!    convfmp.f90) to MPFR calls, for use in the MPFUN-MPFR package. Thus all
!    routines in this module (following !>> below) adhere to this convention:
!    No routine may directly reference the contents of an MP array. Instead,
!    common operations on MP arrays, specific to MPFUN-Fort, are handled by the
!    first few routines below, following !> but prior to !>>.

module mpfune
use mpfuna
use mpfunb
use mpfunc
use mpfund

contains

!>
!  These routines perform simple operations on the MP data structure, specific
!  to MPFUN-Fort. When this module is translated for MPFUN-MPFR, these routines
!  must be modified or replaced by MPFR equivalents.

subroutine mpinitwds (ra, mpnw)
implicit none
integer (mpiknd) ra(0:)
integer mpnw
ra(0) = mpnw + 6
ra(1) = mpnw
ra(2) = 0
ra(3) = 0
ra(4) = 0
return
end subroutine

subroutine mpabs (ra, rb, mpnw)
implicit none
integer (mpiknd) ra(0:), rb(0:)
integer mpnw
call mpeq (ra, rb, mpnw)
rb(2) = min (abs (ra(2)), mpnw)
return
end subroutine

subroutine mpneg (ra, rb, mpnw)
implicit none
integer (mpiknd) ra(0:), rb(0:)
integer mpnw, na
call mpeq (ra, rb, mpnw)
na = min (abs (int (ra(2))), mpnw)
rb(2) = - sign (na, int (ra(2)))
return
end subroutine

function mpsigntr (ra)
implicit none
integer (mpiknd) ra(0:)
integer ia, mpsigntr
ia = ra(2)
if (ia == 0) then
  mpsigntr = 0
elseif (ia > 0) then
  mpsigntr = 1
else
  mpsigntr = -1
endif
return
end function

function mpwprecr (ra)
implicit none
integer (mpiknd) ra(0:)
integer mpwprecr
mpwprecr = ra(1)
return
end function

function mpspacer (ra)
implicit none
integer (mpiknd) ra(0:)
integer mpspacer
mpspacer = ra(0)
return
end function

!>>
!  The following routines compute special functions:


subroutine mpberner (nb1, nb2, berne, mpnw)

!   This returns the array berne, containing Bernoulli numbers indexed 2*k for
!   k = 1 to n, to mpnw words precision. This is done by first computing
!   zeta(2*k), based on the following known formulas:

!   coth (pi*x) = cosh (pi*x) / sinh (pi*x)

!            1      1 + (pi*x)^2/2! + (pi*x)^4/4! + ...
!        =  ---- * -------------------------------------
!           pi*x    1 + (pi*x)^2/3! + (pi*x)^4/5! + ...

!        = 1/(pi*x) * (1 + (pi*x)^2/3 - (pi*x)^4/45 + 2*(pi*x)^6/945 - ...)

!        = 2/(pi*x) * Sum_{k >= 1} (-1)^(k+1) * zeta(2*k) * x^{2*k}

!   The strategy is to calculate the coefficients of the series by polynomial
!   operations. Polynomial division is performed by computing the reciprocal
!   of the denominator polynomial, by a polynomial Newton iteration, as follows. 
!   Let N(x) be the polynomial approximation to the numerator series; let D(x) be
!   a polynomial approximation to the numerator numerator series; and let Q_k(x)
!   be polynomial approximations to R(x) = 1/D(x).  Then iterate: 

!   Q_{k+1} = Q_k(x) + [1 - D(x)*Q_k(x)]*Q_k(x)

!   In these iterations, both the degree of the polynomial Q_k(x) and the
!   precision level in words are initially set to 4. When convergence is 
!   achieved at this precision level, the degree is doubled, and iterations are
!   continued, etc., until the final desired degree is achieved. Then the
!   precision level is doubled and iterations are performed in a similar way,
!   until the final desired precision level is achieved. The reciprocal polynomial
!   R(x) produced by this process is then multiplied by the numerator polynomial
!   N(x) to yield an approximation to the quotient series. The even zeta values
!   are then the coefficients of this series, scaled according to the formula above.

!   Once the even integer zeta values have been computed in this way, the Bernoulli
!   numbers are computed via the formula (for n > 0):

!   B(2*n) = (-1)^(n-1) * 2 * (2*n)! * zeta(2*n) / (2*pi)^(2*n)

!   Note: The notation in the description above is not the same as in the code below.

implicit none
integer, intent(in):: nb1, nb2, mpnw
integer i, i1, ic1, j, kn, mpnw1, nwds, n, nn1
integer:: idb = 0
integer (mpiknd), intent(out):: berne(0:nb1+5,nb2) 
integer (mpiknd) c1(0:mpnw+6,0:nb2), cp2(0:mpnw+6), p1(0:mpnw+6,0:nb2), &
  p2(0:mpnw+6,0:nb2), q(0:mpnw+6,0:nb2), q1(0:mpnw+6), &
  r(0:mpnw+6,0:nb2), s(0:mpnw+6,0:nb2), t1(0:mpnw+6), t2(0:mpnw+6), &
  t3(0:mpnw+6), t4(0:mpnw+6), eps(0:mpnw+6)
real (mprknd) dd1, dd2, dd3
real (mprknd):: alog102 = 0.30102999566398119d0

!  End of declaration.

i1 = 1000000000

do i = 1, nb2
  i1 = min (i1, mpspacer (berne(0:nb1+5,i)))
enddo

if (mpnw < 4 .or. i1 < mpnw + 6) then
  write (mpldb, 2)
2 format ('*** MPBERNER: uninitialized or inadequately sized arrays')
  call mpabrt (62)
endif

n = nb2
mpnw1 = mpnw + 1
nwds = mpnw1

if (idb > 0) write (6, 1) n, mpnw
1 format ('Even Bernoulli number calculation; n, mpnw =',2i6)

call mpinitwds (cp2, nwds)
call mpinitwds (q1, nwds)
call mpinitwds (t1, nwds)
call mpinitwds (t2, nwds)
call mpinitwds (t3, nwds)
call mpinitwds (t4, nwds)
call mpinitwds (eps, nwds)

do i = 0, nb2
  call mpinitwds (c1(0:mpnw1+5,i), nwds)
  call mpinitwds (p1(0:mpnw1+5,i), nwds)
  call mpinitwds (p2(0:mpnw1+5,i), nwds)
  call mpinitwds (q(0:mpnw1+5,i), nwds)
  call mpinitwds (r(0:mpnw1+5,i), nwds)
  call mpinitwds (s(0:mpnw1+5,i), nwds)
enddo

! cp2 = mppi (nwds) ** 2
! c1(0) = mpreal (1.d0, nwds)
! p1(0) = mpreal (1.d0, nwds)
! p2(0) = mpreal (1.d0, nwds)
! q(0) = mpreal (1.d0, nwds)

call mpmul (mppicon, mppicon, cp2, nwds)
call mpdmc (1.d0, 0, c1(0:mpnw1+5,0), nwds)
call mpdmc (1.d0, 0, p1(0:mpnw1+5,0), nwds)
call mpdmc (1.d0, 0, p2(0:mpnw1+5,0), nwds)
call mpdmc (1.d0, 0, q(0:mpnw1+5,0), nwds)

!   Construct numerator and denominator polynomials.

! do i = 1, n
!   c1(i) = mpreal (0.d0, nwds)
!   dd1 = 2.d0 * (i + 1) - 3.d0
!   dd2 = dd1 + 1.d0
!   dd3 = dd2 + 1.d0
!   p1(i) = cp2 * p1(i-1) / (dd1 * dd2)
!   p2(i) = cp2 * p2(i-1) / (dd2 * dd3)
!   q(i) = mpreal (0.d0, nwds)
! enddo

do i = 1, n
  call mpdmc (0.d0, 0, c1(0:mpnw1+5,i), nwds)
  dd1 = 2.d0 * (i + 1) - 3.d0
  dd2 = dd1 + 1.d0
  dd3 = dd2 + 1.d0
  call mpmul (cp2, p1(0:mpnw1+5,i-1), t1, nwds)
  call mpdivd (t1, dd1 * dd2, p1(0:mpnw1+5,i), nwds)
  call mpmul (cp2, p2(0:mpnw1+5,i-1), t1, nwds)
  call mpdivd (t1, dd2 * dd3, p2(0:mpnw1+5,i), nwds)
  call mpdmc (0.d0, 0, q(0:mpnw1+5,i), nwds)
enddo

kn = 4
nwds = 4

! eps = mpreal (2.d0, nwds) ** (70 - nwds * mpnbt)
! call mpdecmd (eps, dd1, nn1)

call mpdmc (2.d0, 0, t1, nwds)
call mpnpwr (t1, 70 - nwds * mpnbt, eps, nwds)
if (idb > 0) then
!   call mpdecmd (eps, dd1, nn1)
  call mpmdc (eps, dd1, nn1, nwds)
  write (6, 4) nwds, nint (alog102*nn1)
4 format ('nwds, log10eps =',2i6)
endif

! q1 = mpreal (0.d0, nwds)

call mpdmc (0.d0, 0, q1, nwds)

!   Perform Newton iterations with dynamic precision levels, using an
!   iteration formula similar to that used to evaluate reciprocals. 

do j = 1, 10000
  if (idb > 0) write (6, 5) j, kn, nwds
5 format ('j, kn, nwds =',3i6)

  call mppolymul (mpnw1, kn, p2, q, r, nwds)
  call mppolysub (mpnw1, kn, c1, r, s, nwds)
  call mppolymul (mpnw1, kn, s, q, r, nwds)
  call mppolyadd (mpnw1, kn, q, r, q, nwds)

!  t1 = q(kn) - q1

  call mpsub (q(0:mpnw1+5,kn), q1, t1, nwds)

  if (idb > 0) then
!     call mpdecmd (t1, dd1, nn1)
    call mpmdc (t1, dd1, nn1, nwds)
    if (dd1 .eq. 0.d0) then
      write (6, 6)
6     format ('Newton error = 0')
    else
      write (6, 7) nint (alog102*nn1) 
7     format ('Newton error = 10^',i6)
    endif
  endif
 
!   if (abs (t1) < eps) then

  call mpabs (t1, t2, nwds)
  call mpcpr (t2, eps, ic1, nwds)
  if (ic1 < 0) then
    if (kn == n .and. nwds == mpnw1) goto 100
    if (kn < n) then
      kn = min (2 * kn, n)
!      q1 = mpreal (0.d0, nwds)
      call mpdmc (0.d0, 0, q1, nwds)
    elseif (nwds < mpnw1) then
      nwds = min (2 * nwds, mpnw1)
!      eps = mpreal (2.d0, nwds) ** (70 - nwds * mpnbt)
!      q1 = mpreal (0.d0, nwds)
!      call mpdecmd (eps, dd1, nn1)

      call mpdmc (2.d0, 0, t1, nwds)
      call mpnpwr (t1, 70 - nwds * mpnbt, eps, nwds)
      call mpdmc (0.d0, 0, q1, nwds)
      if (idb > 0) then
!         call mpdecmd (eps, dd1, nn1)
        call mpmdc (eps, dd1, nn1, nwds)
        write (6, 4) nwds, nint (alog102*nn1)
      endif
    endif
  else
!    q1 = q(kn)
    call mpeq (q(0:mpnw1+5,kn), q1, nwds)
  endif
enddo

write (6, 8)
8 format ('MPBERNER: *** End loop error')
call mpabrt (99)

100 continue

if (idb > 0) write (6, 9)
9 format ('Even zeta computation complete')

!   Multiply numerator polynomial by reciprocal of denominator polynomial.

call mppolymul (mpnw1, n, p1, q, r, nwds)

!  If even zetas are needed, they can be computed with this loop:
! do i = 0, n
!   zev(i) = 0.5d0 * abs (r(i))
! enddo

!   Apply formula to produce Bernoulli numbers.

! t1 = mpreal (-2.d0, nwds)
! t2 = mpreal (1.d0, nwds)

call mpdmc (-2.d0, 0, t1, nwds)
call mpdmc (1.d0, 0, t2, nwds)

do i = 1, n
!   t1 = - dble (2*i-1) * dble (2*i) * t1
!   t2 = 4.d0 * cp2 * t2
!   berne(i) = mpreal (0.5d0 * t1 / t2 * abs (r(i)), mpnw)

  call mpmuld (t1, - dble (2*i-1) * dble (2*i), t3, nwds)
  call mpeq (t3, t1, nwds)
  call mpmuld (cp2, 4.d0, t3, nwds)
  call mpmul (t3, t2, t4, nwds)
  call mpeq (t4, t2, nwds)
  call mpmuld (t1, 0.5d0, t3, nwds)
  call mpdiv (t3, t2, t4, nwds)
  call mpabs (r(0:mpnw1+5,i), t3, nwds)
  call mpmul (t4, t3, berne(0:nb1+5,i), mpnw)
enddo

if (idb > 0) write (6, 10)
10 format ('Bernoulli number computation complete')
return
end subroutine mpberner

subroutine mppolyadd (nd1, n, a, b, c, nwds)

!   This adds two polynomials, as is required by mpberne.
!   The output array C may be the same as A or B.

implicit none
integer, intent(in):: n, nd1, nwds
integer k
integer (mpiknd), intent(in):: a(0:nd1+5,0:n), b(0:nd1+5,0:n)
integer (mpiknd), intent(out):: c(0:nd1+5,0:n)
integer (mpiknd) t1(0:nwds+5), t2(0:nwds+5)

call mpinitwds (t1, nwds)
call mpinitwds (t2, nwds)

do k = 0, n
!   c(k) = mpreal (a(k), nwds) + mpreal (b(k), nwds)

  call mpeq (a(0:nd1+5,k), t1, nwds)
  call mpeq (b(0:nd1+5,k), t2, nwds)
  call mpadd (t1, t2, c(0:nd1+5,k), nwds)
enddo

return
end subroutine mppolyadd

subroutine mppolysub (nd1, n, a, b, c, nwds)

!   This adds two polynomials, as is required by mpberne.
!   The output array C may be the same as A or B.

implicit none
integer, intent(in):: n, nd1, nwds
integer k
integer (mpiknd), intent(in):: a(0:nd1+5,0:n), b(0:nd1+5,0:n)
integer (mpiknd), intent(out):: c(0:nd1+5,0:n)
integer (mpiknd) t1(0:nwds+5), t2(0:nwds+5)

call mpinitwds (t1, nwds)
call mpinitwds (t2, nwds)

do k = 0, n
!   c(k) = mpreal (a(k), nwds) - mpreal (b(k), nwds)

  call mpeq (a(0:nd1+5,k), t1, nwds)
  call mpeq (b(0:nd1+5,k), t2, nwds)
  call mpsub (t1, t2, c(0:nd1+5,k), nwds)
enddo

return
end subroutine mppolysub

subroutine mppolymul (nd1, n, a, b, c, nwds)

!   This adds two polynomials (ignoring high-order terms), as is required by mpberne.
!   The output array C may NOT be the same as A or B.

implicit none
integer, intent(in):: n, nd1, nwds
integer j, k
integer (mpiknd), intent(in):: a(0:nd1+5,0:n), b(0:nd1+5,0:n)
integer (mpiknd), intent(out):: c(0:nd1+5,0:n)
integer (mpiknd) t0(0:nwds+5), t1(0:nwds+5), t2(0:nwds+5), t3(0:nwds+5)

call mpinitwds (t0, nwds)
call mpinitwds (t1, nwds)
call mpinitwds (t2, nwds)
call mpinitwds (t3, nwds)

do k = 0, n
!  t = mpreal (0.d0, nwds)
  call mpdmc (0.d0, 0, t0, nwds)

  do j = 0, k
!     t = t + mpreal (a(j), nwds) * mpreal (b(k-j), nwds)

    call mpeq (a(0:nd1+5,j), t1, nwds)
    call mpeq (b(0:nd1+5,k-j), t2, nwds)
    call mpmul (t1, t2, t3, nwds)
    call mpadd (t0, t3, t2, nwds)
    call mpeq (t2, t0, nwds)
  enddo

!   c(k) = t

  call mpeq (t0, c(0:nd1+5,k), nwds)
enddo

return
end subroutine mppolymul

subroutine mpbesselinr (nu, rr, ss, mpnw)

!   This evaluates the modified Bessel function BesselI (NU,RR).  
!   NU is an integer. The algorithm is DLMF formula 10.25.2.

implicit none
integer, intent(in):: nu, mpnw
integer ic1, itrmax, k, mpnw1, nu1, n1
real (mprknd) dmax, d1
parameter (itrmax = 1000000, dmax = 2000.d0)
integer (mpiknd), intent(in):: rr(0:)
integer (mpiknd), intent(out):: ss(0:)
integer (mpiknd) f1(0:mpnw+6), f2(0:mpnw+6), sum(0:mpnw+6), td(0:mpnw+6), &
  tn(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), &
  tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)
  
! End of declaration

if (mpnw < 4 .or. mpspacer (rr) < mpnw + 4 .or. mpspacer (ss) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPBESSELINR: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

!   Check input rr for proper range.

call mpmdc (rr, d1, n1, mpnw)
d1 = d1 * 2.d0 ** n1
if (mpsigntr (rr) < 0 .or. d1 > dmax) then
  write (mpldb, 2)
2 format ('*** MPBESSELINR: argument is negative or too large')
  call mpabrt (109)
endif

nu1 = abs (nu)
mpnw1 = mpnw + 1
call mpinitwds (f1, mpnw1)
call mpinitwds (f2, mpnw1)
call mpinitwds (sum, mpnw1)
call mpinitwds (td, mpnw1)
call mpinitwds (tn, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)

call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)

! tn = mpreal (1.d0, nwds)
! t1 = 0.25d0 * rr**2
! f1 = mpreal (1.d0, nwds)
! f2 = f1

call mpdmc (1.d0, 0, tn, mpnw1)
call mpdmc (1.d0, 0, f1, mpnw1)
call mpdmc (1.d0, 0, f2, mpnw1)
call mpmul (rr, rr, t2, mpnw1)
call mpmuld (t2, 0.25d0, t1, mpnw1)

! do k = 1, nu1
!   f2 = dble (k) * f2
! enddo

do k = 1, nu1
  call mpmuld (f2, dble (k), t2, mpnw1)
  call mpeq (t2, f2, mpnw1)
enddo

! td = f1 * f2
! t2 = tn / td
! sum = t2

call mpmul (f1, f2, td, mpnw1)
call mpdiv (tn, td, t2, mpnw1)
call mpeq (t2, sum, mpnw1)
  
! do k = 1, itrmax
!   f1 = dble (k) * f1
!   f2 = dble (k + nu1) * f2
!   tn = t1 * tn
!   td = f1 * f2
!   t2 = tn / td
!   sum = sum + t2
!   if (t2%mpr(3) < sum%mpr(3) - nwds) goto 100
! enddo

do k = 1, itrmax
  call mpmuld (f1, dble (k), t2, mpnw1)
  call mpeq (t2, f1, mpnw1)
  call mpmuld (f2, dble (k + nu1), t2, mpnw1)
  call mpeq (t2, f2, mpnw1)
  call mpmul (t1, tn, t2, mpnw1)
  call mpeq (t2, tn, mpnw1)
  call mpmul (f1, f2, td, mpnw1)
  call mpdiv (tn, td, t2, mpnw1)
  call mpadd (sum, t2, t3, mpnw1)
  call mpeq (t3, sum, mpnw1)
  
!  if (t2(0) == 0 .or. t2(3) < sum(3) - mpnw1) goto 100

  call mpabs (t2, tc1, 4)
  call mpmul (eps, sum, tc3, 4)
  call mpabs (tc3, tc2, 4)
  call mpcpr (tc1, tc2, ic1, 4)
  if (ic1 <= 0) goto 100
enddo

write (6, 3) 
3 format ('MPBESSELINR: End loop error')
call mpabrt (101)

100 continue

! besseli = (0.5d0 * rr) ** nu1 * sum

call mpmuld (rr, 0.5d0, t1, mpnw1)
call mpnpwr (t1, nu1, t2, mpnw1)
call mpmul (sum, t2, t3, mpnw1)

call mproun (t3, mpnw)
call mpeq (t3, ss, mpnw)

return
end subroutine mpbesselinr

subroutine mpbesseljnr (nu, rr, ss, mpnw)

!   This evaluates the modified Bessel function BesselJ (NU,RR).  
!   NU is an integer. The algorithm is DLMF formula 10.2.2.
!   This routine scales the working precision, based on the value of rr,
!   up to 2*mpnw words.

implicit none
integer, intent(in):: nu, mpnw
integer ic1, itrmax, k, mpnw1, mpnw2, nu1, n1
real (mprknd) dfact, dmax, d1
parameter (itrmax = 1000000, dfact = 2000.d0, dmax = 2000.d0)
integer (mpiknd), intent(in):: rr(0:)
integer (mpiknd), intent(out):: ss(0:)
integer (mpiknd) f1(0:2*mpnw+6), f2(0:2*mpnw+6), sum(0:2*mpnw+6), &
  td(0:2*mpnw+6), tn(0:2*mpnw+6), t1(0:2*mpnw+6), t2(0:2*mpnw+6), &
  t3(0:2*mpnw+6), t4(0:2*mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)
  
! End of declaration

if (mpnw < 4 .or. mpspacer (rr) < mpnw + 4 .or. mpspacer (ss) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPBESSELJNR: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

!   Check input rr for proper range.

call mpmdc (rr, d1, n1, mpnw)
d1 = d1 * 2.d0 ** n1
if (mpsigntr (rr) < 0 .or. d1 > dmax) then
  write (mpldb, 2)
2 format ('*** MPBESSELJNR: argument is negative or too large')
  call mpabrt (109)
endif

nu1 = abs (nu)
mpnw1 = mpnw + 1
mpnw2 = 2 * mpnw + 1
call mpinitwds (f1, mpnw2)
call mpinitwds (f2, mpnw2)
call mpinitwds (sum, mpnw2)
call mpinitwds (td, mpnw2)
call mpinitwds (tn, mpnw2)
call mpinitwds (t1, mpnw2)
call mpinitwds (t2, mpnw2)
call mpinitwds (t3, mpnw2)
call mpinitwds (t4, mpnw2)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)

! nwds2 = nwds * (1.d0 + dble (r1) / dfact)

mpnw2 = mpnw * (1.d0 + d1 / dfact)

! tn = mpreal (1.d0, nwds)
! t1 = 0.25d0 * rr**2
! f1 = mpreal (1.d0, nwds)
! f2 = f1

call mpdmc (1.d0, 0, tn, mpnw2)
call mpdmc (1.d0, 0, f1, mpnw2)
call mpdmc (1.d0, 0, f2, mpnw2)
call mpmul (rr, rr, t2, mpnw2)
call mpmuld (t2, 0.25d0, t1, mpnw2)

! do k = 1, nu1
!   f2 = dble (k) * f2
! enddo

do k = 1, nu1
  call mpmuld (f2, dble (k), t2, mpnw2)
  call mpeq (t2, f2, mpnw2)
enddo

! td = f1 * f2
! t2 = tn / td
! sum = t2

call mpmul (f1, f2, td, mpnw2)
call mpdiv (tn, td, t2, mpnw2)
call mpeq (t2, sum, mpnw2)
  
! do k = 1, itrmax
!   f1 = dble (k) * f1
!   f2 = dble (k + nu1) * f2
!   tn = - t1 * tn
!   td = f1 * f2
!   t2 = tn / td
!   sum = sum + t2
!   if (t2%mpr(3) < sum%mpr(3) - nwds) goto 100
! enddo

do k = 1, itrmax
  call mpmuld (f1, dble (k), t2, mpnw2)
  call mpeq (t2, f1, mpnw2)
  call mpmuld (f2, dble (k + nu1), t2, mpnw2)
  call mpeq (t2, f2, mpnw2)
  call mpmul (t1, tn, t2, mpnw2)
!  t2(2) = - t2(2)
  call mpneg (t2, tn, mpnw2)
  call mpmul (f1, f2, td, mpnw2)
  call mpdiv (tn, td, t2, mpnw2)
  call mpadd (sum, t2, t3, mpnw2)
  call mpeq (t3, sum, mpnw2)

!  if (t2(0) == 0 .or. t2(3) < sum(3) - mpnw2) goto 100

  call mpabs (t2, tc1, 4)
  call mpmul (eps, sum, tc3, 4)
  call mpabs (tc3, tc2, 4)
  call mpcpr (tc1, tc2, ic1, 4)
  if (ic1 <= 0) goto 100
enddo

write (6, 3) 
3 format ('MPBESSELJNR: End loop error')
call mpabrt (101)

100 continue

! besseli = (0.5d0 * rr) ** nu1 * sum

call mpmuld (rr, 0.5d0, t1, mpnw2)
call mpnpwr (t1, nu1, t2, mpnw2)
call mpmul (sum, t2, t3, mpnw2)

call mproun (t3, mpnw)
call mpeq (t3, ss, mpnw)

return
end subroutine mpbesseljnr

subroutine mpbesselknr (nu, rr, ss, mpnw)

!   This evaluates the modified Bessel function BesselK (NU,RR).  
!   NU is an integer. The algorithm is DLMF formula 10.31.1.

implicit none
integer, intent(in):: nu, mpnw
integer ic1, itrmax, k, mpnw1, nu1, n1
real (mprknd) dmax, d1, egam
parameter (itrmax = 1000000, dmax = 2000.d0, egam = 0.5772156649015328606d0)
integer (mpiknd), intent(in):: rr(0:)
integer (mpiknd), intent(out):: ss(0:)
integer (mpiknd) f1(0:mpnw+6), f2(0:mpnw+6), f3(0:mpnw+6), f4(0:mpnw+6), &
  f5(0:mpnw+6), sum1(0:mpnw+6), sum2(0:mpnw+6), sum3(0:mpnw+6), t1(0:mpnw+6), &
  t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)
  
! End of declaration

if (mpnw < 4 .or. mpspacer (rr) < mpnw + 4 .or. mpspacer (ss) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPBESSELKNR: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

!   Check input rr for proper range.

call mpmdc (rr, d1, n1, mpnw)
d1 = d1 * 2.d0 ** n1
if (mpsigntr (rr) < 0 .or. d1 > dmax) then
  write (mpldb, 2)
2 format ('*** MPBESSELKNR: argument is negative or too large')
  call mpabrt (109)
endif

!   Check if EGAMMA has been precomputed.

call mpmdc (mpegammacon, d1, n1, mpnw)
if (n1 /= -1 .or. abs (d1 * 2.d0**n1 - egam) > mprdfz &
  .or. mpwprecr (mpegammacon) < mpnw) then
  write (mpldb, 3) mpnw
3 format ('*** MPBESSELKNR: EGAMMA must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt (53)
endif

nu1 = abs (nu)
mpnw1 = mpnw + 1
call mpinitwds (f1, mpnw1)
call mpinitwds (f2, mpnw1)
call mpinitwds (f3, mpnw1)
call mpinitwds (f4, mpnw1)
call mpinitwds (f5, mpnw1)
call mpinitwds (sum1, mpnw1)
call mpinitwds (sum2, mpnw1)
call mpinitwds (sum3, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)

! t1 = 0.25d0 * rr**2
! f1 = mpreal (1.d0, nwds)
! f2 = f1
! f3 = f1
! sum1 = mpreal (0.d0, nwds)

call mpmul (rr, rr, t2, mpnw1)
call mpmuld (t2, 0.25d0, t1, mpnw1)
call mpdmc (1.d0, 0, f1, mpnw1)
call mpdmc (1.d0, 0, f2, mpnw1)
call mpdmc (1.d0, 0, f3, mpnw1)
call mpdmc (0.d0, 0,  sum1, mpnw1)

! do k = 1, nu1 - 1
!   f1 = dble (k) * f1
! enddo

do k = 1, nu1 - 1
  call mpmuld (f1, dble (k), t2, mpnw1)
  call mpeq (t2, f1, mpnw1)
enddo

! do k = 0, nu1 - 1
!    if (k > 0) then
!      f1 = f1 / dble (nu1 - k)
!      f2 = - t1 * f2
!      f3 = dble (k) * f3
!    endif
!    t2 = f1 * f2 / f3
!    sum1 = sum1 + t2
! enddo

do k = 0, nu1 - 1
  if (k > 0) then
    call mpdivd (f1, dble (nu1 - k), t2, mpnw1)
    call mpeq (t2, f1, mpnw1)
    call mpmul (t1, f2, t2, mpnw1)
!    t2(2) = - t2(2)
    call mpneg (t2, f2, mpnw1)
    call mpmuld (f3, dble (k), t2, mpnw1)
    call mpeq (t2, f3, mpnw1)
  endif
  call mpmul (f1, f2, t3, mpnw1)
  call mpdiv (t3, f3, t2, mpnw1)
  call mpadd (sum1, t2, t3, mpnw1)
  call mpeq (t3, sum1, mpnw1)
enddo

! sum1 = 0.5d0 * sum1 / (0.5d0 * rr)**nu1

call mpmuld (sum1, 0.5d0, t2, mpnw1)
call mpmuld (rr, 0.5d0, t3, mpnw1)
call mpnpwr (t3, nu1, t4, mpnw1)
call mpdiv (t2, t4, sum1, mpnw1)

! sum2 = (-1.d0) ** (nu1 + 1) * log (0.5d0 * rr) * besseli (nu1, rr)

call mpmuld (rr, 0.5d0, t2, mpnw1)
call mplog (t2, t3, mpnw1)
call mpmuld (t3, (-1.d0) ** (nu1 + 1), t2, mpnw1)
call mpbesselinr (nu1, rr, t3, mpnw1)
call mpmul (t2, t3, sum2, mpnw1)

! f1 = - egam
! f2 = - egam
! f3 = mpreal (1.d0, nwds)
! f4 = f3
! f5 = f3

call mpneg (mpegammacon, f1, mpnw1)
call mpeq (f1, f2, mpnw1)
call mpdmc (1.d0, 0, f3, mpnw1)
call mpdmc (1.d0, 0, f4, mpnw1)
call mpdmc (1.d0, 0, f5, mpnw1)

! do k = 1, nu1
!   f2 = f2 + mpreal (1.d0, nwds) / dble (k)
!   f5 = dble (k) * f5
! enddo

do k = 1, nu1
  call mpdmc (1.d0, 0, t2, mpnw1)
  call mpdivd (t2, dble (k), t3, mpnw1)
  call mpadd (f2, t3, t4, mpnw1)
  call mpeq (t4, f2, mpnw1)
  call mpmuld (f5, dble (k), t2, mpnw1)
  call mpeq (t2, f5, mpnw1)
enddo

! sum3 = (f1 + f2) * f3 / (f4 * f5)

call mpadd (f1, f2, t2, mpnw1)
call mpmul (t2, f3, t3, mpnw1)
call mpmul (f4, f5, t4, mpnw1)
call mpdiv (t3, t4, sum3, mpnw1)

! do k = 1, itrmax
!   f1 = f1 + mpreal (1.d0, nwds) / dble (k)
!   f2 = f2 + mpreal (1.d0, nwds) / dble (nu1 + k)
!   f3 = t1 * f3
!   f4 = dble (k) * f4
!   f5 = dble (nu1 + k) * f5
!   t2 = (f1 + f2) * f3 / (f4 * f5)
!   sum3 = sum3 + t2
!   if (t2%mpr(3) < sum3%mpr(3) - nwds) goto 100
! enddo

do k = 1, itrmax
  call mpdmc (1.d0, 0, t2, mpnw1)
  call mpdivd (t2, dble (k), t3, mpnw1)
  call mpadd (f1, t3, t4, mpnw1)
  call mpeq (t4, f1, mpnw1)
  call mpdivd (t2, dble (nu1 + k), t3, mpnw1)
  call mpadd (f2, t3, t4, mpnw1)
  call mpeq (t4, f2, mpnw1)
  call mpmul (t1, f3, t2, mpnw1)
  call mpeq (t2, f3, mpnw1)
  call mpmuld (f4, dble (k), t2, mpnw1)
  call mpeq (t2, f4, mpnw1)
  call mpmuld (f5, dble (nu1 + k), t2, mpnw1)
  call mpeq (t2, f5, mpnw1)
  call mpadd (f1, f2, t2, mpnw1)
  call mpmul (t2, f3, t3, mpnw1)
  call mpmul (f4, f5, t4, mpnw1)
  call mpdiv (t3, t4, t2, mpnw1)
  call mpadd (sum3, t2, t3, mpnw1)
  call mpeq (t3, sum3, mpnw1)
  
!  if (t2(0) == 0 .or. t2(3) < sum3(3) - mpnw1) goto 100

  call mpabs (t2, tc1, 4)
  call mpmul (eps, sum3, tc3, 4)
  call mpabs (tc3, tc2, 4)
  call mpcpr (tc1, tc2, ic1, 4)
  if (ic1 <= 0) goto 100 
enddo

write (6, 4) 
4 format ('MPBESSELKNR: end loop error')
call mpabrt (101)

100 continue

! sum3 = (-1.d0)**nu1 * 0.5d0 * (0.5d0 * rr)**nu1 * sum3
! besselk = sum1 + sum2 + sum3

call mpmuld (rr, 0.5d0, t2, mpnw1)
call mpnpwr (t2, nu1, t3, mpnw1)
call mpmuld (t3, (-1.d0)**nu1 * 0.5d0, t4, mpnw1)
call mpmul (t4, sum3, t2, mpnw1)
call mpeq (t2, sum3, mpnw1)
call mpadd (sum1, sum2, t2, mpnw1)
call mpadd (t2, sum3, t4, mpnw1)

call mproun (t4, mpnw)
call mpeq (t4, ss, mpnw)

return

end subroutine mpbesselknr

subroutine mpbesselynr (nu, rr, ss, mpnw)

!   This evaluates the modified Bessel function BesselY (NU,RR).  
!   NU is an integer. The algorithm is DLMF formula 10.8.1.

implicit none
integer, intent(in):: nu, mpnw
integer ic1, itrmax, k, mpnw1, mpnw2, nu1, n1
real (mprknd) dfact, dmax, d1, egam, pi
parameter (itrmax = 1000000, dfact = 2000.d0, dmax = 2000.d0, &
  egam = 0.5772156649015328606d0, pi = 3.1415926535897932385d0)
integer (mpiknd), intent(in):: rr(0:)
integer (mpiknd), intent(out):: ss(0:)
integer (mpiknd) f1(0:2*mpnw+6), f2(0:2*mpnw+6), f3(0:2*mpnw+6), &
  f4(0:2*mpnw+6), f5(0:2*mpnw+6), sum1(0:2*mpnw+6), &
  sum2(0:2*mpnw+6), sum3(0:2*mpnw+6), t1(0:2*mpnw+6), t2(0:2*mpnw+6), &
  t3(0:2*mpnw+6), t4(0:2*mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)
  
! End of declaration

if (mpnw < 4 .or. mpspacer (rr) < mpnw + 4 .or. mpspacer (ss) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPBESSELYNR: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

!   Check input rr for proper range.

call mpmdc (rr, d1, n1, mpnw)
d1 = d1 * 2.d0 ** n1
if (mpsigntr (rr) < 0 .or. d1 > dmax) then
  write (mpldb, 2)
2 format ('*** MPBESSELYNR: argument is negative or too large')
  call mpabrt (109)
endif
!   Check if EGAMMA has been precomputed.

call mpmdc (mpegammacon, d1, n1, mpnw)
if (n1 /= -1 .or. abs (d1 * 2.d0**n1 - egam) > mprdfz &
  .or. mpwprecr (mpegammacon) < mpnw) then
  write (mpldb, 3) mpnw
3 format ('*** MPBESSELYNR: EGAMMA must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt (53)
endif

!   Check if PI has been precomputed.

call mpmdc (mppicon, d1, n1, mpnw)
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw) then
  write (mpldb, 4) mpnw
4 format ('*** MPBESSELYNR: PI must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt (53)
endif

nu1 = abs (nu)
mpnw1 = mpnw + 1 
call mpinitwds  (f1, 2*mpnw+1)
call mpinitwds  (f2, 2*mpnw+1)
call mpinitwds  (f3, 2*mpnw+1)
call mpinitwds  (f4, 2*mpnw+1)
call mpinitwds  (f5, 2*mpnw+1)
call mpinitwds  (sum1, 2*mpnw+1)
call mpinitwds  (sum2, 2*mpnw+1)
call mpinitwds  (sum3, 2*mpnw+1)
call mpinitwds  (t1, 2*mpnw+1)
call mpinitwds  (t2, 2*mpnw+1)
call mpinitwds  (t3, 2*mpnw+1)
call mpinitwds  (t4, 2*mpnw+1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)

! nwds2 = nwds * (1.d0 + dble (r1) / dfact)

mpnw2 = mpnw * (1.d0 + d1 / dfact)

! t1 = 0.25d0 * rr**2
! f1 = mpreal (1.d0, nwds)
! f2 = f1
! f3 = f1
! sum1 = mpreal (0.d0, nwds)

call mpmul (rr, rr, t2, mpnw2)
call mpmuld (t2, 0.25d0, t1, mpnw2)
call mpdmc (1.d0, 0, f1, mpnw2)
call mpdmc (1.d0, 0, f2, mpnw2)
call mpdmc (1.d0, 0, f3, mpnw2)
call mpdmc (0.d0, 0,  sum1, mpnw2)

! do k = 1, nu1 - 1
!   f1 = dble (k) * f1
! enddo

do k = 1, nu1 - 1
  call mpmuld (f1, dble (k), t2, mpnw2)
  call mpeq (t2, f1, mpnw2)
enddo

! do k = 0, nu1 - 1
!    if (k > 0) then
!      f1 = f1 / dble (nu1 - k)
!      f2 = t1 * f2
!      f3 = dble (k) * f3
!    endif
!    t2 = f1 * f2 / f3
!    sum1 = sum1 + t2
! enddo

do k = 0, nu1 - 1
  if (k > 0) then
    call mpdivd (f1, dble (nu1 - k), t2, mpnw2)
    call mpeq (t2, f1, mpnw2)
    call mpmul (t1, f2, t2, mpnw2)
    call mpeq (t2, f2, mpnw2)
    call mpmuld (f3, dble (k), t2, mpnw2)
    call mpeq (t2, f3, mpnw2)
  endif
  call mpmul (f1, f2, t3, mpnw2)
  call mpdiv (t3, f3, t2, mpnw2)
  call mpadd (sum1, t2, t3, mpnw2)
  call mpeq (t3, sum1, mpnw2)
enddo

! sum1 = - sum1 / (0.5d0 * rr)**nu1

call mpmuld (rr, 0.5d0, t3, mpnw2)
call mpnpwr (t3, nu1, t4, mpnw2)
call mpdiv (sum1, t4, t3, mpnw2)
call mpneg (t3, sum1, mpnw2)

! sum2 = 2.d0 * log (0.5d0 * rr) * besselj (nu1, rr)

call mpmuld (rr, 0.5d0, t2, mpnw2)
call mplog (t2, t3, mpnw2)
call mpmuld (t3, 2.d0, t2, mpnw2)
call mpbesseljnr (nu1, rr, t3, mpnw2)
call mpmul (t2, t3, sum2, mpnw2)

! f1 = - egam
! f2 = - egam
! f3 = mpreal (1.d0, nwds)
! f4 = f3
! f5 = f3

call mpneg (mpegammacon, f1, mpnw2)
! f1(2) = - f1(2)
call mpeq (f1, f2, mpnw2)
call mpdmc (1.d0, 0, f3, mpnw2)
call mpdmc (1.d0, 0, f4, mpnw2)
call mpdmc (1.d0, 0, f5, mpnw2)

! do k = 1, nu1
!   f2 = f2 + mpreal (1.d0, nwds) / dble (k)
!   f5 = dble (k) * f5
! enddo

do k = 1, nu1
  call mpdmc (1.d0, 0, t2, mpnw2)
  call mpdivd (t2, dble (k), t3, mpnw2)
  call mpadd (f2, t3, t4, mpnw2)
  call mpeq (t4, f2, mpnw2)
  call mpmuld (f5, dble (k), t2, mpnw2)
  call mpeq (t2, f5, mpnw2)
enddo

! sum3 = (f1 + f2) * f3 / (f4 * f5)

call mpadd (f1, f2, t2, mpnw2)
call mpmul (t2, f3, t3, mpnw2)
call mpmul (f4, f5, t4, mpnw2)
call mpdiv (t3, t4, sum3, mpnw2)

! do k = 1, itrmax
!   f1 = f1 + mpreal (1.d0, nwds) / dble (k)
!   f2 = f2 + mpreal (1.d0, nwds) / dble (nu1 + k)
!   f3 = - t1 * f3
!   f4 = dble (k) * f4
!   f5 = dble (nu1 + k) * f5
!   t2 = (f1 + f2) * f3 / (f4 * f5)
!   sum3 = sum3 + t2
!   if (t2%mpr(3) < sum3%mpr(3) - nwds) goto 100
! enddo

do k = 1, itrmax
  call mpdmc (1.d0, 0, t2, mpnw2)
  call mpdivd (t2, dble (k), t3, mpnw2)
  call mpadd (f1, t3, t4, mpnw2)
  call mpeq (t4, f1, mpnw2)
  call mpdivd (t2, dble (nu1 + k), t3, mpnw2)
  call mpadd (f2, t3, t4, mpnw2)
  call mpeq (t4, f2, mpnw2)
  call mpmul (t1, f3, t2, mpnw2)
  call mpneg (t2, f3, mpnw2)
!  f3(2) = - f3(2)
  call mpmuld (f4, dble (k), t2, mpnw2)
  call mpeq (t2, f4, mpnw2)
  call mpmuld (f5, dble (nu1 + k), t2, mpnw2)
  call mpeq (t2, f5, mpnw2)
  call mpadd (f1, f2, t2, mpnw2)
  call mpmul (t2, f3, t3, mpnw2)
  call mpmul (f4, f5, t4, mpnw2)
  call mpdiv (t3, t4, t2, mpnw2)
  call mpadd (sum3, t2, t3, mpnw2)
  call mpeq (t3, sum3, mpnw2)

!   if (t2(0) == 0 .or. t2(3) < sum3(3) - mpnw2) goto 100

  call mpabs (t2, tc1, 4)
  call mpmul (eps, sum3, tc3, 4)
  call mpabs (tc3, tc2, 4)
  call mpcpr (tc1, tc2, ic1, 4)
  if (ic1 <= 0) goto 100
enddo

write (6, 5) 
5 format ('MPBESSELYNR: end loop error')
call mpabrt (101)

100 continue

! sum3 = - (0.5d0 * rr)**nu1 * sum3
! besselk = (sum1 + sum2 + sum3) / mppi (nwds)

call mpmuld (rr, 0.5d0, t2, mpnw2)
call mpnpwr (t2, nu1, t3, mpnw2)
call mpmul (t3, sum3, t2, mpnw2)
call mpneg (t2, sum3, mpnw2)

call mpadd (sum1, sum2, t2, mpnw2)
call mpadd (t2, sum3, t4, mpnw2)
call mpeq (mppicon, t2, mpnw2)
call mpdiv (t4, t2, t3, mpnw2)

call mproun (t3, mpnw)
call mpeq (t3, ss, mpnw)

return
end subroutine mpbesselynr

subroutine mperfr (z, terf, mpnw)

!   This evaluates the erf function, using a combination of two series.
!   In particular, the algorithm is (where B = (mpnw + 1) * mpnbt, and
!   dcon is a constant defined below):

!   if (t == 0) then
!     erf = 0
!   elseif (z > sqrt(B*log(2))) then
!     erf = 1
!   elseif (z < -sqrt(B*log(2))) then
!     erf = -1
!   elseif (abs(z) < B/dcon + 8) then
!     erf = 2 / (sqrt(pi)*exp(z^2)) * Sum_{k>=0} 2^k * z^(2*k+1) 
!             / (1.3....(2*k+1))
!   else
!     erf = 1 - 1 / (sqrt(pi)*exp(z^2)) 
!             * Sum_{k>=0} (-1)^k * (1.3...(2*k-1)) / (2^k * z^(2*k+1))
!   endif

implicit none
integer, intent(in):: mpnw
integer ic1, ic2, ic3, itrmx, k, mpnw1, nbt
real (mprknd) dcon, d1, d2
parameter (dcon = 100.d0, itrmx = 100000)
integer (mpiknd), intent(in):: z(0:)
integer (mpiknd), intent(out):: terf(0:)
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), &
  t5(0:mpnw+6), t6(0:mpnw+6), t7(0:mpnw+6), z2(0:mpnw+6), tc1(0:9), &
  tc2(0:9), tc3(0:9), eps(0:9)

! End of declaration

if (mpnw < 4 .or. mpspacer (z) < mpnw + 4 .or. mpspacer (terf) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPERFR: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

mpnw1 = mpnw + 1
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (t6, mpnw1)
call mpinitwds (t7, mpnw1)
call mpinitwds (z2, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)

nbt = mpnw * mpnbt
d1 = aint (1.d0 + sqrt (nbt * log (2.d0)))
d2 = aint (nbt / dcon + 8.d0)
call mpdmc (d1, 0, t1, mpnw1)
call mpdmc (d2, 0, t2, mpnw1)
call mpcpr (z, t1, ic1, mpnw1)
! t1(2) = - t1(2)
call mpneg (t1, t3, mpnw1)
call mpeq (t3, t1, mpnw1)
call mpcpr (z, t1, ic2, mpnw1)
call mpcpr (z, t2, ic3, mpnw1)

!if (z == 0.d0) then
if (mpsigntr (z) == 0) then

!  terf = mpreal (0.d0, nwds)
  call mpdmc (0.d0, 0, terf, mpnw)

!elseif (z > d2) then
elseif (ic1 > 0) then

!  terf = mpreal (1.d0, nwds)
  call mpdmc (1.d0, 0, terf, mpnw)

!elseif (z < -d2) then
elseif (ic2 < 0) then

!  terf = mpreal (-1.d0, nwds)
  call mpdmc (-1.d0, 0, terf, mpnw)

!elseif (abs (z) < d3) then
elseif (ic3 < 0) then

!  z2 = z**2
  call mpmul (z, z, z2, mpnw1)

!  t1 = mpreal (0.d0, nwds)
  call mpdmc (0.d0, 0, t1, mpnw1)

!  t2 = z
  call mpeq (z, t2, mpnw1)

!  t3 = mpreal (1.d0, nwds)
  call mpdmc (1.d0, 0, t3, mpnw1)

!  t5 = mpreal (1.d10, 4)
  call mpdmc (1.d10, 0, t5, 4)

  do k = 0, itrmx
    if (k > 0) then
!      t2 = 2.d0 * z2 * t2
      call mpmuld (z2, 2.d0, t6, mpnw1)
      call mpmul (t6, t2, t7, mpnw1)
      call mpeq (t7, t2, mpnw1)

!      t3 = (2.d0 * k + 1.d0) * t3
      call mpmuld (t3, 2.d0 * k + 1.d0, t6, mpnw1)
      call mpeq (t6, t3, mpnw1)
    endif

!    t4 = t2 / t3
    call mpdiv (t2, t3, t4, mpnw1)

!    t1 = t1 + t4
    call mpadd (t1, t4, t6, mpnw1)
    call mpeq (t6, t1, mpnw1)
        
!    t6 = abs (mpreal (t4, 4) / mpreal (t1, 4))
    call mpdiv (t4, t1, t6, 4)

!    if (t6 < eps .or. t6 >= t5) goto 120
    call mpcpr (t6, eps, ic1, 4)
    call mpcpr (t6, t5, ic2, 4)
    if (ic1 <= 0 .or. ic2 >= 0) goto 120

!    t5 = t6
    call mpeq (t6, t5, 4)
  enddo
  
write (6, 3) 1, itrmx
3 format ('*** MPERFR: iteration limit exceeded',2i10)
call mpabrt (101)

120 continue

!  terf = 2.d0 * t1 / (sqrt (mppi (nwds)) * exp (z2))
  call mpmuld (t1, 2.d0, t3, mpnw1)
  call mpsqrt (mppicon, t4, mpnw1)
  call mpexp (z2, t5, mpnw1)
  call mpmul (t4, t5, t6, mpnw1)
  call mpdiv (t3, t6, t7, mpnw1)
  call mproun (t7, mpnw)
  call mpeq (t7, terf, mpnw)
else

!  z2 = z ** 2
  call mpmul (z, z, z2, mpnw1)

!  t1 = mpreal (0.d0, nwds)
  call mpdmc (0.d0, 0, t1, mpnw1)

!  t2 = mpreal (1.d0, nwds)
  call mpdmc (1.d0, 0, t2, mpnw1)

!  t3 = abs (z)
  call mpabs (z, t3, mpnw1)

!  t5 = mpreal (1.d10, 4)
  call mpdmc (1.d10, 0, t5, 4)
    
  do k = 0, itrmx
    if (k > 0) then
!      t2 = - (2.d0 * k - 1.d0) * t2 / 2.d0
      call mpmuld (t2, -(2.d0 * k - 1.d0), t6, mpnw1)
      call mpeq (t6, t2, mpnw1)

!      t3 = z2 * t3
      call mpmul (t2, t3, t6, mpnw1)
      call mpeq (t6, t3, mpnw1)
    endif

!    t4 = t2 / t3
    call mpdiv (t2, t3, t4, mpnw1)

!    t1 = t1 + t4
    call mpadd (t1, t4, t6, mpnw1)
    call mpeq (t6, t1, mpnw1)
    
!    t6 = abs (mpreal (t4, 4) / mpreal (t1, 4))
    call mpdiv (t4, t1, t6, 4)

!    if (t6 < eps .or. t6 >= t5) goto 130
    call mpcpr (t6, eps, ic1, 4)
    call mpcpr (t6, t5, ic2, 4)
    if (ic1 <= 0 .or. ic2 >= 0) goto 130

!    t5 = t6
    call mpeq (t6, t5, 4)
  enddo

write (6, 3) 2, itrmx
call mpabrt (101)

130 continue

!  terf = 1.d0 - t1 / (sqrt (mppi (nwds)) * exp (z2))
  call mpdmc (1.d0, 0, t2, mpnw1)
  call mpsqrt (mppicon, t3, mpnw1)
  call mpexp (z2, t4, mpnw1)
  call mpmul (t3, t4, t5, mpnw1)
  call mpdiv (t1, t5, t6, mpnw1)
  call mpsub (t2, t6, t7, mpnw1)
  call mproun (t7, mpnw)
  call mpeq (t7, terf, mpnw)

!  if (z < 0.d0) terf = - terf
!  if (z(2) < 0) terf(2) = - terf(2)
  if (mpsigntr (z) < 0) then
    call mpneg (terf, t6, mpnw)
    call mpeq (t6, terf, mpnw)
  endif
  
endif

return
end subroutine mperfr

subroutine mperfcr (z, terfc, mpnw)

!   This evaluates the erfc function, using a combination of two series.
!   In particular, the algorithm is (where B = (mpnw + 1) * mpnbt, and
!   dcon is a constant defined below):

!   if (t == 0) then
!     erfc = 1
!   elseif (z > sqrt(B*log(2))) then
!     erfc = 0
!   elseif (z < -sqrt(B*log(2))) then
!     erfc = 2
!   elseif (abs(z) < B/dcon + 8) then
!     erfc = 1 - 2 / (sqrt(pi)*exp(z^2)) * Sum_{k>=0} 2^k * z^(2*k+1) 
!               / (1.3....(2*k+1))
!   else
!     erfc = 1 / (sqrt(pi)*exp(z^2)) 
!             * Sum_{k>=0} (-1)^k * (1.3...(2*k-1)) / (2^k * z^(2*k+1))
!   endif

implicit none
integer, intent(in):: mpnw
integer ic1, ic2, ic3, itrmx, k, mpnw1, nbt
real (mprknd) dcon, d1, d2
parameter (dcon = 100.d0, itrmx = 100000)
integer (mpiknd), intent(in):: z(0:)
integer (mpiknd), intent(out):: terfc(0:)
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), &
  t5(0:mpnw+6), t6(0:mpnw+6), t7(0:mpnw+6), z2(0:mpnw+6), tc1(0:9), &
  tc2(0:9), tc3(0:9), eps(0:9)

! End of declaration

if (mpnw < 4 .or. mpspacer (z) < mpnw + 4 .or. mpspacer (terfc) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPERFCR: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

mpnw1 = mpnw + 1
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (t6, mpnw1)
call mpinitwds (t7, mpnw1)
call mpinitwds (z2, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)

nbt = mpnw * mpnbt
d1 = aint (1.d0 + sqrt (nbt * log (2.d0)))
d2 = aint (nbt / dcon + 8.d0)
call mpdmc (d1, 0, t1, mpnw1)
call mpdmc (d2, 0, t2, mpnw1)
call mpcpr (z, t1, ic1, mpnw1)
! t1(2) = - t1(2)
call mpneg (t1, t3, mpnw1)
call mpeq (t3, t1,  mpnw1)
call mpcpr (z, t1, ic2, mpnw1)
call mpcpr (z, t2, ic3, mpnw1)

!if (z == 0.d0) then
if (mpsigntr (z) == 0) then

!  terfc = mpreal (1.d0, nwds)
  call mpdmc (1.d0, 0, terfc, mpnw)

!elseif (z > d2) then
elseif (ic1 > 0) then

!  terfc = mpreal (0.d0, nwds)
  call mpdmc (0.d0, 0, terfc, mpnw)

!elseif (z < -d2) then
elseif (ic2 < 0) then

!  terfc = mpreal (2.d0, nwds)
  call mpdmc (2.d0, 0, terfc, mpnw)

!elseif (abs (z) < d3) then
elseif (ic3 < 0) then

!  z2 = z**2
  call mpmul (z, z, z2, mpnw1)

!  t1 = mpreal (0.d0, nwds)
  call mpdmc (0.d0, 0, t1, mpnw1)

!  t2 = z
  call mpeq (z, t2, mpnw1)

!  t3 = mpreal (1.d0, nwds)
  call mpdmc (1.d0, 0, t3, mpnw1)

!  t5 = mpreal (1.d10, 4)
  call mpdmc (1.d10, 0, t5, 4)

  do k = 0, itrmx
    if (k > 0) then
!      t2 = 2.d0 * z2 * t2
      call mpmuld (z2, 2.d0, t6, mpnw1)
      call mpmul (t6, t2, t7, mpnw1)
      call mpeq (t7, t2, mpnw1)

!      t3 = (2.d0 * k + 1.d0) * t3
      call mpmuld (t3, 2.d0 * k + 1.d0, t6, mpnw1)
      call mpeq (t6, t3, mpnw1)
    endif

!    t4 = t2 / t3
    call mpdiv (t2, t3, t4, mpnw1)

!    t1 = t1 + t4
    call mpadd (t1, t4, t6, mpnw1)
    call mpeq (t6, t1, mpnw1)
        
!    t6 = abs (mpreal (t4, 4) / mpreal (t1, 4))
    call mpdiv (t4, t1, t6, 4)

!    if (t6 < eps .or. t6 >= t5) goto 120
    call mpcpr (t6, eps, ic1, 4)
    call mpcpr (t6, t5, ic2, 4)
    if (ic1 <= 0 .or. ic2 >= 0) goto 120

!    t5 = t6
    call mpeq (t6, t5, 4)
  enddo
  
write (6, 3) 1, itrmx
3 format ('*** MPERFCR: iteration limit exceeded',2i10)
call mpabrt (101)

120 continue

!  terfc = 1.d0 - 2.d0 * t1 / (sqrt (mppi (nwds)) * exp (z2))
  call mpdmc (1.d0, 0, t2, mpnw1)
  call mpmuld (t1, 2.d0, t3, mpnw1)
  call mpsqrt (mppicon, t4, mpnw1)
  call mpexp (z2, t5, mpnw1)
  call mpmul (t4, t5, t6, mpnw1)
  call mpdiv (t3, t6, t7, mpnw1)
  call mpsub (t2, t7, t6, mpnw1)
  call mproun (t6, mpnw)
  call mpeq (t6, terfc, mpnw)
else

!  z2 = z ** 2
  call mpmul (z, z, z2, mpnw1)

!  t1 = mpreal (0.d0, nwds)
  call mpdmc (0.d0, 0, t1, mpnw1)

!  t2 = mpreal (1.d0, nwds)
  call mpdmc (1.d0, 0, t2, mpnw1)

!  t3 = abs (z)
  call mpabs (z, t3, mpnw1)

!  t5 = mpreal (1.d10, 4)
  call mpdmc (1.d10, 0, t5, 4)
    
  do k = 0, itrmx
    if (k > 0) then
!      t2 = - (2.d0 * k - 1.d0) * t2 / 2.d0
      call mpmuld (t2, -(2.d0 * k - 1.d0), t6, mpnw1)
      call mpeq (t6, t2, mpnw1)

!      t3 = z2 * t3
      call mpmul (t2, t3, t6, mpnw1)
      call mpeq (t6, t3, mpnw1)
    endif

!    t4 = t2 / t3
    call mpdiv (t2, t3, t4, mpnw1)

!    t1 = t1 + t4
    call mpadd (t1, t4, t6, mpnw1)
    call mpeq (t6, t1, mpnw1)
    
!    t6 = abs (mpreal (t4, 4) / mpreal (t1, 4))
    call mpdiv (t4, t1, t6, 4)

!    if (t6 < eps .or. t6 >= t5) goto 130
    call mpcpr (t6, eps, ic1, 4)
    call mpcpr (t6, t5, ic2, 4)
    if (ic1 <= 0 .or. ic2 >= 0) goto 130

!    t5 = t6
    call mpeq (t6, t5, 4)
  enddo

write (6, 3) 2, itrmx
call mpabrt (101)

130 continue

!  terfc = t1 / (sqrt (mppi (nwds)) * exp (z2))
  call mpsqrt (mppicon, t3, mpnw1)
  call mpexp (z2, t4, mpnw1)
  call mpmul (t3, t4, t5, mpnw1)
  call mpdiv (t1, t5, t6, mpnw1)

!  if (z < 0.d0) terfc = 2.d0 - terfc
  if (mpsigntr (z) < 0) then
    call mpdmc (2.d0, 0, t2, mpnw1)
    call mpsub (t2, t6, t7, mpnw1)
    call mpeq (t7, t6, mpnw1)
  endif

  call mproun (t6, mpnw)
  call mpeq (t6, terfc, mpnw)
endif

return
end subroutine mperfcr

subroutine mpgammar (t, z, mpnw)

!   This evaluates the gamma function, using an algorithm of R. W. Potter.
!   The argument t must not exceed 10^8 in size (this limit is set below),
!   must not be zero, and if negative must not be integer.

!   In the parameter statement below:
!     itrmx = limit of number of iterations in series; default = 100000.
!     con1 = 1/2 * log (10) to DP accuracy.
!     dmax = maximum size of input argument.

implicit none
integer, intent(in):: mpnw
integer i, ic1, itrmx, j, mpnw1, nt, n1, n2, n3
real (mprknd) alpha, al2, dmax, d1, d2, d3
parameter (al2 = 0.69314718055994530942d0, dmax = 1d8, itrmx = 100000)
integer (mpiknd), intent(in):: t(0:)
integer (mpiknd), intent(out):: z(0:)
integer (mpiknd) f1(0:mpnw+6), sum1(0:mpnw+6), sum2(0:mpnw+6), tn(0:mpnw+6), &
  t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), t5(0:mpnw+6), &
  t6(0:mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)

! End of declaration

if (mpnw < 4 .or. mpspacer (t) < mpnw + 4 .or. mpspacer (z) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPGAMMAR: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

mpnw1 = mpnw + 1
call mpinitwds (f1, mpnw1)
call mpinitwds (sum1, mpnw1)
call mpinitwds (sum2, mpnw1)
call mpinitwds (tn, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (t6, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)

call mpmdc (t, d1, n1, mpnw)
d1 = d1 * 2.d0**n1
call mpnint (t, t1, mpnw)
call mpcpr (t, t1, ic1, mpnw)
if (mpsigntr (t) == 0 .or. d1 > dmax .or. (mpsigntr (t) < 0 .and. ic1 == 0)) then
  write (6, 2) dmax
2 format ('*** MPGAMMAR: input argument must have absolute value <=',f10.0,','/ &
  'must not be zero, and if negative must not be an integer.')
  call mpabrt (65)
endif

call mpdmc (1.d0, 0, f1, mpnw1)

!   Find the integer and fractional parts of t.

call mpinfr (t, t2, t3, mpnw1)

if (mpsigntr (t3) == 0) then

!   If t is a positive integer, then apply the usual factorial recursion.

  call mpmdc (t2, d2, n2, mpnw1)
  nt = d2 * 2.d0 ** n2
  call mpeq (f1, t1, mpnw1)

  do i = 2, nt - 1
    call mpmuld (t1, dble (i), t2, mpnw1)
    call mpeq (t2, t1, mpnw1)
  enddo

  call mproun (t1, mpnw)
  call mpeq (t1, z, mpnw)
  goto 120
elseif (mpsigntr (t) > 0) then

!   Apply the identity Gamma[t+1] = t * Gamma[t] to reduce the input argument
!   to the unit interval.

  call mpmdc (t2, d2, n2, mpnw1)
  nt = d2 * 2.d0 ** n2
  call mpeq (f1, t1, mpnw1)
  call mpeq (t3, tn, mpnw1)
  
  do i = 1, nt
    call mpdmc (dble (i), 0, t4, mpnw1)
    call mpsub (t, t4, t5, mpnw1)
    call mpmul (t1, t5, t6, mpnw1)
    call mpeq (t6, t1, mpnw1)
  enddo
else

!   Apply the gamma identity to reduce a negative argument to the unit interval.

  call mpsub (f1, t, t4, mpnw1)
  call mpinfr (t4, t3, t5, mpnw1)
  call mpmdc (t3, d3, n3, mpnw1)
  nt = d3 * 2.d0 ** n3

  call mpeq (f1, t1, mpnw1)
  call mpsub (f1, t5, t2, mpnw1)
  call mpeq (t2, tn, mpnw1)
    
  do i = 0, nt - 1
!    t1 = t1 / (t + dble (i))
    call mpdmc (dble (i), 0, t4, mpnw1)
    call mpadd (t, t4, t5, mpnw1)
    call mpdiv (t1, t5, t6, mpnw1)
    call mpeq (t6, t1, mpnw1)
  enddo
endif

!   Calculate alpha = bits of precision * log(2) / 2, then take the next integer
!   value mod 4, so that d2 = 0.25 * alpha^2 can be calculated exactly in DP.

alpha = 4.d0 * aint ((0.5d0 * mpnbt * al2 * (mpnw1 + 1)) / 4.d0 + 1.d0)
d2 = 0.25d0 * alpha**2

call mpeq (tn, t2, mpnw1)
call mpdiv (f1, t2, t3, mpnw1)
call mpeq (t3, sum1, mpnw1)

!   Evaluate the series with t.

do j = 1, itrmx
  call mpdmc (dble (j), 0, t6, mpnw1)
  call mpadd (t2, t6, t4, mpnw1)
  call mpmuld (t4, dble (j), t5, mpnw1)
  call mpdiv (t3, t5, t6, mpnw1)
  call mpmuld (t6, d2, t3, mpnw1)
  call mpadd (sum1, t3, t4, mpnw1)
  call mpeq (t4, sum1, mpnw1)

!  if (t3(2) == 0 .or. t3(3) < sum1(3) - mpnw1) goto 100

  call mpabs (t3, tc1, 4)
  call mpmul (eps, sum1, tc3, 4)
  call mpabs (tc3, tc2, 4)
  call mpcpr (tc1, tc2, ic1, 4)
  if (ic1 <= 0) goto 100 
enddo

write (6, 3) 1, itrmx
3 format ('*** MPGAMMAR: iteration limit exceeded',2i10)
call mpabrt (101)

100 continue

call mpneg (tn, t2, mpnw1)
! t2(2) = - t2(2)
call mpdiv (f1, t2, t3, mpnw1)
call mpeq (t3, sum2, mpnw1)

!   Evaluate the same series with -t.

do j = 1, itrmx
  call mpdmc (dble (j), 0, t6, mpnw1)
  call mpadd (t2, t6, t4, mpnw1)
  call mpmuld (t4, dble (j), t5, mpnw1)
  call mpdiv (t3, t5, t6, mpnw1)
  call mpmuld (t6, d2, t3, mpnw1)
  call mpadd (sum2, t3, t4, mpnw1)
  call mpeq (t4, sum2, mpnw1)

!  if (t3(2) == 0 .or. t3(3) < sum2(3) - mpnw1) goto 110

  call mpabs (t3, tc1, 4)
  call mpmul (eps, sum2, tc3, 4)
  call mpabs (tc3, tc2, 4)
  call mpcpr (tc1, tc2, ic1, 4)
  if (ic1 <= 0) goto 110
enddo

write (6, 3) 2, itrmx
call mpabrt (67)

110 continue

!   Compute sqrt (mppic * sum1 / (tn * sin (mppic * tn) * sum2)) 
!   and (alpha/2)^tn terms.

call mpeq (mppicon, t2, mpnw1)
call mpmul (t2, tn, t3, mpnw1)
call mpcssnr (t3, t4, t5, mpnw1)

call mpmul (t5, sum2, t6, mpnw1)
call mpmul (tn, t6, t5, mpnw1)
call mpmul (t2, sum1, t3, mpnw1)
call mpdiv (t3, t5, t6, mpnw1)
! t6(2) = - t6(2)
call mpneg (t6, t4, mpnw1)
call mpeq (t4, t6, mpnw1)
call mpsqrt (t6, t2, mpnw1)
call mpdmc (0.5d0 * alpha, 0, t3, mpnw1)
call mplog (t3, t4, mpnw1)
call mpmul (tn, t4, t5, mpnw1)
call mpexp (t5, t6, mpnw1)
call mpmul (t2, t6, t3, mpnw1)
call mpmul (t1, t3, t4, mpnw1)

!   Round to mpnw words precision.

call mproun (t4, mpnw)
call mpeq (t4, z, mpnw)

120 continue

return
end subroutine mpgammar

subroutine mpincgammar (s, z, g, mpnw)

!  This returns the incomplete gamma function, using a combination of formula
!  8.7.3 of the DLMF (for modest-sized z) and formula 8.11.2 (for large z).

implicit none
integer, intent(in):: mpnw
integer ic1, itrmax, k, mpnw1, n1
real (mprknd) d1, dmax
parameter (dmax = 0.6667d0, itrmax = 1000000)
integer (mpiknd), intent(in):: s(0:), z(0:)
integer (mpiknd), intent(out):: g(0:)
integer (mpiknd) t0(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), &
  t4(0:mpnw+6), t5(0:mpnw+6), f1(0:mpnw+6), tc1(0:9), tc2(0:9), &
  tc3(0:9), eps(0:9)

! End of declaration

if (mpnw < 4 .or. mpspacer (s) < mpnw + 4 .or. mpspacer (z) < mpnw + 4 &
  .or. mpspacer (g) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPINCGAMMAR: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

mpnw1 = mpnw + 1
call mpinitwds (t0, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (f1, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)

call mpdmc (1.d0, 0, f1, mpnw1)

! if (abs (z) < dmax * mpnw * mpnbt) then
!  if (z(3) < 0 .or. (z(3) == 0 .and. z(4) < dmax * mpnw)) then

call mpmdc (z, d1, n1, mpnw1)
d1 = d1 * 2.d0 ** n1

if (abs (d1) < dmax * mpnw1 * mpnbt) then

!  t1 = gamma (s)

  call mpgammar (s, t1, mpnw1)

!  t2 = 1.d0 / (s * t1)

  call mpmul (s, t1, t3, mpnw1)
  call mpdiv (f1, t3, t2, mpnw1)

!   t0 = t2

  call mpeq (t2, t0, mpnw1)
    
  do k = 1, itrmax

!    t2 = t2 * z / (s + dble (k))

    call mpmul (t2, z, t5, mpnw1)
    call mpdmc (dble (k), 0, t3, mpnw1)
    call mpadd (s, t3, t4, mpnw1)
    call mpdiv (t5, t4, t2, mpnw1)

!    t0 = t0 + t2

    call mpadd (t0, t2, t3, mpnw1)
    call mpeq (t3, t0, mpnw1)

!    if (t2(2) == 0 .or. t2(3) < t0(3) - mpnw) goto 100

    call mpabs (t2, tc1, 4)
    call mpmul (eps, t0, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic1, 4)
    if (ic1 <= 0) goto 100 
  enddo
  
  write (mpldb, 2) 1, itrmax
2 format ('*** MPINCGAMMAR: iteration limit exceeded:',2i10)
  call mpabrt (101)

100 continue

!   gammainc = t1 * (1.d0 - z ** s / exp (z) * t0)

  call mppower (z, s, t2, mpnw1)
  call mpexp (z, t3, mpnw1)
  call mpdiv (t2, t3, t4, mpnw1)
  call mpmul (t4, t0, t5, mpnw1)
  call mpsub (f1, t5, t2, mpnw1)
  call mpmul (t1, t2, t3, mpnw1)
  call mpeq (t3, t1, mpnw1)
else
!  t0 = mpreal (1.d0, mpnw)

  call mpdmc (1.d0, 0, t0, mpnw1)
  
!  t1 = mpreal (1.d0, mpnw)

  call mpdmc (1.d0, 0, t1, mpnw1)
  
  do k = 1, itrmax
!    t1 = t1 * (s - dble (k)) / z

    call mpdmc (dble (k), 0, t2, mpnw1)
    call mpsub (s, t2, t3, mpnw1)
    call mpmul (t1, t3, t4, mpnw1)
    call mpdiv (t4, z, t1, mpnw1)

!    t0 = t0 + t1

    call mpadd (t0, t1, t2, mpnw1)
    call mpeq (t2, t0, mpnw1)

!    if (t1(2) == 0 .or. t1(3) < t0(3) - mpnw) goto 110

    call mpabs (t2, tc1, 4)
    call mpmul (eps, t0, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic1, 4)
    if (ic1 <= 0) goto 110 
  enddo

  write (mpldb, 2) 2, itrmax
  call mpabrt (101)

110 continue

!  gammainc = z ** (s - 1.d0) / exp (z) * t0

   call mpsub (s, f1, t2, mpnw1)
   call mppower (z, t2, t3, mpnw1)
   call mpexp (z, t4, mpnw1)
   call mpdiv (t3, t4, t2, mpnw1)
   call mpmul (t2, t0, t1, mpnw1)
endif

200 continue

call mproun (t1, mpnw)
call mpeq (t1, g, mpnw)

return
end subroutine mpincgammar

subroutine mppolylogini (na, nn, arr, mpnw)

!   Initializes the MP array arr with data for mppolylogneeg. 
!   NN must be in the range (-nmax, -1).

implicit none
integer, intent(in):: na, nn, mpnw
integer i1, i2, k, n, nmax, nna, mpnw1
parameter (nmax = 1000)
integer (mpiknd), intent(out):: arr(0:na+5,1:abs(nn))
integer (mpiknd) aa(0:mpnw+6,2,abs(nn)), t1(0:mpnw+6), t2(0:mpnw+6)

!  End of declaration

if (nn < -nmax .or. nn >= 0) then
  write (mpldb, 1) -nmax
1 format ('*** MPPOLYLOGINI: N is <',i6,' or n >= 0.'/ &
  'For n >= 0, call mppolylogpos or polylog_pos.')
  call mpabrt (111)
endif

nna = abs (nn)
i1 = 1000000000

do k = 1, nna
  i1 = min (i1, mpspacer (arr(0:na+5,k)))
enddo

if (mpnw < 4 .or. i1 < mpnw + 6) then
  write (mpldb, 2)
2 format ('*** MPPOLYLOGINI: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

mpnw1 = mpnw + 1
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
i1 = 2
i2 = 1

! aa(1,1) = mpreal (1.d0, nwds)
! aa(2,1) = mpreal (1.d0, nwds)

aa(0,1,1) = mpnw + 7
aa(0,2,1) = mpnw + 7
call mpdmc (1.d0, 0, aa(0:mpnw+6,1,1), mpnw1)
call mpdmc (1.d0, 0, aa(0:mpnw+6,2,1), mpnw1)

do k = 2, nna
!  aa(1,i) = mpreal (0.d0, nwds)
!  aa(2,i) = mpreal (0.d0, nwds)

  aa(0,1,k) = mpnw + 7
  aa(0,2,k) = mpnw + 7
  call mpdmc (0.d0, 0, aa(0:mpnw+6,1,k), mpnw1)
  call mpdmc (0.d0, 0, aa(0:mpnw+6,2,k), mpnw1)
enddo

do n = 2, nna
  i1 = 3 - i1
  i2 = 3 - i1

  do k = 2, n
!    aa(i2,k) = dble (n + 1 - k) * aa(i1,k-1) + dble (k) * aa(i1,k)

    call mpmuld (aa(0:mpnw+6,i1,k-1), dble (n + 1 - k), t1, mpnw1)
    call mpmuld (aa(0:mpnw+6,i1,k), dble (k), t2, mpnw1)
    call mpadd (t1, t2, aa(0:mpnw+6,i2,k), mpnw1)
  enddo
enddo

do k = 1, nna
!  arr(k) = aa(i2,k)

  call mpeq (aa(0:mpnw+6,i2,k), arr(0:na+5,k), mpnw)
enddo

return
end subroutine mppolylogini

subroutine mppolylogneg (na, nn, arr, x, y, mpnw)

!   This returns polylog (nn, x) for the case nn < 0. Before calling this,
!   one must call mppolylognini to initialize the array arr.
!   NN must be in the range (-nmax, -1).
!   The parameter nmxa is the maximum number of additional words of
!   precision needed to overcome cancelation errors when x is negative,
!   for nmax = 1000.

implicit none
integer, intent(in):: na, nn, mpnw
integer i1, i2, k, mpnw1, n1, n2, nna, nmax, nmxa
parameter (nmax = 1000, nmxa = 8525 / mpnbt + 1)
real (mprknd) d1, d2
integer (mpiknd), intent(in):: arr(0:na+5,1:abs(nn)), x(0:)
integer (mpiknd), intent(out):: y(0:)
integer (mpiknd) t1(0:mpnw+6+nmxa), t2(0:mpnw+6+nmxa), t3(0:mpnw+6+nmxa), &
  t4(0:mpnw+6+nmxa)

!  End of declaration

if (nn < -nmax .or. nn >= 0) then
  write (mpldb, 1) -nmax
1 format ('*** MPPOLYLOGNEG: N is <',i6,' or n >= 0.'/ &
  'For n >= 0, call mppolylogpos or polylog_pos.')
  call mpabrt (111)
endif

nna = abs (nn)
i1 = 1000000000
i2 = 1000000000

do k = 1, nna
  i1 = min (i1, mpspacer (arr(0:na+5,k)))
  i2 = min (i2, mpwprecr (arr(0:na+5,k)))
enddo

call mpmdc (arr(0:na+5,1), d1, n1, mpnw)
d1 = d1 * 2.d0 ** n1
call mpmdc (arr(0:na+5,nna), d2, n2, mpnw)
d2 = d2 * 2.d0 ** n2

if (mpnw < 4 .or. i1 < mpnw + 6 .or. i2 < mpnw .or. d1 /= 1.d0 .or. d2 /= 1.d0 &
  .or. mpspacer (x) < mpnw + 4 .or. mpspacer (y) < mpnw + 6) then
  write (mpldb, 2)
2 format ('*** MPPOLYLOGNEG: uninitialized or inadequately sized arrays'/ &
  'Call mppolylogini or polylog_ini to initialize array. See documentation.')
  call mpabrt (99)
endif

mpnw1 = mpnw + 1
call mpinitwds (t1, mpnw + nmxa + 1)
call mpinitwds (t2, mpnw + nmxa + 1)
call mpinitwds (t3, mpnw + nmxa + 1)
call mpinitwds (t4, mpnw + nmxa + 1)

!  na = abs (n)
!  if (x < 0.d0) then
!    call mpbinmd (arr((na+1)/2), d1, n1)
!    nwds1 = nwds + (n1 + 1) / mpnbt + 1
!  else
!    nwds1 = nwds
!  endif

if (mpsigntr (x) < 0) then 
  call mpmdc (arr(0:mpnw+5,(nna+1)/2), d1, n1, mpnw1)
  mpnw1 = mpnw1 + (n1 + 1) / mpnbt + 1
endif

!  t1 = mpreal (x, nwds1)
!  t2 = t1

call mpeq (x, t1, mpnw1)
call mpeq (t1, t2, mpnw1)

!  do k = 2, na
!    t1 = x * t1
!    t2 = t2 + arr(k) * t1
!  enddo

do k = 2, nna
  call mpmul (x, t1, t3, mpnw1)
  call mpeq (t3, t1, mpnw1)
  call mpmul (arr(0:mpnw+5,k), t1, t3, mpnw1)
  call mpadd (t2, t3, t4, mpnw1)
  call mpeq (t4, t2, mpnw1)
enddo

!  polylogr = mpreal (t2 / (1.d0 - x) ** (na + 1), nwds)

call mpdmc (1.d0, 0, t3, mpnw1)
call mpsub (t3, x, t4, mpnw1)
call mpnpwr (t4, nna + 1, t3, mpnw1)
call mpdiv (t2, t3, t4, mpnw1)
call mproun (t4, mpnw)
call mpeq (t4, y, mpnw)

return
end subroutine mppolylogneg

subroutine mppolylogpos (nn, x, y, mpnw)

!   This returns polylog (nn, x) for the case nn >= 0.

implicit none
integer, intent(in):: nn, mpnw
integer ic1, itrmax, k, mpnw1, nn1, n1
parameter (itrmax = 1000000)
integer (mpiknd), intent(in):: x(0:)
integer (mpiknd), intent(out):: y(0:)
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), &
  t5(0:mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps (0:9)

!  End of declaration

if (nn < 0) then
  write (mpldb, 1) nn
1 format ('*** MPPOLYLOGPOS: N is less than zero.'/ &
  'For negative n, call mppolylogneg or polylog_neg. See documentation.')
  call mpabrt (111)
endif

if (mpnw < 4 .or. mpspacer (x) < mpnw + 4 .or. mpspacer (y) < mpnw + 6) then
  write (mpldb, 2)
2 format ('*** MPPOLYLOGPOS: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

mpnw1 = mpnw + 1
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)

if (nn == 0) then
!  polylogr = x / (1.d0 - x)

  call mpdmc (1.d0, 0, t1, mpnw1)
  call mpsub (t1, x, t2, mpnw1)
  call mpdiv (x, t2, t3, mpnw1)
  call mproun (t3, mpnw)
  call mpeq (t3, y, mpnw)
else
!  t1 = x
!  t2 = x

  call mpeq (x, t1, mpnw1)
  call mpeq (x, t2, mpnw1)

  do k = 2, itrmax
!    t2 = x * t2
!    t3 = t2 / mpreal (dble (k), nwds) ** n
!    t1 = t1 + t3
!    if (abs (t3 / x) < eps) goto 100

    call mpmul (x, t2, t3, mpnw1)
    call mpeq (t3, t2, mpnw1)
    call mpdmc (dble (k), 0, t3, mpnw1)
    call mpnpwr (t3, nn, t4, mpnw1)
    call mpdiv (t2, t4, t3, mpnw1)
    call mpadd (t1, t3, t4, mpnw1)
    call mpeq (t4, t1, mpnw1)

!    if (t3(2) == 0 .or. t3(3) < x(3) - mpnw1) goto 100

    call mpabs (t3, tc1, 4)
    call mpmul (eps, t1, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic1, 4)
    if (ic1 <= 0) goto 100 
  enddo

  write (6, 3)
3 format ('MPPOLYLOGPOS: Loop end error')
  call mpabrt (111)

100 continue

!  polylogr = t1

  call mproun (t1, mpnw)
  call mpeq (t1, y, mpnw)
endif

return
end subroutine mppolylogpos

subroutine mpzetar (ss, zz, mpnw)

!   This returns the zeta function of an MPR argument SS using an algorithm
!   due to Peter Borwein.

implicit none
integer, intent(in):: mpnw
integer i, ic1, iss, itrmax, j, mpnw1, n, n1, n2
real (mprknd) dfrac, d1, d2
parameter (itrmax = 1000000, dfrac = 1.d0+ceiling(mpdpw))
integer (mpiknd), intent(in):: ss(0:)
integer (mpiknd), intent(out):: zz(0:)
integer (mpiknd) f1(0:mpnw+6), s(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6), &
  t3(0:mpnw+6), t4(0:mpnw+6), t5(0:mpnw+6), tn(0:mpnw+6), tt(0:mpnw+6), &
  tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)
real (mprknd) sgn

!  End of declaration

if (mpnw < 4 .or. mpspacer (ss) < mpnw + 4 .or. mpspacer (zz) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPZETAR: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

mpnw1 = mpnw + 1
call mpinitwds (f1, mpnw1)
call mpinitwds (s, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (tn, mpnw1)
call mpinitwds (tt, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)

!   Set f1 = 1.

call mpdmc (1.d0, 0, f1, mpnw1)
call mpcpr (ss, f1, ic1, mpnw1)
call mpinfr (ss, t1, t2, mpnw1)

if (ic1 == 0) then
  write (mpldb, 2)
2 format ('*** MPZETAR: argument is 1')
  call mpabrt (63)
elseif (mpsigntr (t2) == 0) then

!   The argument is an integer value. Call mpzetaintr instead.

  call mpmdc (ss, d1, n1, mpnw)
  iss = d1 * 2.d0**n1
  call mpzetaintr (iss, t1, mpnw)
  goto 200
elseif (mpsigntr (ss) < 0) then

!   If arg < 0, compute zeta(1-ss), and later apply Riemann's formula.

  call mpsub (f1, ss, tt, mpnw1)
else
  call mpeq (ss, tt, mpnw1)
endif

!  Check if argument is large enough that computing with definition is faster.

d1 = mpnbt * mpnw * log (2.d0) / log (2.d0 * mpnbt * mpnw / 3.d0)
call mpmdc (tt, d2, n2, mpnw1)
d2 = d2 * 2.d0 ** n2

if (d2 > d1) then

!   Evaluate the infinite series.

!  t1 = mpreal (1.d0, mpnw)

  call mpdmc (1.d0, 0, t1, mpnw1)

  do i = 2, itrmax

!    t2 = mpreal (dble (i), mpnw) ** tt
!    t3 = 1.d0 / t2
!    t1 = t1 + t3

    call mpdmc (dble (i), 0, t4, mpnw1)
    call mppower (t4, tt, t2, mpnw1)
    call mpdiv (f1, t2, t3, mpnw1)
    call mpadd (t1, t3, t2, mpnw1)
    call mpeq (t2, t1, mpnw1)

!    if (t3(2) == 0 .or. t3(3) < - mpnw1) goto 200

    call mpabs (t3, tc1, 4)
    call mpmul (eps, t1, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic1, 4)
    if (ic1 <= 0) goto 200
  enddo
  
  write (mpldb, 3) 1, itrmax
3 format ('*** MPZETAR: iteration limit exceeded',2i10)
  call mpabrt (101)
endif  

n = dfrac * mpnw1

! tn = mpreal (2.d0, mpnw) ** n
! t1 = - tn
! t2 = mpreal (0.d0, mpnw)
! s = mpreal (0.d0, mpnw)

call mpdmc (2.d0, 0, t1, mpnw1)
call mpnpwr (t1, n, tn, mpnw1)
call mpneg (tn, t1, mpnw1)
call mpdmc (0.d0, 0, t2, mpnw1)
call mpdmc (0.d0, 0, s, mpnw1)

sgn = 1.d0

do j = 0, 2 * n - 1
!  t3 = mpreal (dble (j + 1), mpnw) ** tt
!  s = s + sgn * t1 / t3

  call mpdmc (dble (j + 1), 0, t4, mpnw1)
  call mppower (t4, tt, t3, mpnw1)
  call mpdiv (t1, t3, t4, mpnw1)
  call mpmuld (t4, sgn, t5, mpnw1)
  call mpadd (s, t5, t4, mpnw1)
  call mpeq (t4, s, mpnw1)

  sgn = - sgn

  if (j .lt. n - 1) then
    call mpdmc (0.d0, 0, t2, mpnw1)
  elseif (j .eq. n - 1) then
    call mpdmc (1.d0, 0, t2, mpnw1)    
  else
!     t2 = t2 * dble (2 * n - j) / dble (j + 1 - n)

     call mpmuld (t2, dble (2 * n - j), t3, mpnw1)
     call mpdivd (t3, dble (j + 1 - n), t2, mpnw1)
  endif
!  t1 = t1 + t2

  call mpadd (t1, t2, t3, mpnw1)
  call mpeq (t3, t1, mpnw1)
enddo

! t1 = - s / (tn * (1.d0 - mpreal (2.d0, mpnw) ** (1.d0 - tt)))

call mpsub (f1, tt, t3, mpnw1)
call mpdmc (2.d0, 0, t2, mpnw1)
call mppower (t2, t3, t4, mpnw1)
call mpsub (f1, t4, t2, mpnw1)
call mpmul (tn, t2, t3, mpnw1)
call mpdiv (s, t3, t1, mpnw1)
! t1(2) = - t1(2)
call mpneg (t1, t2, mpnw1)
call mpeq (t2, t1, mpnw1)

!   If original argument was negative, apply Riemann's formula.

if (mpsigntr (ss) < 0) then
  call mpgammar (tt, t3, mpnw1)
  call mpmul (t1, t3, t2, mpnw1)
  call mpmul (mppicon, tt, t1, mpnw1)
  call mpmuld (t1, 0.5d0, t3, mpnw1)
  call mpcssnr (t3, t4, t5, mpnw1)
  call mpmul (t2, t4, t1, mpnw1)
  call mpmuld (mppicon, 2.d0, t2, mpnw1)
  call mppower (t2, tt, t3, mpnw1)
  call mpdiv (t1, t3, t2, mpnw1)
  call mpmuld (t2, 2.d0, t1, mpnw1)
endif

200 continue

! zetapbr = t1

call mproun (t1, mpnw)
call mpeq (t1, zz, mpnw)
return
end subroutine mpzetar

subroutine mpzetaintr (iss, zz, mpnw)

!   This returns the zeta function of the integer argument ISS using an algorithm
!   due to Peter Borwein.

implicit none
integer, intent(in):: mpnw
integer i, ic1, itrmax, j, mpnw1, n, n1, n2, itt
real (mprknd) dfrac, d1, d2
parameter (itrmax = 1000000, dfrac = 1.d0+ceiling(mpdpw))
integer, intent(in):: iss
integer (mpiknd), intent(out):: zz(0:)
integer (mpiknd) f1(0:mpnw+6), s(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6), &
  t3(0:mpnw+6), t4(0:mpnw+6), t5(0:mpnw+6), tn(0:mpnw+6), &
  tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)
real (mprknd) sgn

!  End of declaration

if (mpnw < 4 .or. mpspacer (zz) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPZETAINTR: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

mpnw1 = mpnw + 1
call mpinitwds (f1, mpnw1)
call mpinitwds (s, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (tn, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)

!   Set f1 = 1.

call mpdmc (1.d0, 0, f1, mpnw1)

if (iss == 1) then
  write (mpldb, 2)
2 format ('*** MPZETAINTR: argument is 1')
  call mpabrt (63)
elseif (iss == 0) then

!   Argument is zero -- result is -1/2.

  call mpdmc (-0.5d0, 0, t1, mpnw1)
  goto 200
elseif (iss < 0) then

!   If argument is a negative even integer, the result is zero.

  if (mod (iss, 2) == 0) then
    call mpdmc (0.d0, 0, t1, mpnw1)
    goto 200
  endif

!   Otherwise if arg < 0, compute zeta(1-is), and later apply Riemann's formula.

  itt = 1 - iss
else
  itt = iss
endif

!  Check if argument is large enough that computing with definition is faster.

d1 = mpnbt * mpnw * log (2.d0) / log (2.d0 * mpnbt * mpnw / 3.d0)

if (itt > d1) then

!   Evaluate the infinite series.

!  t1 = mpreal (1.d0, mpnw)

  call mpdmc (1.d0, 0, t1, mpnw1)

  do i = 2, itrmax

!    t2 = mpreal (dble (i), mpnw) ** tt
!    t3 = 1.d0 / t2
!    t1 = t1 + t3

    call mpdmc (dble (i), 0, t4, mpnw1)
    call mpnpwr (t4, itt, t2, mpnw1)
    call mpdiv (f1, t2, t3, mpnw1)
    call mpadd (t1, t3, t2, mpnw1)
    call mpeq (t2, t1, mpnw1)

!    if (t3(2) == 0 .or. t3(3) < - mpnw1) goto 200

    call mpabs (t3, tc1, 4)
    call mpmul (eps, t1, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic1, 4)
    if (ic1 <= 0) goto 200
  enddo
  
  write (mpldb, 3) 1, itrmax
3 format ('*** MPZETAINTR: iteration limit exceeded',2i10)
  call mpabrt (101)
endif  

n = dfrac * mpnw1

! tn = mpreal (2.d0, mpnw) ** n
! t1 = - tn
! t2 = mpreal (0.d0, mpnw)
! s = mpreal (0.d0, mpnw)

call mpdmc (2.d0, 0, t1, mpnw1)
call mpnpwr (t1, n, tn, mpnw1)
call mpneg (tn, t1, mpnw1)
call mpdmc (0.d0, 0, t2, mpnw1)
call mpdmc (0.d0, 0, s, mpnw1)

sgn = 1.d0

do j = 0, 2 * n - 1
!  t3 = mpreal (dble (j + 1), mpnw) ** tt
!  s = s + sgn * t1 / t3

  call mpdmc (dble (j + 1), 0, t4, mpnw1)
  call mpnpwr (t4, itt, t3, mpnw1)
  call mpdiv (t1, t3, t4, mpnw1)
  call mpmuld (t4, sgn, t5, mpnw1)
  call mpadd (s, t5, t4, mpnw1)
  call mpeq (t4, s, mpnw1)

  sgn = - sgn

  if (j .lt. n - 1) then
    call mpdmc (0.d0, 0, t2, mpnw1)
  elseif (j .eq. n - 1) then
    call mpdmc (1.d0, 0, t2, mpnw1)    
  else
!     t2 = t2 * dble (2 * n - j) / dble (j + 1 - n)

     call mpmuld (t2, dble (2 * n - j), t3, mpnw1)
     call mpdivd (t3, dble (j + 1 - n), t2, mpnw1)
  endif
!  t1 = t1 + t2

  call mpadd (t1, t2, t3, mpnw1)
  call mpeq (t3, t1, mpnw1)
enddo

! t1 = - s / (tn * (1.d0 - mpreal (2.d0, mpnw) ** (1.d0 - tt)))

call mpdmc (2.d0, 0, t2, mpnw1)
call mpnpwr (t2, 1 - itt, t4, mpnw1)
call mpsub (f1, t4, t2, mpnw1)
call mpmul (tn, t2, t3, mpnw1)
call mpdiv (s, t3, t1, mpnw1)
! t1(2) = - t1(2)
call mpneg (t1, t2, mpnw1)
call mpeq (t2, t1, mpnw1)

!   If original argument was negative, apply Riemann's formula.

if (iss < 0) then
  call mpdmc (1.d0, 0, t3, mpnw1)
  do i = 1, itt - 1
    call mpmuld (t3, dble (i), t4, mpnw1)
    call mpeq (t4, t3, mpnw1)
  enddo

  call mpmul (t1, t3, t2, mpnw1)
  call mpmuld (mppicon, dble (itt), t1, mpnw1)
  call mpmuld (t1, 0.5d0, t3, mpnw1)
  call mpcssnr (t3, t4, t5, mpnw1)
  call mpmul (t2, t4, t1, mpnw1)
  call mpmuld (mppicon, 2.d0, t2, mpnw1)
  call mpnpwr (t2, itt, t3, mpnw1)
  call mpdiv (t1, t3, t2, mpnw1)
  call mpmuld (t2, 2.d0, t1, mpnw1)
endif

200 continue

! zetapbr = t1

call mproun (t1, mpnw)
call mpeq (t1, zz, mpnw)
return
end subroutine mpzetaintr

subroutine mpzetaemr (nb1, nb2, berne, s, z, mpnw)

!  This evaluates the Riemann zeta function, using the combination of 
!  the definition formula (for large s), and an Euler-Maclaurin scheme
!  (see formula 25.2.9 of the DLMF.  The array berne contains precomputed
!  Bernoulli numbers.  Its dimensions must be as shown below. NB2 must be
!  at least as great as the precision level in decimal digits.

implicit none
integer, intent(in):: nb1, nb2, mpnw
integer i, i1, i2, ia, ic1, itrmax, k, mpnw1, n1, n2, na, nn
real (mprknd) dfrac, dlogb, d1, d2
parameter (itrmax = 1000000, dfrac = 8.5d0, dlogb = 33.27106466687737d0)
integer (mpiknd), intent(in):: berne(0:nb1+5,nb2), s(0:)
integer (mpiknd), intent(out):: z(0:)
integer (mpiknd) t0(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), &
  t4(0:mpnw+6), t5(0:mpnw+6), t6(0:mpnw+6), t7(0:mpnw+6), t8(0:mpnw+6), &
  t9(0:mpnw+6), tt(0:mpnw+6), f1(0:mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)

! End of declaration

if (mpnw < 4 .or. mpspacer (s) < mpnw + 4 .or. mpspacer (z) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPZETAEMR: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

i = 0
k = 0
mpnw1 = mpnw + 1
call mpinitwds (t0, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (t6, mpnw1)
call mpinitwds (t7, mpnw1)
call mpinitwds (t8, mpnw1)
call mpinitwds (t9, mpnw1)
call mpinitwds (tt, mpnw1)
call mpinitwds (f1, mpnw1)

call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)

!   Check if argument is 1 -- undefined.

call mpdmc (1.d0, 0, t0, mpnw1)
call mpcpr (s, t0, ic1, mpnw1)
if (ic1 == 0) then
  write (mpldb, 2)
2 format ('*** MPZETAEMR: argument is 1')
  call mpabrt (63)
endif

!   Check if berne array has been initialized.

i1 = 1000000000
i2 = 1000000000

do k = 1, nb2
  i1 = min (i1, mpspacer (berne(0:nb1+5,k)))
  i2 = min (i2, mpwprecr (berne(0:nb1+5,k)))
enddo

call mpmdc (berne(0:na+5,1), d1, n1, mpnw)
d1 = d1 * 2.d0 ** n1

if (i1 < mpnw + 6 .or. i2 < mpnw .or. abs (d1 - 1.d0 / 6.d0) > mprdfz .or. &
  nb2 < int (mpndpw * mpnw)) then
  write (mpldb, 3) int (mpndpw * mpnw)
3 format ('*** MPZETAEMR: Array of even Bernoulli coefficients must be initialized'/ &
   'with at least',i8,' entries.')
  call mpabrt (62)
endif

!   Set f1 = 1.

call mpdmc (1.d0, 0, f1, mpnw1)

!   Check if argument is zero.  If so, result is - 1/2.

if (mpsigntr (s) == 0) then
  call mpdmc (-0.5d0, 0, t1, mpnw)
  goto 200
endif

!   Check if argument is negative.

if (mpsigntr (s) < 0) then

!   Check if argument is a negative even integer.  If so, the result is zero.

  call mpmuld (s, 0.5d0, t1, mpnw1)
  call mpinfr (t1, t2, t3, mpnw1)
  if (mpsigntr (t3) == 0) then
    call mpdmc (0.d0, 0, t1, mpnw1)
    goto 200
  endif

!   Otherwise compute zeta(1-s), and later apply the reflection formula.

  call mpsub (f1, s, tt, mpnw1)
else
  call mpeq (s, tt, mpnw1)
endif

!  Check if argument is large enough that computing with definition is faster.

! if (tt .gt. mpreald (dlogb * mpnw / log (32.d0 * mpnw), mpnw)) then

d1 = dlogb * mpnw1 / log (32.d0 * mpnw1)
call mpmdc (tt, d2, n2, mpnw1)
d2 = d2 * 2.d0**n2
if (d2 > d1) then

!  t1 = mpreal (1.d0, mpnw)

  call mpdmc (1.d0, 0, t1, mpnw1)

  do i = 2, itrmax

!    t2 = mpreal (dble (i), mpnw) ** tt

    call mpdmc (dble (i), 0, t4, mpnw1)
    call mppower (t4, tt, t2, mpnw1)
    
!    t3 = 1.d0 / t2

    call mpdiv (f1, t2, t3, mpnw1)

!    t1 = t1 + t3

    call mpadd (t1, t3, t2, mpnw1)
    call mpeq (t2, t1, mpnw1)

!    if (t3(2) == 0 .or. t3(3) < - mpnw) goto 200

    call mpabs (t3, tc1, 4)
    call mpmul (eps, t1, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic1, 4)
    if (ic1 <= 0) goto 200
  enddo
  
  write (mpldb, 4) 1, itrmax
4 format ('*** MPZETAEMR: iteration limit exceeded',2i10)
  call mpabrt (101)
endif  

! t0 = mpreal (1.d0, mpnw)

call mpdmc (1.d0, 0, t0, mpnw1)

nn = dfrac * mpnw1

do k = 2, nn
!  t1 = mpreal (dble (k), mpnw) ** tt

  call mpdmc (dble (k), 0, t2, mpnw1)
  call mppower (t2, tt, t1, mpnw1)
  
!  t0 = t0 + 1.d0 / t1

  call mpdiv (f1, t1, t2, mpnw1)
  call mpadd (t0, t2, t3, mpnw1)
  call mpeq (t3, t0, mpnw1)
enddo

! t0 = t0 + dble (nn) / (t1 * (tt - 1.d0)) - 0.5d0 / t1

call mpdmc (dble (nn), 0, t2, mpnw1)
call mpsub (tt, f1, t3, mpnw1)
call mpmul (t1, t3, t4, mpnw1)
call mpdiv (t2, t4, t3, mpnw1)
call mpadd (t0, t3, t2, mpnw1)
call mpdmc (0.5d0, 0, t3, mpnw1)
call mpdiv (t3, t1, t4, mpnw1)
call mpsub (t2, t4, t0, mpnw1)

! t3 = tt

call mpeq (tt, t3, mpnw1)

! t2 = t3 / (12.d0 * dble (nn) * t1)

call mpmuld (t1, 12.d0 * dble (nn), t4, mpnw1)
call mpdiv (t3, t4, t2, mpnw1)

! t5 = dble (nn) * t1

call mpmuld (t1, dble (nn), t5, mpnw1)

! t9 = dble (nn) ** 2

call mpdmc (dble (nn), 0, t6, mpnw1)
call mpmul (t6, t6, t9, mpnw1)

do k = 2, min (nb2, itrmax)
!  t3 = t3 * (tt + dble (2*k - 2)) * (tt + dble (2*k - 3)) / &
!         (dble (2 * k - 1) * dble (2 * k - 2))

  call mpdmc (dble (2 * k - 2), 0, t4, mpnw1)
  call mpadd (tt, t4, t6, mpnw1)
  call mpdmc (dble (2 * k - 3), 0, t7, mpnw1)
  call mpadd (tt, t7, t8, mpnw1)
  call mpmul (t6, t8, t7, mpnw1)
  call mpmul (t3, t7, t4, mpnw1)
  call mpdmc (dble (2 * k - 1), 0, t6, mpnw1)
  call mpdmc (dble (2 * k - 2), 0, t7, mpnw1)
  call mpmul (t6, t7, t8, mpnw1)
  call mpdiv (t4, t8, t3, mpnw1)
  
!  t5 = t5 * t9

  call mpmul (t5, t9, t6, mpnw1)
  call mpeq (t6, t5, mpnw1)
  
!  t7 = t3 * berne(k) / (dble (2 * k) * t5)

  call mpmul (t3, berne(0:nb1+5,k), t4, mpnw1)
  call mpmuld (t5, dble (2 * k), t6, mpnw1)
  call mpdiv (t4, t6, t7, mpnw1)

!  t2 = t2 + t7

  call mpadd (t2, t7, t4, mpnw1)
  call mpeq (t4, t2, mpnw1)

!  if (t7(2) == 0 .or. t7(3) < t2(3) - mpnw) goto 110

  call mpabs (t7, tc1, 4)
  call mpmul (eps, t2, tc3, 4)
  call mpabs (tc3, tc2, 4)
  call mpcpr (tc1, tc2, ic1, 4)
  if (ic1 <= 0) goto 110
enddo

write (mpldb, 3) 2, min (nb2, itrmax)
call mpabrt (101)

110 continue

! zetaem = t0 + t2

call mpadd (t0, t2, t1, mpnw1)

!   If original argument was negative, apply the reflection formula.

if (mpsigntr (s) < 0) then
  call mpgammar (tt, t3, mpnw1)
  call mpmul (t1, t3, t2, mpnw1)
  call mpmul (mppicon, tt, t1, mpnw1)
  call mpmuld (t1, 0.5d0, t3, mpnw1)
  call mpcssnr (t3, t4, t5, mpnw1)
  call mpmul (t2, t4, t1, mpnw1)
  call mpmuld (mppicon, 2.d0, t2, mpnw1)
  call mppower (t2, tt, t3, mpnw1)
  call mpdiv (t1, t3, t2, mpnw1)
  call mpmuld (t2, 2.d0, t1, mpnw1)
endif

200 continue

call mproun (t1, mpnw)
call mpeq (t1, z, mpnw)

return
end subroutine mpzetaemr

end module mpfune
