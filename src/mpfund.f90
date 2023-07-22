!*****************************************************************************

!  MPFUN20-Fort: A thread-safe arbitrary precision package with special functions
!  Transcendental function module (module MPFUND)

!  Revision date:  25 Mar 2023

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

!  DESCRIPTION OF THIS MODULE (MPFUND):
!    This module contains subroutines for basic transcendental functions,
!    including routines for cos, sin, inverse cos/sin, hyperbolic cos/sin,
!    and inverse hyperbolic cos, sin, as well as routines to compute pi and log(2).

module mpfund
use mpfuna
use mpfunb
use mpfunc

implicit none

contains

subroutine mpagmr (a, b, c, mpnw)

!   This performs the arithmetic-geometric mean (AGM) iterations on A and B.
!   The AGM algorithm is as follows: Set a_0 = a and b_0 = b, then iterate

!    a_{k+1} = (a_k + b_k)/2
!    b_{k+1} = sqrt (a_k * b_k)

!   until convergence (i.e., until a_k = b_k to available precision).
!   The result is returned in C.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:), b(0:)
integer (mpiknd), intent(out):: c(0:)
integer, parameter:: itrmx = 100
integer j, mpnw1
integer (mpiknd) s0(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6)

!  End of declaration.

if (mpnw < 4 .or. c(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPAGMR: uninitialized or inadequately sized arrays')
  call mpabrt ( 401)
endif

if (mpsgn (a) <= 0 .or. mpsgn (b) <= 0) then
  write (mpldb, 2)
2 format ('*** MPAGMR: argument <= 0')
  call mpabrt ( 402)
endif

s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7
s3(0) = mpnw + 7
mpnw1 = mpnw + 1
call mpeq (a, s1, mpnw1)
call mpeq (b, s2, mpnw1)

do j = 1, itrmx
  call mpadd (s1, s2, s0, mpnw1)
  call mpmuld (s0, 0.5d0, s3, mpnw1)
  call mpmul (s1, s2, s0, mpnw1)
  call mpsqrt (s0, s2, mpnw1)
 call mpeq (s3, s1, mpnw1)

!   Check for convergence.

  call mpsub (s1, s2, s0, mpnw1)
  if (s0(2) == 0 .or. s0(3) < 1 - mpnw1) goto 100
enddo

write (mpldb, 3)
3 format ('*** MPAGMR: Iteration limit exceeded.')
call mpabrt ( 402)

100 continue

call mproun (s1, mpnw)
call mpeq (s1, c, mpnw)

return
end subroutine mpagmr

subroutine mpang (x, y, a, mpnw)

!   This computes the MPR angle A subtended by the MPR pair (X, Y) considered as
!   a point in the x-y plane. This is more useful than an arctan or arcsin
!   routine, since it places the result correctly in the full circle, i.e.
!   -Pi < A <= Pi. Pi must be precomputed to at least MPNW words precision
!   and the stored in the array in module MPFUNA.

!   The Taylor series for Arcsin converges much more slowly than that of Sin,
!   so this routine solves Cos (a) = x or Sin (a) = y using one of the
!   following Newton iterations, both of which converge to a:

!     z_{k+1} = z_k - [x - Cos (z_k)] / Sin (z_k)
!     z_{k+1} = z_k + [y - Sin (z_k)] / Cos (z_k)

!   The first is selected if Abs (x) <= Abs (y); otherwise the second is used.
!   These iterations are performed with a maximum precision level MPNW that
!   is dynamically changed, approximately doubling with each iteration.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: x(0:), y(0:)
integer (mpiknd), intent(out):: a(0:)
integer, parameter:: nit = 3
real (mprknd), parameter:: cl2 = 1.4426950408889633d0, cpi = 3.1415926535897932d0
integer iq, ix, iy, k, kk, mpnw1, mq, nx, ny, n1, n2
real (mprknd) d1, t1, t2, t3
integer (mpiknd) s0(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6), &
  s4(0:mpnw+6), s5(0:mpnw+6)

! End of declaration

if (mpnw < 4 .or. a(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPANG: uninitialized or inadequately sized arrays')
  call mpabrt ( 403)
endif

ix = sign (int (1, mpiknd), x(2))
nx = min (int (abs (x(2))), mpnw)
iy = sign (int (1, mpiknd), y(2))
ny = min (int (abs (y(2))), mpnw)
mpnw1 = mpnw + 1
s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7
s3(0) = mpnw + 7
s4(0) = mpnw + 7
s5(0) = mpnw + 7

!   Check if both X and Y are zero.

if (nx == 0 .and. ny == 0) then
  write (mpldb, 2)
2 format ('*** MPANG: Both arguments are zero.')
  call mpabrt ( 404)
endif

!   Check if Pi has been precomputed.

call mpmdc (mppicon, d1, n1, mpnw1)
d1 = d1 * 2.d0**n1
if (abs (d1 - cpi) > mprdfz .or. mpnw1 > mppicon(1)) then
  write (mpldb, 3) mpnw1
3 format ('*** MPANG: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 405)
endif

!   Check if one of X or Y is zero.

if (nx == 0) then
  call mpeq (mppicon, s0, mpnw1)
  if (iy > 0) then
    call mpmuld (s0, 0.5d0, a, mpnw)
  else
    call mpmuld (s0, -0.5d0, a, mpnw)
  endif
  goto 120
elseif (ny == 0) then
  if (ix > 0) then
    a(1) = mpnw
    a(2) = 0
    a(3) = 0
    a(4) = 0
    a(5) = 0
  else
    call mpeq (s0, a, mpnw)
  endif
  goto 120
endif

!   Determine the least integer MQ such that 2 ^ MQ >= MPNW.

t1 = mpnw1
mq = cl2 * log (t1) + 2.d0 - mprdfz

!   Normalize x and y so that x^2 + y^2 = 1.

call mpmul (x, x, s0, mpnw1)
call mpmul (y, y, s1, mpnw1)
call mpadd (s0, s1, s2, mpnw1)
call mpsqrt (s2, s3, mpnw1)
call mpdiv (x, s3, s1, mpnw1)
call mpdiv (y, s3, s2, mpnw1)

!   Compute initial approximation of the angle.

call mpmdc (s1, t1, n1, mpnw1)
call mpmdc (s2, t2, n2, mpnw1)
n1 = max (n1, -mpnbt)
n2 = max (n2, -mpnbt)
t1 = t1 * 2.d0 ** n1
t2 = t2 * 2.d0 ** n2
t3 = atan2 (t2, t1)
call mpdmc (t3, 0, s5, mpnw1)

!   The smaller of x or y will be used from now on to measure convergence.
!   This selects the Newton iteration (of the two listed above) that has the
!   largest denominator.

if (abs (t1) <= abs (t2)) then
  kk = 1
  call mpeq (s1, s0, mpnw1)
else
  kk = 2
  call mpeq (s2, s0, mpnw1)
endif

mpnw1 = 4
iq = 0

!   Perform the Newton-Raphson iteration described above with a dynamically
!   changing precision level MPNW (one greater than powers of two).

do k = 1, mq
  mpnw1 = min (2 * mpnw1 - 2, mpnw)

100  continue

  call mpcssnr (s5, s1, s2, mpnw1)

  if (kk == 1) then
    call mpsub (s0, s1, s3, mpnw1)
    call mpdiv (s3, s2, s4, mpnw1)
    call mpsub (s5, s4, s1, mpnw1)
  else
    call mpsub (s0, s2, s3, mpnw1)
    call mpdiv (s3, s1, s4, mpnw1)
    call mpadd (s5, s4, s1, mpnw1)
  endif
  call mpeq (s1, s5, mpnw1)

  if (k == mq - nit .and. iq == 0) then
    iq = 1
    goto 100
  endif
enddo

!   Restore original precision level.

call mproun (s5, mpnw)
call mpeq (s5, a, mpnw)

120 continue

return
end subroutine mpang

subroutine mpcagm (a, b, c, mpnw)

!   This performs the arithmetic-geometric mean (AGM) iterations on A and B
!   for MPC arguments A and B.
!   The AGM algorithm is as follows: Set a_0 = a and b_0 = b, then iterate

!    a_{k+1} = (a_k + b_k)/2
!    b_{k+1} = sqrt (a_k * b_k)

!   until convergence (i.e., until a_k = b_k to available precision).
!   The result is returned in C.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:), b(0:)
integer (mpiknd), intent(out):: c(0:)
integer, parameter:: itrmx = 100
integer j, la, lb, lc, mp7, mpnw1
integer (mpiknd) s0(0:2*mpnw+13), s1(0:2*mpnw+13), s2(0:2*mpnw+13), &
  s3(0:2*mpnw+13)

!  End of declaration.

la = a(0)
lb = b(0)
lc = c(0)
if (mpnw < 4 .or.  c(0) < mpnw + 6 .or. c(lc) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCAGM: uninitialized or inadequately sized arrays')
  call mpabrt ( 406)
endif

mp7 = mpnw + 7
s0(0) = mp7
s0(mp7) = mp7
s1(0) = mp7
s1(mp7) = mp7
s2(0) = mp7
s2(mp7) = mp7
s3(0) = mp7
s3(mp7) = mp7
mpnw1 = mpnw + 1
call mpceq (a, s1, mpnw1)
call mpceq (b, s2, mpnw1)

do j = 1, itrmx
  call mpcadd (s1, s2, s0, mpnw1)
  call mpmuld (s0, 0.5d0, s3, mpnw1)
  call mpmuld (s0(mp7:), 0.5d0, s3(mp7:), mpnw1)
  call mpcmul (s1, s2, s0, mpnw1)
  call mpcsqrt (s0, s2, mpnw1)
  call mpceq (s3, s1, mpnw1)
  call mpcsub (s1, s2, s0, mpnw1)

!   Check for convergence.

  if ((s0(2) == 0 .or. s0(3) < 1 - mpnw1) .and. &
    (s0(mp7+2) == 0 .or. s0(mp7+3) < 1 - mpnw1)) goto 100
enddo

write (mpldb, 2)
2 format ('*** MPCAGM: Iteration limit exceeded.')
call mpabrt ( 407)

100 continue

call mproun (s1, mpnw)
call mproun (s1(mp7:), mpnw)
call mpceq (s1, c, mpnw)

return
end subroutine mpcagm

subroutine mpcexp (a, b, mpnw)

!   This computes Exp[A], for MPC A.

!   The formula is:  E^a1 * (Cos[a2] + I * Sin[a2]), where a1 and a2 are
!   the real and imaginary parts of A.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:)
integer (mpiknd), intent(out):: b(0:)
integer la, lb, mpnw1
integer (mpiknd) s0(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6), &
  s4(0:mpnw+6)

! End of declaration

la = a(0)
lb = b(0)
if (mpnw < 4 .or. b(0) < mpnw + 6 .or. b(lb) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCEXP: uninitialized or inadequately sized arrays')
  call mpabrt ( 408)
endif

mpnw1 = mpnw + 1
s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7
s3(0) = mpnw + 7
s4(0) = mpnw + 7

call mpexp (a, s0, mpnw1)
call mpcssnr (a(la:), s1, s2, mpnw1)
call mpmul (s0, s1, s3, mpnw1)
call mpmul (s0, s2, s4, mpnw1)

call mproun (s3, mpnw)
call mproun (s4, mpnw)
call mpeq (s3, b, mpnw)
call mpeq (s4, b(lb:), mpnw)

100 continue

return
end subroutine mpcexp

subroutine mpclog (a, b, mpnw)

!   This computes Log[A], for MPC A.

!   The formula is:  1/2 * Log[r] + I * Theta, where r = a1^2 + a2^2,
!   Theta is the angle corresponding to (a1, a2), and a1 and a2 are the
!   real and imaginary parts of A.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:)
integer (mpiknd), intent(out):: b(0:)
integer la, lb, mpnw1
integer (mpiknd) s0(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6)

! End of declaration

la = a(0)
lb = b(0)
if (mpnw < 4 .or. b(0) < mpnw + 6 .or. b(lb) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCLOG: uninitialized or inadequately sized arrays')
  call mpabrt ( 409)
endif

mpnw1 = mpnw + 1
s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7
s3(0) = mpnw + 7

call mpmul (a, a, s0, mpnw1)
call mpmul (a(la:), a(la:), s1, mpnw1)
call mpadd (s0, s1, s2, mpnw1)
call mplog (s2, s3, mpnw1)
call mpmuld (s3, 0.5d0, s0, mpnw1)
call mpang (a, a(la:), s1, mpnw1)

call mproun (s0, mpnw)
call mproun (s1, mpnw)
call mpeq (s0, b, mpnw)
call mpeq (s1, b(lb:), mpnw)

100 continue

return
end subroutine mpclog

subroutine mpcpowcc (a, b, c, mpnw)

!   This computes A^B, where A and B are MPC.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:), b(0:)
integer (mpiknd), intent(out):: c(0:)
integer la, lb, lc, l3
integer (mpiknd) s1(0:2*mpnw+11), s2(0:2*mpnw+11)

! End of declaration

la = a(0)
lb = b(0)
lc = c(0)
if (mpnw < 4 .or. c(0) < mpnw + 6 .or. c(lc) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCPOWCC: uninitialized or inadequately sized arrays')
  call mpabrt ( 410)
endif

l3 = mpnw + 6
s1(0) = l3
s1(l3) = l3
s2(0) = l3
s2(l3) = l3
call mpclog (a, s1, mpnw)
call mpcmul (s1, b, s2, mpnw)
call mpcexp (s2, c, mpnw)

return
end subroutine mpcpowcc

subroutine mpcpowcr (a, b, c, mpnw)

!   This computes A^B, where A is MPC and B is MPR.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:), b(0:)
integer (mpiknd), intent(out):: c(0:)
integer la, lc, l3
integer (mpiknd) s1(0:2*mpnw+11), s2(0:2*mpnw+11)

! End of declaration

la = a(0)
lc = c(0)
if (mpnw < 4 .or. c(0) < mpnw + 6 .or. c(lc) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCPOWCR: uninitialized or inadequately sized arrays')
  call mpabrt ( 411)
endif

l3 = mpnw + 6
s1(0) = l3
s1(l3) = l3
s2(0) = l3
s2(l3) = l3
call mpclog (a, s1, mpnw)
call mpmul (b, s1, s2, mpnw)
call mpmul (b, s1(l3:), s2(l3:), mpnw)
call mpcexp (s2, c, mpnw)

return
end subroutine mpcpowcr

subroutine mpcpowrc (a, b, c, mpnw)

!   This computes A^B, where A is MPR and and B is MPC.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:), b(0:)
integer (mpiknd), intent(out):: c(0:)
integer lb, lc, l3
integer (mpiknd) s1(0:2*mpnw+11), s2(0:2*mpnw+11)

! End of declaration

lb = b(0)
lc = c(0)
if (mpnw < 4 .or. c(0) < mpnw + 6 .or. c(lc) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCPOWRC: uninitialized or inadequately sized arrays')
  call mpabrt ( 412)
endif

l3 = mpnw + 6
s1(0) = l3
s2(0) = l3
s2(l3) = l3
call mplog (a, s1, mpnw)
call mpmul (s1, b, s2, mpnw)
call mpmul (s1, b(lb:), s2(l3:), mpnw)
call mpcexp (s2, c, mpnw)

return
end subroutine mpcpowrc

subroutine mpcsshr (a, x, y, mpnw)

!   This computes the hyperbolic cosine and sine of the MPR number A and
!   returns the two MPR results in X and Y, respectively. If the argument
!   is very close to zero, a Taylor series is used; otherwise this routine
!   calls mpexp.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:)
integer (mpiknd), intent(out):: x(0:), y(0:)
integer, parameter:: itrmx = 1000000
integer j, mpnw1, mpnw2
real (mprknd) t2
integer (mpiknd) f(0:9), s0(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6)

! End of declaration

if (mpnw < 4 .or. x(0) < mpnw + 6 .or. y(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCSSHR: uninitialized or inadequately sized arrays')
  call mpabrt ( 413)
endif

s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7
s3(0) = mpnw + 7
mpnw1 = mpnw + 1
f(0) = 9
f(1) = mpnw
f(2) = 1
f(3) = 0
f(4) = 1
f(5) = 0
f(6) = 0

!   If argument is very small, compute the sinh using a Taylor series.
!   This avoids accuracy loss that otherwise occurs by using exp.

if (a(3) < -1) then
  call mpeq (a, s0, mpnw1)
  call mpmul (s0, s0, s2, mpnw1)
  mpnw2 =  mpnw1

!   The working precision used to compute each term can be linearly reduced
!   as the computation proceeds.

  do j = 1, itrmx
    t2 = (2.d0 * j) * (2.d0 * j + 1.d0)
    call mpmul (s2, s1, s3, mpnw2)
    call mpdivd (s3, t2, s1, mpnw2)
    call mpadd (s1, s0, s3, mpnw1)
    call mpeq (s3, s0, mpnw1)

!   Check for convergence of the series, and adjust working precision
!   for the next term.

    if (s1(2) == 0 .or. s1(3) < s0(3) - mpnw1) goto 110
    mpnw2 = min (max (mpnw1 + int (s1(3) - s0(3)) + 1, 4), mpnw1)
  enddo

  write (mpldb, 4)
4 format ('*** MPCSSHR: Iteration limit exceeded.')
  call mpabrt ( 414)

110 continue

  call mpmul (s0, s0, s2, mpnw1)
  call mpadd (f, s2, s3, mpnw1)
  call mpsqrt (s3, s1, mpnw1)
  call mproun (s1, mpnw)
  call mpeq (s1, x, mpnw)
  call mproun (s0, mpnw)
  call mpeq (s0, y, mpnw)
else
  call mpexp (a, s0, mpnw1)
  call mpdiv (f, s0, s1, mpnw1)
  call mpadd (s0, s1, s2, mpnw1)
  call mpmuld (s2, 0.5d0, s3, mpnw1)
  call mproun (s3, mpnw)
  call mpeq (s3, x, mpnw)
  call mpsub (s0, s1, s2, mpnw1)
  call mpmuld (s2, 0.5d0, s3, mpnw1)
  call mproun (s3, mpnw)
  call mpeq (s3, y, mpnw)
endif

100 continue

return
end subroutine mpcsshr

subroutine mpcssnr (a, x, y, mpnw)

!   This computes the cosine and sine of the MPR number A and returns the
!   two MPR results in X and Y, respectively. Pi must be precomputed to
!   at least MPNW words precision and the stored in the array MPPICON in
!   module MPFUNA.

!   This routine uses the conventional Taylor series for Sin (s):

!   Sin (s) =  s - s^3 / 3! + s^5 / 5! - s^7 / 7! ...

!   where the argument S has been reduced to (-pi, pi). To further
!   accelerate the series, the reduced argument is divided by 2^NQ, where NQ
!   is computed as int (sqrt (0.5d0 * N)), where N is the precision in bits.
!   After convergence of the series, the double-angle formulas for cos are
!   applied NQ times.

implicit none
integer, intent(in):: mpnw
real (mprknd) d1, t1, t2
integer (mpiknd), intent(in):: a(0:)
integer (mpiknd), intent(out):: x(0:), y(0:)
integer, parameter:: itrmx = 1000000
real (mprknd), parameter:: cpi = 3.1415926535897932d0
integer is, i1, j, mpnw1, mpnw2, na, nq, n1
integer (mpiknd) f1(0:8), f2(0:8),s0(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6), &
  s3(0:mpnw+6), s4(0:mpnw+6), s5(0:mpnw+6)

! End of declaration

if (mpnw < 4 .or. x(0) < mpnw + 6 .or. y(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCSSNR: uninitialized or inadequately sized arrays')
  call mpabrt ( 415)
endif

na = min (int (abs (a(2))), mpnw)
if (na == 0) then
  x(1) = mpnw
  x(2) = 1
  x(3) = 0
  x(4) = 1
  x(5) = 0
  x(6) = 0
  y(1) = mpnw
  y(2) = 0
  y(3) = 0
  y(4) = 0
  y(5) = 0
  goto 120
endif

s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7
s3(0) = mpnw + 7
s4(0) = mpnw + 7
s5(0) = mpnw + 7
mpnw1 = mpnw + 1

!   Set f1 = 1 and f2 = 1/2.

f1(0) = 9
f1(1) = mpnw
f1(2) = 1
f1(3) = 0
f1(4) = 1
f1(5) = 0
f1(6) = 0
f2(0) = 9
f2(1) = mpnw
f2(2) = 1
f2(3) = -1
f2(4) = 0.5d0 * mpbdx
f2(5) = 0
f2(6) = 0

!   Check if Pi and Sqrt(2)/2 have been precomputed in data statements to the
!   requested precision.

call mpmdc (mppicon, d1, n1, mpnw1)
d1 = d1 * 2.d0**n1
if (abs (d1 - cpi) > mprdfz .or. mpnw1 > mppicon(1)) then
  write (mpldb, 2) mpnw1
2 format ('*** MPCSSNR: Pi must be precomputed to precision',i9,' words).'/ &
  'See documentation for details.')
  call mpabrt ( 416)
endif

!   Check if argument is too large to compute meaningful cos/sin values.

call mpmdc (a, t1, n1, mpnw)
if (n1 >= mpnbt * (mpnw - 1)) then
  write (mpldb, 3)
3 format ('*** MPCSSNR: argument is too large to compute cos or sin.')
  call mpabrt ( 417)
endif

!   Reduce to between - Pi and Pi.

call mpmuld (mppicon, 2.d0, s0, mpnw1)
call mpdiv (a, s0, s1, mpnw1)
call mpnint (s1, s2, mpnw1)
call mpmul (s0, s2, s4, mpnw1)
call mpsub (a, s4, s3, mpnw1)

!   Check if reduced argument is zero. If so then cos = 1 and sin = 0.

if (s3(2) == 0) then
  s0(1) = mpnw1
  s0(2) = 1
  s0(3) = 0
  s0(4) = 1
  s0(5) = 0
  s0(6) = 0
  s1(1) = mpnw1
  s1(2) = 0
  s1(3) = 0
  s1(4) = 0
  s1(5) = 0
  goto 115
endif

!   Determine nq to scale reduced argument, then divide by 2^nq.
!   If reduced argument is very close to zero, then nq = 0.

if (s3(3) >= -1) then
  nq = int (sqrt (0.5d0 * mpnw1 * mpnbt))
else
  nq = 0
endif

call mpdivd (s3, 2.d0**nq, s0, mpnw1)
call mpeq (s0, s1, mpnw1)

!   Compute the sin of the reduced argument of s1 using a Taylor series.

call mpmul (s0, s0, s2, mpnw1)
mpnw2 =  mpnw1
is = s0(2)

!   The working precision used to compute each term can be linearly reduced
!   as the computation proceeds.

do i1 = 1, itrmx
  t2 = - (2.d0 * i1) * (2.d0 * i1 + 1.d0)
  call mpmul (s2, s1, s3, mpnw2)
  call mpdivd (s3, t2, s1, mpnw2)
  call mpadd (s1, s0, s3, mpnw1)
  call mpeq (s3, s0, mpnw1)

!   Check for convergence of the series, and adjust working precision
!   for the next term.

  if (s1(2) == 0 .or. s1(3) < s0(3) - mpnw1) goto 110
  mpnw2 = min (max (mpnw1 + int (s1(3) - s0(3)) + 1, 4), mpnw1)
enddo

write (mpldb, 4)
4 format ('*** MPCSSNR: Iteration limit exceeded.')
call mpabrt ( 418)

110 continue

if (nq > 0) then

!   Apply the formula cos(2*x) = 2*cos^2(x) - 1 NQ times to produce
!   the cosine of the reduced argument, except that the first iteration is
!   cos(2*x) = 1 - 2*sin^2(x), since we have computed sin(x) above.
!   Note that these calculations are performed as 2 * (cos^2(x) - 1/2) and
!   2 * (1/2 - sin^2(x)), respectively, to avoid loss of precision.

  call mpmul (s0, s0, s4, mpnw1)
  call mpsub (f2, s4, s5, mpnw1)
  call mpmuld (s5, 2.d0, s0, mpnw1)

  do j = 2, nq
    call mpmul (s0, s0, s4, mpnw1)
    call mpsub (s4, f2, s5, mpnw1)
    call mpmuld (s5, 2.d0, s0, mpnw1)
  enddo

!   Compute sin of result and correct sign.

  call mpmul (s0, s0, s4, mpnw1)
  call mpsub (f1, s4, s5, mpnw1)
  call mpsqrt (s5, s1, mpnw1)
  if (is < 1) s1(2) = - s1(2)
else

!   In case nq = 0, compute cos of result.

  call mpeq (s0, s1, mpnw1)
  call mpmul (s0, s0, s4, mpnw1)
  call mpsub (f1, s4, s5, mpnw1)
  call mpsqrt (s5, s0, mpnw1)
endif

115 continue

!   Restore original precision level.

call mproun (s0, mpnw)
call mproun (s1, mpnw)
call mpeq (s0, x, mpnw)
call mpeq (s1, y, mpnw)

120 continue

return
end subroutine mpcssnr

subroutine mpegammaq (egamma, mpnw)

!   This computes Euler's gamma to available precision (MPNW mantissa words).
!   The algorithm is the following, due to Brent and MacMillan:

!   Given a desired precision D digits, select n1 = 0.6*D and n2 = 2.2*D
!   (approximately). Set A_0 = -log(n1); B_0 = 1; U_0 = A_0; V_0 = 1.
!   Then iterate for k = 1, 2, ..., n2:
!     B_k = B_{k-1} * n1^2 / k^2
!     A_k = (A_{k-1} * n1^2/k + B_k) / k
!     U_k = U_{k-1} + A_k
!     V_k = V_{k-1} + B_k
!   Then set gamma = U_k / V_k.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(out):: egamma(0:)
integer k, mpnw1, ndp, n1, n2
integer (mpiknd) aa(0:mpnw+6), bb(0:mpnw+6), uu(0:mpnw+6), vv(0:mpnw+6), &
  s0(0:mpnw+6), s1(0:mpnw+6)

! End of declaration

if (mpnw < 4 .or. egamma(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPEGAMMAQ: uninitialized or inadequately sized arrays')
  call mpabrt ( 419)
endif

mpnw1 = mpnw + 1
ndp = (mpnw - 1) * mpdpw
n1 = nint (0.6d0 * ndp)
n2 = nint (2.2d0 * ndp)

!   Unless precision is very high, just copy egamma from table.

if (mpnw1 <= mpegammacon(1)) then
  call mpeq (mpegammacon, egamma, mpnw)
  goto 110
endif

!   Initialize aa, bb, uu, vv, s0, s1.

aa(0) = mpnw1 + 6
bb(0) = mpnw1 + 6
uu(0) = mpnw1 + 6
vv(0) = mpnw1 + 6
s0(0) = mpnw1 + 6
s1(0) = mpnw1 + 6

! aa = - log (mpreal (dble (n1), nwds))

call mpdmc (dble (n1), 0, s0, mpnw1)
call mplog (s0, aa, mpnw1)
aa(2) = - aa(2)

! bb = mpreal (1.d0, nwds)
! uu = aa
! vv = bb

call mpdmc (1.d0, 0, bb, mpnw1)
call mpeq (aa, uu, mpnw1)
call mpeq (bb, vv, mpnw1)

do k = 1, n2

!   bb = bb * dble (n1)**2 / dble (k)**2

  call mpmuld (bb, dble (n1)**2, s0, mpnw1)
  call mpdivd (s0, dble (k)**2, bb, mpnw1)

!  aa = (aa * dble (n1)**2 / dble (k) + bb) / dble (k)

  call mpmuld (aa, dble (n1)**2, s0, mpnw1)
  call mpdivd (s0, dble (k), s1, mpnw1)
  call mpadd (s1, bb, s0, mpnw1)
  call mpdivd (s0, dble (k), aa, mpnw1)

!  uu = uu + aa
!  vv = vv + bb

  call mpadd (uu, aa, s0, mpnw1)
  call mpeq (s0, uu, mpnw1)
  call mpadd (vv, bb, s0, mpnw1)
  call mpeq (s0, vv, mpnw1)

!  if (aa%mpr(3) < uu%mpr(3) - nwds .and. bb%mpr(3) < vv%mpr(3) - nwds) &
!    goto 100

  if (aa(3) < uu(3) - mpnw1 .and. bb(3) < vv(3) - mpnw1) goto 100
enddo

100 continue

! egammax = uu / vv

call mpdiv (uu, vv, s0, mpnw1)

!   Restore original precision level.

call mproun (s0, mpnw)
call mpeq (s0, egamma, mpnw)

110 continue

return
end subroutine mpegammaq

subroutine mpexp (a, b, mpnw)

!   This computes the exponential function of the MPR number A and returns
!   the MPR result in B. Log(2) must be precomputed to at least MPNW words
!   precision and the stored in the array MPLOG2CON in module MPFUNA.

!   This routine uses a modification of the Taylor series for Exp (t):

!   Exp (t) =  (1 + r + r^2 / 2! + r^3 / 3! + r^4 / 4! ...) ^ q * 2 ^ n

!   where the argument T has been reduced to within the closest factor of Log(2).
!   To further accelerate the series, the reduced argument is divided by 2^NQ.
!   After convergence of the series, the result is squared NQ times. NQ is
!   set to bb^0.4, where bb is the number of bits of precision.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:)
integer (mpiknd), intent(out):: b(0:)
integer, parameter:: itrmx = 1000000
real (mprknd), parameter:: cl2 = 0.69314718055994530d0
integer j, mpnw1, mpnw2, nq, nz, n1
real (mprknd) d1, t1, t2
integer (mpiknd) f(0:8), s0(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6)

! End of declaration

if (mpnw < 4 .or. b(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPEXP: uninitialized or inadequately sized arrays')
  call mpabrt ( 420)
endif

call mpmdc (a, t1, n1, mpnw)

!   Check for overflows and underflows.

if (n1 > 30) then
  if (t1 > 0.d0) then
    write (mpldb, 2)
2   format ('*** MPEXP: Argument is too large.')
    call mpabrt ( 421)
  else
    b(1) = mpnw
    b(2) = 0
    b(3) = 0
    b(4) = 0
    b(5) = 0
    goto 130
  endif
endif

t1 = t1 * 2.d0 ** n1
if (abs (t1) > 1488522236.d0) then
  if (t1 > 0) then
    write (mpldb, 2)
    call mpabrt ( 422)
  else
    b(1) = mpnw
    b(2) = 0
    b(3) = 0
    b(4) = 0
    b(5) = 0
    goto 130
  endif
endif

s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7
s3(0) = mpnw + 7
mpnw1 = mpnw + 1

!   Set f1 = 1.

f(0) = 9
f(1) = mpnw1
f(2) = 1
f(3) = 0
f(4) = 1
f(5) = 0
f(6) = 0

!   Check if Log(2) has been precomputed.

call mpmdc (mplog2con, d1, n1, mpnw1)
d1 = d1 * 2.d0**n1
if (abs (d1 - cl2) > mprdfz .or. mpnw1 > mplog2con(1)) then
  write (mpldb, 3) mpnw1
3 format ('*** MPEXP: Log(2) must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 423)
endif

!   Compute the reduced argument A' = A - Log(2) * Nint [A / Log(2)]. Save
!   NZ = Nint [A / Log(2)] for correcting the exponent of the final result.

call mpdiv (a, mplog2con, s0, mpnw1)
call mpnint (s0, s1, mpnw1)
call mpmdc (s1, t1, n1, mpnw1)
nz = nint (t1 * 2.d0 ** n1)
call mpmul (mplog2con, s1, s2, mpnw1)
call mpsub (a, s2, s0, mpnw1)

!   Check if the reduced argument is zero.

if (s0(2) == 0) then
  s0(1) = mpnw1
  s0(2) = 1
  s0(3) = 0
  s0(4) = 1
  s0(5) = 0
  s0(6) = 0
  goto 120
endif

!   Divide the reduced argument by 2 ^ NQ.

nq = max (nint (dble (mpnw * mpnbt)** 0.4d0), 1)
call mpdivd (s0, 2.d0**nq, s1, mpnw1)

!   Compute Exp using the usual Taylor series.

call mpeq (f, s2, mpnw1)
call mpeq (f, s3, mpnw1)
mpnw2 =  mpnw1

!   The working precision used to compute each term can be linearly reduced
!   as the computation proceeds.

do j = 1, itrmx
  t2 = dble (j)
  call mpmul (s2, s1, s0, mpnw2)
  call mpdivd (s0, t2, s2, mpnw2)
  call mpadd (s3, s2, s0, mpnw1)
  call mpeq (s0, s3, mpnw1)

!   Check for convergence of the series, and adjust working precision
!   for the next term.

  if (s2(2) == 0 .or. s2(3) < s0(3) - mpnw1) goto 100
  mpnw2 = min (max (mpnw1 + int (s2(3) - s0(3)) + 1, 4), mpnw1)
enddo

write (mpldb, 4)
4 format ('*** MPEXP: Iteration limit exceeded.')
call mpabrt ( 424)

100 continue

!   Raise to the (2 ^ NQ)-th power.

do j = 1, nq
  call mpmul (s0, s0, s1, mpnw1)
  call mpeq (s1, s0, mpnw1)
enddo

!   Multiply by 2 ^ NZ.

120 continue

call mpdmc (1.d0, nz, s2, mpnw1)
call mpmul (s0, s2, s1, mpnw1)

!   Restore original precision level.

call mproun (s1, mpnw)
call mpeq (s1, b, mpnw)

130 continue

return
end subroutine mpexp

subroutine mpinitran (mpnw)

!   This routine computes pi, log(2) and egamma, and stores this data in the
!   proper arrays in module MPFUNA. MPNW is the largest precision level
!   (in words) that will be subsequently required for this run at the user level.

implicit none
integer, intent(in):: mpnw
integer nwds, nwds6

!  End of declaration

!   Add three words to mpnw, since many of the routines in this module
!   increase the working precision level by one word upon entry.

nwds = mpnw + 3

!  Compute pi, log(2) and sqrt(2)/2.

nwds6 = nwds + 6
mplog2con(0) = nwds6
mplog2con(1) = 0
mppicon(0) = nwds6
mppicon(1) = 0
mpegammacon(0) = nwds6
mpegammacon(1) = 0

call mppiq (mppicon, nwds)
call mplog2q (mplog2con, nwds)
call mpegammaq (mpegammacon, nwds - 3)

return
end subroutine mpinitran

subroutine mplog (a, b, mpnw)

!   This computes the natural logarithm of the MPR number A and returns the MPR
!   result in B.

!   The Taylor series for Log converges much more slowly than that of Exp.
!   Thus this routine does not employ Taylor series (except if the argument
!   is extremely close to 1), but instead computes logarithms by solving
!   Exp (b) = a using the following Newton iteration:

!     x_{k+1} = x_k + [a - Exp (x_k)] / Exp (x_k)

!   These iterations are performed with a maximum precision level MPNW that
!   is dynamically changed, approximately doubling with each iteration.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:)
integer (mpiknd), intent(out):: b(0:)
integer, parameter:: itrmax = 1000000, nit = 3
real (mprknd), parameter:: alt = 0.693147180559945309d0, cl2 = 1.4426950408889633d0, &
  rtol = 0.5d0**7
integer ia, iq, is, i1, k, mpnw1, mq, na, n1
real (mprknd) st, tol, t1, t2
integer (mpiknd) f1(0:8), s0(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6), &
  s3(0:mpnw+6), s4(0:mpnw+6)

! End of declaration

if (mpnw < 4 .or. b(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPLOG: uninitialized or inadequately sized arrays')
  call mpabrt ( 425)
endif

ia = sign (int (1, mpiknd), a(2))
na = min (int (abs (a(2))), mpnw)

if (ia < 0 .or. na == 0) then
  write (mpldb, 2)
2 format ('*** MPLOG: Argument is less than or equal to zero.')
  call mpabrt ( 426)
endif

!   Check if input is exactly one.

if (a(2) == 1 .and. a(3) == 0 .and. a(4) == 1) then
  b(1) = mpnw
  b(2) = 0
  b(3) = 0
  b(4) = 0
  b(5) = 0
  goto 130
endif

mpnw1 = mpnw + 1
s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7
s3(0) = mpnw + 7
s4(0) = mpnw + 7

f1(0) = 9
f1(1) = mpnw1
f1(2) = 1
f1(3) = 0
f1(4) = 1
f1(5) = 0
f1(6) = 0

!   If the argument is sufficiently close to 1, employ a Taylor series.

call mpsub (a, f1, s0, mpnw1)

if (s0(2) == 0 .or. s0(3) <= min (-2.d0, - rtol * mpnw1)) then
  call mpeq (s0, s1, mpnw1)
  call mpeq (s1, s2, mpnw1)
  i1 = 1
  is = 1
  tol = s0(3) - mpnw1

  do i1 = 2, itrmax
    is = - is
    st = is * i1
    call mpmul (s1, s2, s3, mpnw1)
    call mpeq (s3, s2, mpnw1)
    call mpdivd (s3, st, s4, mpnw1)
    call mpadd (s0, s4, s3, mpnw1)
    call mpeq (s3, s0, mpnw1)
    if (s4(2) == 0 .or. s4(3) < tol) goto 120
  enddo

  write (mpldb, 3) itrmax
3 format ('*** MPLOG: Iteration limit exceeded',i10)
  call mpabrt ( 427)
endif

!   Determine the least integer MQ such that 2 ^ MQ >= MPNW.

t2 = mpnw
mq = cl2 * log (t2) + 2.d0 - mprdfz

!   Compute initial approximation of Log (A).

call mpmdc (a, t1, n1, mpnw)
t1 = log (t1) + n1 * alt
call mpdmc (t1, 0, s3, mpnw)
mpnw1 = 4
iq = 0

!   Perform the Newton-Raphson iteration described above with a dynamically
!   changing precision level MPNW (one greater than powers of two).

do k = 0, mq
  if (k > 1) mpnw1 = min (2 * mpnw1 - 2, mpnw) + 1

110  continue

  call mpexp (s3, s0, mpnw1)
  call mpsub (a, s0, s1, mpnw1)
  call mpdiv (s1, s0, s2, mpnw1)
  call mpadd (s3, s2, s1, mpnw1)
  call mpeq (s1, s3, mpnw1)
  if (k == mq - nit .and. iq == 0) then
    iq = 1
    goto 110
  endif
enddo

!   Restore original precision level.

120 continue

call mproun (s3, mpnw)
call mpeq (s3, b, mpnw)

130 continue

return
end subroutine mplog

subroutine mplog2q (alog2, mpnw)

!   This computes log(2) to mpnw words precision, using an algorithm due to Salamin
!   and Brent:  Select n > 2^m, where m is the number of bits of desired precision
!   precision in the result. Then

!   Log(2) = Pi / [2 AGM (1, 4/x)]

!   Where AGM (a, b) denotes the arithmetic-geometric mean:  Set a_0 = a and
!   b_0 = b, then iterate
!    a_{k+1} = (a_k + b_k)/2
!    b_{k+1} = sqrt (a_k * b_k)
!   until convergence (i.e., until a_k = b_k to available precision).

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(out):: alog2(0:)
real (mprknd), parameter:: cpi = 3.141592653589793238d0
integer mpnw1, n, n1, n48
real (mprknd) d1
integer (mpiknd) f1(0:8), f4(0:8), s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6), &
  s4(0:mpnw+6)

! End of declaration

if (mpnw < 4 .or. alog2(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPLOG2Q: uninitialized or inadequately sized arrays')
  call mpabrt ( 428)
endif

!   Define sections of the scratch array.

s1(0) = mpnw + 7
s2(0) = mpnw + 7
s3(0) = mpnw + 7
s4(0) = mpnw + 7
mpnw1 = mpnw + 1

!   Unless precision is very high, just copy log2 from table.

if (mpnw1 <= mplog2con(1)) then
  call mpeq (mplog2con, alog2, mpnw)
  goto 100
endif

!   Check if Pi has been precomputed.

call mpmdc (mppicon, d1, n1, mpnw)
d1 = d1 * 2.d0**n1
if (abs (d1 - cpi) > mprdfz .or. mppicon(1) < mpnw) then
  write (mpldb, 2) mpnw
2 format ('*** MPLOG2Q: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 429)
endif

!   Define sections of the scratch array.

s1(0) = mpnw + 7
s2(0) = mpnw + 7
s3(0) = mpnw + 7
s4(0) = mpnw + 7
mpnw1 = mpnw + 1

!   Set f1 = 1.

f1(0) = 9
f1(1) = mpnw1
f1(2) = 1
f1(3) = 0
f1(4) = 1
f1(5) = 0
f1(6) = 0

!   Set f4 = 4.

f4(0) = 9
f4(1) = mpnw1
f4(2) = 1
f4(3) = 0
f4(4) = 4
f4(5) = 0
f4(6) = 0

!   Set s4 to 2^(n/2), where n is the number of bits desired. n48 = n/mpnbt.
!   Note that this value can be directly set in the first few words of s4,
!   avoiding explicit exponentiation.

n = mpnbt * (mpnw1 / 2 + 2)
n48 = n / mpnbt

s4(1) = mpnw1
s4(2) = 1
s4(3) = n48
s4(4) = 1
s4(5) = 0
s4(6) = 0

!   Perform AGM iterations.

call mpeq (f1, s1, mpnw1)
call mpdiv (f4, s4, s2, mpnw1)
call mpagmr (s1, s2, s3, mpnw1)

!   Set Log(2) = Pi / (2 * N * S3), where S3 is the limit of the AGM iterations.

call mpmuld (s3, 2.d0 * n, s1, mpnw1)
call mpdiv (mppicon, s1, s2, mpnw1)
call mproun (s2, mpnw)
call mpeq (s2, alog2, mpnw)

100 continue

return
end subroutine mplog2q

subroutine mppiq (pi, mpnw)

!   This computes Pi to available precision (MPNW mantissa words).
!   The algorithm that is used for computing Pi, which is due to Salamin
!   and Brent, is as follows:

!   Set  A_0 = 1,  B_0 = 1/Sqrt(2)  and  D_0 = Sqrt(2) - 1/2.

!   Then from k = 1 iterate the following operations:

!   A_k = 0.5 * (A_{k-1} + B_{k-1})
!   B_k = Sqrt (A_{k-1} * B_{k-1})
!   D_k = D_{k-1} - 2^k * (A_k - B_k) ^ 2

!   Then  P_k = (A_k + B_k) ^ 2 / D_k  converges quadratically to Pi.
!   In other words, each iteration approximately doubles the number of correct
!   digits, providing all iterations are done with the maximum precision.
!   The constant cl2 (below) = 1 / log(2) (DP approximation).

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(out):: pi(0:)
real (mprknd), parameter:: cl2 = 1.4426950408889633d0
integer k, mpnw1, mq
real (mprknd) t1
integer (mpiknd) f(0:8), s0(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6), &
  s4(0:mpnw+6)

! End of declaration

if (mpnw < 4 .or. pi(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPPIQ: uninitialized or inadequately sized arrays')
  call mpabrt ( 430)
endif

s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7
s3(0) = mpnw + 7
s4(0) = mpnw + 7
mpnw1 = mpnw + 1

!   Unless precision is very high, just copy pi from table.

if (mpnw1 <= mppicon(1)) then
  call mpeq (mppicon, pi, mpnw)
  goto 100
endif

!   Determine the number of iterations required for the given precision level.
!   This formula is good only for this Pi algorithm.

t1 = mpnw1 * log10 (dble (mpbdx))
mq = cl2 * (log (t1) - 1.d0) + 1.d0

!   Initialize as above.

s0(1) = mpnw
s0(2) = 1
s0(3) = 0
s0(4) = 1
s0(5) = 0
s0(6) = 0
f(0) = 9
f(1) = mpnw1
f(2) = 1
f(3) = 0
f(4) = 2
f(5) = 0
f(6) = 0
call mpsqrt (f, s2, mpnw1)
call mpmuld (s2, 0.5d0, s1, mpnw1)
f(3) = -1
f(4) = 0.5d0 * mpbdx
call mpsub (s2, f, s4, mpnw1)

!   Perform iterations as described above.

do k = 1, mq
  call mpadd (s0, s1, s2, mpnw1)
  call mpmul (s0, s1, s3, mpnw1)
  call mpsqrt (s3, s1, mpnw1)
  call mpmuld (s2, 0.5d0, s0, mpnw1)
  call mpsub (s0, s1, s2, mpnw1)
  call mpmul (s2, s2, s3, mpnw1)
  t1 = 2.d0 ** k
  call mpmuld (s3, t1, s2, mpnw1)
  call mpsub (s4, s2, s3, mpnw1)
  call mpeq (s3, s4, mpnw1)
enddo

!   Complete computation.

call mpadd (s0, s1, s2, mpnw1)
call mpmul (s2, s2, s2, mpnw1)
call mpdiv (s2, s4, s2, mpnw1)
call mpeq (s2, s0, mpnw1)

!   Restore original precision level.

call mproun (s0, mpnw)
call mpeq (s0, pi, mpnw)

100 continue

return
end subroutine mppiq

subroutine mppower (a, b, c, mpnw)

!   This computes C = A ^ B, where A, B and C are MPR. It first checks if
!   B is the quotient of two integers up to 10^7 in size, in which case it
!   calls MPNPWR and MPNRTR. Otherwise it calls MPLOG and MPEXP.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:), b(0:)
integer (mpiknd), intent(Out):: c(0:)
real (mprknd), parameter:: mprxx = 5.d-10
integer i, n1
real (mprknd) a1, a2, a3, a4, a5, a6, q1, t0, t1, t2, t3
integer (mpiknd) s0(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6)

!  End of declaration

if (mpnw < 4 .or. c(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPPOWER: uninitialized or inadequately sized arrays')
  call mpabrt ( 431)
endif

!   Check if A < 0 or A = 0 and B < 0 (error); if A = 1 or B = 0 then
!   result = 1; if B = 1 then result = A.

if (a(2) < 0 .or. a(2) == 0 .and. B(2) < 0) then
  write (mpldb, 2)
2 format ('*** MPPOWER: A^B, where A < 0 or A = 0 and B < 0.')
  call mpabrt ( 432)
elseif ((a(2) == 1 .and. a(3) == 0 .and. a(4) == 1) &
  .or. b(2) == 0) then
  c(1) = mpnw
  c(2) = 1
  c(3) = 0
  c(4) = 1
  c(5) = 0
  c(6) = 0
  goto 200
elseif (b(2) == 1 .and. b(3) == 0 .and. b(4) == 1) then
  call mpeq (a, c, mpnw)
  goto 200
endif

s0(0) = mpnw + 6
s1(0) = mpnw + 6
s2(0) = mpnw + 6

!   Check if B is rational using the extended Euclidean algorithm in DP.

call mpmdc (b, t1, n1, mpnw)

if (n1 >= -mpnbt .and. n1 <= mpnbt) then
  t0 = abs (t1 * 2.d0**n1)
  t1 = max (t0, 1.d0)
  t2 = min (t0, 1.d0)
  a1 = 1.d0
  a2 = 0.d0
  a3 = 0.d0
  a4 = 1.d0

  do i = 1, 20
    q1 = aint (t1 / t2)
    a5 = a1 - q1 * a3
    a6 = a2 - q1 * a4
    t3 = t2
    t2 = t1 - q1 * t2
    t1 = t3
    a1 = a3
    a2 = a4
    a3 = a5
    a4 = a6
    if (t2 < mprxx) goto 100
  enddo

  goto 110
endif

100 continue

a3 = abs (a3)
a4 = abs (a4)

!  If b = a3/a4 or a4/a3 (except for sign) or then call mpnpwr and mpnrtr.

if (abs (t0 - a3 / a4) / t0 < mprdfz) then
  a3 = sign (a3, dble (b(2)))
  call mpdmc (a3, 0, s0, mpnw)
  call mpdmc (a4, 0, s1, mpnw)
  call mpdiv (s0, s1, s2, mpnw)
  call mpsub (b, s2, s0, mpnw)
  if (s0(2) == 0 .or. s0(3) < b(3) + 1 - mpnw) then
    call mpnpwr (a, int (a3), s0, mpnw)
    call mpnrtr (s0, int (a4), c, mpnw)
    goto 200
  endif
elseif (abs (t0 - a4 / a3) / t0 < mprdfz) then
  a4 = sign (a4, dble (b(2)))
  call mpdmc (a4, 0, s0, mpnw)
  call mpdmc (a3, 0, s1, mpnw)
  call mpdiv (s0, s1, s2, mpnw)
  call mpsub (b, s2, s0, mpnw)
  if (s0(2) == 0 .or. s0(3) < b(3) + 1 - mpnw) then
    call mpnpwr (a, int (a4), s0, mpnw)
    call mpnrtr (s0, int (a3), c, mpnw)
    goto 200
  endif
endif

110 continue

!  Call mplog and mpexp.

call mplog (a, s0, mpnw)
call mpmul (s0, b, s1, mpnw)
call mpexp (s1, c, mpnw)

200 continue

return
end subroutine mppower

end module mpfund

