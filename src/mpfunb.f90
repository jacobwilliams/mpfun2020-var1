!*****************************************************************************

!  MPFUN20-Fort: A thread-safe arbitrary precision package with special functions
!  Basic function module (module MPFUNB)

!  Revision date:  24 May 2023

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

!  DESCRIPTION OF THIS MODULE (MPFUNB):
!    This module contains routines for: add, subtract, multiply, divide;
!    comparison; double/multi conversion; double/multi multiplication/division;
!    integer and fractional parts; nearest integer; nth power; nth root;
!    rounding and normalization, random numbers, square roots, conversions
!    to/from quad precision (if available), and routines to support FFT-based
!    multiplication and Newton division (used for very high precision). The
!    routines in this package are not intended to be directly called by the
!    user; the high-level language interface modules should be used instead.

module mpfunb
use mpfuna
use mpmask
implicit none

contains

subroutine mpabrt (ier)

!   This routine terminates execution. Users may wish to replace the
!   default STOP with a call to a system routine that provides a traceback.

implicit none
integer, intent(in):: ier

! End of declaration

write (mpldb, 1) ier
1 format ('*** MPABRT: Execution terminated, error code =',i4)
stop
end subroutine mpabrt

subroutine mpabs (ra, rb, mpnw)

!   This routine sets rb = absolute value of ra.

implicit none
integer (mpiknd), intent(in):: ra(0:)
integer (mpiknd), intent(out):: rb(0:)
integer mpnw

! End of declaration

call mpeq (ra, rb, mpnw)
rb(2) = min (int (abs (ra(2))), mpnw)
return
end subroutine mpabs

subroutine mpadd (a, b, c, mpnw)

!   This routine adds MPR numbers A and B to yield C.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:), b(0:)
integer (mpiknd), intent(out):: c(0:)
integer i, ia, ib, idb, ish, ixa, ixb, ixd, m1, m2, m3, m4, m5, na, nb, &
  nd, nsh
integer (mpiknd) d(0:mpnw+6)

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. b(0) < abs (b(2)) + 4 .or. &
  c(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPADD: uninitialized or inadequately sized arrays')
  call mpabrt ( 201)
endif

ia = sign (int (1, mpiknd), a(2))
ib = sign (int (1, mpiknd), b(2))
na = min (int (abs (a(2))), mpnw)
nb = min (int (abs (b(2))), mpnw)

!   Check for zero inputs.

if (na == 0) then

!   A is zero -- the result is B.

  c(1) = mpnw
  c(2) = sign (nb, ib)

  do i = 2, nb + 3
    c(i+1) = b(i+1)
  enddo

  c(nb+4) = 0
  c(nb+5) = 0
  goto 100
elseif (nb == 0) then

!   B is zero -- the result is A.

  c(1) = mpnw
  c(2) = sign (na, ia)

  do i = 2, na + 3
    c(i+1) = a(i+1)
  enddo

  c(na+4) = 0
  c(na+5) = 0
  goto 100
endif

if (ia == ib) then
  idb = 1
else
  idb = -1
endif
ixa = a(3)
ixb = b(3)
ish = ixa - ixb

if (ish >= 0) then

!   A has greater exponent than B, so B must be shifted to the right.

!  m1 = number of initial A words to be copied without adding B.
!  m2 = end index of A words to be added to shifted B words, after end of initial A.
!  m3 = end index of A words to be copied without adding, after end of shifted B section.
!  m4 = end index of zero words after the end of A words.
!  m5 = end index of B words to be copied with a shift, after end of A words.

  m1 = min (na, ish)
  m2 = min (na, nb + ish)
  m3 = na
  m4 = min (max (na, ish), mpnw + 1)
  m5 = min (max (na, nb + ish), mpnw + 1)

  do i = 1, m1
    d(i+3) = a(i+3)
  enddo

  do i = m1 + 1, m2
    d(i+3) = a(i+3) + idb * b(i+2-ish+1)
  enddo

  do i = m2 + 1, m3
    d(i+3) = a(i+3)
  enddo

  do i = m3 + 1, m4
    d(i+3) = 0
  enddo

  do i = m4 + 1, m5
    d(i+3) = idb * b(i+2-ish+1)
  enddo

  nd = m5
  ixd = ixa
  d(nd+4) = 0
  d(nd+5) = 0
else

!   B has greater exponent than A, so A must be shifted to the right.

  nsh = - ish
  m1 = min (nb, nsh)
  m2 = min (nb, na + nsh)
  m3 = nb
  m4 = min (max (nb, nsh), mpnw + 1)
  m5 = min (max (nb, na + nsh), mpnw + 1)

  do i = 1, m1
    d(i+3) = idb * b(i+3)
  enddo

  do i = m1 + 1, m2
    d(i+3) = a(i+2-nsh+1) + idb * b(i+3)
  enddo

  do i = m2 + 1, m3
    d(i+3) = idb * b(i+3)
  enddo

  do i = m3 + 1, m4
    d(i+3) = 0
  enddo

  do i = m4 + 1, m5
    d(i+3) = a(i+2-nsh+1)
  enddo

  nd = m5
  ixd = ixb
  d(nd+4) = 0
  d(nd+5) = 0
endif

!   Call mpnorm to fix up result and store in c.

d(0) = mpnw + 6
d(1) = mpnw
d(2) = sign (nd, ia)
d(3) = ixd

call mpnorm (d, c, mpnw)

100 continue

return
end subroutine mpadd

subroutine mpcabs (a, b, mpnw)

!   This routine returns the absolute value of the MPC argument A (the
!   result is of type MPR).

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:)
integer (mpiknd), intent(out):: b(0:)
integer la, mpnw1
integer (mpiknd) s0(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6)

! End of declaration

la = a(0)
if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. a(la) < abs (a(la+2)) + 4 &
  .or. b(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCABS: uninitialized or inadequately sized arrays')
  call mpabrt ( 202)
endif

mpnw1 = mpnw + 1
s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7
call mpmul (a, a, s0, mpnw1)
call mpmul (a(la:), a(la:), s1, mpnw1)
call mpadd (s0, s1, s2, mpnw1)
call mpsqrt (s2, s0, mpnw1)
call mproun (s0, mpnw)
call mpeq (s0, b, mpnw)

return
end subroutine mpcabs

subroutine mpcadd (a, b, c, mpnw)

!   This routine adds the MPC numbers A and B.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:), b(0:)
integer (mpiknd), intent(out):: c(0:)
integer la, lb, lc

! End of declaration

la = a(0)
lb = b(0)
lc = c(0)
if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. a(la) < abs (a(la+2)) + 4 &
  .or. b(0) < abs (b(2)) + 4 .or. b(lb) < abs (b(lb+2)) + 4 .or. &
  c(0) < mpnw + 6 .or. c(lc) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCADD: uninitialized or inadequately sized arrays')
  call mpabrt ( 203)
endif

call mpadd (a, b, c, mpnw)
call mpadd (a(la:), b(lb:), c(lc:), mpnw)
return
end subroutine mpcadd

subroutine mpcdiv (a, b, c, mpnw)

!   This routine divides the MPC numbers A and B.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:), b(0:)
integer (mpiknd), intent(out):: c(0:)
integer la, lb, lc, mpnw1
integer (mpiknd) s0(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6), &
  s4(0:mpnw+6)

! End of declaration

la = a(0)
lb = b(0)
lc = c(0)
if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. a(la) < abs (a(la+2)) + 4 &
  .or. b(0) < abs (b(2)) + 4 .or. b(lb) < abs (b(lb+2)) + 4 .or. &
  c(0) < mpnw + 6 .or. c(lc) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCDIV: uninitialized or inadequately sized arrays')
  call mpabrt ( 204)
endif

mpnw1 = mpnw + 1
s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7
s3(0) = mpnw + 7
s4(0) = mpnw + 7

call mpmul (a, b, s0, mpnw1)
call mpmul (a(la:), b(lb:), s1, mpnw1)
call mpadd (s0, s1, s2, mpnw1)
call mpmul (a, b(lb:), s0, mpnw1)
s0(2) = - s0(2)
call mpmul (a(la:), b, s1, mpnw1)
call mpadd (s0, s1, s3, mpnw1)

call mpmul (b, b, s0, mpnw1)
call mpmul (b(lb:), b(lb:), s1, mpnw1)
call mpadd (s0, s1, s4, mpnw1)
call mpdiv (s2, s4, s0, mpnw1)
call mpdiv (s3, s4, s1, mpnw1)


call mproun (s0, mpnw)
call mproun (s1, mpnw)
call mpeq (s0, c, mpnw)
call mpeq (s1, c(lc:), mpnw)

return
end subroutine mpcdiv

subroutine mpceq (a, b, mpnw)

!   Sets the MPC number B equal to A.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:)
integer (mpiknd), intent(out):: b(0:)
integer i, ia, la, lb, na

! End of declaration

la = a(0)
lb = b(0)
if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. a(la) < abs (a(la+2)) + 4 &
  .or. b(0) < mpnw + 6 .or. b(lb) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCEQ: uninitialized or inadequately sized arrays')
  call mpabrt ( 205)
endif

ia = sign (int (1, mpiknd), a(2))
na = min (int (abs (a(2))), mpnw)
if (na == 0)  then
  b(1) = mpnw
  b(2) = 0
  b(3) = 0
  b(4) = 0
  b(5) = 0
  goto 110
endif
b(1) = mpnw
b(2) = sign (na, ia)

do i = 2, na + 2
  b(i+1) = a(i+1)
enddo

b(na+4) = 0
b(na+5) = 0

110 continue

ia = sign (int (1, mpiknd), a(la+2))
na = min (int (abs (a(la+2))), mpnw)
if (na == 0)  then
  b(lb+1) = mpnw
  b(lb+2) = 0
  b(lb+3) = 0
  b(lb+4) = 0
  b(lb+5) = 0
  goto 120
endif
b(lb+1) = mpnw
b(lb+2) = sign (na, ia)

do i = 2, na + 2
  b(i+lb+1) = a(i+la+1)
enddo

b(na+lb+4) = 0
b(na+lb+5) = 0

120 continue

return
end subroutine mpceq

subroutine mpcmul (a, b, c, mpnw)

!   This routine multiplies the MPC numbers A and B.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:), b(0:)
integer (mpiknd), intent(out):: c(0:)
integer la, lb, lc, mpnw1
integer (mpiknd) s0(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6)

! End of declaration

la = a(0)
lb = b(0)
lc = c(0)
if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. a(la) < abs (a(la+2)) + 4 &
  .or. b(0) < abs (b(2)) + 4 .or. b(lb) < abs (b(lb+2)) + 4 .or. &
  c(0) < mpnw + 6 .or. c(lc) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCMUL: uninitialized or inadequately sized arrays')
  call mpabrt ( 206)
endif

mpnw1 = mpnw + 1
s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7
s3(0) = mpnw + 7

call mpmul (a, b, s0, mpnw1)
call mpmul (a(la:), b(lb:), s1, mpnw1)
call mpsub (s0, s1, s2, mpnw1)
call mpmul (a, b(lb:), s0, mpnw1)
call mpmul (a(la:), b, s1, mpnw1)
call mpadd (s0, s1, s3, mpnw1)

call mproun (s2, mpnw)
call mproun (s3, mpnw)
call mpeq (s2, c, mpnw)
call mpeq (s3, c(lc:), mpnw)

return
end subroutine mpcmul

subroutine mpcnpwr (a, n, b, mpnw)

!   This computes the N-th power of the MPC number A and returns the MPC result
!   in B. When N is zero, 1 is returned. When N is negative, the reciprocal
!   of A ^ |N| is returned.

implicit none
integer, intent(in):: mpnw, n
integer (mpiknd), intent(in):: a(0:)
integer (mpiknd), intent(out):: b(0:)
real (mprknd), parameter:: cl2 = 1.4426950408889633d0
integer j, kk, kn, la, lb, lc, mn, mpnw1, na, nn
real (mprknd) t1
integer (mpiknd) s0(0:2*mpnw+13), s1(0:2*mpnw+13), s2(0:2*mpnw+13)

! End of declaration

la = a(0)
lb = b(0)
if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. a(la) < abs (a(la+2)) + 4 .or. &
  b(0) < mpnw + 6 .or. b(lb) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCNPWR: uninitialized or inadequately sized arrays')
  call mpabrt ( 207)
endif

na = min (max (int (abs (a(2))), int (abs (a(la+2)))), mpnw)
if (na == 0) then
  if (n >= 0) then
    b(1) = mpnw
    b(2) = 0
    b(3) = 0
    b(4) = 0
    b(5) = 0
    goto 120
  else
    write (mpldb, 2)
2   format ('*** MPCNPWR: Argument is zero and N is negative or zero.')
    call mpabrt ( 208)
  endif
endif

mpnw1 = mpnw + 1
lc = mpnw + 7
s0(0) = mpnw + 7
s0(lc) = mpnw + 7
s1(0) = mpnw + 7
s1(lc) = mpnw + 7
s2(0) = mpnw + 7
s2(lc) = mpnw + 7

nn = abs (n)
if (nn == 0) then
  call mpdmc (1.d0, 0, b, mpnw)
  call mpdmc (0.d0, 0, b(lb:), mpnw)
  goto 120
elseif (nn == 1) then
  call mpceq (a, s2, mpnw1)
  goto 110
elseif (nn == 2) then
  call mpcmul (a, a, s2, mpnw1)
  goto 110
endif

!   Determine the least integer MN such that 2 ^ MN > NN.

t1 = nn
mn = cl2 * log (t1) + 1.d0 + mprdfz
call mpdmc (1.d0, 0, s2, mpnw1)
call mpdmc (0.d0, 0, s2(lc:), mpnw1)
call mpceq (a, s0, mpnw1)
kn = nn

!   Compute B ^ N using the binary rule for exponentiation.

do j = 1, mn
  kk = kn / 2
  if (kn /= 2 * kk) then
    call mpcmul (s2, s0, s1, mpnw1)
    call mpceq (s1, s2, mpnw1)
  endif
  kn = kk
  if (j < mn) then
    call mpcmul (s0, s0, s1, mpnw1)
    call mpceq (s1, s0, mpnw1)
  endif
enddo

!   Compute reciprocal if N is negative.

110 continue

if (n < 0) then
  call mpdmc (1.d0, 0, s1, mpnw1)
  call mpdmc (0.d0, 0, s1(lc:), mpnw1)
  call mpcdiv (s1, s2, s0, mpnw1)
  call mpceq (s0, s2, mpnw1)
endif

!   Restore original precision level.

call mproun (s2, mpnw)
call mproun (s2(lc:), mpnw)
call mpceq (s2, b, mpnw)

120 continue
return
end subroutine mpcnpwr

subroutine mpconjg (a, b, mpnw)

!   This routine returns the conjugate of the MPC argument A.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:)
integer (mpiknd), intent(out):: b(0:)
integer la, lb

! End of declaration

la = a(0)
lb = b(0)
if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. a(la) < abs (a(la+2)) + 4 &
  .or. b(0) < mpnw + 6 .or. b(lb) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCONJ: uninitialized or inadequately sized arrays')
  call mpabrt ( 209)
endif

call mpceq (a, b, mpnw)
b(lb+2) = - b(lb+2)
return
end subroutine mpconjg

subroutine mpcsqrt (a, b, mpnw)

!   This routine returns the square root of the MPC argument A.
!   The formula is:

!   1/Sqrt[2] * (Sqrt[r + a1] + I * a2 / Sqrt[r + a1])  if a1 >= 0, or
!   1/Sqrt[2] * (|a2| / Sqrt[r - a1] + I * Sgn[a2] * Sqrt[r - a1]) if a1 < 0,

!   where r = Sqrt[a1^2 + a2^2], and a1 and a2 are the real and imaginary
!   parts of A.

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
if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. a(la) < abs (a(la+2)) + 4 &
  .or. b(0) < mpnw + 6 .or. b(lb) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCSQRT: uninitialized or inadequately sized arrays')
  call mpabrt ( 210)
endif

mpnw1 = mpnw + 1
s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7
s3(0) = mpnw + 7
s4(0) = mpnw + 7

call mpmul (a, a, s0, mpnw1)
call mpmul (a(la:), a(la:), s1, mpnw1)
call mpadd (s0, s1, s2, mpnw1)
call mpsqrt (s2, s0, mpnw1)

if (a(2) >= 0) then
  call mpadd (s0, a, s1, mpnw1)
  call mpsqrt (s1, s0, mpnw1)
  call mpdiv (a(la:), s0, s1, mpnw1)
else
  call mpsub (s0, a, s2, mpnw1)
  call mpsqrt (s2, s1, mpnw1)
  call mpdiv (a(la:), s1, s0, mpnw1)
  s0(2) = abs (s0(2))
  s1(2) = sign (s1(2), a(la+2))
endif

call mpdmc (0.5d0, 0, s3, mpnw1)
call mpsqrt (s3, s2, mpnw1)
call mpmul (s0, s2, s3, mpnw1)
call mpmul (s1, s2, s4, mpnw1)

call mproun (s3, mpnw)
call mproun (s4, mpnw)
call mpeq (s3, b, mpnw)
call mpeq (s4, b(lb:), mpnw)

return
end subroutine mpcsqrt

subroutine mpcsub (a, b, c, mpnw)

!   This routine subtracts the MPC numbers A and B.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:), b(0:)
integer (mpiknd), intent(out):: c(0:)
integer la, lb, lc

! End of declaration

la = a(0)
lb = b(0)
lc = c(0)
if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. a(la) < abs (a(la+2)) + 4 &
  .or. b(0) < abs (b(2)) + 4 .or. b(lb) < abs (b(lb+2)) + 4 .or. &
  c(0) < mpnw + 6 .or. c(lc) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCSUB: uninitialized or inadequately sized arrays')
  call mpabrt ( 211)
endif

call mpsub (a, b, c, mpnw)
call mpsub (a(la:), b(lb:), c(lc:), mpnw)
return
end subroutine mpcsub

subroutine mpcpr (a, b, ic, mpnw)

!   This routine compares the MPR numbers A and B and returns in IC the value
!   -1, 0, or 1 depending on whether A < B, A = B, or A > B.
!   Note that the first and second words do NOT need to be the same for the
!   result to be "equal".

implicit none
integer (mpiknd), intent(in):: a(0:), b(0:)
integer, intent(in):: mpnw
integer, intent(out):: ic
integer (mpiknd) s0(0:mpnw+5)

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. b(0) < abs (b(2)) + 4) then
  write (mpldb, 1)
1 format ('*** MPCPR: uninitialized or inadequately sized arrays')
  call mpabrt ( 212)
endif

s0(0) = mpnw + 6
call mpsub (a, b, s0, mpnw)
if (s0(2) < 0) then
  ic = -1
elseif (s0(2) == 0) then
  ic = 0
else
  ic = 1
endif

return
end subroutine mpcpr

subroutine mpdiv (a, b, c, mpnw)

!   This divides A by B and returns the result in C.

!   This subroutine employs the following Newton-Raphson iteration, which
!   converges to 1 / B:

!    X_{k+1} = X_k + (1 - X_k * B) * X_k

!   where the multiplication () * X_k is performed with only half of the
!   normal level of precision. These iterations are performed with a
!   working precision level MPNW that is dynamically changed, approximately
!   doubling with each iteration (except that at iteration NIT before the
!   final iteration, the iteration is repeated without doubling the
!   precision, in order to enhance accuracy). The final iteration is
!   performed as follows (this is due to A. Karp):

!    A / B = (A * X_n) + [A - (A * X_n) * B] * X_n  (approx.)

!   where the multiplications A * X_n and [] * X_n are performed with only
!   half of the final level of precision.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:), b(0:)
integer (mpiknd), intent(out):: c(0:)
integer, parameter:: nit = 3
real (mprknd), parameter:: cl2 = 1.4426950408889633d0
integer iq, k, mpnw1, mq, n, na, nb, nw1, nw2
real (mprknd) t1, t2
integer (mpiknd) s0(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6)

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. b(0) < abs (b(2)) + 4 .or. &
  c(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPDIV: uninitialized or inadequately sized arrays')
  call mpabrt ( 213)
endif

na = min (int (abs (a(2))), mpnw)
nb = min (int (abs (b(2))), mpnw)

if (na == 0) then
  c(1) = mpnw
  c(2) = 0
  c(3) = 0
  c(4) = 0
  c(5) = 0
  goto 120
endif
if (nb == 0.d0) then
  write (mpldb, 2)
2 format ('*** MPDIV: Divisor is zero.')
  call mpabrt ( 214)
  return
endif

s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7
s3(0) = mpnw + 7

!   Determine the least integer MQ such that 2 ^ MQ >= MPNW.

t1 = mpnw
mq = cl2 * log (t1) + 1.d0 - mprdfz

!   Compute the initial approximation of 1 / B.

call mpmdc (b, t1, n, mpnw)
t2 = 1.d0 / t1
call mpdmc (t2, -n, s2, mpnw)
call mpdmc (1.d0, 0, s3, mpnw)
mpnw1 = 5
iq = 0
nw1 = mpnw1
nw2 = mpnw1

!   Perform the Newton-Raphson iteration described above with a dynamically
!   changing precision level MPNW (one greater than powers of two).

do k = 1, mq - 1
  if (k > 2) then
    nw1 = mpnw1
    mpnw1 = min (2 * mpnw1 - 2, mpnw) + 1
    nw2 = mpnw1
  endif

100  continue

  call mpmul (b, s2, s1, nw2)
  call mpsub (s3, s1, s0, nw2)
  call mpmul (s2, s0, s1, nw1)
  call mpadd (s2, s1, s0, nw2)
  call mpeq (s0, s2, nw2)

  if (k == mq - nit .and. iq == 0) then
    iq = 1
    goto 100
  endif
enddo

!   Perform last iteration using Karp's trick.

nw1 = mpnw1
mpnw1 = min (2 * mpnw1 - 1, mpnw) + 1
nw2 = mpnw1

call mpmul (a, s2, s0, nw1)
call mpmul (s0, b, s1, nw2)
call mpsub (a, s1, s3, nw2)
call mpmul (s3, s2, s1, nw1)
call mpadd (s0, s1, s2, nw2)

!   Restore original precision level.

call mproun (s2, mpnw)
call mpeq (s2, c, mpnw)

120 continue
return
end subroutine mpdiv

subroutine mpdecmdr (ra, db, ib, mpnw)
implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent (in):: ra(0:)
real (mprknd), intent (out):: db
integer, intent (out):: ib
real (mprknd), parameter:: alg102 = 0.301029995663981195d0
real (mprknd) dt1, dt2
integer i1

call mpmdc (ra, dt1, i1, mpnw)

if (dt1 /= 0.d0) then
  dt2 = alg102 * i1 + log10 (abs (dt1))
  ib = dt2
  if (dt2 < 0.d0) ib = ib - 1
  db = sign (10.d0 ** (dt2 - ib), dt1)
else
  db = 0.d0
  ib = 0
endif

return
end subroutine mpdecmdr

subroutine mpdivd (a, b, c, mpnw)

!   This routine divides the MPR number A by the DP number B to yield C.

!   NOTE however that the product is not fully accurate unless B is an exact
!   binary value.
!   Examples of exact binary values (good): 123456789.d0, 0.25d0, -5.3125d0.
!   Examples of inexact binary values (bad): 0.1d0, 123467.8d0, -3333.3d0.

implicit none
integer, intent(in):: mpnw
real (mprknd), intent(in):: b
integer (mpiknd), intent(in):: a(0:)
integer (mpiknd), intent(out):: c(0:)
integer, parameter:: nbth = mpnbt / 2
real (mprknd), parameter:: bdh = 2.d0**nbth, rdh = 0.5d0**nbth
real (mprknd) bb, bdvd
integer i, ia, ib, j, k, na, n1
integer (mpiknd) cc(0:2*mpnw+10), d(0:2*mpnw+10), ibb, &
  b1, b2, c11, c12, c21, c22, d1, d2, td, t1, t2, t3

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. c(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPDIVD: uninitialized or inadequately sized arrays')
  call mpabrt ( 215)
endif

!   Check for zero inputs.

ia = sign (int (1, mpiknd), a(2))
na = min (int (abs (a(2))), mpnw)
ib = sign (1.d0, b)
if (na == 0) then
  c(1) = mpnw
  c(2) = 0
  c(3) = 0
  c(4) = 0
  c(5) = 0
  goto 140
elseif (b == 0.d0) then
  write (mpldb, 2)
2 format ('*** MPDIVD: Divisor is zero.')
  call mpabrt ( 216)
elseif (b == 1.d0) then
  call mpeq (a, c, mpnw)
  goto 140
endif
bb = abs (b)
n1 = 0

!   Reduce BB to within 1 and MPBDX.

if (bb >= mpbdx) then
  do k = 1, 100
    bb = bb / mpbdx
    if (bb < mpbdx) then
      n1 = n1 + k
      goto 120
    endif
  enddo
elseif (bb < 1.d0) then
  do k = 1, 100
    bb = mpbdx * bb
    if (bb >= 1.d0) then
      n1 = n1 - k
      goto 120
    endif
  enddo
endif

120  continue

ibb = bb

!   If bb is not an integer, call mpdiv instead.

if (bb /= ibb) then
  d(0) = mpnw + 6
  d(1) = mpnw
  call mpdmc (b, 0, d, mpnw)
  call mpdiv (a, d, c, mpnw)
  goto 140
endif

cc(0) = 0
cc(1) = mpnw
cc(2) = sign (mpnw, ia * ib)
cc(3:2*mpnw+10) = 0

!   Split D array into half-word chunks.

d(0:3) = 0

do i = 0, na
  c11 = shifta (a(i+4), nbth)
  c12 = a(i+4) - shiftl (c11, nbth)
  d(2*i+4) = c11
  d(2*i+5) = c12
enddo

d(2*na+6:2*mpnw+10) = 0
b1 = shifta (ibb, nbth)
b2 = ibb - shiftl (b1, nbth)

!   Perform short division algorithm, after splitting inputs.
!   Double precision is employed to find and refine the trial divisor.

do j = 3, 2 * mpnw + 5
  bdvd = bdh * d(j) + d(j+1) + d(j+2) * rdh
  td = floor (bdvd / bb, mpiknd)
  t1 = b1 * td
  c11 = shifta (t1, nbth)
  c12 = t1 - shiftl (c11, nbth)
  t2 = b2 * td
  c21 = shifta (t2, nbth)
  c22 = t2 - shiftl (c21, nbth)
  d1 = c12 + c21 + shiftl (c11, nbth)
  d2 = c22
  d(j) = d(j) - d1
  d(j+1) = d(j+1) - d2 + shiftl (d(j), nbth)
  cc(j+1) = td
enddo

!  Release carries on the full cc vector.

t1 = 0

do i = 2 * mpnw + 5, 3, -1
  t3 = t1 + cc(i+1)
  t1 = shifta (t3, nbth)
  cc(i+1) = t3 - shiftl (t1, nbth)
enddo

cc(3) = cc(3) + t1

!  Recombine half words into full words.

c(1) = mpnw
c(2) = sign (mpnw, ia * ib)
c(3) = cc(3)
c(4) = cc(4)

do i = 0, mpnw + 1
  c(i+4) = shiftl (cc(2*i+4), nbth) + cc(2*i+5)
enddo

!   If c(3) is nonzero, shift the result one cell right.

if (c(3) /= 0) then
  n1 = n1 + 1
  c(2) = sign (abs (c(2)) + 1, c(2))

  do i = mpnw + 4, 3, -1
    c(i+1) = c(i)
  enddo
endif

c(3) = a(3) - n1
call mproun (c, mpnw)

140 continue

return
end subroutine mpdivd

subroutine mpdivd40 (a, b, c, mpnw)

!   This routine divides the MPR number A by the DP number B to yield C.
!   In contrast to mpdivd, this routine only allows 40 significant bits
!   (approximately 12 significant decimal digits) in B. If more nonzero bits
!   are present in B (likely due to inexact binary value), an error is flagged.

!   Examples of exact binary values (good): 123456789.d0, 0.25d0, -5.3125d0.
!   Examples of inexact binary values (bad): 0.1d0, 123467.8d0, -3333.3d0.


implicit none
integer, intent(in):: mpnw
real (mprknd), intent(in):: b
integer (mpiknd), intent(in):: a(0:)
integer (mpiknd), intent(out):: c(0:)
real (mprknd) t2

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. c(0) < mpnw + 6) then
 write (mpldb, 1)
1 format ('*** MPDIVD40: uninitialized or inadequately sized arrays')
  call mpabrt ( 217)
endif

!   Check whether B has more than 40 significant bits (actually whether the
!   trailing 13 bits are zero).

t2 = mpmask13 (b)
if (t2 == abs (b)) then
  call mpdivd (a, b, c, mpnw)
else
  write (mpldb, 2) b
2 format ('*** MPDIVD40: DP value has more than 40 significant bits:', &
  1p,d25.15/'and thus very likely represents an unintended loss of accuracy.'/ &
  'Fix the issue, or else use functions mpprod, mpquot, mpreald or mpcmplxdc.'/ &
  'See documentation for details.')
  call mpabrt ( 218)
endif

return
end subroutine mpdivd40

subroutine mpdmc (a, n, b, mpnw)

!   This routine converts the DP number A * 2^N to MPR form in B.

!   NOTE however that the product is not fully accurate unless A is an exact
!   binary value.
!   Examples of exact binary values (good): 123456789.d0, 0.25d0, -5.3125d0.
!   Examples of inexact binary values (bad): 0.1d0, 123467.8d0, -3333.3d0.

implicit none
integer, intent(in):: mpnw, n
real (mprknd), intent(in):: a
integer (mpiknd), intent(out):: b(0:)
integer i, k, n1, n2
real (mprknd) aa

! End of declaration

if (mpnw < 4 .or. b(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPDMC: uninitialized or inadequately sized arrays')
  call mpabrt ( 219)
endif

!   Check for zero.

if (a == 0.d0) then
  b(1) = mpnw
  b(2) = 0
  b(3) = 0
  b(4) = 0
  b(5) = 0
  goto 150
endif
n1 = n / mpnbt
n2 = n - mpnbt * n1
aa = abs (a) * 2.d0 ** n2

!   Reduce AA to within 1 and MPBDX.

if (aa >= mpbdx) then

  do k = 1, 100
    aa = aa / mpbdx
    if (aa < mpbdx) then
      n1 = n1 + k
      goto 120
    endif
 enddo

elseif (aa < 1.d0) then

  do k = 1, 100
    aa = aa * mpbdx
    if (aa >= 1.d0) then
      n1 = n1 - k
      goto 120
    endif
  enddo

endif

!   Store successive sections of AA into B.

120  continue

b(3) = n1
b(4) = int (aa, mpiknd)
aa = mpbdx * (aa - b(4))
b(5) = int (aa, mpiknd)
b(6) = 0
b(7) = 0
b(8) = 0

do i = 6, 3, -1
  if (b(i+1) /= 0) goto 140
enddo

140  continue

b(1) = mpnw
aa = i - 2
b(2) = sign (aa, a)

150 continue

return
end subroutine mpdmc

subroutine mpdmc40 (a, n, b, mpnw)

!   This routine converts the DP number A * 2^N to MPR form in B.
!   In contrast to mpdmc, this routine only allows 40 significant bits
!   (approximately 12 significant decimal digits) in A. If more nonzero bits
!   are present in A (likely due to inexact binary value), an error is flagged.

!   Examples of exact binary values (good): 123456789.d0, 0.25d0, -5.3125d0.
!   Examples of inexact binary values (bad): 0.1d0, 123467.8d0, -3333.3d0.


implicit none
integer, intent(in):: mpnw, n
real (mprknd), intent(in):: a
integer (mpiknd), intent(out):: b(0:)
real (mprknd) t2

! End of declaration

if (mpnw < 4 .or. b(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPDMC40: uninitialized or inadequately sized arrays')
  call mpabrt ( 220)
endif

!   Check whether A has more than 40 significant bits (actually whether
!   the trailing 13 bits are zero).

t2 = mpmask13 (a)
if (t2 == abs (a)) then
  call mpdmc (a, n, b, mpnw)
else
  write (mpldb, 2) a
2 format ('*** MPDMC40: DP value has more than 40 significant bits:', &
  1p,d25.15/'and thus very likely represents an unintended loss of accuracy.'/ &
  'Fix the issue, or else use functions mpprod, mpquot, mpreald or mprealdm.'/ &
  'See documentation for details.')
  call mpabrt ( 221)
endif

return
end subroutine mpdmc40

subroutine mpeq (a, b, mpnw)

!   Sets the MPR number B equal to the MPR number A.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:)
integer (mpiknd), intent(out):: b(0:)
integer i, ia, na

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. b(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPEQ: uninitialized or inadequately sized arrays')
  call mpabrt ( 222)
endif

ia = sign (int (1, mpiknd), a(2))
na = min (int (abs (a(2))), mpnw)

if (na == 0)  then
  b(1) = mpnw
  b(2) = 0
  b(3) = 0
  b(4) = 0
  b(5) = 0
  goto 110
endif
b(1) = mpnw
b(2) = sign (na, ia)

do i = 2, na + 2
  b(i+1) = a(i+1)
enddo

b(na+4) = 0
b(na+5) = 0

110 continue

return
end subroutine mpeq

subroutine mpinfr (a, b, c, mpnw)

!   Sets B to the integer part of the MPR number A and sets C equal to the
!   fractional part of A. Note this is NOT the quite same as the greatest
!   integer function as often defined in some mathematical books and papers.
!   Examples:  If A = 1.75, then B = 1., C = 0.75.
!     If A = -3.25, then B = -3., C = -0.25.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:)
integer (mpiknd), intent(out):: b(0:), c(0:)
integer i, ia, ma, na, nb, nc

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. b(0) < mpnw + 6 .or. &
  c(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPINFR: uninitialized or inadequately sized arrays')
  call mpabrt ( 223)
endif

!   Check if  A  is zero.

ia = sign (int (1, mpiknd), a(2))
na = min (int (abs (a(2))), mpnw)
ma = a(3)
if (na == 0)  then
  b(1) = mpnw
  b(2) = 0
  b(3) = 0
  b(4) = 0
  b(5) = 0
  c(1) = mpnw
  c(2) = 0
  c(3) = 0
  c(4) = 0
  c(5) = 0
  goto 120
endif

if (ma >= mpnw - 1) then
  write (mpldb, 2)
2 format ('*** MPINFR: Argument is too large.')
  call mpabrt ( 224)
endif

!   Place integer part in  B.

nb = min (max (ma + 1, 0), na)
if (nb == 0) then
  b(1) = mpnw
  b(2) = 0
  b(3) = 0
  b(4) = 0
  b(5) = 0
else
  b(1) = mpnw
  b(2) = sign (nb, ia)
  b(3) = ma
  b(nb+4) = 0
  b(nb+5) = 0

  do i = 3, nb + 2
    b(i+1) = a(i+1)
  enddo
endif

!   Place fractional part in C.

nc = na - nb
if (nc <= 0) then
  c(1) = mpnw
  c(2) = 0
  c(3) = 0
  c(4) = 0
  c(5) = 0
else
  c(1) = mpnw
  c(2) = sign (nc, ia)
  c(3) = ma - nb
  c(nc+4) = 0
  c(nc+5) = 0

  do i = 3, nc + 2
    c(i+1) = a(i+nb+1)
  enddo
endif

!   Fix up results. B may have trailing zeros and C may have leading zeros.

call mproun (b, mpnw)
call mproun (c, mpnw)

120  continue
return
end subroutine mpinfr

subroutine mpinitwds (ra, mpnw)

!   This initializes ra with mpnw+6 space and mpnw working precision.

implicit none
integer (mpiknd), intent(out):: ra(0:)
integer, intent(in):: mpnw

ra(0) = mpnw + 6
ra(1) = mpnw
ra(2) = 0
ra(3) = 0
ra(4) = 0
return
end subroutine mpinitwds

subroutine mpmdc (a, b, n, mpnw)

!   This returns a DP approximation the MPR number A in the form B * 2^n.

implicit none
integer (mpiknd), intent(in):: a(0:)
integer, intent(in):: mpnw
integer, intent(out):: n
real (mprknd), intent(out):: b
integer na
real (mprknd) aa

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4) then
  write (mpldb, 1)
1 format ('*** MPMDC: uninitialized or inadequately sized arrays')
  call mpabrt ( 225)
endif

if (a(2) == 0.d0)  then
  b = 0.d0
  n = 0
  goto 100
endif

na = abs (a(2))
aa = a(4)
if (na >= 2) aa = aa + a(5) / dble (mpbdx)
n = mpnbt * a(3)
b = sign (aa, dble (a(2)))

!   Reduce b to within 1 and 2.

na = log (abs (b)) / log (2.d0) + mprdfz
b = b / 2.d0**na
n = n + na
if (abs (b) < 1.d0) then
  b = 2.d0 * b
  n = n - 1
elseif (abs (b) > 2.d0) then
  b = 0.5d0 * b
  n = n + 1
endif

100  continue
return
end subroutine mpmdc

subroutine mpmul (a, b, c, mpnw)

!   This routine multiplies MPR numbers A and B to yield C.

!   This routine returns up to MPNW mantissa words of the product. If the
!   complete double-long product of A and B is desired (for example in large
!   integer applications), then MPNW must be at least as large as the sum of
!   the mantissa lengths of A and B. In other words, if the precision levels
!   of A and B are both 64 words, then MPNW must be at least 128 words to
!   produce the complete double-long product in C.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:), b(0:)
integer (mpiknd), intent(out):: c(0:)
integer, parameter:: nbth = mpnbt / 2
integer i, ia, ib, j, j3, n2, na, nb, nc
integer (mpiknd) d(0:mpnw+6), a1, a2, b1, b2, c1, c2, c3, dd, t1, t3

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. b(0) < abs (b(2)) + 4 .or. &
  c(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPMUL: uninitialized or inadequately sized arrays')
  call mpabrt ( 226)
endif

ia = sign (int (1, mpiknd), a(2))
ib = sign (int (1, mpiknd), b(2))
na = min (int (abs (a(2))), mpnw)
nb = min (int (abs (b(2))), mpnw)
nc = min (na + nb, mpnw)

if (na == 0 .or. nb == 0) then

!   One of the inputs is zero -- result is zero.

  c(1) = mpnw
  c(2) = 0
  c(3) = 0
  c(4) = 0
  c(5) = 0
  goto 200
endif

if (na == 1 .and. a(4) == 1) then

!   A is 1 or -1 -- result is B or -B.

  c(1) = mpnw
  c(2) = sign (nb, ia * ib)
  c(3) = a(3) + b(3)

  do i = 3, nb + 2
    c(i+1) = b(i+1)
  enddo

  c(nb+4) = 0
  c(nb+5) = 0
  goto 200
elseif (nb == 1 .and. b(4) == 1) then

!   B is 1 or -1 -- result is A or -A.

  c(1) = mpnw
  c(2) = sign (na, ia * ib)
  c(3) = a(3) + b(3)

  do i = 3, na + 2
    c(i+1) = a(i+1)
  enddo

  c(na+4) = 0
  c(na+5) = 0
  goto 200
endif

!   For very high precision, call mpmulx.

if (na > mpmlxm .and. nb > mpmlxm) then
  call mpmulx (a, b, c, mpnw)
  goto 200
endif

dd = a(3) + b(3)
d(0) = mpnw + 6
d(1) = mpnw
d(2) = sign (nc, ia * ib)

do i = 2, nc + 5
  d(i+1) = 0
enddo

!   Perform ordinary long multiplication algorithm, after splitting inputs.
!   Accumulate at most MPNW+2 mantissa words of the product.

do j = 3, na + 2
  j3 = j - 3
  n2 = min (nb + 2, mpnw + 4 - j3)
  a1 = shifta (a(j+1), nbth)
  a2 = a(j+1) - shiftl (a1, nbth)

  do i = 3, n2
    b1 = shifta (b(i+1), nbth)
    b2 = b(i+1) - shiftl (b1, nbth)
    c1 = a1 * b2 + a2 * b1
    c2 = shifta (c1, nbth)
    c3 = c1 - shiftl (c2, nbth)
    d(i+j3) = d(i+j3) + a1 * b1 + c2
    d(i+j3+1) = d(i+j3+1) + a2 * b2 + shiftl (c3, nbth)
  enddo

!  Release carries on the just-computed section of the d vector.

  t1 = 0

  do i = n2, 3, -1
    t3 = t1 + d(i+j3+1)
    t1 = shifta (t3, mpnbt)
    d(i+j3+1) = t3 - shiftl (t1, mpnbt)
  enddo

  d(j3+3) = d(j3+3) + t1
enddo

!  Release carries on the full d vector.

t1 = 0

do i = nc + 1, 1, -1
  t3 = t1 + d(i+3)
  t1 = shifta (t3, mpnbt)
  d(i+3) = t3 - shiftl (t1, mpnbt)
enddo

d(3) = d(3) + t1

!   If d(3) is nonzero, shift the result one cell right.

if (d(3) /= 0) then
  dd = dd + 1
  nc = min (nc + 1, mpnw)
  d(2) = sign (nc, int (d(2)))

  do i = nc + 4, 3, -1
    d(i+1) = d(i)
  enddo
endif

d(3) = dd

do i = 1, nc + 5
  c(i) = d(i)
enddo

call mproun (c, mpnw)

200 continue

return
end subroutine mpmul

subroutine mpmuld (a, b, c, mpnw)

!   This routine multiplies the MPR number A by the DP number B to yield C.

!   NOTE however that the product is not fully accurate unless B is an exact
!   binary value.
!   Examples of exact binary values (good): 123456789.d0, 0.25d0, -5.3125d0.
!   Examples of inexact binary values (bad): 0.1d0, 123467.8d0, -3333.3d0.

implicit none
integer, intent(in):: mpnw
real (mprknd), intent(in):: b
integer (mpiknd), intent(in):: a(0:)
integer (mpiknd), intent(out):: c(0:)
integer, parameter:: nbth = mpnbt / 2
integer i, ia, ib, j, k, na, n1
real (mprknd) bb
integer (mpiknd) d(0:mpnw+6), ibb, a1, a2, b1, b2, c1, c2, c3, t1, t3

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. c(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPMULD: uninitialized or inadequately sized arrays')
  call mpabrt ( 227)
endif

!   Check for zero inputs.

ia = sign (int (1, mpiknd), a(2))
na = min (int (abs (a(2))), mpnw)
ib = sign (1.d0, b)
if (na == 0 .or. b == 0.d0) then
  c(1) = mpnw
  c(2) = 0
  c(3) = 0
  c(4) = 0
  c(5) = 0
  goto 140
elseif (b == 1.d0) then
  call mpeq (a, c, mpnw)
  goto 140
endif
bb = abs (b)
n1 = 0

!   Reduce BB to within 1 and MPBDX.

if (bb >= mpbdx) then
  do k = 1, 100
    bb = bb / mpbdx
    if (bb < mpbdx) then
      n1 = n1 + k
      goto 120
    endif
  enddo
elseif (bb < 1.d0) then
  do k = 1, 100
    bb = mpbdx * bb
    if (bb >= 1.d0) then
      n1 = n1 - k
      goto 120
    endif
  enddo
endif

120  continue

ibb = bb

!   If bb is not an integer, call mpmul instead.

if (bb /= ibb) then
  d(0) = mpnw + 6
  d(1) = mpnw
  call mpdmc (b, 0, d, mpnw)
  call mpmul (a, d, c, mpnw)
  goto 140
endif

d(0) = mpnw + 6
d(1) = mpnw
d(2) = sign (na, ia * ib)

do i = 2, na + 5
  d(i+1) = 0
enddo

b1 = shifta (ibb, nbth)
b2 = ibb - shiftl (b1, nbth)

!   Perform short multiplication algorithm, after splitting inputs.

do j = 3, na + 3
  a1 = shifta (a(j+1), nbth)
  a2 = a(j+1) - shiftl (a1, nbth)
  c1 = a1 * b2 + a2 * b1
  c2 = shifta (c1, nbth)
  c3 = c1 - shiftl (c2, nbth)
  d(j) = d(j) + a1 * b1 + c2
  d(j+1) = d(j+1) + a2 * b2 + shiftl (c3, nbth)
enddo

!  Release carries on the full d vector.

t1 = 0

do i = na + 3, 3, -1
  t3 = t1 + d(i+1)
  t1 = shifta (t3, mpnbt)
  d(i+1) = t3 - shiftl (t1, mpnbt)
enddo

d(3) = d(3) + t1

!   If d(3) is nonzero, shift the result one cell right.

if (d(3) /= 0) then
  n1 = n1 + 1
  d(2) = sign (abs (d(2)) + 1, d(2))

  do i = na + 4, 3, -1
    d(i+1) = d(i)
  enddo
endif

d(3) = a(3) + n1

!   Copy d to c and round.

do i = 1, na + 5
  c(i) = d(i)
enddo

call mproun (c, mpnw)

140 continue

return
end subroutine mpmuld

subroutine mpmuld40 (a, b, c, mpnw)

!   This routine multiples the MP number A by the DP number B to yield C.
!   In contrast to mpmuld, this routine only allows 40 significant bits
!   (approximately 12 significant decimal digits) in B. If more nonzero bits
!   are present in B (likely due to inexact binary value), an error is flagged.

!   Examples of exact binary values (good): 123456789.d0, 0.25d0, -5.3125d0.
!   Examples of inexact binary values (bad): 0.1d0, 123467.8d0, -3333.3d0.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:)
real (mprknd), intent(in):: b
integer (mpiknd), intent(out):: c(0:)
real (mprknd) t2

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. c(0) < mpnw + 6) then
 write (mpldb, 1)
1 format ('*** MPMULD40: uninitialized or inadequately sized arrays')
  call mpabrt ( 228)
endif

!   Check whether B has more than 40 significant bits (actually whether
!   the trailing 13 bits are zero).

t2 = mpmask13 (b)
if (t2 == abs (b)) then
  call mpmuld (a, b, c, mpnw)
else
  write (mpldb, 2) b
2 format ('*** MPMULD40: DP value has more than 40 significant bits:', &
  1p,d25.15/'and thus very likely represents an unintended loss of accuracy.'/ &
  'Fix the issue, or else use functions mpprod, mpquot, mpreald or mpcmplxdc.'/ &
  'See documentation for details.')
  call mpabrt ( 229)
endif

return
end subroutine mpmuld40

subroutine mpneg (ra, rb, mpnw)

!   This routine sets RB = negation of RA.

implicit none
integer (mpiknd), intent(in):: ra(0:)
integer (mpiknd), intent(out):: rb(0:)
integer mpnw, na

! End of declaration

call mpeq (ra, rb, mpnw)
na = min (abs (int (ra(2))), mpnw)
rb(2) = - sign (na, int (ra(2)))
return
end subroutine mpneg

subroutine mpnint (a, b, mpnw)

!   This sets B to the nearest integer to the MPR number A.
!   Examples:  If A = 1.49, B = 1.; if A = 3.5, B = 4; if A = -2.5, B = -3.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:)
integer (mpiknd), intent(out):: b(0:)
integer ia, ma, na
integer (mpiknd) s0(0:mpnw+5), s1(0:mpnw+5)

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. b(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPNINT: uninitialized or inadequately sized arrays')
  call mpabrt ( 230)
endif

ia = sign (int (1, mpiknd), a(2))
na = min (int (abs (a(2))), mpnw)
ma = a(3)
if (na == 0)  then

!   A is zero -- result is zero.

  b(1) = mpnw
  b(2) = 0
  b(3) = 0
  b(4) = 0
  b(5) = 0
  goto 110
endif

if (ma >= mpnw) then

!   A cannot be represented exactly as an integer.

  write (mpldb, 2)
2 format ('*** MPNINT: Argument is too large.')
  call mpabrt ( 231)
endif

!   Add or subtract 1/2 from the input, depending on its sign, then
!   return the greatest integer.

s0(0) = mpnw + 6
s1(0) = mpnw + 6

call mpdmc (0.5d0, 0, s0, mpnw)
if (ia == 1) then
  call mpadd (a, s0, s1, mpnw)
else
  call mpsub (a, s0, s1, mpnw)
endif
call mpinfr (s1, b, s0, mpnw)

110 continue
return
end subroutine mpnint

subroutine mpnorm (d, a, mpnw)

!   This converts the MP number in array D to the standard normalized form
!   in A.

!   MPNORM assumes that two extra mantissa words are input at the end of D.
!   This reduces precision loss when it is necessary to shift the result to
!   the left. All words up to index A(2) + 5 in A *must* have data, even if 0.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(out):: a(0:)
integer (mpiknd), intent(inout):: d(0:)
integer i, ia, na, n4
integer (mpiknd) a2, t1, t3

! End of declaration

if (mpnw < 4 .or. d(0) < abs (d(2)) + 4 .or. a(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPNORM: uninitialized or inadequately sized arrays')
  call mpabrt ( 232)
endif

ia = sign (int (1, mpiknd), d(2))
na = min (int (abs (d(2))), mpnw)
if (na == 0)  then
  a(1) = mpnw
  a(2) = 0
  a(3) = 0
  a(4) = 0
  a(5) = 0
  goto 170
endif
n4 = na + 4
a2 = d(3)
d(3) = 0

110 continue

t1 = 0

do i = n4, 3, -1
  t3 = t1 + d(i+1)
  t1 = shifta (t3, mpnbt)
  d(i+1) = t3 - shiftl (t1, mpnbt)
enddo

d(3) = d(3) + t1

if (d(3) < 0) then

!   D(3) is negative -- negate all words and re-normalize.

  ia = - ia
  d(4) = d(4) + mpbdx * d(3)
  d(3) = 0

  do i = 2, n4
    d(i+1) = - d(i+1)
  enddo

  goto 110
elseif (d(3) > 0) then

!   The fixup loops above "spilled" a nonzero number into D(3). Shift the
!   entire number right one cell. The exponent and length of the result
!   are increased by one.

  do i = n4, 3, -1
    a(i+1) = d(i)
  enddo

  na = min (na + 1, mpnw)
  a2 = a2 + 1
else
  do i = 3, n4
    a(i+1) = d(i+1)
  enddo
endif

!   Perform rounding and truncation.

a(1) = mpnw
a(2) = sign (na, ia)
a(3) = a2

call mproun (a, mpnw)

170 continue

return
end subroutine mpnorm

subroutine mpnpwr (a, n, b, mpnw)

!   This computes the N-th power of the MPR number A and returns the result
!   in B. When N is zero, 1 is returned. When N is negative, the reciprocal
!   of A ^ |N| is returned.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:)
integer (mpiknd), intent(out):: b(0:)
real (mprknd), parameter:: cl2 = 1.4426950408889633d0
integer j, kk, kn, mn, mpnw1, n, na, nn
real (mprknd) t1
integer (mpiknd) s0(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6)

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. b(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPNPWR: uninitialized or inadequately sized arrays')
  call mpabrt ( 233)
endif

na = min (int (abs (a(2))), mpnw)

if (na == 0) then
  if (n >= 0) then
    b(1) = mpnw
    b(2) = 0
    b(3) = 0
    b(4) = 0
    b(5) = 0
    goto 120
  else
    write (mpldb, 2)
2   format ('*** MPNPWR: Argument is zero and N is negative or zero.')
    call mpabrt ( 234)
  endif
endif

mpnw1 = mpnw + 1
s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7

nn = abs (n)
if (nn == 0) then
  call mpdmc (1.d0, 0, b, mpnw)
  goto 120
elseif (nn == 1) then
  call mpeq (a, s2, mpnw1)
  goto 110
elseif (nn == 2) then
  call mpmul (a, a, s2, mpnw1)
  goto 110
endif

!   Determine the least integer MN such that 2 ^ MN > NN.

t1 = nn
mn = cl2 * log (t1) + 1.d0 + mprdfz
call mpdmc (1.d0, 0, s2, mpnw1)
call mpeq (a, s0, mpnw1)
kn = nn

!   Compute B ^ N using the binary rule for exponentiation.

do j = 1, mn
  kk = kn / 2
  if (kn /= 2 * kk) then
    call mpmul (s2, s0, s1, mpnw1)
    call mpeq (s1, s2, mpnw1)
  endif
  kn = kk
  if (j < mn) then
    call mpmul (s0, s0, s1, mpnw1)
    call mpeq (s1, s0, mpnw1)
  endif
enddo

!   Compute reciprocal if N is negative.

110 continue

if (n < 0) then
  call mpdmc (1.d0, 0, s1, mpnw1)
  call mpdiv (s1, s2, s0, mpnw1)
  call mpeq (s0, s2, mpnw1)
endif

!   Restore original precision level.

call mproun (s2, mpnw)
call mpeq (s2, b, mpnw)

120 continue

return
end subroutine mpnpwr

subroutine mpnrtr (a, n, b, mpnw)

!   This computes the N-th root of the MPR number A and returns result inB.
!   N must be at least one and must not exceed 2 ^ 30.

!   This subroutine employs the following Newton-Raphson iteration, which
!   converges to A ^ (-1/N):

!    X_{k+1} = X_k + (X_k / N) * (1 - A * X_k^N)

!   The reciprocal of the final approximation to A ^ (-1/N) is the N-th root.
!   These iterations are performed with a maximum precision level MPNW that
!   is dynamically changed, approximately doubling with each iteration.

!   When N is large and A is very near one, the following binomial series is
!   employed instead of the Newton scheme:

!   (1 + x)^(1/N)  =  1  +  x / N  +  x^2 * (1 - N) / (2! N^2)  +  ...

!   See the comment about the parameter NIT in MPDIVX.

implicit none
integer, intent(in):: mpnw, n
integer (mpiknd), intent(in):: a(0:)
integer (mpiknd), intent(out):: b(0:)
integer, parameter:: nit = 3, n30 = 2 ** 30
real (mprknd), parameter:: alt = 0.693147180559945309d0, cl2 = 1.4426950408889633d0
integer ia, iq, k, mpnw1, mq, na, n1, n2, n3
real (mprknd) t1, t2, tn
integer (mpiknd) f1(0:8), s0(0:mpnw+7), s1(0:mpnw+7), s2(0:mpnw+7), s3(0:mpnw+7)

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. b(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPNRTR: uninitialized or inadequately sized arrays')
  call mpabrt ( 235)
endif

ia = sign (int (1, mpiknd), a(2))
na = min (int (abs (a(2))), mpnw)

if (na == 0) then
  b(1) = mpnw
  b(2) = 0
  b(3) = 0
  b(4) = 0
  b(5) = 0
  goto 140
endif
if (ia < 0) then
  write (mpldb, 2)
2 format ('*** MPNRTR: Argument is negative.')
  call mpabrt ( 236)
endif

if (n <= 0 .or. n > n30) then
  write (mpldb, 3) n
3 format ('*** MPNRTR: Improper value of N',i10)
  call mpabrt ( 237)
endif

!   If N = 1 or 2, call MPEQ or MPSQRT instead.

if (n == 1) then
  call mpeq (a, b, mpnw)
  goto 140
elseif (n == 2) then
  call mpsqrt (a, b, mpnw)
  goto 140
endif

s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7
s3(0) = mpnw + 7
mpnw1 = mpnw + 1

!   Set f1 = 1.

f1(0) = 9
f1(1) = mpnw1
f1(2) = 1
f1(3) = 0
f1(4) = 1
f1(5) = 0
f1(6) = 0

!   Determine the least integer MQ such that 2 ^ MQ >= MPNW.

t1 = mpnw
mq = cl2 * log (t1) + 1.d0 - mprdfz

!   Check how close A is to 1.

call mpsub (a, f1, s0, mpnw1)
if (s0(2) == 0) then
  call mpeq (f1, b, mpnw)
  goto 140
endif
call mpmdc (s0, t1, n1, mpnw1)
n2 = cl2 * log (abs (t1))
t1 = t1 * 0.5d0 ** n2
n1 = n1 + n2

if (n1 <= -30) then
  t2 = n
  n2 = cl2 * log (t2) + 1.d0 + mprdfz
  n3 = - mpnbt * mpnw1 / n1
  if (n3 < 1.25d0 * n2) then

!   A is so close to 1 that it is cheaper to use the binomial series.

    call mpdivd (s0, t2, s1, mpnw1)
    call mpadd (f1, s1, s2, mpnw1)
    k = 0

100 continue

  k = k + 1
    t1 = 1 - k * n
    t2 = (k + 1) * n
    call mpmuld (s1, t1, s3, mpnw1)
    call mpdivd (s3, t2, s1, mpnw1)
    call mpmul (s0, s1, s3, mpnw1)
    call mpeq (s3, s1, mpnw1)
    call mpadd (s1, s2, s3, mpnw1)
    call mpeq (s3, s2, mpnw1)
    if (s1(2) /= 0 .and. s1(3) >= - mpnw1) then
      goto 100
    else
      call mpeq (s2, s1, mpnw1)
      goto 130
    endif
  endif
endif

!   Compute the initial approximation of A ^ (-1/N).

tn = n
call mpmdc (a, t1, n1, mpnw1)
n2 = - n1 / tn
t2 = exp (-1.d0 / tn * (log (t1) + (n1 + tn * n2) * alt))
call mpdmc (t2, n2, s2, mpnw1)
mpnw1 = 5
iq = 0

!   Perform the Newton-Raphson iteration described above with a dynamically
!   changing precision level MPNW (one greater than powers of two).

do k = 1, mq
  if (k > 2) mpnw1 = min (2 * mpnw1 - 2, mpnw) + 1
!  if (k > 2) mpnw1 = min (2 * mpnw1 - 1, mpnw)

110  continue

  call mpnpwr (s2, n, s0, mpnw1)
  call mpmul (a, s0, s1, mpnw1)
  call mpsub (f1, s1, s0, mpnw1)
  call mpmul (s2, s0, s1, mpnw1)
  call mpdivd (s1, tn, s0, mpnw1)
  call mpadd (s2, s0, s1, mpnw1)
  call mpeq (s1, s2, mpnw1)
  if (k == mq - nit .and. iq == 0) then
    iq = 1
    goto 110
  endif
enddo

!   Take the reciprocal to give final result.

call mpdiv (f1, s2, s1, mpnw1)

!   Restore original precision level.

130 continue

call mproun (s1, mpnw)
call mpeq (s1, b, mpnw)

140 continue
return
end subroutine mpnrtr

subroutine mprandr (a, b, mpnw)

!   This returns a pseudorandom number B based the input A.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:)
integer (mpiknd), intent(out):: b(0:)
integer (mpiknd) ia1(0:mpnw+6), ia2(0:mpnw+6), ia3(0:mpnw+6)

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. a(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPRANDR: uninitialized or inadequately sized arrays')
  call mpabrt ( 238)
endif

ia1(0) = mpnw + 7
ia2(0) = mpnw + 7
ia3(0) = mpnw + 7

call mpmuld (a, mprandx, ia1, mpnw + 1)
call mpinfr (ia1, ia2, ia3, mpnw + 1)
call mpeq (ia3, b, mpnw)
return
end subroutine mprandr

subroutine mprealdp (a, b, mpnw)

!   This converts DP argument A to MPR.

implicit none
real (mprknd), intent(in):: a
integer (mpiknd), intent(out):: b(0:)
integer, intent(in):: mpnw

call mpdmc (a, 0, b, mpnw)
return
end subroutine mprealdp

subroutine mprealin (ia, b, mpnw)

!   This converts integer argument IA to MPR.

implicit none
integer, intent(in):: ia
integer (mpiknd), intent(out):: b(0:)
integer, intent(in):: mpnw

call mpdmc (dble (ia), 0, b, mpnw)
return
end subroutine mprealin

subroutine mproun (a, mpnw)

!   This performs rounding and truncation of the MPR number A. It is called
!   by MPNORM, and also by other subroutines when the precision level is
!   modified. It is not intended to be directly called by the user.
!   The parameter MPEXPMX is the absolute value of the largest exponent word
!   allowed for MP numbers (see system parameters at start of this module).

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(inout):: a(0:)
integer i, ia, k, na, n4
integer (mpiknd) a2

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. a(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPROUN: uninitialized or inadequately sized arrays')
  call mpabrt ( 239)
endif

!   Check for initial zeroes.

a2 = a(3)
a(3) = 0
ia = sign (int (1, mpiknd), a(2))
na = min (int (abs (a(2))), mpnw)
n4 = na + 4

if (a(4) == 0) then

!   Find the first nonzero word and shift the entire number left. The length
!   of the result is reduced by the length of the shift.

  do i = 4, n4
    if (a(i+1) /= 0) goto 110
  enddo

  a(2) = 0
  a(3) = 0
  a(4) = 0
  a(5) = 0
  goto 170

110 continue

  k = i - 3

  do i = 3, n4 - k
    a(i+1) = a(i+k+1)
  enddo

  a2 = a2 - k
  na = na - max (k - 2, 0)
  if (k == 2) a(na+4) = 0
endif

!   Perform rounding.

if (na == mpnw) then
  if (a(na+4) >= 0.5d0 * mpbdx) a(na+3) = a(na+3) + 1

!   Release carries as far as necessary due to rounding.

  do i = na + 2, 3, -1
    if (a(i+1) < mpbdx) goto 140
    a(i+1) = a(i+1) - mpbdx
    a(i) = a(i) + 1
  enddo

!   Release of carries due to rounding continued all the way to the start --
!   i.e. number was entirely 9's.

  a(4) = a(3)
  na = 1
  a2 = a2 + 1
endif

140 continue

  if (a(na+3) == 0) then

!   At least the last mantissa word is zero. Find the last nonzero word
!   and adjust the length of the result accordingly.

  do i = na + 2, 3, -1
    if (a(i+1) /= 0) goto 160
  enddo

  a(1) = mpnw
  a(2) = 0
  a(3) = 0
  a(4) = 0
  a(5) = 0
  goto 170

160  continue

  na = i - 2
endif

!   Check for overflow and underflow.

if (a2 < - mpexpmx) then
  write (mpldb, 2)
2 format ('*** MPROUN: Exponent underflow.')
  call mpabrt ( 240)
elseif (a2 > mpexpmx) then
  write (mpldb, 3)
3 format ('*** MPROUN: Exponent overflow.')
  call mpabrt ( 241)
endif

!   Check for zero.

if (a(4) == 0) then
  a(1) = mpnw
  a(2) = 0
  a(3) = 0
  a(4) = 0
  a(5) = 0
else
  a(1) = mpnw
  a(2) = sign (na, ia)
  a(3) = a2
  a(na+4) = 0
  a(na+5) = 0
endif

170  continue

return
end subroutine mproun

integer function mpspacer (ra)

!   This returns the total array space of ra.

implicit none
integer (mpiknd), intent(in):: ra(0:)
mpspacer = ra(0)
return
end function mpspacer

integer function mpsgn (ra)

!   This function returns 1, 0 or -1, depending on whether ra > 0, ra = 0 or ra < 0.

implicit none
integer (mpiknd), intent(in):: ra(0:)
integer ia

! End of declaration

ia = ra(2)
if (ia == 0) then
  mpsgn = 0
elseif (ia > 0) then
  mpsgn = 1
else
  mpsgn = -1
endif
return
end function mpsgn

subroutine mpsqrt (a, b, mpnw)

!   This computes the square root of the MPR number A and returns the result in B.

!   This subroutine employs the following Newton-Raphson iteration, which
!   converges to 1 / Sqrt(A):

!    X_{k+1} = X_k + 0.5 * (1 - X_k^2 * A) * X_k

!   where the multiplication () * X_k is performed with only half of the
!   normal level of precision. These iterations are performed with a
!   working precision level MPNW that is dynamically changed, approximately
!   doubling with each iteration (except that at iteration NIT before the final
!   iteration, the iteration is repeated without doubling the precision, in order
!   to enhance accuracy) . The final iteration is performed as follows
!   (this is due to A. Karp):

!    Sqrt(A) = (A * X_n) + 0.5 * [A - (A * X_n)^2] * X_n  (approx.)

!   where the multiplications A * X_n and [] * X_n are performed with only
!   half of the final level of precision.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:)
integer (mpiknd), intent(out):: b(0:)
integer, parameter:: nit = 3
real (mprknd), parameter:: cl2 = 1.4426950408889633d0
integer ia, iq, k, mpnw1, mq, n, na, nw1, nw2, n2
real (mprknd) t1, t2
integer (mpiknd) s0(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6)

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. b(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPSQRT: uninitialized or inadequately sized arrays')
  call mpabrt ( 242)
endif

ia = sign (int (1, mpiknd), a(2))
na = min (int (abs (a(2))), mpnw)

if (na == 0) then
  b(1) = mpnw
  b(2) = 0
  b(3) = 0
  b(4) = 0
  b(5) = 0
  goto 120
endif
if (ia < 0) then
  write (mpldb, 2)
2 format ('*** MPSQRT: Argument is negative.')
  call mpabrt ( 243)
  return
endif

s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7
s3(0) = mpnw + 7

!   Determine the least integer MQ such that 2 ^ MQ >= MPNW.

t1 = mpnw
mq = cl2 * log (t1) + 1.d0 - mprdfz

!   Compute the initial approximation of 1 / Sqrt(A).

call mpmdc (a, t1, n, mpnw)
n2 = - n / 2
t2 = sqrt (t1 * 2.d0 ** (n + 2 * n2))
t1 = 1.d0 / t2
call mpdmc (t1, n2, s2, mpnw)
call mpdmc (1.d0, 0, s3, mpnw)

mpnw1 = 5
iq = 0
nw1 = mpnw1
nw2 = mpnw1

!   Perform the Newton-Raphson iteration described above with a dynamically
!   changing precision level MPNW (one greater than powers of two).

do k = 1, mq - 1
  if (k > 2) then
    nw1 = mpnw1
    mpnw1 = min (2 * mpnw1 - 2, mpnw) + 1
    nw2 = mpnw1
  endif

100  continue

  call mpmul (s2, s2, s0, nw2)
  call mpmul (a, s0, s1, nw2)
  call mpsub (s3, s1, s0, nw2)
  call mpmul (s2, s0, s1, nw1)
  call mpmuld (s1, 0.5d0, s0, nw1)
  call mpadd (s2, s0, s1, nw2)
  call mpeq (s1, s2, nw2)

  if (k == mq - nit .and. iq == 0) then
    iq = 1
    goto 100
  endif
enddo

!   Perform last iteration using Karp's trick.

nw1 = mpnw1
mpnw1 = min (2 * mpnw1 - 2, mpnw) + 1
nw2 = mpnw1

call mpmul (a, s2, s0, nw1)
call mpmul (s0, s0, s1, nw2)
call mpsub (a, s1, s3, nw2)
call mpmul (s3, s2, s1, nw1)
call mpmuld (s1, 0.5d0, s3, nw1)
call mpadd (s0, s3, s2, nw2)

!   Restore original precision level.

call mproun (s2, mpnw)
call mpeq (s2, b, mpnw)

120 continue

return
end subroutine mpsqrt

subroutine mpsub (a, b, c, mpnw)

!   This routine subtracts MPR numbers A and B to yield C.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:), b(0:)
integer (mpiknd), intent(out):: c(0:)
integer i, nb
integer (mpiknd) s(0:mpnw+5)

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. b(0) < abs (b(2)) + 4 .or. &
  c(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPSUB: uninitialized or inadequately sized arrays')
  call mpabrt ( 244)
endif

nb = min (abs (int (b(2))), mpnw)
s(0) = mpnw + 6
s(1) = mpnw
if (b(2) == 0) then
  s(2) = 0
elseif (b(2) > 0) then
  s(2) = - nb
else
  s(2) = nb
endif

do i = 3, nb + 5
  s(i) = b(i)
enddo

call mpadd (a, s, c, mpnw)

return
end subroutine mpsub

integer function mpwprecr (ra)

!   This returns the working precision of ra.

implicit none
integer (mpiknd), intent(in):: ra(0:)
mpwprecr = ra(1)
return
end function mpwprecr

!   These three subroutines are for real(16) (quad) support:

subroutine mpmqc (a, b, n, mpnw)

!   This returns a quad precision approximation the MPR number A in the form B * 2^n.
!   If IEEE quad floating (128-bit) is not supported on the processor, an error
!   message is output.

implicit none
integer (mpiknd), intent(in):: a(0:)
integer, intent(in):: mpnw
integer, intent(out):: n
integer, parameter:: knd = max (mprknd2, kind(1.0))
real (knd), intent(out):: b
integer na
real (knd) aa

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4) then
  write (mpldb, 1)
1 format ('*** MPMQC: uninitialized or inadequately sized arrays')
  call mpabrt ( 245)
endif

if (mprknd2 < 0) then
  write (mpldb, 2)
2 format ('*** MPMQC: IEEE quad precision is not supported on this processor')
  call mpabrt ( 246)
endif

if (a(2) == 0.d0)  then
  b = 0.d0
  n = 0
  goto 100
endif

na = abs (a(2))
aa = a(4)
if (na >= 2) aa = aa + a(5) / real (mpbdx, knd)
if (na >= 3) aa = aa + a(6) / real (mpbdx, knd)**2

n = mpnbt * a(3)
b = sign (aa, real (a(2), knd))

!   Reduce b to within 1 and 2.

na = log (abs (dble (b))) / log (2.d0) + mprdfz
b = b / 2.d0**na
n = n + na
if (abs (b) < 1.d0) then
  b = 2.d0 * b
  n = n - 1
elseif (abs (b) > 2.d0) then
  b = 0.5d0 * b
  n = n + 1
endif

100  continue
return
end subroutine mpmqc

subroutine mpqmc (a, n, b, mpnw)

!   This routine converts the quad precision number A * 2^N to MPR form in B.
!   If IEEE quad floating (128-bit) is not supported on the processor, an error
!   message is output.

!   NOTE however that the conversion is not fully accurate unless A is an exact
!   binary value.
!   Examples of exact binary values (good): 123456789.q0, 0.25q0, -5.3125q0.
!   Examples of inexact binary values (bad): 0.1q0, 123467.8q0, -3333.3q0.

implicit none
integer, intent(in):: mpnw, n
integer, parameter:: knd = max (mprknd2, kind (1.0))
real (knd), intent(in):: a
integer (mpiknd), intent(out):: b(0:*)
real (knd) aa
integer i, k, n1, n2

! End of declaration

if (mpnw < 4 .or. b(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPQMC: uninitialized or inadequately sized arrays')
  call mpabrt ( 247)
endif

if (mprknd2 < 0) then
  write (mpldb, 2)
2 format ('*** MPQMC: IEEE quad precision is not supported on this processor')
  call mpabrt ( 248)
endif

!   Check for zero.

if (a == 0.d0) then
  b(1) = mpnw
  b(2) = 0
  b(3) = 0
  b(4) = 0
  b(5) = 0
  goto 150
endif
n1 = n / mpnbt
n2 = n - mpnbt * n1
aa = abs (a) * 2.d0 ** n2

!   Reduce AA to within 1 and MPBDX.

if (aa >= mpbdx) then

  do k = 1, 350
    aa = aa / mpbdx
    if (aa < mpbdx) then
      n1 = n1 + k
      goto 120
    endif
 enddo

elseif (aa < 1.d0) then

  do k = 1, 350
    aa = aa * mpbdx
    if (aa >= 1.d0) then
      n1 = n1 - k
      goto 120
    endif
  enddo

endif

!   Store successive sections of AA into B.

120  continue

b(3) = n1
b(4) = int (aa, mpiknd)
aa = mpbdx * (aa - b(3+1))
b(5) = int (aa, mpiknd)
aa = mpbdx * (aa - b(4+1))
b(6) = int (aa, mpiknd)
b(7) = 0
b(8) = 0

do i = 7, 3, -1
  if (b(i+1) /= 0) goto 140
enddo

140  continue

b(1) = mpnw
aa = i - 2
b(2) = sign (aa, a)

150 continue
return
end subroutine mpqmc

subroutine mpqmc90 (a, n, b, mpnw)

!   This routine converts the DP number A * 2^N to MPR form in B.
!   In contrast to mpqmc, this routine only allows 90 significant bits
!   (approximately 27 significant decimal digits) in A. If more nonzero bits
!   are present in A (likely due to inexact binary value), an error is flagged.

!   Examples of exact binary values (good): 123456789.q0, 0.25q0, -5.3125q0.
!   Examples of inexact binary values (bad): 0.1q0, 123467.8q0, -3333.3q0.

implicit none
integer, intent(in):: mpnw, n
integer, parameter:: knd = max (mprknd2, kind (1.0))
real (knd), intent(in):: a
integer (mpiknd), intent(out):: b(0:)
real (knd) t2

! End of declaration

if (mpnw < 4 .or. b(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPQMC40: uninitialized or inadequately sized arrays')
  call mpabrt ( 249)
endif

!   Check whether A has more than 90 significant bits (actually whether
!   the trailing 23 bits are zero).

t2 = mpmask23 (a)
if (t2 == abs (a)) then
  call mpqmc (a, n, b, mpnw)
else
  write (mpldb, 2) a
2 format ('*** MPQMC40: QP value has more than 90 significant bits:'/ &
  1p,d50.35/'and thus very likely represents an unintended loss of accuracy.'/ &
  'Fix the issue, or else use functions mpprod, mpquot, mprealq or mprealqm.'/ &
  'See documentation for details.')
  call mpabrt ( 250)
endif

return
end subroutine mpqmc90

! ***  The following are the extra-high precision multiply routines:

subroutine mpfftcr (is, m, n, nsq, x, y)

!   This performs an N-point complex-to-real FFT, where N = 2^M. X is the
!   double complex input array, and Y is the double precision output array.
!   The array X is used as a scratch array in MPFFT1, and so is overwritten.
!   X and Y must be dimensioned as shown below. IS is the sign of the FFT.
!   The arrays MPUU1 and MPUU2 must have been initialized by calling MPINIFFT.
!   This routine is not intended to be called directly by the user.

implicit none
integer, intent(in):: is, m, n, nsq
real (mprknd), intent(out):: y(n)
complex (mprknd), intent(inout):: x(n/2+nsq*mpnsp1+1)
integer k, ku, mx, n1, n2, n4
complex (mprknd) dc1(n/2), ai, a1, a2, x1, x2

mx = mpuu1(1)

!   Check if input parameters are invalid.

if ((is /= 1 .and. is /= -1) .or. m < 3 .or. m > mx) then
  write (mpldb, 1)  is, m, mx
1 format ('*** MPFFTCR: Either the UU arrays have not been initialized'/ &
  'or else one of the input parameters is invalid', 3i5)
  call mpabrt ( 251)
endif

n1 = 2 ** (m / 2)
n2 = n / 2
n4 = n / 4
ai = cmplx (0.d0, 1.d0, mprknd)

!   Construct the input to MPFFT1.

dc1(1) = 0.5d0 * cmplx (real (x(1) + x(n2+1), mprknd), &
  real (x(1) - x(n2+1), mprknd), mprknd)
if (is == 1) then
  dc1(n4+1) = conjg (x(n4+1))
else
  dc1(n4+1) = x(n4+1)
endif
ku = n2

if (is == 1) then
  do k = 2, n4
    x1 = x(k)
    x2 = conjg (x(n2+2-k))
    a1 = x1 + x2
    a2 = ai * mpuu1(k+ku) * (x1 - x2)
    dc1(k) = 0.5d0 * (a1 + a2)
    dc1(n2+2-k) = 0.5d0 * conjg (a1 - a2)
  enddo
else
  do k = 2, n4
    x1 = x(k)
    x2 = conjg (x(n2+2-k))
    a1 = x1 + x2
    a2 = ai * conjg (mpuu1(k+ku)) * (x1 - x2)
    dc1(k) = 0.5d0 * (a1 + a2)
    dc1(n2+2-k) = 0.5d0 * conjg (a1 - a2)
  enddo
endif

!   Perform a normal N/2-point FFT on DC1.

call mpfft1 (is, m - 1, n1, n2 / n1, dc1, x)

!   Copy DC1 to Y such that DC1(k) = Y(2k-1) + i Y(2k).

do k = 1, n / 2
  y(2*k-1) = real (dc1(k), mprknd)
  y(2*k) = aimag (dc1(k))
enddo

return
end subroutine mpfftcr

subroutine mpfftrc (is, m, n, nsq, x, y)

!   This performs an N-point real-to-complex FFT, where N = 2^M. X is the
!   double precision input array, and Y is the double complex output array.
!   The arrays MPUU1 and MPUU2 must have been initialized by calling MPINIFFT.
!   This routine is not intended to be called directly by the user.

implicit none
integer, intent(in):: is, m, n, nsq
real (mprknd), intent(in):: x(n)
complex (mprknd), intent(out):: y(n/2+nsq*mpnsp1+1)
integer k, ku, mx, n1, n2, n4
complex (mprknd) dc1(n/2), ai, a1, a2, z1, z2

mx = mpuu1(1)

!   Check if input parameters are invalid.

if ((is /= 1 .and. is /= -1) .or. m < 3 .or. m > mx) then
  write (mpldb, 1)  is, m, mx
1 format ('*** MPFFTRC: either the UU arrays have not been initialized'/ &
  'or else one of the input parameters is invalid',3i5)
  call mpabrt ( 252)
endif

n1 = 2 ** (m / 2)
n2 = n / 2
n4 = n / 4
ai = cmplx (0.d0, -1.d0, mprknd)

!   Copy X to DC1 such that DC1(k) = X(2k-1) + i X(2k).

do k = 1, n2
  dc1(k) = cmplx (x(2*k-1), x(2*k), mprknd)
enddo

!   Perform a normal N/2-point FFT on DC1.

call mpfft1 (is, m - 1, n1, n2 / n1, dc1, y)

!   Reconstruct the FFT of X.

y(1) = cmplx (2.d0 * (real (dc1(1), mprknd) + aimag (dc1(1))), &
  0.d0, mprknd)
if (is == 1) then
  y(n4+1) = 2.d0 * dc1(n4+1)
else
  y(n4+1) = 2.d0 * conjg (dc1(n4+1))
endif
y(n2+1) = cmplx (2.d0 * (real (dc1(1), mprknd) - aimag (dc1(1))), &
  0.d0, mprknd)
ku = n2

if (is == 1) then
  do k = 2, n4
    z1 = dc1(k)
    z2 = conjg (dc1(n2+2-k))
    a1 = z1 + z2
    a2 = ai * mpuu1(k+ku) * (z1 - z2)
    y(k) = a1 + a2
    y(n2+2-k) = conjg (a1 - a2)
  enddo
else
  do k = 2, n4
    z1 = dc1(k)
    z2 = conjg (dc1(n2+2-k))
    a1 = z1 + z2
    a2 = ai * conjg (mpuu1(k+ku)) * (z1 - z2)
    y(k) = a1 + a2
    y(n2+2-k) = conjg (a1 - a2)
  enddo
endif

return
end subroutine mpfftrc

subroutine mpfft1 (is, m, n1, n2, x, y)

!   This routine performs a complex-to-complex FFT. IS is the sign of the
!   transform, N = 2^M is the size of the transform. N1 = 2^M1 and N2 = 2^M2,
!   where M1 and M2 are defined as below. X is the input and output array,
!   and Y is a scratch array. X must have at N, and Y at least N + N1*MPNSP1,
!   double complex cells. The arrays MPUU1 and MPUU2 must have been
!   initialized by calling MPINIFFT. This routine is not intended to be called
!   directly by the user.

!   This employs the two-pass variant of the "four-step" FFT. See the
!   article by David H. Bailey in J. of Supercomputing, March 1990, p. 23-35.

implicit none
integer, intent(in):: is, m, n1, n2
complex (mprknd), intent(inout):: x(n1,n2)
complex (mprknd), intent(out):: y(n2+mpnsp1,n1)
integer i, iu, j, j2, k, ku, m1, m2, nr1, nr2
complex (mprknd) z1(mpnrow+mpnsp1,n1), z2(mpnrow+mpnsp1,n1)

m1 = (m + 1) / 2
m2 = m - m1
nr1 = min (n1, mpnrow)
nr2 = min (n2, mpnrow)
ku = mpuu2(m)

do i = 0, n1 - 1, nr1

!   Copy NR1 rows of X (treated as a N1 x N2 complex array) into Z1.

  do j = 1, n2
    do k = 1, nr1
      z1(k,j) = x(i+k,j)
    enddo
  enddo

!   Perform NR1 FFTs, each of length N2.

  call mpfft2 (is, nr1, m2, n2, z1, z2)

!   Multiply the resulting NR1 x N2 complex block by roots of unity and
!   store transposed into the appropriate section of Y.

  iu = i + ku - n1 - 1
  if (is == 1) then
    do j = 1, n2
      do k = 1, nr1
        y(j,i+k) = mpuu2(iu+k+j*n1) * z1(k,j)
      enddo
    enddo
  else
    do j = 1, n2
      do k = 1, nr1
        y(j,i+k) = conjg (mpuu2(iu+k+j*n1)) * z1(k,j)
      enddo
    enddo
  endif
enddo

do i = 0, n2 - 1, nr2

!   Copy NR2 rows of the Y array into Z2.

  do j = 1, n1
    do k = 1, nr2
      z2(k,j) = y(i+k,j)
    enddo
  enddo

!   Perform NR2 FFTs, each of length N1.

  call mpfft2 (is, nr2, m1, n1, z2, z1)

!   Copy NR2 x N1 complex block back into X array. It's a little more
!   complicated if M is odd.

  if (mod (m, 2) == 0) then
    do j = 1, n1
      do k = 1, nr2
        x(i+k,j) = z2(k,j)
      enddo
    enddo
  else
    do j = 1, n1 / 2
      j2 = 2 * j - 1

      do k = 1, nr2
        x(i+k,j) = z2(k,j2)
        x(i+k+n2,j) = z2(k,j2+1)
      enddo
    enddo
  endif
enddo

return
end subroutine mpfft1

subroutine mpfft2 (is, ns, m, n, x, y)

!   This performs NS simultaneous N-point complex-to-complex FFTs, where
!   N = 2^M. X is the input and output array, and Y is a scratch array.
!   The arrays MPUU1 and MPUU2 must have been initialized by calling MPINIFFT.
!   This routine is not intended to be called directly by the user.

implicit none
integer, intent(in):: is, ns, m, n
complex (mprknd), intent(inout):: x(mpnrow+mpnsp1,n)
complex (mprknd), intent(out):: y(mpnrow+mpnsp1,n)
integer i, j, l

!   Perform the second variant of the Stockham FFT.

do l = 1, m, 2
  call mpfft3 (is, l, ns, m, n, x, y)
  if (l == m) goto 100
  call mpfft3 (is, l + 1, ns, m, n, y, x)
enddo

goto 110

!   Copy Y to X.

100 continue

do j = 1, n
  do i = 1, ns
    x(i,j) = y(i,j)
  enddo
enddo

110 continue

return
end subroutine mpfft2

subroutine mpfft3 (is, l, ns, m, n, x, y)

!   This performs the L-th iteration of the second variant of the Stockham FFT
!   on the NS vectors in X. X is input/output, and Y is a scratch array.
!   The arrays MPUU1 and MPUU2 must have been initialized by calling MPINIFFT.
!   This routine is not intended to be called directly by the user.

implicit none
integer, intent(in):: is, l, ns, m, n
complex (mprknd), intent(inout):: x(mpnrow+mpnsp1,n)
complex (mprknd), intent(out):: y(mpnrow+mpnsp1,n)
integer i, i11, i12, i21, i22, j, k, li, lj, lk, ku, n1
complex (mprknd) u1, x1, x2

!   Set initial parameters.

n1 = n / 2
lk = 2 ** (l - 1)
li = 2 ** (m - l)
lj = 2 * lk
ku = li + 1

do i = 0, li - 1
  i11 = i * lk + 1
  i12 = i11 + n1
  i21 = i * lj + 1
  i22 = i21 + lk
  if (is == 1) then
    u1 = mpuu1(i+ku)
  else
    u1 = conjg (mpuu1(i+ku))
  endif

  do k = 0, lk - 1
    do j = 1, ns
      x1 = x(j,i11+k)
      x2 = x(j,i12+k)
      y(j,i21+k) = x1 + x2
      y(j,i22+k) = u1 * (x1 - x2)
    enddo
  enddo
enddo

return
end subroutine mpfft3

subroutine mpinifft (mpnw)

!   This computes the root of unity arrays UU1 and UU2, which are required by
!   the FFT routines, and places this data in the proper arrays defined in
!   module MPFUNA. MPNW is the largest precision level (in words) that will be
!   subsequently used for this run.

implicit none
integer, intent(in):: mpnw
real (mprknd), parameter:: cl2 = 1.4426950408889633d0
integer i, iu, j, k, ku, ln, m, mm, mm1, mm2, mq, nn, nn1, nn2, nq, nwds
real (mprknd) d1
real (mprknd) pi, t1, ti, tpn

!  Determine sizes for FFT arrays. Three words are added to mpnw, since many
!  routines in MPFUND in particular increase the working precision upon entry.

nwds = mpnw + 3
d1 = 2.d0 * (nwds + 1)
m = cl2 * log (d1) + 1.d0 - mprdfz
mq = m + 2
nq = 2 ** mq

if (mq + nq > mplfftx) then
  write (6, 1) mq + nq
1 format ('*** MPINIFFT: Insufficient space for arrays mpuu1 and mpuu2.'/ &
  'At least',i12,' double complex cells must be allocated for each of'/ &
  'these arrays in module mpfuna. See documentation for details.')
  call mpabrt ( 253)
endif

mpuu1(1) = mq
ku = 2
ln = 1
pi = acos (-1.d0)

do j = 1, mq
  t1 = pi / ln

  do i = 0, ln - 1
    ti = i * t1
    mpuu1(i+ku) = cmplx (cos (ti), sin (ti), mprknd)
  enddo

  ku = ku + ln
  ln = 2 * ln
enddo

! write (6, 2) ku - 1
! 2 format ('MPINIFFT: Size of table mpuu1 =',i10)

ku = mq + 1
mpuu2(1) = mq

do k = 2, mq
  mpuu2(k) = cmplx (0.d0, 0.d0, mprknd)
enddo

do k = 2, mq - 1
  mpuu2(k) = ku
  mm = k
  nn = 2 ** mm
  mm1 = (mm + 1) / 2
  mm2 = mm - mm1
  nn1 = 2 ** mm1
  nn2 = 2 ** mm2
  tpn = 2.d0 * pi / nn

  do j = 0, nn2 - 1
    do i = 0, nn1 - 1
      iu = ku + i + j * nn1
      t1 = tpn * i * j
      mpuu2(iu) = cmplx (cos (t1), sin (t1), mprknd)
    enddo
  enddo

  ku = ku + nn
enddo

! write (6, 3) ku - 1
! 3 format ('MPINIFFT: Size of table mpuu2 =',i10)

return
end subroutine mpinifft

subroutine mplconv (iq, n, nsq, a, b, c)

!   This computes the linear convolution of A and B, returning the result
!   in C. If IQ is 1, then it is presumed B = A; if IQ = 2, then A /= B.
!   NSQ is a spacing parameter, which should be set to more than sqrt (3*n).

implicit none
integer, intent(in):: iq, n, nsq
real (mprknd), intent(in):: a(n), b(n)
real (mprknd), intent(out):: c(2*n)
real (mprknd), parameter:: cl2 = 1.4426950408889633d0, ffterrmx = 0.375d0
integer i, m1, m2, n1, n2, n4, nm
real (mprknd) c0
real (mprknd) an, d1(8*n+2), d2(8*n+2), d3(8*n+2), t1, t2
complex (mprknd) dc1(4*n+nsq*mpnsp1+3), dc2(4*n+nsq*mpnsp1+3)

t1 = n
m1 = cl2 * log (t1) + 1.d0 - mprdfz
n1 = 2 ** m1
m2 = m1 + 1
n2 = 2 * n1
n4 = 2 * n2
nm = min (2 * n, n2)

if (abs (iq) == 1) then

!   Compute the square of A -- only one forward FFT is needed.

  do i = 1, n
    d1(i) = a(i)
  enddo

  do i = n + 1, n2
    d1(i) = 0.d0
  enddo

!   Perform a forward real-to-complex FFT on the vector in A.

  call mpfftrc (1, m2, n2, nsq, d1, dc1)

!   Square the resulting complex vector.

  do i = 1, n1 + 1
    dc1(i) = dc1(i) ** 2
  enddo
else

!   Compute the product of A and B -- two forward FFTs are needed.

  do i = 1, n
    d1(i) = a(i)
    d2(i) = b(i)
  enddo

  do i = n + 1, n2
    d1(i) = 0.d0
    d2(i) = 0.d0
  enddo

!   Perform forward real-to-complex FFTs on the vectors in A and B.

  call mpfftrc (1, m2, n2, nsq, d1, dc1)
  call mpfftrc (1, m2, n2, nsq, d2, dc2)

!   Multiply the resulting complex vectors.

  do i = 1, n1 + 1
    dc1(i) = dc1(i) * dc2(i)
  enddo
endif

!   Perform an inverse complex-to-real FFT on the resulting data.

call mpfftcr (-1, m2, n2, nsq, dc1, d3)

!   Divide by N4 and round to nearest whole number.

an = 1.d0 / n4
c0 = 0.d0

do i = 1, nm
  t1 = an * d3(i)
  t2 = anint (t1)
  c(i) = t2
  c0 = max (c0, abs (t2 - t1))
enddo

if (c0 > ffterrmx) then
  write (6, 1) c0
1 format ('*** MPLCONV: excessive rounding error =',f12.6)
  call mpabrt ( 254)
endif

return
end subroutine mplconv

subroutine mpmulx (a, b, c, mpnw)

!   This routine multiplies MP numbers A and B to yield the MP product C,
!   using a FFT-convolution technique. Before calling MPMULX, the arrays
!   UU1 and UU2 must be initialized by calling MPINIFFT. For modest levels
!   of precision, use MPMUL.

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: a(0:), b(0:)
integer (mpiknd), intent(out):: c(0:)
integer i, ia, ib, na, nb, nc, nn, nx
integer (mpiknd) d(0:mpnw+8), i0, i1, i2
real (mprknd) d1(0:4*mpnw+20), d2(0:4*mpnw+20), d3(0:8*mpnw+40)

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. b(0) < abs (b(2)) + 4 .or. &
  c(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPMULX: uninitialized or inadequately sized arrays')
  call mpabrt ( 255)
endif

ia = sign (int (1, mpiknd), a(2))
ib = sign (int (1, mpiknd), b(2))
na = min (int (abs (a(2))), mpnw)
nb = min (int (abs (b(2))), mpnw)
nc = min (na + nb, mpnw)
nn = 4 * max (na, nb)
nx = sqrt (4.d0 * nn) + mprdfz

!   Divide each word of A into four 15-bit chunks.

do i = 0, na - 1
  i1 = a(i+4)
  i2 = shifta (i1, 45)
  d1(4*i) = i2
  i1 = i1 - shiftl (i2, 45)
  i2 = shifta (i1, 30)
  d1(4*i+1) = i2
  i1 = i1 - shiftl (i2, 30)
  i2 = shifta (i1, 15)
  d1(4*i+2) = i2
  i1 = i1 - shiftl (i2, 15)
  d1(4*i+3) = i1
enddo

do i = 4 * na, nn - 1
  d1(i) = 0.d0
enddo

!   Divide each word of B into four 15-bit chunks.

do i = 0, nb - 1
  i1 = b(i+4)
  i2 = shifta (i1, 45)
  d2(4*i) = i2
  i1 = i1 - shiftl (i2, 45)
  i2 = shifta (i1, 30)
  d2(4*i+1) = i2
  i1 = i1 - shiftl (i2, 30)
  i2 = shifta (i1, 15)
  d2(4*i+2) = i2
  i1 = i1 - shiftl (i2, 15)
  d2(4*i+3) = i1
enddo

do i = 4 * nb, nn - 1
  d2(i) = 0.d0
enddo

!   Perform linear convolution.

  call mplconv (2, nn, nx, d1, d2, d3)

!   Release carries.

i0 = 0

do i = min (4 * nc + 16, 2 * nn - 1), 0, -1
  i0 = i0 + d3(i)
  i1 = shifta (i0, 15)
  i2 = i0 - shiftl (i1, 15)
  d1(i) = i2
  i0 = i1
enddo

!  Recombine words, with proper offset.

d(0) = 0
d(1) = 0
d(2) = 0
d(3) = 0
d(4) = shiftl (i0, 45) + shiftl (int (d1(0), mpiknd), 30) &
  + shiftl (int (d1(1), mpiknd), 15) + int (d1(2), mpiknd)

do i = 1, nc + 3
  d(i+4) = shiftl (int (d1(4*i-1), mpiknd), 45) + shiftl (int (d1(4*i), mpiknd), 30) &
    + shiftl (int (d1(4*i+1), mpiknd), 15) + int (d1(4*i+2), mpiknd)
enddo

d(0) = mpnw + 6
d(1) = mpnw
d(2) = sign (nc, ia * ib)
d(3) = a(3) + b(3) + 1

!   Fix up the result.

d1(0) = mpnw + 6
call mpnorm (d, c, mpnw)

190 continue

return
end subroutine mpmulx

end module mpfunb
