!*****************************************************************************

!  program tquadgs

!  Revision date:  24 Sep 2021

!  AUTHOR:
!   David H. Bailey
!   Lawrence Berkeley National Lab (retired) and University of California, Davis
!   Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!   All software in this package (c) 2021 David H. Bailey.
!   By downloading or using this software you agree to the copyright, disclaimer
!   and license agreement in the accompanying file DISCLAIMER.txt.

!  NOTE: This program runs very long (800 seconds with MPFUN-MPFR or 2400 seconds
!   with MPFUN-Fort, on the author's computer).

!  DESCRIPTION:
!   This program demonstrates the quadrature routine 'quadgs', which employs the
!   Gaussian quadrature algorithm for numerical integration. quadgs is suitable to
!   integrate a regular function on a finite interval. It can also be used for an
!   integral on the entire real line, by making a variable transformation or
!   computing a modified set of weights and abscissas, as is done below.

!   Gaussian quadrature is typically faster than tanh-sinh quadrature (see the
!   program tquad.f90 from the author) for entirely regular functions, but it is
!   not effective for functions where the function itself or one of its higher
!   order derivatives has a singularity at either endpoint. Also, the computational
!   cost of generating the abscissa-weight data for Gaussian quadrature is much
!   greater (typically 100-300 times greater) than for tanh-sinh quadrature. On the
!   other hand, the Gaussian abscissa-weight data can be computed once to a certain
!   precision and stored in a file for future use, as is optionally done below. For
!   additional details see:

!   David H. Bailey, Xiaoye S. Li and Karthik Jeyabalan, "A comparison of
!   three high-precision quadrature schemes," Experimental Mathematics, 
!   vol. 14 (2005), no. 3, pg. 317-329, preprint available at
!   http://www.davidhbailey.com/dhbpapers/quadrature-em.pdf.

!   All of these routines are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   Specific instructions:

!   Before calling quadgs, call initqgs to initialize the arrays wkgs, xkgs,
!   wkgu and xkgu. Set x1 and x2 (the limits of integration) in executable
!   statements in the calling program. Then call quadgs with the integrand
!   function name as the first argument, and other arguments as shown below.
!   The problems performed below are a subset of the problems performed in
!   tquad.f90; the problem numbers here correspond to those in tquad.f90. See
!   the initial comments below in initqgs and quadgs for additional details.

!   For problems over the entire real line, call quadgs with wkgu and xkgu as
!   the last two arguments, and set x1 = -one and x2 = one. See problems 15
!   through 18 below (corresponding to functions fun15 through fun18).

!   For some integrand functions, significantly more accurate quadrature results
!   can be achieved by computing x1 and x2 to higher precision (ndp2 digits or
!   nwds2 words) in the calling program. x1 and x2 are used in quadgs for
!   computation of the scaled abscissa, which is in turn passed to the function
!   definition for use in initial subtractions or other numerically sensitive
!   operations involving the input argument. The function evaluation itself
!   should be performed with standard precision (ndp1 digits or nwds1 words)
!   for faster run times.

!   The following integer parameters are set in the parameter statement below;
!   all are default (4-byte) integer:
!   idata  0: Compute abscissas and weights from scratch, but do not write file.
!          1: Compute abscissas and weights from scratch, and write to file
!             before performing the quadrature problems. Default = 1.
!          2: Read precomputed abscissas and weights from a file before
!             performing the quadrature problems.
!          The file name is "gauss-mm-nnnn.dat", where "mm" is the value of
!          nq1 and "nnnn" is the value of ndp1.
!   ndp1   Primary ("standard") precision in digits; this is the the target
!          accuracy of quadrature results; default = 500.
!   ndp2   Secondary ("high") precision in digits; default = 2*ndp1.
!   neps1  Log10 of the primary tolerance;  default = - ndp1.
!   neps2  Log10 of the secondary tolerance; default = -ndp2.
!   nq1    Max number of phases in quadrature routine; adding 1 increases
!          (possibly doubles) the number of accurate digits in the result,
!          but also roughly doubles the run time. nq1 must be at least 3.
!          Default = 11.
!   nq2    Space parameter for wkgs, xkgs, wkgu, xkgu arrays.
!          Increase nq2 if directed by a message produced in initqgs.
!          Default = 6 * 2^nq1.
!   nwds1  standard precision in words; default = int (ndp1 / mpdpw + 2).
!   nwds2  High precision in words; default = int (ndp2 / mpdpw + 2).

!   The functions to be integrated (fun01, fun02, etc) must be declared
!   type (mp_real) and external in the calling program, and must be defined in
!   separate function subprograms. See examples below.

program tquadgs
use mpmodule
implicit none
integer i, i1, i2, idata, ndp1, ndp2, neps1, neps2, nq1, nq2, nwds1, nwds2, n1
parameter (idata = 1, ndp1 = 500, ndp2 = 1000, neps1 = -ndp1, neps2 = -ndp2, &
  nq1 = 11, nq2 = 6 * 2 ** nq1 + 100, nwds1 = int (ndp1 / mpdpw + 2), &
  nwds2 = int (ndp2 / mpdpw + 2))
character*32 chr32
real (mprknd) d1, second, tm0, tm1, tm2
type (mp_real) err, errmx, quadgs, fun01, fun02, fun03, fun04, &
  fun15, fun16, fun17, fun18, one, t1, t2, wkgs(-1:nq2), xkgs(-1:nq2), &
  wkgu(-1:nq2), xkgu(-1:nq2), zero
type (mp_real) mppic, mpl02, x1, x2
external fun01, fun02, fun03, fun04, fun15, fun16, fun17, fun18, &
  quadgs, second

!   Check to see if default precision is high enough.

if (ndp2 > mpipl) then
  write (6, '("Increase default precision in module MPFUNF.")')
  stop
endif

!   Compute pi and log(2) to high precision (nwds2 words).

zero = mpreal (0.d0, nwds2)
one = mpreal (1.d0, nwds2)
mppic = mppi (nwds2)
mpl02 = mplog2 (nwds2)
errmx = zero

write (6, 1) nq1, ndp1, ndp2, neps1, neps2
1 format ('Quadgs test:  Quadlevel =',i6/'Digits1 =',i6,'  Digits2 =',i6, &
  '  Epsilon1 =',i6,'  Epsilon2 =',i6)

if (idata >= 1) then

!    Open file for weights and abscissas.

  write (chr32, '("gauss-",i2.2,"-",i4.4,".dat")') nq1, ndp1
  open (11, file = chr32)
  rewind 11
endif

tm0 = second ()
if (idata <= 1) then

!   Generate abscissas and weights from scratch.

  call initqgs (nq1, nq2, nwds1, wkgs, xkgs, wkgu, xkgu)

  if (idata == 1) then

!   Write abscissas and weights to file.

    n1 = dble (xkgs(-1))

    do i = -1, n1
      write (11, '(2i6)') i, n1
      call mpwrite (11, ndp1 + 30, ndp1 + 10, wkgs(i), xkgs(i))
      call mpwrite (11, ndp1 + 30, ndp1 + 10, wkgu(i), xkgu(i))
    enddo
  endif
elseif (idata == 2) then

!   Read abscissas and weights from file.

  read (11, '(2i6)') i1, i2
  call mpread (11, wkgs(-1), xkgs(-1), nwds1)
  call mpread (11, wkgu(-1), xkgu(-1), nwds1)
  n1 = dble (xkgs(-1))
  if (i1 /= -1 .or. i2 /= n1) stop

  do i = 0, n1
    read (11, '(2i6)') i1, i2
    if (i1 /= i .or. i2 /= n1) stop
    call mpread (11, wkgs(i), xkgs(i), nwds1)
    call mpread (11, wkgu(i), xkgu(i), nwds1)
  enddo
endif
if (idata >= 1) close (11)
tm1 = second ()
tm2 = tm1 - tm0
write (6, 2) tm1 - tm0
2 format ('Quadrature initialization completed: cpu time =',f12.4)

!   Begin quadrature tests.

write (6, 11)
11 format (/'Continuous functions on finite intervals:'//&
  'Problem 1: Int_0^1 t*log(1+t) dt = 1/4'/)
x1 = zero
x2 = one
tm0 = second ()
t1 = quadgs (fun01, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkgs, xkgs)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
3 format ('Quadrature completed: CPU time =',f12.6/'Result =')
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = mpreal (0.25d0, nwds1)
err = t2 - t1
errmx = max (abs (err), errmx)
call mpdecmd (err, d1, n1)
write (6, 4) d1, n1
4 format ('Actual error =',f10.6,'x10^',i6)

write (6, 12)
12 format (/'Problem 2: Int_0^1 t^2*arctan(t) dt = (pi - 2 + 2*log(2))/12'/)
x1 = zero
x2 = one
tm0 = second ()
t1 = quadgs (fun02, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkgs, xkgs)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = (mppic - 2.d0 + 2.d0 * mpl02) / 12.d0
err = t2 - t1
errmx = max (abs (err), errmx)
call mpdecmd (err, d1, n1)
write (6, 4) d1, n1

write (6, 13)
13 format (/'Problem 3: Int_0^(pi/2) e^t*cos(t) dt = 1/2*(e^(pi/2) - 1)'/)
x1 = zero
x2 = 0.5d0 * mppic
tm0 = second ()
t1 = quadgs (fun03, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkgs, xkgs)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = 0.5d0 * (exp (0.5d0 * mppic) - 1.d0)
err = t2 - t1
errmx = max (abs (err), errmx)
call mpdecmd (err, d1, n1)
write (6, 4) d1, n1

write (6, 14)
14 format (/ &
  'Problem 4: Int_0^1 arctan(sqrt(2+t^2))/((1+t^2)sqrt(2+t^2)) dt = 5*Pi^2/96'/)
x1 = zero
x2 = one
tm0 = second ()
t1 = quadgs (fun04, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkgs, xkgs)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = 5.d0 * mppic**2 / 96.d0
err = t2 - t1
errmx = max (abs (err), errmx)
call mpdecmd (err, d1, n1)
write (6, 4) d1, n1

write (6, 25)
25 format (/ &
   'Functions on the entire real line:'// &
   'Problem 15: Int_-inf^inf 1/(1+t^2) dt = Pi'/)
x1 = -one
x2 = one
tm0 = second ()
t1 = quadgs (fun15, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkgu, xkgu)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = mppic
err = t2 - t1
errmx = max (abs (err), errmx)
call mpdecmd (err, d1, n1)
write (6, 4) d1, n1

write (6, 26)
26 format (/'Problem 16: Int_-inf^inf 1/(1+t^4) dt = Pi/Sqrt(2)'/)
x1 = -one
x2 = one
tm0 = second ()
t1 = quadgs (fun16, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkgu, xkgu)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = mppic / sqrt (mpreal (2.d0, nwds1))
err = t2 - t1
errmx = max (abs (err), errmx)
call mpdecmd (err, d1, n1)
write (6, 4) d1, n1

write (6, 27)
27 format (/'Problem 17: Int_-inf^inf e^(-t^2/2) dt = sqrt (2*Pi)'/)
x1 = -one
x2 = one
tm0 = second ()
t1 = quadgs (fun17, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkgu, xkgu)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = sqrt (2.d0 * mppic)
err = t2 - t1
errmx = max (abs (err), errmx)
call mpdecmd (err, d1, n1)
write (6, 4) d1, n1

write (6, 28)
28 format (/'Problem 18: Int_-inf^inf e^(-t^2/2) cos(t) dt = sqrt (2*Pi/e)'/)
x1 = -one
x2 = one
tm0 = second ()
t1 = quadgs (fun18, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkgu, xkgu)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = sqrt (2.d0 * mppic / exp (mpreal (1.d0, nwds1)))
err = t2 - t1
errmx = max (abs (err), errmx)
call mpdecmd (err, d1, n1)
write (6, 4) d1, n1

call mpdecmd (errmx, d1, n1)
write (6, 91) tm2, d1, n1
91 format (/'Total CPU time =',f12.6/'Max abs error =',f10.6,'e',i6)
if (errmx < mpreal (10.d0, nwds1) ** (neps1 + 5)) then
  write (6, '(a)') 'ALL TESTS PASSED'
else
  write (6, '(a)') 'ONE OR MORE TESTS FAILED'
endif

stop
end

!   Function definitions:

function fun01 (t, nwds1, nwds2)

!   fun01(t) = t * log(1+t)

use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) fun01, t1
type (mp_real) t

t1 = mpreal (t, nwds1)
fun01 = t1 * log (1.d0 + t1)
return
end

function fun02 (t, nwds1, nwds2)

!   fun02(t) = t^2 * arctan(t)

use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) fun02, t1
type (mp_real) t

t1 = mpreal (t, nwds1)
fun02 = t1 ** 2 * atan (t1)
return
end

function fun03 (t, nwds1, nwds2)

!   fun03(t) = e^t * cos(t)

use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) fun03, t1
type (mp_real) t

t1 = mpreal (t, nwds1)
fun03 = exp(t1) * cos(t1)
return
end

function fun04 (t, nwds1, nwds2)

!   fun04(t) = arctan(sqrt(2+t^2))/((1+t^2)sqrt(2+t^2))

use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) fun04, t1, t2
type (mp_real) t

t1 = mpreal (t, nwds1)
t2 = sqrt (2.d0 + t1**2)
fun04 = atan(t2) / ((1.d0 + t1**2) * t2)
return
end

function fun15 (t, nwds1, nwds2)

!   fun15(t) = 1 / (1 + t^2)

use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) t, fun15

fun15 = 1.d0 / (1.d0 + t**2)
return
end

function fun16 (t, nwds1, nwds2)

!   fun16(t) = 1 / (1 + t^4)

use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) t, fun16

fun16 = 1.d0 / (1.d0 + t**4)
return
end

function fun17 (t, nwds1, nwds2)

!   fun17(t) = exp(-1/2*t^2)

use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) t, fun17

fun17 = exp (-0.5d0 * t**2)
return
end

function fun18 (t, nwds1, nwds2)

!   fun18(t) = exp(-1/2*t^2) * cos(t)

use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) t, fun18

fun18 = exp (-0.5d0 * t**2) * cos (t)
return
end

!   Gaussian quadrature routines:

subroutine initqgs (nq1, nq2, nwds1, wkgs, xkgs, wkgu, xkgu)

!   David H Bailey    24 Sep 2021

!   This subroutine initializes the quadrature arrays wkgs and xkgs for standard
!   Gaussian quadrature, and also wkgu and xkgu for quadrature over the real line.
!   It employs a Newton iteration scheme with a dynamic precision level that
!   approximately doubles with each iteration.

!   The wkgu and xkgu arrays are computed from wkgs and xkgs based on the
!   transformation t = tan (pi/2 * x), which transforms an integral on
!   (-infinity, infinity) to an integral on (-1, 1). In particular,
!     xkgu(i) = tan (pi/2 * xkgs(i))
!     wkgu(i) = pi/2 * wkgs(i) / cos (pi/2 * xkgs(i))^2

!   Both initqgs and quadgs are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   Input arguments; all are default (4-byte) integer:
!   nq1      Level parameter for data in wkgs, xkgs, wkgu, xkgu arrays.
!   nq2      Space parameter for wkgs, xkgs, wkgu, xkgu arrays.
!   nwds1    Primary precision level, in mantissa words.

!   Output arguments:
!   wkgs     Output array of standard weights; type (mp_real).
!   xkgs     Output array of standard abscissas; type (mp_real).
!   wkgu     Output array of real line weights; type (mp_real).
!   xkgu     Output array of real line abscissas; type (mp_real).

use mpmodule
implicit none
integer i, ik0, iprint, j, j1, k, n, ndebug, nq1, nq2, nwds1, nwp
real (mprknd) dpi
parameter (ik0 = 100, iprint = 1, ndebug = 2, dpi = 3.141592653589793238d0)
type (mp_real) eps, pi2, r, t1, t2, t3, t4, t5, wkgs(-1:nq2), xkgs(-1:nq2), &
  wkgu(-1:nq2), xkgu(-1:nq2), zero

if (ndebug >= 1) then
  write (6, 1)
1 format ('initqgs: Gaussian quadrature initialization')
endif

pi2 = 0.5d0 * mppi (nwds1)
zero = mpreal (0.d0, nwds1)
wkgs(-1) = mpreal (dble (nq1), nwds1)
xkgs(-1) = zero
wkgs(0) = zero
xkgs(0) = zero
wkgs(1) = mpreal (dble (nq1), nwds1)
xkgs(1) = mpreal (dble (ik0), nwds1)
wkgu(-1) = mpreal (dble (nq1), nwds1)
xkgu(-1) = zero
wkgu(0) = zero
xkgu(0) = zero
wkgu(1) = mpreal (dble (nq1), nwds1)
xkgu(1) = mpreal (dble (ik0), nwds1)
i = ik0

do j = 2, ik0
  wkgs(j) = zero
  xkgs(j) = zero
  wkgu(j) = zero
  xkgu(j) = zero
enddo

do k = 1, nq1
  if (ndebug >= 2 .and. mod (k, iprint) == 0) write (6, '(3i8)') k, nq2
  n = 3 * 2 ** (k + 1)

  do j = 1, n / 2

!   Set working precision = 4 words, and compute a DP estimate of the root.

    nwp = 4
    eps = mpreal (2.d0 ** (mpnbt - nwp * mpnbt), nwds1)
    r = mpreald (cos ((dpi * (j - 0.25d0)) / (n + 0.5d0)), nwp)

!   Compute the j-th root of the n-degree Legendre polynomial using Newton's
!   iteration.

100 continue

!   Perform the next 11 lines with working precision nwp words.

    t1 = mpreal (1.d0, nwp)
    t2 = mpreal (0.d0, nwp)
    r = mpreal (r, nwp)

    do j1 = 1, n
      t3 = t2
      t2 = t1
      t1 = ((2 * j1 - 1) * r * t2 - (j1 - 1) * t3) / j1
    enddo

    t4 = n * (r * t1 - t2) / (r ** 2  - 1.d0)
    t5 = r
    r = r - t1 / t4
    
!   Once convergence is achieved at nwp = 4, then start doubling (almost) the
!   working precision level at each iteration until full precision is reached.

    if (nwp == 4) then
      if (abs (r - t5) > eps) goto 100
      nwp = min (2 * nwp - 1, nwds1)
      goto 100
    elseif (nwp < nwds1) then
      nwp = min (2 * nwp - 1, nwds1)
      goto 100
    endif

    i = i + 1
    if (i > nq2) goto 110
    xkgs(i) = r
    t4 = n * (r * t1 - t2) / (r ** 2  - 1.d0)
    wkgs(i) = 2.d0 / ((1.d0 - r ** 2) * t4 ** 2)
    xkgu(i) = tan (pi2 * r)
    wkgu(i) = pi2 * wkgs(i) / cos (pi2 * r)**2
  enddo

!   Save i (starting index for the next block) in the first 100 elements of 
!   xkgs and xkgu.

  xkgs(k+1) = mpreal (dble (i), nwds1)
  xkgu(k+1) = mpreal (dble (i), nwds1)
enddo

!   Save the size of the arrays in index -1 of xkgs and xkgu.

xkgs(-1) = mpreal (dble (i), nwds1)
xkgu(-1) = mpreal (dble (i), nwds1)
if (ndebug >= 2) then
  write (6, 2) i
2 format ('initqgs: Table spaced used =',i8)
endif
goto 130

110 continue

write (6, 3) nq2
3 format ('initqgs: Table space parameter is too small; value =',i8)
stop

130 continue

return
end

function quadgs (fun, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkgs, xkgs)

!   David H Bailey  16 Jul 2021

!   This routine computes the integral of the function fun on the interval
!   (x1, x2) with a target tolerance of 10^neps1. The quadrature level is
!   progressively increased (approximately doubling the work with each level)
!   until level nq1 has been performed or the target tolerance has been met.
!   Prior to calling quadgs, the wkgs and xkgs arrays must first be initialized
!   by calling initqgs. If quadgs outputs the message "Terms too large", adjust
!   nq1 and neps2 as necessary in the call to initqgs.

!   To use quadgs for the entire line, set x1 = -1 and x2 = 1 in executable
!   statements in the calling program, then call quadgs with wkgu and xkgu as
!   the last two arguments. See comments in main program above.

!   For some integrand functions, significantly more accurate quadrature results
!   can be achieved by computing x1 and x2 to higher precision (ndp2 digits or
!   nwds2 words) in the calling program. x1 and x2 are used in quadgs for
!   computation of the scaled abscissa, which is in turn passed to the function
!   definition for use in initial subtractions or other numerically sensitive
!   operations involving the input argument. The function evaluation itself
!   should be performed with standard precision (ndp1 digits or nwds1 words)
!   for faster run times. See the function definitions of fun06, fun07, fun09
!   and fun10 for examples.

!   Both initqgs and quadgs are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   Input arguments; unless indicated otherwise all are default (4-byte) integer:
!   fun      Name of function to be integrated; type (mp_real).
!   x1       Lower limit of integration interval; type (mp_real).
!   x2       Upper limit of integration interval; type (mp_real).
!   nq1      Level parameter for data in wkgs and xkgs arrays.
!   nq2      Space parameter for wkgs, xkgs arrays.
!   nwds1    Primary precision level, in mantissa words.
!   nwds2    Secondary precision level, in mantissa words.
!   neps1    Primary epsilon level (eps1 = 10^neps1).
!   wkgs     Precomputed array of weights; type (mp_real).
!   xkgs     Precomputed array of abscissas; type (mp_real).

!   Output argument:
!   quadgs   Quadrature result; type (mp_real).

use mpmodule
implicit none
integer i, k, j, n, ndebug, nq1, nq2, &
  neps1, nwds1, nwds2
double precision d1, d2, d3, dfrac, dplog10
parameter (dfrac = 100.d0, ndebug = 2)
type (mp_real) ax, bx, c10, eps1, epsilon1, err, fun, &
  quadgs, tsum, s1, s2, s3, t1, t2, tw1, tw2, twmx, wkgs(-1:nq2), xkgs(-1:nq2), &
  x1, x2, xx1, xx2
external fun, dplog10

!  These two lines are performed in high precision (nwds2 words).

ax = 0.5d0 * (x2 - x1)
bx = 0.5d0 * (x2 + x1)

!  The remaining initialization is performed in standard precision (nwds1 words).

epsilon1 = dfrac * mpreal (10.d0, nwds1) ** neps1
s1 = mpreal (0.d0, nwds1)
s2 = mpreal (0.d0, nwds1)
c10 = mpreal (10.d0, nwds1)

if (wkgs(-1) < dble (nq1)) then
  write (6, 1) nq1
1 format ('quadgs: quadrature arrays have not been initialized; nq1 =',i6)
  goto 140
endif

do k = 1, nq1
  n = 3 * 2 ** (k + 1)
  s3 = s2
  s2 = s1
  twmx = mpreal (0.d0, nwds1)
  tsum = mpreal (0.d0, nwds1)
  i = dble (xkgs(k))
  
  do j = 1, n / 2
    i = i + 1

!   These two lines are performed in high precision.

    xx1 = - ax * xkgs(i) + bx
    xx2 = ax * xkgs(i) + bx

    t1 = fun (xx1, nwds1, nwds2)
    tw1 = t1 * wkgs(i)

    if (j + k > 2) then
      t2 = fun (xx2, nwds1, nwds2)
      tw2 = t2 * wkgs(i)
    else
      t2 = mpreal (0.d0, nwds1)
      tw2 = mpreal (0.d0, nwds1)
    endif

    tsum = tsum + tw1 + tw2
    twmx = max (twmx, abs (tw1), abs (tw2))
  enddo

!   Compute s1 = current integral approximation and err = error estimate.
!   Tsum is the sum of all tw1 and tw2 from the loop above.
!   Twmx is the largest absolute value of tw1 and tw2 from the loop above.

  s1 =  mpreal (ax, nwds1) * tsum
  eps1 = twmx * epsilon1
  d1 = dplog10 (abs ((s1 - s2) / s1))
  d2 = dplog10 (abs ((s1 - s3) / s1))
  d3 = dplog10 (eps1) - 1.d0

  if (k <= 2) then
    err = mpreal (1.d0, nwds1)
  elseif (d1 .eq. -999999.d0) then
    err = mpreal (0.d0, nwds1)
  else
    err = c10 ** nint (min (0.d0, max (d1 ** 2 / d2, 2.d0 * d1, d3)))
  endif
  
!   Output current integral approximation and error estimate, to 60 digits.

  if (ndebug >= 2) then
    write (6, 2) k, nq1, nint (dplog10 (abs (err)))
2   format ('quadgs: Iteration',i3,' of',i3,'; est error = 10^',i7, &
      '; approx value =')
    call mpwrite (6, 80, 60, s1)
  endif

  if (k >= 3 .and. err < eps1) then
    write (6, 4) nint (dplog10 (abs (err)))
4   format ('quadgs: Estimated error = 10^',i7)
    goto 140
  endif
enddo

140 continue

quadgs = s1
return
end

function dplog10 (a)

!   For input MPM value a, this routine returns a DP approximation to log10 (a).

use mpmodule
implicit none
integer ia
real (mprknd) da, dplog10
type (mp_real) a

call mpdecmd (a, da, ia)
if (da == 0.d0) then
  dplog10 = -999999.d0
else
  dplog10 = log10 (abs (da)) + ia
endif

100 continue
return
end
