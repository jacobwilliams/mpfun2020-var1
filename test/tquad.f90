!*****************************************************************************

!  program tquad

!  Revision date:  7 Jan 2023

!  AUTHOR:
!   David H. Bailey
!   Lawrence Berkeley National Lab (retired) and University of California, Davis
!   Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!   All software in this package (c) 2023 David H. Bailey.
!   By downloading or using this software you agree to the copyright, disclaimer
!   and license agreement in the accompanying file DISCLAIMER.txt.

!  DESCRIPTION:
!   This program demonstrates three quadrature (numerical integration) routines:

!   quadts: Implements the tanh-sinh quadrature scheme of Takahashi and Mori,
!     for functions on a finite interval such as (0,1).
!   quades: Implements the exp-sinh scheme, a variation of tanh-sinh well-suited
!     for functions on a semi-infinite interval such as (0, +infinity).
!   quadss: Implements the sinh-sinh scheme, a variation of tanh-sinh well-suited
!     for functions on the entire real line.

!   These schemes have highly favorable properties for quadrature calculations,
!   notably the fact that the cost of computing weight-abscissa pairs increases
!   only linearly with the number of evaluation points, instead of quadratically
!   as with Gaussian quadrature. For example, the cost of computing abscissas and
!   weights for a quadrature calculation to, say, 500 digits is approximately 300
!   times faster with one of the above schemes, compared to Gaussian quadrature.
!   Also, quadts and quades usually work well even when the function has a blow-up
!   singularity or infinite derivative at an endpoint, whereas Gaussian quadrature
!   typically performs very poorly for such functions. On the other hand, for
!   entirely regular integrand functions on finite intervals, Gaussian quadrature
!   is usually faster in performing the quadrature itself. For full details, see:

!   David H. Bailey, Xiaoye S. Li and Karthik Jeyabalan, "A comparison of
!   three high-precision quadrature schemes," Experimental Mathematics, 
!   vol. 14 (2005), no. 3, pg. 317-329, preprint available at
!   http://www.davidhbailey.com/dhbpapers/quadrature-em.pdf.

!   All of the routines in this program are 100% THREAD SAFE -- all requisite
!   parameters and arrays are passed through subroutine arguments. 

!   Here are some specific instructions for the individual routines:

!   quadts:  

!   Before calling quadts, call initqts to initialize the arrays wkts and xkts. 
!   Set x1 and x2 (the limits of integration) in executable statements in the
!   calling program. Then call quadts with the integrand function name as the
!   first argument, and other arguments as shown below. quadts is illustrated in
!   problems 1 through 10 below (corresponding to functions fun01 through fun10).
!   See the initial comments below in initqts and quadts for additional details.

!   For some integrand functions, significantly more accurate quadrature results
!   can be achieved by computing x1 and x2 to higher precision (ndp2 digits or
!   nwds2 words) in the calling program. x1 and x2 are used in quadts for
!   computation of the scaled abscissa, which is in turn passed to the function
!   definition for use in initial subtractions or other numerically sensitive
!   operations involving the input argument. The function evaluation itself
!   should be performed with standard precision (ndp1 digits or nwds1 words)
!   for faster run times. See the function definitions of fun06, fun07, fun09
!   and fun10 for examples. 

!   quades:

!   Before calling quades, call initqes to initialize the arrays wkes and xkes. 
!   Set x1 (the lower limit of integration) in an executable statement in the
!   calling program. Then call quades with the integrand function name as the
!   first argument, and other arguments as shown below. quades is illustrated in
!   problems 11 through 14 below (corresponding to functions fun11 through fun14).
!   See the initial comments below in initqes and quades for additional details.

!   The comment above about computing x1 to higher precision also applies here.

!   quadss:

!   Before calling quadss, call initqss to initialize the arrays wkss and xkss. 
!   Then call quadss with the integrand function name as the first argument, and
!   other arguments as shown below. quadss is illustrated in problems 15 through
!   18 below (corresponding to functions fun15 through fun18). See the initial
!   comments below in initqss and quadss for additional details.

!   The following integer parameters are set in the parameter statement below;
!   all are default (4-byte) integer:
!   ndp1   Primary ("standard") precision in digits; this is the target
!          accuracy of quadrature results; default = 500.
!   ndp2   Secondary ("high") precision in digits; default = 2*ndp1.
!   neps1  Log10 of the primary tolerance; default = - ndp1.
!   neps2  Log10 of the secondary tolerance; default = -ndp2.
!   nq1    Max number of phases in quadrature routine; adding 1 increases
!          (possibly doubles) the number of accurate digits in the result,
!          but also roughly doubles the run time. Must be at least 3.
!          Default = 11.
!   nq2    Space parameter for wkts, xkts, wkes, xkes, wkss, xkss arrays.
!          Increase nq2 if directed by a message produced in initqts.
!          Default = 12 * 2^nq1.
!   nwds1  Standard precision in words; default = int (ndp1 / mpdpw + 2).
!   nwds2  High precision in words; default = int (ndp2 / mpdpw + 2).

!   The functions to be integrated (fun01, fun02, etc) must be declared
!   type (mp_real) and external in the calling program, and must be defined in
!   separate function subprograms. See examples below.

program tquad
use mpmodule
implicit none
integer, parameter:: ndp1 = 500, ndp2 = 1000, neps1 = -ndp1, neps2 = -ndp2, &
  nq1 = 11, nq2 = 12 * 2 ** nq1, nwds1 = int (ndp1 / mpdpw + 2), &
  nwds2 = int (ndp2 / mpdpw + 2)
integer n1
real (mprknd) d1, tm0, tm1, tm2
type (mp_real) err, errmx, mppic, mpl02, one, t1, t2, wkes(-1:nq2), x1, x2, &
  xkes(-1:nq2), wkss(-1:nq2), xkss(-1:nq2), wkts(-1:nq2), xkts(-1:nq2), zero
type (mp_real), external:: quades, quadss, quadts, fun01, fun02, fun03, &
  fun04, fun05, fun06, fun07, fun08, fun09, fun10, fun11, fun12, fun13, &
  fun14, fun15, fun16, fun17, fun18
real (mprknd), external:: second

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
1 format ('Quadts test:  Quadlevel =',i6/'Digits1 =',i6,'  Digits2 =',i6, &
  '  Epsilon1 =',i6,'  Epsilon2 =',i6)

!   Initialize tables of weights and abscissas.

tm0 = second ()
call initqts (nq1, nq2, nwds1, neps2, wkts, xkts)
call initqes (nq1, nq2, nwds1, neps2, wkes, xkes)
call initqss (nq1, nq2, nwds1, neps2, wkss, xkss)
tm1 = second ()
tm2 = tm1 - tm0
write (6, 2) tm1 - tm0
2 format ('Quadrature initialization completed: cpu time =',f12.6)

!   Begin quadrature tests.

write (6, 11)
11 format (/'Continuous functions on finite intervals:'//&
  'Problem 1: Int_0^1 t*log(1+t) dt = 1/4'/)
x1 = zero
x2 = one
tm0 = second ()
t1 = quadts (fun01, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkts, xkts)
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
4 format ('Actual error =',f10.6,'e',i6)

write (6, 12)
12 format (/'Problem 2: Int_0^1 t^2*arctan(t) dt = (pi - 2 + 2*log(2))/12'/)
x1 = zero
x2 = one
tm0 = second ()
t1 = quadts (fun02, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkts, xkts)
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
t1 = quadts (fun03, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkts, xkts)
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
t1 = quadts (fun04, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkts, xkts)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = 5.d0 * mppic**2 / 96.d0
err = t2 - t1
errmx = max (abs (err), errmx)
call mpdecmd (err, d1, n1)
write (6, 4) d1, n1

write (6, 15)
15 format (/&
  'Continuous functions on finite intervals, but non-diff at an endpoint:'// &
  'Problem 5: Int_0^1 sqrt(t)*log(t) dt = -4/9'/)
x1 = zero
x2 = one
tm0 = second ()
t1 = quadts (fun05, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkts, xkts)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = mpreal (-4.d0, nwds1) / 9.d0
err = t2 - t1
errmx = max (abs (err), errmx)
call mpdecmd (err, d1, n1)
write (6, 4) d1, n1

write (6, 16)
16 format (/'Problem 6: Int_0^1 sqrt(1-t^2) dt = pi/4'/)
x1 = zero
x2 = one
tm0 = second ()
t1 = quadts (fun06, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkts, xkts)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = 0.25d0 * mppic
err = t2 - t1
errmx = max (abs (err), errmx)
call mpdecmd (err, d1, n1)
write (6, 4) d1, n1

write (6, 17)
17 format (/&
  'Functions on finite intervals with integrable singularity at an endpoint:'//&
  'Problem 7: Int_0^1 sqrt(t)/sqrt(1-t^2) dt = 2*sqrt(pi)*gamma(3/4)/gamma(1/4)'/)
x1 = zero
x2 = one
tm0 = second ()
t1 = quadts (fun07, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkts, xkts)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = 2.d0 * sqrt (mpreal (mppic, nwds1)) * gamma (mpreal (0.75d0, nwds1)) &
  / gamma (mpreal (0.25d0, nwds1))
err = t2 - t1
errmx = max (abs (err), errmx)
call mpdecmd (err, d1, n1)
write (6, 4) d1, n1

write (6, 18)
18 format (/'Problem 8: Int_0^1 log(t)^2 dt = 2'/)
x1 = zero
x2 = one
tm0 = second ()
t1 = quadts (fun08, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkts, xkts)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = mpreal (2.d0, nwds1)
err = t2 - t1
errmx = max (abs (err), errmx)
call mpdecmd (err, d1, n1)
write (6, 4) d1, n1

write (6, 19)
19 format (/'Problem 9: Int_0^(pi/2) log(cos(t)) dt = -pi*log(2)/2'/)
x1 = zero
x2 = 0.5d0 * mppic
tm0 = second ()
t1 = quadts (fun09, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkts, xkts)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = -0.5d0 * mppic * mpl02
err = t2 - t1
errmx = max (abs (err), errmx)
call mpdecmd (err, d1, n1)
write (6, 4) d1, n1

write (6, 20)
20 format (/'Problem 10: Int_0^(pi/2) sqrt(tan(t)) dt = pi*sqrt(2)/2'/)
x1 = zero
x2 = 0.5d0 * mppic
tm0 = second ()
t1 = quadts (fun10, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkts, xkts)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = 0.5d0 * mppic * sqrt (mpreal (2.d0, nwds1))
err = t2 - t1
errmx = max (abs (err), errmx)
call mpdecmd (err, d1, n1)
write (6, 4) d1, n1

write (6, 21)
21 format (/&
  'Functions on a semi-infinite interval:'//&
  'Problem 11: Int_1^inf 1/(1+t^2) dt = pi/4'/)
x1 = one
tm0 = second ()
t1 = quades (fun11, x1, nq1, nq2, nwds1, nwds2, neps1, wkes, xkes)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = 0.25d0 * mppic
err = t2 - t1
errmx = max (abs (err), errmx)
call mpdecmd (err, d1, n1)
write (6, 4) d1, n1

write (6, 22)
22 format (/'Problem 12: Int_0^inf e^(-t)/sqrt(t) dt = sqrt(pi)'/)
x1 = zero
tm0 = second ()
t1 = quades (fun12, x1, nq1, nq2, nwds1, nwds2, neps1, wkes, xkes)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = sqrt (mppic)
err = t2 - t1
errmx = max (abs (err), errmx)
call mpdecmd (err, d1, n1)
write (6, 4) d1, n1

write (6, 23)
23 format (/'Problem 13: Int_0^inf e^(-t^2/2) dt = sqrt(pi/2)'/)
x1 = zero
tm0 = second ()
t1 = quades (fun13, x1, nq1, nq2, nwds1, nwds2, neps1, wkes, xkes)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = sqrt (0.5d0 * mppic)
err = t2 - t1
errmx = max (abs (err), errmx)
call mpdecmd (err, d1, n1)
write (6, 4) d1, n1

write (6, 24)
24 format (/'Problem 14: Int_pi^inf e^(-t)*cos(t) dt = -1/2 * exp(-pi)'/)
x1 = mppic
tm0 = second ()
t1 = quades (fun14, x1, nq1, nq2, nwds1, nwds2, neps1, wkes, xkes)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = -0.5d0 * exp (- mpreal (mppic, nwds1))
err = t2 - t1
errmx = max (abs (err), errmx)
call mpdecmd (err, d1, n1)
write (6, 4) d1, n1

write (6, 25)
25 format (/ &
   'Functions on the entire real line:'// &
   'Problem 15: Int_-inf^inf 1/(1+t^2) dt = Pi'/)
tm0 = second ()
t1 = quadss (fun15, nq1, nq2, nwds1, neps1, wkss, xkss)
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
tm0 = second ()
t1 = quadss (fun16, nq1, nq2, nwds1, neps1, wkss, xkss)
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
tm0 = second ()
t1 = quadss (fun17, nq1, nq2, nwds1, neps1, wkss, xkss)
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
tm0 = second ()
t1 = quadss (fun18, nq1, nq2, nwds1, neps1, wkss, xkss)
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
end program tquad

!   Function definitions:

type (mp_real) function fun01 (t, nwds1, nwds2)

!   fun01(t) = t * log(1+t)

use mpmodule
implicit none
integer, intent(in):: nwds1, nwds2
type (mp_real), intent(in):: t
type (mp_real) t1

t1 = mpreal (t, nwds1)
fun01 = t1 * log (1.d0 + t1)
return
end function fun01

type (mp_real) function fun02 (t, nwds1, nwds2)

!   fun02(t) = t^2 * arctan(t)

use mpmodule
implicit none
integer, intent(in):: nwds1, nwds2
type (mp_real), intent(in):: t
type (mp_real) t1

t1 = mpreal (t, nwds1)
fun02 = t1 ** 2 * atan (t1)
return
end function fun02

type (mp_real) function fun03 (t, nwds1, nwds2)

!   fun03(t) = e^t * cos(t)

use mpmodule
implicit none
integer, intent(in):: nwds1, nwds2
type (mp_real), intent(in):: t
type (mp_real) t1

t1 = mpreal (t, nwds1)
fun03 = exp(t1) * cos(t1)
return
end function fun03

type (mp_real) function fun04 (t, nwds1, nwds2)

!   fun04(t) = arctan(sqrt(2+t^2)) / ((1+t^2) * sqrt(2+t^2))

use mpmodule
implicit none
integer, intent(in):: nwds1, nwds2
type (mp_real), intent(in):: t
type (mp_real) t1, t2

t1 = mpreal (t, nwds1)
t2 = sqrt (2.d0 + t1**2)
fun04 = atan(t2) / ((1.d0 + t1**2) * t2)
return
end function fun04

type (mp_real) function fun05 (t, nwds1, nwds2)

!    fun05(t) = sqrt(t)*log(t)

use mpmodule
implicit none
integer, intent(in):: nwds1, nwds2
type (mp_real), intent(in):: t
type (mp_real) t1

t1 = mpreal (t, nwds1)
fun05 = sqrt (t1) * log (t1)
return
end function fun05

type (mp_real) function fun06 (t, nwds1, nwds2)

!    fun06(t) = sqrt(1-t^2)

use mpmodule
implicit none
integer, intent(in):: nwds1, nwds2
type (mp_real), intent(in):: t
type (mp_real) t1, t2

!   t2 must be calculated using high precision (nwds2 words), but standard
!   precision is fine after the subtraction.

t1 = mpreal (t, nwds1)
t2 = mpreal (1.d0 - t1**2, nwds1)
fun06 = sqrt (t2)
return
end function fun06

type (mp_real) function fun07 (t, nwds1, nwds2)

!   fun07(t) = sqrt(t) / sqrt(1-t^2)

use mpmodule
implicit none
integer, intent(in):: nwds1, nwds2
type (mp_real), intent(in):: t
type (mp_real) t1, t2

!   t2 must be calculated using high precision (nwds2 words), but standard
!   precision is fine after the subtraction.

t1 = mpreal (t, nwds1)
t2 = mpreal (1.d0 - t, nwds1)
fun07 = sqrt (t1) / sqrt (t2 * (1.d0 + t1))
return
end function fun07

type (mp_real) function fun08 (t, nwds1, nwds2)

!   fun08(t) = log(t)^2

use mpmodule
implicit none
integer, intent(in):: nwds1, nwds2
type (mp_real), intent(in):: t
type (mp_real) t1

t1 = mpreal (t, nwds1)
fun08 = log (t1) ** 2
return
end function fun08

type (mp_real) function fun09 (t, nwds1, nwds2)

!   fun09(t) = log(cos(t))

use mpmodule
implicit none
integer, intent(in):: nwds1, nwds2
type (mp_real), intent(in):: t
type (mp_real) pi, t1, t2, t3, t4

!   t2 must be calculated using high precision (nwds2 words), but standard
!   precision is fine after the subtraction.

t1 = mpreal (t, nwds1)
pi = mppi (nwds2)
t3 = mpreal (0.25d0 * pi, nwds1)
t2 = mpreal (0.5d0 * pi - t, nwds1)

if (t1 < t3) then
  t4 = cos (t1)
else
  t4 = sin (t2)
endif
fun09 = log (t4)

return
end function fun09

type (mp_real) function fun10 (t, nwds1, nwds2)

!   fun10(t) = sqrt(tan(t))

use mpmodule
implicit none
integer, intent(in):: nwds1, nwds2
type (mp_real), intent(in):: t
type (mp_real) pi, t1, t2, t3


!   t2 must be calculated using high precision (nwds2 words), but standard
!   precision is fine after the subtraction.

t1 = mpreal (t, nwds1)
pi = mppi (nwds2)
t3 = mpreal (0.25d0 * pi, nwds1)
t2 = mpreal (0.5d0 * pi - t, nwds1)

if (t1 < t3) then
  fun10 = sqrt (tan (t1))
else
  fun10 = 1.d0 / sqrt (tan (t2))
endif
return
end function fun10

type (mp_real) function fun11 (t, nwds1, nwds2)

!   fun11(t) = 1/(1 + t^2)

use mpmodule
implicit none
integer, intent(in):: nwds1, nwds2
type (mp_real), intent(in):: t
type (mp_real) t1

t1 = mpreal (t, nwds1)
fun11 = 1.d0 / (1.d0 + t1 ** 2)
return
end function fun11

type (mp_real) function fun12 (t, nwds1, nwds2)

!   fun12(t) = e^(-t) / sqrt(t)

use mpmodule
implicit none
integer, intent(in):: nwds1, nwds2
type (mp_real), intent(in):: t
type (mp_real) t1

t1 = mpreal (t, nwds1)
fun12 = exp (-t1) / sqrt (t1)
return
end function fun12

type (mp_real) function fun13 (t, nwds1, nwds2)

!   fun13(t) = e^(-t^2/2)

use mpmodule
implicit none
integer, intent(in):: nwds1, nwds2
type (mp_real), intent(in):: t
type (mp_real) t1

t1 = mpreal (t, nwds1)
fun13 = exp (-0.5d0 * t1 ** 2)
return
end function fun13

type (mp_real) function fun14 (t, nwds1, nwds2)

!   fun14(t) = e^(-t) * cos(t)

use mpmodule
implicit none
integer, intent(in):: nwds1, nwds2
type (mp_real), intent(in):: t
type (mp_real) t1

t1 = mpreal (t, nwds1)
fun14 = exp (-t1) * cos (t1)
return
end function fun14

type (mp_real) function fun15 (t, nwds1)

!   fun15(t) = 1 / (1 + t^2)

use mpmodule
implicit none
integer, intent(in):: nwds1
type (mp_real), intent(in):: t
type (mp_real) t1

t1 = mpreal (t, nwds1)
fun15 = 1.d0 / (1.d0 + t1**2)
return
end function fun15

type (mp_real) function fun16 (t, nwds1)

!   fun16(t) = 1 / (1 + t^4)

use mpmodule
implicit none
integer, intent(in):: nwds1
type (mp_real), intent(in):: t
type (mp_real) t1

t1 = mpreal (t, nwds1)
fun16 = 1.d0 / (1.d0 + t1**4)
return
end function fun16

type (mp_real) function fun17 (t, nwds1)

!   fun17(t) = exp(-1/2*t^2)

use mpmodule
implicit none
integer, intent(in):: nwds1
type (mp_real), intent(in):: t
type (mp_real) t1

t1 = mpreal (t, nwds1)
fun17 = exp (-0.5d0 * t1**2)
return
end function fun17

type (mp_real) function fun18 (t, nwds1)

!   fun18(t) = exp(-1/2*t^2) * cos(t)

use mpmodule
implicit none
integer, intent(in):: nwds1
type (mp_real), intent(in):: t
type (mp_real) t1

t1 = mpreal (t, nwds1)
fun18 = exp (-0.5d0 * t1**2) * cos (t1)
return
end function fun18

!   Tanh-sinh routines:

subroutine initqts (nq1, nq2, nwds1, neps2, wkts, xkts)

!   David H Bailey   7 Jan 2023

!   This subroutine initializes the quadrature arrays xkts and wkts for quadts.
!   If initqts outputs the message "Table space parameter is too small", adjust
!   nq2 in calling program. Also, if quadts outputs the message "Terms too large",
!   adjust nq1 and neps2 as necessary in the call to initqts. The argument neps2
!   controls termination of the loop below, which ends when wkts(k) < 10^(neps2).

!   The wkts and xkts arrays are computed based on the transformation
!   t = tanh (pi/2 * sinh (x)).  Note however that xkts contains one minus the
!   conventional abscissas, in order to conserve precision. See comments below.

!   Both initqts and quadts are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   Input arguments; all are default (4-byte) integer:
!   nq1      Level parameter for data in wkts and xkts arrays.
!   nq2      Space parameter for wkts and xkts arrays.
!   nwds1    Primary precision level, in mantissa words.
!   neps2    Secondary epsilon level (eps2 = 10^neps2).

!   Output arguments:
!   wkts     Output array of weights; type (mp_real).
!   xkts     Output array of abscissas; type (mp_real).

use mpmodule
implicit none
integer, intent(in):: nq1, nq2, nwds1, neps2
type (mp_real), intent(out):: wkts(-1:nq2), xkts(-1:nq2)
integer, parameter:: iprint = 1024, ndebug = 2
integer k
type (mp_real) eps2, h, p2, t1, t2, t3, t4, u1, u2

write (6, 1)
1 format ('initqts: Tanh-sinh quadrature initialization')

eps2 = mpreal (10.d0, nwds1) ** neps2
p2 = 0.5d0 * mppi (nwds1)
h = mpreal (0.5d0 ** nq1, nwds1)
wkts(-1) = mpreal (dble (nq1), nwds1)

do k = 0, nq2
  if (ndebug >= 2 .and. mod (k, iprint) == 0) write (6, '(2i8)') k, nq2
  t1 = mpreal (dble (k) * h, nwds1)

!   xkts(k) = 1 - tanh (u1) = 1 /(e^u1 * cosh (u1))
!   wkts(k) = u2 / cosh (u1)^2
!   where u1 = pi/2 * cosh (t1), u2 = pi/2 * sinh (t1)

  t2 = exp (t1)
  u1 = 0.5d0 * p2 * (t2 + 1.d0 / t2)
  u2 = 0.5d0 * p2 * (t2 - 1.d0 / t2)
  t3 = exp (u2)
  t4 = 0.5d0 * (t3 + 1.d0 / t3)
  xkts(k) = 1.d0 / (t3 * t4)
  wkts(k) = u1 / t4 ** 2

  if (wkts(k) < eps2) goto 100
enddo

write (6, 2) nq2
2 format ('initqts: Table space parameter is too small; value =',i8)
stop

100 continue

xkts(-1) = mpreal (dble (k), nwds1)
if (ndebug >= 2) then
  write (6, 3) k
3 format ('initqts: Table spaced used =',i8)
endif

return
end subroutine initqts

type (mp_real) function quadts (fun, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkts, xkts)

!   David H Bailey  7 Jan 2023

!   This routine computes the integral of the function fun on the interval
!   (x1, x2) with a target tolerance of 10^neps1. The quadrature level is
!   progressively increased (approximately doubling the work with each level)
!   until level nq1 has been performed or the target tolerance has been met.
!   Prior to calling quadts, the wkts and xkts arrays must first be initialized
!   by calling initqts. If quadts outputs the message "Terms too large", adjust
!   nq1 and neps2 as necessary in the call to initqts.

!   For some integrand functions, significantly more accurate quadrature results
!   can be achieved by computing x1 and x2 to higher precision (ndp2 digits or
!   nwds2 words) in the calling program. x1 and x2 are used in quadts for
!   computation of the scaled abscissa, which is in turn passed to the function
!   definition for use in initial subtractions or other numerically sensitive
!   operations involving the input argument. The function evaluation itself
!   should be performed with standard precision (ndp1 digits or nwds1 words)
!   for faster run times. See the function definitions of fun06, fun07, fun09
!   and fun10 for examples.

!   Both initqts and quadts are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   Input arguments; unless indicated otherwise all are default (4-byte) integer:
!   fun      Name of function to be integrated; type (mp_real).
!   x1       Lower limit of integration interval; type (mp_real).
!   x2       Upper limit of integration interval; type (mp_real).
!   nq1      Level parameter for data in wkts and xkts arrays.
!   nq2      Space parameter for wkts and xkts arrays.
!   nwds1    Primary precision level, in mantissa words.
!   nwds2    Secondary precision level, in mantissa words.
!   neps1    Primary epsilon level (eps1 = 10^neps1).
!   wkts     Precomputed array of weights; type (mp_real).
!   xkts     Precomputed array of abscissas; type (mp_real).

!   Output argument:
!   quadts   Quadrature result; type (mp_real).

use mpmodule
implicit none
type (mp_real), external:: fun
type (mp_real), intent(in):: x1, x2, wkts(-1:nq2), xkts(-1:nq2)
integer, intent(in):: nq1, nq2, nwds1, nwds2, neps1
integer, parameter:: izx = 5, ndebug = 2
real (mprknd), external:: dplog10
integer i, ip(0:100), iz1, iz2, k, k1, k2, n, nqq1
logical log1, log2
real (mprknd) d1, d2, d3, d4
type (mp_real) c10, eps1, eps2, epsilon1, err, h, &
  tsum, s1, s2, s3, t1, t2, tw1, tw2, twi1, twi2, twmx, &
  ax, bx, xki, xt1, xx1, xx2

!  These two lines are performed in high precision (nwds2 words).

ax = 0.5d0 * (x2 - x1)
bx = 0.5d0 * (x2 + x1)

!  The remaining initialization is performed in standard precision (nwds1 words).

epsilon1 = mpreal (10.d0, nwds1) ** neps1
tsum = mpreal (0.d0, nwds1)
s1 = mpreal (0.d0, nwds1)
s2 = mpreal (0.d0, nwds1)
h = mpreal (1.d0, nwds1)
c10 = mpreal (10.d0, nwds1)

if (wkts(-1) < dble (nq1)) then
  write (6, 1) nq1
1 format ('quadts: Quadrature arrays have not been initialized; nq1 =',i6)
  goto 140
endif
nqq1 = dble (wkts(-1))
n = dble (xkts(-1))

do k = 0, nqq1
  ip(k) = 2 ** k
enddo

do k = 1, nq1
  h = 0.5d0 * h
  s3 = s2
  s2 = s1
  k1 = ip(nqq1-k)
  k2 = ip(nqq1-k+1)
  iz1 = 0
  iz2 = 0
  twmx = mpreal (0.d0, nwds1)

!   Evaluate function at level k in x, avoiding unnecessary computation.

  do i = 0, n, k1  
    if (mod (i, k2) /= 0 .or. k == 1) then

!   These next few lines, which scale the abscissas, must be performed in
!   high precision (nwds2 words) to ensure full accuracy in the quadrature
!   results, even though the abscissas xkts(i) were computed in standard precision.

      xki = xkts(i)
      xt1 = 1.d0 - mpreal (xki, nwds2)
      xx1 = - ax * xt1 + bx
      xx2 = ax * xt1 + bx
      log1 = xx1 > x1
      log2 = xx2 < x2      

!   The remaining computations are performed in standard precision (nwds1 words).

      if (log1 .and. iz1 < izx) then
        t1 = fun (xx1, nwds1, nwds2)
        tw1 = t1 * wkts(i)
        twi1 = abs (tw1)
        if (twi1 < epsilon1) then
          iz1 = iz1 + 1
        else
          iz1 = 0
        endif
      else
        t1 = mpreal (0.d0, nwds1)
        tw1 = mpreal (0.d0, nwds1)
      endif

      if (i > 0 .and. log2 .and. iz2 < izx) then
        t2 = fun (xx2, nwds1, nwds2)
        tw2 = t2 * wkts(i)
        twi2 = abs (tw2)
        if (twi2 < epsilon1) then
          iz2 = iz2 + 1
        else
          iz2 = 0
        endif
      else
        t2 = mpreal (0.d0, nwds1)
        tw2 = mpreal (0.d0, nwds1)
      endif

      tsum = tsum + tw1 + tw2
      twmx = max (twmx, abs (tw1), abs (tw2))
    endif
  enddo

!   Compute s1 = current integral approximation and err = error estimate.
!   Tsum is the sum of all tw1 and tw2 from the loop above.
!   Twmx is the largest absolute value of tw1 and tw2 from the loop above.
!   Twi1 and twi2 are the final nonzero values of abs(tw1) and abs(tw2).

  s1 =  mpreal (ax, nwds1) * h * tsum
  eps1 = twmx * epsilon1
  eps2 = abs (max (twi1, twi2) / s1)
  d1 = dplog10 (abs ((s1 - s2) / s1))
  d2 = dplog10 (abs ((s1 - s3) / s1))
  d3 = dplog10 (eps1) - 1.d0
  d4 = dplog10 (eps2) - 1.d0

  if (k <= 2) then
    err = mpreal (1.d0, nwds1)
  elseif (d1 == -999999.d0) then
    err = mpreal (0.d0, nwds1)
  else
    err = c10 ** nint (min (0.d0, max (d1 ** 2 / d2, 2.d0 * d1, d3, d4)))
  endif

!   Output current integral approximation and error estimate, to 60 digits.

  if (ndebug >= 2) then
    write (6, 2) k, nq1, nint (dplog10 (abs (err)))
2   format ('quadts: Iteration',i3,' of',i3,'; est error = 1.e',i6, &
      '; approx value =')
    call mpwrite (6, 80, 60, s1)
  endif

  if (k >= 3 .and. iz1 == 0 .and. iz2 == 0) then
    write (6, 3)
3   format ('quadts: Terms too large -- adjust neps2 in call to initqts.')
    goto 140
  endif

  if (k >= 3 .and. err < eps1) then
    write (6, 4) nint (dplog10 (abs (err)))
4   format ('quadts: Estimated error = 1.e',i6)
    goto 140
  endif

  if (k >= 3 .and. err < eps2) then
    write (6, 5) nint (dplog10 (abs (err)))
5   format ('quadts: Estimated error = 1.e',i6/&
    'Adjust nq1 and neps2 in initqts for greater accuracy.')
    goto 140
  endif
enddo

140 continue

quadts = s1
return
end function quadts

!   Exp-sinh routines:

subroutine initqes (nq1, nq2, nwds1, neps2, wkes, xkes)

!   David H Bailey   7 Jan 2023

!   This subroutine initializes the quadrature arrays xkes and wkes for quades.
!   If initqes outputs the message "Table space parameter is too small", adjust
!   nq2 in calling program. Also, if quades outputs the message "Terms too large",
!   adjust nq1 and neps2 as necessary in the call to initqes. The argument neps2
!   controls termination of the loop below, which ends when
!   wkes(k) * 10^(neps2) > 1.

!   This subroutine initializes the quadrature arrays xkes and wkes for quades.
!   The argument nq2 is the space allocated for wkes and xkes in the calling
!   program. By default it is set to 12 * 2^nq1. If initqes outputs the message
!   "Table space parameter is too small", adjust nq2. Also, if quades outputs
!   the message "Terms too large", adjust nq1 and neps2 as necessary in the call
!   to initqes. 

!   The wkes and xkes arrays are computed based on the transformation
!   t = exp (pi/2 * sinh (x)). See comments below.

!   Both initqes and quades are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   Input arguments; all are default (4-byte) integer:
!   nq1      Level parameter for data in wkes and xkes arrays.
!   nq2      Space parameter for wkes and xkes arrays.
!   nwds1    Primary precision level, in mantissa words.
!   neps2    Secondary epsilon level (eps2 = 10^neps2).

!   Output arguments:
!   wkes     Output array of weights; type (mp_real).
!   xkes     Output array of abscissas; type (mp_real).

use mpmodule
implicit none
integer, intent(in):: nq1, nq2, nwds1, neps2
type (mp_real), intent(out):: wkes(-1:nq2), xkes(-1:nq2)
integer, parameter:: iprint = 1024, ndebug = 2
integer k
type (mp_real) eps2, h, p2, t1, t2, u1, u2
  
write (6, 1)
1 format ('initqes: Exp-sinh quadrature initialization')

eps2 = mpreal (10.d0, nwds1) ** neps2
p2 = 0.5d0 * mppi (nwds1)
h = mpreal (0.5d0 ** nq1, nwds1)
wkes(-1) = mpreal (dble (nq1), nwds1)

do k = 0, nq2
  if (ndebug >= 2 .and. mod (k, iprint) == 0) write (6, '(2i8)') k, nq2
    t1 = mpreal (dble (k) * h, nwds1)

!   xkes(k) = exp (u1)
!   wkes(k) = exp (u1) * u2
!   where u1 = pi/2 * sinh (t1) and u2 = pi/2 * cosh (t1)

  t2 = exp (t1)
  u1 = 0.5d0 * p2 * (t2 - 1.d0 / t2)
  u2 = 0.5d0 * p2 * (t2 + 1.d0 / t2)
  xkes(k) = exp (u1)
  wkes(k) = xkes(k) * u2

  if (wkes(k) * eps2 > 1.d0) goto 100
enddo

write (6, 2) nq2
2 format ('initqes: Table space parameter is too small; value =',i8)
stop

100 continue

xkes(-1) = mpreal (dble (k), nwds1)
if (ndebug >= 2) then
  write (6, 3) k
3 format ('initqes: Table spaced used =',i8)
endif

return
end subroutine initqes

type (mp_real) function quades (fun, x1, nq1, nq2, nwds1, nwds2, neps1, wkes, xkes)

!   David H Bailey  7 Jan 2023

!   This routine computes the integral of the function fun on the interval
!   (x1, +infinity) with a target tolerance of 10^neps1. The quadrature level is
!   progressively increased (approximately doubling the work with each level)
!   until level nq1 has been performed or the target tolerance has been met.
!   Prior to calling quades, the wkes and xkes arrays must first be initialized
!   by calling initqes. If quades outputs the message "Terms too large", adjust
!   nq1 and neps2 as necessary in the call to initqes.

!   For some integrand functions, significantly more accurate quadrature results
!   can be achieved by computing x1 to higher precision (ndp2 digits or nwds2
!   words) in the calling program. x1 is used in quades for computation of the
!   scaled abscissa, which is in turn passed to the function definition for use
!   in initial subtractions or other numerically sensitive operations involving
!   the input argument. The function evaluation itself should be performed with
!   standard precision (ndp1 digits or nwds1 words) for faster run times.

!   Both initqes and quades are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   Input arguments; unless indicated otherwise all are default (4-byte) integer:
!   fun      Name of function to be integrated; type (mp_real).
!   x1       Lower limit of integration interval; type (mp_real).
!   nq1      Level parameter for data in wkes and xkes arrays.
!   nq2      Space parameter for wkes and xkes arrays.
!   nwds1    Primary precision level, in mantissa words.
!   nwds2    Secondary precision level, in mantissa words.
!   neps1    Primary epsilon level (eps1 = 10^neps1).
!   neps2    Secondary epsilon level (eps2 = 10^neps2).
!   wkes     Precomputed array of weights; type (mp_real).
!   xkes     Precomputed array of abscissas; type (mp_real).

!   Output:
!   quades   Quadrature result; type (mp_real);

use mpmodule
implicit none
type (mp_real), external:: fun
type (mp_real), intent(in):: x1
integer, intent(in):: nq1, nq2, nwds1, nwds2, neps1
type (mp_real), intent(in):: wkes(-1:nq2), xkes(-1:nq2)
real (mprknd), external:: dplog10
integer, parameter:: izx = 5, ndebug = 2
integer i, ip(0:100), iz1, iz2, k, k1, k2, n, nqq1
logical log1
real (mprknd) d1, d2, d3, d4
type (mp_real) c10, eps1, eps2, epsilon1, err, h, &
  tsum, s1, s2, s3, t1, t2, tw1, tw2, twi1, twi2, twmx, xx1, xx2

epsilon1 = mpreal (10.d0, nwds1) ** neps1
tsum = mpreal (0.d0, nwds1)
s1 = mpreal (0.d0, nwds1)
s2 = mpreal (0.d0, nwds1)
h = mpreal (1.d0, nwds1)
c10 = mpreal (10.d0, nwds1)

if (wkes(-1) < dble (nq1)) then
  write (6, 1) nq1
1 format ('quades: quadrature arrays have not been initialized; nq1 =',i6)
  goto 140
endif
nqq1 = dble (wkes(-1))
n = dble (xkes(-1))

do k = 0, nqq1
  ip(k) = 2 ** k
enddo

do k = 1, nq1
  h = 0.5d0 * h
  s3 = s2
  s2 = s1
  k1 = ip(nqq1-k)
  k2 = ip(nqq1-k+1)
  iz1 = 0
  iz2 = 0
  twmx = mpreal (0.d0, nwds1)

!   Evaluate function at level k in x, avoiding unnecessary computation.

  do i = 0, n, k1
    if (mod (i, k2) /= 0 .or. k == 1) then

!   These next few lines, which scale the abscissas, must be performed in
!   high precision (nwds2 words) to ensure full accuracy in the quadrature
!   results, even though the abscissas xkes(i) were computed in regular precision.

      xx1 = x1 + mpreal (xkes(i), nwds2)
      xx2 = x1 + 1.d0 / mpreal (xkes(i), nwds2)
      log1 = xx1 > x1
  
!   The remaining computations are performed in standard precision (nwds1 words).

      if (iz1 < izx) then
        t1 = fun (xx1, nwds1)
        tw1 = t1 * wkes(i)
        twi1 = abs (tw1)
        if (twi1 < epsilon1) then
          iz1 = iz1 + 1
        else
          iz1 = 0
        endif
      else
        t1 = mpreal (0.d0, nwds1)
        tw1 = mpreal (0.d0, nwds1)
      endif

      if (i > 0 .and. log1 .and. iz2 < izx) then
        t2 = fun (xx2, nwds1)
        tw2 = t2 * wkes(i) / xkes(i)**2
        twi2 = abs (tw2)
        if (twi2 < epsilon1) then
          iz2 = iz2 + 1
        else
          iz2 = 0
        endif
      else
        t2 = mpreal (0.d0, nwds1)
        tw2 = mpreal (0.d0, nwds1)
      endif

      tsum = tsum + tw1 + tw2
      twmx = max (twmx, abs (tw1), abs (tw2))
    endif
  enddo

!   Compute s1 = current integral approximation and err = error estimate.
!   Tsum is the sum of all tw1 and tw2 from the loop above.
!   Twmx is the largest absolute value of tw1 and tw2 from the loop above.
!   Twi1 and twi2 are the final nonzero values of abs(tw1) and abs(tw2).

  s1 =  h * tsum
  eps1 = twmx * epsilon1
  eps2 = abs (max (twi1, twi2) / s1)
  d1 = dplog10 (abs ((s1 - s2) / s1))
  d2 = dplog10 (abs ((s1 - s3) / s1))
  d3 = dplog10 (eps1) - 1.d0
  d4 = dplog10 (eps2) - 1.d0

  if (k <= 2) then
    err = mpreal (1.d0, nwds1)
  elseif (d1 == -999999.d0) then
    err = mpreal (0.d0, nwds1)
  else
    err = c10 ** nint (min (0.d0, max (d1 ** 2 / d2, 2.d0 * d1, d3, d4)))
  endif

!   Output current integral approximation and error estimate, to 60 digits.

  if (ndebug >= 2) then
    write (6, 2) k, nq1, nint (dplog10 (abs (err)))
2   format ('quades: Iteration',i3,' of',i3,'; est error = 1.e',i6, &
      '; approx value =')
    call mpwrite (6, 80, 60, s1)
  endif

  if (k >= 3 .and. iz1 == 0 .and. iz2 == 0) then
    write (6, 3)
3   format ('quades: Terms too large -- adjust neps2 in call to initqes.')
    goto 140
  endif

  if (k >= 3 .and. err < eps1) then
    write (6, 4) nint (dplog10 (abs (err)))
4   format ('quades: Estimated error = 1.e',i6)
    goto 140
  endif

  if (k >= 3 .and. err < eps2) then
    write (6, 5) nint (dplog10 (abs (err)))
5   format ('quades: Estimated error = 1.e',i6/&
    'Adjust nq1 and neps2 in initqes for greater accuracy.')
    goto 140
  endif
enddo

140 continue

quades = s1
return
end function quades

!   Sinh-sinh routines:

subroutine initqss (nq1, nq2, nwds1, neps2, wkss, xkss)

!   David H Bailey   7 Jan 2023

!   This subroutine initializes the quadrature arrays xkss and wkss for quadss.
!   If initqss outputs the message "Table space parameter is too small", adjust
!   nq2 in calling program. Also, if quadss outputs the message "Terms too large",
!   adjust nq1 and neps2 as necessary in the call to initqss. The argument neps2
!   controls termination of the loop below, which ends when
!   wkss(k) * 10^(neps2) > 1.

!   The wkss and xkss arrays are computed based on the transformation
!   t = sinh (pi/2 * sinh (x)).  See comments below.

!   Both initqss and quadss are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   Input arguments; all are default (4-byte) integer:
!   nq1      Level parameter for data in wkss and xkss arrays.
!   nq2      Space parameter for wkss and xkss arrays.
!   nwds1    Primary precision level, in mantissa words.
!   neps2    Secondary epsilon level (eps2 = 10^neps2).

!   Output arguments:
!   wkss     Output array of weights; type (mp_real).
!   xkss     Output array of abscissas; type (mp_real).

use mpmodule
implicit none
integer, intent(in):: nq1, nq2, nwds1, neps2
type (mp_real), intent(out):: wkss(-1:nq2), xkss(-1:nq2)
integer, parameter:: iprint = 1024, ndebug = 2
integer k
type (mp_real) eps2, h, p2, t1, t2, t3, u1, u2

write (6, 1)
1 format ('initqss: Sinh-sinh quadrature initialization')

eps2 = mpreal (10.d0, nwds1) ** neps2
p2 = 0.5d0 * mppi (nwds1)
h = mpreal (0.5d0 ** nq1, nwds1)
wkss(-1) = mpreal (dble (nq1), nwds1)

do k = 0, nq2
  if (ndebug >= 2 .and. mod (k, iprint) == 0) write (6, '(2i8)') k, nq2
    t1 = mpreal (dble (k) * h, nwds1)

!   xkss(k) = sinh (u1)
!   wkss(k) = cosh (u1) * u2
!   where u1 = pi/2 * sinh (t1) and u2 = pi/2 * cosh (t1)

  t2 = exp (t1)
  u1 = 0.5d0 * p2 * (t2 - 1.d0 / t2)
  u2 = 0.5d0 * p2 * (t2 + 1.d0 / t2)
  t3 = exp (u1)
  xkss(k) = 0.5d0 * (t3 - 1.d0 / t3)
  wkss(k) = 0.5d0 * (t3 + 1.d0 / t3) * u2

  if (wkss(k) * eps2 > 1.d0) goto 100
enddo

write (6, 2) nq2
2 format ('initqss: Table space parameter is too small; value =',i8)
stop

100 continue

xkss(-1) = mpreal (dble (k), nwds1)
if (ndebug >= 2) then
  write (6, 3) k
3 format ('initqss: Table spaced used =',i8)
endif

return
end subroutine initqss

type (mp_real) function quadss (fun, nq1, nq2, nwds1, neps1, wkss, xkss)

!   David H Bailey  7 Jan 2023

!   This routine computes the integral of the function fun on the interval
!   (-inf, +inf) with a target tolerance of 10^neps1. The quadrature level is
!   progressively increased (approximately doubling the work with each level)
!   until level nq1 has been performed or the target tolerance has been met.
!   Prior to calling quadss, the wkss and xkss arrays must first be initialized
!   by calling initqss. If quadss outputs the message "Terms too large", adjust
!   nq1 and neps2 as necessary in the call to initqss.

!   Both initqss and quadss are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   Input arguments; unless indicated otherwise all are default (4-byte) integer:
!   fun      Name of function to be integrated; type (mp_real).
!   nq1      Level parameter for data in wkss and xkss arrays.
!   nq2      Space parameter for wkss and xkss arrays.
!   nwds1    Primary precision level, in mantissa words.
!   neps1    Primary epsilon level (eps1 = 10^neps1).
!   wkss     Precomputed array of weights; type (mp_real).
!   xkss     Precomputed array of abscissas; type (mp_real).

!   Output argument:
!   quadss   Quadrature result; type (mp_real).

use mpmodule
implicit none
type (mp_real), external:: fun
integer, intent(in):: nq1, nq2, nwds1, neps1
type (mp_real), intent(in):: wkss(-1:nq2), xkss(-1:nq2)
integer, parameter:: izx = 5, ndebug = 2
real (mprknd), external:: dplog10
integer i, ip(0:100), iz1, iz2, k, k1, k2, n, nqq1
real (mprknd) d1, d2, d3, d4
type (mp_real) c10, eps1, eps2, epsilon1, err, h, &
  tsum, s1, s2, s3, t1, t2, tw1, tw2, twi1, twi2, twmx

epsilon1 = mpreal (10.d0, nwds1) ** neps1
tsum = mpreal (0.d0, nwds1)
s1 = mpreal (0.d0, nwds1)
s2 = mpreal (0.d0, nwds1)
h = mpreal (1.d0, nwds1)
c10 = mpreal (10.d0, nwds1)

if (wkss(-1) < dble (nq1)) then
  write (6, 1) nq1
1 format ('quadss: quadrature arrays have not been initialized; nq1 =',i6)
  goto 140
endif
nqq1 = dble (wkss(-1))
n = dble (xkss(-1))

do k = 0, nqq1
  ip(k) = 2 ** k
enddo

do k = 1, nq1
  h = 0.5d0 * h
  s3 = s2
  s2 = s1
  k1 = ip(nqq1-k)
  k2 = ip(nqq1-k+1)
  iz1 = 0
  iz2 = 0
  twmx = mpreal (0.d0, nwds1)

!   Evaluate function at level k in x, avoiding unnecessary computation.

  do i = 0, n, k1
    if (mod (i, k2) /= 0 .or. k == 1) then
      if (iz1 < izx) then
        t1 = fun (xkss(i), nwds1)
        tw1 = t1 * wkss(i)
        twi1 = abs (tw1)
        if (twi1 < epsilon1) then
          iz1 = iz1 + 1
        else
          iz1 = 0
        endif
      else
        t1 = mpreal (0.d0, nwds1)
        tw1 = mpreal (0.d0, nwds1)
      endif

      if (i > 0 .and. iz2 < izx) then
        t2 = fun (-xkss(i), nwds1)
        tw2 = t2 * wkss(i)
        twi2 = abs (tw2)
        if (twi2 < epsilon1) then
          iz2 = iz2 + 1
        else
          iz2 = 0
        endif
      else
        t2 = mpreal (0.d0, nwds1)
        tw2 = mpreal (0.d0, nwds1)
      endif

      tsum = tsum + tw1 + tw2
      twmx = max (twmx, abs (tw1), abs (tw2))
    endif
  enddo

!   Compute s1 = current integral approximation and err = error estimate.
!   Tsum is the sum of all tw1 and tw2 from the loop above.
!   Twmx is the largest absolute value of tw1 and tw2 from the loop above.
!   Twi1 and twi2 are the final nonzero values of abs(tw1) and abs(tw2).

  s1 =  h * tsum
  eps1 = twmx * epsilon1
  eps2 = abs (max (twi1, twi2) / s1)
  d1 = dplog10 (abs ((s1 - s2) / s1))
  d2 = dplog10 (abs ((s1 - s3) / s1))
  d3 = dplog10 (eps1) - 1.d0
  d4 = dplog10 (eps2) - 1.d0

  if (k <= 2) then
    err = mpreal (1.d0, nwds1)
  elseif (d1 == -999999.d0) then
    err = mpreal (0.d0, nwds1)
  else
    err = c10 ** nint (min (0.d0, max (d1 ** 2 / d2, 2.d0 * d1, d3, d4)))
  endif

!   Output current integral approximation and error estimate, to 60 digits.

  if (ndebug >= 2) then
    write (6, 2) k, nq1, nint (dplog10 (abs (err)))
2   format ('quadss: Iteration',i3,' of',i3,'; est error = 1.e',i6, &
      '; approx value =')
    call mpwrite (6, 80, 60, s1)
  endif

  if (k >= 3 .and. iz1 == 0 .and. iz2 == 0) then
    write (6, 3)
3   format ('quadss: Terms too large -- adjust neps2 in call to initqss.')
    goto 140
  endif

  if (k >= 3 .and. err < eps1) then
    write (6, 4) nint (dplog10 (abs (err)))
4   format ('quadss: Estimated error = 1.e',i6)
    goto 140
  endif

  if (k >= 3 .and. err < eps2) then
    write (6, 5) nint (dplog10 (abs (err)))
5   format ('quadss: Estimated error = 1.e',i6/&
    'Adjust nq1 and neps2 in initqss for greater accuracy.')
    goto 140
  endif
enddo

140 continue

quadss = s1
return
end function quadss

real (mprknd) function dplog10 (a)

!   For input MP value a, this routine returns a DP approximation to log10 (a).

use mpmodule
implicit none
type (mp_real), intent(in):: a
integer ia
real (mprknd) da

call mpdecmd (a, da, ia)
if (da == 0.d0) then
  dplog10 = -999999.d0
else
  dplog10 = log10 (abs (da)) + ia
endif

100 continue
return
end function dplog10
