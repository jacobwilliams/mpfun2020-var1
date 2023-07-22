!*****************************************************************************

!  program testmpfun

!  Revision date: 28 Mar 2023

!  AUTHOR:
!   David H. Bailey
!   Lawrence Berkeley National Lab (retired) and University of California, Davis
!   Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!   All software in this package (c) 2023 David H. Bailey.
!   By downloading or using this software you agree to the copyright, disclaimer
!   and license agreement in the accompanying file DISCLAIMER.txt.

!  DESCRIPTION:
!   This briefly tests most individual MPFUN2020 operations and functions
!   (including mixed mode arithmetic, comparison operations, transcendental
!   functions and special functions), by comparing each result with benchmark
!   results in the file testmpfun.ref.txt, which must be present in the same
!   directory. This is not an exhaustive test of all possible scenarios, but it
!   often detects bugs and compiler issues.

program testmpfun
use mpmodule
implicit none
integer, parameter:: nfile = 11, ndp = 500, n1 = ndp + 30, n2 = ndp, np = 2, nq = 3, &
  nbe = 2.d0 * ndp, neps = 5 - ndp, nrr = 10, nwds = ndp / mpdpw + 2
real (mprknd), external:: second
integer i, i1, k
logical l1, l2, l3, l4, l5
character(ndp+30) chrx
character(1) chr1(ndp+30)
real (mprknd) d1, d2, e1, e2, tm0, tm1
complex (mprknd) ec1, ec2
type (mp_real) aa(np), bb(nq), eps, err, errmx, one, be(nbe), rr(nrr), &
  t1, t2, t3, t4, zero
type (mp_complex) z1, z2, z3

!  End of declaration

tm0 = second ()
open (nfile, file = 'testmpfun.ref.txt')
rewind nfile
eps = mpreal (10.d0, nwds) ** neps
one = mpreal (1.d0, nwds)
zero = mpreal (0.d0, nwds)
errmx = zero
aa(1) = mpreal (0.75d0, nwds)
aa(2) = mpreal (1.25d0, nwds)
bb(1) = mpreal (1.d0, nwds)
bb(2) = mpreal (1.5d0, nwds)
bb(3) = mpreal (2.d0, nwds)
call mpberne (nbe, be, nwds)
call polylog_ini (-nrr, rr, 2*nwds)

write (6, '(a)') 'MPUN2020 quick check of operations and functions'

!   Define a few sample data values.

chrx = ' '

t1 = mppi (nwds)
t2 = - mplog2 (nwds)

e1 = 3141.d0 / 8192.d0
e2 = 6931.d0 / 8192.d0
z1 = mpcmplx (0.5d0 * mppi (nwds), exp (mpreal (0.5d0, nwds)), nwds)
z2 = mpcmplx (- gamma (mpreal (0.5d0, nwds)), cos (mpreal (1.d0, nwds)), nwds)
ec1 = cmplx (e1, e2, mprknd)
ec2 = cmplx (-e2, e1, mprknd)
i1 = 5

write (6, '(/a/)') 'Test data:'
write (6, '(a)') 't1 = pi:'
call mpwrite (6, ndp + 20, ndp, t1)
call checkmp (nwds, nfile, 5, t1, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 't2 = -log(2):'
call mpwrite (6, ndp + 20, ndp, t2)
call checkmp (nwds, nfile, 1, t2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'z1 = (0.5*pi, exp(0.5)):'
call mpwrite (6, ndp + 20, ndp, z1)
call checkmpc (nwds, nfile, 1, z1, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'z2 = (-Gamma(0.5), Cos(1)):'
call mpwrite (6, ndp + 20, ndp, z2)
call checkmpc (nwds, nfile, 1, z2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'e1 = 3141/8192:'
write (6, '(1p,d25.15)') e1
write (6, '(a)') 'e2 = 6931/8192:'
write (6, '(1p,d25.15)') e2
write (6, '(a)')  'ec1 = (3141/8192, 6931/8192)'
write (6, '(1p,2d25.15)') ec1
write (6, '(a)') 'ec2 = (-6931/8192, 3141/8192):'
write (6, '(1p,2d25.15)') ec2
write (6, '(a,i4)') 'i1: ', i1

write (6, '(/a/)') 'Real data operations:'

write (6, '(a)') 'addition: t1+t2 ='
t3 = t1 + t2
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 13, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'addition: t1+e2 ='
t3 = t1 + e2
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'addition: e1+t2 ='
t3 = e1 + t2
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'subtraction: t1-t2 ='
t3 = t1 - t2
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'subtraction: t1-e2 ='
t3 = t1 - e2
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'subtraction: e1-t2 ='
t3 = e1 - t2
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'multiplication: t1*t2 ='
t3 = t1 * t2
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'multiplication: t1*e2 ='
t3 = t1 * e2
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'multiplication: e1*t2 ='
t3 = e1 * t2
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'division: t1/t2 ='
t3 = t1 / t2
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'division: t1/e2 ='
t3 = t1 / e2
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'division: e1/t2 ='
t3 = e1 / t2
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'exponentiation: t1**i1 ='
t3 = t1 ** i1
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'exponentiation: t1**t2 ='
t3 = t1 ** t2
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'equal test: t1 == t2, e1 == t2, t1 == e2'
l1 = t1 == t2; l2 = e1 == t2; l3 = t1 == e2
write (6, '(6l4)') l1, l2, l3
if (l1 .eqv. .false. .and. l2 .eqv. .false. .and. l3 .eqv. .false.) then
  err = zero
else
  err = one
endif
errmx = max (err, errmx)

write (6, '(a)') 'not-equal test: t1 /= t2, e1 /= t2, t1 =/ e2'
l1 = t1 /= t2; l2 = e1 /= t2; l3 = t1 /= e2
write (6, '(6l4)') l1, l2, l3
if (l1 .eqv. .true. .and. l2 .eqv. .true. .and. l3 .eqv. .true.) then
  err = zero
else
  err = one
endif
errmx = max (err, errmx)

write (6, '(a)') 'less-than-or-equal test: t1 <= t2, e1 <= t2, t1 <= e2'
l1 = t1 <= t2; l2 = e1 <= t2; l3 = t1 <= e2
write (6, '(6l4)') l1, l2, l3
if (l1 .eqv. .false. .and. l2 .eqv. .false. .and. l3 .eqv. .false.) then
  err = zero
else
  err = one
endif
errmx = max (err, errmx)

write (6, '(a)') 'greater-than-or-equal test: t1 >= t2, e1 >= t2, t1 >= e2'
l1 = t1 >= t2; l2 = e1 >= t2; l3 = t1 >= e2
write (6, '(6l4)') l1,l2, l3
if (l1 .eqv. .true. .and. l2 .eqv. .true. .and. l3 .eqv. .true.) then
  err = zero
else
  err = one
endif
errmx = max (err, errmx)

write (6, '(a)') 'less-than test: t1 < t2, e1 < t2, t1 < e2'
l1 = t1 < t2; l2 = e1 < t2; l3 = t1 < e2
write (6, '(6l4)') l1, l2, l3
if (l1 .eqv. .false. .and. l2 .eqv. .false. .and. l3 .eqv. .false.) then
  err = zero
else
  err = one
endif
errmx = max (err, errmx)

write (6, '(a)') 'greater-than test: t1 > t2, e1 > t2, t1 > e2'
l1 = t1 > t2; l2 = e1 > t2; l3 = t1 > e2
write (6, '(6l4)') l1, l2, l3
if (l1 .eqv. .true. .and. l2 .eqv. .true. .and. l3 .eqv. .true.) then
  err = zero
else
  err = one
endif
errmx = max (err, errmx)

write (6, '(a)') 'abs(t2) ='
t3 = abs (t2)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 13, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'acos(t2) ='
t3 = acos (t2)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'acosh(t1) ='
t3 = acosh (t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'agm(t1,abs(t2)) ='
t3 = agm (t1, abs (t2))
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'aint(t1) ='
t3 = aint (t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)

write (6, '(a)') 'anint(t1) ='
t3 = anint (t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)

write (6, '(a)') 'asin(t2) ='
t3 = asin (t2)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'asinh(t1) ='
t3 = asinh (t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'atan(t1) ='
t3 = atan (t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'atan2(t1,t2) ='
t3 = atan2 (t1,t2)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_i(t1,-t2) ='
t3 = bessel_i (t1, -t2)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_i(t1,200*t1) ='
t3 = bessel_i (t1, 200.d0*t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_i(t1,300*t1) ='
t3 = bessel_i (t1, 300.d0*t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_in(3,t1) ='
t3 = bessel_in (3, t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_in(3,200*t1) ='
t3 = bessel_in (3, 200.d0*t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_in(3,300*t1) ='
t3 = bessel_in (3, 300.d0*t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_j(t1,-t2) ='
t3 = bessel_j (t1, -t2)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_j(t1,200*t1) ='
t3 = bessel_j (t1, 200.d0*t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_j(t1,300*t1) ='
t3 = bessel_j (t1, 300.d0*t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_jn(3,t1) ='
t3 = bessel_jn (3, t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_jn(3,200*t1) ='
t3 = bessel_jn (3, 200.d0*t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_jn(3,300*t1) ='
t3 = bessel_jn (3, 300.d0*t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_k(t1,-t2) ='
t3 = bessel_k (t1, -t2)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_k(t1,200*t1) ='
t3 = bessel_k (t1, 200.d0*t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_k(t1,300*t1) ='
t3 = bessel_k (t1, 300.d0*t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_kn(3,t1) ='
t3 = bessel_kn (3, t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_kn(3,200*t1) ='
t3 = bessel_kn (3, 200.d0*t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_kn(3,300*t1) ='
t3 = bessel_kn (3, 300.d0*t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_y(t1,-t2) ='
t3 = bessel_y (t1, -t2)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_y(t1,200*t1) ='
t3 = bessel_y (t1, 200.d0*t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_y(t1,300*t1) ='
t3 = bessel_y (t1, 300.d0*t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_yn(3,t1) ='
t3 = bessel_yn (3, t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_yn(3,200*t1) ='
t3 = bessel_yn (3, 200.d0*t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_yn(3,300*t1) ='
t3 = bessel_yn (3, 300.d0*t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'cos(t2) ='
t3 = cos (t2)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'cosh(t1) ='
t3 = cosh (t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'digamma_be(nbe,be,t1) ='
t3 = digamma_be (nbe, be, t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'erf(t1) ='
t3 = erf (t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'erfc(t1) ='
t3 = erfc (t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'exp(t1) ='
t3 = exp (t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'expint(t1) ='
t3 = expint (t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'gamma(t1) ='
t3 = gamma (t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'hurwitz_zetan(3,1/t1) ='
t3 = hurwitz_zetan (3,1.d0/t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'hurwitz_zetan_be(nbe,be,5,t1) ='
t3 = hurwitz_zetan_be (nbe, be, 5, t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'hypergeom_pfq(np,nq,aa,bb,t1) ='
t3 = hypergeom_pfq (np, nq, aa, bb, t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'incgamma(t1,-t2) ='
t3 = incgamma (t1, -t2)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'incgamma(-t1,-t2) ='
t3 = incgamma (-t1, -t2)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'incgamma(t1,-300*t2) ='
t3 = incgamma(t1,-300.d0*t2) 
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'incgamma(t1,-2000*t2) ='
t3 = incgamma(t1,-2000.d0*t2) 
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'log(t1) ='
t3 = log (t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'log10(t1) ='
t3 = log10 (t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'max(t1,t2) ='
t3 = max (t1,t2)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'min(t1,t2) ='
t3 = min (t1,t2)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'mod(t1,t2) ='
t3 = mod(t1,t2)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'mpcssh(t1) ='
call mpcssh (t1, t3, t4)
call mpwrite (6, ndp + 20, ndp, t3, t4)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)
call checkmp (nwds, nfile, 0, t4, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'mpcssn(t2) ='
call mpcssn (t2, t3, t4)
call mpwrite (6, ndp + 20, ndp, t3, t4)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)
call checkmp (nwds, nfile, 0, t4, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'mpnrt(t1,i1) ='
t3 = mpnrt (t1,i1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'mpprod (t1,e2) ='
t3 = mpprod (t1, e2)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'mpquot (t1,e2) ='
t3 = mpquot (t1,e2)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'mpeform(t1,n1,n2,chr1) ='
call mpeform (t1, n1, n2, chr1)
write (6, '(80a1,"\")') (chr1(i), i = 1, n1)
do i = 1, n1; chrx(i:i) = chr1(i); enddo; t4 = mpreal (chrx(1:n1), nwds)
call checkmp (nwds, nfile, 1, t4, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'mpfform(t1,n1,n2,chr1) ='
call mpfform (t1, n1, n2, chr1)
write (6, '(80a1,"\")') (chr1(i), i = 1, n1)
do i = 1, n1; chrx(i:i) = chr1(i); enddo; t4 = mpreal (chrx(1:n1), nwds)
call checkmp (nwds, nfile, 1, t4, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'mpreal (chrx, nwds) ='
t3 = mpreal (chrx, nwds)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'mpreal (chr1, n1, nwds) ='
t3 = mpreal (chr1, n1, nwds)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)

write (6, '(a)') 'polygamma(3,1/t1) ='
t3 = polygamma (3,1.d0/t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'polygamma_be(nbe,be,5,t1) ='
t3 = polygamma_be (nbe,be,5,t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'polylog_neg (-10, rr, -t1) ='
t3 = polylog_neg (-10, rr, -t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'polylog_pos (10, 1/t1) ='
t3 = polylog_pos (10, 1.d0/t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'sign(t1,t2) ='
t3 = sign (t1,t2)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'sin(t2) ='
t3 = sin (t2)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'sinh(t1) ='
t3 = sinh (t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'sqrt(t1) ='
t3 = sqrt (t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'struve_hn(10,t1) ='
t3 = struve_hn (10, t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'tan(t1) ='
t3 = tan (t2)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'tanh(t1) ='
t3 = tanh (t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'zeta(t1) ='
t3 = zeta (t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'zeta_be(nbe,be,t1) ='
t3 = zeta_be (nbe, be, t1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'zeta_int(10,nwds) ='
t3 = zeta_int (10,nwds)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(/a/)') 'Complex data operations:'

write (6, '(a)') 'addition: z1+z2 ='
z3 = z1 + z2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 4, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'addition: z1+e2 ='
z3 = z1 + e2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'addition: e1+z2 ='
z3 = e1 + z2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'addition: z1+ec2 ='
z3 = z1 + ec2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'addition: ec1+z2 ='
z3 = ec1 + z2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'addition: z1+t2 ='
z3 = z1 + t2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'addition: t1+z2 ='
z3 = t1 + z2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'subtraction: z1-z2 ='
z3 = z1 - z2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'subtraction: z1-e2 ='
z3 = z1 - e2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'subtraction: e1-z2 ='
z3 = e1 - z2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'subtraction: z1-ec2 ='
z3 = z1 - ec2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'subtraction: ec1-z2 ='
z3 = ec1 - z2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'subtraction: z1-t2 ='
z3 = z1 - t2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'subtraction: t1-z2 ='
z3 = t1 - z2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'multiplication: z1*z2 ='
z3 = z1 * z2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'multiplication: z1*e2 ='
z3 = z1 * e2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'multiplication: e1*z2 ='
z3 = e1 * z2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'multiplication: z1*ec2 ='
z3 = z1 * ec2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'multiplication: ec1*z2 ='
z3 = ec1 * z2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'multiplication: z1*t2 ='
z3 = z1 * t2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'multiplication: t1*z2 ='
z3 = t1 * z2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'multiplication: z1*e2 ='
z3 = z1 * e2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'multiplication: e1*z2 ='
z3 = e1 * z2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'division: z1/z2 ='
z3 = z1 / z2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'division: z1/e2 ='
z3 = z1 / e2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'division: e1/z2 ='
z3 = e1 / z2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'division: z1/ec2 ='
z3 = z1 / ec2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'division: ec1/z2 ='
z3 = ec1 / z2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'division: z1/t2 ='
z3 = z1 / t2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'division: t1/z2 ='
z3 = t1 / z2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'exponentiation: z1**i1 ='
z3 = z1 ** i1
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'exponentiation: z1**z2 ='
z3 = z1 ** z2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'exponentiation: t1**z2 ='
z3 = t1 ** z2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'exponentiation: z1**t2 ='
z3 = z1 ** t2
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'equal test: z1 == z2, e1 == z2, z1 == e2, ec1 == z2, z1 == ec2'
l1 = z1 == z2; l2 = e1 == z2; l3 = z1 == e2; l4 = ec1 == z2; l5 = z1 == ec2
write (6, '(10l4)') l1, l2, l3, l4, l5
if (l1 .eqv. .false. .and. l2 .eqv. .false. .and. l3 .eqv. .false. &
  .and. l4 .eqv. .false. .and. l5 .eqv. .false.) then
  err = zero
else
  err = one
endif
errmx = max (err, errmx)

write (6, '(a)') 'not-equal test: z1 /= z2, e1 /= z2, z1 /= e2, ec1 /= z2, z1 /= ec2'
l1 = z1 /= z2; l2 = e1 /= z2; l3 = z1 /= e2; l4 = ec1 /= z2; l5 = z1 /= ec2
write (6, '(10l4)') l1, l2, l3, l4, l5
if (l1 .eqv. .true. .and. l2 .eqv. .true. .and. l3 .eqv. .true. &
  .and. l4 .eqv. .true. .and. l5 .eqv. .true.) then
  err = zero
else
  err = one
endif
errmx = max (err, errmx)

write (6, '(a)') 'abs(z2) ='
t3 = abs (z2)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 5, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'aimag(z1) ='
t3 = aimag (z1)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'conjg(z1) ='
z3 = conjg (z1)
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'cos(z2) ='
z3 = cos (z2)
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'exp(z1) ='
z3 = exp (z1)
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'log(z1) ='
z3 = log (z1)
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'mpcmplx(t1,t2) ='
z3 = mpcmplx (t1, t2, nwds)
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'mpreal(z1) ='
t3 = mpreal (z1, nwds)
call mpwrite (6, ndp + 20, ndp, t3)
call checkmp (nwds, nfile, 1, t3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'sin(z2) ='
z3 = sin (z2)
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'sqrt(z1) ='
z3 = sqrt (z1)
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'sqrt(z2) ='
z3 = sqrt (z2)
call mpwrite (6, ndp + 20, ndp, z3)
call checkmpc (nwds, nfile, 1, z3, eps, err)
errmx = max (err, errmx)

tm1 = second ()
call mpdecmd (errmx, d1, i1)
write (6, 9) tm1 - tm0, d1, i1
9 format (/'CPU time =',f8.2/'Max relative error =',f9.6,'e',i6)

if (abs (errmx) < eps) then
  write (6, '(a)') 'ALL TESTS PASSED'
else
  write (6, '(a)') 'ONE OR MORE TESTS FAILED'
endif

stop
end program testmpfun

subroutine checkmp (nwds, nfile, i1, t1, eps, err)
use mpmodule
implicit none
type (mp_real) t1, t2, eps, err
integer i, nwds, nfile, i1
character(64) c1

do i = 1, i1
  read (nfile, '(a)') c1
enddo

call mpread (nfile, t2, nwds)
err = abs ((t1 - t2) / t2)

if (abs (err) > eps) then
  write (6, '(a)') 'ERROR:'
  call mpwrite (6, 60, 40, abs (err))
endif

return
end subroutine checkmp

subroutine checkmpc (nwds, nfile, i1, z1, eps, err)
use mpmodule
implicit none
type (mp_complex) z1, z2
type (mp_real) eps, err
integer i, nwds, nfile, i1
character(64) c1

do i = 1, i1
  read (nfile, '(a)') c1
enddo

call mpread (nfile, z2, nwds)
err = abs ((z1 - z2) / z2)

if (abs (err) > eps) then
  write (6, '(a)') 'ERROR:'
  call mpwrite (6, 60, 40, abs (err))
endif

return
end subroutine checkmpc
