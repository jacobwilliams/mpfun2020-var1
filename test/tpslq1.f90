!*****************************************************************************

!  program tpslq1

!  Revision date:  12 Jan 2023

!  AUTHOR:
!   David H. Bailey
!   Lawrence Berkeley National Lab (retired) and University of California, Davis
!   Email: dhbailey@lbl.gov

! COPYRIGHT AND DISCLAIMER:
!   All software in this package (c) 2023 David H. Bailey.
!   By downloading or using this software you agree to the copyright, disclaimer
!   and license agreement in the accompanying file DISCLAIMER.txt.

! DESCRIPTION:
!   This program demonstrates pslq1, which performs the one-level standard
!   PSLQ algorithm on an input vector.  A variety of sample input vectors can
!   be generated as inputs to pslq1, as given in the parameters below. The
!   pslq1 routine is suitable for relations up to degree 25 or so; above this
!   level pslqm2 or pslqm3 should be used for significantly better performance.
!   For additional details, see:

!   David H. Bailey and David J. Broadhurst, "Parallel integer relation
!   detection: Techniques and applications," Mathematics of Computation,
!   vol. 70, no. 236 (Oct 2000), pg. 1719-1736, preprint available at
!   http://www.davidhbailey.com/dhbpapers/ppslq.pdf.

!   The pslq1 routine is 100% THREAD SAFE -- all requisite parameters and
!   arrays are passed through subroutine arguments.

!   PSLQ1 parameters set below; all are default (4-byte) integer:
!     idb   Debug level (0-4); default = 2.
!     n     Integer relation vector length; default = 31.
!           For minimal polynomial problems, n = 1 + polynomial degree.
!     ndp   Full precision level in digits; default = 250.
!           ***Must be <= mpipl in module MPFUNF.
!           ***Must be >= ndpm + precision level required to find relation.
!     ndr   Log10 of the minimum dynamic range in y at detection; default = 20.
!           A detected relation is not deemed reliable unless this is exceeded.
!     nep   Log10 of full precision epsilon for detections; default = 20 - ndp.
!           ***Must not be smaller than the accuracy of input x vector. In other
!           words, if data is accurate say to within 10^(-200), then nep > -200.
!     nrb   Log10 of maximum size (Euclidean norm) of acceptable relation;
!           default = 200. Run will be aborted if this is exceeded.
!     nwds  Full precision level in words; default = int (ndp/mpdpw + 2).

!   TPSLQ1 parameters set below; all are default (4-byte) integer:
!     kq    0: for the algebraic case [1, alpha, alpha^2, ... alpha^(n-1)], where
!             alpha = 3^(1/kr) - 2^(1/ks); this is the default.
!           1: for testing algebraic relations of a number read from a file.
!           2: for testing additive relations of numbers read from a file.
!           3: for testing multiplicative relations of numbers read from a file.
!           4: for custom input.
!     kr    Degree of root of 3 when kq = 0; default = 5.
!     ks    Degree of root of 2 when kq = 0; default = 6.
!     lcx   Size of character*1 array chr1; default = 64.
!           ***Must be a multiple of 64.
!           ***Must be > size in digits of largest result coefficient.

program tpslq1
use mpmodule
implicit none

!  PSLQ1 parameters:

integer, parameter:: kq = 0, kr = 5, ks = 6, lcx = 64

!  TPSLQ1 parameters:

integer, parameter:: idb = 2, n = kr * ks + 1, ndp = 250, ndr = 20, &
  nep = 20 - ndp, nrb = 200, nwds = int (ndp / mpdpw + 2)

integer i, i1, iq, j, j1
real (mprknd) second, tm0, tm1
integer lnm(n)
character(1) chr1(lcx)
character(64) form4, nam(n), namx
type (mp_real) al, r(n), x(n)
external second

!   Check to see if default precision is high enough.

if (ndp > mpipl) then
  write (6, '("Increase default precision in module MPFUNF.")')
  stop
endif

write (6, 1) n, kq, kr, ks, ndr, nrb, ndp, nep
1 format ('PSLQ1 Test Program'/ &
  'n =',i4,3x,'kq =',i2,3x,'kr =',i2,3x,'ks =',i2,3x,'ndr =',i5,3x,'nrb =',i5/ &
  'Full precision level ndp =',i6,' digits'/ &
  'Full precision epsilon level nep = ',i6)

if (kq == 1 .or. kq == 2 .or. kq == 3) then
  open (11, file = 'pslq.inp')
  rewind 11
endif

if (kq == 0) then

!   This code generates al = 3^(1/kr) - 2^(1/ks).  al is algebraic of degree
!   kr * ks.  Set n = kr * ks + 1 to recover the polynomial satisfied by al.

  al = mpnrt (mpreald (3.d0, nwds), kr) - mpnrt (mpreald (2.d0, nwds), ks)
elseif (kq == 1) then

!   Read an algebraic constant from a file.

  call mpread (11, al, nwds)
elseif (kq == 2) then

!   Read constants from a file for additive test.

  do i = 1, n
    call mpread (11, al, nwds)
    x(i) = al
    write (namx, '(''con'',i3.3)') i
    nam(i) = namx(1:6)
    lnm(i) = 6
  enddo
elseif (kq == 3) then

!   Read constants from a file for multiplicative test.

  do i = 1, n
    call mpread (11, al, nwds)
    x(i) = log (al)
    write (namx, '(''log(con'',i3.3,'')'')') i
    nam(i) = namx(1:11)
    lnm(i) = 11
  enddo
elseif (kq == 4) then

!   Produce X vector by a custom scheme.

endif

!   If kq is 0 or 1, generate x = [1, al, al^2, ..., al^(n-1)].

if (kq == 0 .or. kq == 1) then
  x(1) = mpreald (1.d0, nwds)
  nam(1) = '1'
  lnm(1) = 1

  do i = 2, n
    x(i) = al * x(i-1)
    write (namx, '(''al^'',i3)') i - 1
    nam(i) = namx(1:6)
    lnm(i) = 6
  enddo
endif

!   Perform relation search.

tm0 = second ()
call pslq1 (idb, n, nwds, ndr, nrb, nep, x, iq, r)
tm1 = second ()

!   Output relation, if one was found.

if (iq == 1) then
  write (6, 3)
3 format (/'Recovered relation: 0 =')
  form4 = '(64a1'
  i1 = 5

  do i = 65, lcx, 64
    form4(i1+1:i1+9) =  '/64a1'
    i1 = i1 + 5
  enddo

  form4(i1+1:i1+9) = '," * ",a)'
  i1 = i1 + 9
  
  do i = 1, n
    if (r(i) .ne. 0.d0) then
      call mpfform (r(i), lcx, 0, chr1)

      do j = 1, lcx
        if (chr1(j) /= ' ') goto 120
      enddo

120   continue

      j1 = j      
      if (chr1(j1) /= '-') then
        if (j1 > 1) then
          j1 = j1 - 1
          chr1(j1) = '+'
        else
          do j = 1, lcx
            chr1(j) = '*'
          enddo
        endif 
      endif
    
      write (6, form4) (chr1(j), j = 1, lcx), nam(i)(1:lnm(i))
    endif
  enddo
endif

write (6, 6) tm1 - tm0
6 format ('CPU Time =',f12.4)

if (iq == 1) then
  write (6, '(a)') 'TEST PASSED'
else
  write (6, '(a)') 'TEST FAILED'
endif

stop
end program tpslq1

!------------------------------

!   The following code performs the one-level, standard PSLQ algorithm.
!   David H. Bailey   12 Jan 2023

subroutine pslq1 (idb, n, nwds, ndr, nrb, nep, x, iq, r)

!   Arguments are as follows; int means default (4-byte) integer:
!     Name  Type    Description
!     idb   int     Debug flag (0-3); increasing idb produces more output.
!     n     int     Length of input vector x and output relation r.
!     nwds  int     Full precision level in words. This must be sufficient
!                   to recover the relation.
!     ndr   int     Log10 of the min dynamic range at detection, typically 20
!                   or 30. A detected relation is not deemed reliable unless
!                   this is exceeded.
!     nrb   int     Log10 of max size (Euclidean norm) of acceptable relation,
!                   typically 100 or 200. Run will abort if this is exceeded.
!     nep   int     Log10 of tolerance for full precision relation detection.
!     x     mp_real Input mp_real vector.
!     iq    int     Output flag: 0 (unsuccessful) or 1 (successful).
!     r     mp_real Output integer relation vector, if successful.

!   The following parameters are set in this routine:
!     ipi   int     Iteration print interval when idb >= 2; default = 100.
!     ipm   int     Iteration  check interval for MP iterations; default = 10.
!     itm   int     Maximum iteration count; default = 10^5. Run will be aborted
!                   if this is exceeded.
use mpmodule
implicit none
integer, intent(in):: idb, n, nwds, ndr, nrb, nep
integer, parameter:: ipi = 100, ipm = 500, itm = 100000
real (mprknd), external:: bounddp, second
integer i, iq, it, izm, j, j1, n1, n2, n3, n4
real (mprknd) d1, d2, d3, d4, rn, tm0, tm1, times(2)
type (mp_real) eps, b(n,n), h(n,n), r(n), x(n), y(n), t1, t2, t3, t4

!   Initialize.

if (idb >= 2) write (6, 1) n
1 format ('PSLQ1 integer relation detection: n =',i5)
iq = 0
it = 0
rn = mpreal (0.d0, nwds)
eps = mpreal (10.d0, nwds) ** nep

do i = 1, 2
  times(i) = 0.d0
enddo

if (idb >= 2) write (6, 2) it
2 format ('Iteration',i8,2x,'MP initialization')
tm0 = second ()
call initmp (idb, n, nwds, eps, b, h, x, y, izm)
tm1 = second ()
times(1) = tm1 - tm0
if (izm .ne. 0) goto 110

!   MP iterations.

if (idb >= 2) write (6, 3) it
3 format ('Iteration',i8,2x,'Start MP iterations')

100 continue

!   Perform one MP iteration.

it = it + 1
if (idb == 3 .or. idb >= 2 .and. mod (it, ipi) == 0) write (6, 4) it
4 format ('Iteration',i8)
tm0 = second ()
call itermp (idb, it, n, nwds, eps, b, h, y, izm)
tm1 = second ()
times(2) = times(2) + (tm1 - tm0)

!   Test conditions on itermp output flag izm:
!   0: MP update was uneventful; periodically output min, max, norm; continue.
!   1: A small value was found in MP update; go output relation.
!   2: Precision is exhausted; quit.

if (izm == 0) then

!   Periodically output min, max and norm.

  if (mod (it, ipm) == 0) then

!   Find min and max absolute value of y vector.

    call minmax (n, nwds, y, t1, t2)
    if (idb >= 2) then
      call mpdecmd (t1, d1, n1)
      call mpdecmd (t2, d2, n2)
      write (6, 5) it, d1, n1, d2, n2
5     format ('Iteration',i8,2x,'Min, max of y =',f11.6,'e',i6,f11.6,'e',i6)
    endif

!   Compute norm bound.

    d1 = bounddp (n, h)
    if (d1 == -1.d0) goto 120
    rn = max (rn, d1)

    if (idb >= 2) then
      write (6, 6) it, d1, rn
6     format ('Iteration',i8,2x,'Norm bound =',1p,d15.6,4x,'Max. bound =', &
      1p,d15.6)
    endif

!   Check if iteration limit or norm bound limit is exceeded; if so, quit.

    if (it > itm) then
      if (idb >= 1) write (6, 7) itm
7     format ('Iteration limit exceeded',i8)
      goto 120
    endif
    if (log10 (rn) > nrb) then
      if (idb >= 1) write (6, 8) nrb
8     format ('Norm bound limit exceeded.',i5)
      goto 120
    endif
  endif
  goto 100
elseif (izm == 1) then
  goto 110
elseif (izm == 2) then
  goto 120
endif

110 continue

!   A relation has been detected.

tm0 = second ()
t1 = mpreald (1.d300, nwds)
t2 = mpreal (0.d0, nwds)
t3 = mpreal (0.d0, nwds)
t4 = mpreal (0.d0, nwds)

!   Select the relation corresponding to the smallest y entry and compute norm.

do j = 1, n
  if (abs (y(j)) < t1) then
    j1 = j
    t1 = abs (y(j))
  endif
  t2 = max (t2, abs (y(j)))
enddo

t3 = t1 / t2

do i = 1, n
  r(i) = b(i,j1)
  t4 = t4 + r(i) ** 2
enddo

t4 = sqrt (t4)
call mpdecmd (t4, d4, n4)
d4 = d4 * 10.d0 ** n4
d1 = bounddp (n, h)
if (d1 == -1.d0) goto 120
rn = max (rn, d1)

!   Output the final norm bound and other info.

if (idb >= 1) then
  call mpdecmd (t1, d1, n1)
  call mpdecmd (t2, d2, n2)
  call mpdecmd (t3, d3, n3)
  write (6, 9) it, d1, n1, d2, n2, d3, n3, rn
9 format ('Iteration',i8,2x,'Relation detected'/ &
  'Min, max, ratio of y =',0p,f11.6,'e',i5,f11.6,'e',i5,f11.6,'e',i5/ &
  'Max. bound =',1p,d15.6)
  write (6, 10) j1, d4, d1, n1
10 format ('Index of relation =',i4,3x,'Norm =',1p,d15.6,3x, &
  'Residual =',0p,f11.6,'e',i5)
endif

!   If run was successful, set iq = 1.

if (d3 == 0.d0) n3 = nep
if (n4 <= nrb .and. -n3 >= ndr) then
  iq = 1
else
  if (idb >= 2) write (6, 11)
11 format ('Relation is too large or insufficient dynamic range.')
endif

120 continue

!   Output CPU run times and return.

if (idb >= 2) write (6, 12) times
12 format ('CPU times:'/(5f12.2))

return
end subroutine pslq1

!------------------------------

!   First-level subroutines.

subroutine minmax (n, nwds, y, y1, y2)

!   This returns min|y_k| and max|y_k| using MP precision.
!   Input: n, nwds, y.
!   Output: y1, y2.

use mpmodule
implicit none
integer, intent(in):: n, nwds
type (mp_real), intent(in):: y(n)
type (mp_real), intent(out):: y1, y2
integer i
type (mp_real) t1, t2, t3

t1 = mpreald (1.d300, nwds)
t2 = mpreald (0.d0, nwds)

!   Find the min and max absolute value in the y vector.

do i = 1, n
  t3 = abs (mpreal (y(i), nwds))
  t1 = min (t1, t3)
  t2 = max (t2, t3)
enddo

y1 = t1
y2 = t2
return
end subroutine minmax

subroutine initmp (idb, n, nwds, eps, b, h, x, y, izm)

!   This initializes MP arrays at the beginning.
!   This is performed in full precision.
!   Input: idb, n, nwds, eps, x.
!   Output: b, h, y, izm.

use mpmodule
implicit none
integer, intent(in):: idb, n, nwds
type (mp_real), intent(in):: eps, x(n)
type (mp_real), intent(out):: b(n,n), h(n,n), y(n)
integer, intent(out):: izm
integer i, zm, j, k, n1
real (mprknd) d1
type (mp_real) s(n), t1

if (idb >= 3) then
  write (6, 1)
1 format ('initmp: Input x vector:')
  call matoutmd (1, n, x)
endif
izm = 0

!   Set b to the identity matrix.

do j = 1, n
  do i = 1, n
    b(i,j) = mpreal (0.d0, nwds)
  enddo

  b(j,j) = mpreal (1.d0, nwds)
enddo

t1 = mpreal (0.d0, nwds)

!   Compute the x vector, the square root of the partial sum of squares of x,
!   and the y vector, which is the normalized x vector.

do i = n, 1, -1
  t1 = t1 + x(i) ** 2
  s(i) = sqrt (t1)
enddo

t1 = 1.d0 / s(1)

do i = 1, n
  y(i) = t1 * x(i)
  s(i) = t1 * s(i)
enddo

!   Compute the initial h matrix.

do j = 1, n - 1
  do i = 1, j - 1
    h(i,j) = mpreal (0.d0, nwds)
  enddo

  h(j,j) = s(j+1) / s(j)
  t1 = y(j) / (s(j) * s(j+1))

  do i = j + 1, n
!      t1 = y(i) * t1
      h(i,j) = - y(i) * t1
  enddo
enddo

!  Perform reduction on h, updating y, b and h.

do i = 2, n
  do j = i - 1, 1, -1
    t1 = anint (h(i,j) / h(j,j))
    if (t1 .ne. 0.d0) then
      y(j) = y(j) + t1 * y(i)

      do k = i, n
        b(k,j) = b(k,j) + t1 * b(k,i)
      enddo

      do k = 1, j
        h(i,k) = h(i,k) - t1 * h(j,k)
      enddo
    endif
  enddo
enddo

t1 = mpreald (1.d300, nwds)

do i = 1, n
  t1 = min (t1, abs (y(i)))
enddo

if (t1 <= eps) then
  if (idb >= 2) then
    call mpdecmd (t1, d1, n1) 
    write (6, 2) 0, d1, n1
2   format ('Iteration',i8,3x,'initmp: Small value in y =',f11.6,'e',i5)
  endif
  izm = 1
endif

if (idb >= 3) then
  write (6, 3)
3 format ('initmp: Initial y vector:')
  call matoutmd (1, n, y)
  write (6, 4)
4 format ('initmp: Initial h matrix:')
  call matoutmd (n, n - 1, h)
endif

return
end subroutine initmp

subroutine itermp (idb, it, n, nwds, eps, b, h, y, izm)

!   This performs one iteration of the PSLQ algorithm using MP arithmetic.
!   This is performed in medium precision.
!   Input: idb, it, n, nwds, eps, b, h, y.
!   Output: b, h, y, izm.

use mpmodule
implicit none
integer, intent(in):: idb, it, n, nwds
type (mp_real), intent(in)::  eps
type (mp_real), intent(inout):: b(n,n), h(n,n), y(n)
integer, intent(out):: izm
integer, parameter:: ntl = 72
integer i, im, im1, j, j1, k, n1
real (mprknd) d1
type (mp_real) gam, teps, t1, t2, t3, t4

teps = 2.d0 ** ntl * eps
izm = 0
t1 = mpreal (0.d0, nwds)
gam = sqrt (mpreald (4.d0, nwds) / mpreald (3.d0, nwds))

!   Find  im = i  such that  gam^i * |h(i,i)|  is maximized.

do i = 1, n - 1
  t2 = gam ** i * abs (h(i,i))
  if (t2 > t1) then
    im = i
    t1 = t2
  endif
enddo

im1 = im + 1

! Exchange the im and im+1 entries of y, rows of h, and columns of b.

t1 = y(im)
y(im) = y(im1)
y(im1) = t1

do i = 1, n
  t1 = b(i,im)
  b(i,im) = b(i,im1)
  b(i,im1) = t1
enddo

do i = 1, n - 1
  t1 = h(im,i)
  h(im,i) = h(im1,i)
  h(im1,i) = t1
enddo

!  Update h with permutation produced above.

if (im <= n - 2) then
  t1 = h(im,im)
  t2 = h(im,im1)
  t3 = sqrt (t1 ** 2 + t2 ** 2)
  t1 = t1 / t3
  t2 = t2 / t3

  do i = im, n
    t3 = h(i,im)
    t4 = h(i,im1)
    h(i,im) = t1 * t3 + t2 * t4
    h(i,im1) = - t2 * t3 + t1 * t4
  enddo
endif

t2 = mpreal (0.d0, nwds)

!  Perform reduction on h, updating y, a, b and h.

do i = im1, n
  j1 = min (i - 1, im1)

  do j = j1, 1, -1
    t1 = anint (h(i,j) / h(j,j))
    if (t1 .ne. 0.d0) then
      y(j) = y(j) + t1 * y(i)

      do k = 1, n
        b(k,j) = b(k,j) + t1 * b(k,i)
      enddo

      do k = 1, j
        h(i,k) = h(i,k) - t1 * h(j,k)
      enddo
    endif
  enddo
enddo

!   Find the min of y.

t1 = mpreald (1d300, nwds)

do j = 1, n
  if (abs (y(j)) < t1) then
    j1 = j
    t1 = abs (y(j))
  endif
enddo

!  Find the largest entry of b in the same column as the smallest y.

t2 = mpreal (0.d0, nwds)

do i = 1, n
  t2 = max (t2, abs (b(i,j1)))
enddo

if (t1 <= t2 * teps) then
  if (idb >= 2) then
    call mpdecmd (t1, d1, n1) 
    write (6, 1) it, d1, n1
1   format ('Iteration',i8,2x,'itermp: Small value in y =',f11.6,'e',i5)
  endif
  if (t1 <= t2 * eps) then
    izm = 1
  else
    if (idb >= 1) write (6, 2) it
2   format ('Iteration',i8,2x,'itermp: Precision exhausted.')
    izm = 2
  endif
endif

if (idb >= 3) then
  write (6, 4)
4 format ('itermp: Updated y:')
  call matoutmd (1, n, y)
  write (6, 5)
5 format ('itermp: Updated b matrix:')
  call matoutmd (n, n, b)
  write (6, 6)
6 format ('itermp: Updated h matrix:')
  call matoutmd (n, n - 1, h)
endif

return
end subroutine itermp

!------------------------------

!   Second- and third-level subroutines.

real (mprknd) function bounddp (n, h)

!   This computes the norm bound using DP arithmetic.

use mpmodule
implicit none
integer, intent(in):: n
type (mp_real), intent(in):: h(n,n)
integer i, j
real (mprknd) dh(n,n), t1

do j = 1, n - 1
  do i = 1, n
    dh(i,j) = h(i,j)
  enddo
enddo

call lqdp (n, n - 1, dh)
t1 = 0.d0

do i = 1, n - 1
  t1 = max (t1, abs (dh(i,i)))
enddo

if (t1 < 1.d-300) then
  write (6, 1)
1 format ('bounddp: dh matrix too small -- use pslqm3 program instead.')
  bounddp = -1.d0
else
  bounddp = 1.d0 / t1
endif

return
end function bounddp

subroutine lqdp (n, m, dh)

!   This performs an LQ decomposition on the DP matrix dh.  It is a simplified
!   and transposed adaptation of the subroutine dqrdc from Linpack.
!   Input: n, m, dh.
!   Output: dh.

use mpmodule
implicit none
integer, intent(in):: n, m
real (mprknd), intent(inout):: dh(n,m)
integer i, j, l, lup, ml
real (mprknd) nrmxl, one, t, zero

zero = 0.d0
one = 1.d0
lup = min (m,n)

!   Perform the householder reduction of dh.

do l = 1, lup
  if (l == m) go to 280

!   Compute the householder transformation for column l.

  ml = m - l
  t = zero

  do i = 0, ml
    t = t + dh(l,l+i) ** 2
  enddo

  nrmxl = sqrt (t)
  if (nrmxl == zero) go to 270
  if (dh(l,l) /= zero) nrmxl = sign (nrmxl, dh(l,l))
  t = one / nrmxl

  do i = 0, ml
    dh(l,l+i) = t * dh(l,l+i)
  enddo

  dh(l,l) = one + dh(l,l)

!   Apply the transformation to the remaining columns, updating the norms.

  do j = l + 1, n
    t = zero

    do i = 0, ml
      t = t + dh(l,l+i) * dh(j,l+i)
    enddo

    t = - t / dh(l,l)

    do i = 0, ml
      dh(j,l+i) = dh(j,l+i) + t * dh(l,l+i)
    enddo
  enddo

!   Save the transformation.

  dh(l,l) = - nrmxl
270 continue
280 continue
enddo

!   Zero dh above the diagonal.

do j = 1, m
  do i = 1, j - 1
    dh(i,j) = 0.d0
  enddo
enddo

return
end subroutine lqdp

subroutine matoutmd (n1, n2, a)

!   This outputs the MP matrix a as a DP matrix.

use mpmodule
implicit none
integer, intent(in):: n1, n2
type (mp_real), intent(in):: a(n1,n2)
integer i, j, ix(n2)
real (mprknd) dx(n2)

do i = 1, n1
  write (6, 1) i
1 format ('Row',i3)

  do j = 1, n2
    call mpdecmd (a(i,j), dx(j), ix(j))
  enddo

  write (6, 2) (dx(j), ix(j), j = 1, n2)
2 format (4(f13.8,'e',i5))
enddo

return
end subroutine matoutmd

subroutine matoutmp (n1, n2, a)

!   This outputs the MP matrix a.  It may be used in place of calls to matoutmd
!   in the code above if greater accuracy is desired in debug output.

use mpmodule
implicit none
integer, intent(in):: n1, n2
type (mp_real), intent(in):: a(n1,n2)
integer i, j

do i = 1, n1
  write (6, 1) i
1 format ('Row',i3)

  do j = 1, n2
    call mpwrite (6, 80, 60, a(i,j))
  enddo
enddo

return
end subroutine matoutmp

