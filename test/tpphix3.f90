!*****************************************************************************

!  program tpphix3

!  Revision date:  15 Oct 2021

!  AUTHOR:
!   David H. Bailey
!   Lawrence Berkeley National Lab (retired) and University of California, Davis
!   Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!   All software in this package (c) 2021 David H. Bailey.
!   By downloading or using this software you agree to the copyright, disclaimer
!   and license agreement in the accompanying file DISCLAIMER.txt.

!  DESCRIPTION:
!   This computes the Poisson phi function, as described in the paper 
!   David H. Bailey, Jonathan M. Borwein, Jason Kimberley and Watson Ladd, "Computer
!   discovery and analysis of large Poisson polynomials," Experimental Mathematics,
!   27 Aug 2016, vol. 26, pg. 349-363, preprint available at
!   https://www.davidhbailey.com/dhbpapers/poisson-res.pdf.

!   PSLQM3 parameters set below; all are default (4-byte) integer, except fname
!     is character(64):
!     idb   Debug level (0-4); default = 2.
!     n     Integer relation vector length; default = 49.
!           For minimal polynomial problems, n = 1 + polynomial degree.
!     ndp   Full precision level in digits; default = 2500.
!           ***Must be <= mpipl in module MPFUNF.
!           ***Must be >= ndpm + precision level required to find relation.
!     ndpm  Medium precision level in digits; default = max (ndp/10, 100).
!           ***Must be <= mpiplm in module MPFUNF.
!           ***Must exceed the dynamic range of the x vector.
!     ndr   Log10 of the minimum dynamic range in y at detection; default = 30.
!           A detected relation is not deemed reliable unless this is exceeded.
!     nep   Log10 of full precision epsilon for detections; default = 30 - ndp.
!           ***Must not be smaller than the accuracy of input x vector. In other
!           words, if data is accurate to within 10^(-200), then nep > -200.
!     nepm  Log10 of medium precision epsilon; default = 20 - ndpm.
!     nrb   Log10 of maximum size (Euclidean norm) of acceptable relation;
!           default = 200. Run will be aborted if this is exceeded.
!     nrs   0: Restart file is not read or written; default = 0.
!           1: Start new run, and periodically write restart file.
!           2: Read restart file, run and periodically write restart file.
!           The input restart file is overwritten, so make a backup copy.
!     nwds  Full precision level in words; default = int (ndp/mpdpw + 2).
!     nwdsm Medium precision level in words; default = int (ndpm/mpdpw + 2).
!     fname Filename for restart file; default = 'pslqm3.rst'. If more than 
!           one job is run in same directory, use distinct filenames.
!   For some additional details, see comments at start of subroutine pslqm3.

!   TPPHIX3 parameters set below; all are default (4-byte) integer:
!     ipal  0: The palindromic technique (for even kd) is NOT implemented.
!           In this case, the parameter n = 1 + polynomial degree.
!           1: The palindromic technique (for even kd) IS implemented (default).
!           In this case, the parameter n = 1 + half polynomial degree,
!           and the degree of full output polynomial = 2*n-2.
!     kd    Denominator of rational arguments of phi_2; default = 28.
!     kp    Numerator of first argument; default = 1.
!     kq    Numerator of second argument; default = 1.
!     lcx   Size of character(1) array chr1; 
!           default = 64 * int (1 + (2*n + ndp)/(64*n)).
!           ***Must be a multiple of 64. Must be <= lmx.
!           ***Must be > size in digits of largest result coefficient.
!     lmx   Size of line1; default = 131072. 
!           ***Must be > char length of full relation in Mathematica format.
!     mpi   Multiplier of pi; default = 8.
!     nep1  Log10 of full precision epsilon for computation of alpha;
!           default = -ndp.

program tpphix3
use mpmodule
implicit none

!   PSLQM4 parameters:

integer idb, n, ndp, ndpm, ndr, nep, nepm, nrb, nrs, nwds, nwdsm
character(64) fname
parameter (idb = 2, n = 49, ndp = 2500, ndpm = max (ndp / 10, 100), &
  ndr = 30, nep = 30 - ndp, nepm = 20 - ndpm, nrb = 200, nrs = 0, &
  nwds = int (ndp / mpdpw + 2), nwdsm = int (ndpm / mpdpw + 2), &
  fname = 'pslqm3.rst')

!   TPPHIX3 parameters:

integer ipal, kd, kp, kq, lcx, lmx, mpi, nep1
parameter (ipal = 1, kd = 28, kp = 1, kq = 1, &
  lcx = 64 * int (1 + (2*n + ndp)/(64*n)), lmx = 131072, &
  mpi = 8, nep1 = -ndp)

integer i, iq, i1, i2, j, j1, k, l1, nn, n1
real (mprknd) second, tm0, tm1, tm2, tm3
external second
type (mp_real) alpha, beta, eps1, pi, qq, xx, yy
type (mp_complex) zz, c1, c2, c3, c4, theta1, theta2, theta3, theta4
external theta1, theta2, theta3, theta4
integer lnm(n), lnm2(2*n)
character(1) chr1(lcx)
character(64) form4, nam(n), nam2(2*n), namx
character(lmx) line1
type (mp_real) r(n), r2(2*n), x(n)
save

!   Uncomment this line for MPFUN20-Fort when precision > 20,000 digits.

! call mpinit (nwds)

!   Check to see if precision level in MPFUNF is high enough.

if (ndp > mpipl .or. ndpm > mpiplm) then
  write (6, 1) ndp, ndpm
1 format ('Increase the default standard precision in MPFUNF to at least', &
    i8,' digits'/'and the default medium precision in MPFUNF to at least', &
    i8,' digits,'/'then recompile the library.' )
  stop
endif

nn = n

!   Check whether the palindromic option is appropriate.

if (ipal /= 0 .and. mod (kd, 2) == 1) then
  write (6, 2)
2 format ('The palindromic technique may not be used when kd is odd.')
  stop
endif

!   Check whether the parameters lcx and lmx are valid:

if (mod (lcx, 64) /= 0 .or. lcx > lmx) then
  write (6, 21)
21 format ('The parameters lcx and lmx are invalid.')
  stop
endif

write (6, 3) ipal, kd, kp, kq, mpi, nn, ndp, ndpm, nep, nepm, nrs
3 format ('Poisson Phi_2 computation and analysis:'/ &
  'ipal =', i4,'; kd =',i4,'; kp =',i4,'; kq =',i4,'; mpi =',i4,'; nn =',i4/ &
  'ndp =',i6,'; ndpm =',i6,'; nep =',i6,'; nepm =',i6,'; nrs =',i6)

pi = mppi (nwds)
eps1 = mpreald (10.d0, nwds) ** nep1

!   Compute alpha (but no need if nrs = 2).

tm0 = second ()
if (nrs <= 1) then
  qq = exp (- pi)
  xx = mpreal (dble (kp), nwds) / mpreal (dble (kd), nwds)
  yy = mpreal (dble (kq), nwds) / mpreal (dble (kd), nwds)
  zz = 0.5d0 * pi * mpcmplx (yy, xx, nwds)
  c1 = theta1 (zz, qq, eps1, nwds)
  c2 = theta2 (zz, qq, eps1, nwds)
  c3 = theta3 (zz, qq, eps1, nwds)
  c4 = theta4 (zz, qq, eps1, nwds)
  alpha = abs (c2 * c4 / (c1 * c3)) ** (mpi / 2)
else
  alpha = mpreald (1.d0, nwds)
endif
tm1 = second ()
write (6, 4) tm1 - tm0
4 format ('Alpha CPU time =',f12.2/'Alpha =')
call mpwrite (6, ndp+20, ndp, alpha)

!   Construct input x vector for PSLQM3.

if (ipal == 0) then

!   Case: palindromic property is NOT used.

  x(1) = mpreal (1.d0, nwds)

  do i = 2, nn
    x(i) = alpha * x(i-1)
  enddo

  nam(1) = '1'
  lnm(1) = 1

  do i = 2, nn
    if (i <= 10) then
      write (namx, '("al^",i1)') i - 1
      nam(i) = namx(1:4)
      lnm(i) = 4
    elseif (i <= 100) then
      write (namx, '("al^",i2)') i - 1
      nam(i) = namx(1:5)
      lnm(i) = 5
    elseif (i <= 1000) then
      write (namx, '("al^",i3)') i - 1
      nam(i) = namx(1:6)
      lnm(i) = 6
    elseif (i <= 10000) then
      write (namx, '("al^",i4)') i - 1
      nam(i) = namx(1:7)
      lnm(i) = 7
    else
      stop
    endif   
  enddo

  do i = 1, nn
    lnm2(i) = lnm(i)
    nam2(i) = nam(i)
  enddo
elseif (ipal == 1) then

!   Case: palindromic property IS used.

  beta = alpha + 1.d0 / alpha
  x(1) = mpreal (1.d0, nwds)

  do i = 2, nn
    x(i) = beta * x(i-1)
  enddo

  nam(1) = '1'
  lnm(1) = 1

  do i = 2, nn
    if (i <= 10) then
      write (namx, '("beta^",i1)') i - 1
      nam(i) = namx(1:6)
      lnm(i) = 6
    elseif (i <= 100) then
      write (namx, '("beta^",i2)') i - 1
      nam(i) = namx(1:7)
      lnm(i) = 7
    elseif (i <= 1000) then
      write (namx, '("beta^",i3)') i - 1
      nam(i) = namx(1:8)
      lnm(i) = 8
    elseif (i <= 10000) then
      write (namx, '("beta^",i4)') i - 1
      nam(i) = namx(1:9)
      lnm(i) = 9
    else
      stop
    endif   
  enddo

  nam2(1) = '1'
  lnm2(1) = 1

  do i = 2, 2 * nn
    if (i <= 10) then
      write (namx, '("al^",i1)') i - 1
      nam2(i) = namx(1:4)
      lnm2(i) = 4
    elseif (i <= 100) then
      write (namx, '("al^",i2)') i - 1
      nam2(i) = namx(1:5)
      lnm2(i) = 5
    elseif (i <= 1000) then
      write (namx, '("al^",i3)') i - 1
      nam2(i) = namx(1:6)
      lnm2(i) = 6
    elseif (i <= 10000) then
      write (namx, '("al^",i4)') i - 1
      nam2(i) = namx(1:7)
      lnm2(i) = 7
    else
      stop
    endif   
  enddo
endif

!   Perform relation search using pslqm3.

tm2 = second ()
call pslqm3 (idb, nn, nwds, nwdsm, ndr, nrb, nrs, fname, nep, nepm, x, iq, r)
tm3 = second ()
write (6, 5) tm3 - tm2, tm3 - tm0
5 format ('PSLQM3 CPU time =',f12.2/'Total CPU time =',f12.2)

!   Output relation in two formats.

if (iq == 1) then

!   Produce format used below for output.

  form4 = '(64a1'
  i1 = 5

  do i = 65, lcx, 64
    form4(i1+1:i1+9) =  '/64a1'
    i1 = i1 + 5
  enddo

  form4(i1+1:i1+9) = '," * ",a)'
  i1 = i1 + 9

!  If palindromic property is used, output recovered half relation.

  if (ipal == 1) then
    write (6, 6)
6   format (/'Recovered half relation: 0 =')
  
    do i = 1, nn
      if (r(i) .ne. 0.d0) then
        call mpfform (r(i), lcx, 0, chr1)

        do j = 1, lcx
          if (chr1(j) /= ' ') goto 110
        enddo

110     continue

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

!    If palindromic property is not used, we are done; if it is used, then
!    call doublep to expand to full polynomial.

  if (ipal == 0) then
    n1 = nn

    do i = 1, nn
      r2(i) = r(i)
    enddo
  elseif (ipal == 1) then
    n1 = 2 * nn - 1
    call doublep (nwds, 2*nn - 1, r, r2)
  endif

  write (6, 7)
7 format (/'Recovered full relation: 0 =')
  l1 = 0

!   Output full polynomial in Mathematica form.

  do i = 1, n1
    if (r2(i) .ne. 0.d0) then
      call mpfform (r2(i), lcx, 0, chr1)

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
      write (6, form4) (chr1(j), j = 1, lcx), nam2(i)(1:lnm2(i))

      if (l1 + 100 > lmx) then
        write (6, '("Mathematica polynomial is too long; increase lmx =",i8)') lmx
        stop
      endif
      line1(l1+1:l1+1) = ' '
      l1 = l1 + 1
      k = lcx - j1
  
      do j = 1, k
        line1(l1+j:l1+j) = chr1(j+j1-1)
      enddo
  
      l1 = l1 + k
      line1(l1+1:l1+lnm2(i)+1) = '*' // nam2(i)(1:lnm2(i))
      l1 = l1 + lnm2(i) + 1
    endif
  enddo

  i1 = 1
  write (6, 8)
8 format ('Output polynomial in Mathematica notation:')

130 continue

  if (i1 + 64 > l1) goto 140
  i2 = index (line1(i1+65:l1), ' ')
  if (i2 == 0) goto 140
  i2 = i1 + 64 + i2
!  write (6, '(a)') line1(i1:i2) // '\ '
  write (6, '(a)') line1(i1:i2)
  i1 = i2 + 1
  goto 130

140 continue

  write (6, '(a)') line1(i1:l1)
endif

if (iq == 1) then
  write (6, '(a)') 'TEST PASSED'
else
  write (6, '(a)') 'TEST FAILED'
endif

stop
end

subroutine doublep (nwds, n, polyh, poly)

!   When the palindromic option is used, this expands result to full polynomial.

use mpmodule
implicit none
integer i, j, k, n, nwds
type (mp_real) polyh(0:n/2), poly(-n/2:n/2), x(-n/2-1:n/2+1), y(-n/2-1:n/2+1)

do i = -n/2-1, n/2+1
  x(i) = mpreal (0.d0, nwds)
enddo

x(0) = polyh(n/2)

do j = 0, n / 2 - 1
  y(-j-1) = x(-j)
  y(j+1) = x(j)
  
  do k = -j, j
    y(k) = x(k-1) + x(k+1)
  enddo
  
  y(0) = y(0) + polyh(n/2-j-1)

  do k = -j-1, j+1
    x(k) = y(k)
  enddo
enddo

do k = -n/2, n/2
  poly(k) = x(k)
enddo

return
end

function theta1 (zz, qq, eps, nwds)

!   This computes the first of the four theta functions.

use mpmodule
implicit none
integer k1, k2, nwds
real (mprknd) ds
type (mp_real) eps, qq, r1
type (mp_complex) theta1, t0, t2, zz

t0 = mpcmplx (cmplx (0.d0, 0.d0, mprknd), nwds)
r1 = sqrt (sqrt (qq))
ds = -1.d0

do k1 = 1, 10001, 2
  ds = - ds
  k2 = k1 ** 2
  t2 = ds * r1 ** k2 * sin (dble (k1) * zz)
  t0 = t0 + t2
  if (abs (t2) < eps) goto 100
enddo

write (6, 1)
1 format ('theta1: loop end error')
stop

100 continue

theta1 = 2.d0 * t0
write (6, 2) k1
2 format ('theta1 complete: k1 =',i6)
return
end

function theta2 (zz, qq, eps, nwds)

!   This computes the second of the four theta functions.

use mpmodule
implicit none
integer k1, k2, nwds
type (mp_real) eps, qq, r1
type (mp_complex) theta2, t0, t2, zz

t0 = mpcmplx (cmplx (0.d0, 0.d0, mprknd), nwds)
r1 = sqrt (sqrt (qq))

do k1 = 1, 10001, 2
  k2 = k1 ** 2
  t2 = r1 ** k2 * cos (dble (k1) * zz)
  t0 = t0 + t2
  if (abs (t2) < eps) goto 100
enddo

write (6, 1)
1 format ('theta2: loop end error')
stop

100 continue

theta2 = 2.d0 * t0
write (6, 2) k1
2 format ('theta2 complete: k1 =',i6)
return
end

function theta3 (zz, qq, eps, nwds)

!   This computes the third of the four theta functions.

use mpmodule
implicit none
integer k1, k2, nwds
type (mp_real) eps, qq, r1
type (mp_complex) theta3, t0, t2, zz

t0 = mpcmplx (cmplx (0.d0, 0.d0, mprknd), nwds)
r1 = qq

do k1 = 1, 10000
  k2 = k1 ** 2
  t2 = r1 ** k2 * cos (2.d0 * dble (k1) * zz)
  t0 = t0 + t2
  if (abs (t2) < eps) goto 100
enddo

write (6, 1)
1 format ('theta3: loop end error')
stop

100 continue

theta3 = cmplx (1.d0, 0.d0, mprknd) + 2.d0 * t0
write (6, 2) k1
2 format ('theta3 complete: k1 =',i6)
return
end

function theta4 (zz, qq, eps, nwds)

!   This computes the fourth of the four theta functions.

use mpmodule
implicit none
integer k1, k2, nwds
real (mprknd) ds
type (mp_real) eps, qq, r1
type (mp_complex) theta4, t0, t2, zz

t0 = mpcmplx (cmplx (0.d0, 0.d0, mprknd), nwds)
r1 = qq
ds = 1.d0

do k1 = 1, 10000
  ds = - ds
  k2 = k1 ** 2
  t2 = ds * r1 ** k2 * cos (2.d0 * dble (k1) * zz)
  t0 = t0 + t2
  if (abs (t2) < eps) goto 100
enddo

write (6, 1)
1 format ('theta4: loop end error')
stop

100 continue

theta4 = cmplx (1.d0, 0.d0, mprknd) + 2.d0 * t0
write (6, 2) k1
2 format ('theta4 complete: k1 =',i6)
return
end

!------------------------------

!   The following code performs the three-level, multi-pair PSLQ algorithm.
!   David H. Bailey    15 Oct 2021

subroutine pslqm3 (idb, n, nwds, nwdsm, ndr, nrb, nrs, fname, nep, nepm, x, iq, r)

!   Arguments are as follows; int means default (4-byte) integer:
!     Name  Type     Description
!   Input:
!     idb   int      Debug flag (0-4); increasing idb produces more output.
!     n     int      Length of input vector x and output relation vector r.
!     nwds  int      Full precision level in words. This must be sufficient
!                    to recover the relation, plus nwdsm words.
!     nwdsm int      Medium precision level in words, typically nwds / 10 or so.
!     ndr   int      Log10 of the min dynamic range of |y_i| at detection, typically 30.
!                    A detected relation is not deemed reliable unless this is met.
!     nrb   int      Log10 of max size (Euclidean norm) of acceptable relation,
!                    typically 100 or 200. Run will abort if this is exceeded.
!     nrs   int      0: Restart file is not read or written; typically 0.
!                    1: Start new run, and periodically write restart file.
!                    2: Read restart file, run and periodically write restart file.
!                    The input restart file is overwritten, so make a backup copy.
!     fname char(64) Filename for restart; typically 'pslqm3.rst'. If more than 
!                    one job is run in same directory, use distinct filenames.
!     nep   int      Log10 of tolerance for full precision detection.
!     nepm  int      Log10 of tolerance for medium precision detection.
!     x     mp_real  Input x vector.

!   Output:
!     iq    int      Output flag: 0 (unsuccessful) or 1 (successful).
!     r     mp_real  Output integer relation vector, if successful; otherwise 0s.

!   The following parameters are set below in this subroutine:
!     ipi   int      Iteration print interval when idb >= 2; default = 500.
!     ipm   int      Iteration check interval for MPM iterations; default = 10.
!     itm   int      Maximum iteration count; default = 10^7. Run is aborted if
!                    this is exceeded. If itm >= 10^8, change all "i8" in
!                    formats in the subroutines below to "i9" or as needed.
!     ndrm  int      Extra precision beyond wy dynamic range for MPM iterations;
!                    default = 25.
!     nsq   int      Size of tables used in iterdp and itermpw; default = 8.
!     dreps double   Tolerance for DP dynamic range check; default = 1d-10.

use mpmodule
implicit none
integer i, i1, idb, imq, ipi, ipm, iq, it, itm, its, izd, izm, izmm, j, j1, &
  n, n1, n2, n3, n4, ndr, ndrm, nep, nepm, nrb, nrs, nsq, nwds, nwdsm, ndpm2, nwdsm2
character(64) fname
parameter (ipi = 500, ipm = 10, itm = 10000000, ndrm = 25, nsq = 8)
real (mprknd) d1, d2, d3, d4, dreps, dplog10, second, &
  tm0, tm1, times(6)
parameter (dreps = 1.d-10)
real (mprknd) da(n,n), db(n,n), dh(n,n), dsa(n,n), dsb(n,n), &
  dsh(n,n), dsyq(n,nsq), dy(n), dsy(n)
type (mp_real) b(n,n), h(n,n), r(n), x(n), y(n), eps, t1, t2
type (mp_realm) bound, dynrange, dynrangem, epsm, epsm2, wa(n,n), &
  wb(n,n), wh(n,n), wsyq(n,nsq), wy(n), wn, w1, w2
external bound, dplog10, dynrange, dynrangem, second

!   Initialize.

if (idb >= 2) write (6, 1) n
1 format ('PSLQM3 integer relation detection: n =',i5)
iq = 0
it = 0
its = 0
imq = 0
izm = 0
izmm = 0
izd = 0
wn = mprealdm (0.d0, nwdsm)
eps = mpreald (10.d0, nwds) ** nep
epsm = mprealdm (10.d0, nwdsm) ** nepm

if (idb >= 2) write (6, 2) it
2 format ('Iteration',i8,2x,'MP initialization')

!   If nrs >= 1, open restart file; if nrs = 2, read restart file. 

if (nrs >= 1) open (12, file = fname, form = 'unformatted')
if (nrs <= 1) then
  do i = 1, 6
    times(i) = 0.d0
  enddo

  tm0 = second ()
  call initmp (idb, n, nwds, b, h, x, y)
  tm1 = second ()
  times(1) = tm1 - tm0
else
  rewind 12
  read (12) it, izm, izmm, izd
  read (12) times
  read (12) x
  read (12) y
  read (12) b
  read (12) h
  rewind 12

!   Find the min and max absolute value of the y vector.
  
  i1 = 0
  t1 = mpreald (1.d300, nwds)
  t2 = mpreald (0.d0, nwds)

  do i = 1, n
    if (abs (y(i)) < t1) then
      i1 = i
      t1 = abs (y(i))
    endif
    t2 = max (t2, abs (y(i)))
  enddo

  if (idb >= 2) then
    call mpdecmd (t1, d1, n1)
    call mpdecmd (t2, d2, n2)
    write (6, 3) it, d1, n1, d2, n2
3   format ('Iteration',i8,2x,'updtmp: Min, max of y =',f11.6,'e',i6, &
      f11.6,'e',i6)
  endif

!   Find the largest entry of b in the same row as the smallest y.

  t2 = mpreald (0.d0, nwds)

  do i = 1, n
    t2 = max (t2, abs (b(i1,i)))
  enddo

!   Check for small value -- if true, relation has already been detected.
!   This can happen if the previous restart run had too large value of eps
!   and failed to detect a good relation.

  if (t1 <= t2 * eps) then
    if (idb >= 2) then
      call mpdecmd (t1, d1, n1) 
      write (6, 4) it, d1, n1
4     format ('Iteration',i8,2x,'updtmp: Small value in y =',f11.6,'e',i6)
    endif
    goto 150
  endif
endif

!+  Main loop starts here.

100 continue

!   Check the dynamic range of y vector using medium precision.

w1 = dynrange (n, nwdsm, y)
if (idb >= 2) then
  call mpdecmd (w1, d1, n1)
  write (6, 5) it, d1, n1
5 format ('Iteration',i8,2x,'Min/max ratio of y =',f11.6,'e',i6)
endif

!   Initialize MPM arrays from MP arrays.

if (idb >= 3) write (6, 6) it
6 format ('Iteration',i8,2x,'MPM initialization')
tm0 = second ()
call initmpm (idb, n, nsq, nwdsm, wa, wb, wh, wy, wsyq, h, y)
tm1 = second ()
times(2) = times(2) + (tm1 - tm0)

110 continue

!   Check if dynamic range of wy vector is too great for DP or QP iterations
!   (which is often the case at the start of the run), or if the previous
!   DP or QP iteration or was aborted. If so, then go perform MPM iterations.

w1 = dynrangem (n, nwdsm, wy)
if (idb >= 3) then
  call mpdecmd (w1, d1, n1)
  write (6, 7) it, d1, n1
7 format ('Iteration',i8,2x,'Min/max ratio of wy =',f11.6,'e',i6)
endif
if (w1 < mprealdm (dreps, nwdsm) .or. izd == 2) then
    goto 130
endif

!+  Start of DP iteration section.
!   Initialize DP arrays from MPM arrays.

if (idb >= 3) write (6, 8) it
8 format ('Iteration',i8,2x,'Start DP iterations')
call initdp (idb, n, nsq, nwdsm, da, db, dh, dy, dsyq, wh, wy)

!   Save DP arrays and iteration count.

call savedp (n, da, db, dh, dy, dsa, dsb, dsh, dsy)
its = it

!   Perform an LQ decomposition on DH.

call lqdp (n, n - 1, dh)

!+  Start of DP iteration loop.

120 continue

it = it + 1
if (idb >= 4 .or. idb >= 2 .and. mod (it, ipi) == 0) write (6, 9) it
9 format ('Iteration',i8)
tm0 = second ()
call iterdp (idb, it, n, nsq, da, db, dh, dsyq, dy, imq, izd)
tm1 = second ()
times(3) = times(3) + (tm1 - tm0)

!   Test conditions on iterdp output flag izd:
!   0: Iteration was uneventful; periodically save arrays and continue.
!   1: Relation found or DP precision exhausted; perform MPM update.
!   2: Very large value appeared in DA or DB; abort DP iter and do MPM update.

if (izd == 0) then
  if (mod (it - its, ipm) == 0) then
    call savedp (n, da, db, dh, dy, dsa, dsb, dsh, dsy)
    its = it
  endif
  goto 120
else

!   Check if DP iteration was aborted above; if so, then revert to previous data.

  if (izd == 2) then
    it = its
    call savedp (n, dsa, dsb, dsh, dsy, da, db, dh, dy)
  endif

!   Update the MPM arrays from the DP arrays.

  if (idb >= 3) write (6, 10) it
10 format ('Iteration',i8,2x,'MPM update from DP arrays')
  tm0 = second ()
  call updtmpm (idb, it, n, nwdsm, epsm, da, db, dreps, wa, wb, wh, wy, izmm)
  tm1 = second ()
  times(4) = times(4) + (tm1 - tm0)

!   Test conditions on updtmpm output flag izmm:
!   0: MPM update was uneventful; continue.
!   1: MPM update found possible relation, or MPM precision is exhausted, or
!      wy min/max ratio is too small for DP iterations; perform MP update.
!   2: Very large value appeared in wa or wb array; abort run.

  if (izmm == 0 .and. izd /= 2) then
    goto 110
  elseif (izmm == 1 .or. izd == 2) then

!   MPM update found relation, or exhausted MPM precision, or wy min/max ratio
!   too small for DP iterations, or DP iteration produced a very large value.
!   Update the MP arrays from the MPM arrays.

    if (idb >= 2) write (6, 11) it
11  format ('Iteration',i8,2x,'MP update')
    tm0 = second ()
    call updtmp (idb, it, n, nwds, wa, wb, eps, b, h, y, izm)
    tm1 = second ()
    times(5) = times(5) + (tm1 - tm0)

!   If checkpoint-restart is enabled, then save restart data.

    if (nrs >= 1) then
      call saverest (12, n, it, izm, izmm, izd, times, x, y, b, h)
    endif

!   Compute norm bound, using MPM precision.

    w1 = bound (n, nwdsm, wh)
    call mpdecmd (w1, d3, n3)
    wn = max (wn, w1)
    call mpdecmd (wn, d4, n4)
    if (idb >= 2) then
      write (6, 12) it, d3, n3, d4, n4
12    format ('Iteration',i8,2x,'Norm bound =',f11.6,'e',i5,3x, &
        'Max. bound =',f11.6,'e',i5)
    endif
    if (dplog10 (wn) > nrb) then
      if (idb >= 1) write (6, 13) nrb
13    format ('Norm bound limit exceeded.',i5)
      goto 160
    endif
    if (it > itm) then
      if (idb >= 1) write (6, 14) itm
14    format ('Iteration limit exceeded',i8)
      goto 160
    endif

!   Test conditions on updtmp output flag izm:
!   0: MP update was uneventful; go to tag 100 above to initialize MPM.
!   1: A small value was found in MP update; go to output relation.
!   2: Precision is exhausted; abort run.

    if (izm == 0) then
      goto 100
    elseif (izm == 1) then
      goto  150
    elseif (izm == 2) then
      goto 160
    endif
  elseif (izmm == 2) then

!   A very large value was found in wa or wb arrays -- abort run.

    goto 160
  endif
endif

!+  Start of MPM iteration section.

130 continue

!   Perform MPM iterations with a working precision level that is at
!   least ndrm digits greater than the dynamic range of the wy vector, and
!   at least four words.

izd = 0
its = it
w1 = dynrangem (n, nwdsm, wy)
call mpdecmd (w1, d1, n1)
ndpm2 = ndrm - n1
nwdsm2 = max (int (ndpm2 / mpdpw) + 2, 4)
ndpm2 = nint ((nwdsm2 - 1) * mpdpw)
epsm2 = mprealm (2.d0, nwdsm2) ** (70 - nwdsm2 * mpnbt)

if (idb >= 2) write (6, 15) it, nwdsm2, ndpm2
15 format ('Iteration',i8,2x,'Start MPM iter: precision =',i4,' words,',i5,' digits')

call setprec (nwdsm2, n, nsq, wa, wb, wh, wsyq, wy)

!   Perform LQ decomposition using variable MPM precision.

tm0 = second ()
call lqmpm (n, n - 1, nwdsm2, wh)
tm1 = second ()
times(2) = times(2) + (tm1 - tm0)

!+  Start of MPM iteration loop.

140 continue

it = it + 1
if (idb >= 2) write (6, 16) it
16 format ('Iteration',i8)
tm0 = second ()
call itermpm (idb, it, n, nsq, nwdsm2, epsm2, wa, wb, wh, wsyq, wy, imq, izmm)
tm1 = second ()
times(6) = times(6) + (tm1 - tm0)

!   Periodically check if dynamic range of wy vector is least dreps; if so,
!   resume DP iterations; otherwise set flag for an MP update.

if (izmm == 0 .and. mod (it - its, ipm) == 0) then
  w1 = dynrangem (n, nwdsm, wy)
  call mpdecmd (w1, d1, n1)
  if (idb >= 2) write (6, 17) it, d1, n1
17 format ('Iteration',i8,2x,'Min/max ratio of wy =',f11.6,'e',i6)
  if (w1 > mprealdm (dreps, nwdsm2)) then 
    if (idb >= 2) write (6, 18) it
18  format ('Iteration',i8,2x,'Return to DP iterations')
    izmm = 1
  else
    izmm = 1
  endif
endif

!   Test conditions on itermpm output flag izmm:
!   0: Iteration was uneventful; continue.
!   1: A small value was found or MPM precision exhausted; do MP update.
!   2: A very large value appeared in WA or WB; abort run.

if (izmm == 0) then
  goto 140
elseif (izmm == 1) then

!   Update the MP arrays from the MPM arrays.

  if (idb >= 2) write (6, 19) it
19 format ('Iteration',i8,2x,'MP update')
  tm0 = second ()
  call updtmp (idb, it, n, nwds, wa, wb, eps, b, h, y, izm)
  tm1 = second ()
  times(5) = times(5) + (tm1 - tm0)

!   If checkpoint-restart is enabled, then save restart data.

    if (nrs >= 1) then
      call saverest (12, n, it, izm, izmm, izd, times, x, y, b, h)
    endif

!   Test conditions on updtmp output flag izm:
!   0: MP update was uneventful; goto tag 100 above.
!   1: A small value was found in MP update; go output relation.
!   2: Precision is exhausted; abort run.

  if (izm == 0) then
    goto 100
  elseif (izm == 1) then
    goto 150
  elseif (izm == 2) then
    goto 160
  endif
elseif (izmm == 2) then

!   During itermpm, a very large entry was produced in wa or wb; abort run.

  goto 160
endif

!+  Start of detection section.

150 continue

!   A relation has been detected.  Output the final norm bound and other info.

t1 = mpreald (1.d300, nwds)
t2 = mpreald (0.d0, nwds)

!   Select the relation corresponding to the smallest y entry and compute norm.

do j = 1, n
  if (abs (y(j)) < t1) then
    j1 = j
    t1 = abs (y(j))
  endif
  t2 = max (t2, abs (y(j)))
enddo

do i = 1, n
  r(i) = b(j1,i)
enddo

!   The norm and norm bound calculation here are performed in medium precision.

tm0 = second ()
w1 = mprealdm (0.d0, nwdsm)

do i = 1, n
  w1 = w1 + mprealm (r(i), nwdsm) ** 2
enddo

w1 = sqrt (w1)

do j = 1, n
  do i = 1, n
    wh(i,j) = mprealm (h(i,j), nwdsm)
  enddo
enddo
  
w2 = bound (n, nwdsm, wh)
wn = max (wn, w2)
tm1 = second ()
times(3) = times(3) + (tm1 - tm0)

!   Output the final norm bound and other information.

if (idb >= 1) then
  call mpdecmd (t1, d1, n1)
  call mpdecmd (t2, d2, n2)
  call mpdecmd (w1, d3, n3)
  call mpdecmd (wn, d4, n4)
  write (6, 20) it, d1, n1, d2, n2, d4, n4
20 format ('Iteration',i8,2x,'Relation detected'/ &
  'Min, max of y =',f11.6,'e',i6,f11.6,'e',i6/'Max. bound =',f11.6,'e',i6)
  write (6, 21) j1, d3, n3, d1, n1
21 format ('Index of relation =',i4,3x,'Norm =',f11.5,'e',i5,3x, &
  'Residual =',f11.6,'e',i6)
endif

!   If run was successful, set iq = 1; otherwise output message.

if (t1 == 0.d0) n1 = nep
if (n3 <= nrb .and. n2 - n1 >= ndr) then
  iq = 1
else
  if (idb >= 2) write (6, 22)
22 format ('Relation is too large or insufficient dynamic range in y at detection.')
endif

!+  Final section.

160 continue

!   Successful or not, output CPU run times and return.

if (idb >= 2) write (6, 23) times, sum (times)
23 format ('CPU run times:'/(7f12.2))

return
end


!------------------------------

!   First-level subroutines.

subroutine saverest (iu, n, it, izm, izmm, izd, times, x, y, b, h)

!   This saves restart data to the file.

use mpmodule
implicit none
integer iu, n, it, izm, izmm, izd
real (mprknd) times(6)
type (mp_real) b(n,n), h(n,n), x(n), y(n)

rewind (iu)
write (iu) it, izm, izmm, izd
write (iu) times
write (iu) x
write (iu) y
write (iu) b
write (iu) h
flush (iu)
rewind (iu)

return
end

subroutine setprec (nwdsm, n, nsq, wa, wb, wh, wsyq, wy)

!   This sets the working precision of the arrays wa, wb, wh, wsyq and wy
!   to nwdsm words.

use mpmodule
implicit none
integer i, j, n, nsq, nwdsm
type (mp_realm) wa(n,n), wb(n,n), wh(n,n), wsyq(n,nsq), wy(n)

do j = 1, n
  do i = 1, n
    wa(i,j) = mprealm (wa(i,j), nwdsm)
    wb(i,j) = mprealm (wb(i,j), nwdsm)
    wh(i,j) = mprealm (wh(i,j), nwdsm)
  enddo
enddo

do j = 1, nsq
  do i = 1, n
    wsyq(i,j) = mprealm (wsyq(i,j), nwdsm)
  enddo
enddo

do i = 1, n
  wy(i) = mprealm (wy(i), nwdsm)
enddo

return
end

function dynrange (n, nwdsm, y)

!   This returns the dynamic range of y, i.e., ratio of min|y_k| to max|y_k|,
!   normally done using MPM precision. Here y is of type mp_real.

use mpmodule
implicit none
integer i, n, nwdsm
type (mp_realm) dynrange, t1, t2, t3
type (mp_real) y(n)

t1 = mprealdm (1.d300, nwdsm)
t2 = mprealdm (0.d0, nwdsm)

!   Find the min and max absolute value in the y vector.

do i = 1, n
  t3 = abs (mprealm (y(i), nwdsm))
  t1 = min (t1, t3)
  t2 = max (t2, t3)
enddo

dynrange = t1 / t2
return
end

function dynrangem (n, nwdsm, wy)

!   This returns the dynamic range of wy, i.e., ratio of min|wy_k| to max|wy_k|,
!   normally done using MPM precision. Here wy is of type mp_realm.

use mpmodule
implicit none
integer i, n, nwdsm
type (mp_realm) dynrangem, t1, t2, t3, wy(n) 

t1 = mprealdm (1.d300, nwdsm)
t2 = mprealdm (0.d0, nwdsm)

!   Find the min and max absolute value in the wy vector.

do i = 1, n
  t3 = abs (mprealm (wy(i), nwdsm))
  t1 = min (t1, t3)
  t2 = max (t2, t3)
enddo

dynrangem = t1 / t2
return
end

subroutine initdp (idb, n, nsq, nwdsm, da, db, dh, dy, dsyq, wh, wy)

!   This initializes the DP arrays from the MPM arrays.
!   This is performed in medium precision.
!   Input:  idb, n, nsq, nwdsm, wh, wy.
!   Output: da, db, dh, dy, dsyq.

use mpmodule
implicit none
integer i, idb, j, n, nsq, nwdsm
real (mprknd) da(n,n), db(n,n), dh(n,n), dy(n), dsyq(n,nsq)
type (mp_realm) wh(n,n), wy(n), t1, t2

t2 = mprealdm (0.d0, nwdsm)

!   Find the max absolute value in the wy vector.

do i = 1, n
  t2 = max (t2, abs (wy(i)))
enddo

!   Set dy to be the scaled wy vector.

t1 = 1.d0 / t2

do i = 1, n
  dy(i) = t1 * wy(i)
enddo

!   Find the maximum absolute value of the wh matrix diagonals.

t2 = mprealdm (0.d0, nwdsm)

do j = 1, n - 1
  t2 = max (t2, abs (wh(j,j)))
enddo

!   Set dh to be the scaled wh matrix.

t1 = 1.d0 / t2

do j = 1, n - 1
  do i = 1, n
    dh(i,j) = t1 * wh(i,j)
  enddo
enddo

!   Set da and db to the identity.

do j = 1, n
  do i = 1, n
    da(i,j) = 0.d0
    db(i,j) = 0.d0
  enddo

  da(j,j) = 1.d0
  db(j,j) = 1.d0
enddo

!   Zero the dsyq array.

do j = 1, nsq
  do i = 1, n
    dsyq(i,j) = 0.d0
  enddo
enddo

if (idb >= 4) then
  write (6, 2)
2 format ('initdp: Scaled dy vector:')
  call matoutdp (1, n, dy)
  write (6, 3)
3 format ('initdp: Scaled dh matrix:')
  call matoutdp (n, n - 1, dh)
endif

return
end

subroutine initmp (idb, n, nwds, b, h, x, y)

!   Initializes MP arrays at the beginning.
!   This is performed in full precision.
!   Input: idb, n, nwds.
!   Output: b, h, x, y.

use mpmodule
implicit none
integer i, idb, j, n, nwds
integer ix(n)
real (mprknd) dx(n)
type (mp_real) b(n,n), h(n,n), s(n), x(n), y(n), t1

if (idb >= 4) then
  write (6, 1)
1 format ('initmp: Input x vector:')
  call matoutmd (1, n, ix, dx, x)
endif

!   Set b to the identity matrix.

do j = 1, n
  do i = 1, n
    b(i,j) = mpreald (0.d0, nwds)
  enddo

  b(j,j) = mpreald (1.d0, nwds)
enddo

t1 = mpreald (0.d0, nwds)

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
    h(i,j) = mpreald (0.d0, nwds)
  enddo

  h(j,j) = s(j+1) / s(j)
  t1 = y(j) / (s(j) * s(j+1))

  do i = j + 1, n
    h(i,j) = - y(i) * t1
  enddo
enddo

do i = 1, n
  h(i,n) = mpreald (0.d0, nwds)
enddo

if (idb >= 4) then
  write (6, 2)
2 format ('initmp: Initial y vector:')
  call matoutmd (1, n, ix, dx, y)
  write (6, 3)
3 format ('initmp: Initial h matrix:')
  call matoutmd (n, n - 1, ix, dx, h)
endif

return
end

subroutine initmpm (idb, n, nsq, nwdsm, wa, wb, wh, wy, wsyq, h, y)

!   This initializes the MPM arrays from the MP arrays.
!   This is performed in medium precision.
!   Input: idb, n, nsq, nwdsm, h, y.
!   Output: wa, wb, wh, wy, wsyq.

use mpmodule
implicit none
integer i, idb, j, n, nsq, nwdsm
integer ix(n)
real (mprknd) dx(n)
type (mp_realm) wa(n,n), wb(n,n), wh(n,n), wy(n), wsyq(n,nsq), t1, t2, t3
type (mp_real) h(n,n), y(n)

t2 = mprealdm (0.d0, nwdsm)

!   Find the max absolute value in the y vector.

do i = 1, n
  t3 = abs (mprealm (y(i), nwdsm))
  t2 = max (t2, t3)
enddo

!   Set wy to the scaled y vector.

t1 = 1.d0 / t2

do i = 1, n
  wy(i) = t1 * mprealm (y(i), nwdsm)
enddo

!   Set wh to h.

do j = 1, n
  do i = 1, n
    wh(i,j) = mprealm (h(i,j), nwdsm)
  enddo
enddo

!   Set wa and wb to the identity.

do j = 1, n
  do i = 1, n
    wa(i,j) = mprealdm (0.d0, nwdsm)
    wb(i,j) = mprealdm (0.d0, nwdsm)
  enddo

  wa(j,j) = mprealdm (1.d0, nwdsm)
  wb(j,j) = mprealdm (1.d0, nwdsm)
enddo

!   Zero the wsyq array.

do j = 1, nsq
  do i = 1, n
    wsyq(i,j) = mprealdm (0.d0, nwdsm)
  enddo
enddo

if (idb >= 4) then
  write (6, 3)
3 format ('initmpm: wy vector:')
  call matoutmmd (1, n, ix, dx, wy)
  write (6, 4)
4 format ('initmpm: Factored wh matrix:')
  call matoutmmd (n, n - 1, ix, dx, wh)
endif

100 continue

return
end

subroutine iterdp (idb, it, n, nsq, da, db, dh, dsyq, dy, imq, izd)


!   This performs one iteration of the PSLQ algorithm using DP arithmetic.
!   Input: idb, it, n, nsq, da, db, dh, dsyq, dy.
!   Output: da, db, dh, dsyq, dy, imq, izd.

!   NOTE: Parameter tmx2 = 2^52, not 2^53, so as to ensure that values > 2^53
!   never arise, even as intermediate values, in the update loop below.

use mpmodule
implicit none
integer i, idb, ii, ij, im, im1, imq, it, izd, j, j1, j2, k, mpr, mq, n, nsq
real (mprknd) deps, gam, t1, t2, t3, t4, tmx1, tmx2
integer ip(n), ir(n), is(n)
real (mprknd) da(n,n), db(n,n), dh(n,n), dq(n), dsyq(n,nsq), &
  dt(n,n), dy(n)
parameter (tmx1 = 1.d13, tmx2 = 2.d0**52, deps = 1.d-14)

izd = 0
mpr = nint (0.4d0 * n)
gam = sqrt (4.d0 / 3.d0)

!   Compute dq vector = {gam^i * |dh(i,i)|}, then sort in ascending order.

do i = 1, n - 1
  dq(i) = gam ** i * abs (dh(i,i))
enddo

call qsortdp (n - 1, dq, ip)

!   Select up to mpr disjoint pairs of indices (m,m+1), where m is an index
!   from the list of the largest dq(i).

do i = 1, n
  is(i) = 0
enddo

!  If the imq flag is set, perform this iteration with mq = 1.

if (imq == 0) then
  mq = mpr
else
  mq = 1
  imq = 0
endif
ii = n

do i = 1, mq
100 continue
  ii = ii - 1
  if (ii == 0) then
    mq = i - 1
    goto 110
  endif
  j1 = ip(ii)
  j2 = j1 + 1
  if (is(j1) .ne. 0 .or. is(j2) .ne. 0) goto 100
  ir(i) = j1
  is(j1) = 1
  is(j2) = 1
enddo

110 continue

!   Exchange the pairs of entries of dy, and rows of da, db and dh.

do j = 1, mq
  im = ir(j)
  im1 = im + 1
  t1 = dy(im)
  dy(im) = dy(im1)
  dy(im1) = t1

  do i = 1, n
    t1 = da(im,i)
    da(im,i) = da(im1,i)
    da(im1,i) = t1
    t1 = db(im,i)
    db(im,i) = db(im1,i)
    db(im1,i) = t1
  enddo

  do i = 1, n - 1
    t1 = dh(im,i)
    dh(im,i) = dh(im1,i)
    dh(im1,i) = t1
  enddo
enddo

!   Eliminate the "corners" produced by the above permutation in dh.

do j = 1, mq
  im = ir(j)
  im1 = im + 1
  if (im <= n - 2) then
    t1 = dh(im,im)
    t2 = dh(im,im1)
    t3 = sqrt (t1 ** 2 + t2 ** 2)
    t1 = t1 / t3
    t2 = t2 / t3

    do i = im, n
      t3 = dh(i,im)
      t4 = dh(i,im1)
      dh(i,im) = t1 * t3 + t2 * t4
      dh(i,im1) = - t2 * t3 + t1 * t4
    enddo
  endif
enddo

!   Perform reduction on dh, using the diagonal scheme.  Multipliers are
!   saved in the dt array.

do i = 2, n
  do j = 1, n - i + 1
    ij = i + j - 1

    do k = j + 1, ij - 1
      dh(ij,j) = dh(ij,j) - dt(ij,k) * dh(k,j)
    enddo

    dt(ij,j) = anint (dh(ij,j) / dh(j,j))
    dh(ij,j) = dh(ij,j) - dt(ij,j) * dh(j,j)
  enddo
enddo

!   Update dy, using the dt array.  Find min absolute value of dy.

t1 = abs (dy(n))

do j = 1, n - 1
  do i = j + 1, n
    dy(j) = dy(j) + dt(i,j) * dy(i)
  enddo

  t1 = min (t1, abs (dy(j)))
enddo

!   Update da and db, using the dt array.  Find the max absolute value of
!   da and db entries as they are calculated (not merely at the end).

t2 = 0.d0

do k = 1, n
  dq(k) = 0.d0

  do j = 1, n - 1
    do i = j + 1, n
      da(i,k) = da(i,k) - dt(i,j) * da(j,k)
      db(j,k) = db(j,k) + dt(i,j) * db(i,k)
      dq(k) = max (dq(k), abs (da(i,k)), abs (db(j,k)))
    enddo
  enddo
enddo

do k = 1, n
  t2 = max (t2, dq(k))
enddo

if (t1 <= deps) then
  if (idb >= 3) write (6, 1) it, t1
1 format ('Iteration',i8,2x,'iterdp: Small value in dy =',1pd15.6)
  izd = 1
endif

if (t2 > tmx1 .and. t2 <= tmx2) then
  if (idb >= 3) write (6, 2) it, t2
2 format ('Iteration',i8,2x,'iterdp: Large value in da or db =',1pd15.6)
  izd = 1
elseif (t2 > tmx2) then
  if (idb >= 2) write (6, 3) it, t2
3 format ('Iteration',i8,2x,'iterdp: Very large value in da or db =',1pd15.6)
  izd = 2
  return
endif

!   Compare the dy vector with those of recent iterations.  If a duplicate is
!   found, then the next iteration must be performed with mq = 1.

do j = 1, nsq
  t1 = 0.d0

  do i = 1, n
    t1 = max (t1, abs (dy(i) - dsyq(i,j)))
  enddo

  if (t1 <= deps) then
    if (idb >= 2) write (6, 4) it, j
4   format ('Iteration',i8,2x,'iterdp: Duplicate found, j =',i6)
    imq = 1
    goto 120
  endif
enddo

!   Place the vector dy in the table dsyq.

120   continue
k = 1 + mod (it, nsq)

do i = 1, n
  dsyq(i,k) = dy(i)
enddo

if (idb >= 4) then
  write (6, 5)
5 format ('iterdp: Updated dy:')
  call matoutdp (1, n, dy)
  write (6, 6)
6 format ('iterdp: Updated da matrix:')
  call matoutdp (n, n, da)
  write (6, 7)
7 format ('iterdp: Updated db matrix:')
  call matoutdp (n, n, db)
  write (6, 8)
8 format ('iterdp: Updated dh matrix:')
  call matoutdp (n, n - 1, dh)
endif

return
end

subroutine itermpm (idb, it, n, nsq, nwdsm, epsm, wa, wb, wh, wsyq, wy, imq, izmm)

!   This performs one iteration of the PSLQM algorithm using MPM arithmetic.
!   This is performed in medium precision.
!   Input: idb, it, n, nsq, nwdsm, epsm, imq.
!   Output: wa, wb, wh, wsyq, wy, imq, izmm.

use mpmodule
implicit none
integer i, idb, ii, ij, im, im1, imq, it, izmm, j, j1, j2, k, mpr, mq, n, n1, &
  nsq, nwdsm
real (mprknd) d1
type (mp_realm) gam, t1, t2, t3, t4, epsm, tmx1, tmx2
integer ip(n), ir(n), is(n)
real (mprknd) dx(n)
type (mp_realm) wa(n,n), wb(n,n), wh(n,n), wq(n), wsyq(n,nsq), wt(n,n), wy(n)

tmx1 = 1.d0 / mprealm (epsm, nwdsm) / mprealdm (1d20, nwdsm)
tmx2 = mprealdm (2.d0, nwdsm) ** (nwdsm * mpnbt)
izmm = 0
mpr = nint (0.4d0 * n)
gam = sqrt (mprealdm (4.d0, nwdsm) / mprealdm (3.d0, nwdsm))

!   Compute wq vector = {gam^i * |wh(i,i)|}, then sort in ascending order.

do i = 1, n - 1
  wq(i) = gam ** i * mprealm (abs (wh(i,i)), nwdsm)
enddo

call qsortmpm (n - 1, wq, ip)

!   Select up to mpr disjoint pairs of indices (m,m+1), where m is an index
!   from the list of the largest wq(i).

do i = 1, n
  is(i) = 0
enddo

if (imq == 0) then
  mq = mpr
else
  mq = 1
  imq = 0
endif
ii = n

do i = 1, mq
100 continue
  ii = ii - 1
  if (ii == 0) then
    mq = i - 1
    goto 110
  endif
  j1 = ip(ii)
  j2 = j1 + 1
  if (is(j1) .ne. 0 .or. is(j2) .ne. 0) goto 100
  ir(i) = j1
  is(j1) = 1
  is(j2) = 1
enddo

110 continue

!   Exchange the pairs of entries of wy, and rows of wa, wb and wh.

do j = 1, mq
  im = ir(j)
  im1 = im + 1
  t1 = wy(im)
  wy(im) = wy(im1)
  wy(im1) = t1

  do i = 1, n
    t1 = wa(im,i)
    wa(im,i) = wa(im1,i)
    wa(im1,i) = t1
    t1 = wb(im,i)
    wb(im,i) = wb(im1,i)
    wb(im1,i) = t1
  enddo

  do i = 1, n - 1
    t1 = wh(im,i)
    wh(im,i) = wh(im1,i)
    wh(im1,i) = t1
  enddo
enddo

!   Eliminate the "corners" produced by the above permutation in wh.

do j = 1, mq
  im = ir(j)
  im1 = im + 1
  if (im <= n - 2) then
    t1 = wh(im,im)
    t2 = wh(im,im1)
    t3 = sqrt (t1 ** 2 + t2 ** 2)
    t1 = t1 / t3
    t2 = t2 / t3

    do i = im, n
      t3 = wh(i,im)
      t4 = wh(i,im1)
      wh(i,im) = t1 * t3 + t2 * t4
      wh(i,im1) = - t2 * t3 + t1 * t4
    enddo
  endif
enddo

!   Perform reduction on wh, using the diagonal scheme.  Multipliers are
!   saved in the wt array.

do i = 2, n
  do j = 1, n - i + 1
    ij = i + j - 1

    do k = j + 1, ij - 1
      wh(ij,j) = wh(ij,j) - wt(ij,k) * wh(k,j)
    enddo

    wt(ij,j) = anint (wh(ij,j) / wh(j,j))
    wh(ij,j) = wh(ij,j) - wt(ij,j) * wh(j,j)
  enddo
enddo

!   Update wy, using the wt array.  Find min absolute value of wy.

t1 = abs (wy(n))

do j = 1, n - 1
  do i = j + 1, n
    wy(j) = wy(j) + wt(i,j) * wy(i)
  enddo

  t1 = min (t1, abs (wy(j)))
enddo

!   Update wa and wb, using the wt array.  Find the max absolute value of
!   wa and wb entries as they are calculated (not merely at the end).

t2 = mprealdm (0.d0, nwdsm)

do k = 1, n
  wq(k) = mprealdm (0.d0, nwdsm)

  do j = 1, n - 1
    do i = j + 1, n
      wa(i,k) = wa(i,k) - wt(i,j) * wa(j,k)
      wb(j,k) = wb(j,k) + wt(i,j) * wb(i,k)
      wq(k) = max (wq(k), abs (wa(i,k)), abs (wb(j,k)))
    enddo
  enddo
enddo

do k = 1, n
  t2 = max (t2, wq(k))
enddo

if (t1 <= epsm) then
  if (idb >= 2) then
    call mpdecmd (t1, d1, n1) 
    write (6, 1) it, d1, n1
1   format ('Iteration',i8,2x,'itermpm: Small value in wy =',f11.6,'e',i6)
  endif
  izmm = 1
endif

if (t2 > tmx1 .and. t2 <= tmx2) then
  if (idb >= 2) then
    call mpdecmd (t2, d1, n1)
    write (6, 2) it, d1, n1
2   format ('Iteration',i8,2x,'itermpm: Large value in wa or wb =', &
      f11.6,'e',i6)
  endif
  izmm = 1
elseif (t2 > tmx2) then
  if (idb >= 1) then
    call mpdecmd (t2, d1, n1)
    write (6, 3) it, d1, n1
3   format ('Iteration',i8,2x,'itermpm: Very large value in wa or wb =', &
      f11.6,'e',i6/'Run aborted.')
  endif
  izmm = 2
  goto 200
endif

!   Compare the wy vector with those of recent iterations.  If a duplicate is
!   found, then the next iteration must be performed with mq = 1.

do j = 1, nsq
  t1 = mprealdm (0.d0, nwdsm)

  do i = 1, n
    t1 = max (t1, abs (wy(i) - wsyq(i,j)))
  enddo

  if (t1 <= epsm) then
    if (idb >= 2) write (6, 4) it, j
4   format ('Iteration',i8,2x,'itermpm: Duplicate found, j =',i6)
    imq = 1
    goto 120
  endif
enddo

!   Place the vector wy in the table wsyq.

120   continue
k = 1 + mod (it, nsq)

do i = 1, n
  wsyq(i,k) = wy(i)
enddo

if (idb >= 4) then
  write (6, 5)
5 format ('itermpm: Updated wy:')
  call matoutmmd (1, n, ip, dx, wy)
  write (6, 6)
6 format ('itermpm: Updated wa matrix:')
  call matoutmmd (n, n, ip, dx, wa)
  write (6, 7)
7 format ('itermpm: Updated wb matrix:')
  call matoutmmd (n, n, ip, dx, wb)
  write (6, 8)
8 format ('itermpm: Updated wh matrix:')
  call matoutmmd (n, n - 1, ip, dx, wh)
endif

200 continue

return
end

subroutine lqdp (n, m, dh)

!   This performs an LQ decomposition on the DP matrix dh.  It is a simplified
!   and transposed adaptation of the subroutine dqrdc from Linpack.
!   Input: n, m, dh.
!   Output: dh.

use mpmodule
implicit none
integer i, j, l, lup, m, ml, n
real (mprknd) dh(n,m), nrmxl, one, t, zero

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
  if (dh(l,l) .ne. zero) nrmxl = sign (nrmxl, dh(l,l))
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
end

subroutine lqmpm (n, m, nwdsm, h)

!   This performs an LQ decomposition on the MPM matrix h.  It is a simplified
!   and transposed adaptation of the subroutine dqrdc from Linpack.
!   This is normally performed in variable medium precision.
!   Input: n, m, nwdsm, h.
!   Output: h.

use mpmodule
implicit none
integer i, j, l, lup, m, ml, n, nwdsm
type (mp_realm) h(n,m), nrmxl, one, t, zero

zero = mprealdm (0.d0, nwdsm)
one = mprealdm (1.d0, nwdsm)
lup = min (m,n)

!   Perform the householder reduction of h.

do l = 1, lup
  if (l == m) go to 280

!   Compute the householder transformation for column l.

  ml = m - l
  t = zero

  do i = 0, ml
    t = t + h(l,l+i) ** 2
  enddo

  nrmxl = sqrt (t)
  if (nrmxl == zero) go to 270
  if (h(l,l) .ne. zero) nrmxl = sign (nrmxl, h(l,l))
  t = one / nrmxl

  do i = 0, ml
    h(l,l+i) = t * h(l,l+i)
  enddo

  h(l,l) = one + h(l,l)

!   Apply the transformation to the remaining columns, updating the norms.

  do j = l + 1, n
    t = zero

    do i = 0, ml
      t = t + h(l,l+i) * h(j,l+i)
    enddo

    t = - t / h(l,l)

    do i = 0, ml
      h(j,l+i) = h(j,l+i) + t * h(l,l+i)
    enddo
  enddo

!   Save the transformation.

  h(l,l) = - nrmxl
270 continue
280 continue
enddo

!   Zero h above the diagonal.

do j = 1, m
  do i = 1, j - 1
    h(i,j) = zero
  enddo
enddo

return
end

subroutine savedp (n, da, db, dh, dy, dsa, dsb, dsh, dsy)

!   This saves the arrays dy, da, db, dh in case DP iterations must be aborted.
!   A call to the same routine, with (da,db,dh,dy) and (dsa,dsb,dsh,dsy)
!   exchanged, serves to restore these arrays.
!   Input: n, da, db, dh, dy.
!   Output: dsa, dsb, dsh, dsy.

use mpmodule
implicit none
integer i, j, n
real (mprknd) da(n,n), db(n,n), dh(n,n), dy(n), dsa(n,n), dsb(n,n), &
  dsh(n,n), dsy(n)

do i = 1, n
  dsy(i) = dy(i)
enddo

do j = 1, n
  do i = 1, n
    dsa(i,j) = da(i,j)
    dsb(i,j) = db(i,j)
    dsh(i,j) = dh(i,j)
  enddo
enddo

return
end

subroutine updtmp (idb, it, n, nwds, wa, wb, eps, b, h, y, izm)

!   This update the MP arrays from the MPM arrays.
!   This is performed in full precision.
!   Input: idb, it, n, nwds, wa, wb, eps, b, h, y.
!   Output: b, h, y, izm.

use mpmodule
implicit none 
integer i, i1, idb, it, izm, n, n1, n2, nwds
real (mprknd) d1, d2
type (mp_real) eps, t1, t2
integer ix(n)
real (mprknd) dx(n)
type (mp_realm) wa(n,n), wb(n,n)
type (mp_real) b(n,n), h(n,n), y(n)

izm = 0

!   Update y with wb.

call mxm (n, 1, nwds, wb, y)

i1 = 0
t1 = mpreald (1.d300, nwds)
t2 = mpreald (0.d0, nwds)

do i = 1, n
  if (abs (y(i)) < t1) then
    i1 = i
    t1 = abs (y(i))
  endif
  t2 = max (t2, abs (y(i)))
enddo

if (idb >= 2) then
  call mpdecmd (t1, d1, n1)
  call mpdecmd (t2, d2, n2)
  write (6, 1) it, d1, n1, d2, n2
1  format ('Iteration',i8,2x,'updtmp: Min, max of y =',f11.6,'e',i6, &
    f11.6,'e',i6)
endif

!   Update b with wb.

call mxm (n, n, nwds, wb, b)

!   Update h with wa.

call mxm (n, n - 1, nwds, wa, h)

!   Find the largest entry of b in the same row as the smallest y.

t2 = mpreald (0.d0, nwds)

do i = 1, n
  t2 = max (t2, abs (b(i1,i)))
enddo

if (t1 <= t2 * eps) then
  if (idb >= 2) then
    call mpdecmd (t1, d1, n1) 
    write (6, 2) it, d1, n1
2   format ('Iteration',i8,2x,'updtmp: Small value in y =',f11.6,'e',i6)
  endif
  if (t1 <= t2 * eps) then
    izm = 1
  else
    if (idb >= 1) write (6, 3) it
3   format ('Iteration',i8,2x,'updtmp: Precision exhausted.')
    izm = 2
  endif
endif

if (idb >= 4) then
  write (6, 4)
4 format ('updtmp: Updated y:')
  call matoutmd (1, n, ix, dx, y)
  write (6, 5)
5 format ('updtmp: Updated b matrix:')
  call matoutmd (n, n, ix, dx, b)
  write (6, 6)
6 format ('updtmp: Updated h matrix:')
  call matoutmd (n, n - 1, ix, dx, h)
endif

return
end

subroutine updtmpm (idb, it, n, nwdsm, epsm, da, db, dreps, wa, wb, wh, wy, izmm)

!   This updates the MPM arrays from the DP arrays.
!   This is performed in medium precision.
!   Input: idb, it, n, nwdsm, epsm, da, db, dreps, wa, wb, wh, wy.
!   Output: wa, wb, wh, wy, izmm.

use mpmodule
implicit none 
integer i, idb, it, izmm, j, n, n1, n2, nwdsm
real (mprknd) d1, d2
type (mp_realm) t1, t2, epsm, tmx1, tmx2
integer ix(n)
real (mprknd) da(n,n), db(n,n), dx(n), dreps
type (mp_realm) wa(n,n), wb(n,n), wh(n,n), wy(n), w1

tmx1 = 1.d0 / epsm
tmx2 = mprealdm (2.d0, nwdsm) ** (nwdsm * mpnbt)
izmm = 0

!   Update wy with db.

call mxmdm (n, 1, nwdsm, db, wy)

!   Find min/max ratio of wy.

t1 = mprealdm (1.d300, nwdsm)
t2 = mprealdm (0.d0, nwdsm)

do i = 1, n
  t1 = min (t1, abs (wy(i)))
  t2 = max (t2, abs (wy(i)))
enddo

w1 = t1 / t2

if (idb >= 3) then
  call mpdecmd (t1, d1, n1)
  call mpdecmd (t2, d2, n2)
  write (6, 1) it, d1, n1, d2, n2
1 format ('Iteration',i8,2x,'updtmpm: Min, max of wy =',f11.6,'e',i6, &
     f11.6,'e',i6)
endif
if (t1 <= epsm) then
  if (idb >= 2) then
    call mpdecmd (t1, d1, n1)
    write (6, 2) it, d1, n1
2   format ('Iteration',i8,2x,'updtmpm: Small value in wy =',f11.6,'e',i6)
  endif
  izmm = 1
endif

!   Update wa with da.

call mxmdm (n, n, nwdsm, da, wa)

!   Update wb with db.

call mxmdm (n, n, nwdsm, db, wb)
t2 = mprealdm (0.d0, nwdsm)

do j = 1, n
  do i = 1, n
    t2 = max (t2, abs (wa(i,j)), abs (wb(i,j)))
  enddo
enddo

if (t2 > tmx1 .and. t2 <= tmx2) then
  if (idb >= 2) then
    call mpdecmd (t2, d1, n1)
    write (6, 3) it, d1, n1
3   format ('Iteration',i8,2x,'updtmpm: Large value in wa or wb =', &
      f11.6,'e',i6)
  endif
  izmm = 1
elseif (t2 > tmx2) then
  if (idb >= 1) then
    call mpdecmd (t2, d1, n1)
    write (6, 4) it, d1, n1
4   format ('updtmpm: Very large value in wa or wb =',f11.6,'e',i6/ &
      'Run aborted.')
  endif
  izmm = 2
  goto 100
endif

if (w1 < mprealdm (dreps, nwdsm)) then
  if (idb >= 2) then
    call mpdecmd (w1, d1, n1)
    write (6, 9) it, d1, n1
9   format ('Iteration',i8,2x,'updtmpm: Small min/max ratio of wy =',f11.6,'e',i6)
  endif
  izmm = 1
endif

!   Update wh with da.

call mxmdm (n, n - 1, nwdsm, da, wh)

if (idb >= 4) then
  write (6, 5)
5 format ('updtmpm: Updated wy:')
  call matoutmmd (1, n, ix, dx, wy)
  write (6, 6)
6 format ('updtmpm: Updated wa matrix:')
  call matoutmmd (n, n, ix, dx, wa)
  write (6, 7)
7 format ('updtmpm: Updated wb matrix:')
  call matoutmmd (n, n, ix, dx, wb)
  write (6, 8)
8 format ('updtmpm: Updated wh matrix:')
  call matoutmmd (n, n - 1, ix, dx, wh)
endif

100 continue

return
end

!------------------------------

!   Second- and third-level subroutines.

function bound (n, nwdsm, wh)

!   This computes the norm bound, normally done using MPM arithmetic.

use mpmodule
implicit none
integer i, n, nwdsm
type (mp_realm) wh(n,n), bound, t1

call lqmpm (n, n - 1, nwdsm, wh)
t1 = mprealdm (0.d0, nwdsm)

do i = 1, n - 1
  t1 = max (t1, abs (wh(i,i)))
enddo

bound = 1.d0 / t1
return
end

function dplog10 (a)

!   For input MPM value a, this routine returns a DP approximation to log10 (a).

use mpmodule
implicit none
integer ia
real (mprknd) da, dplog10
type (mp_realm) a

call mpdecmd (a, da, ia)
if (da == 0.d0) then
  dplog10 = -999999.d0
else
  dplog10 = log10 (abs (da)) + ia
endif

100 continue
return
end

subroutine matoutdp (n1, n2, a)

!   This outputs the DP matrix a.

use mpmodule
implicit none
integer i, j, n1, n2
real (mprknd) a(n1,n2)

do i = 1, n1
  write (6, 1) i
1 format ('Row',i3)
  write (6, 2) (a(i,j), j = 1, n2)
2 format (1p,5d15.5)
enddo

return
end

subroutine matoutmd (n1, n2, ix, dx, a)

!   This outputs the MP matrix a as a DP matrix.

use mpmodule
implicit none
integer i, j, n1, n2
integer ix(n2)
real (mprknd) dx(n2)
type (mp_real) a(n1,n2)

do i = 1, n1
  write (6, 1) i
1 format ('Row',i3)

  do j = 1, n2
    call mpdecmd (a(i,j), dx(j), ix(j))
  enddo

  write (6, 2) (dx(j), ix(j), j = 1, n2)
2 format (4(f13.8,'e',i6))
enddo

return
end

subroutine matoutmmd (n1, n2, ix, dx, a)

!   This outputs the MPM matrix a as a DP matrix.

use mpmodule
implicit none
integer i, j, n1, n2
integer ix(n2)
real (mprknd) dx(n2)
type (mp_realm) a(n1,n2)

do i = 1, n1
  write (6, 1) i
1 format ('Row',i3)

  do j = 1, n2
    call mpdecmd (a(i,j), dx(j), ix(j))
  enddo

  write (6, 2) (dx(j), ix(j), j = 1, n2)
2 format (4(f13.8,'e',i6))
enddo

return
end

subroutine matoutmp (n1, n2, a)

!   This outputs the MP or MPM matrix a.  It may be used in place of calls to
!   matoutmd in the code above if greater accuracy is desired in debug output.

use mpmodule
implicit none
integer i, j, n1, n2
type (mp_real) a(n1,n2)

do i = 1, n1
  write (6, 1) i
1 format ('Row',i3)

  do j = 1, n2
    call mpwrite (6, 80, 60, a(i,j))
  enddo
enddo

return
end

subroutine mxmdm (n1, n2, nwdsm, a, b)

!   This multiplies the DP square matrix a by the MPM matrix b, and the result
!   is placed in b.  n1, n2 are the matrix dimensions as indicated below.
!   This is performed in medium precision.
!   Input: n1, n2, nwdsm, a, b.
!   Output: b.

use mpmodule
implicit none
integer i, j, k, n1, n2, nwdsm
real (mprknd) a(n1,n1)
type (mp_realm) b(n1,n2), c(n1,n2)

do j = 1, n2
  do i = 1, n1
    c(i,j) = mprealdm (0.d0, nwdsm)
    
    do k = 1, n1
      c(i,j) = c(i,j) + mpprod (b(k,j), a(i,k))
    enddo
  enddo
enddo

do j = 1, n2
  do i = 1, n1
    b(i,j) = c(i,j)
  enddo
enddo

return
end

subroutine mxm (n1, n2, nwds, a, b)

!   This multiplies the MPM square matrix a by the MP matrix b, and the
!   result is returned in b.  n1, n2 are the matrix dimensions as indicated below.
!   This is performed in full precision.
!   Input: n1, n2, nwds, a, b.
!   Output: b.

use mpmodule
implicit none
integer i, j, k, n1, n2, nwds
type (mp_realm) a(n1,n1)
type (mp_real) b(n1,n2), c(n1,n2)

do j = 1, n2
  do i = 1, n1
    c(i,j) = mpreal (0.d0, nwds)

    do k = 1, n1
      c(i,j) = c(i,j) + mpreal (a(i,k), nwds) * b(k,j)
    enddo
  enddo
enddo

do j = 1, n2
  do i = 1, n1
    b(i,j) = c(i,j)
  enddo
enddo

return
end

subroutine qsortdp (n, a, ip)

!   This routine sorts the entries of the N-long DP vector A into ascending
!   order using the quicksort algorithm.  The permutation vector that would
!   sort the vector is returned in IP.
!   Input: n, a.
!   Output: ip.

use mpmodule
implicit none
integer i, iq, it, j, jq, jz, k, l, n
real (mprknd) a(n), s0
integer ip(n), ik(50), jk(50)

do i = 1, n
  ip(i) = i
enddo

if (n == 1) return

k = 1
ik(1) = 1
jk(1) = n

130 continue

i = ik(k)
j = jk(k)
iq = i
jq = j
it = (i + j + 1) / 2
l = ip(j)
ip(j) = ip(it)
ip(it) = l
s0 = a(ip(j))
j = j - 1

140 continue

do l = i, j
  if (s0 < a(ip(l))) goto 160
enddo

i = j
goto 190

160 i = l

do l = j, i, -1
  if (s0 > a(ip(l))) goto 180
enddo

j = i
goto 190

180 continue

j = l
if (i >= j) goto 190
l = ip(i)
ip(i) = ip(j)
ip(j) = l
goto 140

190 continue

if (s0 >= a(ip(i))) goto 200
l = ip(jq)
ip(jq) = ip(i)
ip(i) = l

200 continue

k = k - 1
jz = 0
if (j == iq) goto 210
k = k + 1
jk(k) = j
jz = 1

210 continue

i = i + 1
if (i == jq) goto 220
k = k + 1
ik(k) = i
jk(k) = jq
if (jz == 0) goto 220
if (j - iq >= jq - i) goto 220
ik(k-1) = i
jk(k-1) = jq
ik(k) = iq
jk(k) = j

220 continue

if (k > 0) goto 130

return
end

subroutine qsortmpm (n, a, ip)

!   This routine sorts the entries of the N-long MPM vector A into ascending
!   order using the quicksort algorithm.  The permutation vector that would
!   sort the vector is returned in IP.
!   This is performed in medium precision.
!   Input: n, a, ip.
!   Output: ip.

use mpmodule
implicit none
integer i, iq, it, j, jq, jz, k, l, n
type (mp_realm) a(n), s0
integer ip(n), ik(50), jk(50)

do i = 1, n
  ip(i) = i
enddo

if (n == 1) return

k = 1
ik(1) = 1
jk(1) = n

130 continue

i = ik(k)
j = jk(k)
iq = i
jq = j
it = (i + j + 1) / 2
l = ip(j)
ip(j) = ip(it)
ip(it) = l
s0 = a(ip(j))
j = j - 1

140 continue

do l = i, j
  if (s0 < a(ip(l))) goto 160
enddo

i = j
goto 190

160 continue

i = l

do l = j, i, -1
  if (s0 > a(ip(l))) goto 180
enddo

j = i
goto 190

180 continue

j = l
if (i >= j) goto 190
l = ip(i)
ip(i) = ip(j)
ip(j) = l
goto 140

190 continue

if (s0 >= a(ip(i))) goto 200
l = ip(jq)
ip(jq) = ip(i)
ip(i) = l

200 continue

k = k - 1
jz = 0
if (j == iq) goto 210
k = k + 1
jk(k) = j
jz = 1

210 continue

i = i + 1
if (i == jq) goto 220
k = k + 1
ik(k) = i
jk(k) = jq
if (jz == 0) goto 220
if (j - iq >= jq - i) goto 220
ik(k-1) = i
jk(k-1) = jq
ik(k) = iq
jk(k) = j

220 continue

if (k > 0) goto 130

return
end
