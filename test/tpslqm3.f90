!*****************************************************************************

!  program tpslqm3

!  Revision date:  12 Jan 2023

!  AUTHOR:
!   David H. Bailey
!   Lawrence Berkeley National Lab (retired) and University of California, Davis
!   Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!   All software in this package (c) 2023 David H. Bailey.
!   By downloading or using this software you agree to the copyright, disclaimer
!   and license agreement in the accompanying file DISCLAIMER.txt.

!  DESCRIPTION:
!   This program demonstrates pslqm3, which performs the three-level multipair
!   PSLQ algorithm on an input vector.  A variety of sample input vectors can
!   be generated as inputs to pslqm3, as given in the parameters below.
!   For additional details, see:

!   David H. Bailey and David J. Broadhurst, "Parallel integer relation
!   detection: Techniques and applications," Mathematics of Computation,
!   vol. 70, no. 236 (Oct 2000), pg. 1719-1736, preprint available at
!   http://www.davidhbailey.com/dhbpapers/ppslq.pdf.

!   The pslqm2 routine is 100% THREAD SAFE -- all requisite parameters and
!   arrays are passed through subroutine arguments.

!   PSLQM3 parameters set below; all are default (4-byte) integer, except fname
!     is character(64):
!     idb   Debug level (0-4); default = 2.
!     n     Integer relation vector length; default = 73.
!           For minimal polynomial problems, n = 1 + polynomial degree.
!     ndp   Full precision level in digits; default = 1200.
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

!   TPSLQM3 parameters set below; all are default (4-byte) integer:
!     kq    0: for the algebraic case [1, alpha, alpha^2, ... alpha^(n-1)], where
!             alpha = 3^(1/kr) - 2^(1/ks); this is the default.
!           1: for testing algebraic relations of a number read from a file.
!           2: for testing additive relations of numbers read from a file.
!           3: for testing multiplicative relations of numbers read from a file.
!           4: for custom input.
!     kr    Degree of root of 3 when kq = 0; default = 8.
!     ks    Degree of root of 2 when kq = 0; default = 9.
!     lcx   Size of character*1 array chr1; default = 64.
!           ***Must be a multiple of 64.
!           ***Must be > size in digits of largest result coefficient.

program tpslqm3
use mpmodule
implicit none

!  TPSLQM3 parameters:

integer, parameter:: kq = 0, kr = 8, ks = 9, lcx = 64

!   PSLQM3 parameters:

integer, parameter:: idb = 2, n = kr*ks + 1, ndp = 1200, ndpm = max (ndp / 10, 100), &
  ndr = 30, nep = 30 - ndp, nepm = 20 - ndpm, nrb = 200, nrs = 0, &
  nwds = int (ndp / mpdpw + 2), nwdsm = int (ndpm / mpdpw + 2)
character(64), parameter:: fname= 'pslqm3.rst'

real (mprknd) second, tm0, tm1
integer i, i1, iq, j, j1, nn
character(1) chr1(lcx)
character(64) form4, nam(n), namx
integer lnm(n)
type (mp_real) alpha, r(n), x(n)
external second

!   Uncomment this line for MPFUN20-Fort when precision > 20,000 digits.

! call mpinit (nwds)

!   Check to see if default precision is high enough.

if (ndp > mpipl .or. ndpm > mpiplm) then
  write (6, 1) ndp, ndpm
1 format ('Increase the default standard precision level to at least', &
    i8,' digits'/'and the default medium precision level to at least', &
    i8,' digits. See module MPFUNF.' )
  stop
endif

nn = n
write (6, 2) nn, kq, kr, ks, ndr, nrb, nrs, ndp, ndpm, nep, nepm
2 format ('PSLQM3 Test Program'/ &
  'nn =',i4,3x,'kq =',i2,3x,'kr =',i2,3x,'ks =',i2/ &
  'ndr =',i5,3x,'nrb =',i5,3x,'nrs =',i5/ &
  'Full precision level ndp =',i6,' Intermediate precision ndpm =',i6/ &
  'Full prec. epsilon level nep = ',i6,' Intermediate epsilon nepm =', i6)
 
if (kq == 1 .or. kq == 2 .or. kq == 3) then
  open (11, file = 'pslq.inp')
  rewind 11
endif

if (kq == 0) then

!   This code generates alpha = 3^(1/kr) - 2^(1/ks). alpha is algebraic of degree
!   kr * ks.  Set n = kr * ks + 1 to recover the polynomial satisfied by alpha.

  alpha = mpnrt (mpreald (3.d0, nwds), kr) - mpnrt (mpreald (2.d0, nwds), ks)
elseif (kq == 1) then

!   Read an algebraic constant from a file.

  call mpread (11, alpha, nwds)
elseif (kq == 2) then

!   Read constants from a file for additive test.

  do i = 1, nn
    call mpread (11, alpha, nwds)
    x(i) = alpha
    write (namx, '(''con'',i3.3)') i
    nam(i) = namx(1:6)
    lnm(i) = 6
  enddo
elseif (kq == 3) then

!   Read constants from a file for multiplicative test.

  do i = 1, nn
    call mpread (11, alpha, nwds)
    x(i) = log (alpha)
    write (namx, '(''log(con'',i3.3,'')'')') i
    nam(i) = namx(1:11)
    lnm(i) = 11
  enddo
elseif (kq == 4) then

!   Produce X vector by a custom scheme.

endif

!   If kq is 0 or 1, generate x = [1, alpha, alpha^2, ..., alpha^(n-1)].

if (kq == 0 .or. kq == 1) then
  x(1) = mpreald (1.d0, nwds)
  nam(1) = '1'
  lnm(1) = 1

  do i = 2, nn
    x(i) = alpha * x(i-1)
    write (namx, '(''al^'',i3)') i - 1
    nam(i) = namx(1:6)
    lnm(i) = 6
  enddo
endif

!   Perform relation search.

tm0 = second ()
call pslqm3 (idb, nn, nwds, nwdsm, ndr, nrb, nrs, fname, nep, nepm, x, iq, r)
tm1 = second ()

!   Output relation, if one was found.

if (iq == 1) then
  write (6, 4)
4 format (/'Recovered relation: 0 =')
  form4 = '(64a1'
  i1 = 5

  do i = 65, lcx, 64
    form4(i1+1:i1+9) =  '/64a1'
    i1 = i1 + 5
  enddo

  form4(i1+1:i1+9) = '," * ",a)'
  i1 = i1 + 9
  
  do i = 1, nn
    if (r(i) /= 0.d0) then
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
end program tpslqm3

!------------------------------

!   The following code performs the three-level, multi-pair PSLQ algorithm.
!   David H. Bailey   12 Jan 2023

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
!     nsq   int      Size of tables used in iterdp and itermpm; default = 8.
!     dreps double   Tolerance for DP dynamic range check; default = 1d-10.

use mpmodule
implicit none
integer, intent(in):: idb, n, nwds, nwdsm, ndr, nrb, nrs, nep, nepm
character (64), intent(in):: fname
type (mp_real), intent(inout):: x(n)
integer, intent(out):: iq
type (mp_real), intent(out):: r(n)
real (mprknd), external:: dplog10, second
type (mp_realm), external:: bound, dynrange, dynrangem
real (mprknd), parameter:: dreps = 1.d-10
integer, parameter:: ipi = 500, ipm = 10, itm = 10000000, ndrm = 25, nsq = 8
integer i, i1, imq, it, its, izd, izm, izmm, j, j1, &
  n1, n2, n3, n4, n5, ndpm2, nwdsm2
real (mprknd) da(n,n), db(n,n), dh(n,n), dsyq(n,nsq), dy(n), d1, d2, d3, d4, d5, &
  tm0, tm1, times(6)
type (mp_real) b(n,n), h(n,n), y(n), eps, t1, t2, t3, t4
type (mp_realm) epsm, epsm2, wa(n,n), &
  wb(n,n), wh(n,n), wsyq(n,nsq), wy(n), wn, w1, w2

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
its = it

110 continue

!   Check if dynamic range of the wy vector is too great for DP iterations (which
!   is often the case at the start of the run), or if a very large value appeared
!   in iterdp immediately after an MPM update. If so, then do MPM iterations.

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
!   0: Iteration was uneventful; continue with DP iterations.
!   1: Relation found or DP precision exhausted; perform MPM update from DP arrays.
!   2: Very large value appeared in iterdp; arrays were restored; perform a normal
!      MPM update from DP arrays, as in the case izd = 1. But if this happens
!      immediately after an MPM update, then start MPM iterations.

if (izd == 2 .and. it > its + 1) izd = 1
if (izd == 0) then
  goto 120
else

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
!   too small for DP iterations, or DP iteration produced a very large value
!   immediately after a MPM update. Update the MP arrays from the MPM arrays.

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

    do j = 1, n
      do i = 1, n
        wh(i,j) = mprealm (h(i,j), nwdsm)
      enddo
    enddo
  
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
t3 = mpreald (0.d0, nwds)
t4 = mpreald (0.d0, nwds)

!   Select the relation corresponding to the smallest y entry.

do j = 1, n
  if (abs (y(j)) < t1) then
    j1 = j
    t1 = abs (y(j))
  endif
  t2 = max (t2, abs (y(j)))
enddo

t3 = t1 / t2

do i = 1, n
  r(i) = b(j1,i)
enddo

!   The norm and norm bound calculation here are performed in medium precision.

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

!   Output the final norm bound and other information.

if (idb >= 1) then
  call mpdecmd (t1, d1, n1)
  call mpdecmd (t2, d2, n2)
  call mpdecmd (t3, d3, n3)
  call mpdecmd (w1, d4, n4)
  call mpdecmd (wn, d5, n5)
  write (6, 20) it, d1, n1, d2, n2, d3, n3, d5, n5
20 format ('Iteration',i8,2x,'Relation detected'/ &
  'Min, max, ratio of y =',0p,f11.6,'e',i6,f11.6,'e',i6,f11.6,'e',i6/ &
  'Max. bound =',f11.6,'e',i5)
  write (6, 21) j1, d4, n4, d1, n1
21 format ('Index of relation =',i4,3x,'Norm =',0p,f11.6,'e',i5,3x, &
  'Residual =',0p,f11.6,'e',i6)
endif

!   If run was successful, set iq = 1; otherwise output message.

if (d3 == 0.d0) n3 = -nep
if (n4 <= nrb .and. -n3 >= ndr) then
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
end subroutine pslqm3


!------------------------------

!   First-level subroutines.

subroutine saverest (iu, n, it, izm, izmm, izd, times, x, y, b, h)

!   This saves restart data to the file.

use mpmodule
implicit none
integer, intent(in):: iu, n, it, izm, izmm, izd
real (mprknd), intent(in):: times(6)
type (mp_real), intent(in):: b(n,n), h(n,n), x(n), y(n)

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
end subroutine saverest

subroutine setprec (nwdsm, n, nsq, wa, wb, wh, wsyq, wy)

!   This sets the working precision of the arrays wa, wb, wh, wsyq and wy
!   to nwdsm words.

use mpmodule
implicit none
integer, intent(in):: nwdsm, n, nsq
type (mp_realm), intent(inout):: wa(n,n), wb(n,n), wh(n,n), wsyq(n,nsq), wy(n)
integer i, j

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
end subroutine setprec

type (mp_realm) function dynrange (n, nwdsm, y)

!   This returns the dynamic range of y, i.e., ratio of min|y_k| to max|y_k|,
!   normally done using MPM precision. Here y is of type mp_real.

use mpmodule
implicit none
integer, intent(in):: n, nwdsm
type (mp_real), intent(in):: y(n)
integer i
type (mp_realm) t1, t2, t3

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
end function dynrange

type (mp_realm) function dynrangem (n, nwdsm, wy)

!   This returns the dynamic range of wy, i.e., ratio of min|wy_k| to max|wy_k|,
!   normally done using MPM precision. Here wy is of type mp_realm.

use mpmodule
implicit none
integer, intent(in):: n, nwdsm
type (mp_realm), intent(in):: wy(n)
integer i
type (mp_realm) t1, t2, t3 

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
end function dynrangem

subroutine initdp (idb, n, nsq, nwdsm, da, db, dh, dy, dsyq, wh, wy)

!   This initializes the DP arrays from the MPM arrays.
!   This is performed in medium precision.
!   Input:  idb, n, nsq, nwdsm, wh, wy.
!   Output: da, db, dh, dy, dsyq.

use mpmodule
implicit none
integer, intent(in):: idb, n, nsq, nwdsm
real (mprknd), intent(out):: da(n,n), db(n,n), dh(n,n), dy(n), dsyq(n,nsq)
type (mp_realm), intent(in):: wh(n,n), wy(n)
integer i, j
type (mp_realm) t1, t2

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
end subroutine initdp

subroutine initmp (idb, n, nwds, b, h, x, y)

!   Initializes MP arrays at the beginning.
!   This is performed in full precision.
!   Input: idb, n, nwds.
!   Output: b, h, x, y.

use mpmodule
implicit none
integer, intent(in):: idb, n, nwds
type (mp_real), intent(out):: b(n,n), h(n,n), x(n), y(n)
integer i, j
type (mp_real) s(n), t1

if (idb >= 4) then
  write (6, 1)
1 format ('initmp: Input x vector:')
  call matoutmd (1, n, x)
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
  call matoutmd (1, n, y)
  write (6, 3)
3 format ('initmp: Initial h matrix:')
  call matoutmd (n, n - 1, h)
endif

return
end subroutine initmp

subroutine initmpm (idb, n, nsq, nwdsm, wa, wb, wh, wy, wsyq, h, y)

!   This initializes the MPM arrays from the MP arrays.
!   This is performed in medium precision.
!   Input: idb, n, nsq, nwdsm, h, y.
!   Output: wa, wb, wh, wy, wsyq.

use mpmodule
implicit none
integer, intent(in):: idb, n, nsq, nwdsm
type (mp_real), intent(in):: h(n,n), y(n)
type (mp_realm), intent(out):: wa(n,n), wb(n,n), wh(n,n), wy(n), wsyq(n,nsq)
integer i, j
type (mp_realm) t1, t2, t3

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
  call matoutmmd (1, n, wy)
  write (6, 4)
4 format ('initmpm: Factored wh matrix:')
  call matoutmmd (n, n - 1, wh)
endif

100 continue

return
end subroutine initmpm

subroutine iterdp (idb, it, n, nsq, da, db, dh, dsyq, dy, imq, izd)


!   This performs one iteration of the PSLQ algorithm using DP arithmetic.
!   Input: idb, it, n, nsq, da, db, dh, dsyq, dy, imq.
!   Output: da, db, dh, dsyq, dy, imq, izd.

!   NOTE: Parameter tmx2 = 2^52, not 2^53, so as to ensure that values > 2^53
!   never arise, even as intermediate values, in the update loop below.
!   Gam = sqrt (4/3) = 1.15470053837925153d0 by default.

use mpmodule
implicit none
integer, intent(in):: idb, it, n, nsq
real (mprknd), intent(inout):: da(n,n), db(n,n), dh(n,n), dsyq(n,nsq), dy(n)
integer, intent(inout):: imq
integer, intent(out):: izd
real (mprknd), parameter:: tmx1 = 1.d13, tmx2 = 2.d0**52, deps = 1.d-14, &
  gam = 1.15470053837925153d0
integer i, ii, ij, im, im1, j, j1, j2, k, mpr, mq
integer ip(n), ir(n), is(n)
real (mprknd) dq(n), dt(n,n), dsa(n,n), dsb(n,n), dsh(n,n), dsy(n), &
  t1, t2, t3, t4

!   Save the input DP arrays in case the iteration is aborted below.

call  savedp (n, da, db, dh, dy, dsa, dsb, dsh, dsy)
izd = 0
mpr = nint (0.4d0 * n)

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
  if (is(j1) /= 0 .or. is(j2) /= 0) goto 100
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

!  A small value appeared in dy; set izd = 1.

  if (idb >= 3) write (6, 1) it, t1
1 format ('Iteration',i8,2x,'iterdp: Small value in dy =',1pd15.6)
  izd = 1
endif

if (t2 > tmx1 .and. t2 <= tmx2) then

!  A large (but not very large) value appeared in da or db; set izd = 1.

  if (idb >= 3) write (6, 2) it, t2
2 format ('Iteration',i8,2x,'iterdp: Large value in da or db =',1pd15.6)
  izd = 1
elseif (t2 > tmx2) then

!   A very large value appeared in da or db; restore original arrays, set izd = 2
!   and return.

  if (idb >= 2) write (6, 3) it, t2
3 format ('Iteration',i8,2x,'iterdp: Very large value in da or db =',1pd15.6)
  call savedp (n, dsa, dsb, dsh, dsy, da, db, dh, dy)
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
end subroutine iterdp

subroutine itermpm (idb, it, n, nsq, nwdsm, epsm, wa, wb, wh, wsyq, wy, imq, izmm)

!   This performs one iteration of the PSLQM algorithm using MPM arithmetic.
!   This is performed in medium precision.
!   Input: idb, it, n, nsq, nwdsm, wa, wb, wh, wsyq, wy, epsm, imq.
!   Output: wa, wb, wh, wsyq, wy, imq, izmm.

use mpmodule
implicit none
integer, intent(in):: idb, it, n, nsq, nwdsm
type (mp_realm), intent(inout):: wa(n,n), wb(n,n), wh(n,n), wsyq(n,nsq), wy(n)
integer, intent(inout):: imq
integer, intent(out):: izmm
integer i, ii, ij, im, im1, j, j1, j2, k, mpr, mq, n1, ip(n), ir(n), is(n)
real (mprknd) d1
type (mp_realm) gam, t1, t2, t3, t4, epsm, tmx1, tmx2, wq(n), wt(n,n)

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
  if (is(j1) /= 0 .or. is(j2) /= 0) goto 100
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
  call matoutmmd (1, n, wy)
  write (6, 6)
6 format ('itermpm: Updated wa matrix:')
  call matoutmmd (n, n, wa)
  write (6, 7)
7 format ('itermpm: Updated wb matrix:')
  call matoutmmd (n, n, wb)
  write (6, 8)
8 format ('itermpm: Updated wh matrix:')
  call matoutmmd (n, n - 1, wh)
endif

200 continue

return
end subroutine itermpm

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

subroutine lqmpm (n, m, nwdsm, h)

!   This performs an LQ decomposition on the MPM matrix h.  It is a simplified
!   and transposed adaptation of the subroutine dqrdc from Linpack.
!   This is normally performed in variable medium precision.
!   Input: n, m, nwdsm, h.
!   Output: h.

use mpmodule
implicit none
integer, intent(in):: n, m, nwdsm
type (mp_realm), intent(inout):: h(n,m)
integer i, j, l, lup, ml
type (mp_realm) nrmxl, one, t, zero

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
  if (h(l,l) /= zero) nrmxl = sign (nrmxl, h(l,l))
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
end subroutine lqmpm

subroutine savedp (n, da, db, dh, dy, dsa, dsb, dsh, dsy)

!   This saves the arrays dy, da, db, dh in case DP iterations must be aborted.
!   A call to the same routine, with (da,db,dh,dy) and (dsa,dsb,dsh,dsy)
!   exchanged, serves to restore these arrays.
!   Input: n, da, db, dh, dy.
!   Output: dsa, dsb, dsh, dsy.

use mpmodule
implicit none
integer, intent(in):: n
real (mprknd), intent(in):: da(n,n), db(n,n), dh(n,n), dy(n)
real (mprknd), intent(out):: dsa(n,n), dsb(n,n), dsh(n,n), dsy(n)
integer i, j

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
end subroutine savedp

subroutine updtmp (idb, it, n, nwds, wa, wb, eps, b, h, y, izm)

!   This update the MP arrays from the MPM arrays.
!   This is performed in full precision.
!   Input: idb, it, n, nwds, wa, wb, eps, b, h, y.
!   Output: b, h, y, izm.

use mpmodule
implicit none
integer, intent(in):: idb, it, n, nwds
type (mp_realm), intent(in):: wa(n,n), wb(n,n)
type (mp_real), intent(in):: eps
type (mp_real), intent(inout):: b(n,n), h(n,n), y(n)
integer, intent(out):: izm
integer i, i1, n1, n2
real (mprknd) d1, d2
type (mp_real) t1, t2

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
  call matoutmd (1, n, y)
  write (6, 5)
5 format ('updtmp: Updated b matrix:')
  call matoutmd (n, n, b)
  write (6, 6)
6 format ('updtmp: Updated h matrix:')
  call matoutmd (n, n - 1, h)
endif

return
end subroutine updtmp

subroutine updtmpm (idb, it, n, nwdsm, epsm, da, db, dreps, wa, wb, wh, wy, izmm)

!   This updates the MPM arrays from the DP arrays.
!   This is performed in medium precision.
!   Input: idb, it, n, nwdsm, epsm, da, db, dreps, wa, wb, wh, wy.
!   Output: wa, wb, wh, wy, izmm.

use mpmodule
implicit none
integer, intent(in):: idb, it, n, nwdsm
real (mprknd) da(n,n), db(n,n), dreps
type (mp_realm), intent(inout):: wa(n,n), wb(n,n), wh(n,n), wy(n)
integer, intent(out):: izmm
integer i, j, n1, n2
real (mprknd) d1, d2
type (mp_realm) t1, t2, epsm, tmx1, tmx2, w1

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
  call matoutmmd (1, n, wy)
  write (6, 6)
6 format ('updtmpm: Updated wa matrix:')
  call matoutmmd (n, n, wa)
  write (6, 7)
7 format ('updtmpm: Updated wb matrix:')
  call matoutmmd (n, n, wb)
  write (6, 8)
8 format ('updtmpm: Updated wh matrix:')
  call matoutmmd (n, n - 1, wh)
endif

100 continue

return
end subroutine updtmpm

!------------------------------

!   Second- and third-level subroutines.

type (mp_realm) function bound (n, nwdsm, wh)

!   This computes the norm bound, normally done using MPM arithmetic.

use mpmodule
implicit none
integer, intent(in):: n, nwdsm
type (mp_realm), intent(inout):: wh(n,n)
integer i
type (mp_realm) t1

call lqmpm (n, n - 1, nwdsm, wh)
t1 = mprealdm (0.d0, nwdsm)

do i = 1, n - 1
  t1 = max (t1, abs (wh(i,i)))
enddo

bound = 1.d0 / t1
return
end function bound

real (mprknd) function dplog10 (a)

!   For input MPM value a, this routine returns a DP approximation to log10 (a).

use mpmodule
implicit none
type (mp_realm), intent(in):: a
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

subroutine matoutdp (n1, n2, a)

!   This outputs the DP matrix a.

use mpmodule
implicit none
integer, intent(in):: n1, n2
real (mprknd), intent(in):: a(n1,n2)
integer i, j

do i = 1, n1
  write (6, 1) i
1 format ('Row',i3)
  write (6, 2) (a(i,j), j = 1, n2)
2 format (1p,5d15.5)
enddo

return
end subroutine matoutdp

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
2 format (4(f13.8,'e',i6))
enddo

return
end subroutine matoutmd

subroutine matoutmmd (n1, n2, a)

!   This outputs the MPM matrix a as a DP matrix.

use mpmodule
implicit none
integer, intent(in):: n1, n2
type (mp_realm), intent(in):: a(n1,n2)
integer i, j, ix(n2)
real (mprknd) dx(n2)

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
end subroutine matoutmmd

subroutine matoutmp (n1, n2, a)

!   This outputs the MP or MPM matrix a.  It may be used in place of calls to
!   matoutmd in the code above if greater accuracy is desired in debug output.

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

subroutine mxmdm (n1, n2, nwdsm, a, b)

!   This multiplies the DP square matrix a by the MPM matrix b, and the result
!   is placed in b.  n1, n2 are the matrix dimensions as indicated below.
!   This is performed in medium precision.
!   Input: n1, n2, nwdsm, a, b.
!   Output: b.

use mpmodule
implicit none
integer, intent(in):: n1, n2, nwdsm
real (mprknd), intent(in):: a(n1,n1)
type (mp_realm), intent(inout):: b(n1,n2)
integer i, j, k
type (mp_realm) c(n1,n2)

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
end subroutine mxmdm

subroutine mxm (n1, n2, nwds, a, b)

!   This multiplies the MPM square matrix a by the MP matrix b, and the
!   result is returned in b.  n1, n2 are the matrix dimensions as indicated below.
!   This is performed in full precision.
!   Input: n1, n2, nwds, a, b.
!   Output: b.

use mpmodule
implicit none
integer, intent(in):: n1, n2, nwds
type (mp_realm), intent(in):: a(n1,n1)
type (mp_real), intent(inout):: b(n1,n2)
integer i, j, k
type (mp_real) c(n1,n2)

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
end subroutine mxm

subroutine qsortdp (n, a, ip)

!   This routine sorts the entries of the N-long DP vector A into ascending
!   order using the quicksort algorithm.  The permutation vector that would
!   sort the vector is returned in IP.
!   Input: n, a.
!   Output: ip.

use mpmodule
implicit none
integer, intent(in):: n
real (mprknd), intent(in):: a(n)
integer, intent(out):: ip(n)
integer i, iq, it, j, jq, jz, k, l
real (mprknd) s0
integer ik(100), jk(100)

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
end subroutine qsortdp

subroutine qsortmpm (n, a, ip)

!   This routine sorts the entries of the N-long MPM vector A into ascending
!   order using the quicksort algorithm.  The permutation vector that would
!   sort the vector is returned in IP.
!   This is performed in medium precision.
!   Input: n, a, ip.
!   Output: ip.

use mpmodule
implicit none
integer, intent(in):: n
type (mp_realm), intent(in):: a(n)
integer, intent(out):: ip(n)
integer i, iq, it, j, jq, jz, k, l
type (mp_realm) s0
integer ik(100), jk(100)

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
end subroutine qsortmpm
