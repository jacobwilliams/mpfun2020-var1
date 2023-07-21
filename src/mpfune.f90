
!  module mpfune

!  This is the special special function module (MPFUNE) for the MPFUN20-Fort arbitrary
!  precision library, after conversion using the author's convmp.f90 program, and for the
!  MPFUN20-MPFR library, after conversion using the author's convmpfr.f90 program.

!  Revision date:  10 Jun 2023

!  AUTHOR:
!     David H. Bailey
!     Lawrence Berkeley National Lab (retired)
!     Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2023 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!  DOCUMENTATION:
!    For MPFUN20-Fort, see README-mpfun20-fort.txt file in the mpfun20-fort directory.
!    For MPFUN20-Fort, see README-mpfun20-mpfr.txt file in the mpfun20-mpfr directory.

!  DESCRIPTION OF THIS MODULE (MPFUNE):
!    This module contains all special functions.

!  NOTE:
!    The special comments !mp placed in the file below are annotations to aid in this
!    conversion. See the programs convmp.f90 and convmpfr.f90 for details.

module mpfune
use mpfund
integer, parameter:: mpdk = mpdknd

contains

!   Special functions start here.

subroutine mpberner (nb1, nb2, berne, mpnw)

!   This returns the array berne, containing Bernoulli numbers indexed 2*k for
!   k = 1 to n, using an advanced polynomial Newton iteration scheme as described in
!   the main documentation paper.

implicit none
integer, intent(in):: nb1, mpnw
integer, intent(in):: nb2
!mp dim1="0:nb1+5"
! type (mp_real), intent(out):: berne(1:nb2)
integer (mpiknd), intent(out):: berne(0:nb1+5,1:nb2)
integer, parameter:: itrmax = 10000
integer, parameter:: ibz = 66
integer i, i1, i2, ic1, ic2, j, kn, n, n1, nn1
real (mpdknd) d1, dd1, dd2, dd3
!mp dim1="0:mpnwm+5"
! type (mp_real) eps, tc1, tc2
integer (mpiknd) eps(0:mpnwm+5), tc1(0:mpnwm+5), tc2(0:mpnwm+5)
!mp dim1="0:mpnw+6"
! type (mp_real) c1(0:nb2), p1(0:nb2), p2(0:nb2), q(0:nb2)
integer (mpiknd) c1(0:mpnw+6,0:nb2), p1(0:mpnw+6,0:nb2), p2(0:mpnw+6,0:nb2), q(0:mpnw+6,0:nb2)
! type (mp_real) cp2, q1, r(0:nb2), s(0:nb2)
integer (mpiknd) cp2(0:mpnw+6), q1(0:mpnw+6), r(0:mpnw+6,0:nb2), s(0:mpnw+6,0:nb2)
! type (mp_real) t1, t2, t3, t4, t5
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), t5(0:mpnw+6)
integer mpnw1, mpnw2

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpt3(0:mpnw+6)
integer (mpiknd) mpt4(0:mpnw+6), mpt5(0:mpnw+6), mpt6(0:mpnw+6)

call mpinitwds (mpt1, mpnw+6-5)
call mpinitwds (mpt2, mpnw+6-5)
call mpinitwds (mpt3, mpnw+6-5)
call mpinitwds (mpt4, mpnw+6-5)
call mpinitwds (mpt5, mpnw+6-5)
call mpinitwds (mpt6, mpnw+6-5)
call mpinitwds (eps(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc1(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc2(0:mpnwm+5), mpnwm+5-5)
do mpi1 = 0,nb2
  call mpinitwds (c1(0:mpnw+6,mpi1), mpnw+6-5)
enddo
do mpi1 = 0,nb2
  call mpinitwds (p1(0:mpnw+6,mpi1), mpnw+6-5)
enddo
do mpi1 = 0,nb2
  call mpinitwds (p2(0:mpnw+6,mpi1), mpnw+6-5)
enddo
do mpi1 = 0,nb2
  call mpinitwds (q(0:mpnw+6,mpi1), mpnw+6-5)
enddo
call mpinitwds (cp2(0:mpnw+6), mpnw+6-5)
call mpinitwds (q1(0:mpnw+6), mpnw+6-5)
do mpi1 = 0,nb2
  call mpinitwds (r(0:mpnw+6,mpi1), mpnw+6-5)
enddo
do mpi1 = 0,nb2
  call mpinitwds (s(0:mpnw+6,mpi1), mpnw+6-5)
enddo
call mpinitwds (t1(0:mpnw+6), mpnw+6-5)
call mpinitwds (t2(0:mpnw+6), mpnw+6-5)
call mpinitwds (t3(0:mpnw+6), mpnw+6-5)
call mpinitwds (t4(0:mpnw+6), mpnw+6-5)
call mpinitwds (t5(0:mpnw+6), mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)
mpnw2 = mpnw1

do ic1 = 1, nb2
!   ic2 = mpspacer (berne(ic1))
  ic2 = mpspacer (berne(0:nb1+5,ic1))
  if (mpnw < mpnwm .or. ic2 < mpnw + 6) then
    write (6, 1)
1   format ('*** BERNER: Uninitialized or inadequately sized array')
!     call mpabrt (501)
    call mpabrt (501)
  endif
enddo

! ic1 = mpwprecr (mppicon)
ic1 = mpwprecr (mppicon(0:))
if (ic1 < mpnw2) then
  write (6, 2) mpnw2
2 format ('*** BERNER: Pi must be initialized to',i6,' words by calling mpinit.')
!   call mpabrt (502)
  call mpabrt (502)
endif

!mp prec="mpnw2"

! cp2 = mppicon ** 2
call mpnpwr (mppicon(0:), 2, mpt1(0:mpnw+6), mpnw2)
call mpeq (mpt1(0:mpnw+6), cp2(0:mpnw+6), mpnw2)
! c1(0) = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw2)
call mpeq (mpt1(0:mpnw+6), c1(0:mpnw+6,0), mpnw2)
! p1(0) = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw2)
call mpeq (mpt1(0:mpnw+6), p1(0:mpnw+6,0), mpnw2)
! p2(0) = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw2)
call mpeq (mpt1(0:mpnw+6), p2(0:mpnw+6,0), mpnw2)
! q(0) = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw2)
call mpeq (mpt1(0:mpnw+6), q(0:mpnw+6,0), mpnw2)
n = nb2

!   Construct numerator and denominator polynomials.

do i = 1, n
!   c1(i) = mprealdp (0.e0_mpdk)
  call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw2)
  call mpeq (mpt1(0:mpnw+6), c1(0:mpnw+6,i), mpnw2)
  dd1 = 2.e0_mpdk * real (i+1, mpdk) - 3.e0_mpdk
  dd2 = dd1 + 1.e0_mpdk
  dd3 = dd2 + 1.e0_mpdk
!   t1 = cp2 * p1(i-1)
  call mpmul (cp2(0:mpnw+6), p1(0:mpnw+6,i-1), mpt1(0:mpnw+6), mpnw2)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw2)
!   p1(i) = t1 / dd1 / dd2
  call mpdivd (t1(0:mpnw+6), dd1, mpt1(0:mpnw+6), mpnw2)
  call mpdivd (mpt1(0:mpnw+6), dd2, mpt2(0:mpnw+6), mpnw2)
  call mpeq (mpt2(0:mpnw+6), p1(0:mpnw+6,i), mpnw2)
!   t1 = cp2 * p2(i-1)
  call mpmul (cp2(0:mpnw+6), p2(0:mpnw+6,i-1), mpt1(0:mpnw+6), mpnw2)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw2)
!   p2(i) = t1 / dd2 / dd3
  call mpdivd (t1(0:mpnw+6), dd2, mpt1(0:mpnw+6), mpnw2)
  call mpdivd (mpt1(0:mpnw+6), dd3, mpt2(0:mpnw+6), mpnw2)
  call mpeq (mpt2(0:mpnw+6), p2(0:mpnw+6,i), mpnw2)
!   q(i) = mprealdp (0.e0_mpdk)
  call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw2)
  call mpeq (mpt1(0:mpnw+6), q(0:mpnw+6,i), mpnw2)
enddo

kn = min (mpnwm, nb2)
mpnw2 = mpnwm
! tc2 = mprealdp (2.e0_mpdk)
call mprealdp (2.e0_mpdk, mpt1(0:mpnw+6), mpnw2)
call mpeq (mpt1(0:mpnw+6), tc2(0:mpnwm+5), mpnw2)
ic2 = ibz - mpnw2*mpnbt
! eps = tc2 ** ic2
call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:mpnw+6), mpnw2)
call mpeq (mpt1(0:mpnw+6), eps(0:mpnwm+5), mpnw2)

!mp prec="mpnw2"
! q1 = mprealdp (0.e0_mpdk)
call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw2)
call mpeq (mpt1(0:mpnw+6), q1(0:mpnw+6), mpnw2)

!   Perform Newton iterations with dynamic precision levels, using an
!   iteration formula similar to that used to evaluate reciprocals.

do j = 1, itrmax
!mp prec="mpnw2"
!   call mppolymul (mpnw1, kn, p2, q, r, mpnw2)
  call mppolymul (mpnw1, kn, p2, q, r, mpnw2)
!   call mppolysub (mpnw1, kn, c1, r, s, mpnw2)
  call mppolysub (mpnw1, kn, c1, r, s, mpnw2)
!   call mppolymul (mpnw1, kn, s, q, r, mpnw2)
  call mppolymul (mpnw1, kn, s, q, r, mpnw2)
!   call mppolyadd (mpnw1, kn, q, r, q, mpnw2)
  call mppolyadd (mpnw1, kn, q, r, q, mpnw2)
!   t1 = q(kn) - q1
  call mpsub (q(0:mpnw+6,kn), q1(0:mpnw+6), mpt1(0:mpnw+6), mpnw2)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw2)
!mp prec="mpnwm"
!   tc1 = abs (t1) - eps
  call mpabs (t1(0:mpnw+6), mpt1(0:mpnw+6), mpnwm)
  call mpsub (mpt1(0:mpnw+6), eps(0:mpnwm+5), mpt2(0:mpnw+6), mpnwm)
  call mpeq (mpt2(0:mpnw+6), tc1(0:mpnwm+5), mpnwm)
!   ic1 = mpsgn (tc1)
  mpi1 = mpsgn (tc1(0:mpnwm+5))
  ic1 = mpi1
  if (ic1 < 0) then
    if (kn == n .and. mpnw2 == mpnw1) then
      goto 100
    elseif (kn < n) then
      kn = min (2 * kn, n)
!       q1 = mprealdp (0.e0_mpdk)
      call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnwm)
      call mpeq (mpt1(0:mpnw+6), q1(0:mpnw+6), mpnwm)
    elseif (mpnw2 < mpnw1) then
      mpnw2 = min (2 * mpnw2, mpnw1)
!mp prec="mpnwm"
!       tc2 = mprealdp (2.e0_mpdk)
      call mprealdp (2.e0_mpdk, mpt1(0:mpnw+6), mpnwm)
      call mpeq (mpt1(0:mpnw+6), tc2(0:mpnwm+5), mpnwm)
      ic2 = ibz - mpnw2*mpnbt
!       eps = tc2 ** ic2
      call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:mpnw+6), mpnwm)
      call mpeq (mpt1(0:mpnw+6), eps(0:mpnwm+5), mpnwm)
!mp prec="mpnw2"
!       q1 = mprealdp (0.e0_mpdk)
      call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw2)
      call mpeq (mpt1(0:mpnw+6), q1(0:mpnw+6), mpnw2)
    endif
  else
!     q1 = q(kn)
    call mpeq (q(0:mpnw+6,kn), q1(0:mpnw+6), mpnw2)
  endif
enddo

write (mpldb, 3)
3 format ('*** BERNER: Loop end error')
! call mpabrt (503)
call mpabrt (503)

100 continue

!   Multiply numerator polynomial by reciprocal of denominator polynomial.

! call mppolymul (mpnw1, n, p1, q, r, mpnw2)
call mppolymul (mpnw1, n, p1, q, r, mpnw2)

!   Apply formula to produce Bernoulli numbers.

! t1 = mprealdp (-2.e0_mpdk)
call mprealdp (-2.e0_mpdk, mpt1(0:mpnw+6), mpnw2)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw2)
! t2 = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw2)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw2)

do i = 1, n
!mp prec="mpnw2"
  d1 = - real (2*i-1, mpdk) * real (2*i, mpdk)
!   t1 = t1 * d1
  call mpmuld (t1(0:mpnw+6), d1, mpt1(0:mpnw+6), mpnw2)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw2)
!   t3 = 4.e0_mpdk * cp2
  call mpmuld (cp2(0:mpnw+6), 4.e0_mpdk, mpt1(0:mpnw+6), mpnw2)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw2)
!   t2 = t3 * t2
  call mpmul (t3(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw2)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw2)
!   t4 = 0.5e0_mpdk * t1 / t2
  call mpmuld (t1(0:mpnw+6), 0.5e0_mpdk, mpt1(0:mpnw+6), mpnw2)
  call mpdiv (mpt1(0:mpnw+6), t2(0:mpnw+6), mpt2(0:mpnw+6), mpnw2)
  call mpeq (mpt2(0:mpnw+6), t4(0:mpnw+6), mpnw2)
!   t5 = t4 * abs (r(i))
  call mpabs (r(0:mpnw+6,i), mpt1(0:mpnw+6), mpnw2)
  call mpmul (t4(0:mpnw+6), mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpnw2)
  call mpeq (mpt2(0:mpnw+6), t5(0:mpnw+6), mpnw2)
!mp prec="mpnw"
!   berne(i) = t5
  call mpeq (t5(0:mpnw+6), berne(0:nb1+5,i), mpnw)
enddo

return
end subroutine mpberner

subroutine mppolyadd (mpnw1, n, a, b, c, mpnw2)

!   This adds two polynomials, as is required by mpberne.
!   The output array C may be the same as A or B.

implicit none
integer, intent(in):: mpnw1, mpnw2
integer, intent(in):: n
!mp dim1="0:mpnw1+5"
! type (mp_real), intent(in):: a(0:n), b(0:n)
integer (mpiknd), intent(in):: a(0:mpnw1+5,0:n), b(0:mpnw1+5,0:n)
! type (mp_real), intent(out):: c(0:n)
integer (mpiknd), intent(out):: c(0:mpnw1+5,0:n)
integer k
!mp dim1="0:mpnw2+5"
! type (mp_real) t1, t2
integer (mpiknd) t1(0:mpnw2+5), t2(0:mpnw2+5)

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:mpnw2+5), mpt2(0:mpnw2+5), mpt3(0:mpnw2+5)
integer (mpiknd) mpt4(0:mpnw2+5), mpt5(0:mpnw2+5), mpt6(0:mpnw2+5)

call mpinitwds (mpt1, mpnw2+5-5)
call mpinitwds (mpt2, mpnw2+5-5)
call mpinitwds (mpt3, mpnw2+5-5)
call mpinitwds (mpt4, mpnw2+5-5)
call mpinitwds (mpt5, mpnw2+5-5)
call mpinitwds (mpt6, mpnw2+5-5)
call mpinitwds (t1(0:mpnw2+5), mpnw2+5-5)
call mpinitwds (t2(0:mpnw2+5), mpnw2+5-5)
!mp prec="mpnw2"

do k = 0, n
!   t1 = a(k)
  call mpeq (a(0:mpnw1+5,k), t1(0:mpnw2+5), mpnw2)
!   t2 = b(k)
  call mpeq (b(0:mpnw1+5,k), t2(0:mpnw2+5), mpnw2)
!   c(k) = t1 + t2
  call mpadd (t1(0:mpnw2+5), t2(0:mpnw2+5), mpt1(0:mpnw2+5), mpnw2)
  call mpeq (mpt1(0:mpnw2+5), c(0:mpnw1+5,k), mpnw2)
enddo

return
end subroutine mppolyadd

subroutine mppolysub (mpnw1, n, a, b, c, mpnw2)

!   This adds two polynomials, as is required by mpberne.
!   The output array C may be the same as A or B.

implicit none
integer, intent(in):: mpnw1, mpnw2
integer, intent(in):: n
!mp dim1="0:mpnw1+5"
! type (mp_real), intent(in):: a(0:n), b(0:n)
integer (mpiknd), intent(in):: a(0:mpnw1+5,0:n), b(0:mpnw1+5,0:n)
! type (mp_real), intent(out):: c(0:n)
integer (mpiknd), intent(out):: c(0:mpnw1+5,0:n)
integer k
!mp dim1="0:mpnw2+5"
! type (mp_real) t1, t2
integer (mpiknd) t1(0:mpnw2+5), t2(0:mpnw2+5)

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:mpnw2+5), mpt2(0:mpnw2+5), mpt3(0:mpnw2+5)
integer (mpiknd) mpt4(0:mpnw2+5), mpt5(0:mpnw2+5), mpt6(0:mpnw2+5)

call mpinitwds (mpt1, mpnw2+5-5)
call mpinitwds (mpt2, mpnw2+5-5)
call mpinitwds (mpt3, mpnw2+5-5)
call mpinitwds (mpt4, mpnw2+5-5)
call mpinitwds (mpt5, mpnw2+5-5)
call mpinitwds (mpt6, mpnw2+5-5)
call mpinitwds (t1(0:mpnw2+5), mpnw2+5-5)
call mpinitwds (t2(0:mpnw2+5), mpnw2+5-5)
!mp prec="mpnw2"

do k = 0, n
!   t1 = a(k)
  call mpeq (a(0:mpnw1+5,k), t1(0:mpnw2+5), mpnw2)
!   t2 = b(k)
  call mpeq (b(0:mpnw1+5,k), t2(0:mpnw2+5), mpnw2)
!   c(k) = t1 - t2
  call mpsub (t1(0:mpnw2+5), t2(0:mpnw2+5), mpt1(0:mpnw2+5), mpnw2)
  call mpeq (mpt1(0:mpnw2+5), c(0:mpnw1+5,k), mpnw2)
enddo

return
end subroutine mppolysub

subroutine mppolymul (mpnw1, n, a, b, c, mpnw2)

!   This adds two polynomials (ignoring high-order terms), as is required
!   by mpberne. The output array C may NOT be the same as A or B.

implicit none
integer, intent(in):: mpnw1, mpnw2
integer, intent(in):: n
!mp dim1="0:mpnw1+5"
! type (mp_real), intent(in):: a(0:n), b(0:n)
integer (mpiknd), intent(in):: a(0:mpnw1+5,0:n), b(0:mpnw1+5,0:n)
! type (mp_real), intent(out):: c(0:n)
integer (mpiknd), intent(out):: c(0:mpnw1+5,0:n)
integer j, k
!mp dim1="0:mpnw2+5"
! type (mp_real) t0, t1, t2, t3
integer (mpiknd) t0(0:mpnw2+5), t1(0:mpnw2+5), t2(0:mpnw2+5), t3(0:mpnw2+5)

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:mpnw2+5), mpt2(0:mpnw2+5), mpt3(0:mpnw2+5)
integer (mpiknd) mpt4(0:mpnw2+5), mpt5(0:mpnw2+5), mpt6(0:mpnw2+5)

call mpinitwds (mpt1, mpnw2+5-5)
call mpinitwds (mpt2, mpnw2+5-5)
call mpinitwds (mpt3, mpnw2+5-5)
call mpinitwds (mpt4, mpnw2+5-5)
call mpinitwds (mpt5, mpnw2+5-5)
call mpinitwds (mpt6, mpnw2+5-5)
call mpinitwds (t0(0:mpnw2+5), mpnw2+5-5)
call mpinitwds (t1(0:mpnw2+5), mpnw2+5-5)
call mpinitwds (t2(0:mpnw2+5), mpnw2+5-5)
call mpinitwds (t3(0:mpnw2+5), mpnw2+5-5)
!mp prec="mpnw2"

do k = 0, n
!   t0 = mprealdp (0.e0_mpdk)
  call mprealdp (0.e0_mpdk, mpt1(0:mpnw2+5), mpnw2)
  call mpeq (mpt1(0:mpnw2+5), t0(0:mpnw2+5), mpnw2)

  do j = 0, k
!     t1 = a(j)
    call mpeq (a(0:mpnw1+5,j), t1(0:mpnw2+5), mpnw2)
!     t2 = b(k-j)
    call mpeq (b(0:mpnw1+5,k-j), t2(0:mpnw2+5), mpnw2)
!     t3 = t1 * t2
    call mpmul (t1(0:mpnw2+5), t2(0:mpnw2+5), mpt1(0:mpnw2+5), mpnw2)
    call mpeq (mpt1(0:mpnw2+5), t3(0:mpnw2+5), mpnw2)
!     t0 = t0 + t3
    call mpadd (t0(0:mpnw2+5), t3(0:mpnw2+5), mpt1(0:mpnw2+5), mpnw2)
    call mpeq (mpt1(0:mpnw2+5), t0(0:mpnw2+5), mpnw2)
  enddo

!   c(k) = t0
  call mpeq (t0(0:mpnw2+5), c(0:mpnw1+5,k), mpnw2)
enddo

return
end subroutine mppolymul

subroutine mpbesselinr (nu, rr, ss, mpnw)

!   This evaluates the modified Bessel function BesselI (NU,RR).
!   NU is an integer. The algorithm is DLMF formula 10.25.2 for modest RR,
!   and DLMF 10.40.1 for large RR, relative to precision.

implicit none
integer, intent(in):: mpnw
integer, intent(in):: nu
!mp dim1="0:"
! type (mp_real), intent(in):: rr
integer (mpiknd), intent(in):: rr(0:)
! type (mp_real), intent(out):: ss
integer (mpiknd), intent(out):: ss(0:)
integer, parameter:: itrmax = 1000000
real (mpdknd), parameter:: dfrac = 0.4e0_mpdk
integer i1, i2, ic1, ic2, k, nua, n1
real (mpdknd) d1, d2
!mp dim1="0:mpnwm+5"
! type (mp_real) eps, tc1, tc2
integer (mpiknd) eps(0:mpnwm+5), tc1(0:mpnwm+5), tc2(0:mpnwm+5)
!mp dim1="0:mpnw+6"
! type (mp_real) f1, f2, sum, td, tn
integer (mpiknd) f1(0:mpnw+6), f2(0:mpnw+6), sum(0:mpnw+6), td(0:mpnw+6), tn(0:mpnw+6)
! type (mp_real) t1, t2, t3, t4, rra
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), rra(0:mpnw+6)
integer mpnw1

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpt3(0:mpnw+6)
integer (mpiknd) mpt4(0:mpnw+6), mpt5(0:mpnw+6), mpt6(0:mpnw+6)

call mpinitwds (mpt1, mpnw+6-5)
call mpinitwds (mpt2, mpnw+6-5)
call mpinitwds (mpt3, mpnw+6-5)
call mpinitwds (mpt4, mpnw+6-5)
call mpinitwds (mpt5, mpnw+6-5)
call mpinitwds (mpt6, mpnw+6-5)
call mpinitwds (eps(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc1(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc2(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (f1(0:mpnw+6), mpnw+6-5)
call mpinitwds (f2(0:mpnw+6), mpnw+6-5)
call mpinitwds (sum(0:mpnw+6), mpnw+6-5)
call mpinitwds (td(0:mpnw+6), mpnw+6-5)
call mpinitwds (tn(0:mpnw+6), mpnw+6-5)
call mpinitwds (t1(0:mpnw+6), mpnw+6-5)
call mpinitwds (t2(0:mpnw+6), mpnw+6-5)
call mpinitwds (t3(0:mpnw+6), mpnw+6-5)
call mpinitwds (t4(0:mpnw+6), mpnw+6-5)
call mpinitwds (rra(0:mpnw+6), mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)

! ic1 = mpspacer (ss)
ic1 = mpspacer (ss(0:))
if (mpnw < mpnwm .or. ic1 < mpnw + 6) then
  write (6, 1)
1 format ('*** BESSELINR: Uninitialized or inadequately sized array')
!   call mpabrt (504)
  call mpabrt (504)
endif

! ic1 = mpwprecr (mppicon)
ic1 = mpwprecr (mppicon(0:))
if (ic1 < mpnw1) then
  write (6, 2) mpnw1
2 format ('*** BESSELINR: Pi must be initialized to',i6,' words by calling mpinit.')
!   call mpabrt (505)
  call mpabrt (505)
endif

!mp prec="mpnwm"

! tc2 = mprealdp (2.e0_mpdk)
call mprealdp (2.e0_mpdk, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), tc2(0:mpnwm+5), mpnwm)
ic2 = -mpnw1*mpnbt
! eps = tc2 ** ic2
call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), eps(0:mpnwm+5), mpnwm)

!mp prec="mpnw1"

!   Check for RR = 0.

! ic1 = mpsgn (rr)
mpi1 = mpsgn (rr(0:))
ic1 = mpi1
if (ic1 == 0) then
  write (mpldb, 3)
3 format ('*** BESSELINR: Second argument is zero')
!   call mpabrt (506)
  call mpabrt (506)
endif

nua = abs (nu)
! rra = abs (rr)
call mpabs (rr(0:), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), rra(0:mpnw+6), mpnw1)
! d1 = dpreal (rra)
call mpmdc (rra(0:mpnw+6), mpd1, mpi1, mpnw1)
mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
d1 = mpd1

if (d1 < dfrac * mpnw1 * mpnbt) then
!   tn = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), tn(0:mpnw+6), mpnw1)
!   f1 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), f1(0:mpnw+6), mpnw1)
!   f2 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), f2(0:mpnw+6), mpnw1)
!   t1 = 0.25e0_mpdk * rra ** 2
  call mpnpwr (rra(0:mpnw+6), 2, mpt1(0:mpnw+6), mpnw1)
  call mpmuld (mpt1(0:mpnw+6), 0.25e0_mpdk, mpt2(0:mpnw+6), mpnw1)
  call mpeq (mpt2(0:mpnw+6), t1(0:mpnw+6), mpnw1)

  do k = 1, nua
!     f2 = real (k, mpdk) * f2
    call mpmuld (f2(0:mpnw+6), real (k, mpdk), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), f2(0:mpnw+6), mpnw1)
  enddo

!   td = f1 * f2
  call mpmul (f1(0:mpnw+6), f2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), td(0:mpnw+6), mpnw1)
!   t2 = tn / td
  call mpdiv (tn(0:mpnw+6), td(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   sum = t2
  call mpeq (t2(0:mpnw+6), sum(0:mpnw+6), mpnw1)

  do k = 1, itrmax
!     f1 = f1 * real (k, mpdk)
    call mpmuld (f1(0:mpnw+6), real (k, mpdk), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), f1(0:mpnw+6), mpnw1)
    i1 = k + nua
!     f2 = f2 * real (i1, mpdk)
    call mpmuld (f2(0:mpnw+6), real (i1, mpdk), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), f2(0:mpnw+6), mpnw1)
!     tn = t1 * tn
    call mpmul (t1(0:mpnw+6), tn(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), tn(0:mpnw+6), mpnw1)
!     td = f1 * f2
    call mpmul (f1(0:mpnw+6), f2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), td(0:mpnw+6), mpnw1)
!     t2 = tn / td
    call mpdiv (tn(0:mpnw+6), td(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!     sum = sum + t2
    call mpadd (sum(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), sum(0:mpnw+6), mpnw1)
!mp prec="mpnwm"
!     tc1 = abs (t2) - eps * abs (sum)
    call mpabs (t2(0:mpnw+6), mpt1(0:mpnw+6), mpnwm)
    call mpabs (sum(0:mpnw+6), mpt2(0:mpnw+6), mpnwm)
    call mpmul (eps(0:mpnwm+5), mpt2(0:mpnw+6), mpt3(0:mpnw+6), mpnwm)
    call mpsub (mpt1(0:mpnw+6), mpt3(0:mpnw+6), mpt4(0:mpnw+6), mpnwm)
    call mpeq (mpt4(0:mpnw+6), tc1(0:mpnwm+5), mpnwm)
!     ic1 = mpsgn (tc1)
    mpi1 = mpsgn (tc1(0:mpnwm+5))
    ic1 = mpi1
    if (ic1 < 0) goto 100
  enddo

  write (mpldb, 4)
4 format ('*** BESSELINR: Loop end error 1')
!   call mpabrt (507)
  call mpabrt (507)

100 continue

!mp prec="mpnw1"
!   t1 = 0.5e0_mpdk * rra
  call mpmuld (rra(0:mpnw+6), 0.5e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   t3 = sum * t1 ** nua
  call mpnpwr (t1(0:mpnw+6), nua, mpt1(0:mpnw+6), mpnw1)
  call mpmul (sum(0:mpnw+6), mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
  call mpeq (mpt2(0:mpnw+6), t3(0:mpnw+6), mpnw1)
else
!   sum = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), sum(0:mpnw+6), mpnw1)
  i1 = 2 * nua
  d2 = real (i1, mpdk) ** 2
!   t1 = mprealdp (d2)
  call mprealdp (d2, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   tn = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), tn(0:mpnw+6), mpnw1)
!   td = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), td(0:mpnw+6), mpnw1)

  do k = 1, itrmax
    i1 = 2 * k - 1
!     t2 = t1 - real (i1, mpdk) ** 2
    mpd1 = real (i1, mpdk) ** 2
    call mpdmc (mpd1, 0, mpt2(0:mpnw+6), mpnw1)
    call mpsub (t1(0:mpnw+6), mpt2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!     tn = - tn * t2
    call mpmul (tn(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpneg (mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
    call mpeq (mpt2(0:mpnw+6), tn(0:mpnw+6), mpnw1)
    i2 = 8 * k
!     td = td * real (i2, mpdk) * rra
    call mpmuld (td(0:mpnw+6), real (i2, mpdk), mpt1(0:mpnw+6), mpnw1)
    call mpmul (mpt1(0:mpnw+6), rra(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
    call mpeq (mpt2(0:mpnw+6), td(0:mpnw+6), mpnw1)
!     t4 = tn / td
    call mpdiv (tn(0:mpnw+6), td(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!     sum = sum + t4
    call mpadd (sum(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), sum(0:mpnw+6), mpnw1)
!mp prec="mpnwm"
!     tc2 = abs (t4) - eps * abs (sum)
    call mpabs (t4(0:mpnw+6), mpt1(0:mpnw+6), mpnwm)
    call mpabs (sum(0:mpnw+6), mpt2(0:mpnw+6), mpnwm)
    call mpmul (eps(0:mpnwm+5), mpt2(0:mpnw+6), mpt3(0:mpnw+6), mpnwm)
    call mpsub (mpt1(0:mpnw+6), mpt3(0:mpnw+6), mpt4(0:mpnw+6), mpnwm)
    call mpeq (mpt4(0:mpnw+6), tc2(0:mpnwm+5), mpnwm)
!     ic2 = mpsgn (tc2)
    mpi1 = mpsgn (tc2(0:mpnwm+5))
    ic2 = mpi1
    if (ic2 < 0) goto 110
  enddo

write (mpldb, 5)
5 format ('*** BESSELINR: Loop end error 2')
! call mpabrt (508)
call mpabrt (508)

110 continue

!mp prec="mpnw1"
!   t1 = exp (rra)
  call mpexp (rra(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   t2 = 2.e0_mpdk * mppicon * rra
  call mpmuld (mppicon(0:), 2.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpmul (mpt1(0:mpnw+6), rra(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
  call mpeq (mpt2(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   t4 = sqrt (t2)
  call mpsqrt (t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!   t2 = t1 / t4
  call mpdiv (t1(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   t3 = t2 * sum
  call mpmul (t2(0:mpnw+6), sum(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
endif

!mp prec="mpnw"
! ic1 = mpsgn (rr)
mpi1 = mpsgn (rr(0:))
ic1 = mpi1
if (ic1 < 0 .and. mod (nu, 2) /= 0) then
!   t3 = - t3
  call mpneg (t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw)
endif
! ss = t3
call mpeq (t3(0:mpnw+6), ss(0:), mpnw)

return
end subroutine mpbesselinr

subroutine mpbesselir (qq, rr, ss, mpnw)

!   This evaluates the modified Bessel function BesselI (QQ,RR) for QQ and RR
!   both real. The algorithm is DLMF formula 10.25.2 for modest RR, and
!   DLMF 10.40.1 for large RR, relative to precision.

implicit none
integer, intent(in):: mpnw
!mp dim1="0:"
! type (mp_real), intent(in):: qq, rr
integer (mpiknd), intent(in):: qq(0:), rr(0:)
! type (mp_real), intent(out):: ss
integer (mpiknd), intent(out):: ss(0:)
integer ic1, ic2, i0, i1, i2, k, n1
integer, parameter:: itrmax = 1000000
real (mpdknd), parameter:: dfrac = 0.4e0_mpdk
real (mpdknd) d1, d2
!mp dim1="0:mpnwm+5"
! type (mp_real) eps, tc1, tc2
integer (mpiknd) eps(0:mpnwm+5), tc1(0:mpnwm+5), tc2(0:mpnwm+5)
!mp dim1="0:mpnw+6"
! type (mp_real) f1, f2, sum, td, tn
integer (mpiknd) f1(0:mpnw+6), f2(0:mpnw+6), sum(0:mpnw+6), td(0:mpnw+6), tn(0:mpnw+6)
! type (mp_real) t1, t2, t3, t4, rra
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), rra(0:mpnw+6)
integer mpnw1

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpt3(0:mpnw+6)
integer (mpiknd) mpt4(0:mpnw+6), mpt5(0:mpnw+6), mpt6(0:mpnw+6)

call mpinitwds (mpt1, mpnw+6-5)
call mpinitwds (mpt2, mpnw+6-5)
call mpinitwds (mpt3, mpnw+6-5)
call mpinitwds (mpt4, mpnw+6-5)
call mpinitwds (mpt5, mpnw+6-5)
call mpinitwds (mpt6, mpnw+6-5)
call mpinitwds (eps(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc1(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc2(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (f1(0:mpnw+6), mpnw+6-5)
call mpinitwds (f2(0:mpnw+6), mpnw+6-5)
call mpinitwds (sum(0:mpnw+6), mpnw+6-5)
call mpinitwds (td(0:mpnw+6), mpnw+6-5)
call mpinitwds (tn(0:mpnw+6), mpnw+6-5)
call mpinitwds (t1(0:mpnw+6), mpnw+6-5)
call mpinitwds (t2(0:mpnw+6), mpnw+6-5)
call mpinitwds (t3(0:mpnw+6), mpnw+6-5)
call mpinitwds (t4(0:mpnw+6), mpnw+6-5)
call mpinitwds (rra(0:mpnw+6), mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)

! ic1 = mpspacer (ss)
ic1 = mpspacer (ss(0:))
if (mpnw < mpnwm .or. ic1 < mpnw + 6) then
  write (6, 1)
1 format ('*** BESSELIR: Uninitialized or inadequately sized array')
!   call mpabrt (509)
  call mpabrt (509)
endif

! ic1 = mpwprecr (mppicon)
ic1 = mpwprecr (mppicon(0:))
if (ic1 < mpnw1) then
  write (6, 2) mpnw1
2 format ('*** BESSELIR: Pi must be initialized to',i6,' words by calling mpinit.')
!   call mpabrt (510)
  call mpabrt (510)
endif

!mp prec="mpnwm"

! tc2 = mprealdp (2.e0_mpdk)
call mprealdp (2.e0_mpdk, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), tc2(0:mpnwm+5), mpnwm)
ic2 = -mpnw1*mpnbt
! eps = tc2 ** ic2
call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), eps(0:mpnwm+5), mpnwm)

!mp prec="mpnw1"

!   If QQ is integer, call mpbesselinr; if qq < 0 and rr <= 0, then error.

! t1 = qq - anint (qq)
call mpnint (qq(0:), mpt1(0:mpnw+6), mpnw1)
call mpsub (qq(0:), mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
call mpeq (mpt2(0:mpnw+6), t1(0:mpnw+6), mpnw1)
! i0 = mpsgn (t1)
mpi1 = mpsgn (t1(0:mpnw+6))
i0 = mpi1
! i1 = mpsgn (qq)
mpi1 = mpsgn (qq(0:))
i1 = mpi1
! i2 = mpsgn (rr)
mpi1 = mpsgn (rr(0:))
i2 = mpi1
if (i0 == 0) then
!   d1 = dpreal (qq)
  call mpmdc (qq(0:), mpd1, mpi1, mpnw1)
  mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
  d1 = mpd1
  n1 = nint (d1)
!mp prec="mpnw"
!   call mpbesselinr (n1, rr, t3, mpnw)
  call mpbesselinr (n1, rr, t3, mpnw)
  goto 120
elseif (i1 < 0 .and. i2 <= 0) then
  write (mpldb, 3)
3 format ('*** BESSELIR: First argument < 0 and second argument <= 0')
!   call mpabrt (511)
  call mpabrt (511)
endif

!mp prec="mpnw1"
! rra = abs (rr)
call mpabs (rr(0:), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), rra(0:mpnw+6), mpnw1)
! d1 = dpreal (rra)
call mpmdc (rra(0:mpnw+6), mpd1, mpi1, mpnw1)
mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
d1 = mpd1

if (d1 < dfrac * mpnw1 * mpnbt) then
!   tn = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), tn(0:mpnw+6), mpnw1)
!   f1 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), f1(0:mpnw+6), mpnw1)
!   t1 = qq + f1
  call mpadd (qq(0:), f1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   call mpgammar (t1, f2, mpnw1)
  call mpgammar (t1, f2, mpnw1)
!   t1 = 0.25e0_mpdk * rra ** 2
  call mpnpwr (rra(0:mpnw+6), 2, mpt1(0:mpnw+6), mpnw1)
  call mpmuld (mpt1(0:mpnw+6), 0.25e0_mpdk, mpt2(0:mpnw+6), mpnw1)
  call mpeq (mpt2(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   td = f1 * f2
  call mpmul (f1(0:mpnw+6), f2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), td(0:mpnw+6), mpnw1)
!   sum = tn / td
  call mpdiv (tn(0:mpnw+6), td(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), sum(0:mpnw+6), mpnw1)

  do k = 1, itrmax
!     f1 = f1 * real (k, mpdk)
    call mpmuld (f1(0:mpnw+6), real (k, mpdk), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), f1(0:mpnw+6), mpnw1)
    d2 = real (k, mpdk)
!     t3 = mprealdp (d2)
    call mprealdp (d2, mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!     t4 = qq + t3
    call mpadd (qq(0:), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!     f2 = f2 * t4
    call mpmul (f2(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), f2(0:mpnw+6), mpnw1)
!     tn = t1 * tn
    call mpmul (t1(0:mpnw+6), tn(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), tn(0:mpnw+6), mpnw1)
!     td = f1 * f2
    call mpmul (f1(0:mpnw+6), f2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), td(0:mpnw+6), mpnw1)
!     t2 = tn / td
    call mpdiv (tn(0:mpnw+6), td(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!     sum = sum + t2
    call mpadd (sum(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), sum(0:mpnw+6), mpnw1)
!mp prec="mpnwm"
!     tc1 = abs (t2) - eps * abs (sum)
    call mpabs (t2(0:mpnw+6), mpt1(0:mpnw+6), mpnwm)
    call mpabs (sum(0:mpnw+6), mpt2(0:mpnw+6), mpnwm)
    call mpmul (eps(0:mpnwm+5), mpt2(0:mpnw+6), mpt3(0:mpnw+6), mpnwm)
    call mpsub (mpt1(0:mpnw+6), mpt3(0:mpnw+6), mpt4(0:mpnw+6), mpnwm)
    call mpeq (mpt4(0:mpnw+6), tc1(0:mpnwm+5), mpnwm)
!     ic1 = mpsgn (tc1)
    mpi1 = mpsgn (tc1(0:mpnwm+5))
    ic1 = mpi1
    if (ic1 < 0) goto 100
  enddo

  write (mpldb, 4)
4 format ('*** BESSELIR: Loop end error 1')
!   call mpabrt (512)
  call mpabrt (512)

100 continue

!mp prec="mpnw1"
!   t1 = 0.5e0_mpdk * rr
  call mpmuld (rr(0:), 0.5e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   t2 = t1 ** qq
  call mppower (t1(0:mpnw+6), qq(0:), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   t3 = sum * t2
  call mpmul (sum(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
else
!   sum = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), sum(0:mpnw+6), mpnw1)
!   t1 = 4.e0_mpdk * qq ** 2
  call mpnpwr (qq(0:), 2, mpt1(0:mpnw+6), mpnw1)
  call mpmuld (mpt1(0:mpnw+6), 4.e0_mpdk, mpt2(0:mpnw+6), mpnw1)
  call mpeq (mpt2(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   tn = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), tn(0:mpnw+6), mpnw1)
!   td = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), td(0:mpnw+6), mpnw1)

  do k = 1, itrmax
    i1 = 2 * k - 1
    d2 = real (i1, mpdk)
!     t2 = mprealdp (d2)
    call mprealdp (d2, mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!     t2 = t1 - t2 ** 2
    call mpnpwr (t2(0:mpnw+6), 2, mpt1(0:mpnw+6), mpnw1)
    call mpsub (t1(0:mpnw+6), mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
    call mpeq (mpt2(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!     tn = - tn * t2
    call mpmul (tn(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpneg (mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
    call mpeq (mpt2(0:mpnw+6), tn(0:mpnw+6), mpnw1)
    i2 = 8 * k
!     t2 = real (i2, mpdk) * rra
    call mpmuld (rra(0:mpnw+6), real (i2, mpdk), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!     td = td * t2
    call mpmul (td(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), td(0:mpnw+6), mpnw1)
!     t4 = tn / td
    call mpdiv (tn(0:mpnw+6), td(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!     sum = sum + t4
    call mpadd (sum(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), sum(0:mpnw+6), mpnw1)
!mp prec="mpnwm"
!     tc2 = abs (t4) - eps * abs (sum)
    call mpabs (t4(0:mpnw+6), mpt1(0:mpnw+6), mpnwm)
    call mpabs (sum(0:mpnw+6), mpt2(0:mpnw+6), mpnwm)
    call mpmul (eps(0:mpnwm+5), mpt2(0:mpnw+6), mpt3(0:mpnw+6), mpnwm)
    call mpsub (mpt1(0:mpnw+6), mpt3(0:mpnw+6), mpt4(0:mpnw+6), mpnwm)
    call mpeq (mpt4(0:mpnw+6), tc2(0:mpnwm+5), mpnwm)
!     ic2 = mpsgn (tc2)
    mpi1 = mpsgn (tc2(0:mpnwm+5))
    ic2 = mpi1
    if (ic2 < 0) goto 110
 enddo

write (mpldb, 5)
5 format ('*** BESSELIR: Loop end error 2')
! call mpabrt (513)
call mpabrt (513)

110 continue

!mp prec="mpnw1"
!   t1 = exp (rra)
  call mpexp (rra(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   t2 = 2.e0_mpdk * mppicon * rra
  call mpmuld (mppicon(0:), 2.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpmul (mpt1(0:mpnw+6), rra(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
  call mpeq (mpt2(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   t4 = sqrt (t2)
  call mpsqrt (t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!   t3 = sum * t1 / t4
  call mpmul (sum(0:mpnw+6), t1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpdiv (mpt1(0:mpnw+6), t4(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
  call mpeq (mpt2(0:mpnw+6), t3(0:mpnw+6), mpnw1)
endif

120 continue

!mp prec="mpnw"
! ss = t3
call mpeq (t3(0:mpnw+6), ss(0:), mpnw)

return
end subroutine mpbesselir

subroutine mpbesseljnr (nu, rr, ss, mpnw)

!   This evaluates the modified Bessel function BesselJ (NU,RR).
!   NU is an integer. The algorithm is DLMF formula 10.2.2 for modest RR,
!   and DLMF 10.17.3 for large RR, relative to precision.

implicit none
integer, intent(in):: mpnw
integer, intent(in):: nu
!mp dim1="0:"
! type (mp_real), intent(in):: rr
integer (mpiknd), intent(in):: rr(0:)
! type (mp_real), intent(out):: ss
integer (mpiknd), intent(out):: ss(0:)
integer, parameter:: itrmax = 1000000, ipx = 3
real (mpdknd), parameter:: dfrac1 = 0.4e0_mpdk, dfrac2 = 1.4e0_mpdk
integer i1, ic1, ic2, k, nua, n1
real (mpdknd) d1, d2
!mp dim1="0:mpnwm+5"
! type (mp_real) eps, tc1, tc2
integer (mpiknd) eps(0:mpnwm+5), tc1(0:mpnwm+5), tc2(0:mpnwm+5)
!mp dim1="0:ipx*mpnw+6"
! type (mp_real) f1, f2, sum1, sum2
integer (mpiknd) f1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), sum2(0:ipx*mpnw+6)
! type (mp_real) td1, td2, tn1, tn2
integer (mpiknd) td1(0:ipx*mpnw+6), td2(0:ipx*mpnw+6), tn1(0:ipx*mpnw+6), tn2(0:ipx*mpnw+6)
! type (mp_real) t1, t2, t3, t41
integer (mpiknd) t1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), t41(0:ipx*mpnw+6)
! type (mp_real) t42, t5, rra, rr2
integer (mpiknd) t42(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), rra(0:ipx*mpnw+6), rr2(0:ipx*mpnw+6)
integer mpnw1, mpnw2

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6)
integer (mpiknd) mpt4(0:ipx*mpnw+6), mpt5(0:ipx*mpnw+6), mpt6(0:ipx*mpnw+6)

call mpinitwds (mpt1, ipx*mpnw+6-5)
call mpinitwds (mpt2, ipx*mpnw+6-5)
call mpinitwds (mpt3, ipx*mpnw+6-5)
call mpinitwds (mpt4, ipx*mpnw+6-5)
call mpinitwds (mpt5, ipx*mpnw+6-5)
call mpinitwds (mpt6, ipx*mpnw+6-5)
call mpinitwds (eps(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc1(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc2(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (f1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (f2(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (sum1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (sum2(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (td1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (td2(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (tn1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (tn2(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t2(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t3(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t41(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t42(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t5(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (rra(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (rr2(0:ipx*mpnw+6), ipx*mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)
mpnw2 = min (ipx * mpnw + 1, mpnwx)

! ic1 = mpspacer (ss)
ic1 = mpspacer (ss(0:))
if (mpnw < mpnwm .or. ic1 < mpnw + 6) then
  write (6, 1)
1 format ('*** BESSELJNR: Uninitialized or inadequately sized array')
!   call mpabrt (514)
  call mpabrt (514)
endif

! ic1 = mpwprecr (mppicon)
ic1 = mpwprecr (mppicon(0:))
if (ic1 < mpnw2) then
  write (6, 2) mpnw2
2 format ('*** BESSELJNR: Pi must be initialized to',i6,' words by calling mpinit.')
!   call mpabrt (515)
  call mpabrt (515)
endif

!mp prec="mpnw1"

!   Check for RR = 0.

! ic1 = mpsgn (rr)
mpi1 = mpsgn (rr(0:))
ic1 = mpi1
if (ic1 == 0) then
  write (mpldb, 3)
3 format ('*** BESSELJNR: Second argument is zero')
!   call mpabrt (516)
  call mpabrt (516)
endif

nua = abs (nu)
! d1 = dpreal (rr)
call mpmdc (rr(0:), mpd1, mpi1, mpnw1)
mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
d1 = mpd1

if (d1 < dfrac1 * mpnw1 * mpnbt) then

!   In this case, perform computations with higher precision.

  mpnw2 = min (mpnw + nint (d1 * dfrac2 / mpnbt) + 1, mpnwx)
  if (mpnw2 > ipx*mpnw+1) then
    write (6, 4) mpnw2
4   format ('*** BESSELJNR: Inadequate working precision:',i6)
!     call mpabrt (517)
    call mpabrt (517)
  endif

!mp prec="mpnwm"
!   tc2 = mprealdp (2.e0_mpdk)
  call mprealdp (2.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnwm)
  call mpeq (mpt1(0:ipx*mpnw+6), tc2(0:mpnwm+5), mpnwm)
  ic2 = -mpnw2*mpnbt
!   eps = tc2 ** ic2
  call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:ipx*mpnw+6), mpnwm)
  call mpeq (mpt1(0:ipx*mpnw+6), eps(0:mpnwm+5), mpnwm)
!mp prec="mpnw2"
!   rra = abs (rr)
  call mpabs (rr(0:), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), rra(0:ipx*mpnw+6), mpnw2)
!   tn1 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), tn1(0:ipx*mpnw+6), mpnw2)
!   f1 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), f1(0:ipx*mpnw+6), mpnw2)
!   f2 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpnw2)
!   t1 = 0.25e0_mpdk * rra ** 2
  call mpnpwr (rra(0:ipx*mpnw+6), 2, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpmuld (mpt1(0:ipx*mpnw+6), 0.25e0_mpdk, mpt2(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt2(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw2)

  do k = 1, nua
!     f2 = real (k, mpdk) * f2
    call mpmuld (f2(0:ipx*mpnw+6), real (k, mpdk), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpnw2)
  enddo

!   td1 = f1 * f2
  call mpmul (f1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), td1(0:ipx*mpnw+6), mpnw2)
!   sum1 = tn1 / td1
  call mpdiv (tn1(0:ipx*mpnw+6), td1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), mpnw2)

  do k = 1, itrmax
!     f1 = f1 * real (k, mpdk)
    call mpmuld (f1(0:ipx*mpnw+6), real (k, mpdk), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), f1(0:ipx*mpnw+6), mpnw2)
    i1 = k + nua
!     f2 = f2 * real (i1, mpdk)
    call mpmuld (f2(0:ipx*mpnw+6), real (i1, mpdk), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpnw2)
!     tn1 = - t1 * tn1
    call mpmul (t1(0:ipx*mpnw+6), tn1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpneg (mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt2(0:ipx*mpnw+6), tn1(0:ipx*mpnw+6), mpnw2)
!     td1 = f1 * f2
    call mpmul (f1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), td1(0:ipx*mpnw+6), mpnw2)
!     t2 = tn1 / td1
    call mpdiv (tn1(0:ipx*mpnw+6), td1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!     sum1 = sum1 + t2
    call mpadd (sum1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), mpnw2)
!mp prec="mpnwm"
!     tc1 = abs (t2) - eps * abs (sum1)
    call mpabs (t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnwm)
    call mpabs (sum1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnwm)
    call mpmul (eps(0:mpnwm+5), mpt2(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpnwm)
    call mpsub (mpt1(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpt4(0:ipx*mpnw+6), mpnwm)
    call mpeq (mpt4(0:ipx*mpnw+6), tc1(0:mpnwm+5), mpnwm)
!     ic1 = mpsgn (tc1)
    mpi1 = mpsgn (tc1(0:mpnwm+5))
    ic1 = mpi1
    if (ic1 < 0) goto 100
  enddo

  write (mpldb, 5)
5 format ('*** BESSELJNR: Loop end error 1')
!   call mpabrt (518)
  call mpabrt (518)

100 continue

!mp prec="mpnw2"
!   t1 = 0.5e0_mpdk * rra
  call mpmuld (rra(0:ipx*mpnw+6), 0.5e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw2)
!   t3 = sum1 * t1 ** nua
  call mpnpwr (t1(0:ipx*mpnw+6), nua, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpmul (sum1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt2(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
else
!mp prec="mpnwm"
!   tc2 = mprealdp (2.e0_mpdk)
  call mprealdp (2.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnwm)
  call mpeq (mpt1(0:ipx*mpnw+6), tc2(0:mpnwm+5), mpnwm)
  ic2 = -mpnw1*mpnbt
!   eps = tc2 ** ic2
  call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:ipx*mpnw+6), mpnwm)
  call mpeq (mpt1(0:ipx*mpnw+6), eps(0:mpnwm+5), mpnwm)
!mp prec="mpnw1"
!   rra = abs (rr)
  call mpabs (rr(0:), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), rra(0:ipx*mpnw+6), mpnw1)
!   rr2 = rra ** 2
  call mpnpwr (rra(0:ipx*mpnw+6), 2, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), rr2(0:ipx*mpnw+6), mpnw1)
  i1 = 2 * nua
  d2 = real (i1, mpdk) ** 2
!   t1 = mprealdp (d2)
  call mprealdp (d2, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
!   tn1 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), tn1(0:ipx*mpnw+6), mpnw1)
!   t2 = t1 - tn1
  call mpsub (t1(0:ipx*mpnw+6), tn1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!   tn2 = t2 / 8.e0_mpdk
  call mpdivd (t2(0:ipx*mpnw+6), 8.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), tn2(0:ipx*mpnw+6), mpnw1)
!   td1 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), td1(0:ipx*mpnw+6), mpnw1)
!   td2 = rra
  call mpeq (rra(0:ipx*mpnw+6), td2(0:ipx*mpnw+6), mpnw1)
!   sum1 = tn1 / td1
  call mpdiv (tn1(0:ipx*mpnw+6), td1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), mpnw1)
!   sum2 = tn2 / td2
  call mpdiv (tn2(0:ipx*mpnw+6), td2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), sum2(0:ipx*mpnw+6), mpnw1)

  do k = 1, itrmax
    d1 = real (4*k-3, mpdk) ** 2
    d2 = real (4*k-1, mpdk) ** 2
!     t3 = t1 - d1
    call mpdmc (d1, 0, mpt2(0:ipx*mpnw+6), mpnw1)
    call mpsub (t1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw1)
!     t5 = t1 - d2
    call mpdmc (d2, 0, mpt2(0:ipx*mpnw+6), mpnw1)
    call mpsub (t1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), mpnw1)
!     t2 = t3 * t5
    call mpmul (t3(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!     tn1 = - tn1 * t2
    call mpmul (tn1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpneg (mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt2(0:ipx*mpnw+6), tn1(0:ipx*mpnw+6), mpnw1)
    d1 = real (2*k-1, mpdk) * real (2*k, mpdk) * 64.e0_mpdk
!     t2 = td1 * d1
    call mpmuld (td1(0:ipx*mpnw+6), d1, mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!     td1 = t2 * rr2
    call mpmul (t2(0:ipx*mpnw+6), rr2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), td1(0:ipx*mpnw+6), mpnw1)
!     t41 = tn1 / td1
    call mpdiv (tn1(0:ipx*mpnw+6), td1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t41(0:ipx*mpnw+6), mpnw1)
!     sum1 = sum1 + t41
    call mpadd (sum1(0:ipx*mpnw+6), t41(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), mpnw1)

    d1 = real (2*(2*k)-1, mpdk) ** 2
    d2 = real (2*(2*k+1)-1, mpdk) ** 2
!     t3 = t1 - d1
    call mpdmc (d1, 0, mpt2(0:ipx*mpnw+6), mpnw1)
    call mpsub (t1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw1)
!     t5 = t1 - d2
    call mpdmc (d2, 0, mpt2(0:ipx*mpnw+6), mpnw1)
    call mpsub (t1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), mpnw1)
!     t2 = t3 * t5
    call mpmul (t3(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!     tn2 = - tn2 * t2
    call mpmul (tn2(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpneg (mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt2(0:ipx*mpnw+6), tn2(0:ipx*mpnw+6), mpnw1)
    d1 = real (2*k, mpdk) * real (2*k+1, mpdk) * 64.e0_mpdk
!     t2 = td2 * d1
    call mpmuld (td2(0:ipx*mpnw+6), d1, mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!     td2 = t2 * rr2
    call mpmul (t2(0:ipx*mpnw+6), rr2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), td2(0:ipx*mpnw+6), mpnw1)
!     t42 = tn2 / td2
    call mpdiv (tn2(0:ipx*mpnw+6), td2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t42(0:ipx*mpnw+6), mpnw1)
!     sum2 = sum2 + t42
    call mpadd (sum2(0:ipx*mpnw+6), t42(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), sum2(0:ipx*mpnw+6), mpnw1)
!mp prec="mpnwm"
!     tc1 = abs (t41) - eps * abs (sum1)
    call mpabs (t41(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnwm)
    call mpabs (sum1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnwm)
    call mpmul (eps(0:mpnwm+5), mpt2(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpnwm)
    call mpsub (mpt1(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpt4(0:ipx*mpnw+6), mpnwm)
    call mpeq (mpt4(0:ipx*mpnw+6), tc1(0:mpnwm+5), mpnwm)
!     tc2 = abs (t42) - eps * abs (sum2)
    call mpabs (t42(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnwm)
    call mpabs (sum2(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnwm)
    call mpmul (eps(0:mpnwm+5), mpt2(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpnwm)
    call mpsub (mpt1(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpt4(0:ipx*mpnw+6), mpnwm)
    call mpeq (mpt4(0:ipx*mpnw+6), tc2(0:mpnwm+5), mpnwm)
!     ic1 = mpsgn (tc1)
    mpi1 = mpsgn (tc1(0:mpnwm+5))
    ic1 = mpi1
!     ic2 = mpsgn (tc2)
    mpi1 = mpsgn (tc2(0:mpnwm+5))
    ic2 = mpi1
    if (ic1 < 0 .and. ic2 < 0) goto 110
  enddo

  write (mpldb, 6)
6 format ('*** BESSELJNR: Loop end error 2')
!   call mpabrt (519)
  call mpabrt (519)

110 continue

!mp prec="mpnw1"
!   t1 = mppicon * 0.5e0_mpdk * real (nua, mpdk)
  call mpmuld (mppicon(0:), 0.5e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpmuld (mpt1(0:ipx*mpnw+6), real (nua, mpdk), mpt2(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt2(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
!   t2 = rra - t1
  call mpsub (rra(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!   t1 = mppicon * 0.25e0_mpdk
  call mpmuld (mppicon(0:), 0.25e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
!   t3 = t2 - t1
  call mpsub (t2(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw1)
!   t41 = cos (t3)
  call mpcssnr (t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t41(0:ipx*mpnw+6), mpnw1)
!   t42 = sin (t3)
  call mpcssnr (t3(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t42(0:ipx*mpnw+6), mpnw1)
!   t1 = t41 * sum1
  call mpmul (t41(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
!   t2 = t42 * sum2
  call mpmul (t42(0:ipx*mpnw+6), sum2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!   t5 = t1 - t2
  call mpsub (t1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), mpnw1)
!   t1 = mppicon * rra
  call mpmul (mppicon(0:), rra(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
!   t2 = mprealdp (2.e0_mpdk)
  call mprealdp (2.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!   t3 = t2 / t1
  call mpdiv (t2(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw1)
!   t1 = sqrt (t3)
  call mpsqrt (t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
!   t3 = t1 * t5
  call mpmul (t1(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw1)
endif

if (mod (nu, 2) /= 0) then
!   ic1 = mpsgn (rr)
  mpi1 = mpsgn (rr(0:))
  ic1 = mpi1
  if (nu < 0 .and. ic1 > 0 .or. nu > 0 .and. ic1 < 0) then
!     t3 = - t3
    call mpneg (t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw1)
  endif
endif

!mp prec="mpnw"
! ss = t3
call mpeq (t3(0:ipx*mpnw+6), ss(0:), mpnw)

return
end subroutine mpbesseljnr

subroutine mpbesseljr (qq, rr, ss, mpnw)

!   This evaluates the modified Bessel function BesselJ (QQ,RR) for QQ and RR
!   both MPR. The algorithm is DLMF formula 10.2.2 for modest RR,
!   and DLMF 10.17.3 for large RR, relative to precision.

implicit none
integer, intent(in):: mpnw
!mp dim1="0:"
! type (mp_real), intent(in):: qq, rr
integer (mpiknd), intent(in):: qq(0:), rr(0:)
! type (mp_real), intent(out):: ss
integer (mpiknd), intent(out):: ss(0:)
integer ic1, ic2, i1, i2, k, n1
real (mpdknd) d1, d2
integer, parameter:: itrmax = 1000000, ipx = 3
real (mpdknd), parameter:: dfrac1 = 0.4e0_mpdk, dfrac2 = 1.4e0_mpdk
!mp dim1="0:mpnwm+5"
! type (mp_real) eps, tc1, tc2
integer (mpiknd) eps(0:mpnwm+5), tc1(0:mpnwm+5), tc2(0:mpnwm+5)
!mp dim1="0:ipx*mpnw+6"
! type (mp_real) f1, f2, sum1, sum2
integer (mpiknd) f1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), sum2(0:ipx*mpnw+6)
! type (mp_real) td1, td2, tn1, tn2
integer (mpiknd) td1(0:ipx*mpnw+6), td2(0:ipx*mpnw+6), tn1(0:ipx*mpnw+6), tn2(0:ipx*mpnw+6)
! type (mp_real) t1, t2, t3, t4
integer (mpiknd) t1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), t4(0:ipx*mpnw+6)
! type (mp_real) t41, t42, t5
integer (mpiknd) t41(0:ipx*mpnw+6), t42(0:ipx*mpnw+6), t5(0:ipx*mpnw+6)
! type (mp_real) rra, rr2
integer (mpiknd) rra(0:ipx*mpnw+6), rr2(0:ipx*mpnw+6)
integer mpnw1, mpnw2

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6)
integer (mpiknd) mpt4(0:ipx*mpnw+6), mpt5(0:ipx*mpnw+6), mpt6(0:ipx*mpnw+6)

call mpinitwds (mpt1, ipx*mpnw+6-5)
call mpinitwds (mpt2, ipx*mpnw+6-5)
call mpinitwds (mpt3, ipx*mpnw+6-5)
call mpinitwds (mpt4, ipx*mpnw+6-5)
call mpinitwds (mpt5, ipx*mpnw+6-5)
call mpinitwds (mpt6, ipx*mpnw+6-5)
call mpinitwds (eps(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc1(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc2(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (f1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (f2(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (sum1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (sum2(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (td1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (td2(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (tn1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (tn2(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t2(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t3(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t4(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t41(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t42(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t5(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (rra(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (rr2(0:ipx*mpnw+6), ipx*mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)
mpnw2 = min (ipx*mpnw + 1, mpnwx)

! ic1 = mpspacer (ss)
ic1 = mpspacer (ss(0:))
if (mpnw < mpnwm .or. ic1 < mpnw + 6) then
  write (6, 1)
1 format ('*** BESSELJR: Uninitialized or inadequately sized array')
!   call mpabrt (520)
  call mpabrt (520)
endif

! ic1 = mpwprecr (mppicon)
ic1 = mpwprecr (mppicon(0:))
if (ic1 < mpnw2) then
  write (6, 2) mpnw2
2 format ('*** BESSELJR: Pi must be initialized to',i6,' words by calling mpinit.')
!   call mpabrt (521)
  call mpabrt (521)
endif

!mp prec="mpnw1"

!   If QQ is integer, call mpbesseljnr; if RR <= 0, then error.

! t2 = qq - anint (qq)
call mpnint (qq(0:), mpt1(0:ipx*mpnw+6), mpnw1)
call mpsub (qq(0:), mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnw1)
call mpeq (mpt2(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
! ic1 = mpsgn (t2)
mpi1 = mpsgn (t2(0:ipx*mpnw+6))
ic1 = mpi1
! ic2 = mpsgn (rr)
mpi1 = mpsgn (rr(0:))
ic2 = mpi1
if (ic1 == 0) then
!   d1 = dpreal (qq)
  call mpmdc (qq(0:), mpd1, mpi1, mpnw1)
  mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
  d1 = mpd1
  n1 = nint (d1)
!mp prec="mpnw"
!   call mpbesseljnr (n1, rr, t3, mpnw)
  call mpbesseljnr (n1, rr, t3, mpnw)
  goto 120
elseif (ic2 <= 0) then
  write (mpldb, 3)
3 format ('*** BESSELJR: Second argument <= 0')
!   call mpabrt (522)
  call mpabrt (522)
endif

! d1 = dpreal (rr)
call mpmdc (rr(0:), mpd1, mpi1, mpnw)
mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
d1 = mpd1

if (d1 < dfrac1 * mpnw1 * mpnbt) then

!   In this case, perform computations with higher precision.

  mpnw2 = min (mpnw + nint (d1 * dfrac2 / mpnbt) + 1, mpnwx)
  if (mpnw2 > ipx*mpnw+1) then
    write (6, 4) mpnw2
4   format ('*** BESSELJR: Inadequate working precision:',i6)
!     call mpabrt (523)
    call mpabrt (523)
  endif

!mp prec="mpnwm"
!   tc2 = mprealdp (2.e0_mpdk)
  call mprealdp (2.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnwm)
  call mpeq (mpt1(0:ipx*mpnw+6), tc2(0:mpnwm+5), mpnwm)
  ic2 = -mpnw2*mpnbt
!   eps = tc2 ** ic2
  call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:ipx*mpnw+6), mpnwm)
  call mpeq (mpt1(0:ipx*mpnw+6), eps(0:mpnwm+5), mpnwm)
!mp prec="mpnw2"
!   rra = abs (rr)
  call mpabs (rr(0:), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), rra(0:ipx*mpnw+6), mpnw2)
!   tn1 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), tn1(0:ipx*mpnw+6), mpnw2)
!   f1 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), f1(0:ipx*mpnw+6), mpnw2)
!   t1 = qq + f1
  call mpadd (qq(0:), f1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw2)
!   call mpgammar (t1, f2, mpnw2)
  call mpgammar (t1, f2, mpnw2)
!   t2 = rra ** 2
  call mpnpwr (rra(0:ipx*mpnw+6), 2, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!   t1 = t2 * 0.25e0_mpdk
  call mpmuld (t2(0:ipx*mpnw+6), 0.25e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw2)
!   td1 = f1 * f2
  call mpmul (f1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), td1(0:ipx*mpnw+6), mpnw2)
!   t2 = tn1 / td1
  call mpdiv (tn1(0:ipx*mpnw+6), td1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!   sum1 = t2
  call mpeq (t2(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), mpnw2)

  do k = 1, itrmax
!     f1 = f1 * real (k, mpdk)
    call mpmuld (f1(0:ipx*mpnw+6), real (k, mpdk), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), f1(0:ipx*mpnw+6), mpnw2)
    d2 = real (k, mpdk)
!     t3 = mprealdp (d2)
    call mprealdp (d2, mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!     t4 = qq + t3
    call mpadd (qq(0:), t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpnw2)
!     f2 = f2 * t4
    call mpmul (f2(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpnw2)
!     tn1 = - t1 * tn1
    call mpmul (t1(0:ipx*mpnw+6), tn1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpneg (mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt2(0:ipx*mpnw+6), tn1(0:ipx*mpnw+6), mpnw2)
!     td1 = f1 * f2
    call mpmul (f1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), td1(0:ipx*mpnw+6), mpnw2)
!     t2 = tn1 / td1
    call mpdiv (tn1(0:ipx*mpnw+6), td1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!     sum1 = sum1 + t2
    call mpadd (sum1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), mpnw2)
!mp prec="mpnwm"
!     tc1 = abs (t2) - eps * abs (sum1)
    call mpabs (t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnwm)
    call mpabs (sum1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnwm)
    call mpmul (eps(0:mpnwm+5), mpt2(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpnwm)
    call mpsub (mpt1(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpt4(0:ipx*mpnw+6), mpnwm)
    call mpeq (mpt4(0:ipx*mpnw+6), tc1(0:mpnwm+5), mpnwm)
!     ic1 = mpsgn (tc1)
    mpi1 = mpsgn (tc1(0:mpnwm+5))
    ic1 = mpi1
    if (ic1 < 0) goto 100
  enddo

  write (mpldb, 5)
5 format ('*** BESSELJR: Loop end error 1')
!   call mpabrt (524)
  call mpabrt (524)

100 continue

!mp prec="mpnw2"
!   t1 = rr * 0.5e0_mpdk
  call mpmuld (rr(0:), 0.5e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw2)
!   t2 = t1 ** qq
  call mppower (t1(0:ipx*mpnw+6), qq(0:), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!   t3 = sum1 * t2
  call mpmul (sum1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
else
  mpnw1 = min (mpnw + 1, mpnwx)
!mp prec="mpnwm"
!   tc2 = mprealdp (2.e0_mpdk)
  call mprealdp (2.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnwm)
  call mpeq (mpt1(0:ipx*mpnw+6), tc2(0:mpnwm+5), mpnwm)
  ic2 = -mpnw1*mpnbt
!   eps = tc2 ** ic2
  call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:ipx*mpnw+6), mpnwm)
  call mpeq (mpt1(0:ipx*mpnw+6), eps(0:mpnwm+5), mpnwm)
!mp prec="mpnw1"
!   rra = abs (rr)
  call mpabs (rr(0:), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), rra(0:ipx*mpnw+6), mpnw1)
!   rr2 = rra ** 2
  call mpnpwr (rra(0:ipx*mpnw+6), 2, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), rr2(0:ipx*mpnw+6), mpnw1)
!   t2 = qq ** 2
  call mpnpwr (qq(0:), 2, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!   t1 = t2 * 4.e0_mpdk
  call mpmuld (t2(0:ipx*mpnw+6), 4.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
!   tn1 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), tn1(0:ipx*mpnw+6), mpnw1)
!   t2 = t1 - tn1
  call mpsub (t1(0:ipx*mpnw+6), tn1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!   tn2 = t2 / 8.e0_mpdk
  call mpdivd (t2(0:ipx*mpnw+6), 8.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), tn2(0:ipx*mpnw+6), mpnw1)
!   td1 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), td1(0:ipx*mpnw+6), mpnw1)
!   td2 = rra
  call mpeq (rra(0:ipx*mpnw+6), td2(0:ipx*mpnw+6), mpnw1)
!   sum1 = tn1 / td1
  call mpdiv (tn1(0:ipx*mpnw+6), td1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), mpnw1)
!   sum2 = tn2 / td2
  call mpdiv (tn2(0:ipx*mpnw+6), td2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), sum2(0:ipx*mpnw+6), mpnw1)

  do k = 1, itrmax
    d1 = real (4*k-3, mpdk) ** 2
    d2 = real (4*k-1, mpdk) ** 2
!     t3 = t1 - d1
    call mpdmc (d1, 0, mpt2(0:ipx*mpnw+6), mpnw1)
    call mpsub (t1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw1)
!     t5 = t1 - d2
    call mpdmc (d2, 0, mpt2(0:ipx*mpnw+6), mpnw1)
    call mpsub (t1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), mpnw1)
!     t2 = t3 * t5
    call mpmul (t3(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!     tn1 = - tn1 * t2
    call mpmul (tn1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpneg (mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt2(0:ipx*mpnw+6), tn1(0:ipx*mpnw+6), mpnw1)
    d1 = real (2*k-1, mpdk) * real (2*k, mpdk) * 64.e0_mpdk
!     t2 = td1 * d1
    call mpmuld (td1(0:ipx*mpnw+6), d1, mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!     td1 = t2 * rr2
    call mpmul (t2(0:ipx*mpnw+6), rr2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), td1(0:ipx*mpnw+6), mpnw1)
!     t41 = tn1 / td1
    call mpdiv (tn1(0:ipx*mpnw+6), td1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t41(0:ipx*mpnw+6), mpnw1)
!     sum1 = sum1 + t41
    call mpadd (sum1(0:ipx*mpnw+6), t41(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), mpnw1)

    d1 = real (2*(2*k)-1, mpdk) ** 2
    d2 = real (2*(2*k+1) - 1, mpdk) ** 2
!     t3 = t1 - d1
    call mpdmc (d1, 0, mpt2(0:ipx*mpnw+6), mpnw1)
    call mpsub (t1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw1)
!     t5 = t1 - d2
    call mpdmc (d2, 0, mpt2(0:ipx*mpnw+6), mpnw1)
    call mpsub (t1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), mpnw1)
!     t2 = t3 * t5
    call mpmul (t3(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!     tn2 = - tn2 * t2
    call mpmul (tn2(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpneg (mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt2(0:ipx*mpnw+6), tn2(0:ipx*mpnw+6), mpnw1)
    d1 = real (2*k, mpdk) * real (2*k+1, mpdk) * 64.e0_mpdk
!     t2 = td2 * d1
    call mpmuld (td2(0:ipx*mpnw+6), d1, mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!     td2 = t2 * rr2
    call mpmul (t2(0:ipx*mpnw+6), rr2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), td2(0:ipx*mpnw+6), mpnw1)
!     t42 = tn2 / td2
    call mpdiv (tn2(0:ipx*mpnw+6), td2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t42(0:ipx*mpnw+6), mpnw1)
!     sum2 = sum2 + t42
    call mpadd (sum2(0:ipx*mpnw+6), t42(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), sum2(0:ipx*mpnw+6), mpnw1)
!mp prec="mpnwm"
!     tc1 = abs (t41) - eps * abs (sum1)
    call mpabs (t41(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnwm)
    call mpabs (sum1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnwm)
    call mpmul (eps(0:mpnwm+5), mpt2(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpnwm)
    call mpsub (mpt1(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpt4(0:ipx*mpnw+6), mpnwm)
    call mpeq (mpt4(0:ipx*mpnw+6), tc1(0:mpnwm+5), mpnwm)
!     tc2 = abs (t42) - eps * abs (sum2)
    call mpabs (t42(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnwm)
    call mpabs (sum2(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnwm)
    call mpmul (eps(0:mpnwm+5), mpt2(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpnwm)
    call mpsub (mpt1(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpt4(0:ipx*mpnw+6), mpnwm)
    call mpeq (mpt4(0:ipx*mpnw+6), tc2(0:mpnwm+5), mpnwm)
!     ic1 = mpsgn (tc1)
    mpi1 = mpsgn (tc1(0:mpnwm+5))
    ic1 = mpi1
!     ic2 = mpsgn (tc2)
    mpi1 = mpsgn (tc2(0:mpnwm+5))
    ic2 = mpi1
    if (ic1 < 0 .and. ic2 < 0) goto 110
  enddo

  write (mpldb, 6)
6 format ('*** BESSELJR: Loop end error 2')
!   call mpabrt (525)
  call mpabrt (525)

110 continue

!mp prec="mpnw1"
!   t2 = mppicon * qq
  call mpmul (mppicon(0:), qq(0:), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!   t1 = t2 * 0.5e0_mpdk
  call mpmuld (t2(0:ipx*mpnw+6), 0.5e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
!   t2 = rra - t1
  call mpsub (rra(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!   t1 = mppicon * 0.25e0_mpdk
  call mpmuld (mppicon(0:), 0.25e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
!   t3 = t2 - t1
  call mpsub (t2(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw1)
!   t41 = cos (t3)
  call mpcssnr (t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t41(0:ipx*mpnw+6), mpnw1)
!   t42 = sin (t3)
  call mpcssnr (t3(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t42(0:ipx*mpnw+6), mpnw1)
!   t1 = t41 * sum1
  call mpmul (t41(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
!   t2 = t42 * sum2
  call mpmul (t42(0:ipx*mpnw+6), sum2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!   t5 = t1 - t2
  call mpsub (t1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), mpnw1)
!   t1 = mppicon * rra
  call mpmul (mppicon(0:), rra(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
!   t2 = mprealdp (2.e0_mpdk)
  call mprealdp (2.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!   t3 = t2 / t1
  call mpdiv (t2(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw1)
!   t1 = sqrt (t3)
  call mpsqrt (t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
!   t3 = t1 * t5
  call mpmul (t1(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw1)
endif

120 continue

!mp prec="mpnw"
! ss = t3
call mpeq (t3(0:ipx*mpnw+6), ss(0:), mpnw)

return
end subroutine mpbesseljr

subroutine mpbesselknr (nu, rr, ss, mpnw)

!   This evaluates the modified Bessel function BesselK (NU,RR).
!   NU is an integer. The algorithm is DLMF formula 10.31.1 for modest RR,
!   and DLMF 10.40.2 for large RR, relative to precision.

implicit none
integer, intent(in):: mpnw
integer, intent(in):: nu
!mp dim1="0:"
! type (mp_real), intent(in):: rr
integer (mpiknd), intent(in):: rr(0:)
! type (mp_real), intent(out):: ss
integer (mpiknd), intent(out):: ss(0:)
integer, parameter:: itrmax = 1000000, ipx = 3
integer i1, ic1, ic2, k, nua, n1
real (mpdknd) d1
real (mpdknd), parameter:: dfrac1 = 0.4e0_mpdk, dfrac2 = 2.8e0_mpdk
!mp dim1="0:mpnwm+5"
! type (mp_real) eps, tc1, tc2
integer (mpiknd) eps(0:mpnwm+5), tc1(0:mpnwm+5), tc2(0:mpnwm+5)
!mp dim1="0:ipx*mpnw+6"
! type (mp_real) f1, f2, f3, f4, f5
integer (mpiknd) f1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), f3(0:ipx*mpnw+6), f4(0:ipx*mpnw+6), f5(0:ipx*mpnw+6)
! type (mp_real) sum1, sum2, sum3, td
integer (mpiknd) sum1(0:ipx*mpnw+6), sum2(0:ipx*mpnw+6), sum3(0:ipx*mpnw+6), td(0:ipx*mpnw+6)
! type (mp_real) tn, t1, t2, t3
integer (mpiknd) tn(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), t3(0:ipx*mpnw+6)
! type (mp_real) t4, t5, rra
integer (mpiknd) t4(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), rra(0:ipx*mpnw+6)
integer mpnw1, mpnw2

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6)
integer (mpiknd) mpt4(0:ipx*mpnw+6), mpt5(0:ipx*mpnw+6), mpt6(0:ipx*mpnw+6)

call mpinitwds (mpt1, ipx*mpnw+6-5)
call mpinitwds (mpt2, ipx*mpnw+6-5)
call mpinitwds (mpt3, ipx*mpnw+6-5)
call mpinitwds (mpt4, ipx*mpnw+6-5)
call mpinitwds (mpt5, ipx*mpnw+6-5)
call mpinitwds (mpt6, ipx*mpnw+6-5)
call mpinitwds (eps(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc1(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc2(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (f1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (f2(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (f3(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (f4(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (f5(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (sum1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (sum2(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (sum3(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (td(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (tn(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t2(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t3(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t4(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t5(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (rra(0:ipx*mpnw+6), ipx*mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)
mpnw2 = min (ipx*mpnw + 1, mpnwx)

! ic1 = mpspacer (ss)
ic1 = mpspacer (ss(0:))
if (mpnw < mpnwm .or. ic1 < mpnw + 6) then
  write (6, 1)
1 format ('*** BESSELKNR: Uninitialized or inadequately sized array')
!   call mpabrt (526)
  call mpabrt (526)
endif

! ic1 = mpwprecr (mppicon)
ic1 = mpwprecr (mppicon(0:))
! ic2 = mpwprecr (mpegammacon)
ic2 = mpwprecr (mpegammacon(0:))
if (ic1 < mpnw2 .or. ic2 < mpnw2) then
  write (6, 2) mpnw2
2 format ('*** BESSELKNR: Pi and Gamma must be initialized to',i6,' words by calling mpinit.')
!   call mpabrt (527)
  call mpabrt (527)
endif

!mp prec="mpnw1"

!   Check for RR = 0.

! ic1 = mpsgn (rr)
mpi1 = mpsgn (rr(0:))
ic1 = mpi1
if (ic1 == 0) then
  write (mpldb, 3)
3 format ('*** BESSELKNR: Second argument is zero')
!   call mpabrt (528)
  call mpabrt (528)
endif

nua = abs (nu)
! rra = abs (rr)
call mpabs (rr(0:), mpt1(0:ipx*mpnw+6), mpnw1)
call mpeq (mpt1(0:ipx*mpnw+6), rra(0:ipx*mpnw+6), mpnw1)
! d1 = dpreal (rra)
call mpmdc (rra(0:ipx*mpnw+6), mpd1, mpi1, mpnw1)
mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
d1 = mpd1

if (d1 < dfrac1 * mpnw1 * mpnbt) then

!   This algorithm requires higher working precision.

  mpnw2 = min (mpnw + nint (d1 * dfrac2 / mpnbt) + 1, mpnwx)
  if (mpnw2  > ipx*mpnw+1) then
    write (6, 4) mpnw2
4   format ('*** BESSELKNR: Inadequate working precision:',i6)
!     call mpabrt (529)
    call mpabrt (529)
  endif

!mp prec="mpnwm"
!   tc2 = mprealdp (2.e0_mpdk)
  call mprealdp (2.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnwm)
  call mpeq (mpt1(0:ipx*mpnw+6), tc2(0:mpnwm+5), mpnwm)
  ic2 = -mpnw2*mpnbt
!   eps = tc2 ** ic2
  call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:ipx*mpnw+6), mpnwm)
  call mpeq (mpt1(0:ipx*mpnw+6), eps(0:mpnwm+5), mpnwm)
!mp prec="mpnw2"
!   t2 = rra ** 2
  call mpnpwr (rra(0:ipx*mpnw+6), 2, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!   t1 = t2 * 0.25e0_mpdk
  call mpmuld (t2(0:ipx*mpnw+6), 0.25e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw2)
!   f1 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), f1(0:ipx*mpnw+6), mpnw2)
!   f2 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpnw2)
!   f3 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), f3(0:ipx*mpnw+6), mpnw2)
!   sum1 = mprealdp (0.e0_mpdk)
  call mprealdp (0.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), mpnw2)

  do k = 1, nua - 1
!     f1 = f1 * real (k, mpdk)
    call mpmuld (f1(0:ipx*mpnw+6), real (k, mpdk), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), f1(0:ipx*mpnw+6), mpnw2)
  enddo

  do k = 0, nua - 1
    if (k > 0) then
      i1 = nua - k
!       f1 = f1 / real (i1, mpdk)
      call mpdivd (f1(0:ipx*mpnw+6), real (i1, mpdk), mpt1(0:ipx*mpnw+6), mpnw2)
      call mpeq (mpt1(0:ipx*mpnw+6), f1(0:ipx*mpnw+6), mpnw2)
!       f2 = - t1 * f2
      call mpmul (t1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
      call mpneg (mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnw2)
      call mpeq (mpt2(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpnw2)
!       f3 = f3 * real (k, mpdk)
      call mpmuld (f3(0:ipx*mpnw+6), real (k, mpdk), mpt1(0:ipx*mpnw+6), mpnw2)
      call mpeq (mpt1(0:ipx*mpnw+6), f3(0:ipx*mpnw+6), mpnw2)
    endif
!     t3 = f1 * f2
    call mpmul (f1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!     t2 = t3 / f3
    call mpdiv (t3(0:ipx*mpnw+6), f3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!     sum1 = sum1 + t2
    call mpadd (sum1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), mpnw2)
  enddo

!   t2 = sum1 * 0.5e0_mpdk
  call mpmuld (sum1(0:ipx*mpnw+6), 0.5e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!   t3 = rra * 0.5e0_mpdk
  call mpmuld (rra(0:ipx*mpnw+6), 0.5e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!   t4 = t3 ** nua
  call mpnpwr (t3(0:ipx*mpnw+6), nua, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpnw2)
!   sum1 = t2 / t4
  call mpdiv (t2(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), mpnw2)
!   t5 = rra * 0.5e0_mpdk
  call mpmuld (rra(0:ipx*mpnw+6), 0.5e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), mpnw2)
!   t3 = log (t5)
  call mplog (t5(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
  ic1 = nua + 1
  d1 = (-1.e0_mpdk) ** ic1
!   t2 = t3 * d1
  call mpmuld (t3(0:ipx*mpnw+6), d1, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!   call mpbesselinr (nua, rra, t3, mpnw2)
  call mpbesselinr (nua, rra, t3, mpnw2)
!   sum2 = t2 * t3
  call mpmul (t2(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), sum2(0:ipx*mpnw+6), mpnw2)
!   f1 = - mpegammacon
  call mpneg (mpegammacon(0:), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), f1(0:ipx*mpnw+6), mpnw2)
!   f2 = f1
  call mpeq (f1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpnw2)
!   f3 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), f3(0:ipx*mpnw+6), mpnw2)
!   f4 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), f4(0:ipx*mpnw+6), mpnw2)
!   f5 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), f5(0:ipx*mpnw+6), mpnw2)

  do k = 1, nua
!     t2 = mprealdp (1.e0_mpdk)
    call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!     t3 = t2 / real (k, mpdk)
    call mpdivd (t2(0:ipx*mpnw+6), real (k, mpdk), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!     f2 = f2 + t3
    call mpadd (f2(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpnw2)
!     f5 = f5 * real (k, mpdk)
    call mpmuld (f5(0:ipx*mpnw+6), real (k, mpdk), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), f5(0:ipx*mpnw+6), mpnw2)
  enddo

!   t2 = f1 + f2
  call mpadd (f1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!   t3 = t2 * f3
  call mpmul (t2(0:ipx*mpnw+6), f3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!   t4 = f4 * f5
  call mpmul (f4(0:ipx*mpnw+6), f5(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpnw2)
!   sum3 = t3 / t4
  call mpdiv (t3(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), sum3(0:ipx*mpnw+6), mpnw2)

  do k = 1, itrmax
!     t2 = mprealdp (1.e0_mpdk)
    call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!     t3 = t2 / real (k, mpdk)
    call mpdivd (t2(0:ipx*mpnw+6), real (k, mpdk), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!     f1 = f1 + t3
    call mpadd (f1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), f1(0:ipx*mpnw+6), mpnw2)
    i1 = nua + k
!     t3 = t2 / real (i1, mpdk)
    call mpdivd (t2(0:ipx*mpnw+6), real (i1, mpdk), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!     f2 = f2 + t3
    call mpadd (f2(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpnw2)
!     f3 = t1 * f3
    call mpmul (t1(0:ipx*mpnw+6), f3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), f3(0:ipx*mpnw+6), mpnw2)
!     f4 = f4 * real (k, mpdk)
    call mpmuld (f4(0:ipx*mpnw+6), real (k, mpdk), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), f4(0:ipx*mpnw+6), mpnw2)
!     f5 = f5 * real (i1, mpdk)
    call mpmuld (f5(0:ipx*mpnw+6), real (i1, mpdk), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), f5(0:ipx*mpnw+6), mpnw2)
!     t2 = f1 + f2
    call mpadd (f1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!     t3 = t2 * f3
    call mpmul (t2(0:ipx*mpnw+6), f3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!     t4 = f4 * f5
    call mpmul (f4(0:ipx*mpnw+6), f5(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpnw2)
!     t2 = t3 / t4
    call mpdiv (t3(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!     sum3 = sum3 + t2
    call mpadd (sum3(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), sum3(0:ipx*mpnw+6), mpnw2)
!mp prec="mpnwm"
!     tc1 = abs (t2) - eps * abs (sum3)
    call mpabs (t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnwm)
    call mpabs (sum3(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnwm)
    call mpmul (eps(0:mpnwm+5), mpt2(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpnwm)
    call mpsub (mpt1(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpt4(0:ipx*mpnw+6), mpnwm)
    call mpeq (mpt4(0:ipx*mpnw+6), tc1(0:mpnwm+5), mpnwm)
!     ic1 = mpsgn (tc1)
    mpi1 = mpsgn (tc1(0:mpnwm+5))
    ic1 = mpi1
    if (ic1 < 0) goto 100
  enddo

  write (mpldb, 5)
5 format ('*** BESSELKNR: Loop end error 1')
!   call mpabrt (530)
  call mpabrt (530)

100 continue

!mp prec="mpnw2"
!   t2 = rra * 0.5e0_mpdk
  call mpmuld (rra(0:ipx*mpnw+6), 0.5e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!   t3 = t2 ** nua
  call mpnpwr (t2(0:ipx*mpnw+6), nua, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
  d1 = (-1.e0_mpdk) ** nua * 0.5e0_mpdk
!   t4 = t3 * d1
  call mpmuld (t3(0:ipx*mpnw+6), d1, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpnw2)
!   sum3 = t4 * sum3
  call mpmul (t4(0:ipx*mpnw+6), sum3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), sum3(0:ipx*mpnw+6), mpnw2)
!   t2 = sum1 + sum2
  call mpadd (sum1(0:ipx*mpnw+6), sum2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!   t3 = t2 + sum3
  call mpadd (t2(0:ipx*mpnw+6), sum3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
else
  mpnw1 = mpnw + 1
!mp prec="mpnwm"
!   tc2 = mprealdp (2.e0_mpdk)
  call mprealdp (2.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnwm)
  call mpeq (mpt1(0:ipx*mpnw+6), tc2(0:mpnwm+5), mpnwm)
  ic2 = -mpnw1*mpnbt
!   eps = tc2 ** ic2
  call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:ipx*mpnw+6), mpnwm)
  call mpeq (mpt1(0:ipx*mpnw+6), eps(0:mpnwm+5), mpnwm)
!mp prec="mpnw1"
!   sum1 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), mpnw1)
  d1 = 4.e0_mpdk * real (nua, mpdk) ** 2
!   t1 = mprealdp (d1)
  call mprealdp (d1, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
!   tn = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), tn(0:ipx*mpnw+6), mpnw1)
!   td = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), td(0:ipx*mpnw+6), mpnw1)

  do k = 1, itrmax
    d1 = real (2*k-1, mpdk)
!     t2 = mprealdp (d1)
    call mprealdp (d1, mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!     t3 = t2 ** 2
    call mpnpwr (t2(0:ipx*mpnw+6), 2, mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw1)
!     t2 = t1 - t3
    call mpsub (t1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!     tn = tn * t2
    call mpmul (tn(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), tn(0:ipx*mpnw+6), mpnw1)
    i1 = 8 * k
!     t2 = rra * real (i1, mpdk)
    call mpmuld (rra(0:ipx*mpnw+6), real (i1, mpdk), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!     td = td * t2
    call mpmul (td(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), td(0:ipx*mpnw+6), mpnw1)
!     t4 = tn / td
    call mpdiv (tn(0:ipx*mpnw+6), td(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpnw1)
!     sum1 = sum1 + t4
    call mpadd (sum1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), mpnw1)
!mp prec="mpnwm"
!     tc2 = abs (t4) - eps * abs (sum1)
    call mpabs (t4(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnwm)
    call mpabs (sum1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnwm)
    call mpmul (eps(0:mpnwm+5), mpt2(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpnwm)
    call mpsub (mpt1(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpt4(0:ipx*mpnw+6), mpnwm)
    call mpeq (mpt4(0:ipx*mpnw+6), tc2(0:mpnwm+5), mpnwm)
!     ic2 = mpsgn (tc2)
    mpi1 = mpsgn (tc2(0:mpnwm+5))
    ic2 = mpi1
    if (ic2 < 0) goto 110
  enddo

write (mpldb, 6)
6 format ('*** BESSELKNR: Loop end error 2')
! call mpabrt (531)
call mpabrt (531)

110 continue

!mp prec="mpnw1"
!   t1 = exp (rra)
  call mpexp (rra(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
!   t2 = rra * 2.e0_mpdk
  call mpmuld (rra(0:ipx*mpnw+6), 2.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!   t3 = mppicon / t2
  call mpdiv (mppicon(0:), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw1)
!   t4 = sqrt (t3)
  call mpsqrt (t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpnw1)
!   t2 = t4 / t1
  call mpdiv (t4(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!   t3 = t2 * sum1
  call mpmul (t2(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw1)
endif

!mp prec="mpnw"
! ic1 = mpsgn (rr)
mpi1 = mpsgn (rr(0:))
ic1 = mpi1
if (ic1 < 0 .and. mod (nu, 2) /= 0) then
!   t3 = - t3
  call mpneg (t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw)
endif
! ss = t3
call mpeq (t3(0:ipx*mpnw+6), ss(0:), mpnw)

return
end subroutine mpbesselknr

subroutine mpbesselkr (qq, rr, ss, mpnw)

!   This evaluates the Bessel function BesselK (QQ,RR) for QQ and RR
!   both MPR. This uses DLMF formula 10.27.4.

implicit none
integer, intent(in):: mpnw
!mp dim1="0:"
! type (mp_real), intent(in):: qq, rr
integer (mpiknd), intent(in):: qq(0:), rr(0:)
! type (mp_real), intent(out):: ss
integer (mpiknd), intent(out):: ss(0:)
integer, parameter:: itrmax = 1000000, ipx = 3
real (mpdknd), parameter:: dfrac1 = 0.4e0_mpdk, dfrac2 = 2.8e0_mpdk
integer ic1, ic2, i0, i1, i2, k, n1
real (mpdknd) d1
!mp dim1="0:mpnwm+5"
! type (mp_real) eps, tc1, tc2
integer (mpiknd) eps(0:mpnwm+5), tc1(0:mpnwm+5), tc2(0:mpnwm+5)
!mp dim1="0:ipx*mpnw+6"
! type (mp_real) t1, t2, t3, t4
integer (mpiknd) t1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), t4(0:ipx*mpnw+6)
! type (mp_real) rra, sum1, tn, td
integer (mpiknd) rra(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), tn(0:ipx*mpnw+6), td(0:ipx*mpnw+6)
integer mpnw1, mpnw2

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6)
integer (mpiknd) mpt4(0:ipx*mpnw+6), mpt5(0:ipx*mpnw+6), mpt6(0:ipx*mpnw+6)

call mpinitwds (mpt1, ipx*mpnw+6-5)
call mpinitwds (mpt2, ipx*mpnw+6-5)
call mpinitwds (mpt3, ipx*mpnw+6-5)
call mpinitwds (mpt4, ipx*mpnw+6-5)
call mpinitwds (mpt5, ipx*mpnw+6-5)
call mpinitwds (mpt6, ipx*mpnw+6-5)
call mpinitwds (eps(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc1(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc2(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (t1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t2(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t3(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t4(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (rra(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (sum1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (tn(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (td(0:ipx*mpnw+6), ipx*mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)
mpnw2 = ipx * mpnw + 1

! ic1 = mpspacer (ss)
ic1 = mpspacer (ss(0:))
if (mpnw < mpnwm .or. ic1 < mpnw + 6) then
  write (6, 1)
1 format ('*** BESSELKR: Uninitialized or inadequately sized array')
!   call mpabrt (532)
  call mpabrt (532)
endif

! ic1 = mpwprecr (mppicon)
ic1 = mpwprecr (mppicon(0:))
if (ic1 < mpnw2) then
  write (6, 2) mpnw2
2 format ('*** BESSELKR: Pi must be initialized to',i6,' words by calling mpinit.')
!   call mpabrt (533)
  call mpabrt (533)
endif

!mp prec="mpnw1"

!   If QQ is integer, call mpbesselknr; if qq < 0 and rr <= 0, then error.

! t1 = qq - anint (qq)
call mpnint (qq(0:), mpt1(0:ipx*mpnw+6), mpnw1)
call mpsub (qq(0:), mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnw1)
call mpeq (mpt2(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
! i0 = mpsgn (t1)
mpi1 = mpsgn (t1(0:ipx*mpnw+6))
i0 = mpi1
! i1 = mpsgn (qq)
mpi1 = mpsgn (qq(0:))
i1 = mpi1
! i2 = mpsgn (rr)
mpi1 = mpsgn (rr(0:))
i2 = mpi1
if (i0 == 0) then
!   d1 = dpreal (qq)
  call mpmdc (qq(0:), mpd1, mpi1, mpnw1)
  mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
  d1 = mpd1
  n1 = nint (d1)
!   call mpbesselknr (n1, rr, t1, mpnw)
  call mpbesselknr (n1, rr, t1, mpnw)
  goto 120
elseif (i1 < 0 .and. i2 <= 0) then
  write (mpldb, 3)
3 format ('*** BESSELKR: First argument < 0 and second argument <= 0')
!   call mpabrt (534)
  call mpabrt (534)
endif

! t1 = abs (rr)
call mpabs (rr(0:), mpt1(0:ipx*mpnw+6), mpnw1)
call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
! d1 = dpreal (t1)
call mpmdc (t1(0:ipx*mpnw+6), mpd1, mpi1, mpnw1)
mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
d1 = mpd1

if (d1 < dfrac1 * mpnw * mpnbt) then

!   This algorithm requires higher working precision.

  mpnw2 = min (mpnw + nint (d1 * dfrac2 / mpnbt) + 1, mpnwx)
  if (mpnw2 > ipx*mpnw+1) then
    write (6, 4) mpnw2
4   format ('*** BESSELKR: Inadequate working precision:',i6)
!     call mpabrt (535)
    call mpabrt (535)
  endif

!mp prec="mpnw2"
!   t1 = - qq
  call mpneg (qq(0:), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw2)
!   call mpbesselir (t1, rr, t2, mpnw2)
  call mpbesselir (t1, rr, t2, mpnw2)
!   call mpbesselir (qq, rr, t3, mpnw2)
  call mpbesselir (qq, rr, t3, mpnw2)
!   t4 = t2 - t3
  call mpsub (t2(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpnw2)
!   t1 = qq * mppicon
  call mpmul (qq(0:), mppicon(0:), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw2)
!   t3 = sin (t1)
  call mpcssnr (t1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!   t2 = t4 / t3
  call mpdiv (t4(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!   t3 = mppicon * t2
  call mpmul (mppicon(0:), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!   t1 = t3 * 0.5e0_mpdk
  call mpmuld (t3(0:ipx*mpnw+6), 0.5e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw2)
else
  mpnw1 = mpnw + 1
!mp prec="mpnwm"
!   tc2 = mprealdp (2.e0_mpdk)
  call mprealdp (2.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnwm)
  call mpeq (mpt1(0:ipx*mpnw+6), tc2(0:mpnwm+5), mpnwm)
  ic2 = -mpnw1*mpnbt
!   eps = tc2 ** ic2
  call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:ipx*mpnw+6), mpnwm)
  call mpeq (mpt1(0:ipx*mpnw+6), eps(0:mpnwm+5), mpnwm)
!mp prec="mpnw1"
!   rra = abs (rr)
  call mpabs (rr(0:), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), rra(0:ipx*mpnw+6), mpnw1)
!   sum1 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), mpnw1)
!   t1 = 4.e0_mpdk * qq ** 2
  call mpnpwr (qq(0:), 2, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpmuld (mpt1(0:ipx*mpnw+6), 4.e0_mpdk, mpt2(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt2(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
!   tn = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), tn(0:ipx*mpnw+6), mpnw1)
!   td = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), td(0:ipx*mpnw+6), mpnw1)

  do k = 1, itrmax
    d1 = real (2*k-1, mpdk)
!     t2 = mprealdp (d1)
    call mprealdp (d1, mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!     t3 = t2 ** 2
    call mpnpwr (t2(0:ipx*mpnw+6), 2, mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw1)
!     t2 = t1 - t3
    call mpsub (t1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!     tn = tn * t2
    call mpmul (tn(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), tn(0:ipx*mpnw+6), mpnw1)
    i1 = 8 * k
!     t2 = rra * real (i1, mpdk)
    call mpmuld (rra(0:ipx*mpnw+6), real (i1, mpdk), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!     td = td * t2
    call mpmul (td(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), td(0:ipx*mpnw+6), mpnw1)
!     t4 = tn / td
    call mpdiv (tn(0:ipx*mpnw+6), td(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpnw1)
!     sum1 = sum1 + t4
    call mpadd (sum1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), mpnw1)
!mp prec="mpnwm"
!     tc2 = abs (t4) - eps * abs (sum1)
    call mpabs (t4(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnwm)
    call mpabs (sum1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnwm)
    call mpmul (eps(0:mpnwm+5), mpt2(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpnwm)
    call mpsub (mpt1(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpt4(0:ipx*mpnw+6), mpnwm)
    call mpeq (mpt4(0:ipx*mpnw+6), tc2(0:mpnwm+5), mpnwm)
!     ic2 = mpsgn (tc2)
    mpi1 = mpsgn (tc2(0:mpnwm+5))
    ic2 = mpi1
    if (ic2 < 0) goto 110
  enddo

write (mpldb, 5)
5 format ('*** BESSELKNR: Loop end error')
! call mpabrt (536)
call mpabrt (536)

110 continue

!mp prec="mpnw1"
!   t1 = exp (rra)
  call mpexp (rra(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
!   t2 = rra * 2.e0_mpdk
  call mpmuld (rra(0:ipx*mpnw+6), 2.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!   t3 = mppicon / t2
  call mpdiv (mppicon(0:), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw1)
!   t4 = sqrt (t3)
  call mpsqrt (t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpnw1)
!   t2 = t4 / t1
  call mpdiv (t4(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!   t1 = t2 * sum1
  call mpmul (t2(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
endif

120 continue

!mp prec="mpnw"
! ss = t1
call mpeq (t1(0:ipx*mpnw+6), ss(0:), mpnw)
return
end subroutine mpbesselkr

subroutine mpbesselynr (nu, rr, ss, mpnw)

!   This evaluates the modified Bessel function BesselY (NU,RR).
!   NU is an integer. The algorithm is DLMF formula 10.8.1 for modest RR,
!   and DLMF 10.17.4 for large RR, relative to precision.

implicit none
integer, intent(in):: mpnw
integer, intent(in):: nu
!mp dim1="0:"
! type (mp_real), intent(in):: rr
integer (mpiknd), intent(in):: rr(0:)
! type (mp_real), intent(out):: ss
integer (mpiknd), intent(out):: ss(0:)
integer, parameter:: itrmax = 1000000, ipx = 3
real (mpdknd), parameter:: dfrac1 = 0.4e0_mpdk, dfrac2 = 1.4e0_mpdk
integer i1, ic1, ic2, k, nua, n1
real (mpdknd) d1, d2
!mp dim1="0:mpnwm+5"
! type (mp_real) eps, tc1, tc2
integer (mpiknd) eps(0:mpnwm+5), tc1(0:mpnwm+5), tc2(0:mpnwm+5)
!mp dim1="0:ipx*mpnw+6"
! type (mp_real) f1, f2, f3, f4, f5
integer (mpiknd) f1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), f3(0:ipx*mpnw+6), f4(0:ipx*mpnw+6), f5(0:ipx*mpnw+6)
! type (mp_real) rra, rr2, sum1, sum2
integer (mpiknd) rra(0:ipx*mpnw+6), rr2(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), sum2(0:ipx*mpnw+6)
! type (mp_real) sum3, td1, td2
integer (mpiknd) sum3(0:ipx*mpnw+6), td1(0:ipx*mpnw+6), td2(0:ipx*mpnw+6)
! type (mp_real) tn1, tn2, t1, t2
integer (mpiknd) tn1(0:ipx*mpnw+6), tn2(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6)
! type (mp_real) t3, t4, t41, t42, t5
integer (mpiknd) t3(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), t41(0:ipx*mpnw+6), t42(0:ipx*mpnw+6), t5(0:ipx*mpnw+6)
integer mpnw1, mpnw2

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6)
integer (mpiknd) mpt4(0:ipx*mpnw+6), mpt5(0:ipx*mpnw+6), mpt6(0:ipx*mpnw+6)

call mpinitwds (mpt1, ipx*mpnw+6-5)
call mpinitwds (mpt2, ipx*mpnw+6-5)
call mpinitwds (mpt3, ipx*mpnw+6-5)
call mpinitwds (mpt4, ipx*mpnw+6-5)
call mpinitwds (mpt5, ipx*mpnw+6-5)
call mpinitwds (mpt6, ipx*mpnw+6-5)
call mpinitwds (eps(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc1(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc2(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (f1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (f2(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (f3(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (f4(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (f5(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (rra(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (rr2(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (sum1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (sum2(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (sum3(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (td1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (td2(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (tn1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (tn2(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t2(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t3(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t4(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t41(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t42(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t5(0:ipx*mpnw+6), ipx*mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)
mpnw2 = ipx * mpnw + 1

! ic1 = mpspacer (ss)
ic1 = mpspacer (ss(0:))
if (mpnw < mpnwm .or. ic1 < mpnw + 6) then
  write (6, 1)
1 format ('*** BESSELYNR: Uninitialized or inadequately sized array')
!   call mpabrt (537)
  call mpabrt (537)
endif

! ic1 = mpwprecr (mppicon)
ic1 = mpwprecr (mppicon(0:))
! ic2 = mpwprecr (mpegammacon)
ic2 = mpwprecr (mpegammacon(0:))
if (ic1 < mpnw2 .or. ic2 < mpnw2) then
  write (6, 2) mpnw2
2 format ('*** BESSELYNR: Pi and Gamma must be initialized to',i6,' words by calling mpinit.')
!   call mpabrt (538)
  call mpabrt (538)
endif

!mp prec="mpnw1"

!   Check for RR = 0.

! ic1 = mpsgn (rr)
mpi1 = mpsgn (rr(0:))
ic1 = mpi1
if (ic1 == 0) then
  write (mpldb, 3)
3 format ('*** BESSELYNR: Argument is negative or too large')
!   call mpabrt (539)
  call mpabrt (539)
endif

nua = abs (nu)
! t1 = abs (rr)
call mpabs (rr(0:), mpt1(0:ipx*mpnw+6), mpnw1)
call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
! d1 = dpreal (t1)
call mpmdc (t1(0:ipx*mpnw+6), mpd1, mpi1, mpnw1)
mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
d1 = mpd1

if (d1 < dfrac1 * mpnw1 * mpnbt) then
  mpnw2 = min (mpnw + nint (d1 * dfrac2 / mpnbt) + 1, mpnwx)
  if (mpnw2 > ipx*mpnw+1) then
    write (6, 4) mpnw2
4   format ('*** BESSELYNR: Inadequate working precision:',i6)
!     call mpabrt (540)
    call mpabrt (540)
  endif

!mp prec="mpnwm"
!   tc2 = mprealdp (2.e0_mpdk)
  call mprealdp (2.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnwm)
  call mpeq (mpt1(0:ipx*mpnw+6), tc2(0:mpnwm+5), mpnwm)
  ic2 = -mpnw2*mpnbt
!   eps = tc2 ** ic2
  call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:ipx*mpnw+6), mpnwm)
  call mpeq (mpt1(0:ipx*mpnw+6), eps(0:mpnwm+5), mpnwm)
!mp prec="mpnw2"
!   rra = abs (rr)
  call mpabs (rr(0:), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), rra(0:ipx*mpnw+6), mpnw2)
!   t2 = rra ** 2
  call mpnpwr (rra(0:ipx*mpnw+6), 2, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!   t1 = t2 * 0.25e0_mpdk
  call mpmuld (t2(0:ipx*mpnw+6), 0.25e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw2)
!   f1 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), f1(0:ipx*mpnw+6), mpnw2)
!   f2 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpnw2)
!   f3 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), f3(0:ipx*mpnw+6), mpnw2)
!   sum1 = mprealdp (0.e0_mpdk)
  call mprealdp (0.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), mpnw2)

  do k = 1, nua - 1
!     f1 = f1 * real (k, mpdk)
    call mpmuld (f1(0:ipx*mpnw+6), real (k, mpdk), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), f1(0:ipx*mpnw+6), mpnw2)
  enddo

  do k = 0, nua - 1
    if (k > 0) then
      i1 = nua - k
!       f1 = f1 / real (i1, mpdk)
      call mpdivd (f1(0:ipx*mpnw+6), real (i1, mpdk), mpt1(0:ipx*mpnw+6), mpnw2)
      call mpeq (mpt1(0:ipx*mpnw+6), f1(0:ipx*mpnw+6), mpnw2)
!       f2 = t1 * f2
      call mpmul (t1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
      call mpeq (mpt1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpnw2)
!       f3 = f3 * real (k, mpdk)
      call mpmuld (f3(0:ipx*mpnw+6), real (k, mpdk), mpt1(0:ipx*mpnw+6), mpnw2)
      call mpeq (mpt1(0:ipx*mpnw+6), f3(0:ipx*mpnw+6), mpnw2)
    endif
!     t3 = f1 * f2
    call mpmul (f1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!     t2 = t3 / f3
    call mpdiv (t3(0:ipx*mpnw+6), f3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!     sum1 = sum1 + t2
    call mpadd (sum1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), mpnw2)
  enddo

!   t3 = rra * 0.5e0_mpdk
  call mpmuld (rra(0:ipx*mpnw+6), 0.5e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!   t4 = t3 ** nua
  call mpnpwr (t3(0:ipx*mpnw+6), nua, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpnw2)
!   sum1 = - sum1 / t4
  call mpdiv (sum1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpneg (mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt2(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), mpnw2)
!   t2 = rra * 0.5e0_mpdk
  call mpmuld (rra(0:ipx*mpnw+6), 0.5e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!   t3 = log (t2)
  call mplog (t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!   t2 = t3 * 2.e0_mpdk
  call mpmuld (t3(0:ipx*mpnw+6), 2.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!   call mpbesseljnr (nua, rra, t3, mpnw2)
  call mpbesseljnr (nua, rra, t3, mpnw2)
!   sum2 = t2 * t3
  call mpmul (t2(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), sum2(0:ipx*mpnw+6), mpnw2)

!   f1 = - mpegammacon
  call mpneg (mpegammacon(0:), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), f1(0:ipx*mpnw+6), mpnw2)
!   f2 = f1
  call mpeq (f1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpnw2)
!   f3 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), f3(0:ipx*mpnw+6), mpnw2)
!   f4 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), f4(0:ipx*mpnw+6), mpnw2)
!   f5 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), f5(0:ipx*mpnw+6), mpnw2)

  do k = 1, nua
!     t2 = mprealdp (1.e0_mpdk)
    call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!     t3 = t2 / real (k, mpdk)
    call mpdivd (t2(0:ipx*mpnw+6), real (k, mpdk), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!     f2 = f2 + t3
    call mpadd (f2(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpnw2)
!     f5 = f5 * real (k, mpdk)
    call mpmuld (f5(0:ipx*mpnw+6), real (k, mpdk), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), f5(0:ipx*mpnw+6), mpnw2)
  enddo

!   t2 = f1 + f2
  call mpadd (f1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!   t3 = t2 * f3
  call mpmul (t2(0:ipx*mpnw+6), f3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!   t4 = f4 * f5
  call mpmul (f4(0:ipx*mpnw+6), f5(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpnw2)
!   sum3 = t3 / t4
  call mpdiv (t3(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), sum3(0:ipx*mpnw+6), mpnw2)

  do k = 1, itrmax
!     t2 = mprealdp (1.e0_mpdk)
    call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!     t3 = t2 / real (k, mpdk)
    call mpdivd (t2(0:ipx*mpnw+6), real (k, mpdk), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!     f1 = f1 + t3
    call mpadd (f1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), f1(0:ipx*mpnw+6), mpnw2)
    i1 = nua + k
!     t3 = t2 / real (i1, mpdk)
    call mpdivd (t2(0:ipx*mpnw+6), real (i1, mpdk), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!     f2 = f2 + t3
    call mpadd (f2(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpnw2)
!     f3 = - t1 * f3
    call mpmul (t1(0:ipx*mpnw+6), f3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpneg (mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt2(0:ipx*mpnw+6), f3(0:ipx*mpnw+6), mpnw2)
!     f4 = f4 * real (k, mpdk)
    call mpmuld (f4(0:ipx*mpnw+6), real (k, mpdk), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), f4(0:ipx*mpnw+6), mpnw2)
!     f5 = f5 * real (i1, mpdk)
    call mpmuld (f5(0:ipx*mpnw+6), real (i1, mpdk), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), f5(0:ipx*mpnw+6), mpnw2)
!     t2 = f1 + f2
    call mpadd (f1(0:ipx*mpnw+6), f2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!     t3 = t2 * f3
    call mpmul (t2(0:ipx*mpnw+6), f3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!     t4 = f4 * f5
    call mpmul (f4(0:ipx*mpnw+6), f5(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpnw2)
!     t2 = t3 / t4
    call mpdiv (t3(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!     sum3 = sum3 + t2
    call mpadd (sum3(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), sum3(0:ipx*mpnw+6), mpnw2)
!mp prec="mpnwm"
!     tc1 = abs (t2) - eps * abs (sum3)
    call mpabs (t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnwm)
    call mpabs (sum3(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnwm)
    call mpmul (eps(0:mpnwm+5), mpt2(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpnwm)
    call mpsub (mpt1(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpt4(0:ipx*mpnw+6), mpnwm)
    call mpeq (mpt4(0:ipx*mpnw+6), tc1(0:mpnwm+5), mpnwm)
!     ic1 = mpsgn (tc1)
    mpi1 = mpsgn (tc1(0:mpnwm+5))
    ic1 = mpi1
    if (ic1 < 0) goto 100
  enddo

  write (mpldb, 5)
5 format ('*** BESSELYNR: Loop end error 1')
!   call mpabrt (541)
  call mpabrt (541)

100 continue

!mp prec="mpnw2"
!   t2 = rra * 0.5e0_mpdk
  call mpmuld (rra(0:ipx*mpnw+6), 0.5e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!   t3 = t2 ** nua
  call mpnpwr (t2(0:ipx*mpnw+6), nua, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!   sum3 = - t3 * sum3
  call mpmul (t3(0:ipx*mpnw+6), sum3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpneg (mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt2(0:ipx*mpnw+6), sum3(0:ipx*mpnw+6), mpnw2)
!   t2 = sum1 + sum2
  call mpadd (sum1(0:ipx*mpnw+6), sum2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!   t4 = t2 + sum3
  call mpadd (t2(0:ipx*mpnw+6), sum3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpnw2)
!   t2 = mppicon
  call mpeq (mppicon(0:), t2(0:ipx*mpnw+6), mpnw2)
!   t3 = t4 / t2
  call mpdiv (t4(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
else
  mpnw1 = min (mpnw + 1, mpnwx)
!mp prec="mpnwm"
!   tc2 = mprealdp (2.e0_mpdk)
  call mprealdp (2.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnwm)
  call mpeq (mpt1(0:ipx*mpnw+6), tc2(0:mpnwm+5), mpnwm)
  ic2 = -mpnw1*mpnbt
!   eps = tc2 ** ic2
  call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:ipx*mpnw+6), mpnwm)
  call mpeq (mpt1(0:ipx*mpnw+6), eps(0:mpnwm+5), mpnwm)
!mp prec="mpnw1"
!   rra = abs (rr)
  call mpabs (rr(0:), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), rra(0:ipx*mpnw+6), mpnw1)
!   rr2 = rra ** 2
  call mpnpwr (rra(0:ipx*mpnw+6), 2, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), rr2(0:ipx*mpnw+6), mpnw1)
  d2 = 4.e0_mpdk * real (nua, mpdk) ** 2
!   t1 = mprealdp (d2)
  call mprealdp (d2, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
!   tn1 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), tn1(0:ipx*mpnw+6), mpnw1)
!   t2 = t1 - tn1
  call mpsub (t1(0:ipx*mpnw+6), tn1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!   tn2 = t2 / 8.e0_mpdk
  call mpdivd (t2(0:ipx*mpnw+6), 8.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), tn2(0:ipx*mpnw+6), mpnw1)
!   td1 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), td1(0:ipx*mpnw+6), mpnw1)
!   td2 = rra
  call mpeq (rra(0:ipx*mpnw+6), td2(0:ipx*mpnw+6), mpnw1)
!   sum1 = tn1 / td1
  call mpdiv (tn1(0:ipx*mpnw+6), td1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), mpnw1)
!   sum2 = tn2 / td2
  call mpdiv (tn2(0:ipx*mpnw+6), td2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), sum2(0:ipx*mpnw+6), mpnw1)

  do k = 1, itrmax
    d1 = real (4*k-3, mpdk) ** 2
    d2 = real (4*k-1, mpdk) ** 2
!     t2 = mprealdp (d1)
    call mprealdp (d1, mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!     t3 = t1 - t2
    call mpsub (t1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw1)
!     t2 = mprealdp (d2)
    call mprealdp (d2, mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!     t5 = t1 - t2
    call mpsub (t1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), mpnw1)
!     t2 = t3 * t5
    call mpmul (t3(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!     tn1 = - tn1 * t2
    call mpmul (tn1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpneg (mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt2(0:ipx*mpnw+6), tn1(0:ipx*mpnw+6), mpnw1)
    d1 = real (2*k-1, mpdk) * real (2*k, mpdk) * 64.e0_mpdk
!     t2 = td1 * d1
    call mpmuld (td1(0:ipx*mpnw+6), d1, mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!     td1 = t2 * rr2
    call mpmul (t2(0:ipx*mpnw+6), rr2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), td1(0:ipx*mpnw+6), mpnw1)
!     t41 = tn1 / td1
    call mpdiv (tn1(0:ipx*mpnw+6), td1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t41(0:ipx*mpnw+6), mpnw1)
!     sum1 = sum1 + t41
    call mpadd (sum1(0:ipx*mpnw+6), t41(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), mpnw1)

    d1 = real (2*(2*k)-1, mpdk) ** 2
    d2 = real (2*(2*k+1)-1, mpdk) ** 2
!     t2 = mprealdp (d1)
    call mprealdp (d1, mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!     t3 = t1 - t2
    call mpsub (t1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw1)
!     t2 = mprealdp (d2)
    call mprealdp (d2, mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!     t5 = t1 - t2
    call mpsub (t1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), mpnw1)
!     t2 = t3 * t5
    call mpmul (t3(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!     tn2 = - tn2 * t2
    call mpmul (tn2(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpneg (mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt2(0:ipx*mpnw+6), tn2(0:ipx*mpnw+6), mpnw1)
    d1 = real (2*k, mpdk) * real (2*k+1, mpdk) * 64.e0_mpdk
!     t2 = td2 * d1
    call mpmuld (td2(0:ipx*mpnw+6), d1, mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!     td2 = t2 * rr2
    call mpmul (t2(0:ipx*mpnw+6), rr2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), td2(0:ipx*mpnw+6), mpnw1)
!     t42 = tn2 / td2
    call mpdiv (tn2(0:ipx*mpnw+6), td2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t42(0:ipx*mpnw+6), mpnw1)
!     sum2 = sum2 + t42
    call mpadd (sum2(0:ipx*mpnw+6), t42(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), sum2(0:ipx*mpnw+6), mpnw1)
!mp prec="mpnwm"
!     tc1 = abs (t41) - eps * abs (sum1)
    call mpabs (t41(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnwm)
    call mpabs (sum1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnwm)
    call mpmul (eps(0:mpnwm+5), mpt2(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpnwm)
    call mpsub (mpt1(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpt4(0:ipx*mpnw+6), mpnwm)
    call mpeq (mpt4(0:ipx*mpnw+6), tc1(0:mpnwm+5), mpnwm)
!     tc2 = abs (t42) - eps * abs (sum2)
    call mpabs (t42(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnwm)
    call mpabs (sum2(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnwm)
    call mpmul (eps(0:mpnwm+5), mpt2(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpnwm)
    call mpsub (mpt1(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpt4(0:ipx*mpnw+6), mpnwm)
    call mpeq (mpt4(0:ipx*mpnw+6), tc2(0:mpnwm+5), mpnwm)
!     ic1 = mpsgn (tc1)
    mpi1 = mpsgn (tc1(0:mpnwm+5))
    ic1 = mpi1
!     ic2 = mpsgn (tc2)
    mpi1 = mpsgn (tc2(0:mpnwm+5))
    ic2 = mpi1
    if (ic1 < 0 .and. ic2 < 0) goto 110
  enddo

  write (mpldb, 6)
6 format ('*** BESSELYNR: Loop end error 2')
!   call mpabrt (542)
  call mpabrt (542)

110 continue

!mp prec="mpnw1"
!   t1 = mppicon * 0.5e0_mpdk * real (nua, mpdk)
  call mpmuld (mppicon(0:), 0.5e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpmuld (mpt1(0:ipx*mpnw+6), real (nua, mpdk), mpt2(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt2(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
!   t2 = rra - t1
  call mpsub (rra(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!   t1 = mppicon * 0.25e0_mpdk
  call mpmuld (mppicon(0:), 0.25e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
!   t3 = t2 - t1
  call mpsub (t2(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw1)
!   t41 = cos (t3)
  call mpcssnr (t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t41(0:ipx*mpnw+6), mpnw1)
!   t42 = sin (t3)
  call mpcssnr (t3(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t42(0:ipx*mpnw+6), mpnw1)
!   t1 = t42 * sum1
  call mpmul (t42(0:ipx*mpnw+6), sum1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
!   t2 = t41 * sum2
  call mpmul (t41(0:ipx*mpnw+6), sum2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!   t5 = t1 + t2
  call mpadd (t1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), mpnw1)
!   t1 = mppicon * rra
  call mpmul (mppicon(0:), rra(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
!   t2 = mprealdp (2.e0_mpdk)
  call mprealdp (2.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!   t3 = t2 / t1
  call mpdiv (t2(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw1)
!   t1 = sqrt (t3)
  call mpsqrt (t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
!   t3 = t1 * t5
  call mpmul (t1(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw1)
endif

!mp prec="mpnw"
if (mod (nu, 2) /= 0) then
!   ic1 = mpsgn (rr)
  mpi1 = mpsgn (rr(0:))
  ic1 = mpi1
  if (nu < 0 .and. ic1 > 0 .or. nu > 0 .and. ic1 < 0) then
!     t3 = - t3
    call mpneg (t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw)
  endif
endif

! ss = t3
call mpeq (t3(0:ipx*mpnw+6), ss(0:), mpnw)

return
end subroutine mpbesselynr

subroutine mpbesselyr (qq, rr, ss, mpnw)

!   This evaluates the modified Bessel function BesselY (QQ,RR).
!   NU is an integer. The algorithm is DLMF formula 10.2.2.

implicit none
integer, intent(in):: mpnw
!mp dim1="0:"
! type (mp_real), intent(in):: qq, rr
integer (mpiknd), intent(in):: qq(0:), rr(0:)
! type (mp_real), intent(out):: ss
integer (mpiknd), intent(out):: ss(0:)
integer ic1, i0, i1, i2, n1
real (mpdknd) d1
!mp dim1="0:mpnw+6"
! type (mp_real) t1, t2, t3, t4
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6)
integer mpnw1

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpt3(0:mpnw+6)
integer (mpiknd) mpt4(0:mpnw+6), mpt5(0:mpnw+6), mpt6(0:mpnw+6)

call mpinitwds (mpt1, mpnw+6-5)
call mpinitwds (mpt2, mpnw+6-5)
call mpinitwds (mpt3, mpnw+6-5)
call mpinitwds (mpt4, mpnw+6-5)
call mpinitwds (mpt5, mpnw+6-5)
call mpinitwds (mpt6, mpnw+6-5)
call mpinitwds (t1(0:mpnw+6), mpnw+6-5)
call mpinitwds (t2(0:mpnw+6), mpnw+6-5)
call mpinitwds (t3(0:mpnw+6), mpnw+6-5)
call mpinitwds (t4(0:mpnw+6), mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)

! ic1 = mpspacer (ss)
ic1 = mpspacer (ss(0:))
if (mpnw < mpnwm .or. ic1 < mpnw + 6) then
  write (6, 1)
1 format ('*** BESSELYNR: Uninitialized or inadequately sized array')
!   call mpabrt (543)
  call mpabrt (543)
endif

!mp prec="mpnw1"

!   If QQ is integer, call mpbesselynr; if qq < 0 and rr <= 0, then error.

! t1 = qq - anint (qq)
call mpnint (qq(0:), mpt1(0:mpnw+6), mpnw1)
call mpsub (qq(0:), mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
call mpeq (mpt2(0:mpnw+6), t1(0:mpnw+6), mpnw1)
! i0 = mpsgn (t1)
mpi1 = mpsgn (t1(0:mpnw+6))
i0 = mpi1
! i1 = mpsgn (qq)
mpi1 = mpsgn (qq(0:))
i1 = mpi1
! i2 = mpsgn (rr)
mpi1 = mpsgn (rr(0:))
i2 = mpi1
if (i0 == 0) then
!   d1 = dpreal (qq)
  call mpmdc (qq(0:), mpd1, mpi1, mpnw1)
  mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
  d1 = mpd1
  n1 = nint (d1)
!   call mpbesselynr (n1, rr, t1, mpnw)
  call mpbesselynr (n1, rr, t1, mpnw)
  goto 120
elseif (i1 < 0 .and. i2 <= 0) then
  write (mpldb, 2)
2 format ('*** BESSELYR: First argument < 0 and second argument <= 0')
!   call mpabrt (544)
  call mpabrt (544)
endif

!mp prec="mpnw1"
! t1 = qq * mppicon
call mpmul (qq(0:), mppicon(0:), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
! t2 = cos (t1)
call mpcssnr (t1(0:mpnw+6), mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! t3 = sin (t1)
call mpcssnr (t1(0:mpnw+6), mpt2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
! call mpbesseljr (qq, rr, t4, mpnw1)
call mpbesseljr (qq, rr, t4, mpnw1)
! t1 = t4 * t2
call mpmul (t4(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
! t2 = - qq
call mpneg (qq(0:), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! call mpbesseljr (t2, rr, t4, mpnw1)
call mpbesseljr (t2, rr, t4, mpnw1)
! t2 = t1 - t4
call mpsub (t1(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! t1 = t2 / t3
call mpdiv (t2(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)

120 continue

!mp prec="mpnw"
! ss = t1
call mpeq (t1(0:mpnw+6), ss(0:), mpnw)

return
end subroutine mpbesselyr

subroutine mpdigammabe (nb1, nb2, berne, x, y, mpnw)

!  This evaluates the digamma function, using asymptotic formula DLMF 5.11.2:
!  dig(x) ~ log(x) - 1/(2*x) - Sum_{k=1}^inf B[2k] / (2*k*x^(2*k)).
!  Before using this formula, the recursion dig(x+1) = dig(x) + 1/x is used
!  to shift the argument up by IQQ, where IQQ is set based on mpnw below.
!  The array berne contains precomputed even Bernoulli numbers (see MPBERNER
!  above). Its dimensions must be as shown below. NB2 must be greater than
!  1.4 x precision in decimal digits.

implicit none
integer, intent(in):: nb1, mpnw
integer, intent (in):: nb2
!mp dim1="0:nb1+5"
! type (mp_real), intent(in):: berne(1:nb2)
integer (mpiknd), intent(in):: berne(0:nb1+5,1:nb2)
!mp dim1="0:"
! type (mp_real), intent(in):: x
integer (mpiknd), intent(in):: x(0:)
! type (mp_real), intent(out):: y
integer (mpiknd), intent(out):: y(0:)
real (mpdknd), parameter:: dber = 0.45e0_mpdk, dfrac = 0.12e0_mpdk
integer k, i1, i2, ic1, ic2, iqq, n1
real (mpdknd) d1, d2
!mp dim1="0:mpnwm+5"
! type (mp_real) eps, tc1, tc2
integer (mpiknd) eps(0:mpnwm+5), tc1(0:mpnwm+5), tc2(0:mpnwm+5)
!mp dim1="0:mpnw+6"
! type (mp_real) f1, sum1, sum2, t1
integer (mpiknd) f1(0:mpnw+6), sum1(0:mpnw+6), sum2(0:mpnw+6), t1(0:mpnw+6)
! type (mp_real) t2, t3, t4, t5, xq
integer (mpiknd) t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), t5(0:mpnw+6), xq(0:mpnw+6)
integer mpnw1

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpt3(0:mpnw+6)
integer (mpiknd) mpt4(0:mpnw+6), mpt5(0:mpnw+6), mpt6(0:mpnw+6)

call mpinitwds (mpt1, mpnw+6-5)
call mpinitwds (mpt2, mpnw+6-5)
call mpinitwds (mpt3, mpnw+6-5)
call mpinitwds (mpt4, mpnw+6-5)
call mpinitwds (mpt5, mpnw+6-5)
call mpinitwds (mpt6, mpnw+6-5)
call mpinitwds (eps(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc1(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc2(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (f1(0:mpnw+6), mpnw+6-5)
call mpinitwds (sum1(0:mpnw+6), mpnw+6-5)
call mpinitwds (sum2(0:mpnw+6), mpnw+6-5)
call mpinitwds (t1(0:mpnw+6), mpnw+6-5)
call mpinitwds (t2(0:mpnw+6), mpnw+6-5)
call mpinitwds (t3(0:mpnw+6), mpnw+6-5)
call mpinitwds (t4(0:mpnw+6), mpnw+6-5)
call mpinitwds (t5(0:mpnw+6), mpnw+6-5)
call mpinitwds (xq(0:mpnw+6), mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)

! ic1 = mpspacer (y)
ic1 = mpspacer (y(0:))
if (mpnw < mpnwm .or. ic1 < mpnw + 6) then
  write (6, 1)
1 format ('*** DIGAMMABE: Uninitialized or inadequately sized array')
!   call mpabrt (545)
  call mpabrt (545)
endif


iqq = dfrac * mpnw1 * mpnbt

!mp prec="mpnwm"

! tc2 = mprealdp (2.e0_mpdk)
call mprealdp (2.e0_mpdk, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), tc2(0:mpnwm+5), mpnwm)
ic2 = -mpnw1*mpnbt
! eps = tc2 ** ic2
call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), eps(0:mpnwm+5), mpnwm)

!mp prec="mpnw1"

! f1 = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), f1(0:mpnw+6), mpnw1)

!   Check if argument is less than or equal to 0 -- undefined.

! ic1 = mpsgn (x)
mpi1 = mpsgn (x(0:))
ic1 = mpi1
if (ic1 <= 0) then
  write (mpldb, 2)
2 format ('*** DIGAMMABE: Argument <= 0')
!   call mpabrt (546)
  call mpabrt (546)
endif

!   Check if berne array has been initialized.

! d1 = dpreal (berne(1))
call mpmdc (berne(0:nb1+5,1), mpd1, mpi1, mpnw1)
mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
d1 = mpd1
if (abs (d1 - 1.e0_mpdk / 6.e0_mpdk) > mprdfz .or. nb2 < int (dber * mpnbt * mpnw)) then
  write (mpldb, 3) int (dber * mpnbt * mpnw)
3 format ('*** DIGAMMABE: Array of even Bernoulli coefficients must be initialized'/ &
   'with at least',i8,' entries using BERNE or BERNER.')
!   call mpabrt (547)
  call mpabrt (547)
endif

! sum1 = mprealdp (0.e0_mpdk)
call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), sum1(0:mpnw+6), mpnw1)
! sum2 = mprealdp (0.e0_mpdk)
call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), sum2(0:mpnw+6), mpnw1)
! xq = x + real (iqq, mpdk)
call mpdmc (real (iqq, mpdk), 0, mpt2(0:mpnw+6), mpnw1)
call mpadd (x(0:), mpt2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), xq(0:mpnw+6), mpnw1)

do k = 0, iqq - 1
  d2 = real (k, mpdk)
!   t1 = mprealdp (d2)
  call mprealdp (d2, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   t2 = x + t1
  call mpadd (x(0:), t1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   t3 = f1 / t2
  call mpdiv (f1(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!   sum1 = sum1 + t3
  call mpadd (sum1(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), sum1(0:mpnw+6), mpnw1)
enddo

! t1 = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
! t2 = xq ** 2
call mpnpwr (xq(0:mpnw+6), 2, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)

do k = 1, nb2
!   t1 = t1 * t2
  call mpmul (t1(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   t4 = t1 * 2.e0_mpdk * real (k, mpdk)
  call mpmuld (t1(0:mpnw+6), 2.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpmuld (mpt1(0:mpnw+6), real (k, mpdk), mpt2(0:mpnw+6), mpnw1)
  call mpeq (mpt2(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!   t3 = berne(k) / t4
  call mpdiv (berne(0:nb1+5,k), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!   sum2 = sum2 + t3
  call mpadd (sum2(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), sum2(0:mpnw+6), mpnw1)
!mp prec="mpnwm"
!   tc1 = abs (t3) - eps * abs (sum2)
  call mpabs (t3(0:mpnw+6), mpt1(0:mpnw+6), mpnwm)
  call mpabs (sum2(0:mpnw+6), mpt2(0:mpnw+6), mpnwm)
  call mpmul (eps(0:mpnwm+5), mpt2(0:mpnw+6), mpt3(0:mpnw+6), mpnwm)
  call mpsub (mpt1(0:mpnw+6), mpt3(0:mpnw+6), mpt4(0:mpnw+6), mpnwm)
  call mpeq (mpt4(0:mpnw+6), tc1(0:mpnwm+5), mpnwm)
!   ic1 = mpsgn (tc1)
  mpi1 = mpsgn (tc1(0:mpnwm+5))
  ic1 = mpi1
  if (ic1 < 0) goto 110
enddo

write (mpldb, 4)
4 format ('*** DIGAMMABE: Loop end error: Increase NB2')
! call mpabrt (548)
call mpabrt (548)

110 continue

!mp prec="mpnw1"
! t1 = - sum1
call mpneg (sum1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
! t2 = log (xq)
call mplog (xq(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! t3 = t1 + t2
call mpadd (t1(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
! t4 = xq * 2.e0_mpdk
call mpmuld (xq(0:mpnw+6), 2.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
! t5 = f1 / t4
call mpdiv (f1(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
! t2 = t3 - t5
call mpsub (t3(0:mpnw+6), t5(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! t1 = t2 - sum2
call mpsub (t2(0:mpnw+6), sum2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!mp prec="mpnw"
! y = t1
call mpeq (t1(0:mpnw+6), y(0:), mpnw)

return
end subroutine mpdigammabe

subroutine mperfr (z, terf, mpnw)

!   This evaluates the erf function, using a combination of two series.
!   In particular, the algorithm is (where B = (mpnw + 1) * mpnbt, and
!   dcon is a constant defined below):

!   if (z == 0) then
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
!mp dim1="0:"
! type (mp_real), intent(in):: z
integer (mpiknd), intent(in):: z(0:)
! type (mp_real), intent(out):: terf
integer (mpiknd), intent(out):: terf(0:)
integer, parameter:: itrmx = 100000
real (mpdknd), parameter:: al2 = 0.69314718055994530941723212145817657e0_mpdk
real (mpdknd), parameter:: dcon = 100.e0_mpdk
integer ic1, ic2, ic3, ic4, k, nbt, n1
real (mpdknd) d1, d2
!mp dim1="0:mpnwm+5"
! type (mp_real) eps, tc1, tc2
integer (mpiknd) eps(0:mpnwm+5), tc1(0:mpnwm+5), tc2(0:mpnwm+5)
!mp dim1="0:mpnw+6"
! type (mp_real) t1, t2, t3, t4, t5
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), t5(0:mpnw+6)
! type (mp_real) t6, t7, z2
integer (mpiknd) t6(0:mpnw+6), t7(0:mpnw+6), z2(0:mpnw+6)
integer mpnw1

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpt3(0:mpnw+6)
integer (mpiknd) mpt4(0:mpnw+6), mpt5(0:mpnw+6), mpt6(0:mpnw+6)

call mpinitwds (mpt1, mpnw+6-5)
call mpinitwds (mpt2, mpnw+6-5)
call mpinitwds (mpt3, mpnw+6-5)
call mpinitwds (mpt4, mpnw+6-5)
call mpinitwds (mpt5, mpnw+6-5)
call mpinitwds (mpt6, mpnw+6-5)
call mpinitwds (eps(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc1(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc2(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (t1(0:mpnw+6), mpnw+6-5)
call mpinitwds (t2(0:mpnw+6), mpnw+6-5)
call mpinitwds (t3(0:mpnw+6), mpnw+6-5)
call mpinitwds (t4(0:mpnw+6), mpnw+6-5)
call mpinitwds (t5(0:mpnw+6), mpnw+6-5)
call mpinitwds (t6(0:mpnw+6), mpnw+6-5)
call mpinitwds (t7(0:mpnw+6), mpnw+6-5)
call mpinitwds (z2(0:mpnw+6), mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)

! ic1 = mpspacer (terf)
ic1 = mpspacer (terf(0:))
if (mpnw < mpnwm .or. ic1 < mpnw + 6) then
  write (6, 1)
1 format ('*** ERFR: Uninitialized or inadequately sized array')
!   call mpabrt (549)
  call mpabrt (549)
endif

!mp prec="mpnwm"

! tc2 = mprealdp (2.e0_mpdk)
call mprealdp (2.e0_mpdk, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), tc2(0:mpnwm+5), mpnwm)
ic2 = -mpnw1*mpnbt
! eps = tc2 ** ic2
call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), eps(0:mpnwm+5), mpnwm)

!mp prec="mpnw1"

nbt = mpnbt
d1 = aint (1.e0_mpdk + sqrt (nbt * al2))
d2 = aint (nbt / dcon + 8.e0_mpdk)
! t1 = mprealdp (d1)
call mprealdp (d1, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
! t2 = mprealdp (d2)
call mprealdp (d2, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! t3 = z - t1
call mpsub (z(0:), t1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
! t4 = z + t1
call mpadd (z(0:), t1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
! t5 = z - t2
call mpsub (z(0:), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
! ic1 = mpsgn (z)
mpi1 = mpsgn (z(0:))
ic1 = mpi1
! ic2 = mpsgn (t3)
mpi1 = mpsgn (t3(0:mpnw+6))
ic2 = mpi1
! ic3 = mpsgn (t4)
mpi1 = mpsgn (t4(0:mpnw+6))
ic3 = mpi1
! ic4 = mpsgn (t5)
mpi1 = mpsgn (t5(0:mpnw+6))
ic4 = mpi1

if (ic1 == 0) then
!   terf = mprealdp (0.e0_mpdk)
  call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), terf(0:), mpnw1)
elseif (ic2 > 0) then
!   terf = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), terf(0:), mpnw1)
elseif (ic3 < 0) then
!   terf = mprealdp (-1.e0_mpdk)
  call mprealdp (-1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), terf(0:), mpnw1)
elseif (ic4 < 0) then
!   z2 = z ** 2
  call mpnpwr (z(0:), 2, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), z2(0:mpnw+6), mpnw1)
!   t1 = mprealdp (0.e0_mpdk)
  call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   t2 = z
  call mpeq (z(0:), t2(0:mpnw+6), mpnw1)
!   t3 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!   t5 = mprealdp (1.e2_mpdk)
  call mprealdp (1.e2_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)

  do k = 0, itrmx
    if (k > 0) then
!       t6 = z2 * 2.e0_mpdk
      call mpmuld (z2(0:mpnw+6), 2.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
      call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!       t2 = t6 * t2
      call mpmul (t6(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
      call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
      d1 = real (2*k+1, mpdk)
!       t3 = t3 * d1
      call mpmuld (t3(0:mpnw+6), d1, mpt1(0:mpnw+6), mpnw1)
      call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
    endif

!     t4 = t2 / t3
    call mpdiv (t2(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!     t1 = t1 + t4
    call mpadd (t1(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!     t6 = t4 / t1
    call mpdiv (t4(0:mpnw+6), t1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!mp prec="mpnwm"
!     tc1 = t6 - eps
    call mpsub (t6(0:mpnw+6), eps(0:mpnwm+5), mpt1(0:mpnw+6), mpnwm)
    call mpeq (mpt1(0:mpnw+6), tc1(0:mpnwm+5), mpnwm)
!     tc2 = t5 - t6
    call mpsub (t5(0:mpnw+6), t6(0:mpnw+6), mpt1(0:mpnw+6), mpnwm)
    call mpeq (mpt1(0:mpnw+6), tc2(0:mpnwm+5), mpnwm)
!     ic1 = mpsgn (tc1)
    mpi1 = mpsgn (tc1(0:mpnwm+5))
    ic1 = mpi1
!     ic2 = mpsgn (tc2)
    mpi1 = mpsgn (tc2(0:mpnwm+5))
    ic2 = mpi1
    if (ic1 < 0 .or. ic2 < 0) goto 120
!mp prec="mpnw1"
!     t5 = t6
    call mpeq (t6(0:mpnw+6), t5(0:mpnw+6), mpnw1)
  enddo

write (mpldb, 2)
2 format ('*** ERFR: End loop error 1')
! call mpabrt (550)
call mpabrt (550)

120 continue

!mp prec="mpnw1"
!   t3 = t1 * 2.e0_mpdk
  call mpmuld (t1(0:mpnw+6), 2.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!   t4 = sqrt (mppicon)
  call mpsqrt (mppicon(0:), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!   t5 = exp (z2)
  call mpexp (z2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
!   t6 = t4 * t5
  call mpmul (t4(0:mpnw+6), t5(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!   t7 = t3 / t6
  call mpdiv (t3(0:mpnw+6), t6(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t7(0:mpnw+6), mpnw1)
!mp prec="mpnw"
!   terf = t7
  call mpeq (t7(0:mpnw+6), terf(0:), mpnw)
else
!mp prec="mpnw1"
!   z2 = z ** 2
  call mpnpwr (z(0:), 2, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), z2(0:mpnw+6), mpnw1)
!   t1 = mprealdp (0.e0_mpdk)
  call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   t2 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   t3 = z
  call mpeq (z(0:), t3(0:mpnw+6), mpnw1)
!   t5 = mprealdp (1.e2_mpdk)
  call mprealdp (1.e2_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)

  do k = 0, itrmx
    if (k > 0) then
      d1 = 1.e0_mpdk - 2.e0_mpdk * k
!       t2 = t2 * d1
      call mpmuld (t2(0:mpnw+6), d1, mpt1(0:mpnw+6), mpnw1)
      call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!       t3 = t2 * t3
      call mpmul (t2(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
      call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
    endif

!     t4 = t2 / t3
    call mpdiv (t2(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!     t1 = t1 + t4
    call mpadd (t1(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!     t6 = t4 / t1
    call mpdiv (t4(0:mpnw+6), t1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!mp prec="mpnwm"
!     tc1 = t6 - eps
    call mpsub (t6(0:mpnw+6), eps(0:mpnwm+5), mpt1(0:mpnw+6), mpnwm)
    call mpeq (mpt1(0:mpnw+6), tc1(0:mpnwm+5), mpnwm)
!     tc2 = t5 - t6
    call mpsub (t5(0:mpnw+6), t6(0:mpnw+6), mpt1(0:mpnw+6), mpnwm)
    call mpeq (mpt1(0:mpnw+6), tc2(0:mpnwm+5), mpnwm)
!     ic1 = mpsgn (tc1)
    mpi1 = mpsgn (tc1(0:mpnwm+5))
    ic1 = mpi1
!     ic2 = mpsgn (tc2)
    mpi1 = mpsgn (tc2(0:mpnwm+5))
    ic2 = mpi1
    if (ic1 < 0 .or. ic2 < 0) goto 130
!mp prec="mpnw1"
!     t5 = t6
    call mpeq (t6(0:mpnw+6), t5(0:mpnw+6), mpnw1)
  enddo

write (mpldb, 3)
3 format ('*** ERFR: End loop error 2')
! call mpabrt (551)
call mpabrt (551)

130 continue

!mp prec="mpnw1"
!   t2 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   t3 = sqrt (mppicon)
  call mpsqrt (mppicon(0:), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!   t4 = exp (z2)
  call mpexp (z2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!   t5 = t3 * t4
  call mpmul (t3(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
!   t6 = t1 / t5
  call mpdiv (t1(0:mpnw+6), t5(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!   t7 = t2 - t6
  call mpsub (t2(0:mpnw+6), t6(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t7(0:mpnw+6), mpnw1)
!mp prec="mpnw"
!   terf = t2 - t6
  call mpsub (t2(0:mpnw+6), t6(0:mpnw+6), mpt1(0:mpnw+6), mpnw)
  call mpeq (mpt1(0:mpnw+6), terf(0:), mpnw)
!   ic1 = mpsgn (z)
  mpi1 = mpsgn (z(0:))
  ic1 = mpi1
  if (ic1 < 0) then
!     terf = - terf
    call mpneg (terf(0:), mpt1(0:mpnw+6), mpnw)
    call mpeq (mpt1(0:mpnw+6), terf(0:), mpnw)
  endif
endif

return
end subroutine mperfr

subroutine mperfcr (z, terfc, mpnw)

!   This evaluates the erfc function, using a combination of two series.
!   In particular, the algorithm is (where B = (mpnw + 1) * mpnbt, and
!   dcon is a constant defined below):

!   if (z == 0) then
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
!mp dim1="0:"
! type (mp_real), intent(in):: z
integer (mpiknd), intent(in):: z(0:)
! type (mp_real), intent(out):: terfc
integer (mpiknd), intent(out):: terfc(0:)
integer, parameter:: itrmx = 100000
real (mpdknd), parameter:: dcon = 100.e0_mpdk
real (mpdknd), parameter:: al2 = 0.69314718055994530941723212145817657e0_mpdk
integer ic1, ic2, ic3, ic4, k, nbt, n1
real (mpdknd) d1, d2
!mp dim1="0:mpnwm+5"
! type (mp_real) eps, tc1, tc2
integer (mpiknd) eps(0:mpnwm+5), tc1(0:mpnwm+5), tc2(0:mpnwm+5)
!mp dim1="0:mpnw+6"
! type (mp_real) t1, t2, t3, t4, t5
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), t5(0:mpnw+6)
! type (mp_real) t6, t7, z2
integer (mpiknd) t6(0:mpnw+6), t7(0:mpnw+6), z2(0:mpnw+6)
integer mpnw1

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpt3(0:mpnw+6)
integer (mpiknd) mpt4(0:mpnw+6), mpt5(0:mpnw+6), mpt6(0:mpnw+6)

call mpinitwds (mpt1, mpnw+6-5)
call mpinitwds (mpt2, mpnw+6-5)
call mpinitwds (mpt3, mpnw+6-5)
call mpinitwds (mpt4, mpnw+6-5)
call mpinitwds (mpt5, mpnw+6-5)
call mpinitwds (mpt6, mpnw+6-5)
call mpinitwds (eps(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc1(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc2(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (t1(0:mpnw+6), mpnw+6-5)
call mpinitwds (t2(0:mpnw+6), mpnw+6-5)
call mpinitwds (t3(0:mpnw+6), mpnw+6-5)
call mpinitwds (t4(0:mpnw+6), mpnw+6-5)
call mpinitwds (t5(0:mpnw+6), mpnw+6-5)
call mpinitwds (t6(0:mpnw+6), mpnw+6-5)
call mpinitwds (t7(0:mpnw+6), mpnw+6-5)
call mpinitwds (z2(0:mpnw+6), mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)

! ic1 = mpspacer (terfc)
ic1 = mpspacer (terfc(0:))
if (mpnw < mpnwm .or. ic1 < mpnw + 6) then
  write (6, 1)
1 format ('*** ERFCR: Uninitialized or inadequately sized array')
!   call mpabrt (552)
  call mpabrt (552)
endif

!mp prec="mpnwm"

! tc2 = mprealdp (2.e0_mpdk)
call mprealdp (2.e0_mpdk, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), tc2(0:mpnwm+5), mpnwm)
ic2 = -mpnw1*mpnbt
! eps = tc2 ** ic2
call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), eps(0:mpnwm+5), mpnwm)

!mp prec="mpnw1"

nbt = mpnbt
d1 = aint (1.e0_mpdk + sqrt (nbt * al2))
d2 = aint (nbt / dcon + 8.e0_mpdk)
! t1 = mprealdp (d1)
call mprealdp (d1, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
! t2 = mprealdp (d2)
call mprealdp (d2, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! t3 = z - t1
call mpsub (z(0:), t1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
! t4 = z + t1
call mpadd (z(0:), t1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
! t5 = z - t2
call mpsub (z(0:), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
! ic1 = mpsgn (z)
mpi1 = mpsgn (z(0:))
ic1 = mpi1
! ic2 = mpsgn (t3)
mpi1 = mpsgn (t3(0:mpnw+6))
ic2 = mpi1
! ic3 = mpsgn (t4)
mpi1 = mpsgn (t4(0:mpnw+6))
ic3 = mpi1
! ic4 = mpsgn (t5)
mpi1 = mpsgn (t5(0:mpnw+6))
ic4 = mpi1

if (ic1 == 0) then
!   terfc = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), terfc(0:), mpnw1)
elseif (ic2 > 0) then
!   terfc = mprealdp (0.e0_mpdk)
  call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), terfc(0:), mpnw1)
elseif (ic3 < 0) then
!   terfc = mprealdp (2.e0_mpdk)
  call mprealdp (2.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), terfc(0:), mpnw1)
elseif (ic4 < 0) then
!   z2 = z ** 2
  call mpnpwr (z(0:), 2, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), z2(0:mpnw+6), mpnw1)
!   t1 = mprealdp (0.e0_mpdk)
  call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   t2 = z
  call mpeq (z(0:), t2(0:mpnw+6), mpnw1)
!   t3 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!   t5 = mprealdp (1.e2_mpdk)
  call mprealdp (1.e2_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)

  do k = 0, itrmx
    if (k > 0) then
!       t6 = z2 * 2.e0_mpdk
      call mpmuld (z2(0:mpnw+6), 2.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
      call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!       t2 = t6 * t2
      call mpmul (t6(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
      call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
      d1 = real (2*k+1, mpdk)
!       t3 = t3 * d1
      call mpmuld (t3(0:mpnw+6), d1, mpt1(0:mpnw+6), mpnw1)
      call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
    endif

!     t4 = t2 / t3
    call mpdiv (t2(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!     t1 = t1 + t4
    call mpadd (t1(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!     t6 = t4 / t1
    call mpdiv (t4(0:mpnw+6), t1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!mp prec="mpnwm"
!     tc1 = t6 - eps
    call mpsub (t6(0:mpnw+6), eps(0:mpnwm+5), mpt1(0:mpnw+6), mpnwm)
    call mpeq (mpt1(0:mpnw+6), tc1(0:mpnwm+5), mpnwm)
!     tc2 = t5 - t6
    call mpsub (t5(0:mpnw+6), t6(0:mpnw+6), mpt1(0:mpnw+6), mpnwm)
    call mpeq (mpt1(0:mpnw+6), tc2(0:mpnwm+5), mpnwm)
!     ic1 = mpsgn (tc1)
    mpi1 = mpsgn (tc1(0:mpnwm+5))
    ic1 = mpi1
!     ic2 = mpsgn (tc2)
    mpi1 = mpsgn (tc2(0:mpnwm+5))
    ic2 = mpi1
    if (ic1 < 0 .or. ic2 < 0) goto 120
!mp prec="mpnw1"
!     t5 = t6
    call mpeq (t6(0:mpnw+6), t5(0:mpnw+6), mpnw1)
  enddo

write (mpldb, 2)
2 format ('*** ERFRC: End loop error 1')
! call mpabrt (553)
call mpabrt (553)

120 continue

!mp prec="mpnw1"
!   t2 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   t3 = t1 * 2.e0_mpdk
  call mpmuld (t1(0:mpnw+6), 2.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!   t4 = sqrt (mppicon)
  call mpsqrt (mppicon(0:), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!   t5 = exp (z2)
  call mpexp (z2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
!   t6 = t4 * t5
  call mpmul (t4(0:mpnw+6), t5(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!   t7 = t2 - t3 / t6
  call mpdiv (t3(0:mpnw+6), t6(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpsub (t2(0:mpnw+6), mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
  call mpeq (mpt2(0:mpnw+6), t7(0:mpnw+6), mpnw1)
!mp prec="mpnw"
!   terfc = t7
  call mpeq (t7(0:mpnw+6), terfc(0:), mpnw)
else
!mp prec="mpnw1"
!   z2 = z ** 2
  call mpnpwr (z(0:), 2, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), z2(0:mpnw+6), mpnw1)
!   t1 = mprealdp (0.e0_mpdk)
  call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   t2 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   t3 = z
  call mpeq (z(0:), t3(0:mpnw+6), mpnw1)
!   t5 = mprealdp (1.e2_mpdk)
  call mprealdp (1.e2_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)

  do k = 0, itrmx
    if (k > 0) then
      d1 = 1.e0_mpdk - 2.e0_mpdk * k
!       t2 = t2 * d1
      call mpmuld (t2(0:mpnw+6), d1, mpt1(0:mpnw+6), mpnw1)
      call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!       t3 = t2 * t3
      call mpmul (t2(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
      call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
    endif

!     t4 = t2 / t3
    call mpdiv (t2(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!     t1 = t1 + t4
    call mpadd (t1(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!     t6 = t4 / t1
    call mpdiv (t4(0:mpnw+6), t1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!mp prec="mpnwm"
!     tc1 = t6 - eps
    call mpsub (t6(0:mpnw+6), eps(0:mpnwm+5), mpt1(0:mpnw+6), mpnwm)
    call mpeq (mpt1(0:mpnw+6), tc1(0:mpnwm+5), mpnwm)
!     tc2 = t5 - t6
    call mpsub (t5(0:mpnw+6), t6(0:mpnw+6), mpt1(0:mpnw+6), mpnwm)
    call mpeq (mpt1(0:mpnw+6), tc2(0:mpnwm+5), mpnwm)
!     ic1 = mpsgn (tc1)
    mpi1 = mpsgn (tc1(0:mpnwm+5))
    ic1 = mpi1
!     ic2 = mpsgn (tc2)
    mpi1 = mpsgn (tc2(0:mpnwm+5))
    ic2 = mpi1
    if (ic1 < 0 .or. ic2 < 0) goto 130
!mp prec="mpnw1"
!     t5 = t6
    call mpeq (t6(0:mpnw+6), t5(0:mpnw+6), mpnw1)
  enddo

write (mpldb, 3)
3 format ('*** ERFRC: End loop error 1')
! call mpabrt (554)
call mpabrt (554)

130 continue

!mp prec="mpnw1"
!   t3 = sqrt (mppicon)
  call mpsqrt (mppicon(0:), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!   t4 = exp (z2)
  call mpexp (z2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!   t5 = t3 * t4
  call mpmul (t3(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
!   t6 = t1 / t5
  call mpdiv (t1(0:mpnw+6), t5(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!   ic1 = mpsgn (z)
  mpi1 = mpsgn (z(0:))
  ic1 = mpi1
  if (ic1 < 0) then
!     t6 = 2.e0_mpdk - t6
    call mpdmc (2.e0_mpdk, 0, mpt2(0:mpnw+6), mpnw1)
    call mpsub (mpt2(0:mpnw+6), t6(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
  endif
!mp prec="mpnw"
!   terfc = t6
  call mpeq (t6(0:mpnw+6), terfc(0:), mpnw)
endif

return
end subroutine mperfcr

subroutine mpexpint (x, y, mpnw)

!   This evaluates the exponential integral function Ei(x):
!   Ei(x) = - incgamma (0, -x)

implicit none
integer, intent(in):: mpnw
!mp dim1="0:"
! type (mp_real), intent(in):: x
integer (mpiknd), intent(in):: x(0:)
! type (mp_real), intent(out):: y
integer (mpiknd), intent(out):: y(0:)
integer ic1
!mp dim1="0:mpnw+6"
! type (mp_real) t1, t2, t3
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6)
integer mpnw1

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpt3(0:mpnw+6)
integer (mpiknd) mpt4(0:mpnw+6), mpt5(0:mpnw+6), mpt6(0:mpnw+6)

call mpinitwds (mpt1, mpnw+6-5)
call mpinitwds (mpt2, mpnw+6-5)
call mpinitwds (mpt3, mpnw+6-5)
call mpinitwds (mpt4, mpnw+6-5)
call mpinitwds (mpt5, mpnw+6-5)
call mpinitwds (mpt6, mpnw+6-5)
call mpinitwds (t1(0:mpnw+6), mpnw+6-5)
call mpinitwds (t2(0:mpnw+6), mpnw+6-5)
call mpinitwds (t3(0:mpnw+6), mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)

! ic1 = mpspacer (y)
ic1 = mpspacer (y(0:))
if (mpnw < mpnwm .or. ic1 < mpnw + 6) then
  write (6, 1)
1 format ('*** EXPINT: Uninitialized or inadequately sized array')
!   call mpabrt (555)
  call mpabrt (555)
endif

!mp prec="mpnw1"

! ic1 = mpsgn (x)
mpi1 = mpsgn (x(0:))
ic1 = mpi1
if (ic1 == 0) then
  write (mpldb, 2)
2 format ('*** EXPINT: Argument is zero')
!   call mpabrt (556)
  call mpabrt (556)
endif

! t1 = mprealdp (0.e0_mpdk)
call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
! t2 = - x
call mpneg (x(0:), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! call mpincgammar (t1, t2, t3, mpnw1)
call mpincgammar (t1, t2, t3, mpnw1)
!mp prec="mpnw"
! y = - t3
call mpneg (t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw)
call mpeq (mpt1(0:mpnw+6), y(0:), mpnw)
return
end subroutine mpexpint

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
!mp dim1="0:"
! type (mp_real), intent(in):: t
integer (mpiknd), intent(in):: t(0:)
! type (mp_real), intent(out):: z
integer (mpiknd), intent(out):: z(0:)
integer, parameter:: itrmx = 100000
real (mpdknd), parameter:: al2 = 0.69314718055994530941723212145817657e0_mpdk
real (mpdknd), parameter:: dmax = 1e8_mpdk
integer i, ic1, ic2, i1, i2, j, nt, n1, n2, n3
real (mpdknd) alpha, d1, d2, d3
!mp dim1="0:mpnwm+5"
! type (mp_real) eps, tc1, tc2
integer (mpiknd) eps(0:mpnwm+5), tc1(0:mpnwm+5), tc2(0:mpnwm+5)
!mp dim1="0:mpnw+6"
! type (mp_real) f1, sum1, sum2, tn
integer (mpiknd) f1(0:mpnw+6), sum1(0:mpnw+6), sum2(0:mpnw+6), tn(0:mpnw+6)
! type (mp_real) t1, t2, t3, t4
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6)
! type (mp_real) t5, t6
integer (mpiknd) t5(0:mpnw+6), t6(0:mpnw+6)
integer mpnw1

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpt3(0:mpnw+6)
integer (mpiknd) mpt4(0:mpnw+6), mpt5(0:mpnw+6), mpt6(0:mpnw+6)

call mpinitwds (mpt1, mpnw+6-5)
call mpinitwds (mpt2, mpnw+6-5)
call mpinitwds (mpt3, mpnw+6-5)
call mpinitwds (mpt4, mpnw+6-5)
call mpinitwds (mpt5, mpnw+6-5)
call mpinitwds (mpt6, mpnw+6-5)
call mpinitwds (eps(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc1(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc2(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (f1(0:mpnw+6), mpnw+6-5)
call mpinitwds (sum1(0:mpnw+6), mpnw+6-5)
call mpinitwds (sum2(0:mpnw+6), mpnw+6-5)
call mpinitwds (tn(0:mpnw+6), mpnw+6-5)
call mpinitwds (t1(0:mpnw+6), mpnw+6-5)
call mpinitwds (t2(0:mpnw+6), mpnw+6-5)
call mpinitwds (t3(0:mpnw+6), mpnw+6-5)
call mpinitwds (t4(0:mpnw+6), mpnw+6-5)
call mpinitwds (t5(0:mpnw+6), mpnw+6-5)
call mpinitwds (t6(0:mpnw+6), mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)

! ic1 = mpspacer (z)
ic1 = mpspacer (z(0:))
if (mpnw < mpnwm .or. ic1 < mpnw + 6) then
  write (6, 1)
1 format ('*** GAMMAR: Uninitialized or inadequately sized array')
!   call mpabrt (557)
  call mpabrt (557)
endif

!mp prec="mpnwm"

! tc2 = mprealdp (2.e0_mpdk)
call mprealdp (2.e0_mpdk, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), tc2(0:mpnwm+5), mpnwm)
ic2 = -mpnw1*mpnbt
! eps = tc2 ** ic2
call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), eps(0:mpnwm+5), mpnwm)

!mp prec="mpnw1"

! t1 = t - anint (t)
call mpnint (t(0:), mpt1(0:mpnw+6), mpnw1)
call mpsub (t(0:), mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
call mpeq (mpt2(0:mpnw+6), t1(0:mpnw+6), mpnw1)
! i1 = mpsgn (t)
mpi1 = mpsgn (t(0:))
i1 = mpi1
! i2 = mpsgn (t1)
mpi1 = mpsgn (t1(0:mpnw+6))
i2 = mpi1
if (i1 == 0 .or. d1 > dmax .or. (i1 < 0 .and. i2 == 0)) then
  write (mpldb, 2) dmax
2 format ('*** GAMMAR: Input argument must have absolute value <=',f10.0,','/ &
  'must not be zero, and if negative must not be an integer.')
!   call mpabrt (558)
  call mpabrt (558)
endif

! f1 = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), f1(0:mpnw+6), mpnw1)

!   Find the integer and fractional parts of t.

! t2 = aint (t)
call mpinfr (t(0:), mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! t3 = t - t2
call mpsub (t(0:), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)

! ic1 = mpsgn (t3)
mpi1 = mpsgn (t3(0:mpnw+6))
ic1 = mpi1
! ic2 = mpsgn (t)
mpi1 = mpsgn (t(0:))
ic2 = mpi1

if (ic1 == 0) then

!   If t is a positive integer, then apply the usual factorial recursion.

!   d1 = dpreal (t2)
  call mpmdc (t2(0:mpnw+6), mpd1, mpi1, mpnw1)
  mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
  d1 = mpd1
  nt = int (d1)
!   t1 = f1
  call mpeq (f1(0:mpnw+6), t1(0:mpnw+6), mpnw1)

  do i = 2, nt - 1
!     t1 = t1 * real (i, mpdk)
    call mpmuld (t1(0:mpnw+6), real (i, mpdk), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
  enddo

!mp prec="mpnw"
!   z = t1
  call mpeq (t1(0:mpnw+6), z(0:), mpnw)
  goto 120
elseif (ic2 > 0) then

!   Apply the identity Gamma[t+1] = t * Gamma[t] to reduce the input argument
!   to the unit interval.

!   d1 = dpreal (t2)
  call mpmdc (t2(0:mpnw+6), mpd1, mpi1, mpnw)
  mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
  d1 = mpd1
  nt = nint (d1)
!   t1 = f1
  call mpeq (f1(0:mpnw+6), t1(0:mpnw+6), mpnw)
!   tn = t3
  call mpeq (t3(0:mpnw+6), tn(0:mpnw+6), mpnw)

  do i = 1, nt
    d2 = real (i, mpdk)
!     t4 = mprealdp (d2)
    call mprealdp (d2, mpt1(0:mpnw+6), mpnw)
    call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw)
!     t5 = t - t4
    call mpsub (t(0:), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw)
    call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw)
!     t1 = t1 * t5
    call mpmul (t1(0:mpnw+6), t5(0:mpnw+6), mpt1(0:mpnw+6), mpnw)
    call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw)
  enddo
else

!   Apply the gamma identity to reduce a negative argument to the unit interval.

!   t4 = f1 - t
  call mpsub (f1(0:mpnw+6), t(0:), mpt1(0:mpnw+6), mpnw)
  call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw)
!   t3 = aint (t4)
  call mpinfr (t4(0:mpnw+6), mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpnw)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw)
!   t5 = t4 - t3
  call mpsub (t4(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw)
  call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw)
!   d1 = dpreal (t3)
  call mpmdc (t3(0:mpnw+6), mpd1, mpi1, mpnw)
  mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
  d1 = mpd1
  nt = nint (d1)
!   t1 = f1
  call mpeq (f1(0:mpnw+6), t1(0:mpnw+6), mpnw)
!   t2 = f1 - t5
  call mpsub (f1(0:mpnw+6), t5(0:mpnw+6), mpt1(0:mpnw+6), mpnw)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw)
!   tn = t2
  call mpeq (t2(0:mpnw+6), tn(0:mpnw+6), mpnw)

  do i = 0, nt - 1
    d2 = real (i, mpdk)
!     t4 = mprealdp (d2)
    call mprealdp (d2, mpt1(0:mpnw+6), mpnw)
    call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw)
!     t5 = t + t4
    call mpadd (t(0:), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw)
    call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw)
!     t1 = t1 / t5
    call mpdiv (t1(0:mpnw+6), t5(0:mpnw+6), mpt1(0:mpnw+6), mpnw)
    call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw)
  enddo
endif

!   Calculate alpha = bits of precision * log(2) / 2, then take the next even
!   integer value, so that alpha/2 and alpha^2/4 can be calculated exactly in DP.

alpha = 2.e0_mpdk * aint (0.25e0_mpdk * mpnw1 * mpnbt * al2 + 1.e0_mpdk)
d2 = 0.25e0_mpdk * alpha ** 2
! t2 = tn
call mpeq (tn(0:mpnw+6), t2(0:mpnw+6), mpnw)
! t3 = f1 / t2
call mpdiv (f1(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw)
! sum1 = t3
call mpeq (t3(0:mpnw+6), sum1(0:mpnw+6), mpnw)

!   Evaluate the series with t.

do j = 1, itrmx
  d3 = real (j, mpdk)
!   t6 = mprealdp (d3)
  call mprealdp (d3, mpt1(0:mpnw+6), mpnw)
  call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw)
!   t4 = t2 + t6
  call mpadd (t2(0:mpnw+6), t6(0:mpnw+6), mpt1(0:mpnw+6), mpnw)
  call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw)
!   t5 = t4 * real (j, mpdk)
  call mpmuld (t4(0:mpnw+6), real (j, mpdk), mpt1(0:mpnw+6), mpnw)
  call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw)
!   t6 = t3 / t5
  call mpdiv (t3(0:mpnw+6), t5(0:mpnw+6), mpt1(0:mpnw+6), mpnw)
  call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw)
!   t3 = t6 * d2
  call mpmuld (t6(0:mpnw+6), d2, mpt1(0:mpnw+6), mpnw)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw)
!   sum1 = sum1 + t3
  call mpadd (sum1(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw)
  call mpeq (mpt1(0:mpnw+6), sum1(0:mpnw+6), mpnw)
!mp prec="mpnwm"
!   tc1 = abs (t3) - eps * abs (sum1)
  call mpabs (t3(0:mpnw+6), mpt1(0:mpnw+6), mpnwm)
  call mpabs (sum1(0:mpnw+6), mpt2(0:mpnw+6), mpnwm)
  call mpmul (eps(0:mpnwm+5), mpt2(0:mpnw+6), mpt3(0:mpnw+6), mpnwm)
  call mpsub (mpt1(0:mpnw+6), mpt3(0:mpnw+6), mpt4(0:mpnw+6), mpnwm)
  call mpeq (mpt4(0:mpnw+6), tc1(0:mpnwm+5), mpnwm)
!   ic1 = mpsgn (tc1)
  mpi1 = mpsgn (tc1(0:mpnwm+5))
  ic1 = mpi1
  if (ic1 < 0) goto 100
enddo

write (mpldb, 3)
3 format ('*** GAMMAR: End loop error 1')
! call mpabrt (559)
call mpabrt (559)

100 continue

!mp prec="mpnw1"
! t2 = - tn
call mpneg (tn(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! t3 = f1 / t2
call mpdiv (f1(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
! sum2 = t3
call mpeq (t3(0:mpnw+6), sum2(0:mpnw+6), mpnw1)

!   Evaluate the same series with -t.

do j = 1, itrmx
  d3 = real (j, mpdk)
!   t6 = mprealdp (d3)
  call mprealdp (d3, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!   t4 = t2 + t6
  call mpadd (t2(0:mpnw+6), t6(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!   t5 = t4 * real (j, mpdk)
  call mpmuld (t4(0:mpnw+6), real (j, mpdk), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
!   t6 = t3 / t5
  call mpdiv (t3(0:mpnw+6), t5(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!   t3 = t6 * d2
  call mpmuld (t6(0:mpnw+6), d2, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!   sum2 = sum2 + t3
  call mpadd (sum2(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), sum2(0:mpnw+6), mpnw1)
!mp prec="mpnwm"
!   tc1 = abs (t3) - eps * abs (sum2)
  call mpabs (t3(0:mpnw+6), mpt1(0:mpnw+6), mpnwm)
  call mpabs (sum2(0:mpnw+6), mpt2(0:mpnw+6), mpnwm)
  call mpmul (eps(0:mpnwm+5), mpt2(0:mpnw+6), mpt3(0:mpnw+6), mpnwm)
  call mpsub (mpt1(0:mpnw+6), mpt3(0:mpnw+6), mpt4(0:mpnw+6), mpnwm)
  call mpeq (mpt4(0:mpnw+6), tc1(0:mpnwm+5), mpnwm)
!   ic1 = mpsgn (tc1)
  mpi1 = mpsgn (tc1(0:mpnwm+5))
  ic1 = mpi1
  if (ic1 < 0) goto 110
enddo

write (mpldb, 4)
4 format ('*** GAMMAR: End loop error 2')
! call mpabrt (560)
call mpabrt (560)


110 continue

!   Compute sqrt (pi * sum1 / (tn * sin (pi * tn) * sum2))
!   and (alpha/2)^tn terms. Also, multiply by the factor t1, from the
!   If block above.

!mp prec="mpnw1"
! t2 = mppicon
call mpeq (mppicon(0:), t2(0:mpnw+6), mpnw1)
! t3 = t2 * tn
call mpmul (t2(0:mpnw+6), tn(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
! t4 = cos (t3)
call mpcssnr (t3(0:mpnw+6), mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
! t5 = sin (t3)
call mpcssnr (t3(0:mpnw+6), mpt2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
! t6 = t5 * sum2
call mpmul (t5(0:mpnw+6), sum2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
! t5 = tn * t6
call mpmul (tn(0:mpnw+6), t6(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
! t3 = t2 * sum1
call mpmul (t2(0:mpnw+6), sum1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
! t6 = - t3 / t5
call mpdiv (t3(0:mpnw+6), t5(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpneg (mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
call mpeq (mpt2(0:mpnw+6), t6(0:mpnw+6), mpnw1)
! t2 = sqrt (t6)
call mpsqrt (t6(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
d3 = 0.5e0_mpdk * alpha
! t3 = mprealdp (d3)
call mprealdp (d3, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
! t4 = log (t3)
call mplog (t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
! t5 = tn * t4
call mpmul (tn(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
! t6 = exp (t5)
call mpexp (t5(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
! t3 = t2 * t6
call mpmul (t2(0:mpnw+6), t6(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
! t4 = t1 * t3
call mpmul (t1(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!mp prec="mpnw"
! z = t4
call mpeq (t4(0:mpnw+6), z(0:), mpnw)

120 continue

return
end subroutine mpgammar

subroutine mphurwitzzetan (is, aa, zz, mpnw)

!   This returns the Hurwitz zeta function of IS and AA, using an algorithm from:
!   David H. Bailey and Jonathan M. Borwein, "Crandall's computation of the
!   incomplete gamma function and the Hurwitz zeta function with applications to
!   Dirichlet L-series," Applied Mathematics and Computation, vol. 268C (Oct 2015),
!   pg. 462-477, preprint at:
!   https://www.davidhbailey.com/dhbpapers/lerch.pdf
!   This is limited to IS >= 2 and 0 < AA < 1.

implicit none
integer, intent(in):: mpnw
integer, intent(in):: is
!mp dim1="0:"
! type (mp_real), intent(in):: aa
integer (mpiknd), intent(in):: aa(0:)
! type (mp_real), intent(out):: zz
integer (mpiknd), intent(out):: zz(0:)
integer, parameter:: itrmax = 1000000
integer i1, ic1, ic2, ic3, k, n1
real (mpdknd) d1, d2, dk
!mp dim1="0:mpnwm+5"
! type (mp_real) eps, tc1, tc2, tc3
integer (mpiknd) eps(0:mpnwm+5), tc1(0:mpnwm+5), tc2(0:mpnwm+5), tc3(0:mpnwm+5)
!mp dim1="0:mpnw+6"
! type (mp_real) gs1, gs2, ss, sum1, sum2
integer (mpiknd) gs1(0:mpnw+6), gs2(0:mpnw+6), ss(0:mpnw+6), sum1(0:mpnw+6), sum2(0:mpnw+6)
! type (mp_real) sum3, ss1, ss2, ss3, ss4
integer (mpiknd) sum3(0:mpnw+6), ss1(0:mpnw+6), ss2(0:mpnw+6), ss3(0:mpnw+6), ss4(0:mpnw+6)
! type (mp_real) s1, s2, s3, t1, t2
integer (mpiknd) s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6)
! type (mp_real) t3, t4, t5, t6
integer (mpiknd) t3(0:mpnw+6), t4(0:mpnw+6), t5(0:mpnw+6), t6(0:mpnw+6)
! type (mp_real) t7, t8
integer (mpiknd) t7(0:mpnw+6), t8(0:mpnw+6)
integer mpnw1

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpt3(0:mpnw+6)
integer (mpiknd) mpt4(0:mpnw+6), mpt5(0:mpnw+6), mpt6(0:mpnw+6)

call mpinitwds (mpt1, mpnw+6-5)
call mpinitwds (mpt2, mpnw+6-5)
call mpinitwds (mpt3, mpnw+6-5)
call mpinitwds (mpt4, mpnw+6-5)
call mpinitwds (mpt5, mpnw+6-5)
call mpinitwds (mpt6, mpnw+6-5)
call mpinitwds (eps(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc1(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc2(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc3(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (gs1(0:mpnw+6), mpnw+6-5)
call mpinitwds (gs2(0:mpnw+6), mpnw+6-5)
call mpinitwds (ss(0:mpnw+6), mpnw+6-5)
call mpinitwds (sum1(0:mpnw+6), mpnw+6-5)
call mpinitwds (sum2(0:mpnw+6), mpnw+6-5)
call mpinitwds (sum3(0:mpnw+6), mpnw+6-5)
call mpinitwds (ss1(0:mpnw+6), mpnw+6-5)
call mpinitwds (ss2(0:mpnw+6), mpnw+6-5)
call mpinitwds (ss3(0:mpnw+6), mpnw+6-5)
call mpinitwds (ss4(0:mpnw+6), mpnw+6-5)
call mpinitwds (s1(0:mpnw+6), mpnw+6-5)
call mpinitwds (s2(0:mpnw+6), mpnw+6-5)
call mpinitwds (s3(0:mpnw+6), mpnw+6-5)
call mpinitwds (t1(0:mpnw+6), mpnw+6-5)
call mpinitwds (t2(0:mpnw+6), mpnw+6-5)
call mpinitwds (t3(0:mpnw+6), mpnw+6-5)
call mpinitwds (t4(0:mpnw+6), mpnw+6-5)
call mpinitwds (t5(0:mpnw+6), mpnw+6-5)
call mpinitwds (t6(0:mpnw+6), mpnw+6-5)
call mpinitwds (t7(0:mpnw+6), mpnw+6-5)
call mpinitwds (t8(0:mpnw+6), mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)

! ic1 = mpspacer (zz)
ic1 = mpspacer (zz(0:))
if (mpnw < mpnwm .or. ic1 < mpnw + 6) then
  write (6, 1)
1 format ('*** HURWITZZETAN: Uninitialized or inadequately sized array')
!   call mpabrt (561)
  call mpabrt (561)
endif


! ic1 = mpwprecr (mppicon)
ic1 = mpwprecr (mppicon(0:))
if (ic1 < mpnw1) then
  write (6, 2) mpnw1
2 format ('*** HURWITZZETAN: Pi must be initialized to',i6,' words by calling mpinit.')
!   call mpabrt (562)
  call mpabrt (562)
endif

!mp prec="mpnwm"

! tc2 = mprealdp (2.e0_mpdk)
call mprealdp (2.e0_mpdk, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), tc2(0:mpnwm+5), mpnwm)
ic2 = -mpnw1*mpnbt
! eps = tc2 ** ic2
call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), eps(0:mpnwm+5), mpnwm)

!mp prec="mpnw1"

if (is <= 0) then
  write (mpldb, 3)
3 format ('*** HURWITZZETAN: Integer argument less than or equal to 0:')
!   call mpabrt (563)
  call mpabrt (563)
endif

! t1 = mprealdp (0.e0_mpdk)
call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
! t2 = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! t3 = aa - t1
call mpsub (aa(0:), t1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
! ic1 = mpsgn (t3)
mpi1 = mpsgn (t3(0:mpnw+6))
ic1 = mpi1
! t4 = aa - t2
call mpsub (aa(0:), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
! ic2 = mpsgn (t4)
mpi1 = mpsgn (t4(0:mpnw+6))
ic2 = mpi1
if (ic1 <= 0 .or. ic2 >= 0) then
  write (mpldb, 4)
4 format ('*** HURWITZZETAN: Second argument must be in the range (0, 1)')
!   call mpabrt (564)
  call mpabrt (564)
endif

d2 = real (is, mpdk)
! ss = mprealdp (d2)
call mprealdp (d2, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), ss(0:mpnw+6), mpnw1)
! ss1 = ss * 0.5e0_mpdk
call mpmuld (ss(0:mpnw+6), 0.5e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), ss1(0:mpnw+6), mpnw1)
! t1 = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
! t2 = t1 + ss
call mpadd (t1(0:mpnw+6), ss(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! ss2 = t2 * 0.5e0_mpdk
call mpmuld (t2(0:mpnw+6), 0.5e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), ss2(0:mpnw+6), mpnw1)
! t2 = t1 - ss
call mpsub (t1(0:mpnw+6), ss(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! ss3 = t2 * 0.5e0_mpdk
call mpmuld (t2(0:mpnw+6), 0.5e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), ss3(0:mpnw+6), mpnw1)
! ss4 = t1 - ss1
call mpsub (t1(0:mpnw+6), ss1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), ss4(0:mpnw+6), mpnw1)
! call mpgammar (ss1, gs1, mpnw1)
call mpgammar (ss1, gs1, mpnw1)
! call mpgammar (ss2, gs2, mpnw1)
call mpgammar (ss2, gs2, mpnw1)
! t2 = aa ** 2
call mpnpwr (aa(0:), 2, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! t1 = mppicon * t2
call mpmul (mppicon(0:), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)

! call mpincgammar (ss1, t1, t2, mpnw1)
call mpincgammar (ss1, t1, t2, mpnw1)
! t3 = t2 / gs1
call mpdiv (t2(0:mpnw+6), gs1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
! call mpincgammar (ss2, t1, t2, mpnw1)
call mpincgammar (ss2, t1, t2, mpnw1)
! t4 = t2 / gs2
call mpdiv (t2(0:mpnw+6), gs2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
! t2 = t3 + t4
call mpadd (t3(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! t3 = abs (aa)
call mpabs (aa(0:), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
! t4 = t3 ** is
call mpnpwr (t3(0:mpnw+6), is, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
! sum1 = t2 / t4
call mpdiv (t2(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), sum1(0:mpnw+6), mpnw1)
! sum2 = mprealdp (0.e0_mpdk)
call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), sum2(0:mpnw+6), mpnw1)
! sum3 = mprealdp (0.e0_mpdk)
call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), sum3(0:mpnw+6), mpnw1)

do k = 1, itrmax
  dk = real (k, mpdk)
!   t5 = mprealdp (dk)
  call mprealdp (dk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
!   t6 = t5 + aa
  call mpadd (t5(0:mpnw+6), aa(0:), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!   t5 = t6 ** 2
  call mpnpwr (t6(0:mpnw+6), 2, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
!   t1 = mppicon * t5
  call mpmul (mppicon(0:), t5(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   t5 = - mprealdp (dk)
  call mprealdp (dk, mpt1(0:mpnw+6), mpnw1)
  call mpneg (mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
  call mpeq (mpt2(0:mpnw+6), t5(0:mpnw+6), mpnw1)
!   t6 = t5 + aa
  call mpadd (t5(0:mpnw+6), aa(0:), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!   t7 = t6 ** 2
  call mpnpwr (t6(0:mpnw+6), 2, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t7(0:mpnw+6), mpnw1)
!   t2 = mppicon * t7
  call mpmul (mppicon(0:), t7(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   t6 = t5 ** 2
  call mpnpwr (t5(0:mpnw+6), 2, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!   t3 = mppicon * t6
  call mpmul (mppicon(0:), t6(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!   t5 = mppicon * 2.e0_mpdk * dk
  call mpmuld (mppicon(0:), 2.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpmuld (mpt1(0:mpnw+6), dk, mpt2(0:mpnw+6), mpnw1)
  call mpeq (mpt2(0:mpnw+6), t5(0:mpnw+6), mpnw1)
!   t4 = t5 * aa
  call mpmul (t5(0:mpnw+6), aa(0:), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)

!   call mpincgammar (ss1, t1, t5, mpnw1)
  call mpincgammar (ss1, t1, t5, mpnw1)
!   t6 = t5 / gs1
  call mpdiv (t5(0:mpnw+6), gs1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!   call mpincgammar (ss2, t1, t5, mpnw1)
  call mpincgammar (ss2, t1, t5, mpnw1)
!   t7 = t5 / gs2
  call mpdiv (t5(0:mpnw+6), gs2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t7(0:mpnw+6), mpnw1)
!   t5 = t6 + t7
  call mpadd (t6(0:mpnw+6), t7(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
!   t6 = mprealdp (dk)
  call mprealdp (dk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!   t7 = t6 + aa
  call mpadd (t6(0:mpnw+6), aa(0:), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t7(0:mpnw+6), mpnw1)
!   t6 = abs (t7)
  call mpabs (t7(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!   t7 = t6 ** is
  call mpnpwr (t6(0:mpnw+6), is, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t7(0:mpnw+6), mpnw1)
!   s1 = t5 / t7
  call mpdiv (t5(0:mpnw+6), t7(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), s1(0:mpnw+6), mpnw1)

!   call mpincgammar (ss1, t2, t5, mpnw1)
  call mpincgammar (ss1, t2, t5, mpnw1)
!   t6 = t5 / gs1
  call mpdiv (t5(0:mpnw+6), gs1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!   call mpincgammar (ss2, t2, t5, mpnw1)
  call mpincgammar (ss2, t2, t5, mpnw1)
!   t7 = t5 / gs2
  call mpdiv (t5(0:mpnw+6), gs2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t7(0:mpnw+6), mpnw1)
!   t5 = t6 - t7
  call mpsub (t6(0:mpnw+6), t7(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
!   t6 = - mprealdp (dk)
  call mprealdp (dk, mpt1(0:mpnw+6), mpnw1)
  call mpneg (mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
  call mpeq (mpt2(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!   t7 = t6 + aa
  call mpadd (t6(0:mpnw+6), aa(0:), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t7(0:mpnw+6), mpnw1)
!   t6 = abs (t7)
  call mpabs (t7(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!   t7 = t6 ** is
  call mpnpwr (t6(0:mpnw+6), is, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t7(0:mpnw+6), mpnw1)
!   s2 = t5 / t7
  call mpdiv (t5(0:mpnw+6), t7(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), s2(0:mpnw+6), mpnw1)
!   sum1 = sum1 + s1
  call mpadd (sum1(0:mpnw+6), s1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), sum1(0:mpnw+6), mpnw1)
!   sum2 = sum2 + s2
  call mpadd (sum2(0:mpnw+6), s2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), sum2(0:mpnw+6), mpnw1)

!   call mpincgammar (ss3, t3, t5, mpnw1)
  call mpincgammar (ss3, t3, t5, mpnw1)
!   t6 = cos (t4)
  call mpcssnr (t4(0:mpnw+6), mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!   t7 = sin (t4)
  call mpcssnr (t4(0:mpnw+6), mpt2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t7(0:mpnw+6), mpnw1)
!   t8 = t5 * t6
  call mpmul (t5(0:mpnw+6), t6(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t8(0:mpnw+6), mpnw1)
!   t5 = t8 / gs1
  call mpdiv (t8(0:mpnw+6), gs1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
!   call mpincgammar (ss4, t3, t6, mpnw1)
  call mpincgammar (ss4, t3, t6, mpnw1)
!   t8 = t6 * t7
  call mpmul (t6(0:mpnw+6), t7(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t8(0:mpnw+6), mpnw1)
!   t6 = t8 / gs2
  call mpdiv (t8(0:mpnw+6), gs2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!   t7 = t5 + t6
  call mpadd (t5(0:mpnw+6), t6(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t7(0:mpnw+6), mpnw1)
!   t5 = mprealdp (dk)
  call mprealdp (dk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
  i1 = 1 - is
!   t6 = t5 ** i1
  call mpnpwr (t5(0:mpnw+6), i1, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!   s3 = t7 / t6
  call mpdiv (t7(0:mpnw+6), t6(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), s3(0:mpnw+6), mpnw1)
!   sum3 = sum3 + s3
  call mpadd (sum3(0:mpnw+6), s3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), sum3(0:mpnw+6), mpnw1)
!mp prec="mpnwm"
!   tc1 = abs (s1) - eps * abs (sum1)
  call mpabs (s1(0:mpnw+6), mpt1(0:mpnw+6), mpnwm)
  call mpabs (sum1(0:mpnw+6), mpt2(0:mpnw+6), mpnwm)
  call mpmul (eps(0:mpnwm+5), mpt2(0:mpnw+6), mpt3(0:mpnw+6), mpnwm)
  call mpsub (mpt1(0:mpnw+6), mpt3(0:mpnw+6), mpt4(0:mpnw+6), mpnwm)
  call mpeq (mpt4(0:mpnw+6), tc1(0:mpnwm+5), mpnwm)
!   tc2 = abs (s2) - eps * abs (sum2)
  call mpabs (s2(0:mpnw+6), mpt1(0:mpnw+6), mpnwm)
  call mpabs (sum2(0:mpnw+6), mpt2(0:mpnw+6), mpnwm)
  call mpmul (eps(0:mpnwm+5), mpt2(0:mpnw+6), mpt3(0:mpnw+6), mpnwm)
  call mpsub (mpt1(0:mpnw+6), mpt3(0:mpnw+6), mpt4(0:mpnw+6), mpnwm)
  call mpeq (mpt4(0:mpnw+6), tc2(0:mpnwm+5), mpnwm)
!   tc3 = abs (s3) - eps * abs (sum3)
  call mpabs (s3(0:mpnw+6), mpt1(0:mpnw+6), mpnwm)
  call mpabs (sum3(0:mpnw+6), mpt2(0:mpnw+6), mpnwm)
  call mpmul (eps(0:mpnwm+5), mpt2(0:mpnw+6), mpt3(0:mpnw+6), mpnwm)
  call mpsub (mpt1(0:mpnw+6), mpt3(0:mpnw+6), mpt4(0:mpnw+6), mpnwm)
  call mpeq (mpt4(0:mpnw+6), tc3(0:mpnwm+5), mpnwm)
!   ic1 = mpsgn (tc1)
  mpi1 = mpsgn (tc1(0:mpnwm+5))
  ic1 = mpi1
!   ic2 = mpsgn (tc2)
  mpi1 = mpsgn (tc2(0:mpnwm+5))
  ic2 = mpi1
!   ic3 = mpsgn (tc3)
  mpi1 = mpsgn (tc3(0:mpnwm+5))
  ic3 = mpi1
  if (ic1 < 0 .and. ic2 < 0 .and. ic3 < 0) goto 100
enddo

write (mpldb, 5)
5 format ('*** HURWITZZETAN: Loop end error')
! call mpabrt (565)
call mpabrt (565)

100 continue

!mp prec="mpnw1"
if (mod (is, 2) == 0) then
  i1 = is / 2
!   t2 = mppicon ** i1
  call mpnpwr (mppicon(0:), i1, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   t3 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!   t4 = ss - t3
  call mpsub (ss(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!   call mpgammar (ss1, t5, mpnw1)
  call mpgammar (ss1, t5, mpnw1)
!   t6 = t4 * t5
  call mpmul (t4(0:mpnw+6), t5(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!   t1 = t2 / t6
  call mpdiv (t2(0:mpnw+6), t6(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
else
  i1 = (is - 1) / 2
!   t2 = mppicon ** i1
  call mpnpwr (mppicon(0:), i1, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   t3 = sqrt (mppicon)
  call mpsqrt (mppicon(0:), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!   t4 = t2 * t3
  call mpmul (t2(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!   t2 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   t3 = ss - t2
  call mpsub (ss(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!   call mpgammar (ss1, t5, mpnw1)
  call mpgammar (ss1, t5, mpnw1)
!   t6 = t3 * t5
  call mpmul (t3(0:mpnw+6), t5(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!   t1 = t4 / t6
  call mpdiv (t4(0:mpnw+6), t6(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
endif

! t3 = mppicon ** is
call mpnpwr (mppicon(0:), is, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
! t4 = sqrt (mppicon)
call mpsqrt (mppicon(0:), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
! t2 = t3 / t4
call mpdiv (t3(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! t3 = sum1 * 0.5e0_mpdk
call mpmuld (sum1(0:mpnw+6), 0.5e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
! t4 = sum2 * 0.5e0_mpdk
call mpmuld (sum2(0:mpnw+6), 0.5e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
! t5 = sum3 * t2
call mpmul (sum3(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
! t6 = t1 + t3
call mpadd (t1(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
! t7 = t6 + t4
call mpadd (t6(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t7(0:mpnw+6), mpnw1)
! t1 = t7 + t5
call mpadd (t7(0:mpnw+6), t5(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!mp prec="mpnw"
! zz = t1
call mpeq (t1(0:mpnw+6), zz(0:), mpnw)

return
end subroutine mphurwitzzetan

subroutine mphurwitzzetanbe (nb1, nb2, berne, iss, aa, zz, mpnw)

!  This evaluates the Hurwitz zeta function, using the combination of
!  the definition formula (for large iss), and an Euler-Maclaurin scheme
!  (see formula 25.2.9 of the DLMF). The array berne contains precomputed
!  even Bernoulli numbers (see MPBERNER above). Its dimensions must be as
!  shown below. NB2 must be greater than 1.4 x precision in decimal digits.

implicit none
integer, intent(in):: nb1, mpnw
integer, intent(in):: nb2, iss
!mp dim1="0:nb1+5"
! type (mp_real), intent(in):: berne(1:nb2)
integer (mpiknd), intent(in):: berne(0:nb1+5,1:nb2)
!mp dim1="0:"
! type (mp_real), intent(in):: aa
integer (mpiknd), intent(in):: aa(0:)
! type (mp_real), intent(out):: zz
integer (mpiknd), intent(out):: zz(0:)
integer, parameter:: itrmax = 1000000
real (mpdknd), parameter:: dber = 0.4e0_mpdk, dcon = 0.18e0_mpdk
integer i1, i2, ic1, ic2, iqq, k, n1
real (mpdknd) d1, d2, dp
!mp dim1="0:mpnwm+5"
! type (mp_real) eps, tc1, tc2
integer (mpiknd) eps(0:mpnwm+5), tc1(0:mpnwm+5), tc2(0:mpnwm+5)
!mp dim1="0:mpnw+6"
! type (mp_real) aq, aq2, s1, s2, s3
integer (mpiknd) aq(0:mpnw+6), aq2(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6)
! type (mp_real) s4, t1, t2, t3, t4
integer (mpiknd) s4(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6)
! type (mp_real) t5, t6, f1
integer (mpiknd) t5(0:mpnw+6), t6(0:mpnw+6), f1(0:mpnw+6)
integer mpnw1

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpt3(0:mpnw+6)
integer (mpiknd) mpt4(0:mpnw+6), mpt5(0:mpnw+6), mpt6(0:mpnw+6)

call mpinitwds (mpt1, mpnw+6-5)
call mpinitwds (mpt2, mpnw+6-5)
call mpinitwds (mpt3, mpnw+6-5)
call mpinitwds (mpt4, mpnw+6-5)
call mpinitwds (mpt5, mpnw+6-5)
call mpinitwds (mpt6, mpnw+6-5)
call mpinitwds (eps(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc1(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc2(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (aq(0:mpnw+6), mpnw+6-5)
call mpinitwds (aq2(0:mpnw+6), mpnw+6-5)
call mpinitwds (s1(0:mpnw+6), mpnw+6-5)
call mpinitwds (s2(0:mpnw+6), mpnw+6-5)
call mpinitwds (s3(0:mpnw+6), mpnw+6-5)
call mpinitwds (s4(0:mpnw+6), mpnw+6-5)
call mpinitwds (t1(0:mpnw+6), mpnw+6-5)
call mpinitwds (t2(0:mpnw+6), mpnw+6-5)
call mpinitwds (t3(0:mpnw+6), mpnw+6-5)
call mpinitwds (t4(0:mpnw+6), mpnw+6-5)
call mpinitwds (t5(0:mpnw+6), mpnw+6-5)
call mpinitwds (t6(0:mpnw+6), mpnw+6-5)
call mpinitwds (f1(0:mpnw+6), mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)

! ic1 = mpspacer (zz)
ic1 = mpspacer (zz(0:))
if (mpnw < mpnwm .or. ic1 < mpnw + 6) then
  write (6, 1)
1 format ('*** HURWITZZETANBE: Uninitialized or inadequately sized array')
!   call mpabrt (566)
  call mpabrt (566)
endif

!mp prec="mpnwm"

! tc2 = mprealdp (2.e0_mpdk)
call mprealdp (2.e0_mpdk, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), tc2(0:mpnwm+5), mpnwm)
ic2 = -mpnw1*mpnbt
! eps = tc2 ** ic2
call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), eps(0:mpnwm+5), mpnwm)

!mp prec="mpnw1"

!   Check if berne array has been initialized.

! d1 = dpreal (berne(1))
call mpmdc (berne(0:nb1+5,1), mpd1, mpi1, mpnw1)
mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
d1 = mpd1
if (abs (d1 - 1.e0_mpdk / 6.e0_mpdk) > mprdfz .or. nb2 < int (dber * mpnbt * mpnw1)) then
  write (mpldb, 2) int (dber * mpnbt * mpnw1)
2 format ('*** HURWITZZETANBE: Array of even Bernoulli coefficients must be initialized'/ &
   'with at least',i8,' entries using BERNE or BERNER.')
!   call mpabrt (567)
  call mpabrt (567)
endif

! f1 = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), f1(0:mpnw+6), mpnw1)
! s1 = mprealdp (0.e0_mpdk)
call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), s1(0:mpnw+6), mpnw1)
! s2 = mprealdp (0.e0_mpdk)
call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), s2(0:mpnw+6), mpnw1)
! s3 = mprealdp (0.e0_mpdk)
call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), s3(0:mpnw+6), mpnw1)
! s4 = mprealdp (0.e0_mpdk)
call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), s4(0:mpnw+6), mpnw1)

if (iss <= 0) then
  write (mpldb, 3)
3 format ('*** HURWITZZETANBE: Integer argument <= 0')
!   call mpabrt (568)
  call mpabrt (568)
endif

! ic1 = mpsgn (aa)
mpi1 = mpsgn (aa(0:))
ic1 = mpi1
if (ic1 < 0) then
  write (mpldb, 4)
4 format ('*** HURWITZZETANBE: Argument < 0')
!   call mpabrt (569)
  call mpabrt (569)
endif

!   If iss > a certain value, then use definition formula.

if (iss > mpnbt * mpnw * log (2.e0_mpdk) / log (2.e0_mpdk * mpnbt * mpnw / 3.e0_mpdk)) then
  do k = 0, itrmax
    d2 = real (k, mpdk)
!     t1 = mprealdp (d2)
    call mprealdp (d2, mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!     t2 = aa + t1
    call mpadd (aa(0:), t1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!     t3 = t2 ** iss
    call mpnpwr (t2(0:mpnw+6), iss, mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!     t1 = f1 / t3
    call mpdiv (f1(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!     s1 = s1 + t1
    call mpadd (s1(0:mpnw+6), t1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), s1(0:mpnw+6), mpnw1)
!mp prec="mpnwm"
!     tc1 = abs (t1) - eps * abs (s1)
    call mpabs (t1(0:mpnw+6), mpt1(0:mpnw+6), mpnwm)
    call mpabs (s1(0:mpnw+6), mpt2(0:mpnw+6), mpnwm)
    call mpmul (eps(0:mpnwm+5), mpt2(0:mpnw+6), mpt3(0:mpnw+6), mpnwm)
    call mpsub (mpt1(0:mpnw+6), mpt3(0:mpnw+6), mpt4(0:mpnw+6), mpnwm)
    call mpeq (mpt4(0:mpnw+6), tc1(0:mpnwm+5), mpnwm)
!     ic1 = mpsgn (tc1)
    mpi1 = mpsgn (tc1(0:mpnwm+5))
    ic1 = mpi1
    if (ic1 < 0) goto 110
  enddo

  write (6, 5)
5 format ('*** HURWITZZETANBE: Loop end error 1')
!   call mpabrt (570)
  call mpabrt (570)
endif

!mp prec="mpnw1"
! d1 = dpreal (aa)
call mpmdc (aa(0:), mpd1, mpi1, mpnw1)
mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
d1 = mpd1
iqq = max (dcon * mpnw1 * mpnbt - d1, 0.e0_mpdk)

do k = 0, iqq - 1
  d2 = real (k, mpdk)
!   t1 = mprealdp (d2)
  call mprealdp (d2, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   t2 = aa + t1
  call mpadd (aa(0:), t1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   t3 = t2 ** iss
  call mpnpwr (t2(0:mpnw+6), iss, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!   t1 = f1 / t3
  call mpdiv (f1(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   s1 = s1 + t1
  call mpadd (s1(0:mpnw+6), t1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), s1(0:mpnw+6), mpnw1)
enddo

d2 = real (iqq, mpdk)
! t1 = mprealdp (d2)
call mprealdp (d2, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
! aq = aa + t1
call mpadd (aa(0:), t1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), aq(0:mpnw+6), mpnw1)
i1 = iss - 1
d2 = real (i1, mpdk)
! t1 = mprealdp (d2)
call mprealdp (d2, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
! t2 = aq ** i1
call mpnpwr (aq(0:mpnw+6), i1, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! t3 = t1 * t2
call mpmul (t1(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
! s2 = f1 / t3
call mpdiv (f1(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), s2(0:mpnw+6), mpnw1)
! t1 = aq ** iss
call mpnpwr (aq(0:mpnw+6), iss, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
! t2 = t1 * 2.e0_mpdk
call mpmuld (t1(0:mpnw+6), 2.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! s3 = f1 / t2
call mpdiv (f1(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), s3(0:mpnw+6), mpnw1)

d2 = real (iss, mpdk)
! t1 = mprealdp (d2)
call mprealdp (d2, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
! t2 = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
i1 = iss - 1
! t3 = aq ** i1
call mpnpwr (aq(0:mpnw+6), i1, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
! aq2 = aq ** 2
call mpnpwr (aq(0:mpnw+6), 2, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), aq2(0:mpnw+6), mpnw1)

do k = 1, nb2
  if (k > 1) then
    i1 = iss + 2*k - 3
    i2 = iss + 2*k - 2
!     t5 = t1 * real (i1, mpdk)
    call mpmuld (t1(0:mpnw+6), real (i1, mpdk), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
!     t1 = t5 * real (i2, mpdk)
    call mpmuld (t5(0:mpnw+6), real (i2, mpdk), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
  endif
  i1 = 2*k - 1
  i2 = 2*k
!   t5 = t2 * real (i1, mpdk)
  call mpmuld (t2(0:mpnw+6), real (i1, mpdk), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
!   t2 = t5 * real (i2, mpdk)
  call mpmuld (t5(0:mpnw+6), real (i2, mpdk), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   t3 = t3 * aq2
  call mpmul (t3(0:mpnw+6), aq2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!   t5 = berne(k) * t1
  call mpmul (berne(0:nb1+5,k), t1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
!   t6 = t2 * t3
  call mpmul (t2(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!   t4 = t5 / t6
  call mpdiv (t5(0:mpnw+6), t6(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!   s4 = s4 + t4
  call mpadd (s4(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), s4(0:mpnw+6), mpnw1)
!mp prec="mpnwm"
!   tc2 = abs (t4) - eps * abs (s4)
  call mpabs (t4(0:mpnw+6), mpt1(0:mpnw+6), mpnwm)
  call mpabs (s4(0:mpnw+6), mpt2(0:mpnw+6), mpnwm)
  call mpmul (eps(0:mpnwm+5), mpt2(0:mpnw+6), mpt3(0:mpnw+6), mpnwm)
  call mpsub (mpt1(0:mpnw+6), mpt3(0:mpnw+6), mpt4(0:mpnw+6), mpnwm)
  call mpeq (mpt4(0:mpnw+6), tc2(0:mpnwm+5), mpnwm)
!   ic2 = mpsgn (tc2)
  mpi1 = mpsgn (tc2(0:mpnwm+5))
  ic2 = mpi1
  if (ic2 < 0) goto 110
enddo

write (6, 6)
6 format ('*** HURWITZZETANBE: End loop error 2; call BERNE with larger NB.')
! call mpabrt (571)
call mpabrt (571)

110 continue

!mp prec="mpnw1"
! t5 = s1 + s2
call mpadd (s1(0:mpnw+6), s2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
! t6 = t5 + s3
call mpadd (t5(0:mpnw+6), s3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
! s1 = t6 + s4
call mpadd (t6(0:mpnw+6), s4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), s1(0:mpnw+6), mpnw1)
!mp prec="mpnw"
! zz = s1
call mpeq (s1(0:mpnw+6), zz(0:), mpnw)

return
end subroutine mphurwitzzetanbe

subroutine mphypergeompfq (nw, np, nq, aa, bb, xx, yy, mpnw)

!  This returns the HypergeometricPFQ function, namely the sum of the infinite series

!  Sum_0^infinity poch(aa(1),n)*poch(aa(2),n)*...*poch(aa(np),n) /
!      poch(bb(1),n)*poch(bb(2),n)*...*poch(bb(nq),n) * xx^n / n!

!  This subroutine evaluates the HypergeometricPFQ function directly according to
!  this definitional formula. The arrays aa and bb must be dimensioned as shown below.
!  NP and NQ are limited to [1,10].

implicit none
integer, intent(in):: nw, mpnw
integer, intent(in):: np, nq
!mp dim1="0:nw+5"
! type (mp_real), intent(in):: aa(1:np), bb(1:nq)
integer (mpiknd), intent(in):: aa(0:nw+5,1:np), bb(0:nw+5,1:nq)
!mp dim1="0:"
! type (mp_real), intent(in):: xx
integer (mpiknd), intent(in):: xx(0:)
! type (mp_real), intent(out):: yy
integer (mpiknd), intent(out):: yy(0:)
integer, parameter:: itrmax = 1000000, npq = 10
integer i1, i2, ic1, ic2, j, k
real (mpdknd) d1
!mp dim1="0:mpnwm+5"
! type (mp_real) eps, tc1, tc2
integer (mpiknd) eps(0:mpnwm+5), tc1(0:mpnwm+5), tc2(0:mpnwm+5)
!mp dim1="0:mpnw+6"
! type (mp_real) sum, td, tn, t1, t2
integer (mpiknd) sum(0:mpnw+6), td(0:mpnw+6), tn(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6)
! type (mp_real) t3, t4
integer (mpiknd) t3(0:mpnw+6), t4(0:mpnw+6)
integer mpnw1

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpt3(0:mpnw+6)
integer (mpiknd) mpt4(0:mpnw+6), mpt5(0:mpnw+6), mpt6(0:mpnw+6)

call mpinitwds (mpt1, mpnw+6-5)
call mpinitwds (mpt2, mpnw+6-5)
call mpinitwds (mpt3, mpnw+6-5)
call mpinitwds (mpt4, mpnw+6-5)
call mpinitwds (mpt5, mpnw+6-5)
call mpinitwds (mpt6, mpnw+6-5)
call mpinitwds (eps(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc1(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc2(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (sum(0:mpnw+6), mpnw+6-5)
call mpinitwds (td(0:mpnw+6), mpnw+6-5)
call mpinitwds (tn(0:mpnw+6), mpnw+6-5)
call mpinitwds (t1(0:mpnw+6), mpnw+6-5)
call mpinitwds (t2(0:mpnw+6), mpnw+6-5)
call mpinitwds (t3(0:mpnw+6), mpnw+6-5)
call mpinitwds (t4(0:mpnw+6), mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)

! ic1 = mpspacer (yy)
ic1 = mpspacer (yy(0:))
if (mpnw < mpnwm .or. ic1 < mpnw + 6) then
  write (6, 1)
1 format ('*** HYPERGEOMPFQ: Uninitialized or inadequately sized array')
!   call mpabrt (572)
  call mpabrt (572)
endif

!mp prec="mpnwm"

! tc2 = mprealdp (2.e0_mpdk)
call mprealdp (2.e0_mpdk, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), tc2(0:mpnwm+5), mpnwm)
ic2 = -mpnw1*mpnbt
! eps = tc2 ** ic2
call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), eps(0:mpnwm+5), mpnwm)

if (np < 1 .or. np > npq .or. nq < 1 .or. nq > npq) then
  write (mpldb, 2) npq
2 format ('*** HYPERGEOMPFQ: NP and NQ must be between 1 and',i4)
!   call mpabrt (573)
  call mpabrt (573)
endif

!mp prec="mpnw1"

! sum = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), sum(0:mpnw+6), mpnw1)
! td = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), td(0:mpnw+6), mpnw1)
! tn = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), tn(0:mpnw+6), mpnw1)

do k = 1, itrmax
  i1 = k - 1
  d1 = real (i1, mpdk)
!   t1 = mprealdp (d1)
  call mprealdp (d1, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)

  do j = 1, np
!     t2 = t1 + aa(j)
    call mpadd (t1(0:mpnw+6), aa(0:nw+5,j), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!     tn = tn * t2
    call mpmul (tn(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), tn(0:mpnw+6), mpnw1)
  enddo

  do j = 1, nq
!     t2 = t1 + bb(j)
    call mpadd (t1(0:mpnw+6), bb(0:nw+5,j), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!     td = td * t2
    call mpmul (td(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), td(0:mpnw+6), mpnw1)
  enddo

!   tn = tn * xx
  call mpmul (tn(0:mpnw+6), xx(0:), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), tn(0:mpnw+6), mpnw1)
!   td = td * real (k, mpdk)
  call mpmuld (td(0:mpnw+6), real (k, mpdk), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), td(0:mpnw+6), mpnw1)
!   t1 = tn / td
  call mpdiv (tn(0:mpnw+6), td(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   sum = sum + t1
  call mpadd (sum(0:mpnw+6), t1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), sum(0:mpnw+6), mpnw1)
!mp prec="mpnwm"
!   tc1 = abs (t1) - eps * sum
  call mpabs (t1(0:mpnw+6), mpt1(0:mpnw+6), mpnwm)
  call mpmul (eps(0:mpnwm+5), sum(0:mpnw+6), mpt2(0:mpnw+6), mpnwm)
  call mpsub (mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpt3(0:mpnw+6), mpnwm)
  call mpeq (mpt3(0:mpnw+6), tc1(0:mpnwm+5), mpnwm)
!   ic1 = mpsgn (tc1)
  mpi1 = mpsgn (tc1(0:mpnwm+5))
  ic1 = mpi1
  if (ic1 < 0) goto 100
enddo

write (mpldb, 3) itrmax
3 format ('*** HYPERGEOMPFQ: Loop end error',i10)
! call mpabrt (574)
call mpabrt (574)

100 continue

!mp prec="mpnw"
! yy = sum
call mpeq (sum(0:mpnw+6), yy(0:), mpnw)

return
end subroutine mphypergeompfq

subroutine mpincgammar (s, z, g, mpnw)

!  This returns the incomplete gamma function, using a combination of formula
!  8.7.3 of the DLMF (for modest-sized z), formula 8.11.2 (for large z),
!  a formula from the Wikipedia page for the case S = 0, and another formula
!  from the Wikipedia page for the case S = negative integer. The formula
!  for the case S = 0 requires increased working precision, up to 2.5X normal,
!  depending on the size of Z.

implicit none
integer, intent(in):: mpnw
!mp dim1="0:"
! type (mp_real), intent(in):: s, z
integer (mpiknd), intent(in):: s(0:), z(0:)
! type (mp_real), intent(out):: g
integer (mpiknd), intent(out):: g(0:)
integer, parameter:: itrmax = 1000000, ipx = 3
real (mpdknd), parameter:: dfrac1 = 0.833e0_mpdk, dfrac2 = 1.46e0_mpdk
integer i1, ic1, ic2, ic3, k, nn, n1, n2
real (mpdknd) d1, d2, d3
!mp dim1="0:mpnwm+5"
! type (mp_real) eps, f1, tc1, tc2, tc3
integer (mpiknd) eps(0:mpnwm+5), f1(0:mpnwm+5), tc1(0:mpnwm+5), tc2(0:mpnwm+5), tc3(0:mpnwm+5)
!mp dim1="0:ipx*mpnw+6"
! type (mp_real) t0, t1, t2, t3
integer (mpiknd) t0(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), t3(0:ipx*mpnw+6)
! type (mp_real) t4, t5
integer (mpiknd) t4(0:ipx*mpnw+6), t5(0:ipx*mpnw+6)
integer mpnw1, mpnw2

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6)
integer (mpiknd) mpt4(0:ipx*mpnw+6), mpt5(0:ipx*mpnw+6), mpt6(0:ipx*mpnw+6)

call mpinitwds (mpt1, ipx*mpnw+6-5)
call mpinitwds (mpt2, ipx*mpnw+6-5)
call mpinitwds (mpt3, ipx*mpnw+6-5)
call mpinitwds (mpt4, ipx*mpnw+6-5)
call mpinitwds (mpt5, ipx*mpnw+6-5)
call mpinitwds (mpt6, ipx*mpnw+6-5)
call mpinitwds (eps(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (f1(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc1(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc2(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc3(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (t0(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t2(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t3(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t4(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t5(0:ipx*mpnw+6), ipx*mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)
mpnw2 = min (ipx*mpnw + 1, mpnwx)

! ic1 = mpspacer (g)
ic1 = mpspacer (g(0:))
if (mpnw < mpnwm .or. ic1 < mpnw + 6) then
  write (6, 1)
1 format ('*** INCGAMMAR: Uninitialized or inadequately sized array')
!   call mpabrt (575)
  call mpabrt (575)
endif

!mp prec="mpnwm"

! tc2 = mprealdp (2.e0_mpdk)
call mprealdp (2.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnwm)
call mpeq (mpt1(0:ipx*mpnw+6), tc2(0:mpnwm+5), mpnwm)
ic2 = -mpnw1*mpnbt
! eps = tc2 ** ic2
call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:ipx*mpnw+6), mpnwm)
call mpeq (mpt1(0:ipx*mpnw+6), eps(0:mpnwm+5), mpnwm)
! f1 = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnwm)
call mpeq (mpt1(0:ipx*mpnw+6), f1(0:mpnwm+5), mpnwm)

!mp prec="mpnw1"

! n1 = mpsgn (s)
mpi1 = mpsgn (s(0:))
n1 = mpi1
! n2 = mpsgn (z)
mpi1 = mpsgn (z(0:))
n2 = mpi1
if (n2 == 0 .or. n1 /= 0 .and. n2 < 0) then
  write (mpldb, 2)
2 format ('*** INCGAMMAR: The second argument must not be zero,'/ &
    'and must not be negative unless the first is zero.')
!   call mpabrt (576)
  call mpabrt (576)
endif

! d1 = dpreal (z)
call mpmdc (z(0:), mpd1, mpi1, mpnw1)
mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
d1 = mpd1

if (abs (d1) < dfrac1 * mpnw1 * mpnbt) then

!   This is for modest-sized z.

!   t1 = aint (s)
  call mpinfr (s(0:), mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
!   t2 = s - t1
  call mpsub (s(0:), t1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!   ic1 = mpsgn (t2)
  mpi1 = mpsgn (t2(0:ipx*mpnw+6))
  ic1 = mpi1
!   d2 = dpreal (s)
  call mpmdc (s(0:), mpd1, mpi1, mpnw1)
  mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
  d2 = mpd1
  nn = nint (d2)

  if (ic1 == 0 .and. nn == 1) then
!     t2 = - z
    call mpneg (z(0:), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!     t1 = exp (t2)
    call mpexp (t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
    goto 200
  elseif (ic1 == 0 .and. nn <= 0) then

!    S is zero or a negative integer -- use a different algorithm. In
!    either event, first compute incgamma for S = 0. For large Z, the
!    working precision must be increased.

    mpnw2 = min (mpnw + nint (dfrac2 * d1 / mpnbt) + 1, mpnwx)
    if (mpnw2 > ipx*mpnw+1) then
      write (6, 3) mpnw2
3     format ('*** INCGAMMAR: Inadequate working precision:',i6)
!       call mpabrt (577)
      call mpabrt (577)
    endif

!mp prec="mpnwm"
!     tc2 = mprealdp (2.e0_mpdk)
    call mprealdp (2.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnwm)
    call mpeq (mpt1(0:ipx*mpnw+6), tc2(0:mpnwm+5), mpnwm)
    ic2 = -mpnw2*mpnbt
!     eps = tc2 ** ic2
    call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:ipx*mpnw+6), mpnwm)
    call mpeq (mpt1(0:ipx*mpnw+6), eps(0:mpnwm+5), mpnwm)
!mp prec="mpnw2"
!     t0 = z
    call mpeq (z(0:), t0(0:ipx*mpnw+6), mpnw2)
!     t1 = z
    call mpeq (z(0:), t1(0:ipx*mpnw+6), mpnw2)
!     t2 = mprealdp (1.e0_mpdk)
    call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)

    do k = 2, itrmax
      if (mod (k, 2) == 1) then
        d1 = real (k, mpdk)
!         t3 = f1 / d1
        call mpdivd (f1(0:mpnwm+5), d1, mpt1(0:ipx*mpnw+6), mpnw2)
        call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!         t2 = t2 + t3
        call mpadd (t2(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
        call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
      endif
!       t3 = z * t1
      call mpmul (z(0:), t1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
      call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
      d1 = real (2*k, mpdk)
!       t1 = t3 / d1
      call mpdivd (t3(0:ipx*mpnw+6), d1, mpt1(0:ipx*mpnw+6), mpnw2)
      call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw2)
!       t3 = t1 * t2
      call mpmul (t1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
      call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!       t0 = t0 + t3
      call mpadd (t0(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
      call mpeq (mpt1(0:ipx*mpnw+6), t0(0:ipx*mpnw+6), mpnw2)
!mp prec="mpnwm"
!       tc1 = abs (t3) - eps * abs (t0)
      call mpabs (t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnwm)
      call mpabs (t0(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnwm)
      call mpmul (eps(0:mpnwm+5), mpt2(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpnwm)
      call mpsub (mpt1(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpt4(0:ipx*mpnw+6), mpnwm)
      call mpeq (mpt4(0:ipx*mpnw+6), tc1(0:mpnwm+5), mpnwm)
!       ic1 = mpsgn (tc1)
      mpi1 = mpsgn (tc1(0:mpnwm+5))
      ic1 = mpi1
      if (ic1 < 0) goto 100
    enddo

    write (mpldb, 4)
4   format ('*** INCGAMMAR: Loop end error 1')
!     call mpabrt (578)
    call mpabrt (578)

100  continue

!mp prec="mpnw2"
!     t1 = - mpegammacon
    call mpneg (mpegammacon(0:), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw2)
!     t3 = abs (z)
    call mpabs (z(0:), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!     t2 = log (t3)
    call mplog (t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!     t3 = t1 - t2
    call mpsub (t1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!     t4 = - z * 0.5e0_mpdk
    call mpmuld (z(0:), 0.5e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
    call mpneg (mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt2(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpnw2)
!     t5 = exp (t4)
    call mpexp (t4(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), mpnw2)
!     t4 = t5 * t0
    call mpmul (t5(0:ipx*mpnw+6), t0(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpnw2)
!     t1 = t3 + t4
    call mpadd (t3(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw2)
    if (nn == 0) goto 200

!   S is negative integer (not zero).

    nn = abs (nn)
!     t0 = mprealdp (1.e0_mpdk)
    call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t0(0:ipx*mpnw+6), mpnw2)
!     t2 = t0
    call mpeq (t0(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)

    do k = 1, nn - 1
!       t0 = t0 * real (k, mpdk)
      call mpmuld (t0(0:ipx*mpnw+6), real (k, mpdk), mpt1(0:ipx*mpnw+6), mpnw2)
      call mpeq (mpt1(0:ipx*mpnw+6), t0(0:ipx*mpnw+6), mpnw2)
!       t2 = t0
      call mpeq (t0(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
    enddo

!     t5 = t0 * real (nn, mpdk)
    call mpmuld (t0(0:ipx*mpnw+6), real (nn, mpdk), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), mpnw2)

    do k = 1, nn - 1
!       t3 = t2 * z
      call mpmul (t2(0:ipx*mpnw+6), z(0:), mpt1(0:ipx*mpnw+6), mpnw2)
      call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
      i1 = nn - k
!       t4 = t3 / real (i1, mpdk)
      call mpdivd (t3(0:ipx*mpnw+6), real (i1, mpdk), mpt1(0:ipx*mpnw+6), mpnw2)
      call mpeq (mpt1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpnw2)
!       t2 = - t4
      call mpneg (t4(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
      call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!       t0 = t0 + t2
      call mpadd (t0(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
      call mpeq (mpt1(0:ipx*mpnw+6), t0(0:ipx*mpnw+6), mpnw2)
    enddo

!     t2 = exp (z)
    call mpexp (z(0:), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!     t3 = t0 / t2
    call mpdiv (t0(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!     t4 = z ** nn
    call mpnpwr (z(0:), nn, mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpnw2)
!     t2 = t3 / t4
    call mpdiv (t3(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)

    if (mod (nn, 2) == 0) then
!       t3 = t2 + t1
      call mpadd (t2(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
      call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
    else
!       t3 = t2 - t1
      call mpsub (t2(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
      call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
    endif
!     t1 = t3 / t5
    call mpdiv (t3(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw2)
    goto 200
  endif

  mpnw2 = min (mpnw + nint (dfrac2 * d1 / mpnbt) + 1, mpnwx)
  if (mpnw2 > ipx*mpnw+1) then
    write (6, 5) mpnw2
5   format ('*** INCGAMMAR: Inadequate working precision:',i6)
!     call mpabrt (579)
    call mpabrt (579)
  endif

!mp prec="mpnwm"
!     tc2 = mprealdp (2.e0_mpdk)
    call mprealdp (2.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnwm)
    call mpeq (mpt1(0:ipx*mpnw+6), tc2(0:mpnwm+5), mpnwm)
    ic2 = -mpnw2*mpnbt
!     eps = tc2 ** ic2
    call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:ipx*mpnw+6), mpnwm)
    call mpeq (mpt1(0:ipx*mpnw+6), eps(0:mpnwm+5), mpnwm)
!mp prec="mpnw2"
!   call mpgammar (s, t1, mpnw2)
  call mpgammar (s, t1, mpnw2)
!   t3 = s * t1
  call mpmul (s(0:), t1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!   t2 = f1 / t3
  call mpdiv (f1(0:mpnwm+5), t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!   t0 = t2
  call mpeq (t2(0:ipx*mpnw+6), t0(0:ipx*mpnw+6), mpnw2)

  do k = 1, itrmax
!     t5 = t2  * z
    call mpmul (t2(0:ipx*mpnw+6), z(0:), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), mpnw2)
    d3 = real (k, mpdk)
!     t3 = mprealdp (d3)
    call mprealdp (d3, mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!     t4 = s + t3
    call mpadd (s(0:), t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpnw2)
!     t2 = t5 / t4
    call mpdiv (t5(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!     t0 = t0 + t2
    call mpadd (t0(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
    call mpeq (mpt1(0:ipx*mpnw+6), t0(0:ipx*mpnw+6), mpnw2)
!mp prec="mpnwm"
!     tc2 = abs (t2) - eps * abs (t0)
    call mpabs (t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnwm)
    call mpabs (t0(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnwm)
    call mpmul (eps(0:mpnwm+5), mpt2(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpnwm)
    call mpsub (mpt1(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpt4(0:ipx*mpnw+6), mpnwm)
    call mpeq (mpt4(0:ipx*mpnw+6), tc2(0:mpnwm+5), mpnwm)
!     ic2 = mpsgn (tc2)
    mpi1 = mpsgn (tc2(0:mpnwm+5))
    ic2 = mpi1
    if (ic2 < 0) goto 110
  enddo

  write (mpldb, 6) itrmax
6   format ('*** INCGAMMAR: Loop end error 1')
!   call mpabrt (580)
  call mpabrt (580)

110 continue

!mp prec="mpnw2"
!   t2 = z ** s
  call mppower (z(0:), s(0:), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!   t3 = exp (z)
  call mpexp (z(0:), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!   t4 = t2 / t3
  call mpdiv (t2(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpnw2)
!   t5 = t4 * t0
  call mpmul (t4(0:ipx*mpnw+6), t0(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t5(0:ipx*mpnw+6), mpnw2)
!   t2 = f1 - t5
  call mpsub (f1(0:mpnwm+5), t5(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!   t1 = t1 * t2
  call mpmul (t1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw2)
  goto 200
else

!   This is for large z. Note that if S is a positive integer, this loop
!   is finite.

!mp prec="mpnw1"
!   t0 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t0(0:ipx*mpnw+6), mpnw1)
!   t1 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)

  do k = 1, itrmax
    d3 = real (k, mpdk)
!     t2 = mprealdp (d3)
    call mprealdp (d3, mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!     t3 = s - t2
    call mpsub (s(0:), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw1)
!     t4 = t1 * t3
    call mpmul (t1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpnw1)
!     t1 = t4 / z
    call mpdiv (t4(0:ipx*mpnw+6), z(0:), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
!     t0 = t0 + t1
    call mpadd (t0(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt1(0:ipx*mpnw+6), t0(0:ipx*mpnw+6), mpnw1)
!     tc3 = abs (t1) - eps * abs (t0)
    call mpabs (t1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
    call mpabs (t0(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnw1)
    call mpmul (eps(0:mpnwm+5), mpt2(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpnw1)
    call mpsub (mpt1(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpt4(0:ipx*mpnw+6), mpnw1)
    call mpeq (mpt4(0:ipx*mpnw+6), tc3(0:mpnwm+5), mpnw1)
!     ic3 = mpsgn (tc3)
    mpi1 = mpsgn (tc3(0:mpnwm+5))
    ic3 = mpi1
    if (ic3 < 0) goto 120
  enddo

  write (mpldb, 7)
7 format ('*** INCGAMMAR: Loop end error 3')
!   call mpabrt (581)
  call mpabrt (581)

120 continue

!   t2 = s - f1
  call mpsub (s(0:), f1(0:mpnwm+5), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!   t3 = z ** t2
  call mppower (z(0:), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw1)
!   t4 = exp (z)
  call mpexp (z(0:), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpnw1)
!   t2 = t3 / t4
  call mpdiv (t3(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw1)
!   t1 = t2 * t0
  call mpmul (t2(0:ipx*mpnw+6), t0(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw1)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
  goto 200
endif

200 continue

!mp prec="mpnw"
! g = t1
call mpeq (t1(0:ipx*mpnw+6), g(0:), mpnw)

return
end subroutine mpincgammar

subroutine mppolygamma (nn, x, y, mpnw)

!   This returns polygamma (nn, x) for nn >= 0 and 0 < x < 1, by calling
!   mphurwitzzetan.

implicit none
integer, intent(in):: mpnw
integer, intent(in):: nn
!mp dim1="0:"
! type (mp_real), intent(in):: x
integer (mpiknd), intent(in):: x(0:)
! type (mp_real), intent(out):: y
integer (mpiknd), intent(out):: y(0:)
integer i1, ic1, ic2, k
!mp dim1="0:mpnw+6"
! type (mp_real) t1, t2, t3, t4
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6)
integer mpnw1

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpt3(0:mpnw+6)
integer (mpiknd) mpt4(0:mpnw+6), mpt5(0:mpnw+6), mpt6(0:mpnw+6)

call mpinitwds (mpt1, mpnw+6-5)
call mpinitwds (mpt2, mpnw+6-5)
call mpinitwds (mpt3, mpnw+6-5)
call mpinitwds (mpt4, mpnw+6-5)
call mpinitwds (mpt5, mpnw+6-5)
call mpinitwds (mpt6, mpnw+6-5)
call mpinitwds (t1(0:mpnw+6), mpnw+6-5)
call mpinitwds (t2(0:mpnw+6), mpnw+6-5)
call mpinitwds (t3(0:mpnw+6), mpnw+6-5)
call mpinitwds (t4(0:mpnw+6), mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)

! ic1 = mpspacer (y)
ic1 = mpspacer (y(0:))
if (mpnw < mpnwm .or. ic1 < mpnw + 6) then
  write (6, 1)
1 format ('*** POLYGAMMA: Uninitialized or inadequately sized array')
!   call mpabrt (582)
  call mpabrt (582)
endif

!mp prec="mpnw1"

if (nn <= 0) then
  write (mpldb, 2)
2 format ('*** POLYGAMMA: Integer argument <= 0')
!   call mpabrt (583)
  call mpabrt (583)
endif

! t1 = mprealdp (0.e0_mpdk)
call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
! t2 = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! t3 = x - t1
call mpsub (x(0:), t1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
! ic1 = mpsgn (t3)
mpi1 = mpsgn (t3(0:mpnw+6))
ic1 = mpi1
! t4 = x - t2
call mpsub (x(0:), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
! ic2 = mpsgn (t4)
mpi1 = mpsgn (t4(0:mpnw+6))
ic2 = mpi1

if (ic1 <= 0 .or. ic2 >= 0) then
  write (mpldb, 3)
3 format ('*** POLYGAMMA: Argument must be in the range (0, 1)')
!   call mpabrt (584)
  call mpabrt (584)
endif

! t1 = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)

do k = 1, nn
!   t1 = t1 * real (k, mpdk)
  call mpmuld (t1(0:mpnw+6), real (k, mpdk), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
enddo

if (mod (nn + 1, 2) == 1) then
!   t1 = - t1
  call mpneg (t1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
endif

i1 = nn + 1
! call mphurwitzzetan (i1, x, t2, mpnw1)
call mphurwitzzetan (i1, x, t2, mpnw1)
! t3 = t1 * t2
call mpmul (t1(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!mp prec="mpnw"
! y = t3
call mpeq (t3(0:mpnw+6), y(0:), mpnw)

return
end subroutine mppolygamma

subroutine mppolygammabe (nb1, nb2, berne, nn, x, y, mpnw)

!  This returns polygamma (nn, x) for nn >= 0, by calling mphurwitzzetanbe.
!  The array berne contains precomputed even Bernoulli numbers (see MPBERNER
!  above). Its dimensions must be as shown below. NB2 must be greater than
!  1.4 x precision in decimal digits.

implicit none
integer, intent(in):: nb1, mpnw
integer, intent(in):: nb2, nn
!mp dim1="0:nb1+5"
! type (mp_real), intent(in):: berne(1:nb2)
integer (mpiknd), intent(in):: berne(0:nb1+5,1:nb2)
!mp dim1="0:"
! type (mp_real), intent(in):: x
integer (mpiknd), intent(in):: x(0:)
! type (mp_real), intent(out):: y
integer (mpiknd), intent(out):: y(0:)
real (mpdknd), parameter:: dber = 0.4e0_mpdk
integer ic1, i1, i2, k, n1
real (mpdknd) d1
!mp dim1="0:mpnw+6"
! type (mp_real) t1, t2, t3
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6)
integer mpnw1

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpt3(0:mpnw+6)
integer (mpiknd) mpt4(0:mpnw+6), mpt5(0:mpnw+6), mpt6(0:mpnw+6)

call mpinitwds (mpt1, mpnw+6-5)
call mpinitwds (mpt2, mpnw+6-5)
call mpinitwds (mpt3, mpnw+6-5)
call mpinitwds (mpt4, mpnw+6-5)
call mpinitwds (mpt5, mpnw+6-5)
call mpinitwds (mpt6, mpnw+6-5)
call mpinitwds (t1(0:mpnw+6), mpnw+6-5)
call mpinitwds (t2(0:mpnw+6), mpnw+6-5)
call mpinitwds (t3(0:mpnw+6), mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)

! ic1 = mpspacer (y)
ic1 = mpspacer (y(0:))
if (mpnw < mpnwm .or. ic1 < mpnw + 6) then
  write (6, 1)
1 format ('*** POLYGAMMABE: Uninitialized or inadequately sized array')
!   call mpabrt (585)
  call mpabrt (585)
endif

!mp prec="mpnw1"

!   Check if berne array has been initialized.

! d1 = dpreal (berne(1))
call mpmdc (berne(0:nb1+5,1), mpd1, mpi1, mpnw1)
mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
d1 = mpd1
if (abs (d1 - 1.e0_mpdk / 6.e0_mpdk) > mprdfz .or. nb2 < int (dber * mpnbt * mpnw1)) then
  write (mpldb, 2) int (dber * mpnbt * mpnw1)
2 format ('*** POLYGAMMABE: Array of even Bernoulli coefficients must be initialized'/ &
   'with at least',i8,' entries using BERNE or BERNER.')
!   call mpabrt (586)
  call mpabrt (586)
endif

if (nn <= 0) then
  write (mpldb, 3)
3 format ('*** POLYGAMMABE: Integer argument <= 0')
!   call mpabrt (587)
  call mpabrt (587)
endif

! ic1 = mpsgn (x)
mpi1 = mpsgn (x(0:))
ic1 = mpi1
if (ic1 < 0) then
  write (mpldb, 4)
4 format ('*** POLYGAMMABE: Argument < 0')
!   call mpabrt (588)
  call mpabrt (588)
endif

! t1 = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)

do k = 1, nn
!   t1 = t1 * real (k, mpdk)
  call mpmuld (t1(0:mpnw+6), real (k, mpdk), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
enddo

if (mod (nn + 1, 2) == 1) then
!   t1 = - t1
  call mpneg (t1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
endif

i1 = nn + 1
! call mphurwitzzetanbe (nb1, nb2, berne, i1, x, t2, mpnw1)
call mphurwitzzetanbe (nb1, nb2, berne, i1, x, t2, mpnw1)
! t3 = t1 * t2
call mpmul (t1(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!mp prec="mpnw"
! y = t3
call mpeq (t3(0:mpnw+6), y(0:), mpnw)

return
end subroutine mppolygammabe

subroutine mppolylogini (na, nn, arr, mpnw)

!   Initializes the MP array arr with data for mppolylogneg.
!   NN must be in the range (-nmax, -1).

implicit none
integer, intent(in):: na, mpnw
integer, intent(in):: nn
!mp dim1="0:na+5"
! type (mp_real), intent(out):: arr(1:abs(nn))
integer (mpiknd), intent(out):: arr(0:na+5,1:abs(nn))
integer, parameter:: nmax = 1000
integer ic1, ic2, i1, i2, i3, k, n, nna
!mp dim1="0:mpnw+6"
! type (mp_real) aa(1:2,1:abs(nn)), t1, t2
integer (mpiknd) aa(0:mpnw+6,1:2,1:abs(nn)), t1(0:mpnw+6), t2(0:mpnw+6)
integer mpnw1

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpt3(0:mpnw+6)
integer (mpiknd) mpt4(0:mpnw+6), mpt5(0:mpnw+6), mpt6(0:mpnw+6)

call mpinitwds (mpt1, mpnw+6-5)
call mpinitwds (mpt2, mpnw+6-5)
call mpinitwds (mpt3, mpnw+6-5)
call mpinitwds (mpt4, mpnw+6-5)
call mpinitwds (mpt5, mpnw+6-5)
call mpinitwds (mpt6, mpnw+6-5)
do mpi1 = 1,2; do mpi2 = 1,abs(nn)
  call mpinitwds (aa(0:mpnw+6,mpi1,mpi2), mpnw+6-5)
enddo; enddo
call mpinitwds (t1(0:mpnw+6), mpnw+6-5)
call mpinitwds (t2(0:mpnw+6), mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)
nna = abs (nn)

do ic1 = 1, nna
!   ic2 = mpspacer (arr(ic1))
  ic2 = mpspacer (arr(0:na+5,ic1))
  if (mpnw < mpnwm .or. ic2 < mpnw + 6) then
    write (6, 1)
1   format ('*** POLYLOGINI: Uninitialized or inadequately sized array')
!     call mpabrt (589)
    call mpabrt (589)
  endif
enddo

!mp prec="mpnw1"

i1 = 2
i2 = 1
! aa(1,1) = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), aa(0:mpnw+6,1,1), mpnw1)
! aa(2,1) = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), aa(0:mpnw+6,2,1), mpnw1)

do k = 2, nna
!   aa(1,k) = mprealdp (0.e0_mpdk)
  call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), aa(0:mpnw+6,1,k), mpnw1)
!   aa(2,k) = mprealdp (0.e0_mpdk)
  call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), aa(0:mpnw+6,2,k), mpnw1)
enddo

do n = 2, nna
  i1 = 3 - i1
  i2 = 3 - i1

  do k = 2, n
    i3 = n + 1 - k
!     t1 = aa(i1,k-1) * real (i3, mpdk)
    call mpmuld (aa(0:mpnw+6,i1,k-1), real (i3, mpdk), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!     t2 = aa(i1,k) * real (k, mpdk)
    call mpmuld (aa(0:mpnw+6,i1,k), real (k, mpdk), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!     aa(i2,k) = t1 + t2
    call mpadd (t1(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), aa(0:mpnw+6,i2,k), mpnw1)
  enddo
enddo

!mp prec="mpnw"
do k = 1, nna
!   arr(k) = aa(i2,k)
  call mpeq (aa(0:mpnw+6,i2,k), arr(0:na+5,k), mpnw)
enddo

return
end subroutine mppolylogini

subroutine mppolylogneg (na, nn, arr, x, y, mpnw)

!   This returns polylog (nn, x) for the case nn < 0. Before calling this,
!   one must call mppolylognini to initialize the array arr for this NN.
!   The dimensions of arr must be as shown below.
!   NN must be in the range (-nmax, -1).

implicit none
integer, intent(in):: na, mpnw
integer, intent(in):: nn
!mp dim1="0:na+5"
! type (mp_real), intent(in):: arr(1:abs(nn))
integer (mpiknd), intent(in):: arr(0:na+5,1:abs(nn))
!mp dim1="0:"
! type (mp_real), intent(in):: x
integer (mpiknd), intent(in):: x(0:)
! type (mp_real), intent(out):: y
integer (mpiknd), intent(out):: y(0:)
integer, parameter:: ipx = 3, nmax = 1000
real (mpdknd), parameter:: al2 = 0.69314718055994530941723212145817657e0_mpdk
integer ic1, i1, i2, k, n1, n2, nna
real (mpdknd) d1, d2
!mp dim1="0:ipx*mpnw+6"
! type (mp_real) t1, t2, t3, t4
integer (mpiknd) t1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), t4(0:ipx*mpnw+6)
integer mpnw1, mpnw2

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6)
integer (mpiknd) mpt4(0:ipx*mpnw+6), mpt5(0:ipx*mpnw+6), mpt6(0:ipx*mpnw+6)

call mpinitwds (mpt1, ipx*mpnw+6-5)
call mpinitwds (mpt2, ipx*mpnw+6-5)
call mpinitwds (mpt3, ipx*mpnw+6-5)
call mpinitwds (mpt4, ipx*mpnw+6-5)
call mpinitwds (mpt5, ipx*mpnw+6-5)
call mpinitwds (mpt6, ipx*mpnw+6-5)
call mpinitwds (t1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t2(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t3(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t4(0:ipx*mpnw+6), ipx*mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)
mpnw2 = min (ipx*mpnw + 1, mpnwx)
nna = abs (nn)

! ic1 = mpspacer (y)
ic1 = mpspacer (y(0:))
if (mpnw < mpnwm .or. ic1 < mpnw + 6) then
  write (6, 1)
1 format ('*** POLYLOGNEG: Uninitialized or inadequately sized array')
!   call mpabrt (590)
  call mpabrt (590)
endif

!mp prec="mpnw1"

! d1 = dpreal (arr(1))
call mpmdc (arr(0:na+5,1), mpd1, mpi1, mpnw1)
mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
d1 = mpd1
! d2 = dpreal (arr(nna))
call mpmdc (arr(0:na+5,nna), mpd1, mpi1, mpnw1)
mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
d2 = mpd1

if (d1 /= 1.e0_mpdk .or. d2 /= 1.e0_mpdk) then
  write (mpldb, 2)
2 format ('*** POLYLOGNEG: Uninitialized or inadequately sized arrays'/ &
  'Call polylogini or polylog_ini to initialize array. See documentation.')
!   call mpabrt (591)
  call mpabrt (591)
endif

!   Use higher precision, based on the largest (middle) value in arr.

! ic1 = mpsgn (x)
mpi1 = mpsgn (x(0:))
ic1 = mpi1
if (ic1 < 0) then
  i1 = (nna + 1) / 2
!   d1 = dpreal (arr(i1))
  call mpmdc (arr(0:na+5,i1), mpd1, mpi1, mpnw1)
  mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
  d1 = mpd1
  n1 = log (d1) / al2
  mpnw2 = min (mpnw + (n1 + 1) / mpnbt + 1, ipx*mpnw + 1, mpnwx)
  if (mpnw2 > ipx*mpnw+1) then
    write (6, 3) mpnw2
3   format ('*** POLYLOGNEG: Inadequate working precision:',i6)
!     call mpabrt (592)
    call mpabrt (592)
  endif
endif

!mp prec="mpnw2"
! t1 = x
call mpeq (x(0:), t1(0:ipx*mpnw+6), mpnw2)
! t2 = t1
call mpeq (t1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)

do k = 2, nna
!   t1 = x * t1
  call mpmul (x(0:), t1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw2)
!   t3 = arr(k) * t1
  call mpmul (arr(0:na+5,k), t1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
!   t2 = t2 + t3
  call mpadd (t2(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
enddo

! t3 = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
! t4 = t3 - x
call mpsub (t3(0:ipx*mpnw+6), x(0:), mpt1(0:ipx*mpnw+6), mpnw2)
call mpeq (mpt1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpnw2)
i1 = nna + 1
! t3 = t4 ** i1
call mpnpwr (t4(0:ipx*mpnw+6), i1, mpt1(0:ipx*mpnw+6), mpnw2)
call mpeq (mpt1(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpnw2)
! t4 = t2 / t3
call mpdiv (t2(0:ipx*mpnw+6), t3(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
call mpeq (mpt1(0:ipx*mpnw+6), t4(0:ipx*mpnw+6), mpnw2)
!mp prec="mpnw"
! y = t4
call mpeq (t4(0:ipx*mpnw+6), y(0:), mpnw)

return
end subroutine mppolylogneg

subroutine mppolylogpos (nn, x, y, mpnw)

!   This returns polylog (nn, x) for the case nn >= 0.

implicit none
integer, intent(in):: mpnw
integer, intent(in):: nn
!mp dim1="0:"
! type (mp_real), intent(in):: x
integer (mpiknd), intent(in):: x(0:)
! type (mp_real), intent(out):: y
integer (mpiknd), intent(out):: y(0:)
integer, parameter:: itrmax = 1000000, ipx = 3
integer ic1, ic2, k
real (mpdknd) d1
!mp dim1="0:mpnwm+5"
! type (mp_real) eps, tc1, tc2
integer (mpiknd) eps(0:mpnwm+5), tc1(0:mpnwm+5), tc2(0:mpnwm+5)
!mp dim1="0:mpnw+6"
! type (mp_real) t1, t2, t3, t4, t5
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), t5(0:mpnw+6)
integer mpnw1, mpnw2

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpt3(0:mpnw+6)
integer (mpiknd) mpt4(0:mpnw+6), mpt5(0:mpnw+6), mpt6(0:mpnw+6)

call mpinitwds (mpt1, mpnw+6-5)
call mpinitwds (mpt2, mpnw+6-5)
call mpinitwds (mpt3, mpnw+6-5)
call mpinitwds (mpt4, mpnw+6-5)
call mpinitwds (mpt5, mpnw+6-5)
call mpinitwds (mpt6, mpnw+6-5)
call mpinitwds (eps(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc1(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc2(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (t1(0:mpnw+6), mpnw+6-5)
call mpinitwds (t2(0:mpnw+6), mpnw+6-5)
call mpinitwds (t3(0:mpnw+6), mpnw+6-5)
call mpinitwds (t4(0:mpnw+6), mpnw+6-5)
call mpinitwds (t5(0:mpnw+6), mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)
mpnw2 = min (ipx*mpnw + 1, mpnwx)

! ic1 = mpspacer (y)
ic1 = mpspacer (y(0:))
if (mpnw < mpnwm .or. ic1 < mpnw + 6) then
  write (6, 1)
1 format ('*** polylogpos: Uninitialized or inadequately sized array')
!   call mpabrt (593)
  call mpabrt (593)
endif

!mp prec="mpnwm"

if (nn < 0) then
  write (mpldb, 2)
2 format ('*** POLYLOGPOS: N is less than zero.'/ &
  'For negative n, call POLYLOGNEG or POLYLOG_NEG. See documentation.')
!   call mpabrt (594)
  call mpabrt (594)
endif

! tc2 = mprealdp (2.e0_mpdk)
call mprealdp (2.e0_mpdk, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), tc2(0:mpnwm+5), mpnwm)
ic2 = -mpnw*mpnbt
! eps = tc2 ** ic2
call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), eps(0:mpnwm+5), mpnwm)

!mp prec="mpnw1"

! t1 = abs (x)
call mpabs (x(0:), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
! t2 = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! t3 = t1 - t2
call mpsub (t1(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
! ic1 = mpsgn (t3)
mpi1 = mpsgn (t3(0:mpnw+6))
ic1 = mpi1

if (ic1 >= 0) then
  write (mpldb, 3)
3 format ('*** POLYLOGPOS: Absolute value of argument must be less than one.')
!   call mpabrt (595)
  call mpabrt (595)
endif

if (nn == 0) then
!   t1 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   t2 = t1 - x
  call mpsub (t1(0:mpnw+6), x(0:), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   y = x / t2
  call mpdiv (x(0:), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), y(0:), mpnw1)
else
!   t1 = x
  call mpeq (x(0:), t1(0:mpnw+6), mpnw1)
!   t2 = x
  call mpeq (x(0:), t2(0:mpnw+6), mpnw1)

  do k = 2, itrmax
!     t2 = x * t2
    call mpmul (x(0:), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
    d1 = real (k, mpdk)
!     t3 = mprealdp (d1)
    call mprealdp (d1, mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!     t4 = t3 ** nn
    call mpnpwr (t3(0:mpnw+6), nn, mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!     t3 = t2 / t4
    call mpdiv (t2(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!     t1 = t1 + t3
    call mpadd (t1(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!mp prec="mpnwm"
!     tc1 = abs (t3) - eps * abs (t1)
    call mpabs (t3(0:mpnw+6), mpt1(0:mpnw+6), mpnwm)
    call mpabs (t1(0:mpnw+6), mpt2(0:mpnw+6), mpnwm)
    call mpmul (eps(0:mpnwm+5), mpt2(0:mpnw+6), mpt3(0:mpnw+6), mpnwm)
    call mpsub (mpt1(0:mpnw+6), mpt3(0:mpnw+6), mpt4(0:mpnw+6), mpnwm)
    call mpeq (mpt4(0:mpnw+6), tc1(0:mpnwm+5), mpnwm)
!     ic1 = mpsgn (tc1)
    mpi1 = mpsgn (tc1(0:mpnwm+5))
    ic1 = mpi1
    if (ic1 < 0) goto 100
  enddo

  write (mpldb, 4)
4 format ('*** POLYLOGPOS: Loop end error')
!   call mpabrt (596)
  call mpabrt (596)

100 continue

!mp prec="mpnw"
!   y = t1
  call mpeq (t1(0:mpnw+6), y(0:), mpnw)
endif

return
end subroutine mppolylogpos

subroutine mpstruvehn (nu, ss, zz, mpnw)

!   This returns the StruveH function with integer arg NU and MPFR argument SS.

implicit none
integer, intent(in):: mpnw
integer, intent(in):: nu
!mp dim1="0:"
! type (mp_real), intent(in):: ss
integer (mpiknd), intent(in):: ss(0:)
! type (mp_real), intent(out):: zz
integer (mpiknd), intent(out):: zz(0:)
integer, parameter:: itrmax = 1000000, ipx = 3
real (mpdknd), parameter:: dmax = 1000.e0_mpdk
integer ic1, ic2, k, n1
real (mpdknd) d1, d2
!mp dim1="0:mpnwm+5"
! type (mp_real) eps, tc1, tc2
integer (mpiknd) eps(0:mpnwm+5), tc1(0:mpnwm+5), tc2(0:mpnwm+5)
!mp dim1="0:ipx*mpnw+6"
! type (mp_real) sum, td1, td2, tn1
integer (mpiknd) sum(0:ipx*mpnw+6), td1(0:ipx*mpnw+6), td2(0:ipx*mpnw+6), tn1(0:ipx*mpnw+6)
! type (mp_real) tnm1, t1, t2
integer (mpiknd) tnm1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6)
integer mpnw1, mpnw2

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6)
integer (mpiknd) mpt4(0:ipx*mpnw+6), mpt5(0:ipx*mpnw+6), mpt6(0:ipx*mpnw+6)

call mpinitwds (mpt1, ipx*mpnw+6-5)
call mpinitwds (mpt2, ipx*mpnw+6-5)
call mpinitwds (mpt3, ipx*mpnw+6-5)
call mpinitwds (mpt4, ipx*mpnw+6-5)
call mpinitwds (mpt5, ipx*mpnw+6-5)
call mpinitwds (mpt6, ipx*mpnw+6-5)
call mpinitwds (eps(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc1(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc2(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (sum(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (td1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (td2(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (tn1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (tnm1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t1(0:ipx*mpnw+6), ipx*mpnw+6-5)
call mpinitwds (t2(0:ipx*mpnw+6), ipx*mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)
mpnw2 = min (ipx * mpnw + 1, mpnwx)

! ic1 = mpspacer (zz)
ic1 = mpspacer (zz(0:))
if (mpnw < mpnwm .or. ic1 < mpnw + 6) then
  write (6, 1)
1 format ('*** STRUVEHN: Uninitialized or inadequately sized array')
!   call mpabrt (597)
  call mpabrt (597)
endif

! ic1 = mpwprecr (mppicon)
ic1 = mpwprecr (mppicon(0:))
if (ic1 < mpnw2) then
  write (6, 2) mpnw2
2 format ('*** STRUVEHN: Pi must be initialized to',i6,' words by calling mpinit.')
!   call mpabrt (598)
  call mpabrt (598)
endif

!mp prec="mpnw1"

if (nu < 0) then
  write (mpldb, 3)
3 format ('*** STRUVEHN: Integer argument < 0')
!   call mpabrt (599)
  call mpabrt (599)
endif

! t1 = abs (ss)
call mpabs (ss(0:), mpt1(0:ipx*mpnw+6), mpnw1)
call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw1)
! d1 = dpreal (t1)
call mpmdc (t1(0:ipx*mpnw+6), mpd1, mpi1, mpnw1)
mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
d1 = mpd1
if (d1 > dmax) then
  write (mpldb, 4) dmax
4 format ('*** STRUVEHN: Argument too large >',f8.2)
!   call mpabrt (600)
  call mpabrt (600)
endif

mpnw2 = min (mpnw + nint (d1 * mpnw / dmax) + 1, mpnwx)
if (mpnw2 > ipx*mpnw+1) then
  write (6, 5) mpnw2
5 format ('*** STRUVEHNN: Inadequate working precision:',i6)
!   call mpabrt (601)
  call mpabrt (601)
endif

!mp prec="mpnwm"
! tc2 = mprealdp (2.e0_mpdk)
call mprealdp (2.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnwm)
call mpeq (mpt1(0:ipx*mpnw+6), tc2(0:mpnwm+5), mpnwm)
ic2 = -mpnw2*mpnbt
! eps = tc2 ** ic2
call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:ipx*mpnw+6), mpnwm)
call mpeq (mpt1(0:ipx*mpnw+6), eps(0:mpnwm+5), mpnwm)

!mp prec="mpnw2"
! tn1 = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
call mpeq (mpt1(0:ipx*mpnw+6), tn1(0:ipx*mpnw+6), mpnw2)
! t1 = ss ** 2
call mpnpwr (ss(0:), 2, mpt1(0:ipx*mpnw+6), mpnw2)
call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw2)
! tnm1 = - t1 * 0.25e0_mpdk
call mpmuld (t1(0:ipx*mpnw+6), 0.25e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
call mpneg (mpt1(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnw2)
call mpeq (mpt2(0:ipx*mpnw+6), tnm1(0:ipx*mpnw+6), mpnw2)
! td1 = sqrt (mppicon) * 0.5e0_mpdk
call mpsqrt (mppicon(0:), mpt1(0:ipx*mpnw+6), mpnw2)
call mpmuld (mpt1(0:ipx*mpnw+6), 0.5e0_mpdk, mpt2(0:ipx*mpnw+6), mpnw2)
call mpeq (mpt2(0:ipx*mpnw+6), td1(0:ipx*mpnw+6), mpnw2)
! td2 = td1
call mpeq (td1(0:ipx*mpnw+6), td2(0:ipx*mpnw+6), mpnw2)

do k = 1, nu
  d1 = real (k, mpdk) + 0.5e0_mpdk
!   td2 = td2 * d1
  call mpmuld (td2(0:ipx*mpnw+6), d1, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), td2(0:ipx*mpnw+6), mpnw2)
enddo

! sum = tn1 / td1 / td2
call mpdiv (tn1(0:ipx*mpnw+6), td1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
call mpdiv (mpt1(0:ipx*mpnw+6), td2(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnw2)
call mpeq (mpt2(0:ipx*mpnw+6), sum(0:ipx*mpnw+6), mpnw2)

do k = 1, itrmax
!   tn1 = tn1 * tnm1
  call mpmul (tn1(0:ipx*mpnw+6), tnm1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), tn1(0:ipx*mpnw+6), mpnw2)
  d1 = real (k, mpdk) + 0.5e0_mpdk
!   td1 = td1 * d1
  call mpmuld (td1(0:ipx*mpnw+6), d1, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), td1(0:ipx*mpnw+6), mpnw2)
  d2 = real (nu+k, mpdk) + 0.5e0_mpdk
!   td2 = td2 * d2
  call mpmuld (td2(0:ipx*mpnw+6), d2, mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), td2(0:ipx*mpnw+6), mpnw2)
!   t2 = td1 * td2
  call mpmul (td1(0:ipx*mpnw+6), td2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
!   t1 = tn1 / t2
  call mpdiv (tn1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw2)
!   sum = sum + t1
  call mpadd (sum(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
  call mpeq (mpt1(0:ipx*mpnw+6), sum(0:ipx*mpnw+6), mpnw2)
!mp prec="mpnwm"
!   tc1 = abs (t1) - eps * abs (sum)
  call mpabs (t1(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnwm)
  call mpabs (sum(0:ipx*mpnw+6), mpt2(0:ipx*mpnw+6), mpnwm)
  call mpmul (eps(0:mpnwm+5), mpt2(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpnwm)
  call mpsub (mpt1(0:ipx*mpnw+6), mpt3(0:ipx*mpnw+6), mpt4(0:ipx*mpnw+6), mpnwm)
  call mpeq (mpt4(0:ipx*mpnw+6), tc1(0:mpnwm+5), mpnwm)
!   ic1 = mpsgn (tc1)
  mpi1 = mpsgn (tc1(0:mpnwm+5))
  ic1 = mpi1
  if (ic1 < 0) goto 100
enddo

write (mpldb, 6)
6 format ('*** STRUVEHN: Loop end error')
! call mpabrt (602)
call mpabrt (602)

100 continue

!mp prec="mpnw2"
! t1 = ss * 0.5e0_mpdk
call mpmuld (ss(0:), 0.5e0_mpdk, mpt1(0:ipx*mpnw+6), mpnw2)
call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw2)
ic1 = nu + 1
! t2 = t1 ** ic1
call mpnpwr (t1(0:ipx*mpnw+6), ic1, mpt1(0:ipx*mpnw+6), mpnw2)
call mpeq (mpt1(0:ipx*mpnw+6), t2(0:ipx*mpnw+6), mpnw2)
! t1 = t2 * sum
call mpmul (t2(0:ipx*mpnw+6), sum(0:ipx*mpnw+6), mpt1(0:ipx*mpnw+6), mpnw2)
call mpeq (mpt1(0:ipx*mpnw+6), t1(0:ipx*mpnw+6), mpnw2)
!mp prec="mpnw"
! zz = t1
call mpeq (t1(0:ipx*mpnw+6), zz(0:), mpnw)

return
end subroutine mpstruvehn

subroutine mpzetar (ss, zz, mpnw)

!   This returns the zeta function of an MPR argument SS using an algorithm
!   due to Peter Borwein.

implicit none
integer, intent(in):: mpnw
!mp dim1="0:"
! type (mp_real), intent(in):: ss
integer (mpiknd), intent(in):: ss(0:)
! type (mp_real), intent(out):: zz
integer (mpiknd), intent(out):: zz(0:)
integer, parameter:: itrmax = 1000000
real (mpdknd), parameter:: dcon = 0.31e0_mpdk
real (mpdknd) sgn
integer i, i1, i2, ic1, ic2, ic3, iss, j, n, n1, n2
real (mpdknd) d1, d2, d3
!mp dim1="0:mpnwm+5"
! type (mp_real) eps, tc1, tc2
integer (mpiknd) eps(0:mpnwm+5), tc1(0:mpnwm+5), tc2(0:mpnwm+5)
!mp dim1="0:mpnw+6"
! type (mp_real) f1, s, t1, t2, t3
integer (mpiknd) f1(0:mpnw+6), s(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6)
! type (mp_real) t4, t5, tn, tt
integer (mpiknd) t4(0:mpnw+6), t5(0:mpnw+6), tn(0:mpnw+6), tt(0:mpnw+6)
integer mpnw1

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpt3(0:mpnw+6)
integer (mpiknd) mpt4(0:mpnw+6), mpt5(0:mpnw+6), mpt6(0:mpnw+6)

call mpinitwds (mpt1, mpnw+6-5)
call mpinitwds (mpt2, mpnw+6-5)
call mpinitwds (mpt3, mpnw+6-5)
call mpinitwds (mpt4, mpnw+6-5)
call mpinitwds (mpt5, mpnw+6-5)
call mpinitwds (mpt6, mpnw+6-5)
call mpinitwds (eps(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc1(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc2(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (f1(0:mpnw+6), mpnw+6-5)
call mpinitwds (s(0:mpnw+6), mpnw+6-5)
call mpinitwds (t1(0:mpnw+6), mpnw+6-5)
call mpinitwds (t2(0:mpnw+6), mpnw+6-5)
call mpinitwds (t3(0:mpnw+6), mpnw+6-5)
call mpinitwds (t4(0:mpnw+6), mpnw+6-5)
call mpinitwds (t5(0:mpnw+6), mpnw+6-5)
call mpinitwds (tn(0:mpnw+6), mpnw+6-5)
call mpinitwds (tt(0:mpnw+6), mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)

! ic1 = mpspacer (zz)
ic1 = mpspacer (zz(0:))
if (mpnw < mpnwm .or. ic1 < mpnw + 6) then
  write (6, 1)
1 format ('*** ZETAR: Uninitialized or inadequately sized array')
!   call mpabrt (603)
  call mpabrt (603)
endif

! ic1 = mpwprecr (mppicon)
ic1 = mpwprecr (mppicon(0:))
if (ic1 < mpnw1) then
  write (6, 2) mpnw1
2 format ('*** ZETAR: Pi must be initialized to',i6,' words by calling mpinit.')
!   call mpabrt (604)
  call mpabrt (604)
endif

!mp prec="mpnwm"

! tc2 = mprealdp (2.e0_mpdk)
call mprealdp (2.e0_mpdk, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), tc2(0:mpnwm+5), mpnwm)
ic2 = -mpnw1*mpnbt
! eps = tc2 ** ic2
call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), eps(0:mpnwm+5), mpnwm)

!mp prec="mpnw1"

! f1 = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), f1(0:mpnw+6), mpnw1)
! t1 = ss - f1
call mpsub (ss(0:), f1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
! t2 = ss - anint (ss)
call mpnint (ss(0:), mpt1(0:mpnw+6), mpnw1)
call mpsub (ss(0:), mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
call mpeq (mpt2(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! ic1 = mpsgn (t1)
mpi1 = mpsgn (t1(0:mpnw+6))
ic1 = mpi1
! ic2 = mpsgn (t2)
mpi1 = mpsgn (t2(0:mpnw+6))
ic2 = mpi1
! ic3 = mpsgn (ss)
mpi1 = mpsgn (ss(0:))
ic3 = mpi1

if (ic1 == 0) then
  write (mpldb, 3)
3 format ('*** ZETAR: Argument is 1')
!   call mpabrt (605)
  call mpabrt (605)
elseif (ic2 == 0) then

!   The argument is an integer value. Call mpzetaintr instead.

!   d1 = dpreal (ss)
  call mpmdc (ss(0:), mpd1, mpi1, mpnw1)
  mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
  d1 = mpd1
  iss = nint (d1)
!   call mpzetaintr (iss, zz, mpnw)
  call mpzetaintr (iss, zz, mpnw)
  goto 210
elseif (ic3 < 0) then

!   If arg < 0, compute zeta(1-ss), and later apply Riemann's formula.

!   tt = f1 - ss
  call mpsub (f1(0:mpnw+6), ss(0:), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), tt(0:mpnw+6), mpnw1)
else
!   tt = ss
  call mpeq (ss(0:), tt(0:mpnw+6), mpnw1)
endif

!  Check if argument is large enough that computing with definition is faster.

d3 = log (2.e0_mpdk) / log (2.e0_mpdk * mpnbt * mpnw1 / 3.e0_mpdk)
d1 = mpnbt * mpnw1 * d3
! d2 = dpreal (tt)
call mpmdc (tt(0:mpnw+6), mpd1, mpi1, mpnw1)
mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
d2 = mpd1

if (d2 > d1) then

!   Evaluate the infinite series.

!   t1 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)

  do i = 2, itrmax
    d3 = real (i, mpdk)
!     t4 = mprealdp (d3)
    call mprealdp (d3, mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!     t2 = t4 ** tt
    call mppower (t4(0:mpnw+6), tt(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!     t3 = f1 / t2
    call mpdiv (f1(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!     t1 = t1 + t3
    call mpadd (t1(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!mp prec="mpnwm"
!     tc1 = abs (t3) - eps * abs (t1)
    call mpabs (t3(0:mpnw+6), mpt1(0:mpnw+6), mpnwm)
    call mpabs (t1(0:mpnw+6), mpt2(0:mpnw+6), mpnwm)
    call mpmul (eps(0:mpnwm+5), mpt2(0:mpnw+6), mpt3(0:mpnw+6), mpnwm)
    call mpsub (mpt1(0:mpnw+6), mpt3(0:mpnw+6), mpt4(0:mpnw+6), mpnwm)
    call mpeq (mpt4(0:mpnw+6), tc1(0:mpnwm+5), mpnwm)
!     ic1 = mpsgn (tc1)
    mpi1 = mpsgn (tc1(0:mpnwm+5))
    ic1 = mpi1
    if (ic1 < 0) goto 200
  enddo

  write (mpldb, 4)
4 format ('*** ZETAR: End loop error')
!   call mpabrt (606)
  call mpabrt (606)
endif

!mp prec="mpnw1"
n = nint (dcon * mpnbt * mpnw1)
! t1 = mprealdp (2.e0_mpdk)
call mprealdp (2.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
! tn = t1 ** n
call mpnpwr (t1(0:mpnw+6), n, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), tn(0:mpnw+6), mpnw1)
! t1 = - tn
call mpneg (tn(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
! t2 = mprealdp (0.e0_mpdk)
call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! s = mprealdp (0.e0_mpdk)
call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), s(0:mpnw+6), mpnw1)
sgn = 1.e0_mpdk

do j = 0, 2 * n - 1
  i1 = j + 1
  d3 = real (i1, mpdk)
!   t4 = mprealdp (d3)
  call mprealdp (d3, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!   t3 = t4 ** tt
  call mppower (t4(0:mpnw+6), tt(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!   t4 = t1 / t3
  call mpdiv (t1(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!   t5 = t4 * sgn
  call mpmuld (t4(0:mpnw+6), sgn, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
!   t4 = s + t5
  call mpadd (s(0:mpnw+6), t5(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!   s = t4
  call mpeq (t4(0:mpnw+6), s(0:mpnw+6), mpnw1)
  sgn = - sgn

  if (j < n - 1) then
!     t2 = mprealdp (0.e0_mpdk)
    call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
  elseif (j == n - 1) then
!     t2 = mprealdp (1.e0_mpdk)
    call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
  else
    i1 = 2*n - j
    i2 = j + 1 - n
!     t3 = t2 * real (i1, mpdk)
    call mpmuld (t2(0:mpnw+6), real (i1, mpdk), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!     t2 = t3 / real (i2, mpdk)
    call mpdivd (t3(0:mpnw+6), real (i2, mpdk), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
  endif
!   t1 = t1 + t2
  call mpadd (t1(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
enddo

! t3 = 1.e0_mpdk - tt
call mpdmc (1.e0_mpdk, 0, mpt2(0:mpnw+6), mpnw1)
call mpsub (mpt2(0:mpnw+6), tt(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
! t2 = mprealdp (2.e0_mpdk)
call mprealdp (2.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! t4 = t2 ** t3
call mppower (t2(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
! t2 = f1 - t4
call mpsub (f1(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! t3 = tn * t2
call mpmul (tn(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
! t1 = - s / t3
call mpdiv (s(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpneg (mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
call mpeq (mpt2(0:mpnw+6), t1(0:mpnw+6), mpnw1)

!   If original argument was negative, apply Riemann's formula.

! ic1 = mpsgn (ss)
mpi1 = mpsgn (ss(0:))
ic1 = mpi1
if (ic1 < 0) then
!   call mpgammar (tt, t3, mpnw1)
  call mpgammar (tt, t3, mpnw1)
!   t2 = t1 * t3
  call mpmul (t1(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   t1 = mppicon * tt
  call mpmul (mppicon(0:), tt(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   t3 = t1 * 0.5e0_mpdk
  call mpmuld (t1(0:mpnw+6), 0.5e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!   t4 = cos (t3)
  call mpcssnr (t3(0:mpnw+6), mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!   t5 = sin (t3)
  call mpcssnr (t3(0:mpnw+6), mpt2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
!   t1 = t2 * t4
  call mpmul (t2(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   t2 = mppicon * 2.e0_mpdk
  call mpmuld (mppicon(0:), 2.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   t3 = t2 ** tt
  call mppower (t2(0:mpnw+6), tt(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!   t1 = t1 / t3 * 2.e0_mpdk
  call mpdiv (t1(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpmuld (mpt1(0:mpnw+6), 2.e0_mpdk, mpt2(0:mpnw+6), mpnw1)
  call mpeq (mpt2(0:mpnw+6), t1(0:mpnw+6), mpnw1)
endif

200 continue

!mp prec="mpnw"
! zz = t1
call mpeq (t1(0:mpnw+6), zz(0:), mpnw)

210 continue

return
end subroutine mpzetar

subroutine mpzetaintr (iss, zz, mpnw)

!   This returns the zeta function of the integer argument ISS using an algorithm
!   due to Peter Borwein.

implicit none
integer, intent(in):: mpnw
integer, intent(in):: iss
!mp dim1="0:"
! type (mp_real), intent(out):: zz
integer (mpiknd), intent(out):: zz(0:)
integer, parameter:: itrmax = 1000000
real (mpdknd), parameter:: dcon = 0.31e0_mpdk
integer i, i1, i2, ic1, ic2, j, n, n1, itt
real (mpdknd) d1, d2, d3, sgn
!mp dim1="0:mpnwm+5"
! type (mp_real) eps, tc1, tc2
integer (mpiknd) eps(0:mpnwm+5), tc1(0:mpnwm+5), tc2(0:mpnwm+5)
!mp dim1="0:mpnw+6"
! type (mp_real) f1, s, t1, t2, t3
integer (mpiknd) f1(0:mpnw+6), s(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6)
! type (mp_real) t4, t5, tn
integer (mpiknd) t4(0:mpnw+6), t5(0:mpnw+6), tn(0:mpnw+6)
integer mpnw1, mpnw2

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpt3(0:mpnw+6)
integer (mpiknd) mpt4(0:mpnw+6), mpt5(0:mpnw+6), mpt6(0:mpnw+6)

call mpinitwds (mpt1, mpnw+6-5)
call mpinitwds (mpt2, mpnw+6-5)
call mpinitwds (mpt3, mpnw+6-5)
call mpinitwds (mpt4, mpnw+6-5)
call mpinitwds (mpt5, mpnw+6-5)
call mpinitwds (mpt6, mpnw+6-5)
call mpinitwds (eps(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc1(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc2(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (f1(0:mpnw+6), mpnw+6-5)
call mpinitwds (s(0:mpnw+6), mpnw+6-5)
call mpinitwds (t1(0:mpnw+6), mpnw+6-5)
call mpinitwds (t2(0:mpnw+6), mpnw+6-5)
call mpinitwds (t3(0:mpnw+6), mpnw+6-5)
call mpinitwds (t4(0:mpnw+6), mpnw+6-5)
call mpinitwds (t5(0:mpnw+6), mpnw+6-5)
call mpinitwds (tn(0:mpnw+6), mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)

! ic1 = mpspacer (zz)
ic1 = mpspacer (zz(0:))
if (mpnw < mpnwm .or. ic1 < mpnw + 6) then
  write (6, 1)
1 format ('*** ZETAINTR: Uninitialized or inadequately sized array')
!   call mpabrt (607)
  call mpabrt (607)
endif

! ic1 = mpwprecr (mppicon)
ic1 = mpwprecr (mppicon(0:))
if (ic1 < mpnw1) then
  write (6, 2) mpnw1
2 format ('*** ZETAINTR: Pi must be initialized to',i6,' words by calling mpinit.')
!   call mpabrt (608)
  call mpabrt (608)
endif

!mp prec="mpnwm"

! tc2 = mprealdp (2.e0_mpdk)
call mprealdp (2.e0_mpdk, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), tc2(0:mpnwm+5), mpnwm)
ic2 = -mpnw1*mpnbt
! eps = tc2 ** ic2
call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), eps(0:mpnwm+5), mpnwm)

!mp prec="mpnw1"

! f1 = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), f1(0:mpnw+6), mpnw1)

if (iss == 1) then
  write (mpldb, 3)
3 format ('*** ZETAINTR: Argument is 1')
!   call mpabrt (609)
  call mpabrt (609)
elseif (iss == 0) then

!   Argument is zero -- result is -1/2.

!   t1 = mprealdp (-0.5e0_mpdk)
  call mprealdp (-0.5e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
  goto 200
elseif (iss < 0) then

!   If argument is a negative even integer, the result is zero.

  if (mod (iss, 2) == 0) then
!     t1 = mprealdp (0.e0_mpdk)
    call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
    goto 200
  endif

!   Otherwise if arg < 0, compute zeta(1-is), and later apply Riemann's formula.

  itt = 1 - iss
else
  itt = iss
endif

!  Check if argument is large enough that computing with definition is faster.

d3 = log (2.e0_mpdk) / log (2.e0_mpdk * mpnbt * mpnw1 / 3.e0_mpdk)
d1 = mpnbt * mpnw1 * d3

if (itt > d1) then

!   Evaluate the infinite series.

!   t1 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)

  do i = 2, itrmax
    d3 = real (i, mpdk)
!     t4 = mprealdp (d3)
    call mprealdp (d3, mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!     t2 = t4 ** itt
    call mpnpwr (t4(0:mpnw+6), itt, mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!     t3 = f1 / t2
    call mpdiv (f1(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!     t1 = t1 + t3
    call mpadd (t1(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!mp prec="mpnwm"
!     tc1 = abs (t3) - eps * abs (t1)
    call mpabs (t3(0:mpnw+6), mpt1(0:mpnw+6), mpnwm)
    call mpabs (t1(0:mpnw+6), mpt2(0:mpnw+6), mpnwm)
    call mpmul (eps(0:mpnwm+5), mpt2(0:mpnw+6), mpt3(0:mpnw+6), mpnwm)
    call mpsub (mpt1(0:mpnw+6), mpt3(0:mpnw+6), mpt4(0:mpnw+6), mpnwm)
    call mpeq (mpt4(0:mpnw+6), tc1(0:mpnwm+5), mpnwm)
!     ic1 = mpsgn (tc1)
    mpi1 = mpsgn (tc1(0:mpnwm+5))
    ic1 = mpi1
    if (ic1 < 0) goto 200
  enddo

  write (mpldb, 4)
4 format ('*** ZETAINTR: End loop error')
!   call mpabrt (610)
  call mpabrt (610)
endif

!mp prec="mpnw1"
n = nint (dcon * mpnbt * mpnw1)
! t1 = mprealdp (2.e0_mpdk)
call mprealdp (2.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
! tn = t1 ** n
call mpnpwr (t1(0:mpnw+6), n, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), tn(0:mpnw+6), mpnw1)
! t1 = - tn
call mpneg (tn(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
! t2 = mprealdp (0.e0_mpdk)
call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! s = mprealdp (0.e0_mpdk)
call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), s(0:mpnw+6), mpnw1)
sgn = 1.e0_mpdk

do j = 0, 2 * n - 1
  i1 = j + 1
  d3 = real (i1, mpdk)
!   t4 = mprealdp (d3)
  call mprealdp (d3, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!   t3 = t4 ** itt
  call mpnpwr (t4(0:mpnw+6), itt, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!   t4 = t1 / t3
  call mpdiv (t1(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!   t5 = t4 * sgn
  call mpmuld (t4(0:mpnw+6), sgn, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
!   s = s + t5
  call mpadd (s(0:mpnw+6), t5(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), s(0:mpnw+6), mpnw1)
  sgn = - sgn

  if (j < n - 1) then
!     t2 = mprealdp (0.e0_mpdk)
    call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
  elseif (j == n - 1) then
!     t2 = mprealdp (1.e0_mpdk)
    call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
  else
    i1 = 2*n - j
    i2 = j + 1 - n
!     t3 = t2 * real (i1, mpdk)
    call mpmuld (t2(0:mpnw+6), real (i1, mpdk), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!     t2 = t3 / real (i2, mpdk)
    call mpdivd (t3(0:mpnw+6), real (i2, mpdk), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
  endif

!   t1 = t1 + t2
  call mpadd (t1(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
enddo

! t2 = mprealdp (2.e0_mpdk)
call mprealdp (2.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
ic1 = 1 - itt
! t4 = t2 ** ic1
call mpnpwr (t2(0:mpnw+6), ic1, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
! t2 = f1 - t4
call mpsub (f1(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! t3 = tn * t2
call mpmul (tn(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
! t1 = - s / t3
call mpdiv (s(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpneg (mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
call mpeq (mpt2(0:mpnw+6), t1(0:mpnw+6), mpnw1)

!   If original argument was negative, apply Riemann's formula.

if (iss < 0) then
!   t3 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)

  do i = 1, itt - 1
!     t3 = t3 * real (i, mpdk)
    call mpmuld (t3(0:mpnw+6), real (i, mpdk), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
  enddo

!   t2 = t1 * t3
  call mpmul (t1(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   t1 = mppicon * real (itt, mpdk)
  call mpmuld (mppicon(0:), real (itt, mpdk), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   t3 = t1 * 0.5e0_mpdk
  call mpmuld (t1(0:mpnw+6), 0.5e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!   t4 = cos (t3)
  call mpcssnr (t3(0:mpnw+6), mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!   t5 = sin (t3)
  call mpcssnr (t3(0:mpnw+6), mpt2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
!   t1 = t2 * t4
  call mpmul (t2(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   t2 = mppicon * 2.e0_mpdk
  call mpmuld (mppicon(0:), 2.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   t3 = t2 ** itt
  call mpnpwr (t2(0:mpnw+6), itt, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!   t1 = t1 / t3 * 2.e0_mpdk
  call mpdiv (t1(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpmuld (mpt1(0:mpnw+6), 2.e0_mpdk, mpt2(0:mpnw+6), mpnw1)
  call mpeq (mpt2(0:mpnw+6), t1(0:mpnw+6), mpnw1)
endif

200 continue

!mp prec="mpnw"
! zz = t1
call mpeq (t1(0:mpnw+6), zz(0:), mpnw)

return
end subroutine mpzetaintr

subroutine mpzetabe (nb1, nb2, berne, s, z, mpnw)

!  This evaluates the Riemann zeta function, using the combination of
!  the definition formula (for large s), and an Euler-Maclaurin scheme
!  (see formula 25.2.9 of the DLMF). The array berne contains precomputed
!  even Bernoulli numbers (see MPBERNER above). Its dimensions must be as
!  shown below. NB2 must be greater than 1.4 x precision in decimal digits.

implicit none
integer, intent(in):: nb1, mpnw
integer, intent(in):: nb2
!mp dim1="0:nb1+5"
! type (mp_real), intent(in):: berne(1:nb2)
integer (mpiknd), intent(in):: berne(0:nb1+5,1:nb2)
!mp dim1="0:"
! type (mp_real), intent(in):: s
integer (mpiknd), intent(in):: s(0:)
! type (mp_real), intent(out):: z
integer (mpiknd), intent(out):: z(0:)
integer, parameter:: itrmax = 1000000
real (mpdknd), parameter:: dber = 0.45e0_mpdk, dfrac = 0.18e0_mpdk
integer i, ic1, ic2, ic3, i1, i2, k, n1, n2, nn
real (mpdknd) d1, d2, d3, d4
!mp dim1="0:mpnwm+5"
! type (mp_real) eps, tc1, tc2
integer (mpiknd) eps(0:mpnwm+5), tc1(0:mpnwm+5), tc2(0:mpnwm+5)
!mp dim1="0:mpnw+6"
! type (mp_real) t0, t1, t2, t3, t4
integer (mpiknd) t0(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6)
! type (mp_real) t5, t6, t7, t8, t9
integer (mpiknd) t5(0:mpnw+6), t6(0:mpnw+6), t7(0:mpnw+6), t8(0:mpnw+6), t9(0:mpnw+6)
! type (mp_real) tt, f1
integer (mpiknd) tt(0:mpnw+6), f1(0:mpnw+6)
integer mpnw1

!mp enddec
integer mpi1, mpi2, mpi3, mpi4, mpi5, mpi6
real (mpdknd) mpd1, mpd2, mpd3, mpd4, mpd5, mpd6
integer (mpiknd) mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpt3(0:mpnw+6)
integer (mpiknd) mpt4(0:mpnw+6), mpt5(0:mpnw+6), mpt6(0:mpnw+6)

call mpinitwds (mpt1, mpnw+6-5)
call mpinitwds (mpt2, mpnw+6-5)
call mpinitwds (mpt3, mpnw+6-5)
call mpinitwds (mpt4, mpnw+6-5)
call mpinitwds (mpt5, mpnw+6-5)
call mpinitwds (mpt6, mpnw+6-5)
call mpinitwds (eps(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc1(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (tc2(0:mpnwm+5), mpnwm+5-5)
call mpinitwds (t0(0:mpnw+6), mpnw+6-5)
call mpinitwds (t1(0:mpnw+6), mpnw+6-5)
call mpinitwds (t2(0:mpnw+6), mpnw+6-5)
call mpinitwds (t3(0:mpnw+6), mpnw+6-5)
call mpinitwds (t4(0:mpnw+6), mpnw+6-5)
call mpinitwds (t5(0:mpnw+6), mpnw+6-5)
call mpinitwds (t6(0:mpnw+6), mpnw+6-5)
call mpinitwds (t7(0:mpnw+6), mpnw+6-5)
call mpinitwds (t8(0:mpnw+6), mpnw+6-5)
call mpinitwds (t9(0:mpnw+6), mpnw+6-5)
call mpinitwds (tt(0:mpnw+6), mpnw+6-5)
call mpinitwds (f1(0:mpnw+6), mpnw+6-5)

mpnw1 = min (mpnw + 1, mpnwx)

! ic1 = mpspacer (z)
ic1 = mpspacer (z(0:))
if (mpnw < mpnwm .or. ic1 < mpnw + 6) then
  write (6, 1)
1 format ('*** ZETABE: Uninitialized or inadequately sized array')
!   call mpabrt (611)
  call mpabrt (611)
endif

! ic1 = mpwprecr (mppicon)
ic1 = mpwprecr (mppicon(0:))
if (ic1 < mpnw1) then
  write (6, 2) mpnw1
2 format ('*** ZETABE: Pi must be initialized to',i6,' words by calling mpinit.')
!   call mpabrt (612)
  call mpabrt (612)
endif

!mp prec="mpnw1"

!   Check if berne array has been initialized.

! d1 = dpreal (berne(1))
call mpmdc (berne(0:nb1+5,1), mpd1, mpi1, mpnw1)
mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
d1 = mpd1
if (abs (d1 - 1.e0_mpdk / 6.e0_mpdk) > mprdfz .or. nb2 < int (dber * mpnbt * mpnw1)) then
  write (mpldb, 3) int (dber * mpnbt * mpnw1)
3 format ('*** ZETABE: Array of even Bernoulli coefficients must be initialized'/ &
   'with at least',i8,' entries.')
!   call mpabrt (613)
  call mpabrt (613)
endif

i = 0
k = 0
!mp prec="mpnwm"
! tc2 = mprealdp (2.e0_mpdk)
call mprealdp (2.e0_mpdk, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), tc2(0:mpnwm+5), mpnwm)
ic2 = -mpnw1*mpnbt
! eps = tc2 ** ic2
call mpnpwr (tc2(0:mpnwm+5), ic2, mpt1(0:mpnw+6), mpnwm)
call mpeq (mpt1(0:mpnw+6), eps(0:mpnwm+5), mpnwm)

!   Check if argument is 1 -- undefined.

!mp prec="mpnw1"
! t0 = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t0(0:mpnw+6), mpnw1)
! t1 = s - t0
call mpsub (s(0:), t0(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
! ic1 = mpsgn (t1)
mpi1 = mpsgn (t1(0:mpnw+6))
ic1 = mpi1
! ic2 = mpsgn (s)
mpi1 = mpsgn (s(0:))
ic2 = mpi1

if (ic1 == 0) then
  write (mpldb, 4)
4 format ('*** ZETABE: Argument is 1')
!   call mpabrt (614)
  call mpabrt (614)
endif

! f1 = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), f1(0:mpnw+6), mpnw1)

!   Check if argument is zero. If so, result is - 1/2.

if (ic2 == 0) then
!   t1 = mprealdp (-0.5e0_mpdk)
  call mprealdp (-0.5e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
  goto 200
endif

!   Check if argument is negative.

if (ic2 < 0) then

!   Check if argument is a negative even integer. If so, the result is zero.

!   t1 = s * 0.5e0_mpdk
  call mpmuld (s(0:), 0.5e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   t3 = t1 - anint (t1)
  call mpnint (t1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpsub (t1(0:mpnw+6), mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
  call mpeq (mpt2(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!   ic3 = mpsgn (t3)
  mpi1 = mpsgn (t3(0:mpnw+6))
  ic3 = mpi1
  if (ic3 == 0) then
!     t1 = mprealdp (0.e0_mpdk)
    call mprealdp (0.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
    goto 200
  endif

!   Otherwise compute zeta(1-s), and later apply the reflection formula.

!   tt = f1 - s
  call mpsub (f1(0:mpnw+6), s(0:), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), tt(0:mpnw+6), mpnw1)
else
!   tt = s
  call mpeq (s(0:), tt(0:mpnw+6), mpnw1)
endif

!  Check if argument is large enough that computing with definition is faster.

d1 = mplogb * mpnw1 / log (32.e0_mpdk * mpnw1)
! d2 = dpreal (tt)
call mpmdc (tt(0:mpnw+6), mpd1, mpi1, mpnw1)
mpd1 = mpd1 * 2.e0_mpdknd  ** mpi1
d2 = mpd1

if (d2 > d1) then
!   t1 = mprealdp (1.e0_mpdk)
  call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)

  do i = 2, itrmax
    d3 = real (i, mpdk)
!     t4 = mprealdp (d3)
    call mprealdp (d3, mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!     t2 = t4 ** tt
    call mppower (t4(0:mpnw+6), tt(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!     t3 = f1 / t2
    call mpdiv (f1(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!     t1 = t1 + t3
    call mpadd (t1(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
    call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!mp prec="mpnwm"
!     tc1 = abs (t3) - eps * abs (t1)
    call mpabs (t3(0:mpnw+6), mpt1(0:mpnw+6), mpnwm)
    call mpabs (t1(0:mpnw+6), mpt2(0:mpnw+6), mpnwm)
    call mpmul (eps(0:mpnwm+5), mpt2(0:mpnw+6), mpt3(0:mpnw+6), mpnwm)
    call mpsub (mpt1(0:mpnw+6), mpt3(0:mpnw+6), mpt4(0:mpnw+6), mpnwm)
    call mpeq (mpt4(0:mpnw+6), tc1(0:mpnwm+5), mpnwm)
!     ic1 = mpsgn (tc1)
    mpi1 = mpsgn (tc1(0:mpnwm+5))
    ic1 = mpi1
    if (ic1 < 0) goto 200
  enddo

  write (mpldb, 5)
5 format ('*** ZETABE: End loop error 1')
!   call mpabrt (615)
  call mpabrt (615)
endif

!mp prec="mpnw1"
! t0 = mprealdp (1.e0_mpdk)
call mprealdp (1.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t0(0:mpnw+6), mpnw1)
nn = dfrac * mpnbt * mpnw1

do k = 2, nn
  d3 = real (k, mpdk)
!   t2 = mprealdp (d3)
  call mprealdp (d3, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   t1 = t2 ** tt
  call mppower (t2(0:mpnw+6), tt(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   t2 = f1 / t1
  call mpdiv (f1(0:mpnw+6), t1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   t0 = t0 + t2
  call mpadd (t0(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t0(0:mpnw+6), mpnw1)
enddo

d3 = real (nn, mpdk)
! t2 = mprealdp (d3)
call mprealdp (d3, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! t3 = tt - f1
call mpsub (tt(0:mpnw+6), f1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
! t4 = t1 * t3
call mpmul (t1(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
! t3 = t2 / t4
call mpdiv (t2(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
! t2 = t0 + t3
call mpadd (t0(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! t3 = mprealdp (0.5e0_mpdk)
call mprealdp (0.5e0_mpdk, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
! t4 = t3 / t1
call mpdiv (t3(0:mpnw+6), t1(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
! t0 = t2 - t4
call mpsub (t2(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t0(0:mpnw+6), mpnw1)

! t3 = tt
call mpeq (tt(0:mpnw+6), t3(0:mpnw+6), mpnw1)
d1 = 12.e0_mpdk * real (nn, mpdk)
! t4 = t1 * d1
call mpmuld (t1(0:mpnw+6), d1, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
! t2 = t3 / t4
call mpdiv (t3(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
! t5 = t1 * real (nn, mpdk)
call mpmuld (t1(0:mpnw+6), real (nn, mpdk), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
d3 = real (nn, mpdk)
! t6 = mprealdp (d3)
call mprealdp (d3, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
! t9 = t6 ** 2
call mpnpwr (t6(0:mpnw+6), 2, mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t9(0:mpnw+6), mpnw1)

do k = 2, min (nb2, itrmax)
  i1 = 2*k - 2
  i2 = 2*k - 3
  d3 = real (i1, mpdk)
!   t4 = mprealdp (d3)
  call mprealdp (d3, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!   t6 = tt + t4
  call mpadd (tt(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
  d4 = real (i2, mpdk)
!   t7 = mprealdp (d4)
  call mprealdp (d4, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t7(0:mpnw+6), mpnw1)
!   t8 = tt + t7
  call mpadd (tt(0:mpnw+6), t7(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t8(0:mpnw+6), mpnw1)
!   t7 = t6 * t8
  call mpmul (t6(0:mpnw+6), t8(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t7(0:mpnw+6), mpnw1)
!   t4 = t3 * t7
  call mpmul (t3(0:mpnw+6), t7(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
  i1 = 2*k - 1
  i2 = 2*k - 2
  d3 = real (i1, mpdk)
!   t6 = mprealdp (d3)
  call mprealdp (d3, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
  d4 = real (i2, mpdk)
!   t7 = mprealdp (d4)
  call mprealdp (d4, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t7(0:mpnw+6), mpnw1)
!   t8 = t6 * t7
  call mpmul (t6(0:mpnw+6), t7(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t8(0:mpnw+6), mpnw1)
!   t3 = t4 / t8
  call mpdiv (t4(0:mpnw+6), t8(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!   t5 = t5 * t9
  call mpmul (t5(0:mpnw+6), t9(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
!   t4 = t3 * berne(k)
  call mpmul (t3(0:mpnw+6), berne(0:nb1+5,k), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
  i1 = 2*k
!   t6 = t5 * real (i1, mpdk)
  call mpmuld (t5(0:mpnw+6), real (i1, mpdk), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t6(0:mpnw+6), mpnw1)
!   t7 = t4 / t6
  call mpdiv (t4(0:mpnw+6), t6(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t7(0:mpnw+6), mpnw1)
!   t2 = t2 + t7
  call mpadd (t2(0:mpnw+6), t7(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!mp prec="mpnwm"
!   tc2 = abs (t7) - eps * abs (t2)
  call mpabs (t7(0:mpnw+6), mpt1(0:mpnw+6), mpnwm)
  call mpabs (t2(0:mpnw+6), mpt2(0:mpnw+6), mpnwm)
  call mpmul (eps(0:mpnwm+5), mpt2(0:mpnw+6), mpt3(0:mpnw+6), mpnwm)
  call mpsub (mpt1(0:mpnw+6), mpt3(0:mpnw+6), mpt4(0:mpnw+6), mpnwm)
  call mpeq (mpt4(0:mpnw+6), tc2(0:mpnwm+5), mpnwm)
!   ic2 = mpsgn (tc2)
  mpi1 = mpsgn (tc2(0:mpnwm+5))
  ic2 = mpi1
  if (ic2 < 0) goto 110
enddo

write (mpldb, 6)
6 format ('*** ZETABE: End loop error 2')
! call mpabrt (616)
call mpabrt (616)

110 continue

!mp prec="mpnw1"
! t1 = t0 + t2
call mpadd (t0(0:mpnw+6), t2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)

!   If original argument was negative, apply the reflection formula.

! ic1 = mpsgn (s)
mpi1 = mpsgn (s(0:))
ic1 = mpi1
if (ic1 < 0) then
!   call mpgammar (tt, t3, mpnw1)
  call mpgammar (tt, t3, mpnw1)
!   t2 = t1 * t3
  call mpmul (t1(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   t1 = mppicon * tt
  call mpmul (mppicon(0:), tt(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   t3 = t1 * 0.5e0_mpdk
  call mpmuld (t1(0:mpnw+6), 0.5e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!   t4 = cos (t3)
  call mpcssnr (t3(0:mpnw+6), mpt1(0:mpnw+6), mpt2(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t4(0:mpnw+6), mpnw1)
!   t5 = sin (t3)
  call mpcssnr (t3(0:mpnw+6), mpt2(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t5(0:mpnw+6), mpnw1)
!   t1 = t2 * t4
  call mpmul (t2(0:mpnw+6), t4(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
!   t2 = mppicon * 2.e0_mpdk
  call mpmuld (mppicon(0:), 2.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   t3 = t2 ** t3
  call mppower (t2(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t3(0:mpnw+6), mpnw1)
!   t2 = t1 / t3
  call mpdiv (t1(0:mpnw+6), t3(0:mpnw+6), mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t2(0:mpnw+6), mpnw1)
!   t1 = t2 * 2.e0_mpdk
  call mpmuld (t2(0:mpnw+6), 2.e0_mpdk, mpt1(0:mpnw+6), mpnw1)
  call mpeq (mpt1(0:mpnw+6), t1(0:mpnw+6), mpnw1)
endif

200 continue

!mp prec="mpnw"
! z = t1
call mpeq (t1(0:mpnw+6), z(0:), mpnw)

return
end subroutine mpzetabe

end module mpfune
