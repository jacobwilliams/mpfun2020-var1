!*****************************************************************************

!  MPFUN20: A thread-safe Fortran arbitrary precision computation package

!  Revision date:  9 Jan 2022

!  AUTHOR:
!     David H. Bailey
!     Lawrence Berkeley National Lab (retired) and University of California, Davis
!     Email: dhbailey@lbl.gov

!  All software in this package (c) 2022 David H. Bailey

!  PURPOSE OF THESE ROUTINES:
!    These routines define timing functions for use in the MPFUN20 package.
!
function second ()

!   The permits one to obtain elapsed run time in seconds on most Unix systems.
!   Usage example:

!      double precision second, runtime, tm0, tm1
!      external second
!      tm0 = second ()
!      call sub
!      tm1 = second ()
!      runtime = tm1 - tm0

real (kind (0.d0)) second, timereal

call cpu_time (timereal)
second = timereal
return
end

function secondwc ()
 
!   This routine, when used in conjunction with secondwe below, employs the
!   Fortran 2003 intrinsic system_clock to obtain elapsed wall clock run time.
!   Usage example:

!     double precision runtime, secondwc, secondwe
!     external secondwc, secondwe, tm0, tm1
!     tm0 = secondwc ()
!     call sub
!     tm1 = secondwc ()
!     runtime = secondwe (tm1 - tm0)

real (kind (0.d0)) secondwc
integer (selected_int_kind (18)) itm1, itm2, itm3

call system_clock (itm1, itm2, itm3)
secondwc = itm1 / dble (itm2)
return
end

function secondwe (tm)

!   This corrects the result of secondwc when wrap-around has occurred.
!   The value "2147483.648d0" here should be the same as (itm3+1)/itm2, where
!   itm2 and itm3 are the second and third arguments output by system_clock.
!   See example above.

real (kind (0.d0)) secondwe, tm

if (tm >= 0.d0) then
  secondwe = tm
else
  secondwe = tm + 2147483.648d0
endif

return
end
