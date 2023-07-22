!*****************************************************************************

function mpmask13 (b)

!  Revision date:  9 Jan 2022

!  AUTHOR:
!     David H. Bailey
!     Lawrence Berkeley National Lab (retired) and University of California, Davis
!     Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2022 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!  PURPOSE OF THIS ROUTINE:
!    This convoluted-looking code tests whether the DP value B has more than 40
!    significant bits. It actually returns the absolute value of B, with lower 13
!    bits zeroed out. This function must be compiled separately, with lower
!    optimization, since compiling with -fast with gfortran, for instance, defeats
!    the test.

use mpfuna
implicit none

real (mprknd) b, b13x, mpmask13, t1
parameter (b13x = 2.d0**13)
t1 = b13x * abs (b)
mpmask13 = abs (abs (b) + t1) - abs (t1)
return
end

function mpmask23 (b)

!  PURPOSE OF THIS ROUTINE:
!    This convoluted-looking code tests whether the QP (IEEE quad) value B has more 
!    than 90 significant bits. It actually returns the absolute value of B, with
!    lower 23 bits zeroed out. This function must be compiled separately, with lower
!    optimization, since compiling with -fast with gfortran, for instance, defeats
!    the test.

use mpfuna
implicit none

real (mprknd) b23x
parameter (b23x = 2.d0**23)
real (max (mprknd2, kind (1.0))) b, mpmask23, t1
t1 = b23x * abs (b)
mpmask23 = abs (abs (b) + t1) - abs (t1)
return
end
