!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: ssconv
!
!\Description:
!  Convergence testing for the symmetric Arnoldi eigenvalue routine.
!
!\Usage:
!  call ssconv
!     ( N, RITZ, BOUNDS, TOL, NCONV )
!
!\Arguments
!  N       Integer.  (INPUT)
!          Number of Ritz values to check for convergence.
!
!  RITZ    Real array of length N.  (INPUT)
!          The Ritz values to be checked for convergence.
!
!  BOUNDS  Real array of length N.  (INPUT)
!          Ritz estimates associated with the Ritz values in RITZ.
!
!  TOL     Real scalar.  (INPUT)
!          Desired relative accuracy for a Ritz value to be considered
!          "converged".
!
!  NCONV   Integer scalar.  (OUTPUT)
!          Number of "converged" Ritz values.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Routines called:
!     arscnd  ARPACK utility routine for timing.
!     slamch  LAPACK routine that determines machine constants.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\SCCS Information: @(#)
! FILE: sconv.F   SID: 2.4   DATE OF SID: 4/19/96   RELEASE: 2
!
!\Remarks
!     1. Starting with version 2.4, this routine no longer uses the
!        Parlett strategy using the gap conditions.
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine ssconv (n, ritz, bounds, tol, nconv)
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      integer    n, nconv
      Real
     &           tol
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      Real
     &           ritz(n), bounds(n)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer    i
      Real
     &           temp, eps23
!
!     %-------------------%
!     | External routines |
!     %-------------------%
!
      Real
     &           slamch
      external   slamch

!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      call arscnd (t0)
!
      eps23 = slamch('Epsilon-Machine')
      eps23 = eps23**(2.0E+0 / 3.0E+0)
!
      nconv  = 0
      do 10 i = 1, n
!
!        %-----------------------------------------------------%
!        | The i-th Ritz value is considered "converged"       |
!        | when: bounds(i) .le. TOL*max(eps23, abs(ritz(i)))   |
!        %-----------------------------------------------------%
!
         temp = max( eps23, abs(ritz(i)) )
         if ( bounds(i) .le. tol*temp ) then
            nconv = nconv + 1
         end if
!
   10 continue
!
      call arscnd (t1)
      tsconv = tsconv + (t1 - t0)
!
      return
!
!     %---------------%
!     | End of ssconv |
!     %---------------%
!
      end
