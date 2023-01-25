!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: snconv
!
!\Description:
!  Convergence testing for the nonsymmetric Arnoldi eigenvalue routine.
!
!\Usage:
!  call snconv
!     ( N, RITZR, RITZI, BOUNDS, TOL, NCONV )
!
!\Arguments
!  N       Integer.  (INPUT)
!          Number of Ritz values to check for convergence.
!
!  RITZR,  Real arrays of length N.  (INPUT)
!  RITZI   Real and imaginary parts of the Ritz values to be checked
!          for convergence.

!  BOUNDS  Real array of length N.  (INPUT)
!          Ritz estimates for the Ritz values in RITZR and RITZI.
!
!  TOL     Real scalar.  (INPUT)
!          Desired backward error for a Ritz value to be considered
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
!\Local variables:
!     xxxxxx  real
!
!\Routines called:
!     arscnd  ARPACK utility routine for timing.
!     slamch  LAPACK routine that determines machine constants.
!     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     xx/xx/92: Version ' 2.1'
!
!\SCCS Information: @(#)
! FILE: nconv.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
!
!\Remarks
!     1. xxxx
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine snconv (n, ritzr, ritzi, bounds, tol, nconv)
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

      Real
     &           ritzr(n), ritzi(n), bounds(n)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer    i
      Real
     &           temp, eps23
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real
     &           slapy2, slamch
      external   slapy2, slamch

!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %-------------------------------------------------------------%
!     | Convergence test: unlike in the symmetric code, I am not    |
!     | using things like refined error bounds and gap condition    |
!     | because I don't know the exact equivalent concept.          |
!     |                                                             |
!     | Instead the i-th Ritz value is considered "converged" when: |
!     |                                                             |
!     |     bounds(i) .le. ( TOL * | ritz | )                       |
!     |                                                             |
!     | for some appropriate choice of norm.                        |
!     %-------------------------------------------------------------%
!
      call arscnd (t0)
!
!     %---------------------------------%
!     | Get machine dependent constant. |
!     %---------------------------------%
!
      eps23 = slamch('Epsilon-Machine')
      eps23 = eps23**(2.0E+0 / 3.0E+0)
!
      nconv  = 0
      do 20 i = 1, n
         temp = max( eps23, slapy2( ritzr(i), ritzi(i) ) )
         if (bounds(i) .le. tol*temp)   nconv = nconv + 1
   20 continue
!
      call arscnd (t1)
      tnconv = tnconv + (t1 - t0)
!
      return
!
!     %---------------%
!     | End of snconv |
!     %---------------%
!
      end
