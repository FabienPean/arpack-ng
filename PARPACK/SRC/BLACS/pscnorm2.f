!\BeginDoc
!
!\Name: pscnorm2
!
! Message Passing Layer: BLACS
!
!\Description:
!
!\Usage:
!  call pscnorm2 ( COMM, N, X, INC )
!
!\Arguments
!  COMM    BLACS Communicator for the processor grid.  (INPUT)
!
!\SCCS Information:
! FILE: norm2.F   SID: 1.2   DATE OF SID: 3/6/96
!
!-----------------------------------------------------------------------
!
      Real function pscnorm2 ( comm, n, x, inc )
!
!     %------------------------------%
!     | BLACS Variables and Routines |
!     %------------------------------%
!
      integer    comm
      external   sgsum2d, sgamx2d
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      integer      n, inc
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      Complex
     &             x(n)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      Real
     &             max, buf, zero
      parameter    ( zero = 0.0 )
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    abs, sqrt
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real
     &             scnrm2
      External     scnrm2
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      pscnorm2 = scnrm2( n, x, inc)
!
      max = pscnorm2
      call sgamx2d( comm, 'All', ' ', 1, 1, max, 1, ra, ca,
     &              -1, -1, -1 )
      if ( max .eq. zero ) then
         pscnorm2 = zero
      else
         pscnorm2 = (pscnorm2/max)**2.0
         call sgsum2d( comm, 'All', ' ', 1, 1, pscnorm2, 1, -1, -1 )
         pscnorm2 = max * sqrt(abs(pscnorm2))
      endif
!
!     %-----------------%
!     | End of pscnorm2 |
!     %-----------------%
!
      return
      end
