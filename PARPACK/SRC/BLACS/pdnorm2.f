!\BeginDoc
!
!\Name: pdnorm2
!
! Message Passing Layer: BLACS
!
!\Description:
!
!\Usage:
!  call pdnorm2 ( COMM, N, X, INC )
!
!\Arguments
!  COMM    BLACS Communicator for the processor grid.  (INPUT)
!
!\SCCS Information:
! FILE: norm2.F   SID: 1.2   DATE OF SID: 2/22/96
!
!-----------------------------------------------------------------------
!
      Double precision function pdnorm2 ( comm, n, x, inc )
!
!     %------------------------------%
!     | BLACS Variables and Routines |
!     %------------------------------%
!
      integer    comm
      external   dgsum2d, dgamx2d
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
      Double precision
     &             x(n)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      Double precision
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
      Double precision
     &             dnrm2
      External     dnrm2
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      pdnorm2 = dnrm2( n, x, inc)
!
      max = pdnorm2
      call dgamx2d( comm, 'All', ' ', 1, 1, max, 1, ra, ca,
     &              -1, -1, -1 )
      if ( max .eq. zero ) then
         pdnorm2 = zero
      else
         pdnorm2 = (pdnorm2/max)**2.0
         call dgsum2d( comm, 'All', ' ', 1, 1, pdnorm2, 1, -1, -1 )
         pdnorm2 = max * sqrt(abs(pdnorm2))
      endif
!
!     %----------------%
!     | End of pdnorm2 |
!     %----------------%
!
      return
      end
