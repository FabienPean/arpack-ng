!\BeginDoc
!
!\Name: pslarnv
!
! Message Passing Layer: BLACS
!
!\Description:
!
!  Parallel Version of ARPACK utility routine slarnv
!
!  PSLARNV returns a vector of n (nloc) random real numbers from a uniform or
!  normal distribution. It is assumed that X is distributed across a 1-D array
!  of processors ( nprocs < 1000 )
!
!\Arguments
!  COMM    BLACS Communicator for the processor grid
!
!  IDIST   (input) INTEGER
!          Specifies the distribution of the random numbers:
!          = 1:  uniform (0,1)
!          = 2:  uniform (-1,1)
!          = 3:  normal (0,1)
!
!  ISEED   (input/output) INTEGER array, dimension (4)
!          On entry, the seed of the random number generator; the array
!          elements must be between 0 and 4095, and ISEED(4) must be
!          odd.
!          On exit, the seed is updated.
!
!  N       (input) INTEGER
!          The number of random numbers to be generated.
!
!  X       (output) Real array, dimension (N)
!          The generated random numbers.
!
!\Author: Kristi Maschhoff
!
!\Details
!
!  Simple parallel version of LAPACK auxiliary routine slarnv
!  for X distributed across a 1-D array of processors.
!  This routine calls the auxiliary routine SLARNV to generate random
!  real numbers from a uniform (0,1) distribution. Output is consistent
!  with serial version.
!
!\SCCS Information:
! FILE: larnv.F   SID: 1.4   DATE OF SID: 04/16/99
!
!-----------------------------------------------------------------------
!
      subroutine pslarnv( comm, idist, iseed, n, x )
!
!     .. BLACS VARIABLES AND FUNCTIONS ..
      integer    comm
!
!     .. External Functions ..
      external   BLACS_GRIDINFO
!     ..
!     .. Scalar Arguments ..
      integer			idist, n
!     ..
!     .. Array Arguments ..
      integer			iseed( 4 )
      Real
     &                  x( * )
!     ..
!     .. External Subroutines ..
      external			slarnv
!     ..
!     .. Executable Statements ..
!
      call slarnv ( idist, iseed, n, x )
!
      return
      end
