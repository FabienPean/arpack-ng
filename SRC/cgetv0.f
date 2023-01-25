!\BeginDoc
!
!\Name: cgetv0
!
!\Description:
!  Generate a random initial residual vector for the Arnoldi process.
!  Force the residual vector to be in the range of the operator OP.
!
!\Usage:
!  call cgetv0
!     ( IDO, BMAT, ITRY, INITV, N, J, V, LDV, RESID, RNORM,
!       IPNTR, WORKD, IERR )
!
!\Arguments
!  IDO     Integer.  (INPUT/OUTPUT)
!          Reverse communication flag.  IDO must be zero on the first
!          call to cgetv0.
!          -------------------------------------------------------------
!          IDO =  0: first call to the reverse communication interface
!          IDO = -1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!                    This is for the initialization phase to force the
!                    starting vector into the range of OP.
!          IDO =  2: compute  Y = B * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!          IDO = 99: done
!          -------------------------------------------------------------
!
!  BMAT    Character*1.  (INPUT)
!          BMAT specifies the type of the matrix B in the (generalized)
!          eigenvalue problem A*x = lambda*B*x.
!          B = 'I' -> standard eigenvalue problem A*x = lambda*x
!          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
!
!  ITRY    Integer.  (INPUT)
!          ITRY counts the number of times that cgetv0 is called.
!          It should be set to 1 on the initial call to cgetv0.
!
!  INITV   Logical variable.  (INPUT)
!          .TRUE.  => the initial residual vector is given in RESID.
!          .FALSE. => generate a random initial residual vector.
!
!  N       Integer.  (INPUT)
!          Dimension of the problem.
!
!  J       Integer.  (INPUT)
!          Index of the residual vector to be generated, with respect to
!          the Arnoldi process.  J > 1 in case of a "restart".
!
!  V       Complex N by J array.  (INPUT)
!          The first J-1 columns of V contain the current Arnoldi basis
!          if this is a "restart".
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  RESID   Complex array of length N.  (INPUT/OUTPUT)
!          Initial residual vector to be generated.  If RESID is
!          provided, force RESID into the range of the operator OP.
!
!  RNORM   Real scalar.  (OUTPUT)
!          B-norm of the generated residual.
!
!  IPNTR   Integer array of length 3.  (OUTPUT)
!
!  WORKD   Complex work array of length 2*N.  (REVERSE COMMUNICATION).
!          On exit, WORK(1:N) = B*RESID to be used in SSAITR.
!
!  IERR    Integer.  (OUTPUT)
!          =  0: Normal exit.
!          = -1: Cannot generate a nontrivial restarted residual vector
!                in the range of the operator OP.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  Complex
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!
!\Routines called:
!     arscnd  ARPACK utility routine for timing.
!     cvout   ARPACK utility routine that prints vectors.
!     clarnv  LAPACK routine for generating a random vector.
!     cgemv   Level 2 BLAS routine for matrix vector multiplication.
!     ccopy   Level 1 BLAS that copies one vector to another.
!     cdotc   Level 1 BLAS that computes the scalar product of two vectors.
!     scnrm2  Level 1 BLAS that computes the norm of a vector.
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
! FILE: getv0.F   SID: 2.3   DATE OF SID: 08/27/96   RELEASE: 2
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine cgetv0
     &   ( ido, bmat, itry, initv, n, j, v, ldv, resid, rnorm,
     &     ipntr, workd, ierr )
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
      character  bmat*1
      logical    initv
      integer    ido, ierr, itry, j, ldv, n
      Real
     &           rnorm
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer    ipntr(3)
      Complex
     &           resid(n), v(ldv,j), workd(2*n)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Complex
     &           one, zero
      Real
     &           rzero
      parameter  (one = (1.0E+0, 0.0E+0), zero = (0.0E+0, 0.0E+0),
     &            rzero = 0.0E+0)
!
!     %------------------------%
!     | Local Scalars & Arrays |
!     %------------------------%
!
      logical    first, orth
      integer    idist, iseed(4), iter, msglvl, jj
      Real
     &           rnorm0
      Complex
     &           cnorm
      save       first, iseed, iter, msglvl, orth, rnorm0
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   ccopy, cgemv, clarnv, cvout, arscnd
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real
     &           scnrm2, slapy2
      Complex
     &           ccdotc
      external   ccdotc, scnrm2, slapy2
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!
!     %-----------------------------------%
!     | Initialize the seed of the LAPACK |
!     | random number generator           |
!     %-----------------------------------%
!
      iseed(1) = 1
      iseed(2) = 3
      iseed(3) = 5
      iseed(4) = 7
!
      if (ido .eq.  0) then
!
!        %-------------------------------%
!        | Initialize timing statistics  |
!        | & message level for debugging |
!        %-------------------------------%
!
         call arscnd (t0)
         msglvl = mgetv0
!
         ierr   = 0
         iter   = 0
         first  = .FALSE.
         orth   = .FALSE.
!
!        %-----------------------------------------------------%
!        | Possibly generate a random starting vector in RESID |
!        | Use a LAPACK random number generator used by the    |
!        | matrix generation routines.                         |
!        |    idist = 1: uniform (0,1)  distribution;          |
!        |    idist = 2: uniform (-1,1) distribution;          |
!        |    idist = 3: normal  (0,1)  distribution;          |
!        %-----------------------------------------------------%
!
         if (.not.initv) then
            idist = 2
            call clarnv (idist, iseed, n, resid)
         end if
!
!        %----------------------------------------------------------%
!        | Force the starting vector into the range of OP to handle |
!        | the generalized problem when B is possibly (singular).   |
!        %----------------------------------------------------------%
!
         call arscnd (t2)
         if (itry .eq. 1) then
            nopx = nopx + 1
            ipntr(1) = 1
            ipntr(2) = n + 1
            call ccopy (n, resid, 1, workd, 1)
            ido = -1
            go to 9000
         else if (itry .gt. 1 .and. bmat .eq. 'G') then
            call ccopy (n, resid, 1, workd(n + 1), 1)
         end if
      end if
!
!     %----------------------------------------%
!     | Back from computing OP*(initial-vector) |
!     %----------------------------------------%
!
      if (first) go to 20
!
!     %-----------------------------------------------%
!     | Back from computing OP*(orthogonalized-vector) |
!     %-----------------------------------------------%
!
      if (orth)  go to 40
!
      call arscnd (t3)
      tmvopx = tmvopx + (t3 - t2)
!
!     %------------------------------------------------------%
!     | Starting vector is now in the range of OP; r = OP*r; |
!     | Compute B-norm of starting vector.                   |
!     %------------------------------------------------------%
!
      call arscnd (t2)
      first = .TRUE.
      if (itry .eq. 1) call ccopy (n, workd(n + 1), 1, resid, 1)
      if (bmat .eq. 'G') then
         nbx = nbx + 1
         ipntr(1) = n + 1
         ipntr(2) = 1
         ido = 2
         go to 9000
      else if (bmat .eq. 'I') then
         call ccopy (n, resid, 1, workd, 1)
      end if
!
   20 continue
!
      if (bmat .eq. 'G') then
         call arscnd (t3)
         tmvbx = tmvbx + (t3 - t2)
      end if
!
      first = .FALSE.
      if (bmat .eq. 'G') then
          cnorm  = ccdotc (n, resid, 1, workd, 1)
          rnorm0 = sqrt(slapy2(real(cnorm),aimag(cnorm)))
      else if (bmat .eq. 'I') then
           rnorm0 = scnrm2(n, resid, 1)
      end if
      rnorm  = rnorm0
!
!     %---------------------------------------------%
!     | Exit if this is the very first Arnoldi step |
!     %---------------------------------------------%
!
      if (j .eq. 1) go to 50
!
!     %----------------------------------------------------------------
!     | Otherwise need to B-orthogonalize the starting vector against |
!     | the current Arnoldi basis using Gram-Schmidt with iter. ref.  |
!     | This is the case where an invariant subspace is encountered   |
!     | in the middle of the Arnoldi factorization.                   |
!     |                                                               |
!     |       s = V^{T}*B*r;   r = r - V*s;                           |
!     |                                                               |
!     | Stopping criteria used for iter. ref. is discussed in         |
!     | Parlett's book, page 107 and in Gragg & Reichel TOMS paper.   |
!     %---------------------------------------------------------------%
!
      orth = .TRUE.
   30 continue
!
      call cgemv ('C', n, j-1, one, v, ldv, workd, 1,
     &            zero, workd(n+1), 1)
      call cgemv ('N', n, j-1, -one, v, ldv, workd(n+1), 1,
     &            one, resid, 1)
!
!     %----------------------------------------------------------%
!     | Compute the B-norm of the orthogonalized starting vector |
!     %----------------------------------------------------------%
!
      call arscnd (t2)
      if (bmat .eq. 'G') then
         nbx = nbx + 1
         call ccopy (n, resid, 1, workd(n+1), 1)
         ipntr(1) = n + 1
         ipntr(2) = 1
         ido = 2
         go to 9000
      else if (bmat .eq. 'I') then
         call ccopy (n, resid, 1, workd, 1)
      end if
!
   40 continue
!
      if (bmat .eq. 'G') then
         call arscnd (t3)
         tmvbx = tmvbx + (t3 - t2)
      end if
!
      if (bmat .eq. 'G') then
         cnorm = ccdotc (n, resid, 1, workd, 1)
         rnorm = sqrt(slapy2(real(cnorm),aimag(cnorm)))
      else if (bmat .eq. 'I') then
         rnorm = scnrm2(n, resid, 1)
      end if
!
!     %--------------------------------------%
!     | Check for further orthogonalization. |
!     %--------------------------------------%
!
      if (msglvl .gt. 2) then
          call svout (logfil, 1, [rnorm0], ndigit,
     &                '_getv0: re-orthonalization ; rnorm0 is')
          call svout (logfil, 1, [rnorm], ndigit,
     &                '_getv0: re-orthonalization ; rnorm is')
      end if
!
      if (rnorm .gt. 0.717*rnorm0) go to 50
!
      iter = iter + 1
      if (iter .le. 1) then
!
!        %-----------------------------------%
!        | Perform iterative refinement step |
!        %-----------------------------------%
!
         rnorm0 = rnorm
         go to 30
      else
!
!        %------------------------------------%
!        | Iterative refinement step "failed" |
!        %------------------------------------%
!
         do 45 jj = 1, n
            resid(jj) = zero
   45    continue
         rnorm = rzero
         ierr = -1
      end if
!
   50 continue
!
      if (msglvl .gt. 0) then
         call svout (logfil, 1, [rnorm], ndigit,
     &        '_getv0: B-norm of initial / restarted starting vector')
      end if
      if (msglvl .gt. 2) then
         call cvout (logfil, n, resid, ndigit,
     &        '_getv0: initial / restarted starting vector')
      end if
      ido = 99
!
      call arscnd (t1)
      tgetv0 = tgetv0 + (t1 - t0)
!
 9000 continue
      return
!
!     %---------------%
!     | End of cgetv0 |
!     %---------------%
!
      end
