! \BeginDoc
!
! \Name: snband
!
! \Description:
!
!  This subroutine returns the converged approximations to eigenvalues
!  of A*z = lambda*B*z and (optionally):
!
!      (1) The corresponding approximate eigenvectors;
!
!      (2) An orthonormal basis for the associated approximate
!          invariant subspace;
!
!      (3) Both.
!
!  Matrices A and B are stored in LAPACK-style banded form.
!
!  There is negligible additional cost to obtain eigenvectors.  An orthonormal
!  basis is always computed.  There is an additional storage cost of n*nev
!  if both are requested (in this case a separate array Z must be supplied).
!
!  The approximate eigenvalues and vectors are commonly called Ritz
!  values and Ritz vectors respectively.  They are referred to as such
!  in the comments that follow.  The computed orthonormal basis for the
!  invariant subspace corresponding to these Ritz values is referred to as a
!  Schur basis.
!
!  snband can be called with one of the following modes:
!
!  Mode 1:  A*z = lambda*z.
!           ===> OP = A  and  B = I.
!
!  Mode 2:  A*z = lambda*M*z, M symmetric positive definite
!           ===> OP = inv[M]*A  and  B = M.
!
!  Mode 3:  A*z = lambda*M*z, M symmetric semi-definite
!           ===> OP = Real_Part{ inv[A - sigma*M]*M }  and  B = M.
!           ===> shift-and-invert mode (in real arithmetic)
!           If OP*z = amu*z, then
!           amu = 1/2 * [ 1/(lambda-sigma) + 1/(lambda-conjg(sigma)) ].
!           Note: If sigma is real, i.e. imaginary part of sigma is zero;
!                 Real_Part{ inv[A - sigma*M]*M } == inv[A - sigma*M]*M
!                 amu == 1/(lambda-sigma).
!
!  Mode 4:  A*z = lambda*M*z, M symmetric semi-definite
!           ===> OP = Imaginary_Part{ inv[A - sigma*M]*M }  and  B = M.
!           ===> shift-and-invert mode (in real arithmetic)
!           If OP*z = amu*z, then
!           amu = 1/2i * [ 1/(lambda-sigma) - 1/(lambda-conjg(sigma)) ].
!
!
!  The choice of mode must be specified in IPARAM(7) defined below.
!
! \Usage
!   call snband
!      ( RVEC, HOWMNY, SELECT, DR, DI, Z, LDZ, SIGMAR, SIGMAI,
!        WORKEV, V, N, AB, MB, LDA, RFAC, CFAC, KL, KU, WHICH,
!        BMAT, NEV, TOL, RESID, NCV, V, LDV, IPARAM, WORKD,
!        WORKL, LWORKL, WORKC, IWORK, INFO )
!
! \Arguments
!
!  RVEC    LOGICAL  (INPUT)
!          Specifies whether a basis for the invariant subspace corresponding
!          to the converged Ritz value approximations for the eigenproblem
!          A*z = lambda*B*z is computed.
!
!             RVEC = .FALSE.     Compute Ritz values only.
!
!             RVEC = .TRUE.      Compute the Ritz vectors or Schur vectors.
!                                See Remarks below.
!
!  HOWMNY  Character*1  (INPUT)
!          Specifies the form of the basis for the invariant subspace
!          corresponding to the converged Ritz values that is to be computed.
!
!          = 'A': Compute NEV Ritz vectors;
!          = 'P': Compute NEV Schur vectors;
!          = 'S': compute some of the Ritz vectors, specified
!                 by the logical array SELECT.
!
!  SELECT  Logical array of dimension NCV.  (INPUT)
!          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
!          computed. To select the Ritz vector corresponding to a
!          Ritz value (DR(j), DI(j)), SELECT(j) must be set to .TRUE..
!          If HOWMNY = 'A' or 'P', SELECT is used as internal workspace.
!
!  DR      Real array of dimension NEV+1.  (OUTPUT)
!          On exit, DR contains the real part of the Ritz value approximations
!          to the eigenvalues of A*z = lambda*B*z.
!
!  DI      Real array of dimension NEV+1.  (OUTPUT)
!          On exit, DI contains the imaginary part of the Ritz value
!          approximations to the eigenvalues of A*z = lambda*B*z associated
!          with DR.
!
!          NOTE: When Ritz values are complex, they will come in complex
!                conjugate pairs.  If eigenvectors are requested, the
!                corresponding Ritz vectors will also come in conjugate
!                pairs and the real and imaginary parts of these are
!                represented in two consecutive columns of the array Z
!                (see below).
!
!  Z       Real N by NEV+1 array if RVEC = .TRUE. and HOWMNY = 'A'. (OUTPUT)
!          On exit,
!          if RVEC = .TRUE. and HOWMNY = 'A', then the columns of
!          Z represent approximate eigenvectors (Ritz vectors) corresponding
!          to the NCONV=IPARAM(5) Ritz values for eigensystem
!          A*z = lambda*B*z computed by SNAUPD.
!
!          The complex Ritz vector associated with the Ritz value
!          with positive imaginary part is stored in two consecutive
!          columns.  The first column holds the real part of the Ritz
!          vector and the second column holds the imaginary part.  The
!          Ritz vector associated with the Ritz value with negative
!          imaginary part is simply the complex conjugate of the Ritz vector
!          associated with the positive imaginary part.
!
!          If  RVEC = .FALSE. or HOWMNY = 'P', then Z is not referenced.
!
!          NOTE: If if RVEC = .TRUE. and a Schur basis is not required,
!          the array Z may be set equal to first NEV+1 columns of the Arnoldi
!          basis array V computed by SNAUPD.  In this case the Arnoldi basis
!          will be destroyed and overwritten with the eigenvector basis.
!
!  LDZ     Integer.  (INPUT)
!          The leading dimension of the array Z.  If Ritz vectors are
!          desired, then  LDZ >= max( 1, N ).  In any case,  LDZ >= 1.
!
!  SIGMAR  Real  (INPUT)
!          If IPARAM(7) = 3 or 4, represents the real part of the shift.
!          Not referenced if IPARAM(7) = 1 or 2.
!
!  SIGMAI  Real  (INPUT)
!          If IPARAM(7) = 3 or 4, represents the imaginary part of the
!          shift.
!          Not referenced if IPARAM(7) = 1 or 2.
!
!  WORKEV  Real work array of dimension 3*NCV.  (WORKSPACE)
!
!  N       Integer.  (INPUT)
!          Dimension of the eigenproblem.
!
!  AB      Real array of dimension LDA by N. (INPUT)
!          The matrix A in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!
!  MB      Real array of dimension LDA by N. (INPUT)
!          The matrix M in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of M is stored in the j-th column of the
!          array AB as follows:
!          MB(kl+ku+1+i-j,j) = M(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!          Not referenced if IPARAM(7) = 1
!
!  LDA     Integer. (INPUT)
!          Leading dimension of AB, MB, RFAC and CFAC.
!
!  RFAC    Real array of LDA by N. (WORKSPACE/OUTPUT)
!          RFAC is used to store the LU factors of MB when IPARAM(7) = 2
!          is invoked.  It is used to store the LU factors of
!          (A-sigma*M) when IPARAM(7) = 3 is invoked with a real shift.
!          It is not referenced when IPARAM(7) = 1 or 4.
!
!  CFAC    Complex array of LDA by N. (WORKSPACE/OUTPUT)
!          CFAC is used to store (A-SIGMA*M) and its LU factors
!          when IPARAM(7) = 3 or 4 are used with a complex shift SIGMA.
!          On exit, it contains the LU factors of (A-SIGMA*M).
!          It is not referenced when IPARAM(7) = 1 or 2.
!
!  KL      Integer. (INPUT)
!          Max(number of subdiagonals of A, number of subdiagonals of M)
!
!  KU      Integer. (OUTPUT)
!          Max(number of superdiagonals of A, number of superdiagonals of M)
!
!  WHICH   Character*2.  (INPUT)
!          When IPARAM(7)= 1 or 2,  WHICH can be set to any one of
!          the following.
!
!            'LM' -> want the NEV eigenvalues of largest magnitude.
!            'SM' -> want the NEV eigenvalues of smallest magnitude.
!            'LR' -> want the NEV eigenvalues of largest real part.
!            'SR' -> want the NEV eigenvalues of smallest real part.
!            'LI' -> want the NEV eigenvalues of largest imaginary part.
!            'SI' -> want the NEV eigenvalues of smallest imaginary part.
!
!          When IPARAM(7) = 3 or 4, WHICH should be set to 'LM' only.
!
!  BMAT    Character*1.  (INPUT)
!          BMAT specifies the type of the matrix B that defines the
!          semi-inner product for the operator OP.
!          BMAT = 'I' -> standard eigenvalue problem A*z = lambda*z
!          BMAT = 'G' -> generalized eigenvalue problem A*z = lambda*M*z

!  NEV     Integer. (INPUT)
!          Number of eigenvalues to be computed.
!
!  TOL     Real scalar.  (INPUT)
!          Stopping criteria: the relative accuracy of the Ritz value
!          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I)).
!          If TOL .LE. 0. is passed a default is set:
!          DEFAULT = SLAMCH('EPS')  (machine precision as computed
!                    by the LAPACK auxiliary subroutine SLAMCH).
!
!  RESID   Real array of length N.  (INPUT/OUTPUT)
!          On INPUT:
!          If INFO .EQ. 0, a random initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          On OUTPUT:
!          RESID contains the final residual vector.
!
!  NCV     Integer.  (INPUT)
!          Number of columns of the matrix V (less than or equal to N).
!          Represents the dimension of the Arnoldi basis constructed
!          by snaupd for OP.
!
!  V       Real array N by NCV+1.  (OUTPUT)
!          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
!                       represent approximate Schur vectors that span the
!                       desired invariant subspace.
!          NOTE: The array Z may be set equal to first NEV+1 columns of the
!          Arnoldi basis vector array V computed by SNAUPD. In this case
!          if RVEC = .TRUE. and HOWMNY='A', then the first NCONV=IPARAM(5)
!          are the desired Ritz vectors.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
!          IPARAM(1) = ISHIFT:
!          The shifts selected at each iteration are used to restart
!          the Arnoldi iteration in an implicit fashion.
!          It is set to 1 in this subroutine.  The user do not need
!          to set this parameter.
!           ----------------------------------------------------------
!          ISHIFT = 1: exact shift with respect to the current
!                      Hessenberg matrix H.  This is equivalent to
!                      restarting the iteration from the beginning
!                      after updating the starting vector with a linear
!                      combination of Ritz vectors associated with the
!                      "wanted" eigenvalues.
!          -------------------------------------------------------------
!
!          IPARAM(2) = No longer referenced.
!
!          IPARAM(3) = MXITER
!          On INPUT:  max number of Arnoldi update iterations allowed.
!          On OUTPUT: actual number of Arnoldi update iterations taken.
!
!          IPARAM(4) = NB: blocksize to be used in the recurrence.
!          The code currently works only for NB = 1.
!
!          IPARAM(5) = NCONV: number of "converged" eigenvalues.
!
!          IPARAM(6) = IUPD
!          Not referenced. Implicit restarting is ALWAYS used.
!
!          IPARAM(7) = IPARAM(7):
!          On INPUT determines what type of eigenproblem is being solved.
!          Must be 1,2,3,4; See under \Description of snband for the
!          four modes available.
!
!          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
!          OUTPUT: NUMOP  = total number of OP*z operations,
!                  NUMOPB = total number of B*z operations if BMAT='G',
!                  NUMREO = total number of steps of re-orthogonalization.
!
! WORKD    Real work array of length at least 3*n. (WORKSPACE)
!
! WORKL    Real work array of length LWORKL. (WORKSPACE)
!
! LWORKL   Integer.  (INPUT)
!          LWORKL must be at least 3*NCV**2 + 6*NCV.
!
! WORKC    Complex array of length N. (WORKSPACE)
!          Workspace used when IPARAM(7) = 3 or 4 for storing a temporary
!          complex vector.
!
! IWORK    Integer array of dimension at least N. (WORKSPACE)
!          Used when IPARAM(7)=2,3,4 to store the pivot information in the
!          factorization of M or (A-SIGMA*M).
!
! INFO     Integer.  (INPUT/OUTPUT)
!          Error flag on output.
!          =  0: Normal exit.
!          =  1: The Schur form computed by LAPACK routine slahqr
!                could not be reordered by LAPACK routine strsen.
!                Re-enter subroutine SNEUPD with IPARAM(5)=NCV and
!                increase the size of the arrays DR and DI to have
!                dimension at least NCV and allocate at least NCV
!                columns for Z. NOTE: Not necessary if Z and V share
!                the same space. Please notify the authors.
!
!          = -1: N must be positive.
!          = -2: NEV must be positive.
!          = -3: NCV-NEV >= 2 and less than or equal to N.
!          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
!          = -6: BMAT must be one of 'I' or 'G'.
!          = -7: Length of private work WORKL array is not sufficient.
!          = -8: Error return from calculation of a real Schur form.
!                Informational error from LAPACK routine slahqr.
!          = -9: Error return from calculation of eigenvectors.
!                Informational error from LAPACK routine strevc.
!          = -10: IPARAM(7) must be 1,2,3,4.
!          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
!          = -12: HOWMNY = 'S' not yet implemented
!          = -13: HOWMNY must be one of 'A' or 'P'
!          = -14: SNAUPD did not find any eigenvalues to sufficient
!                 accuracy.
!          = -15: Overflow occurs when we try to transform the Ritz
!                 values returned from SNAUPD to those of the original
!                 problem using Rayleigh Quotient.
!          = -9999: Could not build an Arnoldi factorization.
!                   IPARAM(5) returns the size of the current
!                   Arnoldi factorization.
!
! \EndDoc
!
!------------------------------------------------------------------------
!
!\BeginLib
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
!     Restarted Arnoldi Iteration", Ph.D thesis, TR95-13, Rice Univ,
!     May 1995.
!
!\Routines called:
!     snaupd  ARPACK reverse communication interface routine.
!     sneupd  ARPACK routine that returns Ritz values and (optionally)
!             Ritz vectors.
!     sgbtrf  LAPACK band matrix factorization routine.
!     sgbtrs  LAPACK band linear system solve routine.
!     cgbtrf  LAPACK complex band matrix factorization routine.
!     cgbtrs  LAPACK complex linear system solve routine.
!     slacpy  LAPACK matrix copy routine.
!     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     slamch  LAPACK routine to compute the underflow threshold.
!     scopy   Level 1 BLAS that copies one vector to another.
!     sdot    Level 1 BLAS that computes the dot product of two vectors.
!     snrm2   Level 1 BLAS that computes the norm of a vector.
!     sgbmv   Level 2 BLAS that computes the band matrix vector product.
!
!\Remarks
!
!  1. Currently only HOWMNY = 'A' and 'P' are implemented.
!
!     Let X' denote the transpose of X.
!
!  2. Schur vectors are an orthogonal representation for the basis of
!     Ritz vectors. Thus, their numerical properties are often superior.
!     If RVEC = .TRUE. then the relationship
!             A * V(:,1:IPARAM(5)) = V(:,1:IPARAM(5)) * T, and
!     V(:,1:IPARAM(5))' * V(:,1:IPARAM(5)) = I are approximately satisfied.
!     Here T is the leading submatrix of order IPARAM(5) of the real
!     upper quasi-triangular matrix stored workl(ipntr(12)). That is,
!     T is block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
!     each 2-by-2 diagonal block has its diagonal elements equal and its
!     off-diagonal elements of opposite sign.  Corresponding to each 2-by-2
!     diagonal block is a complex conjugate pair of Ritz values. The real
!     Ritz values are stored on the diagonal of T.
!
!\Author
!     Danny Sorensen
!     Richard Lehoucq
!     Chao Yang
!     Dept. of Computational &
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\SCCS Information: @(#)
! FILE: nband.F   SID: 2.3   DATE OF SID: 10/17/00   RELEASE: 2
!
!\EndLib
!
!---------------------------------------------------------------------
!
      subroutine snband( rvec, howmny, select, dr, di, z, ldz,  sigmar,
     &           sigmai, workev, n, ab, mb, lda, rfac,  cfac, kl, ku,
     &           which, bmat, nev, tol, resid,  ncv, v, ldv,
     &           iparam, workd, workl, lworkl, workc, iwork, info)
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character        which*2, bmat, howmny
      integer          n, lda, kl, ku, nev, ncv, ldv,
     &                 ldz, lworkl, info
      Real
     &                 tol, sigmar, sigmai
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer          iparam(*), iwork(*)
      logical          select(*)
      Real
     &                 dr(*), di(*), resid(*), v(ldv,*), z(ldz,*),
     &                 ab(lda,*), mb(lda,*), rfac(lda,*),
     &                 workd(*), workl(*), workev(*)
      Complex
     &                 cfac(lda,*), workc(*)
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
      integer          ipntr(14)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer          ido, i, j, type, imid, itop, ibot, ierr
      Real
     &                 numr, denr, deni, dmdul, safmin
      logical          rvec, first
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real
     &                  one, zero
      parameter        (one = 1.0E+0, zero = 0.0E+0)
!
!
!     %-----------------------------%
!     | LAPACK & BLAS routines used |
!     %-----------------------------%
!
      Real
     &                 sdot, snrm2, slapy2, slamch
      external         sdot, scopy, sgbmv, cgbtrf, cgbtrs, sgbtrf,
     &                 sgbtrs, snrm2, slapy2, slacpy, slamch
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      Intrinsic        real, aimag, cmplx
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %--------------------------------%
!     | safmin = safe minimum is such  |
!     | that 1/sfmin does not overflow |
!     %--------------------------------%
!
      safmin = slamch('safmin')
!
!     %----------------------------------------------------------------%
!     | Set type of the problem to be solved. Check consistency        |
!     | between BMAT and IPARAM(7).                                    |
!     | type = 1 --> Solving standard problem in regular mode.         |
!     | type = 2 --> Solving standard problem in shift-invert mode.    |
!     | type = 3 --> Solving generalized problem in regular mode.      |
!     | type = 4 --> Solving generalized problem in shift-invert mode. |
!     | type = 5 --> Solving standard problem in shift-invert mode     |
!     |              using iparam(7) = 4 in SNAUPD.                    |
!     | type = 6 --> Solving generalized problem in shift-invert mode. |
!     |              using iparam(7) = 4 in SNAUPD.                    |
!     %----------------------------------------------------------------%
!
      if ( iparam(7) .eq. 1 ) then
         type = 1
      else if ( iparam(7) .eq. 3 .and. bmat .eq. 'I') then
         type = 2
      else if ( iparam(7) .eq. 2 ) then
         type = 3
      else if ( iparam(7) .eq. 3 .and. bmat .eq. 'G') then
         type = 4
      else if ( iparam(7) .eq. 4 .and. bmat .eq. 'I') then
         type = 5
      else if ( iparam(7) .eq. 4 .and. bmat .eq. 'G') then
         type = 6
      else
         print*, ' '
         print*, 'BMAT is inconsistent with IPARAM(7).'
         print*, ' '
         go to 9000
      end if
!
!     %----------------------------------%
!     | When type = 5,6 are used, sigmai |
!     | must be nonzero.                 |
!     %----------------------------------%
!
      if ( type .eq. 5 .or. type .eq. 6 ) then
          if ( sigmai .eq. zero ) then
             print*, ' '
             print*, '_NBAND: sigmai must be nonzero when type 5 or 6
     &                is used. '
             print*, ' '
             go to 9000
          end if
      end if
!
!     %------------------------%
!     | Initialize the reverse |
!     | communication flag.    |
!     %------------------------%
!
      ido   = 0
!
!     %----------------%
!     | Exact shift is |
!     | used.          |
!     %----------------%
!
      iparam(1) = 1
!
!     %-----------------------------------%
!     | Both matrices A and M are stored  |
!     | between rows itop and ibot.  Imid |
!     | is the index of the row that      |
!     | stores the diagonal elements.     |
!     %-----------------------------------%
!
      itop = kl + 1
      imid = kl + ku + 1
      ibot = 2*kl + ku + 1
!
      if ( type .eq. 2 .or. type .eq. 5 ) then
!
!         %-------------------------------%
!         | Solving a standard eigenvalue |
!         | problem in shift-invert mode. |
!         | Factor (A-sigma*I).           |
!         %-------------------------------%
!
          if (sigmai .eq. zero) then
!
!            %-----------------------------------%
!            | Construct (A-sigmar*I) and factor |
!            | in real arithmetic.               |
!            %-----------------------------------%
!
             call slacpy ('A', ibot, n, ab, lda, rfac, lda )
             do 10 j = 1, n
                rfac(imid,j) =  ab(imid,j) - sigmar
  10         continue
             call sgbtrf(n, n, kl, ku, rfac, lda, iwork, ierr )
             if (ierr .ne. 0) then
                print*, ' '
                print*, ' _NBAND: Error with _gbtrf. '
                print*, ' '
                go to  9000
             end if
!
          else
!
!            %-----------------------------------%
!            | Construct (A-sigmar*I) and factor |
!            | in COMPLEX arithmetic.            |
!            %-----------------------------------%
!
             do 30 j = 1, n
                do 20 i = itop, ibot
                   cfac(i,j) = cmplx(ab(i,j))
  20            continue
  30         continue
!
             do 40 j = 1, n
                cfac(imid,j) = cfac(imid,j)
     $                         - cmplx(sigmar, sigmai)
  40         continue
!
             call cgbtrf(n, n, kl, ku, cfac, lda, iwork, ierr )
             if ( ierr .ne. 0) then
                print*, ' '
                print*, ' _NBAND: Error with _gbtrf. '
                print*, ' '
                go to  9000
             end if
!
          end if

      else if ( type .eq. 3 ) then
!
!        %-----------------------------------------------%
!        | Solving generalized eigenvalue problem in     |
!        | regular mode. Copy M to rfac, and call LAPACK |
!        | routine sgbtrf to factor M.                   |
!        %-----------------------------------------------%
!
         call slacpy ('A', ibot, n, mb, lda, rfac, lda )
         call sgbtrf(n, n, kl, ku, rfac, lda, iwork, ierr)
         if (ierr .ne. 0) then
             print*, ' '
             print*,'_NBAND:  Error with _gbtrf.'
             print*, ' '
             go to 9000
         end if
!
      else if ( type .eq. 4 .or. type .eq. 6 ) then
!
!        %-------------------------------------------%
!        | Solving generalized eigenvalue problem in |
!        | shift-invert mode.                        |
!        %-------------------------------------------%
!
         if ( sigmai .eq. zero ) then
!
!            %--------------------------------------------%
!            | Construct (A - sigma*M) and factor in real |
!            | arithmetic.                                |
!            %--------------------------------------------%
!
             do 60 j = 1,n
                do 50 i = itop, ibot
                   rfac(i,j) = ab(i,j) - sigmar*mb(i,j)
  50            continue
  60         continue
!
             call sgbtrf(n, n, kl, ku, rfac, lda, iwork, ierr)
             if ( ierr .ne. 0 )  then
                 print*, ' '
                 print*, '_NBAND: Error with _gbtrf.'
                 print*, ' '
                 go to 9000
             end if
!
         else
!
!            %-----------------------------------------------%
!            | Construct (A - sigma*M) and factor in complex |
!            | arithmetic.                                   |
!            %-----------------------------------------------%
!
             do 80 j = 1,n
                do 70 i = itop, ibot
                   cfac(i,j) = cmplx( ab(i,j)-sigmar*mb(i,j),
     &                         -sigmai*mb(i,j) )
  70            continue
  80         continue
!
             call cgbtrf(n, n, kl, ku, cfac, lda, iwork, ierr)
             if ( ierr .NE. 0 )  then
                print*, ' '
                print*, '_NBAND: Error with _gbtrf.'
                print*, ' '
                go to 9000
             end if
!
         end if
!
      end if
!
!     %--------------------------------------------%
!     |  M A I N   L O O P (reverse communication) |
!     %--------------------------------------------%
!
  90  continue
!
      call snaupd ( ido, bmat, n, which, nev, tol, resid, ncv,
     &              v, ldv, iparam, ipntr, workd, workl, lworkl,
     &              info )
!
      if (ido .eq. -1) then
!
         if ( type .eq. 1) then
!
!           %----------------------------%
!           | Perform  y <--- OP*x = A*x |
!           %----------------------------%
!
            call sgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1),
     &                 lda, workd(ipntr(1)), 1, zero,
     &                 workd(ipntr(2)), 1)
!
         else if ( type .eq. 2 ) then
!
            if (sigmai .eq. zero) then
!
!              %----------------------------------%
!              | Shift is real.  Perform          |
!              | y <--- OP*x = inv[A-sigmar*I]*x  |
!              | to force the starting vector     |
!              | into the range of OP.            |
!              %----------------------------------%
!
               call scopy (n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
               call sgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
     &                       iwork, workd(ipntr(2)), n, ierr)
               if ( ierr .ne. 0 ) then
                  print*, ' '
                  print*, ' _NBAND: Error with _bgtrs. '
                  print*, ' '
                  go to 9000
               end if
!
            else
!
!              %--------------------------------------------%
!              | Shift is COMPLEX. Perform                  |
!              | y <--- OP*x = Real_Part{inv[A-sigma*I]*x}  |
!              | to force the starting vector into the      |
!              | range of OP.                               |
!              %--------------------------------------------%
!
               do 100 j = 1, n
                  workc(j) = cmplx(workd(ipntr(1)+j-1))
  100          continue
!
               call cgbtrs ('Notranspose', n, kl, ku, 1, cfac, lda,
     &                       iwork, workc, n, ierr)
               if ( ierr .ne. 0 ) then
                  print*, ' '
                  print*, ' _NBAND: Error with _gbtrs. '
                  print*, ' '
                  go to 9000
               end if
!
               do 110 j = 1, n
                  workd(ipntr(2)+j-1) = real(workc(j))
  110          continue
!
            end if
!
         else if ( type .eq. 3 ) then
!
!           %-----------------------------------%
!           | Perform  y <--- OP*x = inv[M]*A*x |
!           | to force the starting vector into |
!           | the range of OP.                  |
!           %-----------------------------------%
!
            call sgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1),
     &                  lda, workd(ipntr(1)), 1, zero,
     &                  workd(ipntr(2)), 1)
!
            call sgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
     &                    iwork, workd(ipntr(2)), n, ierr)
            if (ierr .ne. 0) then
               print*, ' '
               print*, '_NBAND: Error with _bgtrs.'
               print*, ' '
               go to 9000
            end if
!
         else if ( type .eq. 4 ) then
!
!           %-----------------------------------------%
!           | Perform y <-- OP*x                      |
!           |         = Real_part{inv[A-SIGMA*M]*M}*x |
!           | to force the starting vector into the   |
!           | range of OP.                            |
!           %-----------------------------------------%
!
            call sgbmv('Notranspose', n, n, kl, ku, one, mb(itop,1),
     &                 lda, workd(ipntr(1)), 1, zero,
     &                 workd(ipntr(2)), 1)
!
            if ( sigmai .eq. zero ) then
!
!              %---------------------%
!              | Shift is real, stay |
!              | in real arithmetic. |
!              %---------------------%
!
               call sgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
     &                      iwork, workd(ipntr(2)), n, ierr)
               if (ierr .ne. 0) then
                  print*, ' '
                  print*, '_NBAND: Error with _gbtrs.'
                  print*, ' '
                  go to 9000
               end if
!
            else
!
!              %--------------------------%
!              | Goto complex arithmetic. |
!              %--------------------------%
!
               do 120 i = 1,n
                  workc(i) = cmplx(workd(ipntr(2)+i-1))
  120           continue
!
               call cgbtrs ('Notranspose', n, kl, ku, 1, cfac, lda,
     &                      iwork, workc, n, ierr)
               if (ierr .ne. 0) then
                  print*, ' '
                  print*, '_NBAND: Error with _gbtrs.'
                  print*, ' '
                  go to 9000
               end if
!
               do  130 i = 1, n
                  workd(ipntr(2)+i-1) = real(workc(i))
  130          continue
!
            end if
!
         else if ( type .eq. 5) then
!
!           %---------------------------------------%
!           | Perform y <-- OP*x                    |
!           |    = Imaginary_part{inv[A-SIGMA*I]}*x |
!           | to force the starting vector into the |
!           | range of OP.                          |
!           %---------------------------------------%
!
            do 140 j = 1, n
                  workc(j) = cmplx(workd(ipntr(1)+j-1))
  140       continue
!
            call cgbtrs ('Notranspose', n, kl, ku, 1, cfac, lda,
     &                    iwork, workc, n, ierr)
            if ( ierr .ne. 0 ) then
               print*, ' '
               print*, ' _NBAND: Error with _gbtrs. '
               print*, ' '
               go to 9000
            end if
!
            do 150 j = 1, n
               workd(ipntr(2)+j-1) = aimag(workc(j))
  150       continue
!
         else if ( type .eq. 6 ) then
!
!           %----------------------------------------%
!           | Perform y <-- OP*x                     |
!           |       Imaginary_part{inv[A-SIGMA*M]*M} |
!           | to force the starting vector into the  |
!           | range of OP.                           |
!           %----------------------------------------%
!
            call sgbmv('Notranspose', n, n, kl, ku, one, mb(itop,1),
     &                 lda, workd(ipntr(1)), 1, zero,
     &                 workd(ipntr(2)), 1)
!
            do 160 i = 1,n
               workc(i) = cmplx(workd(ipntr(2)+i-1))
  160       continue
!
            call cgbtrs ('Notranspose', n, kl, ku, 1, cfac, lda,
     &                   iwork, workc, n, ierr)
            if (ierr .ne. 0) then
               print*, ' '
               print*, '_NBAND: Error with _gbtrs.'
               print*, ' '
               go to 9000
            end if
!
            do  170 i = 1, n
               workd(ipntr(2)+i-1) = aimag(workc(i))
  170       continue
!
         end if
!
      else if (ido .eq. 1) then
!
         if ( type .eq. 1) then
!
!           %----------------------------%
!           | Perform  y <--- OP*x = A*x |
!           %----------------------------%
!
            call sgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1),
     &                 lda, workd(ipntr(1)), 1, zero,
     &                 workd(ipntr(2)), 1)
!
         else if ( type .eq. 2) then
!
            if ( sigmai .eq. zero) then
!
!              %----------------------------------%
!              | Shift is real.  Perform          |
!              | y <--- OP*x = inv[A-sigmar*I]*x. |
!              %----------------------------------%
!
               call scopy (n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
               call sgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
     &                       iwork, workd(ipntr(2)), n, ierr)
            else
!
!              %------------------------------------------%
!              | Shift is COMPLEX. Perform                |
!              | y <-- OP*x = Real_Part{inv[A-sigma*I]*x} |
!              | in COMPLEX arithmetic.                   |
!              %------------------------------------------%
!
               do 180 j = 1, n
                  workc(j) = cmplx(workd(ipntr(1)+j-1))
  180          continue
!
               call cgbtrs ('Notranspose', n, kl, ku, 1, cfac, lda,
     &                       iwork, workc, n, ierr)
               if ( ierr .ne. 0 ) then
                  print*, ' '
                  print*, '_NBAND: Error with _gbtrs.'
                  print*, ' '
                  go to 9000
               end if
!
               do 190 j = 1, n
                  workd(ipntr(2)+j-1) = real(workc(j))
  190          continue
!
            end if
!
         else if ( type .eq. 3 ) then
!
!           %-----------------------------------%
!           | Perform  y <--- OP*x = inv[M]*A*x |
!           %-----------------------------------%
!
            call sgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1),
     &                  lda, workd(ipntr(1)), 1, zero,
     &                  workd(ipntr(2)), 1)
!
            call sgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
     &                    iwork, workd(ipntr(2)), n, ierr)
            if (ierr .ne. 0) then
               print*, ' '
               print*, '_NBAND: Error with _bgtrs.'
               print*, ' '
               go to 9000
            end if
!
         else if ( type .eq. 4 ) then
!
!           %--------------------------------------%
!           | Perform  y <-- inv(A-sigma*M)*(M*x). |
!           | (M*x) has been computed and stored   |
!           | in workd(ipntr(3)).                  |
!           %--------------------------------------%
!
            if ( sigmai .eq. zero ) then
!
!              %------------------------%
!              | Shift is real, stay in |
!              | real arithmetic.       |
!              %------------------------%
!
               call scopy(n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)
               call sgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
     &                       iwork, workd(ipntr(2)), n, ierr)
               if (ierr .ne. 0) then
                  print*, ' '
                  print*, '_NBAND: Error with _gbtrs.'
                  print*, ' '
                  go to 9000
               end if
!
            else
!
!              %---------------------------%
!              | Go to COMPLEX arithmetic. |
!              %---------------------------%
!
               do 200 i = 1,n
                  workc(i) = cmplx(workd(ipntr(3)+i-1))
  200          continue
!
               call cgbtrs ('Notranspose', n, kl, ku, 1, cfac, lda,
     &                       iwork, workc, n, ierr)
               if (ierr .ne. 0) then
                  print*, ' '
                  print*, '_NBAND: Error in _gbtrs.'
                  print*, ' '
                  go to 9000
               end if
!
               do 210 i = 1,n
                  workd(ipntr(2)+i-1) = real(workc(i))
  210          continue
!
            end if
!
         else if ( type .eq. 5 ) then
!
!           %---------------------------------------%
!           | Perform y <-- OP*x                    |
!           |    = Imaginary_part{inv[A-SIGMA*I]*x} |
!           %---------------------------------------%
!
            do 220 j = 1, n
                  workc(j) = cmplx(workd(ipntr(1)+j-1))
  220       continue
!
            call cgbtrs ('Notranspose', n, kl, ku, 1, cfac, lda,
     &                    iwork, workc, n, ierr)
            if ( ierr .ne. 0 ) then
               print*, ' '
               print*, ' _NBAND: Error with _gbtrs. '
               print*, ' '
               go to 9000
            end if
!
            do 230 j = 1, n
               workd(ipntr(2)+j-1) = aimag(workc(j))
  230       continue
!
         else if ( type .eq. 6) then
!
!           %-----------------------------------------%
!           | Perform y <-- OP*x                      |
!           |   = Imaginary_part{inv[A-SIGMA*M]*M}*x. |
!           %-----------------------------------------%
!
            do 240 i = 1,n
               workc(i) = cmplx(workd(ipntr(3)+i-1))
  240       continue
!
            call cgbtrs ('Notranspose', n, kl, ku, 1, cfac, lda,
     &                   iwork, workc, n, ierr)
            if (ierr .ne. 0) then
               print*, ' '
               print*, '_NBAND: Error with _gbtrs.'
               print*, ' '
               go to 9000
            end if
!
            do  250 i = 1, n
               workd(ipntr(2)+i-1) = aimag(workc(i))
  250       continue
!
         end if
!
      else if (ido .eq. 2) then
!
!        %--------------------%
!        | Perform y <-- M*x  |
!        | Not used when      |
!        | type = 1,2.        |
!        %--------------------%
!
          call sgbmv('Notranspose', n, n, kl, ku, one, mb(itop,1),
     &                lda, workd(ipntr(1)), 1, zero,
     &                workd(ipntr(2)), 1)
!
      else
!
!        %-----------------------------------------%
!        | Either we have convergence, or there is |
!        | error.                                  |
!        %-----------------------------------------%
!
         if ( info .lt. 0) then
!
!           %--------------------------%
!           | Error message, check the |
!           | documentation in SNAUPD  |
!           %--------------------------%
!
            print *, ' '
            print *, ' Error with _naupd info = ',info
            print *, ' Check the documentation of _naupd '
            print *, ' '
            go to 9000
!
         else
!
            if ( info .eq. 1) then
               print *, ' '
               print *, ' Maximum number of iterations reached.'
               print *, ' '
            else if ( info .eq. 3) then
               print *, ' '
               print *, ' No shifts could be applied during implicit',
     &                  ' Arnoldi update, try increasing NCV.'
               print *, ' '
            end if
!
            if (iparam(5) .gt. 0) then
!
               call sneupd ( rvec, 'A', select, dr, di, z, ldz,
     &                 sigmar, sigmai, workev, bmat, n, which,
     &                 nev, tol, resid, ncv, v, ldv, iparam,
     &                 ipntr, workd, workl, lworkl, info )
!
               if ( info .ne. 0) then
!
!                 %------------------------------------%
!                 | Check the documentation of SNEUPD. |
!                 %------------------------------------%
!
                  print *, ' '
                  print *, ' Error with _neupd = ', info
                  print *, ' Check the documentation of _neupd '
                  print *, ' '
                  go to 9000
!
               else if ( sigmai .ne. zero ) then
!
                  if ( type .eq. 4 .or. type .eq. 6 ) then
!
                     first = .true.
                     do 270 j = 1, iparam(5)
!
!                    %----------------------------------%
!                    | Use Rayleigh Quotient to recover |
!                    | eigenvalues of the original      |
!                    | generalized eigenvalue problem.  |
!                    %----------------------------------%
!
                     if ( di(j) .eq. zero ) then
!
!                       %--------------------------------------%
!                       | Eigenvalue is real. Compute          |
!                       | d = (x'*inv[A-sigma*M]*M*x) / (x'*x) |
!                       %--------------------------------------%
!
                        call sgbmv('Nontranspose', n, n, kl, ku, one,
     $                     mb(itop,1), lda, z(1,j), 1, zero,
     $                     workd, 1)
                        do i = 1, n
                           workc(i) = cmplx(workd(i))
                        end do
                        call cgbtrs ('Notranspose', n, kl, ku, 1,
     $                        cfac, lda, iwork, workc, n, info)
                        do i = 1, n
                           workd(i) = real(workc(i))
                           workd(i+n) = aimag(workc(i))
                        end do
                        denr = sdot(n, z(1,j), 1, workd, 1)
                        deni = sdot(n, z(1,j), 1, workd(n+1), 1)
                        numr  = snrm2(n, z(1,j), 1)**2
                        dmdul = slapy2(denr,deni)**2
                        if ( dmdul .ge. safmin ) then
                           dr(j) = sigmar + numr*denr / dmdul
                        else
!
!                          %---------------------%
!                          | dmdul is too small. |
!                          | Exit to avoid       |
!                          | overflow.           |
!                          %---------------------%
!
                           info = -15
                           go to 9000
                        end if
!
                     else if (first) then
!
!                       %------------------------%
!                       | Eigenvalue is complex. |
!                       | Compute the first one  |
!                       | of the conjugate pair. |
!                       %------------------------%
!
!                       %-------------%
!                       | Compute M*x |
!                       %-------------%
!
                        call sgbmv('Nontranspose', n, n, kl, ku,
     $                      one, mb(itop,1), lda, z(1,j), 1, zero,
     $                      workd, 1)
                        call sgbmv('Nontranspose', n, n, kl, ku,
     $                       one, mb(itop,1), lda, z(1,j+1), 1,
     $                       zero, workd(n+1), 1)
                        do i = 1, n
                           workc(i) = cmplx(workd(i),workd(i+n))
                        end do
!
!                       %----------------------------%
!                       | Compute inv(A-sigma*M)*M*x |
!                       %----------------------------%
!
                        call cgbtrs('Notranspose',n,kl,ku,1,cfac,
     $                     lda, iwork, workc, n, info)
!
!                       %-------------------------------%
!                       | Compute x'*inv(A-sigma*M)*M*x |
!                       %-------------------------------%
!
                        do i = 1, n
                           workd(i) = real(workc(i))
                           workd(i+n) = aimag(workc(i))
                        end do
                        denr = sdot(n,z(1,j),1,workd,1)
                        denr = denr+sdot(n,z(1,j+1),1,workd(n+1),1)
                        deni = sdot(n,z(1,j),1,workd(n+1),1)
                        deni = deni - sdot(n,z(1,j+1),1,workd,1)
!
!                       %----------------%
!                       | Compute (x'*x) |
!                       %----------------%
!
                        numr = slapy2( snrm2(n, z(1,j), 1),
     &                         snrm2(n, z(1, j+1), 1) )**2
!
!                       %----------------------------------------%
!                       | Compute (x'x) / (x'*inv(A-sigma*M)*Mx) |
!                       %----------------------------------------%
!
                        dmdul = slapy2(denr,deni)**2
                        if ( dmdul .ge. safmin ) then
                           dr(j) = sigmar+numr*denr / dmdul
                           di(j) = sigmai-numr*deni / dmdul
                           first = .false.
                        else
!
!                          %---------------------%
!                          | dmdul is too small. |
!                          | Exit to avoid       |
!                          | overflow.           |
!                          %---------------------%
!
                           info = -15
                           go to 9000
!
                        end if
!
                     else
!
!                       %---------------------------%
!                       | Get the second eigenvalue |
!                       | of the conjugate pair by  |
!                       | taking the conjugate of   |
!                       | previous one.             |
!                       %---------------------------%
!
                        dr(j) = dr(j-1)
                        di(j) = -di(j-1)
                        first = .true.
!
                     end if
!
  270                continue
!
                  else if ( type .eq. 2 .or. type .eq. 5) then
!
                     first = .true.
                     do 280 j = 1, iparam(5)
!
!                    %----------------------------------%
!                    | Use Rayleigh Quotient to recover |
!                    | eigenvalues of the original      |
!                    | standard eigenvalue problem.     |
!                    %----------------------------------%
!
                     if ( di(j) .eq. zero ) then
!
!                       %-------------------------------------%
!                       | Eigenvalue is real. Compute         |
!                       | d = (x'*inv[A-sigma*I]*x) / (x'*x). |
!                       %-------------------------------------%
!
                        do i = 1, n
                           workc(i) = cmplx(z(i,j))
                        end do
                        call cgbtrs ('Notranspose', n, kl, ku, 1,
     $                        cfac, lda, iwork, workc, n, info)
                        do i = 1, n
                           workd(i) = real(workc(i))
                           workd(i+n) = aimag(workc(i))
                        end do
                        denr = sdot(n,z(1,j),1,workd,1)
                        deni = sdot(n,z(1,j),1,workd(n+1),1)
                        numr  = snrm2(n, z(1,j), 1)**2
                        dmdul = slapy2(denr,deni)**2
                        if ( dmdul .ge. safmin ) then
                           dr(j) = sigmar + numr*denr / dmdul
                        else
!
!                          %---------------------%
!                          | dmdul is too small. |
!                          | Exit to avoid       |
!                          | overflow.           |
!                          %---------------------%
!
                           info = -15
                           go to 9000
!
                        end if
!
                     else if (first) then
!
!                       %------------------------%
!                       | Eigenvalue is complex. |
!                       | Compute the first one  |
!                       | of the conjugate pair. |
!                       %------------------------%
!
                        do i = 1, n
                           workc(i) = cmplx( z(i,j), z(i,j+1) )
                        end do
!
!                       %---------------------------%
!                       | Compute inv[A-sigma*I]*x. |
!                       %---------------------------%
!
                        call cgbtrs('Notranspose',n,kl,ku,1,cfac,
     $                       lda, iwork, workc, n, info)
!
!                       %-----------------------------%
!                       | Compute x'*inv(A-sigma*I)*x |
!                       %-----------------------------%
!
                        do i = 1, n
                           workd(i) = real(workc(i))
                           workd(i+n) = aimag(workc(i))
                        end do
                        denr = sdot(n,z(1,j),1,workd,1)
                        denr = denr+sdot(n,z(1,j+1),1,workd(n+1),1)
                        deni = sdot(n,z(1,j),1,workd(n+1),1)
                        deni = deni - sdot(n,z(1,j+1),1,workd,1)
!
!                       %----------------%
!                       | Compute (x'*x) |
!                       %----------------%
!
                        numr = slapy2( snrm2(n, z(1,j), 1),
     &                         snrm2(n, z(1,j+1), 1))**2
!
!                       %----------------------------------------%
!                       | Compute (x'x) / (x'*inv(A-sigma*I)*x). |
!                       %----------------------------------------%
!
                        dmdul = slapy2(denr,deni)**2
                        if (dmdul .ge. safmin) then
                           dr(j) = sigmar+numr*denr / dmdul
                           di(j) = sigmai-numr*deni / dmdul
                           first = .false.
                        else
!
!                          %---------------------%
!                          | dmdul is too small. |
!                          | Exit to avoid       |
!                          | overflow.           |
!                          %---------------------%
!
                           info = -15
                           go to 9000
                        end if
!
                     else
!
!                       %---------------------------%
!                       | Get the second eigenvalue |
!                       | of the conjugate pair by  |
!                       | taking the conjugate of   |
!                       | previous one.             |
!                       %---------------------------%
!
                        dr(j) = dr(j-1)
                        di(j) = -di(j-1)
                        first = .true.
!
                     end if
!
  280                continue
!
                  end if
!
               end if
!
            end if
!
         end if
!
         go to 9000
!
      end if
!
!     %----------------------------------------%
!     | L O O P  B A C K to call SNAUPD again. |
!     %----------------------------------------%
!
      go to 90
!
 9000 continue
!
      end
