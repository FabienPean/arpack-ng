      program dsbdr3
!
!     ... Construct the matrix A in LAPACK-style band form.
!         The matrix A is the 1-dimensional discrete Laplacian on [0,1]
!         with zero Dirichlet boundary condition, M is the mass
!         formed by using piecewise linear elements on [0,1].
!
!     ... Call DSBAND  with regular mode to find eigenvalues LAMBDA
!         such that
!                          A*x = LAMBDA*M*x.
!
!     ... Use mode 2 of DSAUPD .
!
!\BeginLib
!
!\Routines called:
!     dsband   ARPACK banded eigenproblem solver.
!     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     dlaset   LAPACK routine to initialize a matrix to zero.
!     daxpy    Level 1 BLAS that computes y <- alpha*x+y.
!     dnrm2    Level 1 BLAS that computes the norm of a vector.
!     dgbmv    Level 2 BLAS that computes the band matrix vector product
!\Author
!     Richard Lehoucq
!     Danny Sorensen
!     Chao Yang
!     Dept. of Computational &
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\SCCS Information: @(#)
! FILE: sbdr3.F   SID: 2.5   DATE OF SID: 08/26/96   RELEASE: 2
!
!\Remarks
!     1. None
!
!\EndLib
!
!----------------------------------------------------------------------
!
!     %-------------------------------------%
!     | Define leading dimensions for all   |
!     | arrays.                             |
!     | MAXN   - Maximum size of the matrix |
!     | MAXNEV - Maximum number of          |
!     |          eigenvalues to be computed |
!     | MAXNCV - Maximum number of Arnoldi  |
!     |          vectors stored             |
!     | MAXBDW - Maximum bandwidth          |
!     %-------------------------------------%
!
      integer          maxn, maxnev, maxncv, maxbdw, lda,
     &                 lworkl, ldv
      parameter        ( maxn = 1000, maxnev = 25, maxncv=50,
     &                   maxbdw=50, lda = maxbdw, ldv = maxn )
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
      integer          iparam(11), iwork(maxn)
      logical          select(maxncv)
      Double precision
     &                 a(lda,maxn), m(lda,maxn), rfac(lda,maxn),
     &                 workl(maxncv*maxncv+8*maxncv), workd(3*maxn),
     &                 v(ldv, maxncv), resid(maxn), d(maxncv, 2),
     &                 ax(maxn), mx(maxn)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character        which*2, bmat
      integer          nev, ncv, ku, kl, info, j, ido,
     &                 n, isub, isup, idiag, maxitr, mode, nconv
      Double precision
     &                 tol, h, sigma, r1, r2
      logical          rvec
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Double precision
     &                 one, zero, two, four, six
      parameter        (one = 1.0D+0 , zero = 0.0D+0 , two = 2.0D+0 ,
     &                  four = 4.0D+0 , six = 6.0D+0 )
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      Double precision
     &                  dlapy2 , dnrm2
      external          dlapy2 , dnrm2 , daxpy , dgbmv
!
!     %--------------------%
!     | Intrinsic function |
!     %--------------------%
!
      intrinsic         abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %-------------------------------------------------%
!     | The number N is the dimension of the matrix.  A |
!     | generalized eigenvalue problem is solved        |
!     | (BMAT = 'G').  NEV is the number of eigenvalues |
!     | to be approximated. The user can modify N, NEV, |
!     | NCV and WHICH to solve problems of different    |
!     | sizes, and to get different parts the spectrum. |
!     | However, the following conditions must be       |
!     | satisfied:                                      |
!     |                   N <= MAXN                     |
!     |                 NEV <= MAXNEV                   |
!     |           NEV + 1 <= NCV <= MAXNCV              |
!     %-------------------------------------------------%
!
      n    = 100
      nev  = 4
      ncv  = 10
      if ( n .gt. maxn ) then
         print *, ' ERROR with _SBDR3: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _SBDR3: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _SBDR3: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat = 'G'
      which = 'LM'
!
!     %-----------------------------------------------------%
!     | The work array WORKL is used in DSAUPD  as           |
!     | workspace.  Its dimension LWORKL is set as          |
!     | illustrated below.  The parameter TOL determines    |
!     | the stopping criterion. If TOL<=0, machine          |
!     | precision is used.  The variable IDO is used for    |
!     | reverse communication, and is initially set to 0.   |
!     | Setting INFO=0 indicates that a random vector is    |
!     | generated in DSAUPD  to start the Arnoldi iteration. |
!     %-----------------------------------------------------%
!
      lworkl  = ncv**2+8*ncv
      tol  = zero
      ido  = 0
      info = 0
!
!     %---------------------------------------------------%
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 2 of DSAUPD  is used     |
!     | (IPARAM(7) = 2). All these options can be changed |
!     | by the user. For details see the documentation in |
!     | DSBAND .                                            |
!     %---------------------------------------------------%
!
      maxitr = 300
      mode   = 2
!
      iparam(3) = maxitr
      iparam(7) = mode
!
!     %----------------------------------------%
!     | Construct the matrix A in LAPACK-style |
!     | banded form.                           |
!     %----------------------------------------%
!
!     %---------------------------------------------%
!     | Zero out the workspace for banded matrices. |
!     %---------------------------------------------%
!
      call dlaset ('A', lda, n, zero, zero, a, lda)
      call dlaset ('A', lda, n, zero, zero, m, lda)
      call dlaset ('A', lda, n, zero, zero, rfac, lda)
!
!     %-------------------------------------%
!     | KU, KL are number of superdiagonals |
!     | and subdiagonals within the band of |
!     | matrices A and M.                   |
!     %-------------------------------------%
!
      kl   = 1
      ku   = 1
!
!     %---------------%
!     | Main diagonal |
!     %---------------%
!
      h = one / dble (n+1)
      r1 = four / six
      idiag = kl+ku+1
      do 30 j = 1, n
         a(idiag,j) = two / h
         m(idiag,j) = r1 * h
  30  continue
!
!     %-------------------------------------%
!     | First subdiagonal and superdiagonal |
!     %-------------------------------------%
!
      r2 = one / six
      isup = kl+ku
      isub = kl+ku+2
      do 60 j = 1, n-1
         a(isup,j+1) = -one / h
         a(isub,j) = -one / h
         m(isup,j+1) = r2 * h
         m(isub,j) = r2 * h
  60  continue
!
!     %-------------------------------------%
!     | Call DSBAND  to find eigenvalues and |
!     | eigenvectors.  Eigenvalues are      |
!     | returned in the first column of D.  |
!     | Eigenvectors are returned in the    |
!     | first NCONV (=IPARAM(5)) columns of |
!     | V.                                  |
!     %-------------------------------------%
!
      rvec = .true.
      call dsband ( rvec, 'A', select, d, v, ldv, sigma, n, a, m, lda,
     &             rfac, kl, ku, which, bmat, nev, tol,
     &             resid, ncv, v, ldv, iparam, workd, workl, lworkl,
     &             iwork, info)
!
      if ( info .eq. 0) then
!
         nconv = iparam(5)
!
!        %-----------------------------------%
!        | Print out convergence information |
!        %-----------------------------------%
!
         print *, ' '
         print *, ' _SBDR3 '
         print *, ' ====== '
         print *, ' '
         print *, ' The size of the matrix is ', n
         print *, ' Number of eigenvalue requested is ', nev
         print *, ' The number of Lanczos vectors generated',
     &            ' (NCV) is ', ncv
         print *, ' The number of converged Ritz values is ',
     &              nconv
         print *, ' What portion of the spectrum ', which
         print *, ' The number of Implicit Arnoldi',
     &              ' update taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence tolerance is ', tol
         print *, ' '
!
!        %----------------------------%
!        | Compute the residual norm. |
!        |    ||  A*x - lambda*x ||   |
!        %----------------------------%
!
         do 90 j = 1, nconv
            call dgbmv ('Notranspose', n, n, kl, ku, one,
     &                 a(kl+1,1), lda, v(1,j), 1, zero,
     &                 ax, 1)
            call dgbmv ('Notranspose', n, n, kl, ku, one,
     &                 m(kl+1,1), lda, v(1,j), 1, zero,
     &                 mx, 1)
            call daxpy (n, -d(j,1), mx, 1, ax, 1)
            d(j,2) = dnrm2 (n, ax, 1)
            d(j,2) = d(j,2) / abs(d(j,1))
!
 90      continue

         call dmout (6, nconv, 2, d, maxncv, -6,
     &             'Ritz values and relative residuals')
      else
!
!        %-------------------------------------%
!        | Either convergence failed, or there |
!        | is error.  Check the documentation  |
!        | for DSBAND .                         |
!        %-------------------------------------%
!
          print *, ' '
          print *, ' Error with _sband, info= ', info
          print *, ' Check the documentation of _sband '
          print *, ' '
!
      end if
!
 9000 end
