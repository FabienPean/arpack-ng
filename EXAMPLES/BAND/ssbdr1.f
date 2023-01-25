      program ssbdr1
!
!     ... Construct the matrix A in LAPACK-style band form.
!         The matrix A is derived from the discretization of
!         the 2-dimensional Laplacian on the unit square with
!         zero Dirichlet boundary condition using standard
!         central difference.
!
!     ... Call SSBAND to find eigenvalues LAMBDA such that
!                          A*x = x*LAMBDA.
!
!     ... Use mode 1 of SSAUPD.
!
!\BeginLib
!
!\Routines called:
!     ssband  ARPACK banded eigenproblem solver.
!     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     slaset  LAPACK routine to initialize a matrix to zero.
!     saxpy   Level 1 BLAS that computes y <- alpha*x+y.
!     snrm2   Level 1 BLAS that computes the norm of a vector.
!     sgbmv   Level 2 BLAS that computes the band matrix vector product
!
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
! FILE: sbdr1.F   SID: 2.5   DATE OF SID: 08/26/96   RELEASE: 2
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
      Real
     &                 a(lda,maxn), m(lda,maxn), rfac(lda,maxn),
     &                 workl(maxncv*maxncv+8*maxncv), workd(3*maxn),
     &                 v(ldv, maxncv), resid(maxn), d(maxncv, 2),
     &                 ax(maxn)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character        which*2, bmat
      integer          nev, ncv, ku, kl, info, i, j, ido,
     &                 n, nx, lo, isub, isup, idiag, maxitr, mode,
     &                 nconv
      Real
     &                 tol, sigma, h2
      logical          rvec
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real
     &                 one, zero, two
      parameter        (one = 1.0E+0 , zero = 0.0E+0 , two = 2.0E+0 )
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      Real
     &                  slapy2, snrm2
      external          slapy2, snrm2, sgbmv, saxpy
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
!     | The number NX is the number of interior points  |
!     | in the discretization of the 2-dimensional      |
!     | Laplacian operator on the unit square with zero |
!     | Dirichlet boundary condition. The number        |
!     | N(=NX*NX) is the dimension of the matrix.  A    |
!     | standard eigenvalue problem is solved           |
!     | (BMAT = 'I').  NEV is the number of eigenvalues |
!     | to be approximated. The user can modify NX,NEV, |
!     | NCV and WHICH to solve problems of different    |
!     | sizes, and to get different parts the spectrum. |
!     | However, the following conditions must be       |
!     | satisfied:                                      |
!     |                   N <= MAXN                     |
!     |                 NEV <= MAXNEV                   |
!     |           NEV + 1 <= NCV <= MAXNCV              |
!     %-------------------------------------------------%
!
      nx  = 10
      n    = nx*nx
      nev  = 4
      ncv  = 10
      if ( n .gt. maxn ) then
         print *, ' ERROR with _SBDR1: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _SBDR1: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _SBDR1: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat = 'I'
      which = 'LM'
!
!     %-----------------------------------------------------%
!     | The work array WORKL is used in SSAUPD as           |
!     | workspace.  Its dimension LWORKL is set as          |
!     | illustrated below.  The parameter TOL determines    |
!     | the stopping criterion. If TOL<=0, machine          |
!     | precision is used.  The variable IDO is used for    |
!     | reverse communication, and is initially set to 0.   |
!     | Setting INFO=0 indicates that a random vector is    |
!     | generated in SSAUPD to start the Arnoldi iteration. |
!     %-----------------------------------------------------%
!
      lworkl  = ncv**2+8*ncv
      tol  = zero
      ido  = 0
      info = 0
!
!     %---------------------------------------------------%
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of SSAUPD is used     |
!     | (IPARAM(7) = 1). All these options can be changed |
!     | by the user. For details see the documentation in |
!     | SSBAND.                                           |
!     %---------------------------------------------------%
!
      maxitr  = 300
      mode = 1
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
      call slaset('A', lda, n, zero, zero, a, lda)
      call slaset('A', lda, n, zero, zero, m, lda)
      call slaset('A', lda, n, zero, zero, rfac, lda)
!
!     %-------------------------------------%
!     | KU, KL are number of superdiagonals |
!     | and subdiagonals within the band of |
!     | matrices A and M.                   |
!     %-------------------------------------%
!
      kl   = nx
      ku   = nx
!
!     %---------------%
!     | Main diagonal |
!     %---------------%
!
      h2 = one / ((nx+1)*(nx+1))
      idiag = kl+ku+1
      do 30 j = 1, n
         a(idiag,j) = 4.0E+0  / h2
  30  continue
!
!     %-------------------------------------%
!     | First subdiagonal and superdiagonal |
!     %-------------------------------------%
!
      isup = kl+ku
      isub = kl+ku+2
      do 50 i = 1, nx
        lo = (i-1)*nx
        do 40 j = lo+1, lo+nx-1
           a(isup,j+1) = -one / h2
           a(isub,j) = -one / h2
  40    continue
  50  continue
!
!     %------------------------------------%
!     | KL-th subdiagonal and KU-th super- |
!     | diagonal.                          |
!     %------------------------------------%
!
      isup = kl+1
      isub = 2*kl+ku+1
      do 80 i = 1, nx-1
         lo = (i-1)*nx
         do 70 j = lo+1, lo+nx
            a(isup,nx+j)  = -one / h2
            a(isub,j) = -one / h2
 70      continue
 80   continue
!
!     %-------------------------------------%
!     | Call SSBAND to find eigenvalues and |
!     | eigenvectors.  Eigenvalues are      |
!     | returned in the first column of D.  |
!     | Eigenvectors are returned in the    |
!     | first NCONV (=IPARAM(5)) columns of |
!     | V.                                  |
!     %-------------------------------------%
!
      rvec = .true.
      call ssband( rvec, 'A', select, d, v, ldv, sigma, n, a, m, lda,
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
         print *, ' _SBDR1 '
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
            call sgbmv('Notranspose', n, n, kl, ku, one,
     &                 a(kl+1,1), lda, v(1,j), 1, zero,
     &                 ax, 1)
            call saxpy(n, -d(j,1), v(1,j), 1, ax, 1)
            d(j,2) = snrm2(n, ax, 1)
            d(j,2) = d(j,2) / abs(d(j,1))
!
 90      continue

         call smout(6, nconv, 2, d, maxncv, -6,
     &             'Ritz values and relative residuals')
      else
!
!        %-------------------------------------%
!        | Either convergence failed, or there |
!        | is error.  Check the documentation  |
!        | for SSBAND.                         |
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
