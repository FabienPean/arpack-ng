      program snbdr5
!
!     ... Construct matrices A and M in LAPACK-style band form.
!         The matrix A is a block tridiagonal matrix.  Each
!         diagonal block is a tridiagonal matrix with
!         4 on the diagonal, 1-rho*h/2 on the subdiagonal and
!         1+rho*h/2 on the superdiagonal.  Each off-diagonal block
!         of A is an identity matrices.
!
!     ... Define COMPLEX shift SIGMA = (SIGMAR,SIGMAI), SIGMAI .ne. 0.
!
!     ... Call SNBAND to find eigenvalues LAMBDA closest to SIGMA
!         such that
!                       A*x = LAMBDA*x.
!
!     ... Use mode 4 of SNAUPD.
!
!\BeginLib
!
!\Routines called:
!     snband  ARPACK banded eigenproblem solver.
!     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     slaset  LAPACK routine to initialize a matrix to zero.
!     saxpy   Level 1 BLAS that computes y <- alpha*x+y.
!     snrm2   Level 1 BLAS that computes the norm of a vector.
!     sgbmv   Level 2 BLAS that computes the band matrix vector product.
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
! FILE: nbdr5.F   SID: 2.5   DATE OF SID: 08/26/96   RELEASE: 2
!
!\Remarks
!     1. None
!
!\EndLib
!
!-------------------------------------------------------------------------
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
     &                 workl(3*maxncv*maxncv+6*maxncv), workd(3*maxn),
     &                 workev(3*maxncv), v(ldv, maxncv),
     &                 resid(maxn), d(maxncv, 3), ax(maxn)
      Complex
     &                 cfac(lda, maxn), workc(maxn)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character        which*2, bmat
      integer          nev, ncv, ku, kl, info, i, j, ido,
     &                 n, nx, lo, idiag, isup, isub, mode, maxitr,
     &                 nconv
      logical          rvec, first
      Real
     &                 tol, rho, h, sigmar, sigmai
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real
     &                 one, zero, two
      parameter        (one = 1.0E+0 , zero = 0.0E+0 ,
     &                  two = 2.0E+0 )
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
!     %--------------------------------------------------%
!     | The number NX is the size of each block diagonal |
!     | of A. The number N(=NX*NX) is the dimension of   |
!     | the matrix.  A standard eigenvalue problem is    |
!     | solved (BMAT = 'I').  NEV is the number of       |
!     | eigenvalues (closest to (SIGMAR,SIGMAI)) to be   |
!     | approximated. Since the shift-invert moded is    |
!     | used, WHICH is set to 'LM'. The user can modify  |
!     | NX, NEV, NCV, SIGMAR, SIGMAI to solve problems   |
!     | of different sizes, and to get different parts   |
!     | the spectrum. However, The following conditions  |
!     | must be satisfied:                               |
!     |                   N <= MAXN                      |
!     |                 NEV <= MAXNEV                    |
!     |           NEV + 2 <= NCV <= MAXNCV               |
!     %--------------------------------------------------%
!
      nx   = 10
      n    = nx*nx
      nev  = 4
      ncv  = 10
      if ( n .gt. maxn ) then
         print *, ' ERROR with _NBDR5: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NBDR5: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NBDR5: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat = 'I'
      which = 'LM'
      sigmar = 4.0E-1
      sigmai = 6.0E-1
!
!     %-----------------------------------------------------%
!     | The work array WORKL is used in SNAUPD as           |
!     | workspace.  Its dimension LWORKL is set as          |
!     | illustrated below.  The parameter TOL determines    |
!     | the stopping criterion. If TOL<=0, machine          |
!     | precision is used.  The variable IDO is used for    |
!     | reverse communication, and is initially set to 0.   |
!     | Setting INFO=0 indicates that a random vector is    |
!     | generated in SNAUPD to start the Arnoldi iteration. |
!     %-----------------------------------------------------%
!
      lworkl  = 3*ncv**2+6*ncv
      tol  = zero
      ido  = 0
      info = 0
!
!     %---------------------------------------------------%
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 4 of SNAUPD is used     |
!     | (IPARAM(7) = 4). All these options can be changed |
!     | by the user. For details, see the documentation   |
!     | in SNBAND.                                        |
!     %---------------------------------------------------%
!
      maxitr = 300
      mode = 4
!
      iparam(3) = maxitr
      iparam(7) = mode
!
!     %--------------------------------------------%
!     | Construct matrices A and M in LAPACK-style |
!     | banded form.                               |
!     %--------------------------------------------%
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
      idiag = kl+ku+1
      do 30 j = 1, n
         a(idiag,j) = 4.0E+0
         m(idiag,j) = 4.0E+0
  30  continue
!
!     %-------------------------------------%
!     | First subdiagonal and superdiagonal |
!     %-------------------------------------%
!
      isup = kl+ku
      isub = kl+kl+2
      h = one / real (nx+1)
      rho = 1.0E+2
      do 50 i = 1, nx
        lo = (i-1)*nx
        do 40 j = lo+1, lo+nx-1
           a(isup,j+1) = -one+h*rho/two
           a(isub,j) = -one-h*rho/two
  40    continue
  50  continue
!
      do 60 j = 1, n-1
         m(isup,j+1) = one
         m(isub,j) = one
  60  continue
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
            a(isup,nx+j)  = -one
            a(isub,j) = -one
 70      continue
 80   continue
!
!     %------------------------------------------------%
!     | Call ARPACK banded solver to find eigenvalues  |
!     | and eigenvectors. The real parts of the        |
!     | eigenvalues are returned in the first column   |
!     | of D, the imaginary parts are returned in the  |
!     | second column of D.  Eigenvectors are returned |
!     | in the first NCONV (=IPARAM(5)) columns of V.  |
!     %------------------------------------------------%
!
      rvec = .true.
      call snband(rvec, 'A', select, d, d(1,2), v, ldv, sigmar,
     &     sigmai, workev, n, a, m, lda, rfac, cfac, ku, kl,
     &     which, bmat, nev, tol, resid, ncv, v, ldv, iparam,
     &     workd, workl, lworkl, workc, iwork, info)
!
      if ( info .eq. 0) then
!
!        %-----------------------------------%
!        | Print out convergence information |
!        %-----------------------------------%
!
         nconv = iparam(5)
!
         print *, ' '
         print *, ' _NBDR5 '
         print *, ' ====== '
         print *, ' '
         print *, ' The size of the matrix is ', n
         print *, ' Number of eigenvalue requested is ', nev
         print *, ' The number of Arnoldi vectors generated',
     &            ' (NCV) is ', ncv
         print *, ' The number of converged Ritz values is ',
     &              nconv
         print *, ' What portion of the spectrum ', which
         print *, ' The number of Implicit Arnoldi ',
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
         first = .true.
         do 90 j = 1, nconv
!
            if ( d(j,2) .eq. zero ) then
!
!              %--------------------%
!              | Ritz value is real |
!              %--------------------%
!
               call sgbmv('Notranspose', n, n, kl, ku, one,
     &                    a(kl+1,1), lda, v(1,j), 1, zero,
     &                    ax, 1)
               call saxpy(n, -d(j,1), v(1,j), 1, ax, 1)
               d(j,3) = snrm2(n, ax, 1)
               d(j,3) = d(j,3) / abs(d(j,1))
!
            else if ( first ) then
!
!              %------------------------%
!              | Ritz value is complex  |
!              | Residual of one Ritz   |
!              | value of the conjugate |
!              | pair is computed.      |
!              %------------------------%
!
               call sgbmv('Notranspose', n, n, kl, ku, one,
     &                    a(kl+1,1), lda, v(1,j), 1, zero,
     &                    ax, 1)
               call saxpy(n, -d(j,1), v(1,j), 1, ax, 1)
               call saxpy(n, d(j,2), v(1,j+1), 1, ax, 1)
               d(j,3) = snrm2(n, ax, 1)
               call sgbmv('Notranspose', n, n, kl, ku, one,
     &                    a(kl+1,1), lda, v(1,j+1), 1, zero,
     &                    ax, 1)
               call saxpy(n, -d(j,1), v(1,j+1), 1, ax, 1)
               call saxpy(n, -d(j,2), v(1,j), 1, ax, 1)
               d(j,3) = slapy2( d(j,3), snrm2(n, ax, 1) )
               d(j,3) = d(j,3) / slapy2(d(j,1),d(j,2))
               d(j+1,3) = d(j,3)
               first = .false.
            else
               first = .true.
            end if
!
 90      continue

         call smout(6, nconv, 3, d, maxncv, -6,
     &             'Ritz values (Real,Imag) and relative residuals')
      else
!
!        %-------------------------------------%
!        | Either convergence failed, or there |
!        | is error.  Check the documentation  |
!        | for SNBAND.                         |
!        %-------------------------------------%
!
          print *, ' '
          print *, ' Error with _nband, info= ', info
          print *, ' Check the documentation of _nband '
          print *, ' '
!
      end if
!
 9000 end
