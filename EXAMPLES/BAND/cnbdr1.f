      program cnbdr1
!
!     ... Construct the matrix A in LAPACK-style band form.
!         The matrix A is derived from the discretization of
!         the 2-d convection-diffusion operator
!
!              -Laplacian(u) + rho*partial(u)/partial(x).
!
!         on the unit square with zero Dirichlet boundary condition
!         using standard central difference.
!
!     ... Call CNBAND to find eigenvalues LAMBDA such that
!                          A*x = x*LAMBDA.
!
!     ... Use mode 1 of CNAUPD.
!
!\BeginLib
!
!     cnband  ARPACK banded eigenproblem solver.
!     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     claset  LAPACK routine to initialize a matrix to zero.
!     caxpy   Level 1 BLAS that computes y <- alpha*x+y.
!     scnrm2  Level 1 BLAS that computes the norm of a vector.
!     cgbmv   Level 2 BLAS that computes the band matrix vector product
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
! FILE: nbdr1.F   SID: 2.3   DATE OF SID: 08/26/96   RELEASE: 2
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
      Complex
     &                 a(lda,maxn), m(lda,maxn), fac(lda,maxn),
     &                 workl(3*maxncv*maxncv+5*maxncv), workd(3*maxn),
     &                 workev(2*maxncv), v(ldv, maxncv),
     &                 resid(maxn), d(maxncv), ax(maxn)
      Real
     &                 rwork(maxn), rd(maxncv,3)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character        which*2, bmat
      integer          nev, ncv, kl, ku, info, i, j,
     &                 n, nx, lo, isub, isup, idiag, maxitr, mode,
     &                 nconv
      logical          rvec
      Real
     &                 tol
      Complex
     &                 rho, h, h2, sigma
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Complex
     &                 one, zero, two
      parameter        ( one = (1.0E+0, 0.0E+0) ,
     &                   zero = (0.0E+0, 0.0E+0) ,
     &                   two = (2.0E+0, 0.0E+0)  )
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      Real
     &                  scnrm2, slapy2
      external          scnrm2, cgbmv, caxpy, slapy2, claset
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %-------------------------------------------------%
!     | The number NX is the number of interior points  |
!     | in the discretization of the 2-dimensional      |
!     | convection-diffusion operator on the unit       |
!     | square with zero Dirichlet boundary condition.  |
!     | The number N(=NX*NX) is the dimension of the    |
!     | matrix.  A standard eigenvalue problem is       |
!     | solved (BMAT = 'I').  NEV is the number of      |
!     | eigenvalues to be approximated. The user can    |
!     | modify NX, NEV, NCV and WHICH to solve problems |
!     | of different sizes, and to get different parts  |
!     | the spectrum.  However, the following           |
!     | conditions must be satisfied:                   |
!     |                   N <= MAXN                     |
!     |                 NEV <= MAXNEV                   |
!     |           NEV + 2 <= NCV <= MAXNCV              |
!     %-------------------------------------------------%
!
      nx  = 10
      n    = nx*nx
      nev  = 4
      ncv  = 10
      if ( n .gt. maxn ) then
         print *, ' ERROR with _NBDR1: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NBDR1: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NBDR1: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat = 'I'
      which = 'LM'
!
!     %-----------------------------------------------------%
!     | The work array WORKL is used in CNAUPD as           |
!     | workspace.  Its dimension LWORKL is set as          |
!     | illustrated below.  The parameter TOL determines    |
!     | the stopping criterion. If TOL<=0, machine          |
!     | precision is used.  Setting INFO=0 indicates that a |
!     | random vector is generated in CNAUPD to start the   |
!     | Arnoldi iteration.                                  |
!     %-----------------------------------------------------%
!
      lworkl  = 3*ncv**2+5*ncv
      tol  = 0.0
      info = 0
!
!     %---------------------------------------------------%
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of CNAUPD is used     |
!     | (IPARAM(7) = 1). All these options can be changed |
!     | by the user. For details, see the documentation   |
!     | in cnband.                                        |
!     %---------------------------------------------------%
!
      maxitr = 300
      mode   = 1
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
      call claset('A', lda, n, zero, zero, a, lda)
      call claset('A', lda, n, zero, zero, m, lda)
      call claset('A', lda, n, zero, zero, fac, lda)
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
      h = one / cmplx(nx+1)
      h2 = h*h
!
      idiag = kl+ku+1
      do 30 j = 1, n
         a(idiag,j) = (4.0E+0, 0.0E+0)  / h2
  30  continue
!
!     %-------------------------------------%
!     | First subdiagonal and superdiagonal |
!     %-------------------------------------%
!
      rho = (1.0E+2, 0.0E+0)
      isup = kl+ku
      isub = kl+ku+2
      do 50 i = 1, nx
        lo = (i-1)*nx
        do 40 j = lo+1, lo+nx-1
           a(isup,j+1) = -one/h2 + rho/two/h
           a(isub,j) = -one/h2 - rho/two/h
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
!     %-----------------------------------------------%
!     | Call ARPACK banded solver to find eigenvalues |
!     | and eigenvectors. Eigenvalues are returned in |
!     | the one dimensional array D.  Eigenvectors    |
!     | are returned in the first NCONV (=IPARAM(5))  |
!     | columns of V.                                 |
!     %-----------------------------------------------%
!
      rvec = .true.
      call cnband(rvec, 'A', select, d, v, ldv, sigma,
     &           workev, n, a, m, lda, fac, kl, ku, which,
     &           bmat, nev, tol, resid, ncv, v, ldv, iparam,
     &           workd, workl, lworkl, rwork, iwork, info)
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
         print *, '_NBDR1 '
         print *, '====== '
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
         do 90 j = 1, nconv
!
!           %---------------------------%
!           | Compute the residual norm |
!           |   ||  A*x - lambda*x ||   |
!           %---------------------------%
!
            call cgbmv('Notranspose', n, n, kl, ku, one,
     &                 a(kl+1,1), lda, v(1,j), 1, zero,
     &                 ax, 1)
            call caxpy(n, -d(j), v(1,j), 1, ax, 1)
            rd(j,1) = real (d(j))
            rd(j,2) = aimag(d(j))
            rd(j,3) = scnrm2(n, ax, 1)
            rd(j,3) = rd(j,3) / slapy2(rd(j,1),rd(j,2))
 90      continue

         call smout(6, nconv, 3, rd, maxncv, -6,
     &             'Ritz values (Real,Imag) and relative residuals')
      else
!
!        %-------------------------------------%
!        | Either convergence failed, or there |
!        | is error.  Check the documentation  |
!        | for cnband.                         |
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
