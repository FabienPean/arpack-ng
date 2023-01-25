      program cnbdr3
!
!     ... Construct matrices A and M in LAPACK-style band form.
!         Matrices A and M are derived from the finite
!         element discretization of the 1-dimensional
!         convection-diffusion operator
!                         (d^2u/dx^2) + rho*(du/dx)
!         on the interval [0,1] with zero boundary condition using
!         piecewise linear elements.
!
!     ... Call CNBAND to find eigenvalues LAMBDA such that
!                    A*x = M*x*LAMBDA.
!
!     ... Eigenvalues with largest real parts are sought.
!
!     ... Use mode 2 of CNAUPD.
!
!\BeginLib
!
!\Routines called:
!     cnband  ARPACK banded eigenproblem solver.
!     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     claset  LAPACK routine to initialize a matrix to zero.
!     caxpy   Level 1 BLAS that computes y <- alpha*x+y.
!     scnrm2  Level 1 BLAS that computes the norm of a vector.
!     cgbmv   Level 2 BLAS that computes the band matrix vector product.
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
! FILE: nbdr3.F   SID: 2.4   DATE OF SID: 10/20/00   RELEASE: 2
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
     &                   maxbdw=50, lda = maxbdw, ldv = maxn)
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
     &                 resid(maxn), d(maxncv), ax(maxn), mx(maxn)
      Real
     &                 rwork(maxn), rd(maxncv,3)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character        which*2, bmat
      integer          nev, ncv, ku, kl, info, j,
     &                 n, idiag, isup, isub, maxitr,
     &                 mode, nconv
      logical          rvec
      Real
     &                 tol
      Complex
     &                 rho, h, sigma
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Complex
     &                 one, zero, two
      parameter        (one  = (1.0E+0, 0.0E+0) ,
     &                  zero = (0.0E+0, 0.0E+0) ,
     &                  two  = (2.0E+0, 0.0E+0) )
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
!     | matrix.  A generalized eigenvalue problem is    |
!     | solved (BMAT = 'G').  NEV is the number of      |
!     | eigenvalues to be approximated. The user can    |
!     | modify NX, NEV, NCV and WHICH to solve problems |
!     | of different sizes, and to get different parts  |
!     | the spectrum.  However, The following           |
!     | conditions must be satisfied:                   |
!     |                   N <= MAXN                     |
!     |                 NEV <= MAXNEV                   |
!     |           NEV + 2 <= NCV <= MAXNCV              |
!     %-------------------------------------------------%
!
      n    = 100
      nev  = 4
      ncv  = 10
      if ( n .gt. maxn ) then
         print *, ' ERROR with _NBDR3: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NBDR3: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NBDR3: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat = 'G'
      which = 'LM'
      sigma = zero
!
!     %----------------------------------------------------%
!     | The work array WORKL is used in CNAUPD as          |
!     | workspace.  Its dimension LWORKL has to be set as  |
!     | illustrated below.  The parameter TOL determines   |
!     | the stopping criterion. If TOL<=0, machine machine |
!     | precision is used.  Setting INFO=0 indicates that  |
!     | using a randomly generated vector to start the     |
!     | the ARNOLDI process.                               |
!     %----------------------------------------------------%
!
      lworkl  = 3*ncv**2+5*ncv
      info = 0
      tol  = 0.0
!
!     %---------------------------------------------------%
!     | IPARAm(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 2 of CNAUPD is used     |
!     | (IPARAm(7) = 2). All these options can be changed |
!     | by the user. For details, see the documentation   |
!     | in cnband.                                        |
!     %---------------------------------------------------%
!
      maxitr = 300
      mode   = 2
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
      kl   = 1
      ku   = 1
!
!     %---------------%
!     | Main diagonal |
!     %---------------%
!
      h = one / cmplx(n+1)
!
      idiag = kl+ku+1
      do 30 j = 1, n
         a(idiag,j) = (2.0E+0, 0.0E+0)  / h
         m(idiag,j) = (4.0E+0, 0.0E+0)  * h
  30  continue
!
!     %-------------------------------------%
!     | First subdiagonal and superdiagonal |
!     %-------------------------------------%
!
      isup = kl+ku
      isub = kl+ku+2
      rho = (1.0E+1, 0.0E+0)
      do 40 j = 1, n-1
           a(isup,j+1) = -one/h + rho/two
           a(isub,j) = -one/h - rho/two
           m(isup,j+1) = one*h
           m(isub,j) = one*h
  40  continue
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
         print *, '_NBDR3 '
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
         do 50 j = 1, nconv
!
!           %----------------------------%
!           | Compute the residual norm. |
!           |    ||  A*x - lambda*x ||   |
!           %----------------------------%
!
            call cgbmv('Notranspose', n, n, kl, ku, one,
     &                 a(kl+1,1), lda, v(1,j), 1, zero,
     &                 ax, 1)
            call cgbmv('Notranspose', n, n, kl, ku, one,
     &                 m(kl+1,1), lda, v(1,j), 1, zero,
     &                 mx, 1)
            call caxpy(n, -d(j), mx, 1, ax, 1)
            rd(j,1) = real (d(j))
            rd(j,2) = aimag(d(j))
            rd(j,3) = scnrm2(n, ax, 1)
            rd(j,3) = rd(j,3) / slapy2(rd(j,1), rd(j,2))
 50      continue

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
          print *, ' Error with _band, info= ', info
          print *, ' Check the documentation of _band '
          print *, ' '
!
      end if
!
 9000 end
