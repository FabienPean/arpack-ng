      program pdsdrv1
!
!     Message Passing Layer: BLACS
!
!     Simple program to illustrate the idea of reverse communication
!     in regular mode for a standard symmetric eigenvalue problem.
!
!     We implement example one of ex-sym.doc in SRC directory
!
!\Example-1
!     ... Suppose we want to solve A*x = lambda*x in regular mode,
!         where A is derived from the central difference discretization
!         of the 2-dimensional Laplacian on the unit square with
!         zero Dirichlet boundary condition.
!     ... OP = A  and  B = I.
!     ... Assume "call av (n,x,y)" computes y = A*x.
!     ... Use mode 1 of DSAUPD.
!
!\BeginLib
!
!\Routines called:
!     pdsaupd  Parallel ARPACK reverse communication interface routine.
!     pdseupd  Parallel ARPACK routine that returns Ritz values and (optionally)
!              Ritz vectors.
!     pdnorm2  Parallel version of Level 1 BLAS that computes the norm of a vector.
!     daxpy    Level 1 BLAS that computes y <- alpha*x+y.
!     av       Matrix vector multiplication routine that computes A*x.
!     tv       Matrix vector multiplication routine that computes T*x,
!              where T is a tridiagonal matrix.  It is used in routine
!              av.
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
!\Parallel Modifications
!     Kristi Maschhoff
!
!\Revision history:
!     Starting Point: Serial Code FILE: sdrv1.F   SID: 2.2
!
! FILE: sdrv1.F   SID: 1.4   DATE OF SID: 3/19/97   RELEASE: 1
!
!\Remarks
!     1. None
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      include 'debug.h'
      include 'stat.h'

!     %-----------------%
!     | BLACS INTERFACE |
!     %-----------------%
!
      integer           comm, iam, nprocs, nloc,
     &                  nprow, npcol, myprow, mypcol
!
      external          BLACS_PINFO, BLACS_SETUP, BLACS_GET,
     &                  BLACS_GRIDINIT, BLACS_GRIDINFO
!
!     %-----------------------------%
!     | Define leading dimensions   |
!     | for all arrays.             |
!     | MAXN:   Maximum dimension   |
!     |         of the A allowed.   |
!     | MAXNEV: Maximum NEV allowed |
!     | MAXNCV: Maximum NCV allowed |
!     %-----------------------------%
!
      integer          maxnloc, maxnev, maxncv, ldv
      parameter       (maxnloc=256, maxnev=10, maxncv=25,
     &                 ldv=maxnloc )
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
      Double precision
     &                 v(ldv,maxncv), workl(maxncv*(maxncv+8)),
     &                 workd(3*maxnloc), d(maxncv,2), resid(maxnloc),
     &                 ax(maxnloc)
      logical          select(maxncv)
      integer          iparam(11), ipntr(11)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character        bmat*1, which*2
      integer          ido, n, nev, ncv, lworkl, info, ierr, j,
     &                 nx, nconv, maxitr, mode, ishfts
      logical          rvec
      Double precision
     &                 tol, sigma
!
!     %----------------------------------------------%
!     | Local Buffers needed for BLACS communication |
!     %----------------------------------------------%
!
      Double precision
     &                  mv_buf(maxnloc)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Double precision
     &                 zero
      parameter        (zero = 0.0)
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      Double precision
     &                 pdnorm2
      external         pdnorm2, daxpy
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic         abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      call BLACS_PINFO( iam, nprocs )
!
!     If in PVM, create virtual machine if it doesn't exist
!
      if (nprocs .lt. 1) then
         if (iam .eq. 0) then
              write(*,1000)
              read(*, 2000) nprocs
         endif
         call BLACS_SETUP( iam, nprocs )
      endif
!
1000  format('How many processes in machine?')
2000  format(I3)
!
!     Set up processors in 1D Grid
!
      nprow = nprocs
      npcol = 1
!
!     Get default system context, and define grid
!
      call BLACS_GET( 0, 0, comm )
      call BLACS_GRIDINIT( comm, 'Row', nprow, npcol )
      call BLACS_GRIDINFO( comm, nprow, npcol, myprow, mypcol )
!
!     If I'm not in grid, go to end of program
!
      if ( (myprow .ge. nprow) .or. (mypcol .ge. npcol) ) goto 9000
!
      ndigit = -3
      logfil = 6
      msaupd = 1
!
!     %----------------------------------------------------%
!     | The number NX is the number of interior points     |
!     | in the discretization of the 2-dimensional         |
!     | Laplacian on the unit square with zero Dirichlet   |
!     | boundary condition.  The number N(=NX*NX) is the   |
!     | dimension of the matrix.  A standard eigenvalue    |
!     | problem is solved (BMAT = 'I'). NEV is the number  |
!     | of eigenvalues to be approximated.  The user can   |
!     | modify NEV, NCV, WHICH to solve problems of        |
!     | different sizes, and to get different parts of the |
!     | spectrum.  However, The following conditions must  |
!     | be satisfied:                                      |
!     |                   N <= MAXN,                       |
!     |                 NEV <= MAXNEV,                     |
!     |             NEV + 2 <= NCV <= MAXNCV               |
!     %----------------------------------------------------%
!
      nx = 10
      n = nx*nx
      nev =  4
      ncv =  20
!
!     %--------------------------------------%
!     | Set up distribution of data to nodes |
!     %--------------------------------------%
!
      nloc = (nx / nprocs)*nx
      if ( mod(nx, nprocs) .gt. myprow ) nloc = nloc + nx
!
      if ( nloc .gt. maxnloc ) then
         print *, ' ERROR with _SDRV1: NLOC is greater than MAXNLOC '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _SDRV1: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _SDRV1: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat = 'I'
      which = 'SM'
!
!     %--------------------------------------------------%
!     | The work array WORKL is used in PSSAUPD as       |
!     | workspace.  Its dimension LWORKL is set as       |
!     | illustrated below.  The parameter TOL determines |
!     | the stopping criterion.  If TOL<=0, machine      |
!     | precision is used.  The variable IDO is used for |
!     | reverse communication and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is |
!     | generated in PSSAUPD to start the Arnoldi        |
!     | iteration.                                       |
!     %--------------------------------------------------%
!
      lworkl = ncv*(ncv+8)
      tol = zero
      info = 0
      ido = 0
!
!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of PSSAUPD is used    |
!     | (IPARAM(7) = 1).  All these options may be        |
!     | changed by the user. For details, see the         |
!     | documentation in PSSAUPD.                         |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 300
      mode   = 1
!
      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%
!
 10   continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine PSSAUPD and take|
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         call pdsaupd ( comm, ido, bmat, nloc, which, nev, tol, resid,
     &                 ncv, v, ldv, iparam, ipntr, workd, workl,
     &                 lworkl, info )
!
         if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %--------------------------------------%
!           | Perform matrix vector multiplication |
!           |              y <--- OP*x             |
!           | The user should supply his/her own   |
!           | matrix vector multiplication routine |
!           | here that takes workd(ipntr(1)) as   |
!           | the input, and return the result to  |
!           | workd(ipntr(2)).                     |
!           %--------------------------------------%
!
            call av ( comm, nloc, nx, mv_buf,
     &               workd(ipntr(1)), workd(ipntr(2)))
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call PSSAUPD again.|
!           %-----------------------------------------%
!
            go to 10
!
         end if
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
      if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message. Check the |
!        | documentation in PSSAUPD.|
!        %--------------------------%
!
         if ( myprow .eq. 0 ) then
            print *, ' '
            print *, ' Error with _saupd, info = ', info
            print *, ' Check documentation in _saupd '
            print *, ' '
         endif
!
      else
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using PSSEUPD.               |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        %-------------------------------------------%
!
         rvec = .true.
!
         call pdseupd ( comm, rvec, 'All', select,
     &        d, v, ldv, sigma,
     &        bmat, nloc, which, nev, tol, resid, ncv, v, ldv,
     &        iparam, ipntr, workd, workl, lworkl, ierr )
!        %----------------------------------------------%
!        | Eigenvalues are returned in the first column |
!        | of the two dimensional array D and the       |
!        | corresponding eigenvectors are returned in   |
!        | the first NEV columns of the two dimensional |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
!
         if ( ierr .ne. 0) then
!
!            %------------------------------------%
!            | Error condition:                   |
!            | Check the documentation of PSSEUPD.|
!            %------------------------------------%
!
!
            if ( myprow .eq. 0 ) then
             	print *, ' '
             	print *, ' Error with _seupd, info = ', ierr
             	print *, ' Check the documentation of _seupd. '
             	print *, ' '
            endif
!
         else
!
             nconv =  iparam(5)
             do 20 j=1, nconv
!
!               %---------------------------%
!               | Compute the residual norm |
!               |                           |
!               |   ||  A*x - lambda*x ||   |
!               |                           |
!               | for the NCONV accurately  |
!               | computed eigenvalues and  |
!               | eigenvectors.  (iparam(5) |
!               | indicates how many are    |
!               | accurate to the requested |
!               | tolerance)                |
!               %---------------------------%
!
                call av(comm, nloc, nx, mv_buf, v(1,j), ax)
                call daxpy(nloc, -d(j,1), v(1,j), 1, ax, 1)
                d(j,2) = pdnorm2( comm, nloc, ax, 1 )
!
 20          continue
!
!            %-------------------------------%
!            | Display computed residuals    |
!            %-------------------------------%
!
             call pdmout(comm, 6, nconv, 2, d, maxncv, -6,
     &            'Ritz values and direct residuals')
         end if
!
!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%
!
         if (myprow .eq. 0)then
         if ( info .eq. 1) then
            print *, ' '
            print *, ' Maximum number of iterations reached.'
            print *, ' '
         else if ( info .eq. 3) then
            print *, ' '
            print *, ' No shifts could be applied during implicit
     &                 Arnoldi update, try increasing NCV.'
            print *, ' '
         end if
!
         print *, ' '
         print *, '_SDRV1 '
         print *, '====== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of processors is ', nprocs
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated',
     &            ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ',
     &              nconv
         print *, ' The number of Implicit Arnoldi update',
     &            ' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
         endif
!
      end if
!
!     %---------------------------%
!     | Done with program pdsdrv1.|
!     %---------------------------%
!
 9000 continue
!
      call BLACS_GRIDEXIT ( comm )
      call BLACS_EXIT(0)
!
      end
!
! ------------------------------------------------------------------
!     parallel matrix vector subroutine
!
!     The matrix used is the 2 dimensional discrete Laplacian on unit
!     square with zero Dirichlet boundary condition.
!
!     Computes w <--- OP*v, where OP is the nx*nx by nx*nx block
!     tridiagonal matrix
!
!                  | T -I          |
!                  |-I  T -I       |
!             OP = |   -I  T       |
!                  |        ...  -I|
!                  |           -I T|
!
!     The subroutine TV is called to computed y<---T*x.
!-------------------------------------------------------------------
!
      subroutine av (comm, nloc, nx, mv_buf, v, w)
!
!     .. BLACS Declarations ...
      integer           comm, nprow, npcol, myprow, mypcol
      external          BLACS_GRIDINFO, DGESD2D, DGERV2D
      integer           nloc, nx, np, j, lo, next, prev
      Double precision
     &                  v(nloc), w(nloc), mv_buf(nx), one
      parameter         (one = 1.0 )
      external          daxpy

      call BLACS_GRIDINFO( comm, nprow, npcol, myprow, mypcol )
!
      np = nloc/nx
      call tv(nx,v(1),w(1))
      call daxpy(nx, -one, v(nx+1), 1, w(1), 1)
!
      if ( np .gt. 2) then
         do 10 j = 2, np-1
            lo = (j-1)*nx
            call tv(nx, v(lo+1), w(lo+1))
            call daxpy(nx, -one, v(lo-nx+1), 1, w(lo+1), 1)
            call daxpy(nx, -one, v(lo+nx+1), 1, w(lo+1), 1)
  10     continue
      end if
!
      if ( np .gt. 1) then
         lo = (np-1)*nx
         call tv(nx, v(lo+1), w(lo+1))
         call daxpy(nx, -one, v(lo-nx+1), 1, w(lo+1), 1)
      end if
!
      next = myprow + 1
      prev = myprow - 1
      if ( myprow .lt. nprow-1 ) then
         call dgesd2d( comm, nx, 1, v((np-1)*nx+1), nx, next, mypcol)
      endif
      if ( myprow .gt. 0 ) then
         call dgerv2d( comm, nx, 1, mv_buf, nx, prev, mypcol )
         call daxpy( nx, -one, mv_buf, 1, w(1), 1 )
      endif
!
      if ( myprow .gt. 0 ) then
         call dgesd2d( comm, nx, 1, v(1), nx, prev, mypcol)
      endif
      if ( myprow .lt. nprow-1 ) then
         call dgerv2d( comm, nx, 1, mv_buf, nx, next, mypcol )
         call daxpy( nx, -one, mv_buf, 1, w(lo+1), 1 )
      endif
!
      return
      end
!=========================================================================
      subroutine tv (nx, x, y)
!
      integer           nx, j
      Double precision
     &                  x(nx), y(nx), dd, dl, du
!
      Double precision
     &                 one
      parameter        (one = 1.0 )
!
!     Compute the matrix vector multiplication y<---T*x
!     where T is a nx by nx tridiagonal matrix with DD on the
!     diagonal, DL on the subdiagonal, and DU on the superdiagonal.
!
!
      dd  = 4.0
      dl  = -one
      du  = -one
!
      y(1) =  dd*x(1) + du*x(2)
      do 10 j = 2,nx-1
         y(j) = dl*x(j-1) + dd*x(j) + du*x(j+1)
 10   continue
      y(nx) =  dl*x(nx-1) + dd*x(nx)
      return
      end
