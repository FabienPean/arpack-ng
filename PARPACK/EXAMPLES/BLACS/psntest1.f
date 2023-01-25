      program psntest1
!
!     Message Passing Layer: BLACS
!
!     Example program to illustrate the idea of reverse communication
!     for a standard nonsymmetric eigenvalue problem.
!
!     We implement example one of ex-nonsym.doc in DOCUMENTS directory
!
!\Test-1
!     ... Suppose we want to solve A*x = lambda*x in regular mode,
!         where A is random diagonal matrix with 4 separated eigenvalues.
!     ... OP = A  and  B = I.
!     ... Assume "call av ( nloc, diag, x, y)" computes y = A*x.
!     ... Use mode 1 of PSNAUPD.
!
!\BeginLib
!
!\Routines called:
!     psnaupd  Parallel ARPACK reverse communication interface routine.
!     psneupd  Parallel ARPACK routine that returns Ritz values and (optionally)
!              Ritz vectors.
!     slapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     saxpy    Level 1 BLAS that computes y <- alpha*x+y.
!     psnorm2  Parallel version of Level 1 BLAS that computes the norm of a vector.
!     av       Distributed matrix vector multiplication routine that computes A*x.
!
!\Author
!     Kristi Maschhoff
!     Dept. of Computational &
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\SCCS Information:
! FILE: %M%   SID: %I%   DATE OF SID: %G%   RELEASE: %R%
!
!\Remarks
!     1. None
!
!\EndLib
!---------------------------------------------------------------------------
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
!     | Define maximum dimensions   |
!     | for all arrays.             |
!     | MAXN:   Maximum dimension   |
!     |         of the distributed  |
!     |         block of A allowed. |
!     | MAXNEV: Maximum NEV allowed |
!     | MAXNCV: Maximum NCV allowed |
!     %-----------------------------%
!
      integer           maxn, maxnev, maxncv, ldv
      parameter         (maxn=256, maxnev=12, maxncv=30, ldv=maxn)
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
      integer           iparam(11), ipntr(14), iseed(4)
      logical           select(maxncv)
      Real
     &                  ax(maxn), d(maxncv,3), resid(maxn), diag(maxn),
     &                  v(ldv,maxncv), workd(3*maxn),
     &                  workev(3*maxncv),
     &                  workl(3*maxncv*maxncv+6*maxncv)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character         bmat*1, which*2
      integer           ido, n, nx, nev, ncv, lworkl, info, j,
     &                  ierr, nconv, maxitr, ishfts, mode, idist
      Real
     &                  tol, sigmar, sigmai
      logical           first, rvec
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real
     &                  zero, one
      parameter         (zero = 0.0, one = 1.0)
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      Real
     &                  slapy2, psnorm2
      external          slapy2, saxpy, psnorm2, slarnv
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic         abs, sqrt
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
      mnaupd = 1
!
      n     = 10000*nprocs
      nev   = 4
      ncv   = 20
!
!     %--------------------------------------%
!     | Set up distribution of data to nodes |
!     %--------------------------------------%
!
      nloc = 10000
!
      if ( nloc .gt. maxn ) then
         print *, ' ERROR with _NDRV1: NLOC is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NDRV1: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NDRV1: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat  = 'I'
      which = 'SM'
!
!     %-----------------------------------%
!     | Generate random diagonal matrix   |
!     | Isolate 4 extreamal eigenvalues   |
!     %-----------------------------------%
!
      idist = 1
      iseed(1) = 15
      iseed(2) = 35
      iseed(3) = 52
      iseed(4) = 7
      call slarnv ( idist, iseed, nloc, diag )
      diag(1) = diag(1) + 1.01
      diag(2) = diag(2) + 1.01
      diag(3) = diag(3) + 1.01
      diag(4) = diag(4) + 1.01
!
!     %-----------------------------------------------------%
!     | The work array WORKL is used in PSNAUPD as          |
!     | workspace.  Its dimension LWORKL is set as          |
!     | illustrated below.  The parameter TOL determines    |
!     | the stopping criterion. If TOL<=0, machine          |
!     | precision is used.  The variable IDO is used for    |
!     | reverse communication, and is initially set to 0.   |
!     | Setting INFO=0 indicates that a random vector is    |
!     | generated in PSNAUPD to start the Arnoldi iteration.|
!     %-----------------------------------------------------%
!
      lworkl  = 3*ncv**2+6*ncv
      tol    = zero
      ido    = 0
      info   = 1
      do 50 j=1,nloc
         resid(j) = one
 50   continue
!
!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of PSNAUPD is used    |
!     | (IPARAM(7) = 1). All these options can be changed |
!     | by the user. For details see the documentation in |
!     | PSNAUPD.                                          |
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
!        | Repeatedly call the routine PSNAUPD and take|
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         call psnaupd(comm, ido, bmat, nloc, which, nev, tol, resid,
     &        ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info )
!
         if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %-------------------------------------------%
!           | Perform matrix vector multiplication      |
!           |                y <--- OP*x                |
!           | The user should supply his/her own        |
!           | matrix vector multiplication routine here |
!           | that takes workd(ipntr(1)) as the input   |
!           | vector, and return the matrix vector      |
!           | product to workd(ipntr(2)).               |
!           %-------------------------------------------%
!
            call av ( nloc, diag, workd(ipntr(1)), workd(ipntr(2)))
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call PSNAUPD again.|
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
!        | Error message, check the |
!        | documentation in PSNAUPD.|
!        %--------------------------%
!
         if ( myprow .eq. 0 ) then
            print *, ' '
            print *, ' Error with _naupd, info = ', info
            print *, ' Check the documentation of _naupd'
            print *, ' '
         endif
!
      else
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using PSNEUPD.               |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        %-------------------------------------------%
!
         rvec = .true.
!
         call psneupd ( comm, rvec, 'A', select, d, d(1,2), v, ldv,
     &        sigmar, sigmai, workev, bmat, nloc, which, nev, tol,
     &        resid, ncv, v, ldv, iparam, ipntr, workd, workl,
     &        lworkl, ierr )
!
!        %-----------------------------------------------%
!        | The real part of the eigenvalue is returned   |
!        | in the first column of the two dimensional    |
!        | array D, and the imaginary part is returned   |
!        | in the second column of D.  The corresponding |
!        | eigenvectors are returned in the first NEV    |
!        | columns of the two dimensional array V if     |
!        | requested.  Otherwise, an orthogonal basis    |
!        | for the invariant subspace corresponding to   |
!        | the eigenvalues in D is returned in V.        |
!        %-----------------------------------------------%
!
         if ( ierr .ne. 0) then
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of PSNEUPD.|
!           %------------------------------------%
!
         	if ( myprow .eq. 0 ) then
             	print *, ' '
             	print *, ' Error with _neupd, info = ', ierr
             	print *, ' Check the documentation of _neupd. '
             	print *, ' '
            endif
!
         else
!
             first  = .true.
             nconv  = iparam(5)
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
                if (d(j,2) .eq. zero)  then
!
!                  %--------------------%
!                  | Ritz value is real |
!                  %--------------------%
!
                   call av( nloc, diag, v(1,j), ax)
                   call saxpy(nloc, -d(j,1), v(1,j), 1, ax, 1)
                   d(j,3) = psnorm2( comm, nloc, ax, 1)
!
                else if (first) then
!
!                  %------------------------%
!                  | Ritz value is complex. |
!                  | Residual of one Ritz   |
!                  | value of the conjugate |
!                  | pair is computed.      |
!                  %------------------------%
!
                   call av( nloc, diag, v(1,j), ax)
                   call saxpy(nloc, -d(j,1), v(1,j), 1, ax, 1)
                   call saxpy(nloc, d(j,2), v(1,j+1), 1, ax, 1)
                   d(j,3) = psnorm2( comm, nloc, ax, 1)
                   call av( nloc, diag, v(1,j+1), ax)
                   call saxpy(nloc, -d(j,2), v(1,j), 1, ax, 1)
                   call saxpy(nloc, -d(j,1), v(1,j+1), 1, ax, 1)
                   d(j,3) = slapy2(d(j,3), psnorm2(comm,nloc,ax,1) )
                   d(j+1,3) = d(j,3)
                   first = .false.
                else
                   first = .true.
                end if
!
 20          continue
!
!            %-----------------------------%
!            | Display computed residuals. |
!            %-----------------------------%
!
             call psmout(comm, 6, nconv, 3, d, maxncv, -6,
     &            'Ritz values (Real,Imag) and direct residuals')
          end if
!
!        %-------------------------------------------%
!        | Print additional convergence information. |
!        %-------------------------------------------%
!
         if (myprow .eq. 0)then
         if ( info .eq. 1) then
             print *, ' '
             print *, ' Maximum number of iterations reached.'
             print *, ' '
         else if ( info .eq. 3) then
             print *, ' '
             print *, ' No shifts could be applied during implicit
     &                  Arnoldi update, try increasing NCV.'
             print *, ' '
         end if
!
         print *, ' '
         print *, '_NDRV1 '
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
!
         endif
      end if
!
!     %---------------------------%
!     | Done with program pdndrv1.|
!     %---------------------------%
!
 9000 continue
!
!     %-------------------------%
!     | Release resources BLACS |
!     %-------------------------%
!
      call BLACS_GRIDEXIT ( comm )
      call BLACS_EXIT(0)
!
      end
!
!==========================================================================
!
!     parallel matrix vector subroutine
!
      subroutine av (n, diag, v, w)
      integer           n, j
      Real
     &                  v(n), w(n), diag(n)
!
      do 10 j = 1, n
         w(j) = diag(j)*v(j)
  10  continue
!
      return
      end
