      subroutine mmread(iunit,rep,field,symm,rows,cols,nnz,nnzmax,
     *                 indx,jndx,ival,rval,cval)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! This routine will read data from a matrix market formatted file.
! The data may be either sparse coordinate format, or dense array format.
!
! The unit iunit must be open, and the file will be rewound on return.
!
! 20-Sept-96  Karin A. Remington, NIST ACMD (karin@cam.nist.gov)
! 18-Oct-96   Change in routine name to match C and Matlab routines.
! 30-Oct-96   Bug fixes in mmio.f:
!                  -looping for comment lines
!                  -fixed non-ansi zero stringlength
!                  -incorrect size calculation for skew-symmetric arrays
! 	      Other changes in mmio.f:
!                  -added integer value parameter to calling sequences
!                  -enforced proper count in size info line
!                  -added routine to count words in string (countwd)
!            (Thanks to G.P.Leendetse and H.Oudshoom for their review
!             of the initial version and suggested fixes.)
! 15-Oct-08  fixed illegal attempt of mimicking "do while" construct
!            by redifing limits inside loop. (lines 443-450)
!            (Thanks to Geraldo Veiga for his comments.)
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!   Arguments:
!
!   name     type      in/out description
!   ---------------------------------------------------------------
!
!   iunit    integer     in   Unit identifier for the file
!                             containing the data to be read.
!                             Must be open prior to call.
!                             Will be rewound on return.
!
!   rep     character*10 out  Matrix Market 'representation'
!                             indicator. On return:
!
!                                coordinate   (for sparse data)
!                                array        (for dense data)
!                                elemental    (to be added)
!
!   field   character*7  out  Matrix Market 'field'. On return:
!
!                                real
!                                complex
!                                integer
!                                pattern
!
!   symm    character*19 out  Matrix Market 'field'. On return:
!
!                                symmetric
!                                hermitian
!                                skew-symmetric
!                                general
!
!   rows     integer     out  Number of rows in matrix.
!
!   cols     integer     out  Number of columns in matrix.
!
!   nnz      integer     out  Number of nonzero entries required to
!                             store matrix.
!
!   nnzmax   integer     in   Maximum dimension of data arrays.
!
!   indx     integer(nnz)out  Row indices for coordinate format.
!                             Undefined for array format.
!
!   jndx     integer(nnz)out  Column indices for coordinate format.
!                             Undefined for array format.
!
!   ival     integer(nnz) out Integer data (if applicable, see 'field')
!
!   rval     double(nnz) out  Real data (if applicable, see 'field')
!
!   cval     complex(nnz)out  Complex data (if applicable, see 'field')
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! Declarations:
!
      integer ival(*)
      double precision rval(*)
      complex cval(*)
      double precision rpart,ipart
      integer indx(*)
      integer jndx(*)
      integer i, rows, cols, nnz, nnzreq, nnzmax, iunit
      integer count
      character mmhead*15
      character mmtype*6
      character rep*10
      character field*7
      character symm*19
      character tmp1*1024
      character tmp2*2
!
! Read header line and check validity:
!
      read (iunit,end=1000,fmt=5) tmp1
 5    format(1024A)
      call getwd(mmhead,tmp1,1024,1,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(mmtype,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(rep,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(field,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(symm,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      if ( mmhead .ne. '%%MatrixMarket' ) go to 5000
!
! Convert type code to lower case for easier comparisons:
!
      call lowerc(mmtype,1,6)
      if ( mmtype .ne. 'matrix' ) then
         print *,'Invalid matrix type: ',mmtype
         print *,'This reader only understands type ''matrix''.'
         stop
      else
         call lowerc(rep,1,10)
         call lowerc(field,1,7)
         call lowerc(symm,1,19)
      endif
!
! Test input qualifiers:
!
      if (rep .ne. 'coordinate' .and. rep .ne. 'array' )
     *   go to 6000
      if (rep .eq. 'coordinate' .and. field .ne. 'integer' .and.
     *    field .ne. 'real' .and. field .ne. 'complex' .and.
     *    field .ne. 'pattern') go to 7000
      if (rep .eq. 'array' .and. field .ne. 'integer' .and.
     *    field .ne. 'real' .and. field .ne. 'complex' ) go to 8000
      if (symm .ne. 'general' .and. symm .ne. 'symmetric' .and.
     *    symm .ne. 'hermitian' .and. symm .ne. 'skew-symmetric')
     *   go to 9000
!
! Read through comment lines, ignoring content:
!
      read (iunit,end=2000,fmt=200) tmp2
 200  format(1a)
 10   continue
        if ( tmp2(1:1) .ne. '%' ) then
           go to 20
        endif
        read (iunit,end=2000,fmt=200) tmp2
        go to 10
 20   continue
!
! Just read a non-comment.
!   Now, back up a line, and read for first int, and back up
!   again. This will set pointer to just before apparent size
!   info line.
!   Before continuing with free form input, count the number of
!   words on the size info line to ensure there is the right amount
!   of info (2 words for array matrices, 3 for coordinate matrices).
!
      backspace (iunit)
      read (iunit,end=1000,fmt=5) tmp1
      call countwd(tmp1,1024,1,count)
      if ( rep .eq. 'array' .and. count .ne. 2 ) go to 3000
      if ( rep .eq. 'coordinate' .and. count .ne. 3 ) go to 3500
!
!   Correct number of words are present, now back up and read them.
!
      backspace (iunit)
!
      if ( rep .eq. 'coordinate' ) then
!
! Read matrix in sparse coordinate format
!
        read (iunit,fmt=*) rows,cols,nnz
!
! Check to ensure adequate storage is available
!
        if ( nnz .gt. nnzmax ) then
          print *,'insufficent array lengths for matrix of ',nnz,
     *            ' nonzeros.'
          print *,'resize nnzmax to at least ',nnz,'. (currently ',
     *            nnzmax,')'
          stop
        endif
!
! Read data according to data type (real,integer,complex, or pattern)
!
        if ( field .eq. 'integer' ) then
          do 30 i=1,nnz
            read (iunit,fmt=*,end=4000) indx(i),jndx(i),ival(i)
 30       continue
        elseif ( field .eq. 'real' ) then
          do 35 i=1,nnz
            read (iunit,fmt=*,end=4000) indx(i),jndx(i),rval(i)
 35       continue
        elseif ( field .eq. 'complex' ) then
          do 40 i=1,nnz
            read (iunit,fmt=*,end=4000) indx(i),jndx(i),rpart,ipart
            cval(i) = cmplx(rpart,ipart)
 40       continue
        elseif ( field .eq. 'pattern' ) then
          do 50 i=1,nnz
            read (iunit,fmt=*,end=4000) indx(i),jndx(i)
 50       continue
        else
           print *,'''',field,''' data type not recognized.'
           stop
        endif
        rewind(iunit)
        return
!
      elseif ( rep .eq. 'array' ) then
!
! Read matrix in dense column-oriented array format
!
        read (iunit,fmt=*) rows,cols
!
! Check to ensure adequate storage is available
!
        if ( symm .eq. 'symmetric' .or. symm .eq. 'hermitian' ) then
          nnzreq = (rows*cols - rows)/2 + rows
          nnz = nnzreq
        elseif ( symm .eq. 'skew-symmetric' ) then
          nnzreq = (rows*cols - rows)/2
          nnz = nnzreq
        else
          nnzreq = rows*cols
          nnz = nnzreq
        endif
        if ( nnzreq .gt. nnzmax ) then
          print *,'insufficent array length for ',rows, ' by ',
     *             cols,' dense ',symm,' matrix.'
          print *,'resize nnzmax to at least ',nnzreq,'. (currently ',
     *             nnzmax,')'
          stop
        endif
!
! Read data according to data type (real,integer,complex, or pattern)
!
        if ( field .eq. 'integer' ) then
          do 60 i=1,nnzreq
            read (iunit,fmt=*,end=4000) ival(i)
 60      continue
        elseif ( field .eq. 'real' ) then
          do 65 i=1,nnzreq
            read (iunit,fmt=*,end=4000) rval(i)
 65      continue
        elseif ( field .eq. 'complex' ) then
          do 70 i=1,nnzreq
            read (iunit,fmt=*,end=4000) rpart,ipart
            cval(i) = cmplx(rpart,ipart)
 70      continue
        else
           print *,'''pattern'' data not consistent with type ''array'''
           stop
        endif
        rewind(iunit)
        return
      else
        print *,'''',rep,''' representation not recognized.'
        print *, 'Recognized representations:'
        print *, '   array'
        print *, '   coordinate'
        stop
      endif
!
! Various error conditions:
!
 1000 print *,'Premature end-of-file.'
      print *,'No lines found.'
      stop
 2000 print *,'Premature end-of-file.'
      print *,'No data lines found.'
      stop
 3000 print *,'Size info inconsistent with representation.'
      print *,'Array matrices need exactly 2 size descriptors.'
      print *, count,' were found.'
      stop
 3500 print *,'Size info inconsistent with representation.'
      print *,'Coordinate matrices need exactly 3 size descriptors.'
      print *, count,' were found.'
      stop
 4000 print *,'Premature end-of-file.'
      print *,'Check that the data file contains ',nnz,
     *        ' lines of  i,j,[val] data.'
      print *,'(it appears there are only ',i,' such lines.)'
      stop
 5000 print *,'Invalid matrix header: ',tmp1
      print *,'Correct header format:'
      print *,'%%MatrixMarket type representation field symmetry'
      print *
      print *,'Check specification and try again.'
 6000 print *,'''',rep,''' representation not recognized.'
      print *, 'Recognized representations:'
      print *, '   array'
      print *, '   coordinate'
      stop
 7000 print *,'''',field,''' field is not recognized.'
      print *, 'Recognized fields:'
      print *, '   real'
      print *, '   complex'
      print *, '   integer'
      print *, '   pattern'
      stop
 8000 print *,'''',field,''' arrays are not recognized.'
      print *, 'Recognized fields:'
      print *, '   real'
      print *, '   complex'
      print *, '   integer'
      stop
 9000 print *,'''',symm,''' symmetry is not recognized.'
      print *, 'Recognized symmetries:'
      print *, '   general'
      print *, '   symmetric'
      print *, '   hermitian'
      print *, '   skew-symmetric'
      stop
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! End of subroutine mmread
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine mminfo(iunit,rep,field,symm,rows,cols,nnz)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! This routine will read header information from a Matrix Market
! formatted file.
!
! The unit iunit must be open, and the file will be rewound on return.
!
! 20-Sept-96  Karin A. Remington, NIST ACMD (karin@cam.nist.gov)
! 18-Oct-96   Change in routine name to match C and Matlab routines.
! 30-Oct-96   Bug fixes in mmio.f:
!                  -looping for comment lines
!                  -fixed non-ansi zero stringlength
!                  -incorrect size calculation for skew-symmetric arrays
! 	      Other changes in mmio.f:
!                  -added integer value parameter to calling sequences
!                  -enforced proper count in size info line
!                  -added routine to count words in string (countwd)
!            (Thanks to G.P.Leendetse and H.Oudshoom for their review
!             of the initial version and suggested fixes.)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!   Arguments:
!
!   name     type      in/out description
!   ---------------------------------------------------------------
!
!   iunit  integer     in   Unit identifier for the open file
!                             containing the data to be read.
!
!   rep     character*10 out  Matrix Market 'representation'
!                             indicator. On return:
!
!                                coordinate   (for sparse data)
!                                array        (for dense data)
!                                elemental    (to be added)
!
!   field   character*7  out  Matrix Market 'field'. On return:
!
!                                real
!                                complex
!                                integer
!                                pattern
!
!   symm    character*19 out  Matrix Market 'field'. On return:
!
!                                symmetric
!                                hermitian
!                                skew-symmetric
!                                general
!
!   rows     integer     out  Number of rows in matrix.
!
!   cols     integer     out  Number of columns in matrix.
!
!   nnz      integer     out  Number of nonzero entries required to store
!                             the matrix.
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! Declarations:
!
      integer i, rows, cols, nnz, iunit
      integer count
      character mmhead*14
      character mmtype*6
      character rep*10
      character field*7
      character symm*19
      character tmp1*1024
      character tmp2*2
!
! Read header line and check validity:
!
      read (iunit,end=1000,fmt=5) tmp1
 5    format(1024A)
!
! Parse words from header line:
!
      call getwd(mmhead,tmp1,1024,1,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(mmtype,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(rep,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(field,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(symm,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      if ( mmhead .ne. '%%MatrixMarket' ) go to 5000
!
! Convert type code to upper case for easier comparisons:
!
      call lowerc(mmtype,1,6)
      if ( mmtype .ne. 'matrix' ) then
         print *,'Invalid matrix type: ',mmtype
         print *,'This reader only understands type ''matrix''.'
        stop
      else
         call lowerc(rep,1,10)
         call lowerc(field,1,7)
         call lowerc(symm,1,19)
      endif
!
! Test input qualifiers:
!
      if (rep .ne. 'coordinate' .and. rep .ne. 'array' )
     *   go to 6000
      if (rep .eq. 'coordinate' .and. field .ne. 'integer' .and.
     *    field .ne. 'real' .and. field .ne. 'complex' .and.
     *    field .ne. 'pattern') go to 7000
      if (rep .eq. 'array' .and. field .ne. 'integer' .and.
     *    field .ne. 'real' .and. field .ne. 'complex' ) go to 8000
      if (symm .ne. 'general' .and. symm .ne. 'symmetric' .and.
     *    symm .ne. 'hermitian' .and. symm .ne. 'skew-symmetric')
     *   go to 9000
!
! Read through comment lines, ignoring content:
!
      read (iunit,end=2000,fmt=200) tmp2
 200  format(1a)
 10   continue
        if ( tmp2(1:1) .ne. '%' ) then
           go to 20
        endif
        read (iunit,end=2000,fmt=200) tmp2
        go to 10
 20   continue
!
! Just read a non-comment.
!   Now, back up a line, and read for first int, and back up
!   again. This will set pointer to just before apparent size
!   info line.
!   Before continuing with free form input, count the number of
!   words on the size info line to ensure there is the right amount
!   of info (2 words for array matrices, 3 for coordinate matrices).
!
      backspace (iunit)
      read (iunit,end=1000,fmt=5) tmp1
      call countwd(tmp1,1024,1,count)
      if ( rep .eq. 'array' .and. count .ne. 2 ) go to 3000
      if ( rep .eq. 'coordinate' .and. count .ne. 3 ) go to 3500
!
!   Correct number of words are present, now back up and read them.
!
      backspace (iunit)
!
      if ( rep .eq. 'coordinate' ) then
!
! Read matrix in sparse coordinate format
!
        read (iunit,fmt=*) rows,cols,nnz
!
! Rewind before returning
!
        rewind(iunit)
        return
!
      elseif ( rep .eq. 'array' ) then
!
! Read matrix in dense column-oriented array format
!
        read (iunit,fmt=*) rows,cols
        if ( symm .eq. 'symmetric' .or. symm .eq. 'hermitian' ) then
          nnz = (rows*cols - rows)/2 + rows
        elseif ( symm .eq. 'skew-symmetric' ) then
          nnz = (rows*cols - rows)/2
        else
          nnz = rows*cols
        endif
!
! Rewind before returning
!
        rewind(iunit)
        return
      else
        print *,'''',rep,''' representation not recognized.'
        print *, 'Recognized representations:'
        print *, '   array'
        print *, '   coordinate'
        stop
      endif
!
! Various error conditions:
!
 1000 print *,'Premature end-of-file.'
      print *,'No lines found.'
      stop
 2000 print *,'Premature end-of-file.'
      print *,'No data found.'
      stop
 3000 print *,'Size info inconsistent with representation.'
      print *,'Array matrices need exactly 2 size descriptors.'
      print *, count,' were found.'
      stop
 3500 print *,'Size info inconsistent with representation.'
      print *,'Coordinate matrices need exactly 3 size descriptors.'
      print *, count,' were found.'
      stop
 5000 print *,'Invalid matrix header: ',tmp1
      print *,'Correct header format:'
      print *,'%%MatrixMarket type representation field symmetry'
      print *
      print *,'Check specification and try again.'
      stop
 6000 print *,'''',rep,''' representation not recognized.'
      print *, 'Recognized representations:'
      print *, '   array'
      print *, '   coordinate'
      stop
 7000 print *,'''',field,''' field is not recognized.'
      print *, 'Recognized fields:'
      print *, '   real'
      print *, '   complex'
      print *, '   integer'
      print *, '   pattern'
      stop
 8000 print *,'''',field,''' arrays are not recognized.'
      print *, 'Recognized fields:'
      print *, '   real'
      print *, '   complex'
      print *, '   integer'
      stop
 9000 print *,'''',symm,''' symmetry is not recognized.'
      print *, 'Recognized symmetries:'
      print *, '   general'
      print *, '   symmetric'
      print *, '   hermitian'
      print *, '   skew-symmetric'
      stop
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! End of subroutine mmread
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine mmwrite(ounit,rep,field,symm,rows,cols,nnz,
     *                    indx,jndx,ival,rval,cval)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! This routine will write data to a matrix market formatted file.
! The data may be either sparse coordinate format, or dense array format.
!
! The unit ounit must be open.
!
! 20-Sept-96  Karin A. Remington, NIST ACMD (karin@cam.nist.gov)
! 18-Oct-96   Change in routine name to match C and Matlab routines.
! 30-Oct-96   Bug fixes in mmio.f:
!                  -looping for comment lines
!                  -fixed non-ansi zero stringlength
!                  -incorrect size calculation for skew-symmetric arrays
! 	      Other changes in mmio.f:
!                  -added integer value parameter to calling sequences
!                  -enforced proper count in size info line
!                  -added routine to count words in string (countwd)
!            (Thanks to G.P.Leendetse and H.Oudshoom for their review
!             of the initial version and suggested fixes.)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!   Arguments:
!
!   name     type      in/out description
!   ---------------------------------------------------------------
!
!   ounit  integer     in   Unit identifier for the file
!                             to which the data will be written.
!                             Must be open prior to call.
!
!   rep     character*   in   Matrix Market 'representation'
!                             indicator. Valid inputs:
!
!                                coordinate   (for sparse data)
!                                array        (for dense data)
!                               *elemental*    (to be added)
!
!   field   character*   in   Matrix Market 'field'. Valid inputs:
!
!                                real
!                                complex
!                                integer
!                                pattern (not valid for dense arrays)
!
!   symm    character*   in   Matrix Market 'field'. Valid inputs:
!
!                                symmetric
!                                hermitian
!                                skew-symmetric
!                                general
!
!   rows     integer     in   Number of rows in matrix.
!
!   cols     integer     in   Number of columns in matrix.
!
!   nnz      integer     in   Number of nonzero entries in matrix.
!                             (rows*cols for array matrices)
!
!   indx     integer(nnz)in   Row indices for coordinate format.
!                             Undefined for array format.
!
!   jndx     integer(nnz)in   Column indices for coordinate format.
!                             Undefined for array format.
!
!   ival     integer(nnz) in  Integer data (if applicable, see 'field')
!
!   rval     double(nnz) in   Real data (if applicable, see 'field')
!
!   cval     complex(nnz)in   Complex data (if applicable, see 'field')
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! Declarations:
!
      integer ival(*)
      double precision rval(*)
      complex cval(*)
      integer indx(*)
      integer jndx(*)
      integer i, rows, cols, nnz, nnzreq, ounit
      character*(*)rep,field,symm
!
! Test input qualifiers:
!
      if (rep .ne. 'coordinate' .and. rep .ne. 'array' )
     *   go to 1000
      if (rep .eq. 'coordinate' .and. field .ne. 'integer' .and.
     *    field .ne. 'real' .and. field .ne. 'complex' .and.
     *    field .ne. 'pattern') go to 2000
      if (rep .eq. 'array' .and. field .ne. 'integer' .and.
     *    field .ne. 'real' .and. field .ne. 'complex' ) go to 3000
      if (symm .ne. 'general' .and. symm .ne. 'symmetric' .and.
     *    symm .ne. 'hermitian' .and. symm .ne. 'skew-symmetric')
     *   go to 4000
!
! Write header line:
!
      write(unit=ounit,fmt=5)rep,' ',field,' ',symm
 5    format('%%MatrixMarket matrix ',11A,1A,8A,1A,20A)
!
! Write size information:
!
      if ( rep .eq. 'coordinate' ) then
         nnzreq=nnz
         write(unit=ounit,fmt=*) rows,cols,nnz
         if ( field .eq. 'integer' ) then
            do 10 i=1,nnzreq
               write(unit=ounit,fmt=*)indx(i),jndx(i),ival(i)
 10         continue
         elseif ( field .eq. 'real' ) then
            do 20 i=1,nnzreq
               write(unit=ounit,fmt=*)indx(i),jndx(i),rval(i)
 20         continue
         elseif ( field .eq. 'complex' ) then
            do 30 i=1,nnzreq
               write(unit=ounit,fmt=*)indx(i),jndx(i),
     *                                real(cval(i)),aimag(cval(i))
 30         continue
         else
!        field .eq. 'pattern'
            do 40 i=1,nnzreq
               write(unit=ounit,fmt=*)indx(i),jndx(i)
 40         continue
         endif
      else
!        rep .eq. 'array'
         if ( symm .eq. 'general' ) then
           nnzreq = rows*cols
         elseif ( symm .eq. 'symmetric' .or.
     *            symm .eq. 'hermitian' ) then
           nnzreq = (rows*cols - rows)/2 + rows
         else
!        symm .eq. 'skew-symmetric'
           nnzreq = (rows*cols - rows)/2
         endif
         write(unit=ounit,fmt=*)rows,cols
         if ( field .eq. 'integer' ) then
            do 50 i=1,nnzreq
               write(unit=ounit,fmt=*)ival(i)
 50         continue
         elseif ( field .eq. 'real' ) then
            do 60 i=1,nnzreq
               write(unit=ounit,fmt=*)rval(i)
 60         continue
         else
!        field .eq. 'complex'
            do 70 i=1,nnzreq
               write(unit=ounit,fmt=*)real(cval(i)),aimag(cval(i))
 70         continue
         endif
      endif
      return
!
! Various errors
!
 1000 print *,'''',rep,''' representation not recognized.'
      print *, 'Recognized representations:'
      print *, '   array'
      print *, '   coordinate'
      stop
 2000 print *,'''',field,''' field is not recognized.'
      print *, 'Recognized fields:'
      print *, '   real'
      print *, '   complex'
      print *, '   integer'
      print *, '   pattern'
      stop
 3000 print *,'''',field,''' arrays are not recognized.'
      print *, 'Recognized fields:'
      print *, '   real'
      print *, '   complex'
      print *, '   integer'
      stop
 4000 print *,'''',symm,''' symmetry is not recognized.'
      print *, 'Recognized symmetries:'
      print *, '   general'
      print *, '   symmetric'
      print *, '   hermitian'
      print *, '   skew-symmetric'
      stop
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! End of subroutine mmwrite
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine lowerc(string,pos,len)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Convert uppercase letters to lowercase letters in string with
! starting position pos and length len.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer pos, len
      character*(*) string

      character*26 lcase, ucase
      save lcase,ucase
      data lcase/'abcdefghijklmnopqrstuvwxyz'/
      data ucase/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      do 10 i=pos,len
        k = index(ucase,string(i:i))
        if (k.ne.0) string(i:i) = lcase(k:k)
 10   continue
      return
      end

      subroutine getwd(word,string,slen,start,next,wlen)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Getwd extracts the first  word from string starting
!     at position start.  On return, next is the position
!     of the blank which terminates the word in string.
!     If the found word is longer than the allocated space
!     for the word in the calling program, the word will be
!     truncated to fit.
!     Count is set to the length of the word found.
!
! 30-Oct-96   Bug fix: fixed non-ansi zero stringlength
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer slen, start, next, begin, space, wlen
      character*(*) word
      character*(*) string

      begin = start
      do 5 i=start,slen
         space = index(string(i:slen),' ')
         if ( space .gt. 1) then
            next = i+space-1
            go to 100
         endif
         begin=begin+1
 5    continue
 100  continue
      wlen=next-begin
      if ( wlen .le. 0 ) then
        wlen = 0
        word = ' '
        return
      endif
      word=string(begin:begin+wlen)
      return
      end

      subroutine countwd(string,slen,start,count)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Countwd counts the number of words in string starting
!     at position start.  On return, count is the number of words.
! 30-Oct-96   Routine added
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      character*(*) string
      integer slen, start, next, wordlength, count
      character tmp2*2

      count = 0
      next = 1
 10   call getwd(tmp2,string,1024,next,next,wordlength)
      if ( wordlength .gt. 0 ) then
         count = count + 1
         go to 10
      endif
      return
      end
