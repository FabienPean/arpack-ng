      complex function ccdotc(n,zx,incx,zy,incy)
!
!     forms the dot product of a vector.
!     jack dongarra, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      complex zx(*),zy(*),ztemp
      integer i,incx,incy,ix,iy,n
      ztemp = (0.0d0,0.0d0)
      ccdotc = (0.0d0,0.0d0)
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp = ztemp + conjg(zx(ix))*zy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ccdotc = ztemp
      return
!
!        code for both increments equal to 1
!
   20 do 30 i = 1,n
        ztemp = ztemp + conjg(zx(i))*zy(i)
   30 continue
      ccdotc = ztemp
      return
      end
