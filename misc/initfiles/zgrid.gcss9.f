      program main
c
c   compute zt and zw (layer mid points and boundaries) using
c   sin^2 stretching in BL, fixed dz above, increasing aloft 
c
      integer, parameter :: nz = 96

      real zt(nz)
      real zw(nz+1)

c     ... input parameters

c     ... top of domain (m)
      domain_H = 1500.

c     ... minimum dz (m)
      dzmin = 5.0 

c     ... thickness of layers in normalized (bar) units
      dz = domain_H / nz

c     ... heights of dz minimum aloft and overlying dz minimum (m)
      zi = 800.
      zf = 925.

c     ... possibly adjust dzmin to fit evenly between zi and zf
      if ( zi < zf ) then
        nzmin = nint((zf - zi)/dzmin)
        dzmin = (zf - zi)/nzmin
      endif

c     ... fraction of grid below zi
      fzi = 0.55

c     ... indices of zi and zf
      ki = nint(fzi*nz)
      kf = ki + (zf - zi)/dzmin

c     ... normalized zi and zf 
      zibar = dz*ki
      zfbar = dz*kf

c     ... ratio of minimum dz to H/nz
      fmin = dzmin / dz

c     ... derived parameters
      const_PI = 2*acos(0.)

      b = (2.0)*( (zi/zibar) - fmin )
      d = (3.0/((domain_H - zfbar)**2) )*(
     &     (domain_H - zf)/(domain_H - zfbar) - fmin )

c     ... impossible grids
      if ( kf < ki .or. kf > nz+1 ) then
         print *, 'Error :: ki, kf, nz = ', ki, kf, nz
         stop
      endif
      if ( b < 0 ) then
         print *, 'Error :: b less than zero, b = ', b
         stop
      endif
      if ( d < 0 ) then
         print *, 'Error :: d less than zero, d = ', d
         stop
      endif

c     ... vertical positions of cell boundaries
      do k = 1, nz+1
         zbar = dz*(k-1)
         if ( k <= ki+1 ) then
            u  = const_PI*zbar/zibar
            zz = fmin*zbar+ (b*zibar/const_PI)*
     &            (u/2 - sin( 2.0*u )/4.0 )
         elseif ( k <= kf+1 ) then
            zz = zz + dzmin
         else
            dd = zbar - zfbar
            zz = zf + fmin*dd + d*dd*dd*dd/3.0
         endif
         zw(k) = zz
      enddo

c     ... layer mid-points
      do k = nz, 1, -1
         zt(k)  = 0.5*zw(k) + 0.5*zw(k+1)
      enddo

c     ... print them out
      write(*,'("  k     zw     zt    dz")')
      write(*,'(i4,2f7.1)') nz+1, zw(nz+1)
      do k = nz, 1, -1
         write(*,'(i4,3f7.1)') k, zw(k), zt(k), zw(k+1) - zw(k)
      enddo
      write(*,'("  k     zw     zt    dz")')


c     ... prepare files for the UCLA LES input
      open (1,file='zm_grid_in',status='new',form='formatted')
      do k = 1, nz+1
        write(1,'(f7.1)') zw(k)
      end do
      close (1)

      open (2,file='zt_grid_in',status='new',form='formatted')
      zt0 = zt(1) - (zt(2)-zt(1))
      write(2,'(f7.1)') zt0
      do k = 1, nz
        write(2,'(f7.1)') zt(k)
      end do
      close (2)

      end
