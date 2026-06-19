module ellip

USE decomp_2d_constants
USE decomp_2d_mpi
USE decomp_2d
USE variables
USE param

IMPLICIT NONE

integer :: FS
character(len=100) :: fileformat
character(len=1),parameter :: NL=char(10) !new line character

PRIVATE ! All functions/subroutines private by default
PUBLIC :: init_ellip, boundary_conditions_ellip, postprocess_ellip, &
            geomcomplex_ellip, visu_ellip, visu_ellip_init, update_ellipsoid, &
            check_body_proximity, update_ellipsoid_cv, set_ellipsoid_cv_bounds, &
            init_body_dat, ellipsoid_bc_diagnostic, ellipsoid_pressure_correction_diagnostic, &
            ellipsoid_projection_slip_correction

contains

subroutine geomcomplex_ellip(epsi,nxi,nxf,ny,nyi,nyf,nzi,nzf,dx,yp,dz,remp)

    use param, only : one, two, ten
    use ibm_param
    use ellipsoid_utils, only: NormalizeQuaternion, EllipsoidalRadius, EllipsoidalRadius_debug
    use complex_geometry, only: nraf,nyraf

    implicit none

    integer                    :: nxi,nxf,ny,nyi,nyf,nzi,nzf
    real(mytype),dimension(nxi:nxf,nyi:nyf,nzi:nzf) :: epsi
    real(mytype),dimension(ny) :: yp
    real(mytype)               :: dx,dz
    real(mytype)               :: remp
    integer                    :: i,j,k, i_body
    real(mytype)               :: xm,ym,zm,r,rads2,kcon
    real(mytype)               :: zeromach
    real(mytype)               :: cexx,ceyy,cezz,dist_axi
    real(mytype)               :: point(3)
    logical                    :: is_inside

    zeromach=one
    do while ((one + zeromach / two) .gt. one)
        zeromach = zeromach/two
    end do
    zeromach = ten*zeromach
    is_inside=.false.
    !  orientation=[oriw, orii, orij, orik]
    do i = 1,nbody
        call NormalizeQuaternion(orientation(i,:))
    enddo
    !  shape=[shx, shy, shz]
    !  write(*,*) shape, 'SHAPE'


    ! Intitialise epsi
    epsi(:,:,:)=zero



    ! Update center of moving ellipsoid
    ! if (t.ne.0.) then
    !     cexx=cex+lvx*(t-ifirst*dt)
    !     ceyy=cey+lvy*(t-ifirst*dt)
    !     cezz=cez+lvz*(t-ifirst*dt)
    ! else
    !     cexx=cex
    !     ceyy=cey
    !     cezz=cez
    ! endif
    ! position=[cexx,ceyy,cezz]
    !  write(*,*) position
    !  ce=[cexx, ceyy, cezz]
    !
    ! Define adjusted smoothing constant
!    kcon = log((one-0.0001)/0.0001)/(smoopar*0.5*dx) ! 0.0001 is the y-value, smoopar: desired number of affected points
!   write(*,*) nzi, nzf
    do k=nzi,nzf
    zm=(real(k-1,mytype))*dz
    ! write(*,*) k, zm
        do j=nyi,nyf
        ! ym=(real(j-1,mytype))*dy
            ym=yp(j)
            do i=nxi,nxf
               xm=real(i-1,mytype)*dx
               point=[xm, ym, zm]
            ! call EllipsoidalRadius(point, position, orientation, shape, r)
               do i_body = 1,nbody
                if (cube_flag.eq.0) then
                    call EllipsoidalRadius(point,position(i_body,:),orientation(i_body,:),shape(i_body,:),r)
                    is_inside = (r-ra(i_body)).lt.zeromach
                    ! if (is_inside) then
                    !     call EllipsoidalRadius_debug(point,position(i_body,:),orientation(i_body,:),shape(i_body,:),r)
                    ! endif
                    if (ra(i_body) /= ra(i_body)) then
                        write(*,*) "Nrank = ", nrank
                        write(*,*) "Point = ", point
                    endif
                else if (cube_flag.eq.1) then
                    is_inside = (abs(xm-position(i_body,1)).lt.ra(i_body)).and.(abs(ym-position(i_body,2)).lt.ra(i_body)).and.(abs(zm-position(i_body,3)).lt.ra(i_body))
                endif
                !  r=sqrt_prec((xm-cexx)**two+(ym-ceyy)**two+(zm-cezz)**two)
                !  r=sqrt_prec((xm-cexx)**two+(ym-ceyy)**two)
                if (is_inside) then
                    !  write(*,*) i, j, k
                    epsi(i,j,k)=remp
                    cycle
                endif
            enddo
            ! write(*,*) is_inside

            !  write(*,*) i, j, k, zm
            ! epsi(i,j,k)=remp
            !  write(*,*) remp
        enddo
        enddo
    enddo

    return
end subroutine geomcomplex_ellip

!********************************************************************
subroutine boundary_conditions_ellip (ux,uy,uz,phi)

    USE param
    USE variables
    USE decomp_2d

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    call inflow (phi)
    call outflow (ux,uy,uz,phi)

    return
end subroutine boundary_conditions_ellip
!********************************************************************
subroutine inflow (phi)

    USE param
    USE variables
    USE decomp_2d
    USE ibm_param

    implicit none

    integer  :: i,j,k,is
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    if ((shear_flow_ybc.eq.1).or.(shear_flow_zbc.eq.1)) then
        u1 = 0.0_mytype
        u2 = 0.0_mytype
    endif


    !call random_number(bxo)
    !call random_number(byo)
    !call random_number(bzo)
    do k=1,xsize(3)
        do j=1,xsize(2)
        bxx1(j,k)=u1+bxo(j,k)*inflow_noise
        bxy1(j,k)=zero+byo(j,k)*inflow_noise
        bxz1(j,k)=zero+bzo(j,k)*inflow_noise
        enddo
    enddo

    if (shear_flow_ybc.eq.1) then
        do k=1,xsize(3)
            do i=1,xsize(1)
                byxn(i,k)=+shear_velocity
            enddo
        enddo
        do k=1,xsize(3)
            do i=1,xsize(1)
                byx1(i,k)=-shear_velocity
            enddo
        enddo
    endif

    if (shear_flow_zbc.eq.1) then
        do j=1,xsize(2)
            do i=1,xsize(1)
                bzxn(i,j)=+shear_velocity
            enddo
        enddo
        do j=1,xsize(2)
            do i=1,xsize(1)
                bzx1(i,j)=-shear_velocity
            enddo
        enddo
    endif


    if (iscalar.eq.1) then
        do is=1, numscalar
        do k=1,xsize(3)
            do j=1,xsize(2)
                phi(1,j,k,is)=cp(is)
            enddo
        enddo
        enddo
    endif

    return
end subroutine inflow
!********************************************************************
subroutine outflow (ux,uy,uz,phi)

    USE param
    USE variables
    USE decomp_2d
    USE MPI
    USE ibm_param

    implicit none

    integer :: j,k,code
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cx,uxmin,uxmax,uxmin1,uxmax1

    udx=one/dx; udy=one/dy; udz=one/dz; uddx=half/dx; uddy=half/dy; uddz=half/dz

    uxmax=-1609._mytype
    uxmin=1609._mytype
    do k=1,xsize(3)
        do j=1,xsize(2)
        if (ux(nx-1,j,k).gt.uxmax) uxmax=ux(nx-1,j,k)
        if (ux(nx-1,j,k).lt.uxmin) uxmin=ux(nx-1,j,k)
        enddo
    enddo

    call MPI_ALLREDUCE(uxmax,uxmax1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(uxmin,uxmin1,1,real_type,MPI_MIN,MPI_COMM_WORLD,code)

    if (u1 == zero) then
        cx=(half*(uxmax1+uxmin1))*gdt(itr)*udx
    elseif (u1 == one) then
        cx=uxmax1*gdt(itr)*udx
    elseif (u1 == two) then
        cx=u2*gdt(itr)*udx    !works better
    else
        cx=(half*(u1+u2))*gdt(itr)*udx
    endif

    do k=1,xsize(3)
        do j=1,xsize(2)
        bxxn(j,k)=ux(nx,j,k)-cx*(ux(nx,j,k)-ux(nx-1,j,k))
        bxyn(j,k)=uy(nx,j,k)-cx*(uy(nx,j,k)-uy(nx-1,j,k))
        bxzn(j,k)=uz(nx,j,k)-cx*(uz(nx,j,k)-uz(nx-1,j,k))
        enddo
    enddo

    if (iscalar==1) then
        if (u2==zero) then
        cx=(half*(uxmax1+uxmin1))*gdt(itr)*udx
        elseif (u2==one) then
        cx=uxmax1*gdt(itr)*udx
        elseif (u2==two) then
        cx=u2*gdt(itr)*udx    !works better
        else
        stop
        endif

        do k=1,xsize(3)
        do j=1,xsize(2)
            phi(nx,j,k,:)=phi(nx,j,k,:)-cx*(phi(nx,j,k,:)-phi(nx-1,j,k,:))
        enddo
        enddo
    endif

    if (nrank==0.and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime == ilast)) &
        write(*,*) "Outflow velocity ux nx=n min max=",real(uxmin1,4),real(uxmax1,4)

    return
end subroutine outflow
!********************************************************************
subroutine init_ellip (ux1,uy1,uz1,phi1)

    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI
    USE ibm_param
    use ellipsoid_utils, only: NormalizeQuaternion,ellipInertiaCalculate,ellipMassCalculate


    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype) :: y,um,eqr,ym
    integer :: k,j,i,ii,is,code,jj

    ! write(*,*) 'INSIDE INIT ELLIP'

    ! eqr=(sh(1)*sh(2)*sh(3))**(1.0/3.0)
    ! shape=sh(:)/eqr

    ! orientation=ori
    ! call NormalizeQuaternion(orientation)
    ! position=ce
    ! linearVelocity=lv
    ! angularVelocity=[zero, av(1), av(2), av(3)]
    ! call ellipInertiaCalculate(shape,rho_s,inertia)
    ! call ellipMassCalculate(shape,rho_s,ellip_m)

    ! if (nrank==0) then
    !     write(*,*) 'set shape             = ', shape
    !     write(*,*) 'set orientation       = ', orientation
    !     write(*,*) 'set position          = ', position
    !     write(*,*) 'set linear velocity   = ', linearVelocity
    !     write(*,*) 'set angular velocity  = ', angularVelocity
    !     write(*,*) 'set moment of inertia = ', inertia
    !     write(*,*) 'density of solid      = ', rho_s
    ! end if

    if (iscalar==1) then

        phi1(:,:,:,:) = zero !change as much as you want

    endif
    ! if (shear_flow_ybc.eq.1) then
    !     do i=1,xsize(1)
    !         do j=1,xsize(2)
    !             jj=j+xstart(2)-1
    !             ym=real(jj)*dy
    !             do k=1,xsize(3)
    !                 ux1(i,j,k)=real((jj-(ny/2)))/(yly/2.0)*shear_velocity
    !             enddo
    !         enddo
    !     enddo
    ! else
        ux1=zero;
    ! endif


    uy1=zero; uz1=zero

    if (iin.ne.0) then
        call system_clock(count=code)
        if (iin.eq.2) code=0
        call random_seed(size = ii)
        call random_seed(put = code+63946*(nrank+1)*(/ (i - 1, i = 1, ii) /))

        call random_number(ux1)
        call random_number(uy1)
        call random_number(uz1)

        do k=1,xsize(3)
        do j=1,xsize(2)
            do i=1,xsize(1)
                ux1(i,j,k)=init_noise*(ux1(i,j,k)-0.5)
                uy1(i,j,k)=init_noise*(uy1(i,j,k)-0.5)
                uz1(i,j,k)=init_noise*(uz1(i,j,k)-0.5)
            enddo
        enddo
        enddo

        !modulation of the random noise
        do k=1,xsize(3)
        do j=1,xsize(2)
            if (istret.eq.0) y=(j+xstart(2)-1-1)*dy-yly/2.
            if (istret.ne.0) y=yp(j+xstart(2)-1)-yly/2.
            um=exp(-zptwo*y*y)
            do i=1,xsize(1)
                ux1(i,j,k)=um*ux1(i,j,k)
                uy1(i,j,k)=um*uy1(i,j,k)
                uz1(i,j,k)=um*uz1(i,j,k)
            enddo
        enddo
        enddo
    endif

    !INIT FOR G AND U=MEAN FLOW + NOISE
    do k=1,xsize(3)
        do j=1,xsize(2)
        do i=1,xsize(1)
            ux1(i,j,k)=ux1(i,j,k)+u1
            uy1(i,j,k)=uy1(i,j,k)
            if (shear_flow_ybc.eq.1) then
                ux1(i,j,k)=ux1(i,j,k)+((j+xstart(2)-1-1)*dy-yly/2.)/(yly/2.0)*shear_velocity
            endif
            if (shear_flow_zbc.eq.1) then
                ux1(i,j,k)=ux1(i,j,k)+((k+xstart(3)-1-1)*dz-zlz/2.)/(zlz/2.0)*shear_velocity
            endif
            uz1(i,j,k)=uz1(i,j,k)
        enddo
        enddo
    enddo

    if (ellipsoid_init_potential.eq.1) then
        call init_potential_sphere(ux1, uy1, uz1)
    endif

#ifdef DEBG
    if (nrank .eq. 0) write(*,*) '# init end ok'
#endif

    call init_body_dat()

    return
end subroutine init_ellip

subroutine sphere_potential_velocity(point, velocity)
    !! Exact unbounded-domain inviscid velocity around a single sphere.

    use param
    use ibm_param

    implicit none

    real(mytype), dimension(3), intent(in) :: point
    real(mytype), dimension(3), intent(out) :: velocity
    real(mytype) :: rel(3), ufar(3), urel(3), ub(3)
    real(mytype) :: radius, r2, rmag, r3, r5, udotr, correction(3)

    radius = ra(1) * shape(1,1)
    ufar = [u1, zero, zero]
    ub = linearVelocity(1,:)
    urel = ufar - ub

    rel = point - position(1,:)
    r2 = sum(rel * rel)
    rmag = sqrt(r2)

    if (rmag.le.radius) then
        velocity = ub
    else
        r3 = r2 * rmag
        r5 = r3 * r2
        udotr = sum(urel * rel)
        correction = (radius**3 / two) * (urel / r3 - three * udotr * rel / r5)
        velocity = ub + urel + correction
    endif
end subroutine sphere_potential_velocity

subroutine init_potential_sphere(ux1, uy1, uz1)
    !! Initialise the exact unbounded-domain potential-flow velocity around a
    !! single sphere in a uniform stream. This is intended for inviscid
    !! no-penetration validation; it is not the ellipsoid potential solution.

    use decomp_2d
    use variables
    use param
    use ibm_param

    implicit none

    real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(out) :: ux1, uy1, uz1
    real(mytype) :: point(3), velocity(3), ufar(3), ub(3)
    real(mytype) :: radius
    real(mytype) :: shape_tol
    integer :: i, j, k

    if (nbody.ne.1 .and. nrank.eq.0) then
        write(*,*) "WARNING: ellipsoid_init_potential currently uses body 1 only."
    endif

    shape_tol = 1.0e-10_mytype
    if ((abs(shape(1,1)-shape(1,2)).gt.shape_tol .or. &
         abs(shape(1,1)-shape(1,3)).gt.shape_tol) .and. nrank.eq.0) then
        write(*,*) "WARNING: ellipsoid_init_potential is analytic only for a sphere."
        write(*,*) "         Current body 1 shape = ", shape(1,:)
    endif

    radius = ra(1) * shape(1,1)
    ufar = [u1, zero, zero]
    ub = linearVelocity(1,:)

    if (nrank.eq.0) then
        write(*,*) "Initialising analytic inviscid sphere potential flow"
        write(*,*) "  centre = ", position(1,:)
        write(*,*) "  radius = ", radius
        write(*,*) "  ufar   = ", ufar
        write(*,*) "  ub     = ", ub
    endif

    do k=1,xsize(3)
        point(3) = real(k + xstart(3) - 2, mytype) * dz
        do j=1,xsize(2)
            if (istret.eq.0) then
                point(2) = real(j + xstart(2) - 2, mytype) * dy
            else
                point(2) = yp(j + xstart(2) - 1)
            endif
            do i=1,xsize(1)
                point(1) = real(i + xstart(1) - 2, mytype) * dx
                call sphere_potential_velocity(point, velocity)
                ux1(i,j,k) = velocity(1)
                uy1(i,j,k) = velocity(2)
                uz1(i,j,k) = velocity(3)
            enddo
        enddo
    enddo
end subroutine init_potential_sphere

!############################################################################
subroutine init_body_dat()
  !! Open body.dat{N} output files. On a fresh run they are truncated;
  !! on restart they are trimmed back to t0 so re-restarting from the
  !! same checkpoint never leaves stale or duplicate rows.
  use ibm_param, only : nbody
  implicit none
  integer :: i
  character(len=30) :: filename
  if (nrank /= 0) return
  do i = 1, nbody
     write(filename,"('body.dat',I1.1)") i
     if (irestart == 0) then
        open(unit=11+i, file=filename, status='replace', form='formatted')
     else
        call open_body_dat(11+i, filename, t0)
     endif
  enddo
end subroutine init_body_dat

!############################################################################
subroutine open_body_dat(iunit, filename, t_restart)
  !! Open body.dat for append on restart, trimming any entries with
  !! t > t_restart so that re-restarting from the same checkpoint never
  !! leaves stale or duplicate data in the file.
  implicit none
  integer,          intent(in) :: iunit
  character(len=*), intent(in) :: filename
  real(mytype),     intent(in) :: t_restart

  integer            :: ios, tmp_unit
  character(len=512) :: line
  real(mytype)       :: t_val
  logical            :: exists
  character(len=64)  :: tmpfile

  tmp_unit = iunit + 50

  inquire(file=filename, exist=exists)
  if (.not. exists) then
     open(unit=iunit, file=filename, status='new', form='formatted')
     return
  endif

  ! Pass 1: copy lines with t <= t_restart to a temp file
  write(tmpfile,"('body_tmp.dat',I1.1)") iunit - 11
  open(unit=iunit,    file=filename, status='old',     form='formatted', action='read')
  open(unit=tmp_unit, file=tmpfile,  status='replace', form='formatted')
  do
     read(iunit, '(A)', iostat=ios) line
     if (ios /= 0) exit
     read(line, *, iostat=ios) t_val
     if (ios /= 0) cycle
     if (t_val <= t_restart + epsilon(t_restart)) write(tmp_unit, '(A)') trim(line)
  enddo
  close(iunit)
  close(tmp_unit)

  ! Pass 2: overwrite original with trimmed content; leave open at end for append
  open(unit=iunit,    file=filename, status='replace', form='formatted')
  open(unit=tmp_unit, file=tmpfile,  status='old',     form='formatted', action='read')
  do
     read(tmp_unit, '(A)', iostat=ios) line
     if (ios /= 0) exit
     write(iunit, '(A)') trim(line)
  enddo
  close(tmp_unit, status='delete')
  ! iunit remains open, positioned at end of file, ready for append

end subroutine open_body_dat

!********************************************************************

!############################################################################
subroutine postprocess_ellip(ux1,uy1,uz1,ep1)

    USE MPI
    USE decomp_2d
    USE decomp_2d_io
    USE var, only : uvisu
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    USE ibm_param

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1

end subroutine postprocess_ellip

subroutine visu_ellip_init (visu_initialised)

    use decomp_2d
    use decomp_2d_io, only : decomp_2d_register_variable
    use visu, only : io_name, output2D

    implicit none

    logical, intent(out) :: visu_initialised

    call decomp_2d_register_variable(io_name, "vort", 1, 0, output2D, mytype)
    call decomp_2d_register_variable(io_name, "critq", 1, 0, output2D, mytype)

    visu_initialised = .true.

end subroutine visu_ellip_init
!############################################################################
!!
!!  SUBROUTINE: visu_ellip
!!      AUTHOR: FS
!! DESCRIPTION: Performs ellipinder-specific visualization
!!
!############################################################################
subroutine visu_ellip(ux1, uy1, uz1, pp3, phi1, ep1, num)

    use var, only : ux2, uy2, uz2, ux3, uy3, uz3
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    use var, ONLY : nxmsize, nymsize, nzmsize
    use visu, only : write_field
    use ibm_param, only : ubcx,ubcy,ubcz,inviscid_output

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ep1
    integer, intent(in) :: num

    ! Write vorticity as an example of post processing

    ! Perform communications if needed
    if (sync_vel_needed) then
    call transpose_x_to_y(ux1,ux2)
    call transpose_x_to_y(uy1,uy2)
    call transpose_x_to_y(uz1,uz2)
    call transpose_y_to_z(ux2,ux3)
    call transpose_y_to_z(uy2,uy3)
    call transpose_y_to_z(uz2,uz3)
    sync_vel_needed = .false.
    endif

    !x-derivatives
    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,1) !ubcx is 1. etc
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,2)
    call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,3)
    !y-derivatives
    call dery (ta2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,1)
    call dery (tb2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,2)
    call dery (tc2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,3)
    !!z-derivatives
    call derz (ta3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,1)
    call derz (tb3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,2)
    call derz (tc3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,3)
    !!all back to x-pencils
    call transpose_z_to_y(ta3,td2)
    call transpose_z_to_y(tb3,te2)
    call transpose_z_to_y(tc3,tf2)
    call transpose_y_to_x(td2,tg1)
    call transpose_y_to_x(te2,th1)
    call transpose_y_to_x(tf2,ti1)
    call transpose_y_to_x(ta2,td1)
    call transpose_y_to_x(tb2,te1)
    call transpose_y_to_x(tc2,tf1)
    !du/dx=ta1 du/dy=td1 and du/dz=tg1
    !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
    !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1
    !VORTICITY FIELD
    di1 = zero
    di1(:,:,:)=sqrt(  (tf1(:,:,:)-th1(:,:,:))**2 &
                    + (tg1(:,:,:)-tc1(:,:,:))**2 &
                    + (tb1(:,:,:)-td1(:,:,:))**2)

    if (inviscid_output.eq.0) then

    call write_field(di1, ".", "vort", num, flush = .true.) ! Reusing temporary array, force flush

    !Q=-0.5*(ta1**2+te1**2+ti1**2)-td1*tb1-tg1*tc1-th1*tf1
    di1 = zero
    di1(:,:,: ) = - half*(ta1(:,:,:)**2+te1(:,:,:)**2+ti1(:,:,:)**2) &
                - td1(:,:,:)*tb1(:,:,:) &
                - tg1(:,:,:)*tc1(:,:,:) &
                - th1(:,:,:)*tf1(:,:,:)
    call write_field(di1, ".", "critq", num, flush = .true.) ! Reusing temporary array, force flush
    endif
end subroutine visu_ellip

subroutine set_ellipsoid_cv_bounds()
  !! Set xld/xrd/yld/yud/zld/zrd from the current body position.
  !! Call this before init_forces() so icvlf is computed from the correct
  !! position (not the input-file ForceCVs defaults).
  use forces, only : iforces, xld, xrd, yld, yud, zld, zrd, nvol
  use ibm_param
  use complex_geometry, only : nobjmax
  implicit none
  real(mytype) :: maxrad
  integer :: i
  do i = 1, min(nobjmax, nvol)
     if (iforces .eq. 1) then
        maxrad = max(shape(i,1), shape(i,2), shape(i,3))
        xld(i) = position(i,1) - maxrad * ra(i) * cvl_scalar
        xrd(i) = position(i,1) + maxrad * ra(i) * cvl_scalar
        yld(i) = position(i,2) - maxrad * ra(i) * cvl_scalar
        yud(i) = position(i,2) + maxrad * ra(i) * cvl_scalar
        zld(i) = position(i,3) - maxrad * ra(i) * cvl_scalar
        zrd(i) = position(i,3) + maxrad * ra(i) * cvl_scalar
     endif
  enddo
end subroutine set_ellipsoid_cv_bounds

subroutine update_ellipsoid_cv()
  !! Update the CV bounds (xld/xrd/etc.) from the current body position/shape,
  !! check for boundary violations, run proximity detection, and recompute
  !! integer CV indices via init_forces/update_forces.
  use forces, only : iforces, init_forces, update_forces, xld, xrd, yld, yud, zld, zrd, nvol
  use ibm_param
  use complex_geometry, only : nobjmax
  use MPI
  implicit none

  real(mytype) :: maxrad
  integer :: i, code, ierror

  do i = 1, min(nobjmax, nvol)
     maxrad = max(shape(i,1), shape(i,2), shape(i,3))
     if (iforces .eq. 1) then
        xld(i) = position(i,1) - maxrad * ra(i) * cvl_scalar
        xrd(i) = position(i,1) + maxrad * ra(i) * cvl_scalar
        yld(i) = position(i,2) - maxrad * ra(i) * cvl_scalar
        yud(i) = position(i,2) + maxrad * ra(i) * cvl_scalar
        zld(i) = position(i,3) - maxrad * ra(i) * cvl_scalar
        zrd(i) = position(i,3) + maxrad * ra(i) * cvl_scalar
        if ((xld(i).lt.0).or.(xrd(i).gt.xlx).or.(yld(i).lt.0).or.(yud(i).gt.yly)) then
           write(*,*) "Body is too close to boundary!"
           call MPI_ABORT(MPI_COMM_WORLD, code, ierror)
        endif
        if ((zld(i).lt.0).or.(zrd(i).gt.zlz)) then
           write(*,*) "Body is too close to boundary!"
           call MPI_ABORT(MPI_COMM_WORLD, code, ierror)
        endif
     endif
  enddo

  call check_body_proximity()

  if (iforces .eq. 1) then
     if (itime.eq.ifirst .and. irestart.eq.0) then
        call init_forces()
     else
        call update_forces()
     endif
  endif

end subroutine update_ellipsoid_cv

subroutine update_ellipsoid(ux1, uy1, uz1, ep1)

    use forces, only : force, torque_calc, nvol, iforces
    use ellipsoid_utils, only : lin_step, ang_step
    use ibm_param
    use param, only : zero, one, dt, dx, dy, dz, zpfive
    use variables, only : ilist
    use var, only : itime, t
    use decomp_2d_mpi, only : nrank
    use MPI

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1

    real(mytype) :: drag(10), lift(10), lat(10)
    real(mytype) :: grav_effx(10), grav_effy(10), grav_effz(10)
    real(mytype) :: xtorq(10), ytorq(10), ztorq(10)
    real(mytype) :: eek
    integer :: i, code

    ! Body dynamics need forces; if force calculation is disabled, skip the
    ! whole update (body stays at initial position/velocity).
    if (iforces.ne.1) return

    xtorq = zero; ytorq = zero; ztorq = zero

    call force(ux1, uy1, uz1, ep1, drag, lift, lat, 1)

    grav_effx = grav_x * (rho_s - 1.0_mytype)
    grav_effy = grav_y * (rho_s - 1.0_mytype)
    grav_effz = grav_z * (rho_s - 1.0_mytype)
    do i = 1, nbody
       linearForce(i,:) = [drag(i)-grav_effx(i), lift(i)-grav_effy(i), lat(i)-grav_effz(i)]
    enddo

    if (torques_flag.eq.1) then
       call torque_calc(ux1, uy1, uz1, ep1, xtorq, ytorq, ztorq, 1)
    endif
    if (orientations_free.eq.1) then
       do i = 1, nvol
          torque(i,:) = [xtorq(i), ytorq(i), ztorq(i)]
       enddo
       if (ztorq_only.eq.1) then
          torque(:,1) = zero
          torque(:,2) = zero
       endif
    else
       torque(:,:) = zero
    endif

    if (bodies_fixed.ne.1) then
       do i = 1, nvol
          call lin_step(position(i,:), linearVelocity(i,:), linearForce(i,:), ellip_m(i), ellip_m_added(i,:), orientation(i,:), dt, position_1, linearVelocity_1)
          call ang_step(orientation(i,:), angularVelocity(i,:), torque(i,:), inertia(i,:,:), inertia_rot_added(i,:), dt, orientation_1, angularVelocity_1)
          position(i,:) = position_1
          linearVelocity(i,:) = linearVelocity_1
          orientation(i,:) = orientation_1
          angularVelocity(i,:) = angularVelocity_1
       enddo
    endif

    if (mod(itime,ilist)==0) then
       ! All ranks contribute to the masked fluid kinetic energy
       eek = sum(zpfive * (one - ep1) * (ux1**2 + uy1**2 + uz1**2)) * dx * dy * dz
       call MPI_Allreduce(MPI_IN_PLACE, eek, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
       call ellipsoid_bc_diagnostic(ux1, uy1, uz1)

       if (nrank==0) then
          write(*,*) "Kinetic Energy =   ", eek
          do i = 1, nbody
             write(11+i,*) t, position(i,1), position(i,2), position(i,3), &
                  orientation(i,1), orientation(i,2), orientation(i,3), orientation(i,4), &
                  linearVelocity(i,1), linearVelocity(i,2), linearVelocity(i,3), &
                  angularVelocity(i,2), angularVelocity(i,3), angularVelocity(i,4), &
                  linearForce(i,1), linearForce(i,2), linearForce(i,3), &
                  torque(i,1), torque(i,2), torque(i,3), &
                  eek
             flush(11+i)
             write(*,*) "Body", i
             write(*,*) "Position =         ", position(i,:)
             write(*,*) "Orientation =      ", orientation(i,:)
             write(*,*) "Linear velocity =  ", linearVelocity(i,:)
             write(*,*) "Angular velocity = ", angularVelocity(i,:)
             write(*,*) "Linear Force =     ", linearForce(i,:)
             write(*,*) "Torque =           ", torque(i,:)
          enddo
       endif
    endif

end subroutine update_ellipsoid

!********************************************************************
subroutine ellipsoid_bc_diagnostic(ux1, uy1, uz1, stage)

    use complex_geometry, only : nobjx, nobjy, nobjz, xi, xf, yi, yf, zi, zf
    use param, only : zero, dx, dy, dz, xnu, xlx, yly, zlz, izap
    use variables, only : yp
    use var, only : ux2, uy2, uz2, ux3, uy3, uz3, t
    use decomp_2d_mpi, only : nrank
    use MPI

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    character(len=*), intent(in), optional :: stage

    integer :: i, j, k, iobj, ix, jy, kz, ig, jg, kg, code, iunit
    integer :: local_count, global_count
    logical :: file_exists
    logical :: projection_stage
    real(mytype) :: xm, ym, zm
    real(mytype) :: point(3), nearestVelocity(3)
    real(mytype) :: local_normal_max, local_full_max
    real(mytype) :: local_normal_sum2, local_full_sum2
    real(mytype) :: global_normal_max, global_full_max
    real(mytype) :: global_normal_sum2, global_full_sum2
    real(mytype) :: normal_rms, full_rms

    projection_stage = present(stage)
    local_count = 0
    local_normal_max = zero
    local_full_max = zero
    local_normal_sum2 = zero
    local_full_sum2 = zero

    call transpose_x_to_y(ux1, ux2)
    call transpose_x_to_y(uy1, uy2)
    call transpose_x_to_y(uz1, uz2)
    call transpose_y_to_z(ux2, ux3)
    call transpose_y_to_z(uy2, uy3)
    call transpose_y_to_z(uz2, uz3)

    do k = 1, xsize(3)
       zm = real(xstart(3)+k-1, mytype) * dz
       do j = 1, xsize(2)
          ym = real(xstart(2)+j-1, mytype) * dy
          do iobj = 1, nobjx(j,k)
             if (xi(iobj,j,k) .gt. zero) then
                ix = xi(iobj,j,k) / dx + 1
                if (izap .eq. 1) then
                   nearestVelocity = [ux1(ix-1,j,k), uy1(ix-1,j,k), uz1(ix-1,j,k)]
                else
                   nearestVelocity = [ux1(ix,j,k), uy1(ix,j,k), uz1(ix,j,k)]
                endif
                point = [xi(iobj,j,k), ym, zm]
                call accumulate_ellipsoid_bc_sample(point, nearestVelocity, local_normal_max, &
                     local_normal_sum2, local_full_max, local_full_sum2, local_count)
             endif

             if (xf(iobj,j,k) .lt. xlx) then
                ix = (xf(iobj,j,k) + dx) / dx + 1
                if (izap .eq. 1) then
                   nearestVelocity = [ux1(ix+1,j,k), uy1(ix+1,j,k), uz1(ix+1,j,k)]
                else
                   nearestVelocity = [ux1(ix,j,k), uy1(ix,j,k), uz1(ix,j,k)]
                endif
                point = [xf(iobj,j,k), ym, zm]
                call accumulate_ellipsoid_bc_sample(point, nearestVelocity, local_normal_max, &
                     local_normal_sum2, local_full_max, local_full_sum2, local_count)
             endif
          enddo
       enddo
    enddo

    do k = 1, ysize(3)
       zm = real(ystart(3)+k-1, mytype) * dz
       kg = ystart(3) + k - 1
       do i = 1, ysize(1)
          xm = real(ystart(1)+i-1, mytype) * dx
          ig = ystart(1) + i - 1
          do j = 1, nobjy(i,k)
             if (yi(j,i,k) .gt. zero) then
                jy = 1
                do while (yp(jy) .lt. yi(j,i,k))
                   jy = jy + 1
                enddo
                jy = jy - 1
                if (izap .eq. 1) then
                   nearestVelocity = [ux2(i,jy-1,k), uy2(i,jy-1,k), uz2(i,jy-1,k)]
                else
                   nearestVelocity = [ux2(i,jy,k), uy2(i,jy,k), uz2(i,jy,k)]
                endif
                point = [xm, yi(j,i,k), zm]
                call accumulate_ellipsoid_bc_sample(point, nearestVelocity, local_normal_max, &
                     local_normal_sum2, local_full_max, local_full_sum2, local_count)
             endif

             if (yf(j,i,k) .lt. yly) then
                jy = 1
                do while (yp(jy) .lt. yf(j,i,k))
                   jy = jy + 1
                enddo
                if (izap .eq. 1) then
                   nearestVelocity = [ux2(i,jy+1,k), uy2(i,jy+1,k), uz2(i,jy+1,k)]
                else
                   nearestVelocity = [ux2(i,jy,k), uy2(i,jy,k), uz2(i,jy,k)]
                endif
                point = [xm, yf(j,i,k), zm]
                call accumulate_ellipsoid_bc_sample(point, nearestVelocity, local_normal_max, &
                     local_normal_sum2, local_full_max, local_full_sum2, local_count)
             endif
          enddo
       enddo
    enddo

    do j = 1, zsize(2)
       ym = real(zstart(2)+j-1, mytype) * dy
       jg = zstart(2) + j - 1
       do i = 1, zsize(1)
          xm = real(zstart(1)+i-1, mytype) * dx
          ig = zstart(1) + i - 1
          do k = 1, nobjz(i,j)
             if (zi(k,i,j) .gt. zero) then
                kz = zi(k,i,j) / dz + 1
                if (izap .eq. 1) then
                   nearestVelocity = [ux3(i,j,kz-1), uy3(i,j,kz-1), uz3(i,j,kz-1)]
                else
                   nearestVelocity = [ux3(i,j,kz), uy3(i,j,kz), uz3(i,j,kz)]
                endif
                point = [xm, ym, zi(k,i,j)]
                call accumulate_ellipsoid_bc_sample(point, nearestVelocity, local_normal_max, &
                     local_normal_sum2, local_full_max, local_full_sum2, local_count)
             endif

             if (zf(k,i,j) .lt. zlz) then
                kz = (zf(k,i,j) + dz) / dz + 1
                if (izap .eq. 1) then
                   nearestVelocity = [ux3(i,j,kz+1), uy3(i,j,kz+1), uz3(i,j,kz+1)]
                else
                   nearestVelocity = [ux3(i,j,kz), uy3(i,j,kz), uz3(i,j,kz)]
                endif
                point = [xm, ym, zf(k,i,j)]
                call accumulate_ellipsoid_bc_sample(point, nearestVelocity, local_normal_max, &
                     local_normal_sum2, local_full_max, local_full_sum2, local_count)
             endif
          enddo
       enddo
    enddo

    global_count = local_count
    global_normal_max = local_normal_max
    global_full_max = local_full_max
    global_normal_sum2 = local_normal_sum2
    global_full_sum2 = local_full_sum2

    call MPI_Allreduce(MPI_IN_PLACE, global_count, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, code)
    call MPI_Allreduce(MPI_IN_PLACE, global_normal_max, 1, real_type, MPI_MAX, MPI_COMM_WORLD, code)
    call MPI_Allreduce(MPI_IN_PLACE, global_full_max, 1, real_type, MPI_MAX, MPI_COMM_WORLD, code)
    call MPI_Allreduce(MPI_IN_PLACE, global_normal_sum2, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
    call MPI_Allreduce(MPI_IN_PLACE, global_full_sum2, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)

    normal_rms = zero
    full_rms = zero
    if (global_count .gt. 0) then
       normal_rms = sqrt(global_normal_sum2 / real(global_count, mytype))
       full_rms = sqrt(global_full_sum2 / real(global_count, mytype))
    endif

    if (nrank .eq. 0) then
       if (projection_stage) then
          write(*,*) "Ellipsoid projection BC diagnostic ", trim(stage), &
               ": samples=", global_count, &
               " normal max/rms=", global_normal_max, normal_rms, &
               " full max/rms=", global_full_max, full_rms

          inquire(file="ellipsoid_projection_bc_error.dat", exist=file_exists)
          open(newunit=iunit, file="ellipsoid_projection_bc_error.dat", status="unknown", position="append")
          if (.not. file_exists) then
             write(iunit,*) "# t itime itr stage samples normal_max normal_rms full_max full_rms xnu"
          endif
          write(iunit,*) t, itime, itr, trim(stage), global_count, &
               global_normal_max, normal_rms, global_full_max, full_rms, xnu
          close(iunit)
       else
          if (xnu .eq. zero) then
             write(*,*) "Ellipsoid slip nearest-cell diagnostic: samples=", global_count, &
                  " normal max/rms=", global_normal_max, normal_rms, &
                  " full max/rms=", global_full_max, full_rms
          else
             write(*,*) "Ellipsoid no-slip nearest-cell diagnostic: samples=", global_count, &
                  " full max/rms=", global_full_max, full_rms, &
                  " normal max/rms=", global_normal_max, normal_rms
          endif

          inquire(file="ellipsoid_bc_error.dat", exist=file_exists)
          open(newunit=iunit, file="ellipsoid_bc_error.dat", status="unknown", position="append")
          if (.not. file_exists) then
             write(iunit,*) "# t samples normal_max normal_rms full_max full_rms xnu"
          endif
          write(iunit,*) t, global_count, global_normal_max, normal_rms, &
               global_full_max, full_rms, xnu
          close(iunit)
       endif
    endif

end subroutine ellipsoid_bc_diagnostic

!********************************************************************
subroutine ellipsoid_pressure_correction_diagnostic(ux1, uy1, uz1, px1, py1, pz1, stage)

    use complex_geometry, only : nobjx, nobjy, nobjz, xi, xf, yi, yf, zi, zf
    use param, only : zero, dx, dy, dz, xlx, yly, zlz, izap, xnu
    use variables, only : yp
    use var, only : ux2, uy2, uz2, ux3, uy3, uz3, t
    use decomp_2d_mpi, only : nrank
    use MPI

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: px1, py1, pz1
    character(len=*), intent(in) :: stage

    integer :: i, j, k, iobj, ix, jy, kz, code, iunit
    integer :: local_count, global_count
    logical :: file_exists
    real(mytype) :: xm, ym, zm
    real(mytype) :: point(3), nearestVelocity(3), pressureCorrection(3)
    real(mytype) :: local_required_max, local_required_sum2
    real(mytype) :: local_applied_max, local_applied_sum2
    real(mytype) :: local_error_max, local_error_sum2
    real(mytype) :: global_required_max, global_required_sum2
    real(mytype) :: global_applied_max, global_applied_sum2
    real(mytype) :: global_error_max, global_error_sum2
    real(mytype) :: required_rms, applied_rms, error_rms
    real(mytype), allocatable, dimension(:,:,:) :: px2, py2, pz2
    real(mytype), allocatable, dimension(:,:,:) :: px3, py3, pz3

    local_count = 0
    local_required_max = zero
    local_required_sum2 = zero
    local_applied_max = zero
    local_applied_sum2 = zero
    local_error_max = zero
    local_error_sum2 = zero

    allocate(px2(ysize(1), ysize(2), ysize(3)))
    allocate(py2(ysize(1), ysize(2), ysize(3)))
    allocate(pz2(ysize(1), ysize(2), ysize(3)))
    allocate(px3(zsize(1), zsize(2), zsize(3)))
    allocate(py3(zsize(1), zsize(2), zsize(3)))
    allocate(pz3(zsize(1), zsize(2), zsize(3)))

    call transpose_x_to_y(ux1, ux2)
    call transpose_x_to_y(uy1, uy2)
    call transpose_x_to_y(uz1, uz2)
    call transpose_x_to_y(px1, px2)
    call transpose_x_to_y(py1, py2)
    call transpose_x_to_y(pz1, pz2)
    call transpose_y_to_z(ux2, ux3)
    call transpose_y_to_z(uy2, uy3)
    call transpose_y_to_z(uz2, uz3)
    call transpose_y_to_z(px2, px3)
    call transpose_y_to_z(py2, py3)
    call transpose_y_to_z(pz2, pz3)

    do k = 1, xsize(3)
       zm = real(xstart(3)+k-1, mytype) * dz
       do j = 1, xsize(2)
          ym = real(xstart(2)+j-1, mytype) * dy
          do iobj = 1, nobjx(j,k)
             if (xi(iobj,j,k) .gt. zero) then
                ix = xi(iobj,j,k) / dx + 1
                if (izap .eq. 1) then
                   nearestVelocity = [ux1(ix-1,j,k), uy1(ix-1,j,k), uz1(ix-1,j,k)]
                   pressureCorrection = [px1(ix-1,j,k), py1(ix-1,j,k), pz1(ix-1,j,k)]
                else
                   nearestVelocity = [ux1(ix,j,k), uy1(ix,j,k), uz1(ix,j,k)]
                   pressureCorrection = [px1(ix,j,k), py1(ix,j,k), pz1(ix,j,k)]
                endif
                point = [xi(iobj,j,k), ym, zm]
                call accumulate_ellipsoid_pressure_sample(point, nearestVelocity, pressureCorrection, &
                     local_required_max, local_required_sum2, local_applied_max, local_applied_sum2, &
                     local_error_max, local_error_sum2, local_count)
             endif

             if (xf(iobj,j,k) .lt. xlx) then
                ix = (xf(iobj,j,k) + dx) / dx + 1
                if (izap .eq. 1) then
                   nearestVelocity = [ux1(ix+1,j,k), uy1(ix+1,j,k), uz1(ix+1,j,k)]
                   pressureCorrection = [px1(ix+1,j,k), py1(ix+1,j,k), pz1(ix+1,j,k)]
                else
                   nearestVelocity = [ux1(ix,j,k), uy1(ix,j,k), uz1(ix,j,k)]
                   pressureCorrection = [px1(ix,j,k), py1(ix,j,k), pz1(ix,j,k)]
                endif
                point = [xf(iobj,j,k), ym, zm]
                call accumulate_ellipsoid_pressure_sample(point, nearestVelocity, pressureCorrection, &
                     local_required_max, local_required_sum2, local_applied_max, local_applied_sum2, &
                     local_error_max, local_error_sum2, local_count)
             endif
          enddo
       enddo
    enddo

    do k = 1, ysize(3)
       zm = real(ystart(3)+k-1, mytype) * dz
       do i = 1, ysize(1)
          xm = real(ystart(1)+i-1, mytype) * dx
          do j = 1, nobjy(i,k)
             if (yi(j,i,k) .gt. zero) then
                jy = 1
                do while (yp(jy) .lt. yi(j,i,k))
                   jy = jy + 1
                enddo
                jy = jy - 1
                if (izap .eq. 1) then
                   nearestVelocity = [ux2(i,jy-1,k), uy2(i,jy-1,k), uz2(i,jy-1,k)]
                   pressureCorrection = [px2(i,jy-1,k), py2(i,jy-1,k), pz2(i,jy-1,k)]
                else
                   nearestVelocity = [ux2(i,jy,k), uy2(i,jy,k), uz2(i,jy,k)]
                   pressureCorrection = [px2(i,jy,k), py2(i,jy,k), pz2(i,jy,k)]
                endif
                point = [xm, yi(j,i,k), zm]
                call accumulate_ellipsoid_pressure_sample(point, nearestVelocity, pressureCorrection, &
                     local_required_max, local_required_sum2, local_applied_max, local_applied_sum2, &
                     local_error_max, local_error_sum2, local_count)
             endif

             if (yf(j,i,k) .lt. yly) then
                jy = 1
                do while (yp(jy) .lt. yf(j,i,k))
                   jy = jy + 1
                enddo
                if (izap .eq. 1) then
                   nearestVelocity = [ux2(i,jy+1,k), uy2(i,jy+1,k), uz2(i,jy+1,k)]
                   pressureCorrection = [px2(i,jy+1,k), py2(i,jy+1,k), pz2(i,jy+1,k)]
                else
                   nearestVelocity = [ux2(i,jy,k), uy2(i,jy,k), uz2(i,jy,k)]
                   pressureCorrection = [px2(i,jy,k), py2(i,jy,k), pz2(i,jy,k)]
                endif
                point = [xm, yf(j,i,k), zm]
                call accumulate_ellipsoid_pressure_sample(point, nearestVelocity, pressureCorrection, &
                     local_required_max, local_required_sum2, local_applied_max, local_applied_sum2, &
                     local_error_max, local_error_sum2, local_count)
             endif
          enddo
       enddo
    enddo

    do j = 1, zsize(2)
       ym = real(zstart(2)+j-1, mytype) * dy
       do i = 1, zsize(1)
          xm = real(zstart(1)+i-1, mytype) * dx
          do k = 1, nobjz(i,j)
             if (zi(k,i,j) .gt. zero) then
                kz = zi(k,i,j) / dz + 1
                if (izap .eq. 1) then
                   nearestVelocity = [ux3(i,j,kz-1), uy3(i,j,kz-1), uz3(i,j,kz-1)]
                   pressureCorrection = [px3(i,j,kz-1), py3(i,j,kz-1), pz3(i,j,kz-1)]
                else
                   nearestVelocity = [ux3(i,j,kz), uy3(i,j,kz), uz3(i,j,kz)]
                   pressureCorrection = [px3(i,j,kz), py3(i,j,kz), pz3(i,j,kz)]
                endif
                point = [xm, ym, zi(k,i,j)]
                call accumulate_ellipsoid_pressure_sample(point, nearestVelocity, pressureCorrection, &
                     local_required_max, local_required_sum2, local_applied_max, local_applied_sum2, &
                     local_error_max, local_error_sum2, local_count)
             endif

             if (zf(k,i,j) .lt. zlz) then
                kz = (zf(k,i,j) + dz) / dz + 1
                if (izap .eq. 1) then
                   nearestVelocity = [ux3(i,j,kz+1), uy3(i,j,kz+1), uz3(i,j,kz+1)]
                   pressureCorrection = [px3(i,j,kz+1), py3(i,j,kz+1), pz3(i,j,kz+1)]
                else
                   nearestVelocity = [ux3(i,j,kz), uy3(i,j,kz), uz3(i,j,kz)]
                   pressureCorrection = [px3(i,j,kz), py3(i,j,kz), pz3(i,j,kz)]
                endif
                point = [xm, ym, zf(k,i,j)]
                call accumulate_ellipsoid_pressure_sample(point, nearestVelocity, pressureCorrection, &
                     local_required_max, local_required_sum2, local_applied_max, local_applied_sum2, &
                     local_error_max, local_error_sum2, local_count)
             endif
          enddo
       enddo
    enddo

    deallocate(px2, py2, pz2, px3, py3, pz3)

    global_count = local_count
    global_required_max = local_required_max
    global_required_sum2 = local_required_sum2
    global_applied_max = local_applied_max
    global_applied_sum2 = local_applied_sum2
    global_error_max = local_error_max
    global_error_sum2 = local_error_sum2

    call MPI_Allreduce(MPI_IN_PLACE, global_count, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, code)
    call MPI_Allreduce(MPI_IN_PLACE, global_required_max, 1, real_type, MPI_MAX, MPI_COMM_WORLD, code)
    call MPI_Allreduce(MPI_IN_PLACE, global_required_sum2, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
    call MPI_Allreduce(MPI_IN_PLACE, global_applied_max, 1, real_type, MPI_MAX, MPI_COMM_WORLD, code)
    call MPI_Allreduce(MPI_IN_PLACE, global_applied_sum2, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
    call MPI_Allreduce(MPI_IN_PLACE, global_error_max, 1, real_type, MPI_MAX, MPI_COMM_WORLD, code)
    call MPI_Allreduce(MPI_IN_PLACE, global_error_sum2, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)

    required_rms = zero
    applied_rms = zero
    error_rms = zero
    if (global_count .gt. 0) then
       required_rms = sqrt(global_required_sum2 / real(global_count, mytype))
       applied_rms = sqrt(global_applied_sum2 / real(global_count, mytype))
       error_rms = sqrt(global_error_sum2 / real(global_count, mytype))
    endif

    if (nrank .eq. 0) then
       write(*,*) "Ellipsoid pressure correction diagnostic ", trim(stage), &
            ": samples=", global_count, &
            " required max/rms=", global_required_max, required_rms, &
            " applied max/rms=", global_applied_max, applied_rms, &
            " error max/rms=", global_error_max, error_rms

       inquire(file="ellipsoid_pressure_correction_error.dat", exist=file_exists)
       open(newunit=iunit, file="ellipsoid_pressure_correction_error.dat", status="unknown", position="append")
       if (.not. file_exists) then
          write(iunit,*) "# t itime itr stage samples required_max required_rms applied_max applied_rms error_max error_rms xnu"
       endif
       write(iunit,*) t, itime, itr, trim(stage), global_count, &
            global_required_max, required_rms, global_applied_max, applied_rms, &
            global_error_max, error_rms, xnu
       close(iunit)
    endif

end subroutine ellipsoid_pressure_correction_diagnostic

!********************************************************************
subroutine ellipsoid_projection_slip_correction(ux1, uy1, uz1, stage)

    use complex_geometry, only : nobjx, nobjy, nobjz, xi, xf, yi, yf, zi, zf
    use param, only : zero, dx, dy, dz, xlx, yly, zlz, izap, xnu
    use variables, only : yp
    use var, only : ux2, uy2, uz2, ux3, uy3, uz3, t
    use decomp_2d_mpi, only : nrank
    use MPI

    implicit none

    real(mytype), intent(inout), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    character(len=*), intent(in) :: stage

    integer :: i, j, k, iobj, ix, jy, kz, code, iunit
    integer :: local_count, global_count
    logical :: file_exists
    real(mytype) :: xm, ym, zm
    real(mytype) :: point(3), velocity(3)
    real(mytype) :: local_before_max, local_before_sum2
    real(mytype) :: local_after_max, local_after_sum2
    real(mytype) :: local_correction_max, local_correction_sum2
    real(mytype) :: global_before_max, global_before_sum2
    real(mytype) :: global_after_max, global_after_sum2
    real(mytype) :: global_correction_max, global_correction_sum2
    real(mytype) :: before_rms, after_rms, correction_rms

    local_count = 0
    local_before_max = zero
    local_before_sum2 = zero
    local_after_max = zero
    local_after_sum2 = zero
    local_correction_max = zero
    local_correction_sum2 = zero

    do k = 1, xsize(3)
       zm = real(xstart(3)+k-1, mytype) * dz
       do j = 1, xsize(2)
          ym = real(xstart(2)+j-1, mytype) * dy
          do iobj = 1, nobjx(j,k)
             if (xi(iobj,j,k) .gt. zero) then
                ix = xi(iobj,j,k) / dx + 1
                if (izap .eq. 1) ix = ix - 1
                if (ix.ge.1 .and. ix.le.xsize(1)) then
                   point = [xi(iobj,j,k), ym, zm]
                   velocity = [ux1(ix,j,k), uy1(ix,j,k), uz1(ix,j,k)]
                   call project_ellipsoid_slip_sample(point, velocity, local_before_max, &
                        local_before_sum2, local_after_max, local_after_sum2, &
                        local_correction_max, local_correction_sum2, local_count)
                   ux1(ix,j,k) = velocity(1)
                   uy1(ix,j,k) = velocity(2)
                   uz1(ix,j,k) = velocity(3)
                endif
             endif

             if (xf(iobj,j,k) .lt. xlx) then
                ix = (xf(iobj,j,k) + dx) / dx + 1
                if (izap .eq. 1) ix = ix + 1
                if (ix.ge.1 .and. ix.le.xsize(1)) then
                   point = [xf(iobj,j,k), ym, zm]
                   velocity = [ux1(ix,j,k), uy1(ix,j,k), uz1(ix,j,k)]
                   call project_ellipsoid_slip_sample(point, velocity, local_before_max, &
                        local_before_sum2, local_after_max, local_after_sum2, &
                        local_correction_max, local_correction_sum2, local_count)
                   ux1(ix,j,k) = velocity(1)
                   uy1(ix,j,k) = velocity(2)
                   uz1(ix,j,k) = velocity(3)
                endif
             endif
          enddo
       enddo
    enddo

    call transpose_x_to_y(ux1, ux2)
    call transpose_x_to_y(uy1, uy2)
    call transpose_x_to_y(uz1, uz2)

    do k = 1, ysize(3)
       zm = real(ystart(3)+k-1, mytype) * dz
       do i = 1, ysize(1)
          xm = real(ystart(1)+i-1, mytype) * dx
          do j = 1, nobjy(i,k)
             if (yi(j,i,k) .gt. zero) then
                jy = 1
                do while (yp(jy) .lt. yi(j,i,k))
                   jy = jy + 1
                enddo
                jy = jy - 1
                if (izap .eq. 1) jy = jy - 1
                if (jy.ge.1 .and. jy.le.ysize(2)) then
                   point = [xm, yi(j,i,k), zm]
                   velocity = [ux2(i,jy,k), uy2(i,jy,k), uz2(i,jy,k)]
                   call project_ellipsoid_slip_sample(point, velocity, local_before_max, &
                        local_before_sum2, local_after_max, local_after_sum2, &
                        local_correction_max, local_correction_sum2, local_count)
                   ux2(i,jy,k) = velocity(1)
                   uy2(i,jy,k) = velocity(2)
                   uz2(i,jy,k) = velocity(3)
                endif
             endif

             if (yf(j,i,k) .lt. yly) then
                jy = 1
                do while (yp(jy) .lt. yf(j,i,k))
                   jy = jy + 1
                enddo
                if (izap .eq. 1) jy = jy + 1
                if (jy.ge.1 .and. jy.le.ysize(2)) then
                   point = [xm, yf(j,i,k), zm]
                   velocity = [ux2(i,jy,k), uy2(i,jy,k), uz2(i,jy,k)]
                   call project_ellipsoid_slip_sample(point, velocity, local_before_max, &
                        local_before_sum2, local_after_max, local_after_sum2, &
                        local_correction_max, local_correction_sum2, local_count)
                   ux2(i,jy,k) = velocity(1)
                   uy2(i,jy,k) = velocity(2)
                   uz2(i,jy,k) = velocity(3)
                endif
             endif
          enddo
       enddo
    enddo

    call transpose_y_to_z(ux2, ux3)
    call transpose_y_to_z(uy2, uy3)
    call transpose_y_to_z(uz2, uz3)

    do j = 1, zsize(2)
       ym = real(zstart(2)+j-1, mytype) * dy
       do i = 1, zsize(1)
          xm = real(zstart(1)+i-1, mytype) * dx
          do k = 1, nobjz(i,j)
             if (zi(k,i,j) .gt. zero) then
                kz = zi(k,i,j) / dz + 1
                if (izap .eq. 1) kz = kz - 1
                if (kz.ge.1 .and. kz.le.zsize(3)) then
                   point = [xm, ym, zi(k,i,j)]
                   velocity = [ux3(i,j,kz), uy3(i,j,kz), uz3(i,j,kz)]
                   call project_ellipsoid_slip_sample(point, velocity, local_before_max, &
                        local_before_sum2, local_after_max, local_after_sum2, &
                        local_correction_max, local_correction_sum2, local_count)
                   ux3(i,j,kz) = velocity(1)
                   uy3(i,j,kz) = velocity(2)
                   uz3(i,j,kz) = velocity(3)
                endif
             endif

             if (zf(k,i,j) .lt. zlz) then
                kz = (zf(k,i,j) + dz) / dz + 1
                if (izap .eq. 1) kz = kz + 1
                if (kz.ge.1 .and. kz.le.zsize(3)) then
                   point = [xm, ym, zf(k,i,j)]
                   velocity = [ux3(i,j,kz), uy3(i,j,kz), uz3(i,j,kz)]
                   call project_ellipsoid_slip_sample(point, velocity, local_before_max, &
                        local_before_sum2, local_after_max, local_after_sum2, &
                        local_correction_max, local_correction_sum2, local_count)
                   ux3(i,j,kz) = velocity(1)
                   uy3(i,j,kz) = velocity(2)
                   uz3(i,j,kz) = velocity(3)
                endif
             endif
          enddo
       enddo
    enddo

    call transpose_z_to_y(ux3, ux2)
    call transpose_z_to_y(uy3, uy2)
    call transpose_z_to_y(uz3, uz2)
    call transpose_y_to_x(ux2, ux1)
    call transpose_y_to_x(uy2, uy1)
    call transpose_y_to_x(uz2, uz1)

    global_count = local_count
    global_before_max = local_before_max
    global_before_sum2 = local_before_sum2
    global_after_max = local_after_max
    global_after_sum2 = local_after_sum2
    global_correction_max = local_correction_max
    global_correction_sum2 = local_correction_sum2

    call MPI_Allreduce(MPI_IN_PLACE, global_count, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, code)
    call MPI_Allreduce(MPI_IN_PLACE, global_before_max, 1, real_type, MPI_MAX, MPI_COMM_WORLD, code)
    call MPI_Allreduce(MPI_IN_PLACE, global_before_sum2, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
    call MPI_Allreduce(MPI_IN_PLACE, global_after_max, 1, real_type, MPI_MAX, MPI_COMM_WORLD, code)
    call MPI_Allreduce(MPI_IN_PLACE, global_after_sum2, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
    call MPI_Allreduce(MPI_IN_PLACE, global_correction_max, 1, real_type, MPI_MAX, MPI_COMM_WORLD, code)
    call MPI_Allreduce(MPI_IN_PLACE, global_correction_sum2, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)

    before_rms = zero
    after_rms = zero
    correction_rms = zero
    if (global_count .gt. 0) then
       before_rms = sqrt(global_before_sum2 / real(global_count, mytype))
       after_rms = sqrt(global_after_sum2 / real(global_count, mytype))
       correction_rms = sqrt(global_correction_sum2 / real(global_count, mytype))
    endif

    if (nrank .eq. 0) then
       write(*,*) "Ellipsoid projection slip correction ", trim(stage), &
            ": samples=", global_count, &
            " before max/rms=", global_before_max, before_rms, &
            " after max/rms=", global_after_max, after_rms, &
            " correction max/rms=", global_correction_max, correction_rms

       inquire(file="ellipsoid_projection_slip_correction.dat", exist=file_exists)
       open(newunit=iunit, file="ellipsoid_projection_slip_correction.dat", status="unknown", position="append")
       if (.not. file_exists) then
          write(iunit,*) "# t itime itr stage samples before_max before_rms after_max after_rms correction_max correction_rms xnu"
       endif
       write(iunit,*) t, itime, itr, trim(stage), global_count, &
            global_before_max, before_rms, global_after_max, after_rms, &
            global_correction_max, correction_rms, xnu
       close(iunit)
    endif

end subroutine ellipsoid_projection_slip_correction

!********************************************************************
subroutine accumulate_ellipsoid_bc_sample(point, nearestVelocity, normal_max, &
     normal_sum2, full_max, full_sum2, count)

    use ellipsoid_utils, only : CalculatePointVelocity_Multi, EllipsoidNormal_Multi

    implicit none

    real(mytype), intent(in) :: point(3), nearestVelocity(3)
    real(mytype), intent(inout) :: normal_max, normal_sum2, full_max, full_sum2
    integer, intent(inout) :: count

    real(mytype) :: bodyVelocity(3), normal(3), delta(3)
    real(mytype) :: normal_err, full_err

    call CalculatePointVelocity_Multi(point, bodyVelocity)
    call EllipsoidNormal_Multi(point, normal)

    delta = nearestVelocity - bodyVelocity
    normal_err = abs(sum(delta * normal))
    full_err = sqrt(sum(delta * delta))

    normal_max = max(normal_max, normal_err)
    full_max = max(full_max, full_err)
    normal_sum2 = normal_sum2 + normal_err * normal_err
    full_sum2 = full_sum2 + full_err * full_err
    count = count + 1

end subroutine accumulate_ellipsoid_bc_sample

!********************************************************************
subroutine accumulate_ellipsoid_pressure_sample(point, nearestVelocity, pressureCorrection, &
     required_max, required_sum2, applied_max, applied_sum2, error_max, error_sum2, count)

    use ellipsoid_utils, only : CalculatePointVelocity_Multi, EllipsoidNormal_Multi

    implicit none

    real(mytype), intent(in) :: point(3), nearestVelocity(3), pressureCorrection(3)
    real(mytype), intent(inout) :: required_max, required_sum2
    real(mytype), intent(inout) :: applied_max, applied_sum2
    real(mytype), intent(inout) :: error_max, error_sum2
    integer, intent(inout) :: count

    real(mytype) :: bodyVelocity(3), normal(3)
    real(mytype) :: required, applied, error

    call CalculatePointVelocity_Multi(point, bodyVelocity)
    call EllipsoidNormal_Multi(point, normal)

    required = sum((nearestVelocity - bodyVelocity) * normal)
    applied = sum(pressureCorrection * normal)
    error = applied - required

    required_max = max(required_max, abs(required))
    applied_max = max(applied_max, abs(applied))
    error_max = max(error_max, abs(error))
    required_sum2 = required_sum2 + required * required
    applied_sum2 = applied_sum2 + applied * applied
    error_sum2 = error_sum2 + error * error
    count = count + 1

end subroutine accumulate_ellipsoid_pressure_sample

!********************************************************************
subroutine project_ellipsoid_slip_sample(point, velocity, before_max, before_sum2, &
     after_max, after_sum2, correction_max, correction_sum2, count)

    use param, only : zero
    use ellipsoid_utils, only : CalculatePointVelocity_Multi, EllipsoidNormal_Multi

    implicit none

    real(mytype), intent(in) :: point(3)
    real(mytype), intent(inout) :: velocity(3)
    real(mytype), intent(inout) :: before_max, before_sum2
    real(mytype), intent(inout) :: after_max, after_sum2
    real(mytype), intent(inout) :: correction_max, correction_sum2
    integer, intent(inout) :: count

    real(mytype) :: bodyVelocity(3), normal(3), correction(3)
    real(mytype) :: before_err, after_err, correction_mag

    call CalculatePointVelocity_Multi(point, bodyVelocity)
    call EllipsoidNormal_Multi(point, normal)

    before_err = sum((velocity - bodyVelocity) * normal)
    correction = -before_err * normal
    velocity = velocity + correction
    after_err = sum((velocity - bodyVelocity) * normal)
    correction_mag = sqrt(sum(correction * correction))

    before_max = max(before_max, abs(before_err))
    after_max = max(after_max, abs(after_err))
    correction_max = max(correction_max, correction_mag)
    before_sum2 = before_sum2 + before_err * before_err
    after_sum2 = after_sum2 + after_err * after_err
    correction_sum2 = correction_sum2 + correction_mag * correction_mag
    count = count + 1

end subroutine project_ellipsoid_slip_sample

!********************************************************************
! check_body_proximity
!
! Aborts the simulation if any pair of bodies comes within a separation
! at which their force control volumes would overlap. At that point the
! force integration starts double-counting fluid cells and any further
! physics is invalid, so it's cleaner to stop with a clear message than
! to continue silently into garbage.
!
! Threshold: cvl_scalar * (max(shape(i,:))*ra(i) + max(shape(j,:))*ra(j))
! This is the same scaling used for the body-vs-domain-boundary check
! in xcompact3d.f90 — it's the centre-to-centre distance at which the
! two force CVs first touch.
!********************************************************************
subroutine check_body_proximity()

    use param, only : itype, itype_cyl
    use ibm_param, only : nbody, position, shape, ra, cvl_scalar
    use decomp_2d_mpi, only : nrank
    use MPI

    implicit none

    integer :: i, j, code, ierror
    real(mytype) :: maxrad_i, maxrad_j, dist, min_sep

    if (nbody <= 1) return

    do i = 1, nbody - 1
       if (itype .eq. itype_cyl) then
          maxrad_i = max(shape(i,1), shape(i,2))
       else
          maxrad_i = max(shape(i,1), shape(i,2), shape(i,3))
       endif
       do j = i + 1, nbody
          if (itype .eq. itype_cyl) then
             maxrad_j = max(shape(j,1), shape(j,2))
          else
             maxrad_j = max(shape(j,1), shape(j,2), shape(j,3))
          endif
          dist = sqrt(sum((position(i,:) - position(j,:))**2))
          min_sep = cvl_scalar * (maxrad_i*ra(i) + maxrad_j*ra(j))
          if (dist < min_sep) then
             if (nrank == 0) then
                write(*,*) "Bodies", i, "and", j, "are too close!"
                write(*,*) "  centre-to-centre distance =", dist
                write(*,*) "  minimum allowed separation =", min_sep
             endif
             call MPI_ABORT(MPI_COMM_WORLD, 1, ierror)
          endif
       enddo
    enddo

end subroutine check_body_proximity

end module ellip
