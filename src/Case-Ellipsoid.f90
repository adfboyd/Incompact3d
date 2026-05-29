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
            check_body_proximity

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

#ifdef DEBG
    if (nrank .eq. 0) write(*,*) '# init end ok'
#endif

    call init_body_dat()

    return
end subroutine init_ellip

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

subroutine update_ellipsoid(ux1, uy1, uz1, ep1)

    use forces, only : force, torque_calc, nvol, iforces
    use ellipsoid_utils, only : lin_step, ang_step
    use ibm_param
    use param, only : zero, dt
    use variables, only : ilist
    use var, only : itime, t
    use decomp_2d_mpi, only : nrank

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1

    real(mytype) :: drag(10), lift(10), lat(10)
    real(mytype) :: grav_effx(10), grav_effy(10), grav_effz(10)
    real(mytype) :: xtorq(10), ytorq(10), ztorq(10)
    integer :: i

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

    do i = 1, nvol
       call lin_step(position(i,:), linearVelocity(i,:), linearForce(i,:), ellip_m(i), ellip_m_added(i,:), orientation(i,:), dt, position_1, linearVelocity_1)
       call ang_step(orientation(i,:), angularVelocity(i,:), torque(i,:), inertia(i,:,:), inertia_rot_added(i,:), dt, orientation_1, angularVelocity_1)
       position(i,:) = position_1
       linearVelocity(i,:) = linearVelocity_1
       orientation(i,:) = orientation_1
       angularVelocity(i,:) = angularVelocity_1
    enddo

    if (nrank==0 .and. mod(itime,ilist)==0) then
       do i = 1, nbody
          write(11+i,*) t, position(i,1), position(i,2), position(i,3), &
               orientation(i,1), orientation(i,2), orientation(i,3), orientation(i,4), &
               linearVelocity(i,1), linearVelocity(i,2), linearVelocity(i,3), &
               angularVelocity(i,2), angularVelocity(i,3), angularVelocity(i,4), &
               linearForce(i,1), linearForce(i,2), linearForce(i,3), &
               torque(i,1), torque(i,2), torque(i,3)
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

end subroutine update_ellipsoid

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
  