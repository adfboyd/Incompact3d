!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module ellipsoid_utils

    use decomp_2d, only: mytype
    use param, only: zero, one, two
    use dbg_schemes, only: sqrt_prec, cos_prec, exp_prec, sin_prec
    
    implicit none
    ! public QuatRot, cross, IsoKernel, AnIsoKernel, int2str

contains

    ! !*******************************************************************************
    ! !
    ! real(mytype) function trilinear_interpolation(x0,y0,z0, &
    !                                               x1,y1,z1, &
    !                                               x,y,z, &
    !                                               u000,u100,u001,u101, &
    !                                               u010,u110,u011,u111)
    ! !
    ! !*******************************************************************************

    !   implicit none
    !   real(mytype),intent(in) :: x0,y0,z0,x1,y1,z1,x,y,z,u000,u100,u001,u101,u010,u110,u011,u111
    !   real(mytype) :: c00,c01,c10,c11,c0,c1,xd,yd,zd

    !   if (x1/=x0) then
    !      xd=(x-x0)/(x1-x0)
    !   else
    !      xd=zero
    !   endif

    !   if (y1/=y0) then
    !      yd=(y-y0)/(y1-y0)
    !   else
    !      yd=zero
    !   endif

    !   if (z1/=z0) then
    !      zd=(z-z0)/(z1-z0)
    !   else
    !      zd=zero
    !   endif

    !   ! Interpolate along X
    !   c00=u000*(one-xd)+u100*xd
    !   c01=u001*(one-xd)+u101*xd
    !   c10=u010*(one-xd)+u110*xd
    !   c11=u011*(one-xd)+u111*xd

    !   ! Interpolate along Y
    !   c0 = c00*(one-yd)+c10*yd
    !   c1 = c01*(one-yd)+c11*yd

    !   ! Interpolate along Z
    !   trilinear_interpolation=c0*(one-zd)+c1*zd

    !   return

    ! end function trilinear_interpolation

    ! !*******************************************************************************
    ! !
    ! subroutine cross(ax,ay,az,bx,by,bz,cx,cy,cz)
    ! !
    ! !*******************************************************************************

    !   real(mytype) :: ax,ay,az,bx,by,bz,cx,cy,cz

    !   cx = ay*bz - az*by
    !   cy = az*bx - ax*bz
    !   cz = ax*by - ay*bx

    ! end subroutine cross

    ! !*******************************************************************************
    ! !
    ! subroutine QuatRot(vx,vy,vz,Theta,Rx,Ry,Rz,Ox,Oy,Oz,vRx,vRy,vRz)
    ! !
    ! !*******************************************************************************

    !   implicit none
    !   real(mytype), intent(in) :: vx,vy,vz,Theta,Rx,Ry,Rz,Ox,Oy,Oz
    !   real(mytype),intent(inout) :: vRx,vRy,vRz
    !   real(mytype) :: nRx,nRy,nRz
    !   real(mytype) :: p(4,1), pR(4,1), q(4), qbar(4), RMag, vOx, vOy, vOz
    !   real(mytype) :: QL(4,4), QbarR(4,4)

    !   ! Perform rotation of vector v around normal vector nR using the quaternion machinery.
    !   ! v: input vector
    !   ! Theta: rotation angle (rad)
    !   ! nR: normal vector around which to rotate
    !   ! Origin: origin point of rotation
    !   ! vR: Rotated vector

    !   ! Force normalize nR
    !   RMag=sqrt_prec(Rx**2.0+Ry**2.0+Rz**2.0)
    !   nRx=Rx/RMag
    !   nRy=Ry/RMag
    !   nRz=Rz/RMag

    !   ! Quaternion form of v
    !   vOx=vx-Ox
    !   vOy=vy-Oy
    !   vOz=vz-Oz
    !   p=reshape([zero,vOx,vOy,vOz],[4,1])

    !   ! Rotation quaternion and conjugate
    !   q=(/cos_prec(Theta/2),nRx*sin_prec(Theta/2),nRy*sin_prec(Theta/2),nRz*sin_prec(Theta/2)/)
    !   qbar=(/q(1),-q(2),-q(3),-q(4)/)

    !   QL=transpose(reshape((/q(1), -q(2), -q(3), -q(4), &
    !                          q(2),  q(1), -q(4),  q(3), &
    !                          q(3),  q(4),  q(1), -q(2), &
    !                          q(4), -q(3),  q(2),  q(1)/),(/4,4/)))

    !   QbarR=transpose(reshape((/qbar(1), -qbar(2), -qbar(3), -qbar(4), &
    !                             qbar(2),  qbar(1),  qbar(4), -qbar(3), &
    !                             qbar(3), -qbar(4),  qbar(1),  qbar(2), &
    !                             qbar(4),  qbar(3), -qbar(2),  qbar(1)/),(/4,4/)))

    !   ! Rotate p
    !   pR=matmul(matmul(QbarR,QL),p)
    !   vRx=pR(2,1)+Ox
    !   vRy=pR(3,1)+Oy
    !   vRz=pR(4,1)+Oz

    ! end subroutine QuatRot



    subroutine NormalizeQuaternion(quaternion) 
      real(mytype), intent(in) :: quaternion(4)
      real(mytype) :: normalizedQuaternion(4)
    
      ! Compute the magnitude of the quaternion
      real(mytype) :: magnitude
      magnitude = sqrt(quaternion(1)**2 + quaternion(2)**2 + quaternion(3)**2 + quaternion(4)**2)
    
      ! Normalize the quaternion
      normalizedQuaternion = quaternion / magnitude
    
    end subroutine NormalizeQuaternion

    subroutine QuaternionConjugate(q, q_c)
      real(mytype), intent(in) :: q(4)
      real(mytype), intent(out) :: q_c(4)

      q_c = [q(1), -q(2), -q(3), -q(4)]
    end subroutine


    subroutine QuaternionMultiply(q1, q2, result)
      real(mytype), intent(in) :: q1(4), q2(4)
      real(mytype), intent(out) :: result(4)
      
      result(1) = q1(1) * q2(1) - q1(2) * q2(2) - q1(3) * q2(3) - q1(4) * q2(4)
      result(2) = q1(1) * q2(2) + q1(2) * q2(1) + q1(3) * q2(4) - q1(4) * q2(3)
      result(3) = q1(1) * q2(3) - q1(2) * q2(4) + q1(3) * q2(1) + q1(4) * q2(2)
      result(4) = q1(1) * q2(4) + q1(2) * q2(3) - q1(3) * q2(2) + q1(4) * q2(1)
    end subroutine QuaternionMultiply


    subroutine RotatePoint(point, quaternion, rotatedPoint)
      real(mytype), intent(in) :: point(3), quaternion(4)
      real(mytype), intent(out) :: rotatedPoint(3)
      real(mytype) :: conjugateQuaternion(4)
      real(mytype) :: resultQuaternion(4)
      real(mytype) :: rotatedPointQuaternion(4)

      ! Convert the point to a quaternion
      real(mytype) :: pointQuaternion(4)
      pointQuaternion(1) = 0.0D0
      pointQuaternion(2:4) = point(:)
    
      ! Perform the rotation
      
      conjugateQuaternion = [quaternion(1), -quaternion(2), -quaternion(3), -quaternion(4)]
    
      call QuaternionMultiply(quaternion, pointQuaternion, resultQuaternion)
      call QuaternionMultiply(resultQuaternion, conjugateQuaternion, rotatedPointQuaternion)
    
      ! Convert the rotated quaternion back to a 3D point
      rotatedPoint = rotatedPointQuaternion(2:4)
    end subroutine RotatePoint

    subroutine EllipsoidalRadius(point, centre, orientation, shape, radius)
      real(mytype), intent(in) :: point(3), centre(3), orientation(4), shape(3)
      real(mytype), intent(out) :: radius
      real(mytype) :: trans_point(3),rotated_point(3),scaled_point(3), orientation_c(4)
      integer :: i

      !translate point to body frame
      trans_point = point-centre

      call QuaternionConjugate(orientation, orientation_c)

      !rotate point into body frame (using inverse(conjugate) of orientation)
      call RotatePoint(trans_point, orientation, rotated_point)

      do i = 1,3
         scaled_point(i)=rotated_point(i)/shape(i)
      end do 

      radius=sqrt_prec(scaled_point(1)**two+scaled_point(2)**two+scaled_point(3)**two)

   end subroutine

    
    
    

    subroutine CrossProduct(a, b, result)
      real(mytype), intent(in) :: a(3), b(3)
      real(mytype), intent(inout) :: result(3)
    
      result(1) = a(2) * b(3) - a(3) * b(2)
      result(2) = a(3) * b(1) - a(1) * b(3)
      result(3) = a(1) * b(2) - a(2) * b(1)
    end subroutine CrossProduct
    
    subroutine CalculatePointVelocity(point, center, linearVelocity, angularVelocity, pointVelocity)
      real(mytype), intent(in) :: point(3), center(3), linearVelocity(3), angularVelocity(3)
      real(mytype), intent(out) :: pointVelocity(3)
      real(mytype) :: crossed(3)
      ! Compute the distance vector from the center to the point
      real(mytype) :: distance(3)
      distance = point - center
    
      ! Compute the cross product of angular velocity and distance vector
      
      call CrossProduct(angularVelocity, distance, crossed)
    
      ! Calculate the velocity at the point
      pointVelocity = crossed + linearVelocity
    end subroutine CalculatePointVelocity

    subroutine navierFieldGen(center, linearVelocity, angularVelocity, ep1, ep1_x, ep1_y, ep1_z)
      use param
      use decomp_2d
      real(mytype), intent(in) :: center(3), linearVelocity(3), angularVelocity(3)
      real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ep1
      real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(out) :: ep1_x, ep1_y, ep1_z
      real(mytype) :: xm, ym, zm, point(3), x_pv, y_pv, z_pv, pointVelocity(3)
      integer :: i,j,k

      do k = 1,xsize(3)
        zm=real(k+xstart(3)-2, mytype)*dz
        do j = 1,xsize(2)
          ym=real(j+xstart(2)-2, mytype)*dy
          do i = 1,xsize(1)
            xm=real(i+xstart(1)-2, mytype)*dx
            point=[xm,ym,zm]
            if (ep1(i,j,k).eq.1) then 
              call CalculatePointVelocity(point, center, linearVelocity, angularVelocity, pointVelocity)
              x_pv=pointVelocity(1)
              y_pv=pointVelocity(2)
              z_pv=pointVelocity(3)
            else
              x_pv=0
              y_pv=0
              z_pv=0
            endif
            ep1_x(i,j,k)=x_pv
            ep1_y(i,j,k)=y_pv
            ep1_z(i,j,k)=z_pv
          end do
        end do
      end do

      end subroutine navierFieldGen

    

   !  subroutine CalculatePointVelocity(point, center, angularVelocity, velocity)
   !    real(mytype), intent(in) :: point(3), center(3), angularVelocity(3)
   !    real(mytype), intent(out) :: velocity(3)
   !    real(mytype) :: distance(3)
   !    real(mytype) :: crossProduct(3)

    
   !    ! Compute the distance vector from the center to the point
   !    distance = point - center
    
   !    ! Compute the cross product of angular velocity and distance vector
   !    call CrossProduct(angularVelocity, distance, crossProduct)
    
   !    ! Calculate the velocity at the point
   !    velocity = crossProduct
   !  end subroutine CalculatePointVelocity
    

    ! !*******************************************************************************
    ! !
    ! subroutine IDW(Ncol,Xcol,Ycol,Zcol,Fxcol,Fycol,Fzcol,p,Xmesh,Ymesh,Zmesh,Fxmesh,Fymesh,Fzmesh)
    ! !
    ! !*******************************************************************************

    !   implicit none
    !   integer, intent(in) :: Ncol
    !   real(mytype), dimension(Ncol), intent(in) :: Xcol,Ycol,Zcol,Fxcol,Fycol,Fzcol
    !   real(mytype), intent(in) :: Xmesh,Ymesh,Zmesh
    !   integer,intent(in) :: p
    !   real(mytype), intent(inout) :: Fxmesh,Fymesh,Fzmesh
    !   real(mytype), dimension(Ncol) :: d(Ncol), w(Ncol)
    !   real(mytype) ::  wsum
    !   integer :: i,imin

    !   wsum=zero
    !   do i=1,Ncol
    !      d(i)=sqrt_prec((Xcol(i)-Xmesh)**2+(Ycol(i)-Ymesh)**2+(Zcol(i)-Zmesh)**2)
    !      w(i)=one/d(i)**p
    !      wsum=wsum+w(i)
    !   enddo

    !   if (minval(d)<0.001_mytype) then
    !      imin=minloc(d,1)
    !      Fxmesh=Fxcol(imin)
    !      Fymesh=Fycol(imin)
    !      Fzmesh=Fzcol(imin)
    !   else
    !      Fxmesh=zero
    !      Fymesh=zero
    !      Fzmesh=zero
    !      do i=1,Ncol
    !         Fxmesh=Fxmesh+w(i)*Fxcol(i)/wsum
    !         Fymesh=Fymesh+w(i)*Fycol(i)/wsum
    !         Fzmesh=Fzmesh+w(i)*Fzcol(i)/wsum
    !      enddo
    !   endif

    ! end subroutine IDW

    ! !*******************************************************************************
    ! !
    ! real(mytype) function IsoKernel(dr,epsilon_par,dim)
    ! !
    ! !*******************************************************************************

    !   use constants

    !   implicit none
    !   integer, intent(in) :: dim
    !   real(mytype), intent(in) :: dr, epsilon_par

    !   if (dim==2) then
    !      IsoKernel = one/(epsilon_par**2*pi)*exp_prec(-(dr/epsilon_par)**2.0)
    !   else if (dim==3) then
    !      IsoKernel = one/(epsilon_par**3.0*pi**1.5)*exp_prec(-(dr/epsilon_par)**2.0)
    !   else
    !      write(*,*) "1D source not implemented"
    !      stop
    !   endif

    ! end function IsoKernel

    ! !*******************************************************************************
    ! !
    ! real(mytype) function AnIsoKernel(dx,dy,dz,nx,ny,nz,tx,ty,tz,sx,sy,sz,ec,et,es)
    ! !
    ! !*******************************************************************************

    !   use constants

    !   implicit none
    !   real(mytype), intent(in) :: dx,dy,dz,nx,ny,nz,tx,ty,tz,sx,sy,sz,ec,et,es
    !   real(mytype) :: n,t,s

    !   n=dx*nx+dy*ny+dz*nz ! normal projection
    !   t=dx*tx+dy*ty+dz*tz ! Chordwise projection
    !   s=dx*sx+dy*sy+dz*sz ! Spanwise projection

    !   if (abs(s)<=es) then
    !      AnIsoKernel = exp_prec(-((n/et)**2.0+(t/ec)**2.0))/(ec*et*pi)
    !   else
    !      AnIsoKernel = zero
    !   endif

    ! end function AnIsoKernel

    ! !*******************************************************************************
    ! !
    ! integer function FindMinimum(x,Start,End)
    ! !
    ! !*******************************************************************************

    !   implicit none
    !   integer, dimension(1:), intent(in) :: x
    !   integer, intent(in) :: Start, End
    !   integer :: Minimum
    !   integer :: Location
    !   integer :: i

    !   minimum = x(start)
    !   Location = Start
    !   do i=start+1,End
    !      if (x(i) < Minimum) then
    !         Minimum = x(i)
    !         Location = i
    !      endif
    !   enddo
    !   FindMinimum = Location

    ! end function FindMinimum

    ! !*******************************************************************************
    ! !
    ! subroutine swap(a,b)
    ! !
    ! !*******************************************************************************

    !   implicit none
    !   integer, intent(inout) :: a,b
    !   integer :: Temp

    !   Temp = a
    !   a = b
    !   b = Temp

    ! end subroutine swap

    ! !*******************************************************************************
    ! !
    ! subroutine sort(x,size)
    ! !
    ! !*******************************************************************************

    !   implicit none
    !   integer, dimension(1:), intent(inout) :: x
    !   integer, intent(in) :: size
    !   integer :: i
    !   integer :: Location

    !   do i=1,Size-1
    !      location=FindMinimum(x,i,size)
    !      call swap(x(i),x(Location))
    !   enddo

    ! end subroutine sort

    ! !*******************************************************************************
    ! !
    ! function dirname(number)
    ! !
    ! !*******************************************************************************

    !   integer, intent(in) :: number
    !   character(len=6) :: dirname

    !   ! Cast the (rounded) number to string using 6 digits and leading zeros
    !   write (dirname, '(I6.1)')  number
    !   ! This is the same w/o leading zeros
    !   !write (dirname, '(I6)')  nint(number)
    !   ! This is for one digit (no rounding)
    !   !write (dirname, '(F4.1)')  number

    ! end function dirname

    ! !*******************************************************************************
    ! !
    ! function outdirname(number)
    ! !
    ! !*******************************************************************************

    !   integer, intent(in) :: number
    !   character(len=6) :: outdirname

    !   ! Cast the (rounded) number to string using 6 digits and leading zeros
    !   write (outdirname, '(I6.1)')  number
    !   ! This is the same w/o leading zeros
    !   !write (dirname, '(I6)')  nint(number)
    !   ! This is for one digit (no rounding)
    !   !write (dirname, '(F4.1)')  number
 
    ! end function outdirname

    ! !*******************************************************************************
    ! !
    ! character(20) function int2str(num)
    ! !
    ! !*******************************************************************************

    !   integer, intent(in) ::num
    !   character(20) :: str

    !   ! convert integer to string using formatted write
    !   write(str, '(i20)') num
    !   int2str = adjustl(str)

    ! end function int2str

end module ellipsoid_utils
