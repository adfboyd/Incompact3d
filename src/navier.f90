!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module navier

  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d

  implicit none

  private
  logical, save :: suppress_divergence_output = .false.

  public :: solve_poisson, divergence,divergence2, calc_divu_constraint
  public :: pre_correc, cor_vel
  public :: ellipsoid_schur_projection
  public :: lmn_t_to_rho_trans, momentum_to_velocity, velocity_to_momentum
  public :: gradp, tbl_flrt

contains
  !############################################################################
  !!  SUBROUTINE: solve_poisson
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Takes the intermediate momentum field as input,
  !!              computes div and solves pressure-Poisson equation.
  !############################################################################
  SUBROUTINE solve_poisson(div_visu_var, pp3, px1, py1, pz1, rho1, ux1, uy1, uz1, ep1, drho1, divu3, &
       save_bc_grad, quiet)

    USE decomp_2d_poisson, ONLY : poisson
    USE var, ONLY : nzmsize
    USE var, ONLY : dv3
    USE param, ONLY : ntime, nrhotime, npress
    USE param, ONLY : ilmn, ivarcoeff, zero, one 
    USE mpi
    
    implicit none

    !! Inputs
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ux1, uy1, uz1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ep1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime), INTENT(IN) :: rho1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), ntime), INTENT(IN) :: drho1
    REAL(mytype), DIMENSION(zsize(1), zsize(2), zsize(3)), INTENT(IN) :: divu3
    LOGICAL, INTENT(IN), OPTIONAL :: save_bc_grad
    LOGICAL, INTENT(IN), OPTIONAL :: quiet

    !! Outputs
    REAL(mytype), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3, div_visu_var
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: px1, py1, pz1

    !! Locals
    INTEGER :: nlock, poissiter
    LOGICAL :: converged
    LOGICAL :: save_bc_grad_local, quiet_local, old_suppress_divergence_output
    REAL(mytype) :: atol, rtol, rho0, divup3norm
#ifdef DEBG
    real(mytype) :: dep
    integer :: code
#endif

    nlock = 1 !! Corresponds to computing div(u*)
    converged = .FALSE.
    poissiter = 0
    save_bc_grad_local = .TRUE.
    IF (PRESENT(save_bc_grad)) save_bc_grad_local = save_bc_grad
    quiet_local = .FALSE.
    IF (PRESENT(quiet)) quiet_local = quiet
    call calc_rho0(rho1(:,:,:,1), rho0)
#ifdef DOUBLE_PREC
    atol = 1.0e-14_mytype !! Absolute tolerance for Poisson solver
    rtol = 1.0e-14_mytype !! Relative tolerance for Poisson solver
#else
    atol = 1.0e-9_mytype !! Absolute tolerance for Poisson solver
    rtol = 1.0e-9_mytype !! Relative tolerance for Poisson solver
#endif

    IF (ilmn.AND.ivarcoeff) THEN
       !! Variable-coefficient Poisson solver works on div(u), not div(rho u)
       !! rho u -> u
       CALL momentum_to_velocity(rho1, ux1, uy1, uz1)
    ENDIF

    old_suppress_divergence_output = suppress_divergence_output
    if (quiet_local) suppress_divergence_output = .true.
    call divergence(div_visu_var(:,:,:,1),pp3(:,:,:,1),rho1,ux1,uy1,uz1,ep1,drho1,divu3,nlock)
    suppress_divergence_output = old_suppress_divergence_output
    IF (ilmn.AND.ivarcoeff) THEN
       dv3(:,:,:) = pp3(:,:,:,1)
    ENDIF

    do while(.not.converged)
       if (ivarcoeff) then
          !! Test convergence
          CALL test_varcoeff(converged, divup3norm, pp3, dv3, atol, rtol, poissiter)

          IF (.NOT.converged) THEN
             !! Evaluate additional RHS terms
             CALL calc_varcoeff_rhs(pp3(:,:,:,1), rho1, px1, py1, pz1, dv3, drho1, ep1, divu3, rho0)
          ENDIF

       ENDIF

       IF (.NOT.converged) THEN
#ifdef DEBG
          dep=maxval(abs(pp3(:,:,:,1)))
          call MPI_ALLREDUCE(MPI_IN_PLACE,dep,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
          if (nrank == 0) write(*,*)'## Solve Poisson before1 pp3', dep
#endif
          CALL poisson(pp3(:,:,:,1))
#ifdef DEBG
          dep=maxval(abs(pp3(:,:,:,1)))
          call MPI_ALLREDUCE(MPI_IN_PLACE,dep,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
          if (nrank == 0) write(*,*)'## Solve Poisson after call  pp3', dep
#endif

          !! Need to update pressure gradient here for varcoeff
          CALL gradp(px1,py1,pz1,pp3(:,:,:,1),save_bc_grad_local)
#ifdef DEBG
          dep=maxval(abs(pp3(:,:,:,1)))
          call MPI_ALLREDUCE(MPI_IN_PLACE,dep,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
          if (nrank == 0) write(*,*)'## Solve Poisson pp3', dep
          dep=maxval(abs(px1))
          call MPI_ALLREDUCE(MPI_IN_PLACE,dep,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
          if (nrank == 0) write(*,*)'## Solve Poisson px', dep
          dep=maxval(abs(py1))
          call MPI_ALLREDUCE(MPI_IN_PLACE,dep,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
          if (nrank == 0) write(*,*)'## Solve Poisson py', dep
          dep=maxval(abs(pz1))
          call MPI_ALLREDUCE(MPI_IN_PLACE,dep,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
          if (nrank == 0) write(*,*)'## Solve Poisson pz', dep
#endif
         

          IF ((.NOT.ilmn).OR.(.NOT.ivarcoeff)) THEN
             !! Once-through solver
             !! - Incompressible flow
             !! - LMN - constant-coefficient solver
             converged = .TRUE.
          ENDIF
       ENDIF

       poissiter = poissiter + 1
    enddo

    IF (ilmn.AND.ivarcoeff) THEN
       !! Variable-coefficient Poisson solver works on div(u), not div(rho u)
       !! u -> rho u
       CALL velocity_to_momentum(rho1, ux1, uy1, uz1)
    ENDIF

  END SUBROUTINE solve_poisson
  !############################################################################
  !!
  !!  SUBROUTINE: lmn_t_to_rho_trans
  !! DESCRIPTION: Converts the temperature transient to the density transient
  !!              term. This is achieved by application of EOS and chain rule.
  !!      INPUTS: dtemp1 - the RHS of the temperature equation.
  !!                rho1 - the density field.
  !!     OUTPUTS:  drho1 - the RHS of the density equation.
  !!      AUTHOR: Paul Bartholomew
  !!
  !############################################################################
  SUBROUTINE lmn_t_to_rho_trans(drho1, dtemp1, rho1, dphi1, phi1)

    USE param, ONLY : zero
    USE param, ONLY : imultispecies, massfrac, mol_weight
    USE param, ONLY : ntime
    USE var, ONLY : numscalar
    USE var, ONLY : ta1, tb1

    IMPLICIT NONE

    !! INPUTS
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3)) :: dtemp1, rho1
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), numscalar) :: phi1
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), ntime, numscalar) :: dphi1

    !! OUTPUTS
    REAL(mytype), INTENT(OUT), DIMENSION(xsize(1), xsize(2), xsize(3)) :: drho1

    !! LOCALS
    INTEGER :: is

    drho1(:,:,:) = zero

    IF (imultispecies) THEN
       DO is = 1, numscalar
          IF (massfrac(is)) THEN
             drho1(:,:,:) = drho1(:,:,:) - dphi1(:,:,:,1,is) / mol_weight(is)
          ENDIF
       ENDDO

       ta1(:,:,:) = zero !! Mean molecular weight
       DO is = 1, numscalar
          IF (massfrac(is)) THEN
             ta1(:,:,:) = ta1(:,:,:) + phi1(:,:,:,is) / mol_weight(is)
          ENDIF
       ENDDO
       drho1(:,:,:) = drho1(:,:,:) / ta1(:,:,:)  !! XXX ta1 is the inverse molecular weight
    ENDIF

    CALL calc_temp_eos(ta1, rho1, phi1, tb1, xsize(1), xsize(2), xsize(3))
    drho1(:,:,:) = drho1(:,:,:) - dtemp1(:,:,:) / ta1(:,:,:)

    drho1(:,:,:) = rho1(:,:,:) * drho1(:,:,:)

  ENDSUBROUTINE lmn_t_to_rho_trans
  !############################################################################
  !subroutine COR_VEL
  !Correction of u* by the pressure gradient to get a divergence free
  !field
  ! input : px,py,pz
  ! output : ux,uy,uz
  !written by SL 2018
  !############################################################################
  subroutine cor_vel (ux,uy,uz,px,py,pz)

    USE variables
    USE param
    USE mpi

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: px,py,pz
#ifdef DEBG
    real(mytype) :: dep
    integer :: code
#endif

#ifdef DEBG
        dep=maxval(abs(ux))
    call MPI_ALLREDUCE(MPI_IN_PLACE,dep,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## Cor Vel ux', dep
        dep=maxval(abs(uy))
    call MPI_ALLREDUCE(MPI_IN_PLACE,dep,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## Cor Vel uy', dep
        dep=maxval(abs(uz))
    call MPI_ALLREDUCE(MPI_IN_PLACE,dep,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## Cor Vel uz', dep
        dep=maxval(abs(px))
    call MPI_ALLREDUCE(MPI_IN_PLACE,dep,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## Cor Vel px', dep
        dep=maxval(abs(py))
    call MPI_ALLREDUCE(MPI_IN_PLACE,dep,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## Cor Vel py', dep
        dep=maxval(abs(pz))
    call MPI_ALLREDUCE(MPI_IN_PLACE,dep,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## Cor Vel pz', dep
#endif

    ux(:,:,:)=ux(:,:,:)-px(:,:,:)
    uy(:,:,:)=uy(:,:,:)-py(:,:,:)
    uz(:,:,:)=uz(:,:,:)-pz(:,:,:)

    sync_vel_needed = .true.

    return
  end subroutine cor_vel

  !############################################################################
  ! Matrix-free Schur complement correction for inviscid immersed ellipsoids.
  !
  ! This solves, approximately, for a divergence-free velocity correction whose
  ! normal component reduces the immersed no-penetration residual after the
  ! usual pressure projection:
  !
  !     C P E^T lambda = U_body.n - C u
  !
  ! C samples normal velocity at immersed boundary samples, E^T spreads scalar
  ! boundary multipliers back to the nearest velocity cells along the sample
  ! normals, and P is the standard pressure projection.
  !############################################################################
  subroutine ellipsoid_schur_projection(div_visu_var, pp3, px1, py1, pz1, rho1, ux1, uy1, uz1, ep1, drho1, divu3)

    use param, only : ntime, nrhotime, npress, zero, one, dx, dy, dz, xlx, yly, zlz, izap, &
         xnu, ifirst, ilast, itime, itr, istret, sync_vel_needed
    use variables, only : yp, ilist
    use complex_geometry, only : nobjx, nobjy, nobjz, xi, xf, yi, yf, zi, zf
    use ibm_param, only : ellipsoid_schur_projection_iters, ellipsoid_schur_projection_relax, &
         ellipsoid_schur_projection_tol, ellipsoid_schur_projection_passes
    use ellipsoid_utils, only : CalculatePointVelocity_Multi, EllipsoidNormal_Multi
    use var, only : t, nzmsize
    use MPI

    implicit none

    real(mytype), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress), intent(inout) :: pp3
    real(mytype), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress), intent(inout) :: div_visu_var
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(inout) :: ux1, uy1, uz1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(inout) :: px1, py1, pz1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ep1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), nrhotime), intent(in) :: rho1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime), intent(in) :: drho1
    real(mytype), dimension(zsize(1), zsize(2), zsize(3)), intent(in) :: divu3

    integer, allocatable :: ix_sx(:), jy_sx(:), kz_sx(:)
    integer, allocatable :: ix_sy(:), jy_sy(:), kz_sy(:)
    integer, allocatable :: ix_sz(:), jy_sz(:), kz_sz(:)
    real(mytype), allocatable :: nx_sx(:), ny_sx(:), nz_sx(:), body_sx(:), weight_sx(:)
    real(mytype), allocatable :: nx_sy(:), ny_sy(:), nz_sy(:), body_sy(:), weight_sy(:)
    real(mytype), allocatable :: nx_sz(:), ny_sz(:), nz_sz(:), body_sz(:), weight_sz(:)
    real(mytype), allocatable :: sx_b(:), sy_b(:), sz_b(:)
    real(mytype), allocatable :: sx_r(:), sy_r(:), sz_r(:)
    real(mytype), allocatable :: sx_p(:), sy_p(:), sz_p(:)
    real(mytype), allocatable :: sx_ap(:), sy_ap(:), sz_ap(:)
    real(mytype), allocatable :: sx_lambda(:), sy_lambda(:), sz_lambda(:)
    real(mytype), allocatable :: wx(:,:,:), wy(:,:,:), wz(:,:,:)
    integer :: capx, capy, capz, nsx, nsy, nsz, sample_count_local, sample_count_global
    integer :: iter, iter_done, code, iunit, pass, npass
    logical :: file_exists
    real(mytype) :: rr, rr_new, pap, alpha, beta, initial_rms, final_rms, correction_rms
    real(mytype) :: pass1_rms
    real(mytype) :: relax, tol, denom

    if (ellipsoid_schur_projection_iters.le.0) return
    if (xnu.ne.zero) return

    relax = ellipsoid_schur_projection_relax
    tol = ellipsoid_schur_projection_tol
    if (relax.le.zero) return

    capx = max(1, two_int_sum(nobjx))
    capy = max(1, two_int_sum(nobjy))
    capz = max(1, two_int_sum(nobjz))

    call allocate_samples(capx, ix_sx, jy_sx, kz_sx, nx_sx, ny_sx, nz_sx, body_sx, weight_sx)
    call allocate_samples(capy, ix_sy, jy_sy, kz_sy, nx_sy, ny_sy, nz_sy, body_sy, weight_sy)
    call allocate_samples(capz, ix_sz, jy_sz, kz_sz, nx_sz, ny_sz, nz_sz, body_sz, weight_sz)

    nsx = 0
    nsy = 0
    nsz = 0
    call build_schur_samples()

    sample_count_local = nsx + nsy + nsz
    sample_count_global = sample_count_local
    call MPI_Allreduce(MPI_IN_PLACE, sample_count_global, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, code)
    if (sample_count_global.le.0) return

    call allocate_vectors(max(1, nsx), sx_b, sx_r, sx_p, sx_ap, sx_lambda)
    call allocate_vectors(max(1, nsy), sy_b, sy_r, sy_p, sy_ap, sy_lambda)
    call allocate_vectors(max(1, nsz), sz_b, sz_r, sz_p, sz_ap, sz_lambda)

    ! Outer passes: re-sample the CURRENT (already-corrected-by-previous-pass)
    ! velocity field, re-solve, and re-apply, all within this one timestep's
    ! call. Sample GEOMETRY (build_schur_samples above) is body-position-only
    ! and does not change between passes; only the sampled velocity residual
    ! does. With passes=1 this reduces exactly to the original single-pass
    ! behaviour. Motivation: a single small-relax pass per timestep lets
    ! advection re-introduce a fresh (and front/back-asymmetric) boundary
    ! defect faster than one relaxed correction can remove it, so the
    ! correction chases a moving target across many timesteps rather than
    ! converging within one. Iterating the full correction here, before the
    ! next advection step gets a chance to perturb things further, aims to
    ! satisfy the constraint more completely per timestep instead.
    allocate(wx(xsize(1),xsize(2),xsize(3)))
    allocate(wy(xsize(1),xsize(2),xsize(3)))
    allocate(wz(xsize(1),xsize(2),xsize(3)))

    npass = max(1, ellipsoid_schur_projection_passes)
    iter_done = 0
    pass1_rms = zero
    final_rms = zero
    correction_rms = zero

    do pass = 1, npass
       sx_b = zero
       sy_b = zero
       sz_b = zero
       call sample_field(ux1, uy1, uz1, sx_ap, sy_ap, sz_ap)
       if (nsx.gt.0) sx_b(1:nsx) = body_sx(1:nsx) - sx_ap(1:nsx)
       if (nsy.gt.0) sy_b(1:nsy) = body_sy(1:nsy) - sy_ap(1:nsy)
       if (nsz.gt.0) sz_b(1:nsz) = body_sz(1:nsz) - sz_ap(1:nsz)

       sx_r = sx_b
       sy_r = sy_b
       sz_r = sz_b
       sx_p = sx_r
       sy_p = sy_r
       sz_p = sz_r
       sx_lambda = zero
       sy_lambda = zero
       sz_lambda = zero

       rr = schur_dot(sx_r, sy_r, sz_r, sx_r, sy_r, sz_r)
       initial_rms = sqrt(rr / real(sample_count_global, mytype))
       if (pass.eq.1) pass1_rms = initial_rms
       final_rms = initial_rms

       if (initial_rms.gt.zero) then
          do iter = 1, ellipsoid_schur_projection_iters
             call schur_matvec(sx_p, sy_p, sz_p, sx_ap, sy_ap, sz_ap)
             pap = schur_dot(sx_p, sy_p, sz_p, sx_ap, sy_ap, sz_ap)
             if (pap.le.epsilon(one)) exit

             alpha = rr / pap
             if (nsx.gt.0) sx_lambda(1:nsx) = sx_lambda(1:nsx) + alpha * sx_p(1:nsx)
             if (nsy.gt.0) sy_lambda(1:nsy) = sy_lambda(1:nsy) + alpha * sy_p(1:nsy)
             if (nsz.gt.0) sz_lambda(1:nsz) = sz_lambda(1:nsz) + alpha * sz_p(1:nsz)

             if (nsx.gt.0) sx_r(1:nsx) = sx_r(1:nsx) - alpha * sx_ap(1:nsx)
             if (nsy.gt.0) sy_r(1:nsy) = sy_r(1:nsy) - alpha * sy_ap(1:nsy)
             if (nsz.gt.0) sz_r(1:nsz) = sz_r(1:nsz) - alpha * sz_ap(1:nsz)

             rr_new = schur_dot(sx_r, sy_r, sz_r, sx_r, sy_r, sz_r)
             final_rms = sqrt(rr_new / real(sample_count_global, mytype))
             iter_done = iter
             if (final_rms.le.tol * max(initial_rms, epsilon(one))) exit

             denom = max(rr, epsilon(one))
             beta = rr_new / denom
             if (nsx.gt.0) sx_p(1:nsx) = sx_r(1:nsx) + beta * sx_p(1:nsx)
             if (nsy.gt.0) sy_p(1:nsy) = sy_r(1:nsy) + beta * sy_p(1:nsy)
             if (nsz.gt.0) sz_p(1:nsz) = sz_r(1:nsz) + beta * sz_p(1:nsz)
             rr = rr_new
          enddo
       endif

       call spread_vector(sx_lambda, sy_lambda, sz_lambda, wx, wy, wz)
       call project_correction(wx, wy, wz)
       call sample_field(wx, wy, wz, sx_ap, sy_ap, sz_ap)

       if (nsx.gt.0) sx_r(1:nsx) = sx_b(1:nsx) - relax * sx_ap(1:nsx)
       if (nsy.gt.0) sy_r(1:nsy) = sy_b(1:nsy) - relax * sy_ap(1:nsy)
       if (nsz.gt.0) sz_r(1:nsz) = sz_b(1:nsz) - relax * sz_ap(1:nsz)
       final_rms = sqrt(schur_dot(sx_r, sy_r, sz_r, sx_r, sy_r, sz_r) / real(sample_count_global, mytype))
       correction_rms = sqrt(schur_dot(sx_ap, sy_ap, sz_ap, sx_ap, sy_ap, sz_ap) / &
            real(sample_count_global, mytype))

       ux1(:,:,:) = ux1(:,:,:) + relax * wx(:,:,:)
       uy1(:,:,:) = uy1(:,:,:) + relax * wy(:,:,:)
       uz1(:,:,:) = uz1(:,:,:) + relax * wz(:,:,:)
       sync_vel_needed = .true.
    enddo

    initial_rms = pass1_rms

    deallocate(wx, wy, wz)

    if (nrank.eq.0 .and. (mod(itime,ilist).eq.0 .or. itime.eq.ifirst .or. itime.eq.ilast)) then
       write(*,*) "Ellipsoid Schur projection: samples=", sample_count_global, &
            " passes=", npass, " iterations=", iter_done, " residual rms initial/final=", initial_rms, final_rms, &
            " correction rms=", correction_rms, " relax=", relax

       inquire(file="ellipsoid_schur_projection.dat", exist=file_exists)
       open(newunit=iunit, file="ellipsoid_schur_projection.dat", status="unknown", position="append")
       if (.not. file_exists) then
          write(iunit,*) "# t itime itr samples iterations initial_rms final_rms correction_rms relax tol xnu"
       endif
       write(iunit,*) t, itime, itr, sample_count_global, iter_done, initial_rms, final_rms, &
            correction_rms, relax, tol, xnu
       close(iunit)
    endif

  contains

    integer function two_int_sum(values)
      integer, intent(in), dimension(:,:) :: values
      two_int_sum = 2 * sum(values) + 1
    end function two_int_sum

    subroutine allocate_samples(capacity, ix, jy, kz, nxn, nyn, nzn, body, weight)
      integer, intent(in) :: capacity
      integer, allocatable, intent(out) :: ix(:), jy(:), kz(:)
      real(mytype), allocatable, intent(out) :: nxn(:), nyn(:), nzn(:), body(:), weight(:)

      allocate(ix(capacity), jy(capacity), kz(capacity))
      allocate(nxn(capacity), nyn(capacity), nzn(capacity), body(capacity), weight(capacity))
      ix = 0
      jy = 0
      kz = 0
      nxn = zero
      nyn = zero
      nzn = zero
      body = zero
      weight = one
    end subroutine allocate_samples

    subroutine allocate_vectors(n, b, r, p, ap, lambda)
      integer, intent(in) :: n
      real(mytype), allocatable, intent(out) :: b(:), r(:), p(:), ap(:), lambda(:)

      allocate(b(n), r(n), p(n), ap(n), lambda(n))
      b = zero
      r = zero
      p = zero
      ap = zero
      lambda = zero
    end subroutine allocate_vectors

    subroutine build_schur_samples()
      integer :: i, j, k, iobj, ix, jy, kz
      real(mytype) :: point(3)

      do k = 1, xsize(3)
         do j = 1, xsize(2)
            do iobj = 1, nobjx(j,k)
               if (xi(iobj,j,k).gt.zero) then
                  ix = int(xi(iobj,j,k) / dx) + 1
                  if (izap.eq.1) ix = ix - 1
                  if (ix.ge.1 .and. ix.le.xsize(1)) then
                     point = [xi(iobj,j,k), y_x_coord(j), z_x_coord(k)]
                     call add_sample_x(point, ix, j, k, dy*dz)
                  endif
               endif

               if (xf(iobj,j,k).lt.xlx) then
                  ix = int((xf(iobj,j,k) + dx) / dx) + 1
                  if (izap.eq.1) ix = ix + 1
                  if (ix.ge.1 .and. ix.le.xsize(1)) then
                     point = [xf(iobj,j,k), y_x_coord(j), z_x_coord(k)]
                     call add_sample_x(point, ix, j, k, dy*dz)
                  endif
               endif
            enddo
         enddo
      enddo

      do k = 1, ysize(3)
         do i = 1, ysize(1)
            do j = 1, nobjy(i,k)
               if (yi(j,i,k).gt.zero) then
                  jy = lower_velocity_y_index(yi(j,i,k))
                  if (izap.eq.1) jy = jy - 1
                  if (jy.ge.1 .and. jy.le.ysize(2)) then
                     point = [x_y_coord(i), yi(j,i,k), z_y_coord(k)]
                     call add_sample_y(point, i, jy, k, dx*dz)
                  endif
               endif

               if (yf(j,i,k).lt.yly) then
                  jy = upper_velocity_y_index(yf(j,i,k))
                  if (izap.eq.1) jy = jy + 1
                  if (jy.ge.1 .and. jy.le.ysize(2)) then
                     point = [x_y_coord(i), yf(j,i,k), z_y_coord(k)]
                     call add_sample_y(point, i, jy, k, dx*dz)
                  endif
               endif
            enddo
         enddo
      enddo

      do j = 1, zsize(2)
         do i = 1, zsize(1)
            do k = 1, nobjz(i,j)
               if (zi(k,i,j).gt.zero) then
                  kz = int(zi(k,i,j) / dz) + 1
                  if (izap.eq.1) kz = kz - 1
                  if (kz.ge.1 .and. kz.le.zsize(3)) then
                     point = [x_z_coord(i), y_z_coord(j), zi(k,i,j)]
                     call add_sample_z(point, i, j, kz, dx*dy)
                  endif
               endif

               if (zf(k,i,j).lt.zlz) then
                  kz = int((zf(k,i,j) + dz) / dz) + 1
                  if (izap.eq.1) kz = kz + 1
                  if (kz.ge.1 .and. kz.le.zsize(3)) then
                     point = [x_z_coord(i), y_z_coord(j), zf(k,i,j)]
                     call add_sample_z(point, i, j, kz, dx*dy)
                  endif
               endif
            enddo
         enddo
      enddo
    end subroutine build_schur_samples

    subroutine add_sample_x(point, ix, jy, kz, projected_area)
      real(mytype), intent(in) :: point(3), projected_area
      integer, intent(in) :: ix, jy, kz
      real(mytype) :: normal(3), body_velocity(3), normal_axis_abs

      call CalculatePointVelocity_Multi(point, body_velocity)
      call EllipsoidNormal_Multi(point, normal)
      if (dominant_normal_axis(normal).ne.1) return
      normal_axis_abs = abs(normal(1))
      if (normal_axis_abs.le.zero) return
      nsx = nsx + 1
      ix_sx(nsx) = ix
      jy_sx(nsx) = jy
      kz_sx(nsx) = kz
      nx_sx(nsx) = normal(1)
      ny_sx(nsx) = normal(2)
      nz_sx(nsx) = normal(3)
      body_sx(nsx) = sum(body_velocity * normal)
      weight_sx(nsx) = projected_area / normal_axis_abs
    end subroutine add_sample_x

    subroutine add_sample_y(point, ix, jy, kz, projected_area)
      real(mytype), intent(in) :: point(3), projected_area
      integer, intent(in) :: ix, jy, kz
      real(mytype) :: normal(3), body_velocity(3), normal_axis_abs

      call CalculatePointVelocity_Multi(point, body_velocity)
      call EllipsoidNormal_Multi(point, normal)
      if (dominant_normal_axis(normal).ne.2) return
      normal_axis_abs = abs(normal(2))
      if (normal_axis_abs.le.zero) return
      nsy = nsy + 1
      ix_sy(nsy) = ix
      jy_sy(nsy) = jy
      kz_sy(nsy) = kz
      nx_sy(nsy) = normal(1)
      ny_sy(nsy) = normal(2)
      nz_sy(nsy) = normal(3)
      body_sy(nsy) = sum(body_velocity * normal)
      weight_sy(nsy) = projected_area / normal_axis_abs
    end subroutine add_sample_y

    subroutine add_sample_z(point, ix, jy, kz, projected_area)
      real(mytype), intent(in) :: point(3), projected_area
      integer, intent(in) :: ix, jy, kz
      real(mytype) :: normal(3), body_velocity(3), normal_axis_abs

      call CalculatePointVelocity_Multi(point, body_velocity)
      call EllipsoidNormal_Multi(point, normal)
      if (dominant_normal_axis(normal).ne.3) return
      normal_axis_abs = abs(normal(3))
      if (normal_axis_abs.le.zero) return
      nsz = nsz + 1
      ix_sz(nsz) = ix
      jy_sz(nsz) = jy
      kz_sz(nsz) = kz
      nx_sz(nsz) = normal(1)
      ny_sz(nsz) = normal(2)
      nz_sz(nsz) = normal(3)
      body_sz(nsz) = sum(body_velocity * normal)
      weight_sz(nsz) = projected_area / normal_axis_abs
    end subroutine add_sample_z

    integer function dominant_normal_axis(normal)
      real(mytype), intent(in) :: normal(3)
      real(mytype) :: ax, ay, az

      ax = abs(normal(1))
      ay = abs(normal(2))
      az = abs(normal(3))
      if (ax.ge.ay .and. ax.ge.az) then
         dominant_normal_axis = 1
      elseif (ay.ge.ax .and. ay.ge.az) then
         dominant_normal_axis = 2
      else
         dominant_normal_axis = 3
      endif
    end function dominant_normal_axis

    subroutine spread_vector(vx, vy, vz, qx, qy, qz)
      real(mytype), intent(in) :: vx(:), vy(:), vz(:)
      real(mytype), intent(out), dimension(xsize(1),xsize(2),xsize(3)) :: qx, qy, qz
      real(mytype), allocatable :: qx2(:,:,:), qy2(:,:,:), qz2(:,:,:)
      real(mytype), allocatable :: qx3(:,:,:), qy3(:,:,:), qz3(:,:,:)
      integer :: s, d, ii, jj, kk
      real(mytype) :: w

      ! Each sample is spread with a 3-point (1/4,1/2,1/4) hat smoothing along
      ! its dominant-normal axis (the boundary's fractional direction, contiguous
      ! within that sample group's pencil). This regularises the single-cell
      ! delta so the projected correction is smooth and full-strength enforcement
      ! does not create near-body velocity spikes. sample_field is the adjoint.
      qx = zero
      qy = zero
      qz = zero
      do s = 1, nsx
         do d = -1, 1
            ii = ix_sx(s) + d
            if (ii.lt.1 .or. ii.gt.xsize(1)) cycle
            w = merge(0.5_mytype, 0.25_mytype, d.eq.0)
            qx(ii,jy_sx(s),kz_sx(s)) = qx(ii,jy_sx(s),kz_sx(s)) + w * vx(s) * nx_sx(s)
            qy(ii,jy_sx(s),kz_sx(s)) = qy(ii,jy_sx(s),kz_sx(s)) + w * vx(s) * ny_sx(s)
            qz(ii,jy_sx(s),kz_sx(s)) = qz(ii,jy_sx(s),kz_sx(s)) + w * vx(s) * nz_sx(s)
         enddo
      enddo

      allocate(qx2(ysize(1),ysize(2),ysize(3)), qy2(ysize(1),ysize(2),ysize(3)), qz2(ysize(1),ysize(2),ysize(3)))
      call transpose_x_to_y(qx, qx2)
      call transpose_x_to_y(qy, qy2)
      call transpose_x_to_y(qz, qz2)
      do s = 1, nsy
         do d = -1, 1
            jj = jy_sy(s) + d
            if (jj.lt.1 .or. jj.gt.ysize(2)) cycle
            w = merge(0.5_mytype, 0.25_mytype, d.eq.0)
            qx2(ix_sy(s),jj,kz_sy(s)) = qx2(ix_sy(s),jj,kz_sy(s)) + w * vy(s) * nx_sy(s)
            qy2(ix_sy(s),jj,kz_sy(s)) = qy2(ix_sy(s),jj,kz_sy(s)) + w * vy(s) * ny_sy(s)
            qz2(ix_sy(s),jj,kz_sy(s)) = qz2(ix_sy(s),jj,kz_sy(s)) + w * vy(s) * nz_sy(s)
         enddo
      enddo

      allocate(qx3(zsize(1),zsize(2),zsize(3)), qy3(zsize(1),zsize(2),zsize(3)), qz3(zsize(1),zsize(2),zsize(3)))
      call transpose_y_to_z(qx2, qx3)
      call transpose_y_to_z(qy2, qy3)
      call transpose_y_to_z(qz2, qz3)
      deallocate(qx2, qy2, qz2)

      do s = 1, nsz
         do d = -1, 1
            kk = kz_sz(s) + d
            if (kk.lt.1 .or. kk.gt.zsize(3)) cycle
            w = merge(0.5_mytype, 0.25_mytype, d.eq.0)
            qx3(ix_sz(s),jy_sz(s),kk) = qx3(ix_sz(s),jy_sz(s),kk) + w * vz(s) * nx_sz(s)
            qy3(ix_sz(s),jy_sz(s),kk) = qy3(ix_sz(s),jy_sz(s),kk) + w * vz(s) * ny_sz(s)
            qz3(ix_sz(s),jy_sz(s),kk) = qz3(ix_sz(s),jy_sz(s),kk) + w * vz(s) * nz_sz(s)
         enddo
      enddo

      allocate(qx2(ysize(1),ysize(2),ysize(3)), qy2(ysize(1),ysize(2),ysize(3)), qz2(ysize(1),ysize(2),ysize(3)))
      call transpose_z_to_y(qx3, qx2)
      call transpose_z_to_y(qy3, qy2)
      call transpose_z_to_y(qz3, qz2)
      deallocate(qx3, qy3, qz3)
      call transpose_y_to_x(qx2, qx)
      call transpose_y_to_x(qy2, qy)
      call transpose_y_to_x(qz2, qz)
      deallocate(qx2, qy2, qz2)
    end subroutine spread_vector

    subroutine sample_field(qx, qy, qz, vx, vy, vz)
      real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: qx, qy, qz
      real(mytype), intent(out) :: vx(:), vy(:), vz(:)
      real(mytype), allocatable :: qx2(:,:,:), qy2(:,:,:), qz2(:,:,:)
      real(mytype), allocatable :: qx3(:,:,:), qy3(:,:,:), qz3(:,:,:)
      integer :: s, d, ii, jj, kk
      real(mytype) :: w, acc

      ! Adjoint of spread_vector: 3-point (1/4,1/2,1/4) weighted sample of the
      ! normal velocity along each sample's dominant axis.
      vx = zero
      vy = zero
      vz = zero
      do s = 1, nsx
         acc = zero
         do d = -1, 1
            ii = ix_sx(s) + d
            if (ii.lt.1 .or. ii.gt.xsize(1)) cycle
            w = merge(0.5_mytype, 0.25_mytype, d.eq.0)
            acc = acc + w * (qx(ii,jy_sx(s),kz_sx(s)) * nx_sx(s) + &
                 qy(ii,jy_sx(s),kz_sx(s)) * ny_sx(s) + &
                 qz(ii,jy_sx(s),kz_sx(s)) * nz_sx(s))
         enddo
         vx(s) = acc
      enddo

      allocate(qx2(ysize(1),ysize(2),ysize(3)), qy2(ysize(1),ysize(2),ysize(3)), qz2(ysize(1),ysize(2),ysize(3)))
      call transpose_x_to_y(qx, qx2)
      call transpose_x_to_y(qy, qy2)
      call transpose_x_to_y(qz, qz2)
      do s = 1, nsy
         acc = zero
         do d = -1, 1
            jj = jy_sy(s) + d
            if (jj.lt.1 .or. jj.gt.ysize(2)) cycle
            w = merge(0.5_mytype, 0.25_mytype, d.eq.0)
            acc = acc + w * (qx2(ix_sy(s),jj,kz_sy(s)) * nx_sy(s) + &
                 qy2(ix_sy(s),jj,kz_sy(s)) * ny_sy(s) + &
                 qz2(ix_sy(s),jj,kz_sy(s)) * nz_sy(s))
         enddo
         vy(s) = acc
      enddo

      allocate(qx3(zsize(1),zsize(2),zsize(3)), qy3(zsize(1),zsize(2),zsize(3)), qz3(zsize(1),zsize(2),zsize(3)))
      call transpose_y_to_z(qx2, qx3)
      call transpose_y_to_z(qy2, qy3)
      call transpose_y_to_z(qz2, qz3)
      deallocate(qx2, qy2, qz2)
      do s = 1, nsz
         acc = zero
         do d = -1, 1
            kk = kz_sz(s) + d
            if (kk.lt.1 .or. kk.gt.zsize(3)) cycle
            w = merge(0.5_mytype, 0.25_mytype, d.eq.0)
            acc = acc + w * (qx3(ix_sz(s),jy_sz(s),kk) * nx_sz(s) + &
                 qy3(ix_sz(s),jy_sz(s),kk) * ny_sz(s) + &
                 qz3(ix_sz(s),jy_sz(s),kk) * nz_sz(s))
         enddo
         vz(s) = acc
      enddo
      deallocate(qx3, qy3, qz3)
    end subroutine sample_field

    subroutine schur_matvec(vx, vy, vz, ax, ay, az)
      real(mytype), intent(in) :: vx(:), vy(:), vz(:)
      real(mytype), intent(out) :: ax(:), ay(:), az(:)
      real(mytype), allocatable :: qx(:,:,:), qy(:,:,:), qz(:,:,:)

      allocate(qx(xsize(1),xsize(2),xsize(3)))
      allocate(qy(xsize(1),xsize(2),xsize(3)))
      allocate(qz(xsize(1),xsize(2),xsize(3)))
      call spread_vector(vx, vy, vz, qx, qy, qz)
      call project_correction(qx, qy, qz)
      call sample_field(qx, qy, qz, ax, ay, az)
      deallocate(qx, qy, qz)
    end subroutine schur_matvec

    subroutine project_correction(qx, qy, qz)
      real(mytype), intent(inout), dimension(xsize(1),xsize(2),xsize(3)) :: qx, qy, qz

      call solve_poisson(div_visu_var, pp3, px1, py1, pz1, rho1, qx, qy, qz, ep1, drho1, divu3, &
           save_bc_grad=.false., quiet=.true.)
      call cor_vel(qx, qy, qz, px1, py1, pz1)
    end subroutine project_correction

    real(mytype) function schur_dot(ax, ay, az, bx, by, bz)
      real(mytype), intent(in) :: ax(:), ay(:), az(:), bx(:), by(:), bz(:)
      integer :: code_dot

      schur_dot = zero
      if (nsx.gt.0) schur_dot = schur_dot + sum(ax(1:nsx) * bx(1:nsx))
      if (nsy.gt.0) schur_dot = schur_dot + sum(ay(1:nsy) * by(1:nsy))
      if (nsz.gt.0) schur_dot = schur_dot + sum(az(1:nsz) * bz(1:nsz))
      call MPI_Allreduce(MPI_IN_PLACE, schur_dot, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code_dot)
    end function schur_dot

    real(mytype) function y_x_coord(jloc)
      integer, intent(in) :: jloc
      integer :: jglob
      jglob = xstart(2) + jloc - 1
      if (istret.eq.0) then
         y_x_coord = real(jglob - 1, mytype) * dy
      else
         y_x_coord = yp(jglob)
      endif
    end function y_x_coord

    real(mytype) function z_x_coord(kloc)
      integer, intent(in) :: kloc
      z_x_coord = real(xstart(3) + kloc - 2, mytype) * dz
    end function z_x_coord

    real(mytype) function x_y_coord(iloc)
      integer, intent(in) :: iloc
      x_y_coord = real(ystart(1) + iloc - 2, mytype) * dx
    end function x_y_coord

    real(mytype) function z_y_coord(kloc)
      integer, intent(in) :: kloc
      z_y_coord = real(ystart(3) + kloc - 2, mytype) * dz
    end function z_y_coord

    real(mytype) function x_z_coord(iloc)
      integer, intent(in) :: iloc
      x_z_coord = real(zstart(1) + iloc - 2, mytype) * dx
    end function x_z_coord

    real(mytype) function y_z_coord(jloc)
      integer, intent(in) :: jloc
      integer :: jglob
      jglob = zstart(2) + jloc - 1
      if (istret.eq.0) then
         y_z_coord = real(jglob - 1, mytype) * dy
      else
         y_z_coord = yp(jglob)
      endif
    end function y_z_coord

    integer function lower_velocity_y_index(ypos)
      real(mytype), intent(in) :: ypos
      integer :: jj

      lower_velocity_y_index = 1
      do jj = 1, size(yp)
         if (yp(jj).lt.ypos) lower_velocity_y_index = jj
      enddo
    end function lower_velocity_y_index

    integer function upper_velocity_y_index(ypos)
      real(mytype), intent(in) :: ypos
      integer :: jj

      upper_velocity_y_index = size(yp)
      do jj = 1, size(yp)
         if (yp(jj).gt.ypos) then
            upper_velocity_y_index = jj
            return
         endif
      enddo
    end function upper_velocity_y_index

  end subroutine ellipsoid_schur_projection
  !############################################################################
  !subroutine DIVERGENCe
  !Calculation of div u* for nlock=1 and of div u^{n+1} for nlock=2
  ! input : ux1,uy1,uz1,ep1 (on velocity mesh)
  ! output : pp3 (on pressure mesh)
  !written by SL 2018
  !############################################################################
  subroutine divergence (div_visu_var,pp3,rho1,ux1,uy1,uz1,ep1,drho1,divu3,nlock)

    USE param
    USE variables
    USE var, ONLY: ta1, tb1, tc1, pp1, pgy1, pgz1, di1, &
         duxdxp2, uyp2, uzp2, duydypi2, upi2, ta2, dipp2, &
         duxydxyp3, uzp3, po3, dipp3, nxmsize, nymsize, nzmsize, &
         ux2, uy2, uz2, ux3, uy3, uz3
    USE MPI
    USE ibm_param
    USE ellipsoid_utils, ONLY: navierFieldGen
    USE ellip, ONLY: ellipsoid_projection_flux_rhs_x, ellipsoid_projection_flux_rhs_y, &
         ellipsoid_projection_flux_rhs_z

    implicit none

    !  TYPE(DECOMP_INFO) :: ph1,ph3,ph4

    !X PENCILS NX NY NZ  -->NXM NY NZ
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime),intent(in) :: drho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime),intent(in) :: rho1
    !Z PENCILS NXM NYM NZ  -->NXM NYM NZM
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)),intent(in) :: divu3
    real(mytype),dimension(xsize(1),xsize(2),xsize(3))     :: ep1_ux,ep1_uy,ep1_uz
    real(mytype),dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize),intent(out) :: pp3,div_visu_var
   !  real(mytype),dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize) :: div_visu_var


    integer :: nvect3,i,j,k,nlock
    integer :: code
    real(mytype) :: tmax,tmoy,pres_ref

    nvect3=(ph1%zen(1)-ph1%zst(1)+1)*(ph1%zen(2)-ph1%zst(2)+1)*nzmsize

    if (iibm.eq.0) then
       ta1(:,:,:) = ux1(:,:,:)
       tb1(:,:,:) = uy1(:,:,:)
       tc1(:,:,:) = uz1(:,:,:)
    else if (itype.eq.itype_ellip) then
       call navierFieldGen(ep1, ep1_ux, ep1_uy, ep1_uz, ux1, uy1, uz1)
       ta1(:,:,:) = (one - ep1(:,:,:)) * ux1(:,:,:) + ep1(:,:,:)*ep1_ux(:,:,:)
       tb1(:,:,:) = (one - ep1(:,:,:)) * uy1(:,:,:) + ep1(:,:,:)*ep1_uy(:,:,:)
       tc1(:,:,:) = (one - ep1(:,:,:)) * uz1(:,:,:) + ep1(:,:,:)*ep1_uz(:,:,:)
    else
       ta1(:,:,:) = (one - ep1(:,:,:)) * ux1(:,:,:) + ep1(:,:,:)*ubcx
       tb1(:,:,:) = (one - ep1(:,:,:)) * uy1(:,:,:) + ep1(:,:,:)*ubcy
       tc1(:,:,:) = (one - ep1(:,:,:)) * uz1(:,:,:) + ep1(:,:,:)*ubcz
    endif

    !WORK X-PENCILS

    call derxvp(pp1,ta1,di1,sx,cfx6,csx6,cwx6,xsize(1),nxmsize,xsize(2),xsize(3),0)

    if (itype.eq.itype_ellip .and. xnu.eq.zero .and. ellipsoid_projection_flux_fix.gt.0) then
       call ellipsoid_projection_flux_rhs_x(pp1, ux1, uy1, uz1)
    endif

    if (ilmn.and.(nlock.gt.0)) then
       if ((nlock.eq.1).and.(.not.ivarcoeff)) then
          !! Approximate -div(rho u) using ddt(rho)
          call extrapol_drhodt(ta1, rho1, drho1)
       elseif ((nlock.eq.2).or.ivarcoeff) then
          !! Need to check our error against divu constraint
          !! Or else we are solving the variable-coefficient Poisson equation
          call transpose_z_to_y(-divu3, ta2)
          call transpose_y_to_x(ta2, ta1)
       endif
       call interxvp(pgy1,ta1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)
       pp1(:,:,:) = pp1(:,:,:) + pgy1(:,:,:)
    endif

    call interxvp(pgy1,tb1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)
    call interxvp(pgz1,tc1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)

    call transpose_x_to_y(pp1,duxdxp2,ph4)!->NXM NY NZ
    call transpose_x_to_y(pgy1,uyp2,ph4)
    call transpose_x_to_y(pgz1,uzp2,ph4)

    !WORK Y-PENCILS
    call interyvp(upi2,duxdxp2,dipp2,sy,cifyp6,cisyp6,ciwyp6,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),1)
    call deryvp(duydypi2,uyp2,dipp2,sy,cfy6,csy6,cwy6,ppyi,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),0)

    !! Compute sum dudx + dvdy
    duydypi2(:,:,:) = duydypi2(:,:,:) + upi2(:,:,:)

    if (itype.eq.itype_ellip .and. xnu.eq.zero .and. ellipsoid_projection_flux_fix.gt.0) then
       call transpose_x_to_y(ux1, ux2)
       call transpose_x_to_y(uy1, uy2)
       call transpose_x_to_y(uz1, uz2)
       call ellipsoid_projection_flux_rhs_y(duydypi2, ux2, uy2, uz2)
    endif

    call interyvp(upi2,uzp2,dipp2,sy,cifyp6,cisyp6,ciwyp6,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),1)

    call transpose_y_to_z(duydypi2,duxydxyp3,ph3)!->NXM NYM NZ
    call transpose_y_to_z(upi2,uzp3,ph3)

    !WORK Z-PENCILS
    call interzvp(pp3,duxydxyp3,dipp3,sz,cifzp6,ciszp6,ciwzp6,(ph1%zen(1)-ph1%zst(1)+1),&
         (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzmsize,1)
    call derzvp(po3,uzp3,dipp3,sz,cfz6,csz6,cwz6,(ph1%zen(1)-ph1%zst(1)+1),&
         (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzmsize,0)

    !! Compute sum dudx + dvdy + dwdz
    pp3(:,:,:) = pp3(:,:,:) + po3(:,:,:)

    if (itype.eq.itype_ellip .and. xnu.eq.zero .and. ellipsoid_projection_flux_fix.gt.0) then
       call transpose_y_to_z(ux2, ux3)
       call transpose_y_to_z(uy2, uy3)
       call transpose_y_to_z(uz2, uz3)
       call ellipsoid_projection_flux_rhs_z(pp3, ux3, uy3, uz3)
    endif

    if (nlock==2) then
       ! Line below sometimes generates issues with Intel
       !pp3(:,:,:)=pp3(:,:,:)-pp3(ph1%zst(1),ph1%zst(2),nzmsize)
       ! Using a tmp variable seems to sort the issue
       pres_ref = pp3(ph1%zst(1),ph1%zst(2),nzmsize)
       pp3(:,:,:)=pp3(:,:,:)-pres_ref
    endif

    tmax=-1609._mytype
    tmoy=zero
    do k=1,nzmsize
       do j=ph1%zst(2),ph1%zen(2)
          do i=ph1%zst(1),ph1%zen(1)
             if (pp3(i,j,k).gt.tmax) tmax=pp3(i,j,k)
             tmoy=tmoy+abs(pp3(i,j,k))
          enddo
       enddo
    enddo
    tmoy=tmoy/real(nvect3,mytype)
    call MPI_ALLREDUCE(MPI_IN_PLACE, tmoy, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
    call MPI_ALLREDUCE(MPI_IN_PLACE, tmax, 1, real_type, MPI_MAX, MPI_COMM_WORLD, code)

    if ((nrank == 0) .and. (nlock > 0) .and. (.not. suppress_divergence_output) .and. &
         (mod(itime, ilist) == 0 .or. itime == ifirst .or. itime==ilast)) then
       if (nlock == 2) then
          write(*,*) 'DIV U  max mean=',real(tmax,mytype),real(tmoy/real(nproc),mytype)
       else
          write(*,*) 'DIV U* max mean=',real(tmax,mytype),real(tmoy/real(nproc),mytype)
       endif
    endif

    div_visu_var(:,:,:)=pp3(:,:,:)

    return
  end subroutine divergence

  subroutine divergence2(pp3,rho1,ux1,uy1,uz1,ep1,drho1,divu3,nlock)

   USE param
   USE variables
   USE var, ONLY: ta1, tb1, tc1, pp1, pgy1, pgz1, di1, &
        duxdxp2, uyp2, uzp2, duydypi2, upi2, ta2, dipp2, &
        duxydxyp3, uzp3, po3, dipp3, nxmsize, nymsize, nzmsize
   USE MPI
   USE ibm_param
   USE ellipsoid_utils, ONLY: navierFieldGen

   implicit none

   !  TYPE(DECOMP_INFO) :: ph1,ph3,ph4

   !X PENCILS NX NY NZ  -->NXM NY NZ
   real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1,uz1,ep1
   real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime),intent(in) :: drho1
   real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime),intent(in) :: rho1
   !Z PENCILS NXM NYM NZ  -->NXM NYM NZM
   real(mytype),dimension(zsize(1),zsize(2),zsize(3)),intent(in) :: divu3
   real(mytype),dimension(xsize(1),xsize(2),xsize(3))     :: ep1_ux,ep1_uy,ep1_uz
   real(mytype),dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize),intent(out) :: pp3
   real(mytype),dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize) :: div_visu_var


   integer :: nvect3,i,j,k,nlock
   integer :: code
   real(mytype) :: tmax,tmoy,pres_ref

   nvect3=(ph1%zen(1)-ph1%zst(1)+1)*(ph1%zen(2)-ph1%zst(2)+1)*nzmsize

   if (iibm.eq.0) then
      ta1(:,:,:) = ux1(:,:,:)
      tb1(:,:,:) = uy1(:,:,:)
      tc1(:,:,:) = uz1(:,:,:)
   else if (itype.eq.itype_ellip) then
      call navierFieldGen(ep1, ep1_ux, ep1_uy, ep1_uz, ux1, uy1, uz1)
      ta1(:,:,:) = (one - ep1(:,:,:)) * ux1(:,:,:) + ep1(:,:,:)*ep1_ux(:,:,:)
      tb1(:,:,:) = (one - ep1(:,:,:)) * uy1(:,:,:) + ep1(:,:,:)*ep1_uy(:,:,:)
      tc1(:,:,:) = (one - ep1(:,:,:)) * uz1(:,:,:) + ep1(:,:,:)*ep1_uz(:,:,:)
   else
      ta1(:,:,:) = (one - ep1(:,:,:)) * ux1(:,:,:) + ep1(:,:,:)*ubcx
      tb1(:,:,:) = (one - ep1(:,:,:)) * uy1(:,:,:) + ep1(:,:,:)*ubcy
      tc1(:,:,:) = (one - ep1(:,:,:)) * uz1(:,:,:) + ep1(:,:,:)*ubcz
   endif

   !WORK X-PENCILS

   call derxvp(pp1,ta1,di1,sx,cfx6,csx6,cwx6,xsize(1),nxmsize,xsize(2),xsize(3),0)

   if (ilmn.and.(nlock.gt.0)) then
      if ((nlock.eq.1).and.(.not.ivarcoeff)) then
         !! Approximate -div(rho u) using ddt(rho)
         call extrapol_drhodt(ta1, rho1, drho1)
      elseif ((nlock.eq.2).or.ivarcoeff) then
         !! Need to check our error against divu constraint
         !! Or else we are solving the variable-coefficient Poisson equation
         call transpose_z_to_y(-divu3, ta2)
         call transpose_y_to_x(ta2, ta1)
      endif
      call interxvp(pgy1,ta1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)
      pp1(:,:,:) = pp1(:,:,:) + pgy1(:,:,:)
   endif

   call interxvp(pgy1,tb1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)
   call interxvp(pgz1,tc1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)

   call transpose_x_to_y(pp1,duxdxp2,ph4)!->NXM NY NZ
   call transpose_x_to_y(pgy1,uyp2,ph4)
   call transpose_x_to_y(pgz1,uzp2,ph4)

   !WORK Y-PENCILS
   call interyvp(upi2,duxdxp2,dipp2,sy,cifyp6,cisyp6,ciwyp6,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),1)
   call deryvp(duydypi2,uyp2,dipp2,sy,cfy6,csy6,cwy6,ppyi,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),0)

   !! Compute sum dudx + dvdy
   duydypi2(:,:,:) = duydypi2(:,:,:) + upi2(:,:,:)

   call interyvp(upi2,uzp2,dipp2,sy,cifyp6,cisyp6,ciwyp6,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),1)

   call transpose_y_to_z(duydypi2,duxydxyp3,ph3)!->NXM NYM NZ
   call transpose_y_to_z(upi2,uzp3,ph3)

   !WORK Z-PENCILS
   call interzvp(pp3,duxydxyp3,dipp3,sz,cifzp6,ciszp6,ciwzp6,(ph1%zen(1)-ph1%zst(1)+1),&
        (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzmsize,1)
   call derzvp(po3,uzp3,dipp3,sz,cfz6,csz6,cwz6,(ph1%zen(1)-ph1%zst(1)+1),&
        (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzmsize,0)

   !! Compute sum dudx + dvdy + dwdz
   pp3(:,:,:) = pp3(:,:,:) + po3(:,:,:)

   if (nlock==2) then
      ! Line below sometimes generates issues with Intel
      !pp3(:,:,:)=pp3(:,:,:)-pp3(ph1%zst(1),ph1%zst(2),nzmsize)
      ! Using a tmp variable seems to sort the issue
      pres_ref = pp3(ph1%zst(1),ph1%zst(2),nzmsize)
      pp3(:,:,:)=pp3(:,:,:)-pres_ref
   endif

   tmax=-1609._mytype
   tmoy=zero
   do k=1,nzmsize
      do j=ph1%zst(2),ph1%zen(2)
         do i=ph1%zst(1),ph1%zen(1)
            if (pp3(i,j,k).gt.tmax) tmax=pp3(i,j,k)
            tmoy=tmoy+abs(pp3(i,j,k))
         enddo
      enddo
   enddo
   tmoy=tmoy/real(nvect3,mytype)
   call MPI_ALLREDUCE(MPI_IN_PLACE, tmoy, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
   call MPI_ALLREDUCE(MPI_IN_PLACE, tmax, 1, real_type, MPI_MAX, MPI_COMM_WORLD, code)

   if ((nrank == 0) .and. (nlock > 0) .and. (.not. suppress_divergence_output) .and. &
        (mod(itime, ilist) == 0 .or. itime == ifirst .or. itime==ilast)) then
      if (nlock == 2) then
         write(*,*) 'DIV U  max mean=',real(tmax,mytype),real(tmoy/real(nproc),mytype)
      else
         write(*,*) 'DIV U* max mean=',real(tmax,mytype),real(tmoy/real(nproc),mytype)
      endif
   endif

   div_visu_var(:,:,:)=pp3(:,:,:)

   return
 end subroutine divergence2
  !############################################################################
  !subroutine GRADP
  !Computation of the pressure gradient from the pressure mesh to the
  !velocity mesh
  !Saving pressure gradients on boundaries for correct imposition of
  !BCs on u* via the fractional step methodi (it is not possible to
  !impose BC after correction by pressure gradient otherwise lost of
  !incompressibility--> BCs are imposed on u*
  !
  ! input: pp3 - pressure field (on pressure mesh)
  ! output: px1, py1, pz1 - pressure gradients (on velocity mesh)
  !written by SL 2018
  !############################################################################
  subroutine gradp(px1,py1,pz1,pp3,save_bc_grad)

    USE param
    USE variables
    USE MPI
    USE var, only: pp1,pgy1,pgz1,di1,pp2,ppi2,pgy2,pgz2,pgzi2,dip2,&
         pgz3,ppi3,dip3,nxmsize,nymsize,nzmsize

    USE forces, only : iforces, ppi1

    implicit none

    logical, intent(in), optional :: save_bc_grad
    integer :: i,j,k
    logical :: save_bc_grad_local

    real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize) :: pp3
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: px1,py1,pz1

    save_bc_grad_local = .true.
    if (present(save_bc_grad)) save_bc_grad_local = save_bc_grad

    !WORK Z-PENCILS
    call interzpv(ppi3,pp3,dip3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
         (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
    call derzpv(pgz3,pp3,dip3,sz,cfip6z,csip6z,cwip6z,cfz6,csz6,cwz6,&
         (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)

    !WORK Y-PENCILS
    call transpose_z_to_y(pgz3,pgz2,ph3) !nxm nym nz
    call transpose_z_to_y(ppi3,pp2,ph3)

    call interypv(ppi2,pp2,dip2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
         (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
    call derypv(pgy2,pp2,dip2,sy,cfip6y,csip6y,cwip6y,cfy6,csy6,cwy6,ppy,&
         (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
    call interypv(pgzi2,pgz2,dip2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
         (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)

    !WORK X-PENCILS

    call transpose_y_to_x(ppi2,pp1,ph2) !nxm ny nz
    call transpose_y_to_x(pgy2,pgy1,ph2)
    call transpose_y_to_x(pgzi2,pgz1,ph2)

    call derxpv(px1,pp1,di1,sx,cfip6,csip6,cwip6,cfx6,csx6,cwx6,&
         nxmsize,xsize(1),xsize(2),xsize(3),1)
    call interxpv(py1,pgy1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
         nxmsize,xsize(1),xsize(2),xsize(3),1)
    call interxpv(pz1,pgz1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
         nxmsize,xsize(1),xsize(2),xsize(3),1)

    if (.not. save_bc_grad_local) return

    if (iforces.eq.1) then
       call interxpv(ppi1,pp1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
            nxmsize,xsize(1),xsize(2),xsize(3),1)
    endif

    !we are in X pencils:
    if (nclx1.eq.2) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             dpdyx1(j,k)=py1(1,j,k)/gdt(itr)
             dpdzx1(j,k)=pz1(1,j,k)/gdt(itr)
          enddo
       enddo
    endif
    if (nclxn.eq.2) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             dpdyxn(j,k)=py1(nx,j,k)/gdt(itr)
             dpdzxn(j,k)=pz1(nx,j,k)/gdt(itr)
          enddo
       enddo
    endif

    if (ncly1.eq.2) then
       if (xstart(2)==1) then
          do k=1,xsize(3)
             do i=1,xsize(1)
                dpdxy1(i,k)=px1(i,1,k)/gdt(itr)
                dpdzy1(i,k)=pz1(i,1,k)/gdt(itr)
             enddo
          enddo
       endif
    endif
    if (nclyn.eq.2) then
       if (xend(2)==ny) then
          do k=1,xsize(3)
             do i=1,xsize(1)
                dpdxyn(i,k)=px1(i,xsize(2),k)/gdt(itr)
                dpdzyn(i,k)=pz1(i,xsize(2),k)/gdt(itr)
             enddo
          enddo
       endif
    endif

    if (nclz1.eq.2) then
       if (xstart(3)==1) then
          do j=1,xsize(2)
             do i=1,xsize(1)
                dpdxz1(i,j)=py1(i,j,1)/gdt(itr)
                dpdyz1(i,j)=pz1(i,j,1)/gdt(itr)
             enddo
          enddo
       endif
    endif
    if (nclzn.eq.2) then
       if (xend(3)==nz) then
          do j=1,xsize(2)
             do i=1,xsize(1)
                dpdxzn(i,j)=py1(i,j,xsize(3))/gdt(itr)
                dpdyzn(i,j)=pz1(i,j,xsize(3))/gdt(itr)
             enddo
          enddo
       endif
    endif

    return
  end subroutine gradp
  !############################################################################
  !############################################################################
  subroutine pre_correc(ux,uy,uz,ep)

    USE variables
    USE param
    USE var
    USE MPI
    use ibm, only : corgp_ibm, body

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,ep
    integer :: i,j,k,is
    real(mytype) :: ut,ut1

    integer :: code
    integer, dimension(2) :: dims, dummy_coords
    logical, dimension(2) :: dummy_periods
#ifdef DEBG
    real(mytype) dep
#endif

    call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, dummy_periods, dummy_coords, code)

    !********NCLX==0*************************************
    !we are in X pencils:
    if ((itype.eq.itype_pipe).and.(nclx1==0.and.nclxn==0)) then

        !Correct bulk velocity and temperature
        call pipe_bulk(ux,uy,uz,ep)

    endif

    !********NCLX==2*************************************
    !we are in X pencils:
    if ((itype.eq.itype_channel.or.itype.eq.itype_uniform.or.itype.eq.itype_abl).and.(nclx1==2.and.nclxn==2)) then

       !Computation of the flow rate Inflow/Outflow
       ut1=zero
       do k=1,xsize(3)
          do j=1,xsize(2)
             ut1=ut1+bxx1(j,k)
          enddo
       enddo
       call MPI_ALLREDUCE(MPI_IN_PLACE,ut1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       ut1=ut1/(real(ny*nz,mytype))
       ut=zero
       do k=1,xsize(3)
          do j=1,xsize(2)
             ut=ut+bxxn(j,k)
          enddo
       enddo
       call MPI_ALLREDUCE(MPI_IN_PLACE,ut,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       ut=ut/(real(ny*nz,mytype))
       if ((nrank==0).and.(mod(itime,ilist)==0)) &
          write(*,*) 'Flow rate x I/O/O-I',real(ut1,4),real(ut,4),real(ut-ut1,4)
       do k=1,xsize(3)
          do j=1,xsize(2)
             bxxn(j,k)=bxxn(j,k)-ut+ut1
          enddo
       enddo

    endif

    if (itype.eq.itype_tbl) call tbl_flrt(ux,uy,uz)

    if (nclx1==2) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             dpdyx1(j,k)=dpdyx1(j,k)*gdt(itr)
             dpdzx1(j,k)=dpdzx1(j,k)*gdt(itr)
          enddo
       enddo
       do k=1,xsize(3)
          do j=1,xsize(2)
             ux(1 ,j,k)=bxx1(j,k)
             uy(1 ,j,k)=bxy1(j,k)+dpdyx1(j,k)
             uz(1 ,j,k)=bxz1(j,k)+dpdzx1(j,k)
          enddo
       enddo
    endif
    if (nclxn==2) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             dpdyxn(j,k)=dpdyxn(j,k)*gdt(itr)
             dpdzxn(j,k)=dpdzxn(j,k)*gdt(itr)
          enddo
       enddo
       do k=1,xsize(3)
          do j=1,xsize(2)
             ux(nx,j,k)=bxxn(j,k)
             uy(nx,j,k)=bxyn(j,k)+dpdyxn(j,k)
             uz(nx,j,k)=bxzn(j,k)+dpdzxn(j,k)
          enddo
       enddo
    endif


    !********NCLX==1*************************************
    if (nclx1==1) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             ux(1 ,j,k)=zero
          enddo
       enddo
    endif
    if (nclxn==1) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             ux(nx,j,k)=zero
          enddo
       enddo
    endif

    !********NCLY==2*************************************
    if (ncly1==2) then
       if (xstart(2)==1) then
          do k=1,xsize(3)
             do i=1,xsize(1)
                dpdxy1(i,k)=dpdxy1(i,k)*gdt(itr)
                dpdzy1(i,k)=dpdzy1(i,k)*gdt(itr)
             enddo
          enddo
          if (mhd_active) then
             do k=1,xsize(3)
                do i=1,xsize(1)
                   ux(i,1,k)=byx1(i,k)
                   uy(i,1,k)=byy1(i,k)
                   uz(i,1,k)=byz1(i,k)
                enddo
             enddo
          else
             do k=1,xsize(3)
                do i=1,xsize(1)
                   ux(i,1,k)=byx1(i,k) + dpdxy1(i,k)
                   uy(i,1,k)=byy1(i,k)
                   uz(i,1,k)=byz1(i,k) + dpdzy1(i,k)
                enddo
             enddo
          endif
       endif
    endif

    if (nclyn==2) then
       if (xend(2)==ny) then
          do k=1,xsize(3)
             do i=1,xsize(1)
                dpdxyn(i,k)=dpdxyn(i,k)*gdt(itr)
                dpdzyn(i,k)=dpdzyn(i,k)*gdt(itr)
             enddo
          enddo
       endif
       if (dims(1)==1) then
          if (mhd_active) then
             do k=1,xsize(3)
                do i=1,xsize(1)
                   ux(i,xsize(2),k)=byxn(i,k)
                   uy(i,xsize(2),k)=byyn(i,k)
                   uz(i,xsize(2),k)=byzn(i,k)
                enddo
             enddo
          else
             do k=1,xsize(3)
                do i=1,xsize(1)
                   ux(i,xsize(2),k)=byxn(i,k) + dpdxyn(i,k)
                   uy(i,xsize(2),k)=byyn(i,k)
                   uz(i,xsize(2),k)=byzn(i,k) + dpdzyn(i,k)
                enddo
             enddo
          endif
       elseif (ny - (nym / dims(1)) == xstart(2)) then
          if (mhd_active) then
             do k=1,xsize(3)
                do i=1,xsize(1)
                   ux(i,xsize(2),k)=byxn(i,k)
                   uy(i,xsize(2),k)=byyn(i,k)
                   uz(i,xsize(2),k)=byzn(i,k)
                enddo
             enddo
          else
             do k=1,xsize(3)
                do i=1,xsize(1)
                   ux(i,xsize(2),k)=byxn(i,k) + dpdxyn(i,k)
                   uy(i,xsize(2),k)=byyn(i,k)
                   uz(i,xsize(2),k)=byzn(i,k) + dpdzyn(i,k)
                enddo
             enddo
          endif
       endif
    endif

    !********NCLY==1*************************************
    if (ncly1==1) then
       if (xstart(2)==1) then
          do k=1,xsize(3)
             do i=1,xsize(1)
                uy(i,1,k)=zero
             enddo
          enddo
       endif
    endif

    if (nclyn==1) then
       if (xend(2)==ny) then
          do k=1,xsize(3)
             do i=1,xsize(1)
                uy(i,xsize(2),k)=zero
             enddo
          enddo
       endif
    endif


    !********NCLZ==2*************************************
    if (nclz1==2) then
       if (xstart(3)==1) then
          do j=1,xsize(2)
             do i=1,xsize(1)
                dpdxz1(i,j)=dpdxz1(i,j)*gdt(itr)
                dpdyz1(i,j)=dpdyz1(i,j)*gdt(itr)
             enddo
          enddo
          do j=1,xsize(2)
             do i=1,xsize(1)
                ux(i,j,1)=bzx1(i,j)+dpdxz1(i,j)
                uy(i,j,1)=bzy1(i,j)+dpdyz1(i,j)
                uz(i,j,1)=bzz1(i,j)
             enddo
          enddo
       endif
    endif

    if (nclzn==2) then
       if (xend(3)==nz) then
          do j=1,xsize(2)
             do i=1,xsize(1)
                dpdxzn(i,j)=dpdxzn(i,j)*gdt(itr)
                dpdyzn(i,j)=dpdyzn(i,j)*gdt(itr)
             enddo
          enddo
          do j=1,xsize(2)
             do i=1,xsize(1)
                ux(i,j,xsize(3))=bzxn(i,j)+dpdxzn(i,j)
                uy(i,j,xsize(3))=bzyn(i,j)+dpdyzn(i,j)
                uz(i,j,xsize(3))=bzzn(i,j)
             enddo
          enddo
       endif
    endif
    !********NCLZ==1************************************* !just to reforce free-slip condition
    if (nclz1==1) then
       if (xstart(3)==1) then
          do j=1,xsize(2)
             do i=1,xsize(1)
                uz(i,j,1)=zero
             enddo
          enddo
       endif
    endif

    if (nclzn==1) then
       if (xend(3)==nz) then
          do j=1,xsize(2)
             do i=1,xsize(1)
                uz(i,j,xsize(3))=zero
             enddo
          enddo
       endif
    endif
#ifdef DEBG
    dep=maxval(abs(ux))
    call MPI_ALLREDUCE(MPI_IN_PLACE,dep,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## Pres corr ux ', dep
    dep=maxval(abs(uy))
    call MPI_ALLREDUCE(MPI_IN_PLACE,dep,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## Pres corr uy ', dep
    dep=maxval(abs(uz))
    call MPI_ALLREDUCE(MPI_IN_PLACE,dep,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## Pres corr uz ', dep
#endif

    if (iibm==1) then !solid body old school
       call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,1)
       call body(ux1,uy1,uz1,ep1)
       call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,2)
    endif

    return
  end subroutine pre_correc
  !############################################################################
  !############################################################################
  !! Convert to/from conserved/primary variables
  SUBROUTINE primary_to_conserved(rho1, var1)

    USE param, ONLY : nrhotime

    IMPLICIT NONE

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: var1

    var1(:,:,:) = rho1(:,:,:,1) * var1(:,:,:)

  ENDSUBROUTINE primary_to_conserved
  !############################################################################
  !############################################################################
  SUBROUTINE velocity_to_momentum (rho1, ux1, uy1, uz1)

    USE param, ONLY : nrhotime
    USE var, ONLY : ilmn

    IMPLICIT NONE

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1

    IF (.NOT.ilmn) THEN
       RETURN
    ENDIF

    CALL primary_to_conserved(rho1, ux1)
    CALL primary_to_conserved(rho1, uy1)
    CALL primary_to_conserved(rho1, uz1)

  ENDSUBROUTINE velocity_to_momentum
  !############################################################################
  !############################################################################
  SUBROUTINE conserved_to_primary(rho1, var1)

    USE param, ONLY : nrhotime

    IMPLICIT NONE

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: var1

    var1(:,:,:) = var1(:,:,:) / rho1(:,:,:,1)

  ENDSUBROUTINE conserved_to_primary
  !############################################################################
  !############################################################################
  SUBROUTINE momentum_to_velocity (rho1, ux1, uy1, uz1)

    USE param, ONLY : nrhotime
    USE var, ONLY : ilmn

    IMPLICIT NONE

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1

    IF (.NOT.ilmn) THEN
       RETURN
    ENDIF

    CALL conserved_to_primary(rho1, ux1)
    CALL conserved_to_primary(rho1, uy1)
    CALL conserved_to_primary(rho1, uz1)

  ENDSUBROUTINE momentum_to_velocity
  !############################################################################
  !############################################################################
  !! Calculate velocity-divergence constraint
  SUBROUTINE calc_divu_constraint(divu3, rho1, phi1)

    USE param, ONLY : nrhotime, zero, ilmn, pressure0, imultispecies, massfrac, mol_weight
    USE param, ONLY : ibirman_eos
    USE param, ONLY : xnu, prandtl
    USE param, ONLY : one
    USE param, ONLY : iimplicit
    USE variables

    USE var, ONLY : ta1, tb1, tc1, td1, di1
    USE var, ONLY : phi2, ta2, tb2, tc2, td2, te2, di2
    USE var, ONLY : phi3, ta3, tb3, tc3, td3, rho3, di3
    USE param, only : zero
    IMPLICIT NONE

    INTEGER :: is, tmp

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), numscalar) :: phi1
    REAL(mytype), INTENT(OUT), DIMENSION(zsize(1), zsize(2), zsize(3)) :: divu3

    IF (ilmn.and.(.not.ibirman_eos)) THEN
       !!------------------------------------------------------------------------------
       !! X-pencil

       !! We need temperature
       CALL calc_temp_eos(ta1, rho1(:,:,:,1), phi1, tb1, xsize(1), xsize(2), xsize(3))

       CALL derxx (tb1, ta1, di1, sx, sfxp, ssxp, swxp, xsize(1), xsize(2), xsize(3), 1, 0)
       IF (imultispecies) THEN
          tb1(:,:,:) = (xnu / prandtl) * tb1(:,:,:) / ta1(:,:,:)

          !! Calc mean molecular weight
          td1(:,:,:) = zero
          DO is = 1, numscalar
             IF (massfrac(is)) THEN
                td1(:,:,:) = td1(:,:,:) + phi1(:,:,:,is) / mol_weight(is)
             ENDIF
          ENDDO
          td1(:,:,:) = one / td1(:,:,:)

          DO is = 1, numscalar
             IF (massfrac(is)) THEN
                CALL derxx (tc1, phi1(:,:,:,is), di1, sx, sfxp, ssxp, swxp, xsize(1), xsize(2), xsize(3), 1, 0)
                tb1(:,:,:) = tb1(:,:,:) + (xnu / sc(is)) * (td1(:,:,:) / mol_weight(is)) * tc1(:,:,:)
             ENDIF
          ENDDO
       ENDIF

       CALL transpose_x_to_y(ta1, ta2)        !! Temperature
       CALL transpose_x_to_y(tb1, tb2)        !! d2Tdx2
       IF (imultispecies) THEN
          DO is = 1, numscalar
             IF (massfrac(is)) THEN
                CALL transpose_x_to_y(phi1(:,:,:,is), phi2(:,:,:,is))
             ENDIF
          ENDDO
       ENDIF

       !!------------------------------------------------------------------------------
       !! Y-pencil
       tmp = iimplicit
       iimplicit = 0
       CALL deryy (tc2, ta2, di2, sy, sfyp, ssyp, swyp, ysize(1), ysize(2), ysize(3), 1, 0)
       iimplicit = tmp
       IF (imultispecies) THEN
          tc2(:,:,:) = (xnu / prandtl) * tc2(:,:,:) / ta2(:,:,:)

          !! Calc mean molecular weight
          te2(:,:,:) = zero
          DO is = 1, numscalar
             IF (massfrac(is)) THEN
                te2(:,:,:) = te2(:,:,:) + phi2(:,:,:,is) / mol_weight(is)
             ENDIF
          ENDDO
          te2(:,:,:) = one / te2(:,:,:)

          DO is = 1, numscalar
             IF (massfrac(is)) THEN
                tmp = iimplicit
                iimplicit = 0
                CALL deryy (td2, phi2(:,:,:,is), di2, sy, sfyp, ssyp, swyp, ysize(1), ysize(2), ysize(3), 1, 0)
                iimplicit = tmp
                tc2(:,:,:) = tc2(:,:,:) + (xnu / sc(is)) * (te2(:,:,:) / mol_weight(is)) * td2(:,:,:)
             ENDIF
          ENDDO
       ENDIF
       tb2(:,:,:) = tb2(:,:,:) + tc2(:,:,:)

       CALL transpose_y_to_z(ta2, ta3)        !! Temperature
       CALL transpose_y_to_z(tb2, tb3)        !! d2Tdx2 + d2Tdy2
       IF (imultispecies) THEN
          DO is = 1, numscalar
             IF (massfrac(is)) THEN
                CALL transpose_y_to_z(phi2(:,:,:,is), phi3(:,:,:,is))
             ENDIF
          ENDDO
       ENDIF

       !!------------------------------------------------------------------------------
       !! Z-pencil
       CALL derzz (divu3, ta3, di3, sz, sfzp, sszp, swzp, zsize(1), zsize(2), zsize(3), 1, 0)
       IF (imultispecies) THEN
          divu3(:,:,:) = (xnu / prandtl) * divu3(:,:,:) / ta3(:,:,:)

          !! Calc mean molecular weight
          td3(:,:,:) = zero
          DO is = 1, numscalar
             IF (massfrac(is)) THEN
                td3(:,:,:) = td3(:,:,:) + phi3(:,:,:,is) / mol_weight(is)
             ENDIF
          ENDDO
          td3(:,:,:) = one / td3(:,:,:)

          DO is = 1, numscalar
             IF (massfrac(is)) THEN
                CALL derzz (tc3, phi3(:,:,:,is), di3, sz, sfzp, sszp, swzp, zsize(1), zsize(2), zsize(3), 1, 0)
                divu3(:,:,:) = divu3(:,:,:) + (xnu / sc(is)) * (td3(:,:,:) / mol_weight(is)) * tc3(:,:,:)
             ENDIF
          ENDDO
       ENDIF
       divu3(:,:,:) = divu3(:,:,:) + tb3(:,:,:)

       IF (imultispecies) THEN
          !! Thus far we have computed rho * divu, want divu
          CALL calc_rho_eos(rho3, ta3, phi3, tb3, zsize(1), zsize(2), zsize(3))
          divu3(:,:,:) = divu3(:,:,:) / rho3(:,:,:)
       ELSE
          divu3(:,:,:) = (xnu / prandtl) * divu3(:,:,:) / pressure0
       ENDIF
    ELSE
       divu3(:,:,:) = zero
    ENDIF

  ENDSUBROUTINE calc_divu_constraint

  !############################################################################
  !############################################################################
  ! Calculate extrapolation drhodt 
  SUBROUTINE extrapol_drhodt(drhodt1_next, rho1, drho1)

    USE param, ONLY : ntime, nrhotime, itime, itimescheme, itr, dt, gdt, irestart
    USE param, ONLY : half, three, four
    USE param, ONLY : ibirman_eos

    IMPLICIT NONE

    INTEGER :: subitr

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), ntime) :: drho1
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), INTENT(OUT), DIMENSION(xsize(1), xsize(2), xsize(3)) :: drhodt1_next

    IF (itimescheme.EQ.1) THEN
       !! EULER
       drhodt1_next(:,:,:) = (rho1(:,:,:,1) - rho1(:,:,:,2)) / dt
    ELSEIF (itimescheme.EQ.2) THEN
       !! AB2
       IF ((itime.EQ.1).AND.(irestart.EQ.0)) THEN
          drhodt1_next(:,:,:) = (rho1(:,:,:,1) - rho1(:,:,:,2)) / dt
       ELSE
          drhodt1_next(:,:,:) = three * rho1(:,:,:,1) - four * rho1(:,:,:,2) + rho1(:,:,:,3)
          drhodt1_next(:,:,:) = half * drhodt1_next(:,:,:) / dt
       ENDIF
       ! ELSEIF (itimescheme.EQ.3) THEN
       !    !! AB3
       ! ELSEIF (itimescheme.EQ.4) THEN
       !    !! AB4
    ELSEIF (itimescheme.EQ.5) THEN
       !! RK3
       IF (itime.GT.1) THEN
          drhodt1_next(:,:,:) = rho1(:,:,:,2)
          DO subitr = 1, itr
             drhodt1_next(:,:,:) = drhodt1_next(:,:,:) + (gdt(subitr) / dt) &
                  * (rho1(:,:,:,2) - rho1(:,:,:,3))
          ENDDO
       ELSE
          drhodt1_next(:,:,:) = drho1(:,:,:,1)
       ENDIF
    ELSE
       IF (nrank.EQ.0) THEN
          PRINT *, "Extrapolating drhodt not implemented for timescheme:", itimescheme
          STOP
       ENDIF
    ENDIF

    IF (ibirman_eos) THEN
       CALL birman_drhodt_corr(drhodt1_next, rho1)
    ENDIF

  ENDSUBROUTINE extrapol_drhodt
  !############################################################################
  !!
  !! Subroutine : birman_drhodt_corr
  !! Author     :
  !! Description: Calculate extrapolation drhodt correction
  !!
  !############################################################################

  SUBROUTINE birman_drhodt_corr(drhodt1_next, rho1)

    USE variables, ONLY : derxx, deryy, derzz
    USE param, ONLY : nrhotime
    USE param, ONLY : xnu, prandtl
    USE param, ONLY : iimplicit

    USE var, ONLY : td1, te1, di1, sx, sfxp, ssxp, swxp
    USE var, ONLY : rho2, ta2, tb2, di2, sy, sfyp, ssyp, swyp
    USE var, ONLY : rho3, ta3, di3, sz, sfzp, sszp, swzp
    USE param, only : zero
    IMPLICIT NONE

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: drhodt1_next

    REAL(mytype) :: invpe

    invpe = xnu / prandtl

    CALL transpose_x_to_y(rho1(:,:,:,1), rho2)
    CALL transpose_y_to_z(rho2, rho3)

    !! Diffusion term
    CALL derzz (ta3,rho3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1, 0)
    CALL transpose_z_to_y(ta3, tb2)

    iimplicit = -iimplicit
    CALL deryy (ta2,rho2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1, 0)
    iimplicit = -iimplicit
    ta2(:,:,:) = ta2(:,:,:) + tb2(:,:,:)
    CALL transpose_y_to_x(ta2, te1)

    CALL derxx (td1,rho1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1, 0)
    td1(:,:,:) = td1(:,:,:) + te1(:,:,:)

    drhodt1_next(:,:,:) = drhodt1_next(:,:,:) - invpe * td1(:,:,:)

  ENDSUBROUTINE birman_drhodt_corr
  !############################################################################
  !!
  !!  SUBROUTINE: test_varcoeff
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Tests convergence of the variable-coefficient Poisson solver
  !!
  !############################################################################
  SUBROUTINE test_varcoeff(converged, divup3norm, pp3, dv3, atol, rtol, poissiter)

    USE MPI
    USE var, ONLY : nzmsize
    USE param, ONLY : npress, itime
    USE variables, ONLY : nxm, nym, nzm, ilist

    IMPLICIT NONE

    !! INPUTS
    REAL(mytype), INTENT(INOUT), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3
    REAL(mytype), INTENT(IN), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize) :: dv3
    REAL(mytype), INTENT(IN) :: atol, rtol
    INTEGER, INTENT(IN) :: poissiter

    !! OUTPUTS
    LOGICAL, INTENT(OUT) :: converged
    REAL(mytype) :: divup3norm

    !! LOCALS
    INTEGER :: ierr
    REAL(mytype) :: errglob

    IF (poissiter.EQ.0) THEN
       divup3norm = SUM(dv3**2)
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,divup3norm,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
       divup3norm = SQRT(divup3norm / nxm / nym / nzm)

       if (nrank.eq.0.and.mod(itime, ilist) == 0) then
          write(*,*)  "solving variable-coefficient poisson equation:"
          write(*,*)  "+ rms div(u*) - div(u): ", divup3norm
       endif
    else
       !! Compute RMS change
       errglob = SUM((pp3(:,:,:,1) - pp3(:,:,:,2))**2)
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,errglob,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
       errglob = SQRT(errglob / nxm / nym / nzm)

       if (nrank.eq.0.and.mod(itime, ilist) == 0) then
          write(*,*)  "+ RMS change in pressure: ", errglob
       endif

       if (errglob.le.atol) then
          converged = .true.
          if (nrank.eq.0.and.mod(itime, ilist) == 0) then
             write(*,*)  "- Converged: atol"
          endif
       endif

       !! Compare RMS change to size of |div(u*) - div(u)|
       if (errglob.lt.(rtol * divup3norm)) then
          converged = .true.
          if (nrank.eq.0.and.mod(itime, ilist) == 0) then
             write(*,*)  "- Converged: rtol"
          endif
       endif

       if (.not.converged) then
          pp3(:,:,:,2) = pp3(:,:,:,1)
       endif
    endif

  ENDSUBROUTINE test_varcoeff
  !############################################################################
  !!
  !! SUBROUTINE: calc_rho0
  !! DESCRIPTION: Computes the density used to normalise the variable-coefficient Poisson equation.
  subroutine calc_rho0(rho1, rho0)

    use mpi

    use param, only: ilmn, ivarcoeff

    implicit none

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: rho1 ! The density field
    real(mytype), intent(out) :: rho0                                         ! The normalising density (globally constant)

    integer :: ierr

    if (ilmn .and. ivarcoeff) then
       !! Compute rho0
       rho0 = minval(rho1(:,:,:))

       call mpi_allreduce(MPI_IN_PLACE,rho0,1,real_type,mpi_min,mpi_comm_world,ierr)
    else
       ! Ensure a default value
       rho0 = 1.0_mytype
    end if

  end subroutine calc_rho0
  !############################################################################
  !############################################################################
  !!
  !!  SUBROUTINE: calc_varcoeff_rhs
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Computes RHS of the variable-coefficient Poisson solver
  !!
  !############################################################################
  SUBROUTINE calc_varcoeff_rhs(pp3, rho1, px1, py1, pz1, dv3, drho1, ep1, divu3, rho0)

    USE MPI

    USE param, ONLY : nrhotime, ntime
    USE param, ONLY : one

    USE var, ONLY : ta1, tb1, tc1
    USE var, ONLY : nzmsize

    IMPLICIT NONE

    !! INPUTS
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3)) :: px1, py1, pz1
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), ntime) :: drho1
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ep1
    REAL(mytype), INTENT(IN), DIMENSION(zsize(1), zsize(2), zsize(3)) :: divu3
    REAL(mytype), INTENT(IN), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize) :: dv3
    real(mytype), intent(in) :: rho0

    !! OUTPUTS
    REAL(mytype), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize) :: pp3,div_visu_var

    !! LOCALS
    INTEGER :: nlock, ierr

    ta1(:,:,:) = (one - rho0 / rho1(:,:,:,1)) * px1(:,:,:)
    tb1(:,:,:) = (one - rho0 / rho1(:,:,:,1)) * py1(:,:,:)
    tc1(:,:,:) = (one - rho0 / rho1(:,:,:,1)) * pz1(:,:,:)

    nlock = -1 !! Don't do any funny business with LMN
    call divergence(div_visu_var,pp3,rho1,ta1,tb1,tc1,ep1,drho1,divu3,nlock)

    !! lapl(p) = div((1 - rhomin/rho) grad(p)) + rhomin(div(u*) - div(u))
    !! dv3 contains div(u*) - div(u)
    pp3(:,:,:) = pp3(:,:,:) + rho0 * dv3(:,:,:)

  ENDSUBROUTINE calc_varcoeff_rhs
  !############################################################################
  !********************************************************************
  !
  subroutine tbl_flrt (ux1,uy1,uz1)
  !
  !********************************************************************

    USE decomp_2d_poisson
    USE variables
    USE param
    USE MPI

    implicit none
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2

    integer :: j,i,k,code
    real(mytype) :: can,ut1,ut2,ut3,ut4,udif

    ux1(1,:,:)=bxx1(:,:)
    ux1(nx,:,:)=bxxn(:,:)

    call transpose_x_to_y(ux1,ux2)
    call transpose_x_to_y(uy1,uy2)
    ! Flow rate at the inlet
    ut1=zero
    if (ystart(1)==1) then !! CPUs at the inletdiv_vi
      do k=1,ysize(3)
        do j=1,ysize(2)-1
          ut1=ut1+(yp(j+1)-yp(j))*(ux2(1,j+1,k)-half*(ux2(1,j+1,k)-ux2(1,j,k)))
        enddo
      enddo
      ! ut1=ut1/real(ysize(3),mytype)
    endif
    call MPI_ALLREDUCE(MPI_IN_PLACE,ut1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    ut1=ut1/real(nz,mytype) !! Volume flow rate per unit spanwise dist
    ! Flow rate at the outlet
    ut2=zero
    if (yend(1)==nx) then !! CPUs at the outlet
      do k=1,ysize(3)
        do j=1,ysize(2)-1
          ut2=ut2+(yp(j+1)-yp(j))*(ux2(ysize(1),j+1,k)-half*(ux2(ysize(1),j+1,k)-ux2(ysize(1),j,k)))
        enddo
      enddo
      ! ut2=ut2/real(ysize(3),mytype)
    endif

    call MPI_ALLREDUCE(MPI_IN_PLACE,ut2,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    ut2=ut2/real(nz,mytype) !! Volume flow rate per unit spanwise dist

    ! Flow rate at the top and bottom
    ut3=zero
    ut4=zero
    do k=1,ysize(3)
      do i=1,ysize(1)
        ut3=ut3+uy2(i,1,k)
        ut4=ut4+uy2(i,ny,k)
      enddo
    enddo
    call MPI_ALLREDUCE(MPI_IN_PLACE,ut3,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(MPI_IN_PLACE,ut4,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    ut3=ut3/(real(nx*nz,mytype))*xlx  !!! Volume flow rate per unit spanwise dist
    ut4=ut4/(real(nx*nz,mytype))*xlx  !!! Volume flow rate per unit spanwise dist

    !! velocity correction
    udif=(ut1-ut2+ut3-ut4)/yly
    if ((nrank==0).and.(mod(itime,ilist)==0)) then
      write(*,"(' Mass balance: L-BC, R-BC,',2f12.6)") ut1,ut2
      write(*,"(' Mass balance: B-BC, T-BC, Crr-Vel',3f11.5)") ut3,ut4,udif
    endif
    ! do k=1,xsize(3)
    !   do j=1,xsize(2)
    !     ux1(nx,i,k)=ux1(nx,i,k)+udif
    !   enddo
    ! enddo
    do k=1,xsize(3)
      do j=1,xsize(2)
        bxxn(j,k)=bxxn(j,k)+udif
      enddo
    enddo

  end subroutine tbl_flrt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  subroutine: pipe_bulk / pipe_bulk_u / pipe_bulk_phi
  !!      AUTHOR: Rodrigo Vicente Cruz
  !! DESCRIPTION: Correction of pipe's bulk velocity (constant
  !!              flow rate) and bulk temperature.
  !!              See Thesis Vicente Cruz 2021 for help.
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !********************************************************************
  !
  subroutine pipe_bulk(ux,uy,uz,ep)
  !
  !********************************************************************

    use param, only: one
    use variables, only: numscalar
    use var, only: ta1, phi1

    implicit none
    real(mytype),dimension(xsize(1),xsize(2),xsize(3))  :: ux,uy,uz,ep
    real(mytype),save                                   :: ncount = -one
    integer                                             :: is

    !Compute the number of cells inside the pipe at the beginning
    if (ncount < 0) then
       ta1(:,:,:) = one
       call pipe_volume_avg(ta1,ncount,ep,one)
    endif

    !Bulk velocity correction
    call pipe_bulk_u(ux,uy,uz,ep,one,ncount)

    !Bulk temperature correction
    if (numscalar.ne.0) then
        do is=1,numscalar
            call pipe_bulk_phi(phi1(:,:,:,is),ux,ep,is,one,ncount)
        enddo
    endif

  end subroutine pipe_bulk
  !********************************************************************
  !
  subroutine pipe_bulk_u(ux,uy,uz,ep,ub_constant,ncount)
  !
  !********************************************************************

    use variables
    use param
    use var
    use ibm_param, only: rai
    use MPI

    implicit none
    !INPUTS
    real(mytype),intent(inout),dimension(xsize(1),xsize(2),xsize(3))    :: ux,uy,uz
    real(mytype),intent(in   ),dimension(xsize(1),xsize(2),xsize(3))    :: ep
    real(mytype),intent(in   )                                          :: ub_constant !bulk velocity value
    real(mytype),intent(in   )                                          :: ncount !numer of cells inside the pipe
    !LOCALS
    real(mytype)                                                        :: qm      !flow rate
    real(mytype)                                                        :: ym,zm,yc,zc,r
    integer                                                             :: j,i,k
    integer, save                                                       :: local_io_unit=-1

    yc = yly/two
    zc = zlz/two

    if (local_io_unit.eq.-1 .and. nrank.eq.0) then
       open(newunit=local_io_unit,file='Ub.dat',status='unknown')
    endif

    !--------------------------- Bulk Velocity ---------------------------
    !Calculate loss of streamwise mean pressure gradient
    call pipe_volume_avg(ux,qm,ep,ncount)
    if (nrank==0) then
       if (mod(itime, ilist)==0) print *,'Velocity:'
       if (mod(itime, ilist)==0) print *,'    Bulk velocity before',qm
       write(local_io_unit,*) real((itime-1)*dt,mytype), (ub_constant-qm) !write pressure drop
    endif

    !Correction
    do k=1,xsize(3)
        zm=dz*real(xstart(3)-1+k-1,mytype)-zc
        do j=1,xsize(2)
            if (istret.eq.0) ym=real(j+xstart(2)-1-1,mytype)*dy-yc
            if (istret.ne.0) ym=yp(j+xstart(2)-1)-yc
            r=sqrt(ym*ym+zm*zm)
            do i=1,xsize(1)
                if (r.le.rai.and.ep(i,j,k).eq.0) then
                    ux(i,j,k)=ux(i,j,k)+(ub_constant-qm)
                else
                    !Cancel solid zone (rai <= r <= rao)
                    !and buffer zone (r > rao)
                    ux(i,j,k)=zero
                    uy(i,j,k)=zero
                    uz(i,j,k)=zero
                endif
            enddo
        enddo
    enddo

    !Check new bulk velocity
    if (mod(itime, ilist)==0) then
        call pipe_volume_avg(ux,qm,ep,ncount)
        if (nrank==0) print *,'    Bulk velocity  after',qm
    endif
    !
    return
    !
  end subroutine pipe_bulk_u
  !********************************************************************
  !
  subroutine pipe_bulk_phi(phi,ux,ep,is,phib_constant,ncount)
  !
  !********************************************************************

    use decomp_2d_poisson
    use variables
    use param
    use ibm_param, only: rai
    use MPI

    implicit none
    !INPUTS
    real(mytype),intent(inout),dimension(xsize(1),xsize(2),xsize(3))  :: phi,ux
    real(mytype),intent(in   ),dimension(xsize(1),xsize(2),xsize(3))  :: ep
    real(mytype),intent(in   )                                        :: phib_constant !bulk temperature value
    real(mytype),intent(in   )                                        :: ncount !numer of cells inside the pipe
    !LOCALS
    real(mytype)                                                      :: qv,qm !volumetric averaged values
    real(mytype)                                                      :: ym,zm,yc,zc,r
    real(mytype)                                                      :: phi_out
    integer                                                           :: is,j,i,k
    character(len=30)                                                 :: filename
    integer, allocatable, save, dimension(:)                          :: local_io_unit

    if (iscalar.eq.0) return

255 format('Tb_',I2.2,'.dat')
256 format(' Scalar:                       #',I2)

    ! Safety check
    if (is<1.or.is>numscalar) return

    if (.not.allocated(local_io_unit)) then
        allocate(local_io_unit(numscalar))
        local_io_unit(:)=-1
    endif

    if (local_io_unit(is).eq.-1.and.nrank.eq.0) then
        write(filename,255) is
        open(newunit=local_io_unit(is),file=filename,status='unknown')
    endif

    yc = yly / two
    zc = zlz / two
    phi_out =zero !for reconstruction smoothness

    !--------------------------- Bulk Temperature ---------------------------
    !                  with corrected streamwise velocity
    call pipe_volume_avg(ux*phi,qm,ep,ncount)
    call pipe_volume_avg(ux*ux ,qv,ep,ncount)
    if (nrank.eq.0) then
        if (mod(itime, ilist)==0) write(*,256) is
        if (mod(itime, ilist)==0) print *,'         Bulk phi before',qm
        write(local_io_unit(is),*) real((itime-1)*dt,mytype),(phib_constant-qm)/(qv)
    endif

    !Correction
    do k=1,xsize(3)
        zm=dz*real(xstart(3)-1+k-1,mytype)-zc
        do j=1,xsize(2)
            if (istret.eq.0) ym=real(j+xstart(2)-1-1,mytype)*dy-yc
            if (istret.ne.0) ym=yp(j+xstart(2)-1)-yc
            r=sqrt(ym*ym+zm*zm)
            do i=1,xsize(1)
                if (r.le.rai.and.ep(i,j,k).eq.0) then
                    phi(i,j,k)=phi(i,j,k)+ux(i,j,k)*((phib_constant-qm)/qv)
                else !smoothness for reconstruction
                    phi(i,j,k)=phi_out
                endif
            enddo
        enddo
    enddo

    !Check new bulk temperature
    if (mod(itime, ilist)==0) then
        call pipe_volume_avg(ux*phi,qm,ep,ncount)
        if (nrank==0) print *,'          Bulk phi after',qm
    endif
    !
    return

  end subroutine pipe_bulk_phi
  !********************************************************************
  !
  subroutine pipe_volume_avg(var,qm,ep,ncount)
  !
  !********************************************************************

    use param
    use variables
    use MPI
    use ibm_param, only: rai

    implicit none

    !INPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3))  :: var,ep
    real(mytype),intent(out)                            :: qm
    real(mytype),intent(in)                             :: ncount
    !LOCALS
    real(mytype)                                        :: ym,yc,zm,zc,r
    integer                                             :: i,j,k,code

    !Compute volumetric average of var
    !in the inner fluid zone r <= rai
    yc=yly/two
    zc=zlz/two
    qm=zero
    do k=1,xsize(3)
        zm=dz*real(xstart(3)-1+k-1,mytype)-zc
        do j=1,xsize(2)
            if (istret.eq.0) ym=real(j+xstart(2)-1-1,mytype)*dy-yc
            if (istret.ne.0) ym=yp(j+xstart(2)-1)-yc
            r=sqrt(ym*ym+zm*zm)
            do i=1,xsize(1)
                if (r.le.rai.and.ep(i,j,k).eq.0) then
                    qm=qm+var(i,j,k)
                endif
            enddo
        enddo
    enddo
    call MPI_ALLREDUCE(MPI_IN_PLACE,qm,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    qm=qm/ncount
    !
    return

  end subroutine pipe_volume_avg

endmodule navier
