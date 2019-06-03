#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "SpaceMask.h"
#include "GRHydro_Macros.h"


subroutine Conservative2PrimitiveHot(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  integer :: i, j, k, itracer, nx, ny, nz, myproc
  CCTK_REAL :: uxx, uxy, uxz, uyy, uyz, uzz, sdet, pmin, epsmin, dummy1, dummy2
  CCTK_REAL :: vlowx, vlowy, vlowz
  logical :: epsnegative
  character*512 :: warnline

  CCTK_REAL :: local_min_tracer
  CCTK_REAL :: local_perc_ptol
  integer :: reflevel

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup
  pointer (pvup,vup)

! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1), reset_to_atmo
  CCTK_REAL :: xpress(1),xeps(1),xtemp(1),xye(1),xrho(1)
  n = 1;keytemp = 0;anyerr = 0;keyerr(1) = 0
  xpress = 0.0d0;xeps = 0.0d0;xtemp = 0.0d0;xye = 0.0d0
! end EOS Omni vars

  myproc = CCTK_MyProc(cctkGH)
  reflevel = int(log10(dble(cctk_levfac(1)))/log10(2.0d0))

  ! save memory when MP is not used
  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
    pvup = loc(lvel)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
    pvup = loc(vel)
  end if

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  if (use_min_tracer .ne. 0) then
    local_min_tracer = min_tracer
  else
    local_min_tracer = 0.0d0
  end if

  ! this is a poly call
  xrho(1) = GRHydro_rho_min
  call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
       xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
  call EOS_Omni_EpsFromPress(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
       xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
  pmin = xpress(1)
  epsmin = xeps(1)

  !$OMP PARALLEL DO PRIVATE(i,j,k,itracer,&
  !$OMP uxx, uxy, uxz, uyy, uyz, uzz, sdet, epsnegative, anyerr, keyerr, keytemp,&
  !$OMP warnline, dummy1, dummy2,reset_to_atmo)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx


        !do not compute if in atmosphere
        if (atmosphere_mask(i,j,k) .gt. 0) cycle

        epsnegative = .false.

        sdet = sdetg(i,j,k)

        call UpperMetric(uxx,uxy,uxz,uyy,uyz,uzz,sdet*sdet,&
             g11(i,j,k),g12(i,j,k),g13(i,j,k),g22(i,j,k),&
             g23(i,j,k),g33(i,j,k))

        if (evolve_tracer .ne. 0) then
          do itracer=1,number_of_tracers
            call Con2Prim_ptTracer(cons_tracer(i,j,k,itracer), tracer(i,j,k,itracer), &
                   dens(i,j,k))

            if (use_min_tracer .ne. 0) then
              if (tracer(i,j,k,itracer) .le. local_min_tracer) then
                tracer(i,j,k,itracer) = local_min_tracer
              end if
            end if
          end do
        end if

        if(evolve_Y_e.ne.0) then
           Y_e(i,j,k) = max(min(Y_e_con(i,j,k) / dens(i,j,k),GRHydro_Y_e_max),&
                GRHydro_Y_e_min)
        endif

        reset_to_atmo = 0
        IF_BELOW_ATMO(dens(i,j,k), sdet*GRHydro_rho_min, GRHydro_atmo_tolerance, r(i,j,k)) then
           reset_to_atmo = 1
        endif

        if (reset_to_atmo .gt. 0 .or. (GRHydro_enable_internal_excision /= 0 .and. hydro_excision_mask(i,j,k) .gt. 0)) then
          SET_ATMO_MIN(dens(i,j,k), sdet*GRHydro_rho_min, r(i,j,k))
          SET_ATMO_MIN(rho(i,j,k), GRHydro_rho_min, r(i,j,k))
          scon(i,j,k,:) = 0.d0
          vup(i,j,k,:) = 0.d0
          w_lorentz(i,j,k) = 1.d0

          ! set hot atmosphere values
          temperature(i,j,k) = grhydro_hot_atmo_temp
          y_e(i,j,k) = grhydro_hot_atmo_Y_e
          y_e_con(i,j,k) = y_e(i,j,k) * dens(i,j,k)
          keytemp = 1
          call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
               rho(i,j,k),eps(i,j,k),temperature(i,j,k),y_e(i,j,k),&
               press(i,j,k),keyerr,anyerr)
          keytemp = 0
          ! w_lorentz=1, so the expression for tau reduces to:
          tau(i,j,k)  = sdet * (rho(i,j,k)+rho(i,j,k)*eps(i,j,k)) - dens(i,j,k)
          cycle
        end if
        
        call Con2Prim_pt_hot3(int(cctk_iteration,ik),myproc,int(i,ik),int(j,ik),int(k,ik),GRHydro_eos_handle,&
             dens(i,j,k),scon(i,j,k,1),&
             scon(i,j,k,2),scon(i,j,k,3),tau(i,j,k),Y_e_con(i,j,k),rho(i,j,k),vup(i,j,k,1),&
             vup(i,j,k,2),vup(i,j,k,3),eps(i,j,k),temperature(i,j,k),y_e(i,j,k),&
             press(i,j,k),w_lorentz(i,j,k), &
             uxx,uxy,uxz,uyy,uyz,uzz,sdet,x(i,j,k),y(i,j,k), &
             z(i,j,k),r(i,j,k),epsnegative,GRHydro_rho_min,pmin, epsmin, &
             int(reflevel,ik), GRHydro_C2P_failed(i,j,k), GRHydro_perc_ptol)

        if(temperature(i,j,k).gt.GRHydro_max_temp) then
          !$OMP CRITICAL
          call CCTK_WARN(1,"C2P: Temperature too high")
          write(warnline,"(i8)") cctk_iteration
          call CCTK_WARN(1,warnline)
          write(warnline,"(4i5,1P10E15.6)") GRHydro_Reflevel,i,j,k,x(i,j,k),y(i,j,k),z(i,j,k)
          call CCTK_WARN(1,warnline)
          write(warnline,"(1P10E15.6)") dens(i,j,k),scon(i,j,k,1:3),&
                tau(i,j,k),w_lorentz(i,j,k)
          call CCTK_WARN(1,warnline)
          write(warnline,"(1P10E15.6)") rho(i,j,k),eps(i,j,k),&
                temperature(i,j,k),Y_e(i,j,k)
          call CCTK_WARN(1,warnline)
          write(warnline,"(A7,i8)") "code: ",keyerr(1)
          call CCTK_WARN(1,warnline)
          write(warnline,"(A10,i5)") "reflevel: ", reflevel
          call CCTK_WARN(1,warnline)
          call CCTK_ERROR("Aborting!!!")
          STOP
          !$OMP END CRITICAL
        endif


        if(abs(GRHydro_C2P_failed(i,j,k)-2.0d0) .lt. 1.0d-10) then
          ! this means c2p did not converge.
          ! In this case, we attempt to call c2p with a reduced
          ! accuracy requirement; if it fails again, we abort
          GRHydro_C2P_failed(i,j,k) = 0
          if(temperature(i,j,k) .gt. GRHydro_hot_atmo_temp) then
            local_perc_ptol = GRHydro_perc_ptol*100.0d0
          else
            ! If we are in the extrapolation regime for the EOS,
            ! we accept a larger c2p error, since we will be resetting
            ! the temperature anyway. The error scale here is chosen
            ! somewhat arbitrarily to be 0.01% in pressnew/pressold.
            local_perc_ptol = 1.0d-4
          endif
          call Con2Prim_pt_hot3(int(cctk_iteration,ik),myproc,int(i,ik),int(j,ik),int(k,ik),GRHydro_eos_handle,&
                dens(i,j,k),scon(i,j,k,1),&
                scon(i,j,k,2),scon(i,j,k,3),tau(i,j,k),Y_e_con(i,j,k),rho(i,j,k),&
                vup(i,j,k,1),vup(i,j,k,2), vup(i,j,k,3),eps(i,j,k),&
                temperature(i,j,k),y_e(i,j,k),press(i,j,k),w_lorentz(i,j,k), &
                uxx,uxy,uxz,uyy,uyz,uzz,sdet,x(i,j,k),y(i,j,k), &
                z(i,j,k),r(i,j,k),epsnegative,GRHydro_rho_min,pmin, epsmin, &
                int(reflevel,ik), GRHydro_C2P_failed(i,j,k), local_perc_ptol)
          if(abs(GRHydro_C2P_failed(i,j,k)-2.0d0) .lt. 1.0d-10) then
            !$OMP CRITICAL
            if (reflevel.ge.GRHydro_c2p_warn_from_reflevel) then
              call CCTK_WARN(1,"Convergence problem in c2p")
              write(warnline,"(A10,i5)") "reflevel: ",reflevel
              call CCTK_WARN(1,warnline)
              write(warnline,"(A10,i10)") "iteration: ",cctk_iteration
              call CCTK_WARN(1,warnline)
              write(warnline,"(A10,F5.2)") "error: ", GRHydro_C2P_failed(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,"(3i5,1P10E15.6)") i,j,k,x(i,j,k),y(i,j,k),z(i,j,k),r(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,"(1P10E15.6)") rho(i,j,k),dens(i,j,k),eps(i,j,k),&
                    temperature(i,j,k),y_e(i,j,k),press(i,j,k)
              call CCTK_WARN(1,warnline)
              call EOS_Omni_press(GRHydro_eos_handle,1,GRHydro_eos_rf_prec,1,&
                    rho(i,j,k),eps(i,j,k),temperature(i,j,k),&
                    y_e(i,j,k),press(i,j,k),keyerr,anyerr)
              write(warnline,"(1P10E15.6)") eps(i,j,k),press(i,j,k)
              call CCTK_WARN(1,warnline)
              call CCTK_ERROR("Aborting!!!")
              STOP
            endif
            !$OMP END CRITICAL
          endif
        endif
        
        if(abs(GRHydro_C2P_failed(i,j,k)-3.0d0) .lt. 1.0d-10 .or. &
           temperature(i,j,k).lt.GRHydro_hot_atmo_temp) then
          ! dropped out off the EOS table in temperature or below
          ! the temperature of the atmosphere.
          ! Now reset this point to minimum temperature.
          GRHydro_C2P_failed(i,j,k) = 0
          if (reflevel.ge.GRHydro_c2p_warn_from_reflevel) then
            !$OMP CRITICAL
            write(warnline,"(A10,i7,4i5,1P10E15.6)") "reset T:",&
                  cctk_iteration,GRHydro_Reflevel,&
                  i,j,k,rho(i,j,k),&
                  eps(i,j,k),temperature(i,j,k),y_e(i,j,k)
            call CCTK_WARN(1,warnline)
            write(warnline,"(A10,i7,4i5,1P10E15.6)") "reset T:",&
                  cctk_iteration,GRHydro_Reflevel,&
                  i,j,k,x(i,j,k),y(i,j,k),z(i,j,k),r(i,j,k)
            call CCTK_WARN(1,warnline)
            !$OMP END CRITICAL
          endif
          temperature(i,j,k) = GRHydro_hot_atmo_temp
          keytemp = 1
          call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                rho(i,j,k),eps(i,j,k),temperature(i,j,k),y_e(i,j,k),&
                press(i,j,k),keyerr,anyerr)
          keytemp = 0
          if(anyerr.ne.0) then
            !$OMP CRITICAL
            call CCTK_ERROR("EOS Problem in C2P hot!!!")
            STOP
            !$OMP END CRITICAL
          endif
          !prim2con
          w_lorentz(i,j,k) = &
            1.d0 / sqrt(1.d0 - (gxx(i,j,k)*vel(i,j,k,1)**2 &
            + gyy(i,j,k)*vel(i,j,k,2)**2 &
            + gzz(i,j,k) *vel(i,j,k,3)**2 &
            + 2.0d0*gxy(i,j,k)*vel(i,j,k,1)*vel(i,j,k,2) &
            + 2.0d0*gxz(i,j,k)*vel(i,j,k,1)*vel(i,j,k,3)  &
            + 2.0d0*gyz(i,j,k)*vel(i,j,k,2)*vel(i,j,k,3)))
          vlowx = gxx(i,j,k)*vel(i,j,k,1) &
            + gxy(i,j,k)*vel(i,j,k,2)  &
            + gxz(i,j,k)*vel(i,j,k,3)
          vlowy = gxy(i,j,k)*vel(i,j,k,1) &
            + gyy(i,j,k)*vel(i,j,k,2)  &
            + gyz(i,j,k)*vel(i,j,k,3)
          vlowz = gxz(i,j,k)*vel(i,j,k,1) &
            + gyz(i,j,k)*vel(i,j,k,2)  &
            + gzz(i,j,k)*vel(i,j,k,3)
          scon(i,j,k,1) = sdet * (rho(i,j,k)*(1.0d0+eps(i,j,k)) &
            + press(i,j,k))*w_lorentz(i,j,k)**2 * vlowx
          scon(i,j,k,2) = sdet * (rho(i,j,k)*(1.0d0+eps(i,j,k)) &
            +press(i,j,k))*w_lorentz(i,j,k)**2 * vlowy
          scon(i,j,k,3) = sdet * (rho(i,j,k)*(1.0d0+eps(i,j,k)) &
            + press(i,j,k))*w_lorentz(i,j,k)**2 * vlowz
          tau(i,j,k) = sdet * ((rho(i,j,k)*(1.0d0+eps(i,j,k)) &
            + press(i,j,k))*w_lorentz(i,j,k)**2 - press(i,j,k)) &
            - dens(i,j,k)
        endif
      end do
    end do
  end do
  !$OMP END PARALLEL DO

  return

end subroutine Conservative2PrimitiveHot


subroutine Con2Prim_3DNR_1DD_hot(cctk_iteration, myproc, ii,jj,kk,handle, dens, &
    sx, sy, sz, tau, ye_con, rho, velx, vely, &
    velz, epsilon, temp, ye, press, w_lorentz, uxx, uxy, uxz, uyy, &
    uyz, uzz, sdet, x, y, z, r, epsnegative, GRHydro_rho_min, pmin, epsmin, &
    GRHydro_reflevel, GRHydro_C2P_failed, local_perc_ptol, W_guess, z_guess, T_guess)

  implicit none

  DECLARE_CCTK_PARAMETERS
  !DECLARE_CCTK_FUNCTIONS

  CCTK_REAL dens, sx, sy, sz, tau, ye_con, rho, velx, vely, velz, epsilon, &
      press, uxx, uxy, uxz, uyy, uyz, uzz, invsdet, sdet, w_lorentz, x, &
      y, z, r, GRHydro_rho_min
  CCTK_REAL temp, ye
  CCTK_REAL s2, f, df, vlowx, vlowy, vlowz
  CCTK_INT cctk_iteration, ii,jj,kk,count, i, handle, GRHydro_reflevel
  CCTK_REAL GRHydro_C2P_failed
  CCTK_REAL udens, usx, usy, usz, utau, uye_con, pold, pnew, &
           temp1, drhobydpress, depsbydpress, dpressbydeps, dpressbydrho, pmin, epsmin
  CCTK_REAL pminl,plow,tmp, dummy1, dummy2
  CCTK_REAL local_perc_ptol

  CCTK_REAL :: W_guess, z_guess, T_guess

  CCTK_REAL dpf

  integer :: myproc

  character(len=256) warnline
  logical epsnegative, mustbisect

  integer :: failwarnmode
  integer :: failinfomode

  !begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress,temp0
  CCTK_INT  :: nfudgemax,nf
  n=1;keytemp=0;anyerr=0;keyerr(1)=0
  temp0 = 0.0d0;xpress = 0.0d0
  nf=0;nfudgemax=30
  !end EOS Omni vars

  logical :: NR_failed, Dekker_failed

  mustbisect = .false.
  
  failinfomode = 1
  failwarnmode = 0

  if(con2prim_oct_hack.ne.0.and.&
      (x .lt. -1.0d-12 .or.&
       y .lt. -1.0d-12 .or.&
       z .lt. -1.0d-12)) then
    failwarnmode = 2
    failinfomode = 2
  end if

  !!$Undensitize the variables
  invsdet = 1.0d0/sdet
  udens = dens * invsdet
  usx = sx * invsdet
  usy = sy * invsdet
  usz = sz * invsdet
  utau = tau * invsdet
  uye_con = ye_con * invsdet
  s2 = usx*usx*uxx + usy*usy*uyy + usz*usz*uzz + 2.*usx*usy*uxy + &
      2.*usx*usz*uxz + 2.*usy*usz*uyz

  !Try 3D Newton Raphson in W,z,T
  NR_failed = .true.

  call Con2Prim_3DNR_hot(W_guess, z_guess, T_guess, W_result, z_result, T_result,&
                         udens, usx, usy, usz, s2, utau, uye_con)
  xrho = udens/W_result
  xYe = uye_con/udens
  call EOS_Omni_press(GRHydro_eos_handle,1,GRHydro_eos_rf_prec,1,&
          xrho,xeps,T_result,xYe,xpress,keyerr,anyerr)
  xs2 = (z_result**2) * (1.0d0 - (W_result**(-2.0d0)))
  xh = 1.0d0 + xeps + xpress/xrho
  xtau = xrho * xh * (W_result**2.0d0) - xpress - udens
  
  error_s2  = 2.0d0*abs(s2-xs2)/(s2+xs2)
  error_tau = 2.0d0*abs(tau-xtau)/(tau+xtau)

  if ((error_s2.gt.local_perc_ptol) .or. (error_tau.gt.local_perc_ptol)) then
    NR_failed = .true.
  end if

  if(NR_failed) then
    !If NR fails, try 1D Dekker in h*W
    Dekker_failed = .true.
    call Con2Prim_1DD_hot(dens, sx, sy, sz, tau, ye_con)

    if(Dekker_failed) then
      !If Dekker also failed, abort.
      !$OMP CRITICAL
      write(warnline,"Dekker fallback failed, aborting.")
      call CCTK_WARN(1,warnline)
      call CCTK_ERROR("C2P hot failed, aborting.")
      STOP
      !$OMP END CRITICAL
    end if
  end if
    

  
end subroutine Con2Prim_3DNR_1DD_hot

subroutine Con2Prim_3DNR_hot(W_guess, z_guess, T_guess, W_result, z_result, T_result, &
                             udens, usx, usy, usz, s2, utau, uye_con, success)
  
  implicit none

  DECLARE_CCTK_PARAMETERS
  
  CCTK_REAL, intent(in) :: W_guess, z_guess, T_guess
  CCTK_REAL, intent(out) :: W_result, z_result, T_result
  CCTK_REAL, intent(in) :: udens, usx, usy, usz, s2, utau, uye_con
  
  CCTK_REAL :: W_current, z_current, T_current
  CCTK_REAL :: Ye_current
  CCTK_REAL :: cons1, cons2, cons3

  CCTK_INT :: anyerr = 0
  CCTK_INT :: keyerr = 0

  CCTK_REAL :: press_current, eps_current
  CCTK_REAL :: dpdW, dpdT, dedW, dedT
  CCTK_REAL :: J1W,J1z,J1T,J2W,J2z,J2T,J3W,J3z,J3T
  CCTK_REAL :: iJ1W,iJ1z,iJ1T,iJ2W,iJ2z,iJ2T,iJ3W,iJ3z,iJ3T
  CCTK_REAL :: delta_W, delta_z, delta_T
  CCTK_REAL :: delta_W_old, delta_z_old, delta_T_old

  logical, intent(out) :: success
  logical :: done

  integer :: count

  done = .false.
  success = .false.

  Ye_current = uye_con/udens

  W_current = W_guess
  z_current = z_guess
  T_current = T_guess

  count = 0
  !loop
  do while(done .eq. .false.)
    !Get press and eps
    call EOS_Omni_press(GRHydro_eos_handle,1,GRHydro_eos_rf_prec,1,&
          rho_current,eps_current,T_current,Ye_current,press_current,keyerr,anyerr)

    !Calculate constraints
    call Con2Prim_3DNR_hot_cons1(W_current,z_current,udens,utau,press_current,cons1)
    call Con2Prim_3DNR_hot_cons1(W_current,z_current,s2,cons2)
    call Con2Prim_3DNR_hot_cons1(W_current,z_current,udens,press_current,eps_current,cons3)

    !Calculate Jacobian
    dpdW = !TODO
    dpdT = !TODO
    dedW = !TODO
    dedT = !TODO

    J1W = 2.0d0*W_current*(utau + udens - z_current + press_current) + (W_current**2.0d0)*dpdW
    J1z = -(W_current**2.0d0)
    J1T = (W_current**2.0d0)*dpdT

    J2W = 2.0d0*W_current*((z_current**2.0d0) - s2)
    J2z = 2.0d0*z_current*((W_current**2.0d0) - 1.0d0)
    J2T = 0.0d0

    J3W = - (z_current/(udens*(W_current**2))) - ((press_current + W*dpdW)/udens) - dedW
    J3z = 1.0d0/(udens*W_current)
    J3T = - dedT - ((W_current/udens)*dpdT)

    !Invert Jacobian
    call Con2Prim_3DNR_InvertJacobian(J1W,J1z,J1T,&
                                    J2W,J2z,J2T,&
                                    J3W,J3z,J3T,&
                                    iJ1W,iJ2W,iJ3W,&
                                    iJ1z,iJ2z,iJ3z,&
                                    iJ1T,iJ2T,iJ3T))
  
    !Calculate deltas
    delta_W = cons1*iJW1 + cons2*iJW2 + cons3*iJW3
    delta_z = cons1*iJz1 + cons2*iJz2 + cons3*iJz3
    delta_T = cons1*iJT1 + cons2*iJT2 + cons3*iJT3
    
    !Update WzT
    W_current = W_current - delta_W
    z_current = z_current - delta_z
    T_current = T_current - delta_T

    !Check bounds and oscillation
    if (W_current.lt.1.0d0) then
      W_current = 1.0d0
    else if ((delta_W*delta_W_old.lt.0.0d0) .and. &
      (delta_W+delta_W_old.lt.1.0d-6*abs(delta_W))) then
      W_current = W_current + 0.5d0*delta_W
      delta_W = 0.5d0*delta_W
    end if

    if ((delta_z*delta_z_old.lt.0.0d0) .and. &
      (delta_z+delta_z_old.lt.1.0d-6*abs(delta_z))) then
      z_current = z_current + 0.5d0*delta_z
      delta_z = 0.5d0*delta_z
    end if

    if (T_current.lt.eos_compose_temp_min) then
      T_current = eos_compose_temp_min
    else if (T_current.gt.eos_compose_temp_max) then
      T_current = eos_compose_temp_max
    else if ((delta_T*delta_T_old.lt.0.0d0) .and. &
      (delta_T+delta_T_old.lt.1.0d-6*abs(delta_T))) then
      T_current = T_current + 0.5d0*delta_T
      delta_T = 0.5d0*delta_T
    end if
    
    delta_W_old = delta_W
    delta_z_old = delta_z
    delta_T_old = delta_T

    !Update counter
    count = count + 1
    !Check whether done
    if(((abs(delta_W)/W_current).lt.local_perc_ptol) .and. &
       ((abs(delta_z)/z_current).lt.local_perc_ptol) .and. &
       ((abs(delta_T)/T_current).lt.local_perc_ptol)) then
      done = .true.
      success = .true.

    else if(count.gt.GRHydro_countmax) then
      done = .true.
      success = .false.

    end if

  end do

  W_result = W_current
  z_result = z_current
  T_result = T_current

  return
end subroutine Con2Prim_3DNR_hot

subroutine Con2Prim_3DNR_hot_cons1(W,z,D,tau,press,result)
  implicit none

  CCTK_REAL, intent(in)  :: W, z, D, tau, press
  CCTK_REAL, intent(out) :: result

  result = (tau + D - z + press)*(W**2)

end subroutine Con2Prim_3DNR_hot_cons1

subroutine Con2Prim_3DNR_hot_cons2(W,z,s2,result)

  implicit none
  
  CCTK_REAL, intent(in)  :: W, z, s2
  CCTK_REAL, intent(out) :: result

  result = (((z**2) - s2)*(W**2)) - (z**2)

end subroutine Con2Prim_3DNR_hot_cons2

subroutine Con2Prim_3DNR_hot_cons3(W,z,D,press,eps,result)

  implicit none

  CCTK_REAL, intent(in)  :: W, z, D, press, eps
  CCTK_REAL, intent(out) :: result

  result = ((z - D*W - press*(W**2))/(D*W)) - eps

end subroutine Con2Prim_3DNR_hot_cons3

subroutine Con2Prim_3DNR_InvertJacobian(J1W,J1z,J1T,J2W,J2z,J2T,J3W,J3z,J3T,&
            iJ1W,iJ2W,iJ3W,iJ1z,iJ2z,iJ3z,iJ1T,iJ2T,iJ3T)

  implicit none

  CCTK_REAL :: J1W,J1z,J1T
  CCTK_REAL :: J2W,J2z,J2T
  CCTK_REAL :: J3W,J3z,J3T

  CCTK_REAL :: iJ1W,iJ2W,iJ3W
  CCTK_REAL :: iJ1z,iJ2z,iJ3z
  CCTK_REAL :: iJ1T,iJ2T,iJ3T

  CCTK_REAL :: A,B,C,D,E,F,G,H,I
  CCTK_REAL :: invDetMat

  A =  (J2z*J3T - J2T*J3z)
  B = -(J2W*J3T - J2T*J3W)
  C =  (J2W*J3z - J2z*J3W)
  D = -(J1z*J3T - J1T*J3z)
  E =  (J1W*J3T - J1T*J3W)
  F = -(J1W*J3z - J1z*J3W)
  G =  (J1z*J2T - J1T*J2z)
  H = -(J1W*J2T - J1T*J2W)
  I =  (J2W*J3z - J2z*J3W)

  invDetMat = 1.0d0/(J1W*A + J2z*B + J3T*C)

  iJ1W = A*invDetMat
  iJ2W = D*invDetMat
  iJ3W = G*invDetMat

  iJ1z = B*invDetMat
  iJ2z = E*invDetMat
  iJ3z = H*invDetMat

  iJ1T = C*invDetMat
  iJ2T = F*invDetMat
  iJ3T = I*invDetMat

  return
end subroutine Con2Prim_3DNR_InvertJacobian

subroutine Con2Prim_1DD_hot(dens, sx, sy, sz, tau, ye_con)
  
  implicit none

end subroutine Con2Prim_1DD_hot