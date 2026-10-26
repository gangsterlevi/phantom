!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setfixedentropycore
!
! This module replaces the core of a MESA stellar profile with a flat-
! entropy profile that is in hydrostatic equilibrium with an added sink
! particle.
!
! :References:
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: dim, eos, io, kernel, physcon, table_utils
!
 implicit none
 integer :: ientropy
 public :: set_fixedS_softened_core,calc_mass_from_rho,gcore

 private
 integer, parameter :: ierr_pres=1,ierr_rho=2,ierr_mass=3,ierr_convergence=4

contains

!-----------------------------------------------------------------------
!+
!  Main subroutine that calculates the constant entropy softened profile
!+
!-----------------------------------------------------------------------
subroutine set_fixedS_softened_core(eos_type,mcore,rcore,rho,r,pres,m,Xcore,Ycore,ierr)
 use dim,         only:do_radiation
 use physcon,     only:pi,gg,solarm,solarr
 use table_utils, only:interpolator
 use io,          only:fatal
 integer, intent(in)    :: eos_type
 real,    intent(inout) :: r(:),rho(:),m(:),pres(:),mcore
 real,    intent(in)    :: rcore,Xcore,Ycore
 integer, intent(out)   :: ierr
 real, allocatable    :: r_alloc(:),rho_alloc(:),pres_alloc(:),m_alloc(:)
 real                 :: mc,msoft,rc
 integer              :: icore,iverbose

 ierr = 0
 if (mcore < 0.) then
    ierr = ierr_convergence
    return
 endif
 rc = rcore*solarr  ! convert to cm
 mc = mcore*solarm  ! convert to g
 call interpolator(r,rc,icore)  ! find index in r closest to rc
 msoft = m(icore) - mc
 if (msoft<=0.) call fatal('setup','mcore cannot exceed m(r=h)')

 if (do_radiation) then
    ientropy = 2
 else
    select case(eos_type)
    case(2)
       ientropy = 1
    case(10,12,20)
       ientropy = 2
    case default
       call fatal('setfixedentropycore',&
                   'eos_type not one of 2 (ideal gas), 12 (ideal plus rad.), 10 (MESA), or 20 (gas+rad+recombination)')
    end select
 endif

 ! Make allocatable copies, see instructions of calc_rho_and_pres
 allocate(r_alloc(0:icore+1))
 r_alloc(0) = 0.
 r_alloc(1:icore+1) = r(1:icore+1)
 allocate(rho_alloc(0:icore))
 rho_alloc(icore) = rho(icore)
 allocate(pres_alloc(0:icore+1))
 pres_alloc(icore:icore+1) = pres(icore:icore+1)
 iverbose = 0
 call calc_rho_and_pres(r_alloc,mc,m(icore),rho_alloc,pres_alloc,Xcore,Ycore,iverbose,ierr)
 if (ierr /= 0) return
 mcore = mc / solarm
 write(*,'(1x,a,f12.5,a)') 'Obtained core mass of ',mcore,' Msun'
 write(*,'(1x,a,f12.5,a)') 'Softened mass is ',m(icore)/solarm-mcore,' Msun'
 rho(1:icore)  = rho_alloc(1:icore)
 pres(1:icore) = pres_alloc(1:icore)
 ! Include the central shell, using the same quadrature as the shooting solve.
 allocate(m_alloc(icore+1))
 call calc_mass_from_rho(r_alloc(0:icore),rho_alloc(0:icore),m_alloc)
 m(1:icore) = m_alloc(2:icore+1) + mc

end subroutine set_fixedS_softened_core

!-----------------------------------------------------------------------
!+
!  Returns softened core profile with fixed entropy
!+
!-----------------------------------------------------------------------
subroutine calc_rho_and_pres(r,mcore,mh,rho,pres,Xcore,Ycore,iverbose,ierr_out)
 use eos, only:entropy,get_mean_molecular_weight
 use, intrinsic :: ieee_arithmetic, only:ieee_is_finite
 real, allocatable, intent(in)    :: r(:)
 integer,           intent(in)    :: iverbose
 real,              intent(in)    :: mh,Xcore,Ycore
 real,              intent(inout) :: mcore
 real, allocatable, intent(inout) :: rho(:),pres(:)
 integer, intent(out)                           :: ierr_out
 integer                                        :: Nmax,it,ierr
 real                                           :: Sc,mass,msoft,mu,mc_lo,mc_hi,mc_next
 integer, parameter :: it_max = 100
 real, parameter :: tol_mass = 1.e-10

! INSTRUCTIONS

! Input variables should be given in the following format:

! r(0:Nmax+1): Array of radial grid to be softened, satisfying r(0)=0 and r(Nmax)=rcore
! mcore:       Core particle mass, need to provide initial guess
! mh:          Mass coordinate at rcore, m(r=rcore)
! rho(0:Nmax): Give rho(Nmax)=(rho at rcore) as input. Outputs density profile.
! p(0:Nmax+1): Give p(Nmax:Nmax+1)=(p at r(Nmax:Nmax+1)) as input. Outputs pressure profile.

 msoft = mh - mcore
 Nmax  = size(rho)-1 ! Index corresponding to r = h
 mu = get_mean_molecular_weight(Xcore,1.-Xcore-Ycore)
 Sc = entropy(rho(Nmax),pres(Nmax),mu,ientropy)
 rho(0:Nmax-1) = rho(Nmax)
 pres(0:Nmax-1) = pres(Nmax)

 ! Bisect the physical core-mass interval. An exhausted gas mass means
 ! the trial core is too massive; a positive central residual means too small.
 ! Unlike multiplicative updates, this remains bounded and cannot stall at zero.
 mc_lo = 0.
 mc_hi = mh
 ierr_out = ierr_convergence
 do it = 1,it_max
    msoft = mh - mcore
    ierr = 0
    call one_shot(Sc,r,mcore,msoft,mu,rho,pres,mass,iverbose,ierr)
    if (.not.ieee_is_finite(mass) .or. any(.not.ieee_is_finite(rho)) .or. &
        any(.not.ieee_is_finite(pres))) exit
    ! A density inversion is diagnostic, as in the original shooting method;
    ! pressure failure or an incomplete integration must never count as closure.
    if ((ierr==0 .or. ierr==ierr_rho) .and. abs(mass) <= tol_mass*msoft) then
       ierr_out = 0
       write(*,'(1x,a,i0,a,es12.5)') 'Core mass converged in ',it, &
                                   ' iterations; central mass / softened mass = ',mass/msoft
       return
    endif
    if (ierr==ierr_pres) then
       ierr_out = ierr
       exit
    endif
    if (mass < 0.) then
       mc_hi = mcore
    else
       mc_lo = mcore
    endif
    mc_next = mc_lo + 0.5*(mc_hi-mc_lo)
    if (mc_next==mcore .or. mc_next==mc_lo .or. mc_next==mc_hi) exit
    mcore = mc_next
    if (iverbose > 0) write(*,'(1x,i5,2(a,e12.5))') it,' m(r=0) = ',mass,' mcore = ',mcore
 enddo
 write(*,'(/,1x,a,i0,a,es12.5)') 'WARNING: Constant-entropy core solve failed (status ', &
                               ierr_out,'); central mass / softened mass = ',mass/msoft

end subroutine calc_rho_and_pres

!-----------------------------------------------------------------------
!+
!  Calculate a hydrostatic structure for a given entropy
!+
!-----------------------------------------------------------------------
subroutine one_shot(Sc,r,mcore,msoft,mu,rho,pres,mass,iverbose,ierr)
 use physcon, only:gg,pi,solarm
 use eos,     only:get_rho_from_p_s
 real,              intent(in)    :: Sc,mcore,msoft,mu
 integer,           intent(in)    :: iverbose
 real, allocatable, intent(in)    :: r(:)
 real, allocatable, intent(inout) :: rho(:),pres(:)
 real,              intent(out)   :: mass
 integer,           intent(out)   :: ierr
 integer                                        :: i,Nmax
 real                                           :: rcore,rhoguess
 real, allocatable :: dr(:),dvol(:)

 Nmax = size(rho)-1
 ierr = 0
 allocate(dr(1:Nmax+1),dvol(1:Nmax+1))

 ! Pre-fill arrays
 do i = 1,Nmax+1
    dr(i) = r(i)-r(i-1)
    dvol(i) = 4.*pi/3. * (r(i)**3 - r(i-1)**3)
 enddo

 rcore = r(Nmax)
 mass  = msoft

 do i = Nmax, 1, -1
    pres(i-1) = ( dr(i) * dr(i+1) * sum(dr(i:i+1)) &
                * rho(i) * gg * (mass/r(i)**2 + mcore * gcore(r(i),rcore)) &
                + dr(i)**2 * pres(i+1) &
                + ( dr(i+1)**2 - dr(i)**2) * pres(i) ) / dr(i+1)**2
    rhoguess = rho(i)
    call get_rho_from_p_s(pres(i-1),Sc,rho(i-1),mu,rhoguess,ientropy)
    mass = mass - 0.5*(rho(i)+rho(i-1)) * dvol(i)

    if (iverbose > 2) print*,Nmax-i+1,pres(i-1),rhoguess,rho(i-1),mass
    if (mass < 0.) then ! m(r) < 0 encountered, exit and decrease mcore
       if (iverbose > 1) print*,'WARNING: Negative mass reached at i = ',i, 'm = ',mass/solarm
       ierr = ierr_mass
       return
    endif
    if (rho(i-1)<rho(i)) then
       if (iverbose > 1) then
          print*,'WARNING: Density inversion at i = ',i, 'm = ',mass/solarm
          write(*,'(i5,2x,e12.4,2x,e12.4,2x,e12.4)') i,rho(i),rho(i-1),mass
       endif
       ierr = ierr_rho
    endif
    if (pres(i-1)<pres(i)) then
       if (iverbose > 1) then
          print*,'WARNING: Pressure inversion at i = ',i, 'm = ',mass/solarm
          write(*,'(i5,2x,e12.4,2x,e12.4,2x,e12.4)') i,pres(i-1),rho(i),mass
       endif
       ierr = ierr_pres
       return
    endif
 enddo

end subroutine one_shot

!-----------------------------------------------------------------------
!+
!  Acceleration from softened potential of stellar core
!+
!-----------------------------------------------------------------------
function gcore(r,rcore)
 use kernel, only:kernel_softening
 real, intent(in) :: r,rcore
 real             :: gcore,dum
 real             :: hsoft,q

 hsoft = 0.5 * rcore
 q = r / hsoft
 call kernel_softening(q**2,q,dum,gcore)
 gcore = gcore / hsoft**2 ! Note: gcore is not multiplied by G or mcore yet.

end function gcore

subroutine calc_mass_from_rho(r,rho,m)
 use physcon, only:pi
 real, intent(in)    :: r(:),rho(:)
 real, intent(inout) :: m(:)
 real                :: densi
 integer             :: i

 m(1) = 0
 do i = 1,size(r)-1
    densi = 0.5*(rho(i) + rho(i+1))
    m(i+1) = m(i) + 4./3. * pi * (r(i+1)**3 - r(i)**3) * densi
 enddo

end subroutine calc_mass_from_rho

end module setfixedentropycore
