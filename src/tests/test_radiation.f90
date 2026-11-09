!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testradiation
!
! Unit tests for radiation hydro
!
! :References:
!   - Whitehouse & Bate (2004), 353, 1078
!   - Whitehouse, Bate & Monaghan (2005), 364, 1367
!   - Biriukov (2019), PhD thesis, Monash Univ.
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, deriv, dim, eos, implicit, io, kernel, mpidomain,
!   mpiutils, neighkdtree, options, part, physcon, radiation_implicit,
!   radiation_utils, readwrite_dumps, sortutils, step_lf_global, testutils,
!   timestep, unifdis, units
!
 use part,      only:ithick,iradxi,ifluxx,ifluxy,ifluxz,ikappa
 use io,        only:id,master
 use testutils, only:checkval,update_test_scores,checkvalbuf,checkvalbuf_end
 implicit none

 public :: test_radiation
 private

contains
!-------------------------------------------------
!+
!  unit tests of radiation hydrodynamics
!+
!-------------------------------------------------
subroutine test_radiation(ntests,npass)
 use physcon, only:solarm,au
 use units,   only:set_units
 use dim,     only:do_radiation,periodic
 use io,      only:iverbose
 use radiation_utils, only:implicit_radiation,tol_rad
 integer, intent(inout) :: ntests,npass

 if (.not.do_radiation) then
    if (id==master) write(*,"(/,a,/)") '--> SKIPPING RADIATION TEST (NEED RADIATION=yes)'
    return
 endif
 if (id==master) write(*,"(/,a,/)") '--> TESTING RADIATION MODULE'

 call test_implicit_workspace_growth(ntests,npass)
 call set_units(dist=au,mass=solarm,G=1.d0)
 call test_implicit_thin_diffusion(ntests,npass)
 call test_exchange_terms(ntests,npass,use_implicit=.false.)
 call test_exchange_terms(ntests,npass,use_implicit=.true.)

 if (.not.periodic) then
    if (id==master) write(*,"(/,a)") '--> SKIPPING TEST OF RADIATION DERIVS (need -DPERIODIC)'
 else
    call test_implicit_matches_explicit(ntests,npass)

    iverbose = -1
    implicit_radiation = .false.
    call test_radiation_diffusion(ntests,npass)

    implicit_radiation = .true.
    tol_rad = 1.e-6
    call test_radiation_diffusion(ntests,npass)
    call test_radiation_diffusion(ntests,npass,hierarchy=.true.)
 endif

 if (id==master) write(*,"(/,a)") '<-- RADIATION TEST COMPLETE'

end subroutine test_radiation

! Setup adds particles as each star is constructed. The implicit workspace
! must grow between solves and also return a defined status on reuse.
subroutine test_implicit_workspace_growth(ntests,npass)
 use implicit, only:allocate_memory_implicit,ivar,vari,varinew,mask,ijvar,varij,varij2,icompactmax
 integer, intent(inout) :: ntests,npass
 integer :: ierr,nold,nnew,nerr(6)

 call allocate_memory_implicit(8,2.,1.2,ierr)
 nold = size(ivar,2)
 nnew = nold+8
 call allocate_memory_implicit(nnew,2.,1.2,ierr)
 call checkval(ierr,0,0,nerr(1),'implicit workspace grows successfully')
 call checkval(size(ivar,2),nnew,0,nerr(2),'implicit particle index capacity grows')
 call checkval(size(vari,2),nnew,0,nerr(3),'implicit per-particle state capacity grows')
 call checkval(size(varinew,2),nnew,0,nerr(4),'implicit new state capacity grows')
 call checkval(size(mask),nnew,0,nerr(5),'implicit particle mask capacity grows')
 call checkval(min(size(ijvar),size(varij,2),size(varij2,2)),icompactmax,0,nerr(6), &
               'implicit neighbour storage grows consistently')
 call update_test_scores(ntests,nerr,npass)
 ierr = -1
 call allocate_memory_implicit(4,2.,1.2,ierr)
 call checkval(ierr,0,0,nerr(1),'implicit workspace reuse returns success')
 call checkval(size(ivar,2),nnew,0,nerr(2),'implicit workspace retains sufficient capacity')
 call update_test_scores(ntests,nerr(1:2),npass)

end subroutine test_implicit_workspace_growth

! An optically thin particle between fixed radiation reservoirs can cycle
! as its flux limiter changes. Check convergence against the diffusion
! equation, including when the gas has no energy to exchange.
subroutine test_implicit_thin_diffusion(ntests,npass)
 use part, only:init_part,npart,npartoftype,massoftype,hfact,xyzh,vxyzu,rho,rad,radprop,drad,&
                iphase,isetphase,igas,iboundary,gradh,luminosity,init_rho_from_h
 use boundary, only:set_boundary
 use kernel, only:grkern,cnormk,hfact_default
 use units, only:unit_opacity,get_c_code
 use eos, only:iopacity_type,get_cv,gamma,gmw
 use radiation_utils, only:kappa_cgs,cv_type,tol_rad,itsmax_rad,exchange_radiation_energy,&
                          limit_radiation_flux,get_rad_R
 use radiation_implicit, only:do_radiation_implicit,rad_errorE
 use neighkdtree, only:build_tree
 use mpiutils, only:reduceall_mpi
 integer, intent(inout) :: ntests,npass
 integer :: ierr,j,nerr(5),opacity_prev,cv_prev,limit_prev
 real :: kappa_prev,tol_prev,gamma_prev,gmw_prev,dt,c_code,gradient(3),lambda,coeff,residual,diagonal,dW,r,xi
 real, parameter :: kappa_code = 1.e-7,xi0 = 1.0001
 logical :: exchange_prev,limiter_prev

 opacity_prev = iopacity_type
 cv_prev = cv_type
 kappa_prev = kappa_cgs
 gamma_prev = gamma
 gmw_prev = gmw
 tol_prev = tol_rad
 limit_prev = itsmax_rad
 exchange_prev = exchange_radiation_energy
 limiter_prev = limit_radiation_flux
 iopacity_type = 2
 cv_type = 0
 kappa_cgs = kappa_code*unit_opacity
 gamma = 5./3.
 gmw = 2.
 tol_rad = 1.e-6
 itsmax_rad = 32
 exchange_radiation_energy = .false.
 limit_radiation_flux = .true.

 call init_part()
 call set_boundary(-4.,4.,-4.,4.,-4.,4.)
 hfact = hfact_default
 massoftype(igas) = 1./hfact**3
 massoftype(iboundary) = massoftype(igas)
 if (id==master) then
    npart = 3
    npartoftype(igas) = 1
    npartoftype(iboundary) = 2
    iphase(1) = isetphase(igas,iactive=.true.)
    iphase(2:3) = isetphase(iboundary,iactive=.true.)
    xyzh(:,1:3) = 0.
    xyzh(1,2:3) = [1.,1.5]
    xyzh(4,1) = 1.
    xyzh(4,2:3) = 0.05
 endif
 call init_rho_from_h()
 gradh = 0.
 luminosity = 0.
 vxyzu(:,1:npart) = 0.
 vxyzu(4,1:npart) = 3000.*get_cv(0)
 if (id==master) rad(iradxi,1:3) = [xi0,1.0004/rho(2),0.9996/rho(3)]
 call build_tree(npart,npart,xyzh,vxyzu)
 c_code = get_c_code()
 dt = 2./c_code
 call do_radiation_implicit(dt,dt,npart,rad,xyzh,vxyzu,radprop,drad,ierr)
 call checkval(ierr,0,0,nerr(1),'optically thin diffusion converges without gas exchange')
 call checkval(rad_errorE,0.,tol_rad,nerr(2),'full radiation proposal satisfies tolerance')

 ! Assemble the backwards Euler defect for the gas particle independently.
 ! The reservoirs have no neighbours inside their own kernels, so lambda=1/3.
 residual = 0.
 xi = 0.
 if (id==master) then
    xi = rad(iradxi,1)
    gradient = 0.
    do j=2,3
       r = xyzh(1,j)
       dW = grkern(r*r,r)*cnormk
       gradient(1) = gradient(1) - massoftype(igas)*dW*(rho(j)*rad(iradxi,j)-rho(1)*xi)/rho(1)
    enddo
    r = get_rad_R(rho(1),xi,gradient,kappa_code)
    lambda = (2.+r)/(6.+3.*r+r*r)
    residual = xi0-xi
    diagonal = 1.
    do j=2,3
       r = xyzh(1,j)
       dW = grkern(r*r,r)*cnormk
       coeff = -0.5*dt*c_code*massoftype(igas)*dW/(r*rho(1)*rho(j))* &
               (lambda/(kappa_code*rho(1)) + 1./(3.*kappa_code*rho(j)))
       residual = residual + coeff*(rho(j)*rad(iradxi,j)-rho(1)*xi)
       diagonal = diagonal + coeff*rho(1)
    enddo
    residual = abs(residual)/(diagonal*xi)
 endif
 residual = reduceall_mpi('max',residual)
 call checkval(residual,0.,tol_rad,nerr(3),'optically thin backwards Euler defect')
 xi = reduceall_mpi('max',xi)
 call checkval(xi,1.000258953,1.e-6,nerr(4),'optically thin implicit diffusion solution')
 residual = 0.
 if (npart > 0) residual = maxval(abs(vxyzu(4,1:npart)/(3000.*get_cv(0))-1.))
 residual = reduceall_mpi('max',residual)
 call checkval(residual,0.,1.e-12,nerr(5),'diffusion without exchange preserves gas energy')
 call update_test_scores(ntests,nerr,npass)

 iopacity_type = opacity_prev
 cv_type = cv_prev
 kappa_cgs = kappa_prev
 gamma = gamma_prev
 gmw = gmw_prev
 tol_rad = tol_prev
 itsmax_rad = limit_prev
 exchange_radiation_energy = exchange_prev
 limit_radiation_flux = limiter_prev

end subroutine test_implicit_thin_diffusion

!----------------------------------------------------
!+
!  unit tests of gas-radiation energy exchange terms
!+
!----------------------------------------------------
subroutine test_exchange_terms(ntests,npass,use_implicit)
 use radiation_utils, only:update_radenergy,kappa_cgs,exchange_radiation_energy
 use units,      only:set_units,unit_ergg,unit_density,unit_opacity,utime
 use physcon,    only:au,solarm,seconds
 use dim,        only:maxp,periodic
 use io,         only:iverbose
 use apr,        only:sync_aprmassoftype
 use part,       only:init_part,npart,rho,xyzh,fxyzu,vxyzu,massoftype,igas,&
                      iphase,maxphase,isetphase,drad,&
                      npartoftype,rad,radprop,maxvxyzu,luminosity,init_rho_from_h,itlastrad
 use kernel,     only:hfact_default
 use unifdis,    only:set_unifdis
 use eos,        only:gmw,gamma,polyk,iopacity_type
 use boundary,   only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use mpiutils,   only:reduceall_mpi
 use mpidomain,  only:i_belong
 use radiation_implicit, only:do_radiation_implicit
 use neighkdtree,        only:build_tree
 integer, intent(inout) :: ntests,npass
 logical, intent(in)    :: use_implicit
 real :: psep,hfact
 real :: pmassi,rhozero,totmass
 real :: dt,t,physrho,rhoi,maxt,laste
 integer :: i,nerr(1),ndiff(1),ncheck,ierrmax,ierr,itest,N_implicit_steps
 integer(kind=8) :: nptot
 logical, parameter :: write_output = .false.
 character(len=12) :: string,filestr
 real :: step,log_start,u_before,xi_before

 call init_part()
 iverbose = 0
 exchange_radiation_energy = .false.
 iopacity_type = -1  ! preserve the opacity value
 kappa_cgs = 0.4

 if (use_implicit) then
    string = ' implicit'
    iopacity_type = 2
 else
    string = ' explicit'
 endif

 if (id==master) write(*,"(/,a)") '--> checking radiation exchange terms'//trim(string)

 psep = 1./16.
 hfact = hfact_default
 npart = 0
 call set_boundary(-0.5,0.5,-0.5,0.5,-0.5,0.5)
 call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                  psep,hfact,npart,xyzh,periodic,mask=i_belong)
 rhozero = 1.e-7/unit_density  ! 1e-7 g/cm^3
 totmass = rhozero*(dxbound*dybound*dzbound)
 nptot = reduceall_mpi('+',npart)
 massoftype(igas) = totmass/nptot
 call sync_aprmassoftype
 gamma = 5./3.
 gmw = 2.0
 polyk = 0.
 if (maxphase==maxp) iphase(:) = isetphase(igas,iactive=.true.)
 npartoftype(:) = 0
 npartoftype(1) = npart
 pmassi = massoftype(igas)
 call init_rho_from_h()

 if (use_implicit) call build_tree(npart,npart,xyzh,vxyzu)
 !
 ! first version of the test: set gas temperature high and radiation temperature low
 ! so that gas cools towards radiation temperature (itest=1)
 !
 ! second version: gas heats towards radiation temperature (itest=2)
 !
 do itest = 1,2
    do i=1,npart
       rhoi              = rho(i)
       rad(iradxi,i)     = 1e12/(unit_ergg*unit_density)/rhoi
       radprop(ikappa,i) = kappa_cgs/unit_opacity
       if (itest==2) then
          vxyzu(4,i)     = 1e2/(unit_ergg*unit_density)
       else
          vxyzu(4,i)     = 1e10/(unit_ergg*unit_density)
       endif
       vxyzu(4,i)        = vxyzu(4,i)/rhoi
       fxyzu(4,i)        = 0.
       luminosity(i)     = 0.
    enddo

    if (use_implicit) then
       ! With uniform radiation and no work terms, disabling exchange must
       ! preserve both energies even when the gas and radiation are not in LTE.
       exchange_radiation_energy = .false.
       u_before = vxyzu(4,1)
       xi_before = rad(iradxi,1)
       radprop(itlastrad,1:npart) = 0.
       call do_radiation_implicit(1.e-12*seconds/utime,1.e-12*seconds/utime,npart,rad,xyzh,vxyzu,radprop,drad,ierr)
       call checkval(ierr,0,0,nerr(1),'implicit exchange disabled: solver status')
       call update_test_scores(ntests,nerr,npass)
       call checkval(maxval(abs(vxyzu(4,1:npart)/u_before-1.)),0.,1.e-12,nerr(1),&
                     'implicit exchange disabled: gas energy unchanged')
       call update_test_scores(ntests,nerr,npass)
       call checkval(maxval(abs(rad(iradxi,1:npart)/xi_before-1.)),0.,1.e-12,nerr(1),&
                     'implicit exchange disabled: radiation energy unchanged')
       call update_test_scores(ntests,nerr,npass)
       exchange_radiation_energy = .true.
    endif

    if (write_output) then
       write(filestr,"(i2.2)") itest
       open(unit=23+itest,file='radiation-test'//trim(filestr)//trim(adjustl(string))//'.ev')
       write(23+itest,"(a)") '# time [s],Egas [erg/cm^3],Erad [erg/cm^3]'
    endif
    maxt = 5e-7*seconds
    t = 0.
    rhoi    = rho(1)
    physrho = rhoi*unit_density
    i = 0
    ndiff = 0
    ncheck = 0
    ierrmax = 0

    ! logarithmically spaced time steps for implicit
    if (use_implicit) then
       N_implicit_steps = 40
       log_start = log10(1.d-18*seconds/utime)
       step = (log10(maxt/utime) - log_start) / real(N_implicit_steps - 1)
    endif

    do while(t < maxt/utime)
       dt = max(1d-18*seconds/utime,0.05d0*t)
       if (t + dt > maxt/utime) dt = maxt/utime - t
       if (use_implicit) then
          if (i >1) dt = min(0.05*maxt/utime, 10.d0**(log_start + (i - 1) * step) - t)
          if (t + dt > maxt/utime) dt = maxt/utime - t  ! take last step to maxt
          ! integrate all particles from t to t+dt (also with individual timesteps)
          radprop(itlastrad,1:npart) = t
          call do_radiation_implicit(dt,t+dt,npart,rad,xyzh,vxyzu,radprop,drad,ierr)
          call checkvalbuf(ierr,0,0,'no errors from implicit solver',ndiff(1),ncheck,ierrmax)
       else
          if (t + dt > maxt/utime) dt = maxt/utime - t  ! take last step to maxt
          call update_radenergy(1,fxyzu,vxyzu,rho,rad,radprop,dt)
       endif
       t = t + dt
       if (mod(i,10)==0 .or. use_implicit) then
          laste = (vxyzu(4,1)*unit_ergg)*physrho
          if (write_output) write(23+itest,*) t*utime, laste,(rad(iradxi,1)*unit_ergg)*physrho
       endif
       i = i + 1
    enddo
    if (ncheck > 0) then
       call checkvalbuf_end('no errors from implicit solver',ncheck,ndiff(1),ierrmax,0)
       call update_test_scores(ntests,ndiff,npass)
    endif
    if (itest==2) then
       call checkval(laste,21144463.0313597,1.e-14,nerr(1),'energy exchange for gas heating'//trim(string))
    else
       call checkval(laste,21197127.9406196,1.e-14,nerr(1),'energy exchange for gas cooling'//trim(string))
    endif
    call update_test_scores(ntests,nerr,npass)

    if (write_output) close(23+itest)
 enddo

 if (use_implicit) call test_exchange_equilibrium_composition(ntests,npass)

end subroutine test_exchange_terms

! A uniform radiation field in LTE must remain in equilibrium even when
! particles have different compositions from the global EoS defaults.
subroutine test_exchange_equilibrium_composition(ntests,npass)
 use part, only:npart,xyzh,vxyzu,rad,radprop,drad,rho,eos_vars,iX,iZ,dvdx,luminosity,itlastrad
 use radiation_utils, only:cv_type,radxi_from_Trad
 use radiation_implicit, only:do_radiation_implicit
 use eos, only:use_var_comp
 use eos_gasradrec, only:init_eos_gasradrec,calc_uP_from_rhoT_gasradrec
 use units, only:unit_density,unit_ergg,utime
 integer, intent(inout) :: ntests,npass
 integer :: i,itest,ierr,nerr(3),cv_prev
 logical :: var_comp_prev
 real :: temp,ui,pi,imui
 real, allocatable :: u0(:),xi0(:)

 allocate(u0(npart),xi0(npart))
 cv_prev = cv_type
 var_comp_prev = use_var_comp
 cv_type = 20
 use_var_comp = .true.
 call init_eos_gasradrec(ierr)
 dvdx = 0.
 luminosity = 0.
 do itest=1,2
    temp = 10.**(3+itest)
    do i=1,npart
       eos_vars(iX,i) = 0.
       if (mod(i,2)==0) eos_vars(iX,i) = 0.35
       eos_vars(iZ,i) = 0.02
       call calc_uP_from_rhoT_gasradrec(rho(i)*unit_density,temp,eos_vars(iX,i),&
                                      1.-eos_vars(iX,i)-eos_vars(iZ,i),ui,pi,imui,do_radiation=.true.)
       vxyzu(4,i) = ui/unit_ergg
       rad(iradxi,i) = radxi_from_Trad(rho(i),temp)
    enddo
    u0 = vxyzu(4,1:npart)
    xi0 = rad(iradxi,1:npart)
    radprop(itlastrad,1:npart) = 0.
    call do_radiation_implicit(1.e-3/utime,1.e-3/utime,npart,rad,xyzh,vxyzu,radprop,drad,ierr)
    call checkval(ierr,0,0,nerr(1),'composition-dependent LTE: solver status')
    call checkval(maxval(abs(vxyzu(4,1:npart)/u0-1.)),0.,1.e-8,nerr(2),&
                  'composition-dependent LTE: gas energy unchanged')
    call checkval(maxval(abs(rad(iradxi,1:npart)/xi0-1.)),0.,1.e-8,nerr(3),&
                  'composition-dependent LTE: radiation energy unchanged')
    call update_test_scores(ntests,nerr,npass)
 enddo
 call test_implicit_iteration_limit(ntests,npass)
 cv_type = cv_prev
 use_var_comp = var_comp_prev
 deallocate(u0,xi0)

end subroutine test_exchange_equilibrium_composition

! Reaching the cap must return on every OpenMP thread, including when one
! thread reaches the next iteration before another checks whether to stop.
subroutine test_implicit_iteration_limit(ntests,npass)
 use part, only:npart,xyzh,vxyzu,rad,radprop,drad,itlastrad
 use radiation_utils, only:itsmax_rad,tol_rad
 use radiation_implicit, only:do_radiation_implicit,ierr_failed_to_converge,its_global
 use units, only:utime
 integer, intent(inout) :: ntests,npass
 integer :: limit,ierr,nerr(2),limit_prev
 real :: tol_prev
 real, allocatable :: u0(:),xi0(:)

 allocate(u0(npart),xi0(npart))
 u0 = vxyzu(4,1:npart)
 xi0 = rad(iradxi,1:npart)
 limit_prev = itsmax_rad
 tol_prev = tol_rad
 tol_rad = 1.e-14
 do limit=2,3
    vxyzu(4,1:npart) = u0
    rad(iradxi,1:npart) = xi0
    ! Drive a local gas-radiation imbalance requiring more than a few sweeps.
    rad(iradxi,1) = 2.*xi0(1)
    radprop(itlastrad,1:npart) = 0.
    itsmax_rad = limit
    call do_radiation_implicit(1.e-3/utime,1.e-3/utime,npart,rad,xyzh,vxyzu,radprop,drad,ierr)
    call checkval(ierr,ierr_failed_to_converge,0,nerr(1),'implicit iteration cap returns nonconvergence status')
    call checkval(its_global,limit,0,nerr(2),'all threads stop at the radiation iteration cap')
    call update_test_scores(ntests,nerr,npass)
 enddo
 vxyzu(4,1:npart) = u0
 rad(iradxi,1:npart) = xi0
 itsmax_rad = limit_prev
 tol_rad = tol_prev

end subroutine test_implicit_iteration_limit

!---------------------------------------------------------
!+
!  unit tests to check that derivatives computed
!  in the radiation_implicit routines match those
!  computed in the regular density/force routines
!  in an explicit calculation
!+
!---------------------------------------------------------
subroutine test_implicit_matches_explicit(ntests,npass)
 use part,     only:npart,xyzh,vxyzu,rad,radprop,drad,iorig,&
                    xyzh_label,init_part,radprop_label,apr_level
 use sortutils, only:indexx
 use dim,      only:use_apr
 use boundary, only:xmin,xmax,ymin,ymax,zmin,zmax
 use options,  only:tolh
 use physcon,  only:pi
 use deriv,    only:get_derivs_global
 use io,       only:iverbose
 use timestep, only:dtmax
 use radiation_utils,    only:implicit_radiation
 use radiation_implicit, only:do_radiation_implicit
 integer, intent(inout) :: ntests,npass
 real(kind=kind(radprop)), allocatable :: flux_explicit(:,:),flux_implicit(:,:)
 integer(kind=8), allocatable :: iorig_explicit(:)
 integer, allocatable :: inow(:),iexpl(:)
 real :: kappa_code,c_code,xi0,rho0,errmax_e,tol_e,tolh_old,pmassi !,exact(9)
 integer :: i,j,n,itry,nerr_e(9),ierr,imass,nmass
 character(len=24) :: label

 ! with APR, repeat with neighbouring particles of different mass
 nmass = 1
 if (use_apr) nmass = 2
 do imass=1,nmass
 label = ''
 if (imass == 2) label = ' (APR, unequal masses)'
 if (id==master) write(*,"(/,a)") '--> checking implicit routine matches explicit for flux terms'//trim(label)

 implicit_radiation = .false.
 iverbose = 0
 call init_part()
 tolh_old = tolh
 tolh = 1.e-8
 dtmax = 5.e-22
 !--set up once, so that both passes see the same particles
 call setup_radiation_diffusion_problem_sinusoid(kappa_code,c_code,xi0,rho0,pmassi)
 if (imass == 2) then
    do i=1,npart
       if (mod(i,2)==0) apr_level(i) = 2
    enddo
 endif

 ! set some non-zero velocities
 do i=1,npart
    vxyzu(1,i) = 0.1*sin((xyzh(1,i)-xmin)*2.*pi/(xmax-xmin))
    vxyzu(2,i) = 0.1*cos((xyzh(2,i)-ymin)*2.*pi/(ymax-ymin))
    vxyzu(3,i) = 0.1*sin((xyzh(3,i)-zmin)*2.*pi/(zmax-zmin))
 enddo
 do itry = 1,2
    if (itry == 2) implicit_radiation = .true.
    call get_derivs_global()
    if (itry==1) then
       call get_derivs_global()  ! twice to get density on neighbours correct

       !--allocate and copy flux, and the particle IDs
       allocate(flux_explicit(3,npart),iorig_explicit(npart))
       flux_explicit  = radprop(ifluxx:ifluxz,1:npart)
       iorig_explicit = iorig(1:npart)
    else
       radprop = 0.  ! also sets the time of the last update to zero
       call do_radiation_implicit(1.e-24,1.e-24,npart,rad,xyzh,vxyzu,radprop,drad,ierr)
    endif
 enddo

 ! now check that things match, matching the particles by their IDs
 ! since the particles may have been reordered in between
 tol_e = 1.e-15
 nerr_e = 0
 errmax_e = 0.
 ! (every MPI task makes the same checks, since the checks combine results across tasks)
 call checkval(npart,size(iorig_explicit),0,nerr_e(4),'number of particles')
 n = min(npart,size(iorig_explicit))
 allocate(inow(npart),iexpl(size(iorig_explicit)),flux_implicit(3,n))
 call indexx(npart,iorig,inow)
 call indexx(size(iorig_explicit),iorig_explicit,iexpl)
 call checkval(n,iorig(inow(1:n)),iorig_explicit(iexpl(1:n)),0_8,nerr_e(5),'same particle IDs')
 flux_explicit = flux_explicit(:,iexpl)
 flux_implicit = radprop(ifluxx:ifluxz,inow(1:n))
 do j=1,3
    call checkval(n,flux_implicit(j,:),flux_explicit(j,1:n),tol_e,nerr_e(j),radprop_label(ifluxx+j-1))
 enddo
 call update_test_scores(ntests,nerr_e,npass)
 deallocate(flux_explicit,iorig_explicit,flux_implicit,inow,iexpl)

 vxyzu(1:3,:) = 0.  ! zero velocities again
 tolh = tolh_old
 if (use_apr) apr_level(1:npart) = 1
 enddo

end subroutine test_implicit_matches_explicit

!---------------------------------------------------------
!+
!  unit tests of radiation derivatives: grad E and div F,
!  then diffusion of a sine wave, evolved either with step
!  or (hierarchy=.true.) by calling the implicit solver
!  directly, on a timestep hierarchy if ind_timesteps
!+
!---------------------------------------------------------
subroutine test_radiation_diffusion(ntests,npass,hierarchy)
 use part,            only:npart,xyzh,rho,vxyzu,&
                           radprop,drad,ifluxx,maxvxyzu,fxyzu,init_part
 use boundary,        only:xmin,xmax
 use physcon,         only:pi
 use readwrite_dumps, only:write_fulldump
 use deriv,           only:get_derivs_global
 use step_lf_global,  only:init_step,step
 use timestep,        only:dtmax
 use radiation_utils, only:implicit_radiation,limit_radiation_flux,implicit_radiation_store_drad
 use mpiutils,        only:reduceall_mpi
 integer, intent(inout) :: ntests,npass
 logical, intent(in), optional :: hierarchy
 real :: rhoi,dtext,pmassi,dt,t,kappa_code
 real :: xi0,D0,rho0,l0,c_code
 real :: dtnew,tmax
 real :: exact_grE,exact_DgrF
 real :: errmax_e,errmax_f,tol_e,tol_f,tol_xi,de,dekin,degas,derad
 integer :: i
 integer :: nactive,nerr_e(1),ncheck_e,nerr_f(1),ncheck_f
 character(len=20) :: string,filename
 logical, parameter :: write_output = .false.
 logical :: use_hierarchy

 use_hierarchy = .false.
 if (present(hierarchy)) use_hierarchy = hierarchy
 if (id==master) write(*,"(/,a,l1)") '--> checking radiation diffusion on sine function'//&
                                      get_tag(implicit_radiation)//', direct calls to solver=',use_hierarchy

 call init_part()
 call setup_radiation_diffusion_problem_sinusoid(kappa_code,c_code,xi0,rho0,pmassi)

 tmax = 5.e-22
 dtmax = tmax
 implicit_radiation_store_drad = .true.
 limit_radiation_flux = .false.
 do i = 1,2
    call get_derivs_global(dt_new=dtnew,dt=1.e-25)
 enddo

 nerr_e = 0
 ncheck_e = 0
 errmax_e = 0.
 tol_e = 1e-10

 nerr_f = 0
 ncheck_f = 0
 errmax_f = 0.
 tol_f = 2e-2

 l0 = 2*pi/(xmax-xmin)
 do i=1,npart
    rhoi = rho(i)
    D0  = c_code*(1./3)/kappa_code/rhoi
    exact_grE  =  xi0*rho0*0.1*l0   *cos(xyzh(1,i)*l0)
    exact_DgrF = -xi0*D0  *0.1*l0*l0*sin(xyzh(1,i)*l0)
    call checkvalbuf(radprop(ifluxx,i),exact_grE,tol_e, '  grad{E}',nerr_e(1),ncheck_e,errmax_e)
    !  this test only works with fixed lambda = 1/3
    if (.not.limit_radiation_flux) then
       call checkvalbuf(drad(iradxi,i),exact_DgrF,tol_f,'dxi/dt = -div F = div(D grad E)',nerr_f(1),ncheck_f,errmax_f)
    endif
 enddo
 call checkvalbuf_end('  grad{E}',ncheck_e,nerr_e(1),errmax_e,tol_e)
 call checkvalbuf_end('D*grad{F}',ncheck_f,nerr_f(1),errmax_f,tol_f)
 call update_test_scores(ntests,nerr_e,npass)
 call update_test_scores(ntests,nerr_f,npass)
 !
 ! check that energy is conserved (i.e. dEtot/dt = 0)
 !
 de = 0.; degas = 0.; derad = 0.; dekin = 0.
 do i=1,npart
    dekin = dekin + dot_product(vxyzu(1:3,i),fxyzu(1:3,i))  ! v.dv/dt
    degas = degas + fxyzu(4,i)       ! du/dt
    derad = derad + drad(iradxi,i)   ! dxi/dt
 enddo
 de = pmassi*(dekin + degas + derad)/xi0
 de = reduceall_mpi('+',de)
 call checkval(de,0.,4.e-9,nerr_e(1),'dE/dt = 0')
 call update_test_scores(ntests,nerr_e,npass)
 !
 ! now solve diffusion as a function of time
 !
 implicit_radiation_store_drad = .false.
 if (use_hierarchy) then
    call evolve_diffusion_hierarchy(ntests,npass,kappa_code,c_code,xi0)
 else
    t  = 0.
    if (write_output) then
       write (filename,'(a5,i3.3)') 'rad_test_',0
       call write_fulldump(t,filename)
    endif
    dt = dtnew
    dtext = dt
    if (implicit_radiation) dt = 6.7e-24 ! force the explicit timestep

    call init_step(npart,t,dtmax)
    i = 0
    D0  = c_code*(1./3)/kappa_code/rho0
    if (write_output) print "(/,a,3(es10.3,a))", ' exact solution: ',xi0,'*(1 + 0.1*sin(x*',l0,')*exp(',-l0*l0*D0,'*t))'
    do while(t < tmax)
       dtext = dt
       call step(npart,nactive,t,dt,dtext,dtnew)  ! step from t to t+dt
       t = t + dt
       dt = dtnew
       if (implicit_radiation) dt = 6.7e-24 ! force the explicit timestep

       i = i + 1

       if (mod(i,10) == 0) then
          if (implicit_radiation) then
             tol_xi = 5.e-3
          else
             tol_xi = 3.5e-4
          endif
          write (string,"(a,i3.3,a)") 'xi(t_', i, ')'
          call check_xi_exact(ntests,npass,t,kappa_code,c_code,xi0,tol_xi,trim(string))
       endif
       if (write_output) then
          write (filename,'(a5,i3.3)') 'rad_test_', i
          call write_fulldump(t,filename)
       endif
    enddo
 endif

 ! reset various things
 call init_part()
 limit_radiation_flux = .true.
 drad = 0.

end subroutine test_radiation_diffusion

!---------------------------------------------------------
!+
!  evolve the diffusion test by calling the implicit solver
!  directly. With individual timesteps, particles with x <= 0.1
!  are only active every second substep, so the solver sees
!  inactive particles in the tree and as neighbours
!+
!---------------------------------------------------------
subroutine evolve_diffusion_hierarchy(ntests,npass,kappa_code,c_code,xi0)
 use dim,                only:maxphase,maxp,ind_timesteps
 use part,               only:npart,xyzh,vxyzu,rad,radprop,drad,iphase,isetphase,igas,itlastrad
 use neighkdtree,        only:build_tree
 use mpiutils,           only:reduceall_mpi
 use radiation_implicit, only:do_radiation_implicit
 integer, intent(inout) :: ntests,npass
 real,    intent(in)    :: kappa_code,c_code,xi0
 integer, parameter :: nsub = 64     ! number of substeps (even, so all particles finish together)
 real,    parameter :: dt = 6.7e-24  ! shortest timestep
 real    :: t,etot0,etot1
 integer :: i,k,nactive,ierr,nfail,nerr(1)
 logical :: active

 radprop(itlastrad,1:npart) = 0.
 etot0 = reduceall_mpi('+',sum(vxyzu(4,1:npart) + rad(iradxi,1:npart)))

 nfail = 0
 do k = 1,nsub
    t = k*dt
    !--0.1 is not a tree split, so some leaves hold both active and inactive particles
    nactive = 0
    do i = 1,npart
       active = (.not.ind_timesteps .or. xyzh(1,i) > 0.1 .or. mod(k,2)==0)
       if (maxphase==maxp) iphase(i) = isetphase(igas,iactive=active)
       if (active) nactive = nactive + 1
    enddo
    call build_tree(npart,nactive,xyzh,vxyzu)
    call do_radiation_implicit(dt,t,npart,rad,xyzh,vxyzu,radprop,drad,ierr)
    if (ierr /= 0) nfail = nfail + 1
 enddo
 nerr = 0
 call checkval(nfail,0,0,nerr(1),'implicit solve succeeded on every substep')
 call update_test_scores(ntests,nerr,npass)

 !--twice the measured error (about 1e-3), as for the global-timestep test
 call check_xi_exact(ntests,npass,t,kappa_code,c_code,xi0,2.e-3,'xi(t)')

 !--total energy (gas plus radiation, equal masses), to twice the measured error. With
 !  individual timesteps the exchange with inactive neighbours is not symmetric, so energy
 !  is conserved only to first order in dt (about 3e-4 here, against 5e-6 without)
 etot1 = reduceall_mpi('+',sum(vxyzu(4,1:npart) + rad(iradxi,1:npart)))
 nerr = 0
 call checkval(etot1,etot0,merge(1.e-3,1.e-5,ind_timesteps),nerr(1),'total energy conserved')
 call update_test_scores(ntests,nerr,npass)

end subroutine evolve_diffusion_hierarchy

!---------------------------------------------------------
!+
!  check xi against the exact decay of the sine wave at time t
!+
!---------------------------------------------------------
subroutine check_xi_exact(ntests,npass,t,kappa_code,c_code,xi0,tol_xi,label)
 use part,     only:npart,xyzh,rho,rad
 use boundary, only:xmin,xmax
 use physcon,  only:pi
 integer,          intent(inout) :: ntests,npass
 real,             intent(in)    :: t,kappa_code,c_code,xi0,tol_xi
 character(len=*), intent(in)    :: label
 real    :: l0,D0,exact_xi,errmax_xi
 integer :: j,nerr_xi(1),ncheck_xi

 nerr_xi = 0
 ncheck_xi = 0
 errmax_xi = 0.
 l0 = 2.*pi/(xmax-xmin)
 do j = 1,npart
    D0 = c_code*(1./3)/kappa_code/rho(j)
    exact_xi = xi0*(1.+0.1*sin(xyzh(1,j)*l0)*exp(-l0*l0*t*D0))
    call checkvalbuf(rad(iradxi,j),exact_xi,tol_xi,label,nerr_xi(1),ncheck_xi,errmax_xi)
 enddo
 call checkvalbuf_end(label,ncheck_xi,nerr_xi(1),errmax_xi,tol_xi)
 call update_test_scores(ntests,nerr_xi,npass)

end subroutine check_xi_exact

!---------------------------------------------------------
!+
!  subroutine to setup sinusoidal diffusion problem
!  (split so this can be called by different subtests)
!+
!---------------------------------------------------------
subroutine setup_radiation_diffusion_problem_sinusoid(kappa_code,c_code,xi0,rho0,pmassi)
 use dim,             only:maxp
 use io,              only:id,master
 use eos,             only:gamma,gmw,iopacity_type
 use apr,             only:sync_aprmassoftype
 use part,            only:npart,hfact,xyzh,vxyzu,massoftype,igas,periodic,&
                           iphase,maxphase,isetphase,npartoftype,&
                           rad,radprop,ifluxx,maxvxyzu
 use units,           only:unit_opacity,get_c_code,unit_ergg,get_radconst_code
 use physcon,         only:Rg,pi,seconds
 use mpidomain,       only:i_belong
 use mpiutils,        only:reduceall_mpi
 use boundary,        only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use kernel,          only:hfact_default
 use unifdis,         only:set_unifdis
 use radiation_utils, only:kappa_cgs
 real, intent(out) :: kappa_code,c_code,xi0,rho0,pmassi
 real :: psep,a,cv1,Tref
 integer :: i,nactive
 integer(kind=8) :: nptot

 psep = 1./32.
 hfact = hfact_default
 npart = 0
 call set_boundary(-0.5,0.5,-0.1,0.1,-0.1,0.1)
 call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep,&
                  hfact,npart,xyzh,periodic,mask=i_belong)
 nptot = reduceall_mpi('+',npart)

 rho0 = 2.5e-24
 massoftype(igas) = rho0*dxbound*dybound*dzbound/nptot
 call sync_aprmassoftype
 pmassi = massoftype(igas)
 if (maxphase==maxp) iphase(1:npart) = isetphase(igas,iactive=.true.)
 npartoftype(:) = 0
 npartoftype(igas) = npart
 nactive = npart

 iopacity_type = 2  ! constant opacity value
 c_code = get_c_code()
 gamma = 5./3.
 gmw = 2.0
 cv1 = (gamma-1.)*gmw/Rg*unit_ergg
 a   = get_radconst_code()
 pmassi = massoftype(igas)
 kappa_cgs = 1.0
 kappa_code = kappa_cgs/unit_opacity
 Tref = 100.

 xi0 = a*Tref**4/rho0
 do i=1,npart
    vxyzu(4,i) = (Tref/cv1)
    radprop(ikappa,i) = kappa_code
    rad(iradxi,i)     = xi0*(1. + 1e-1*sin(xyzh(1,i)*2.*pi/(xmax-xmin)))
    radprop(ithick,i) = 1.
    !etot = vxyzu(4,i) + rad(iradxi,i)
    !Tgas = vxyzu(4,i)*cv1
    !rhoi = rhoh(xyzh(4,i),pmassi)
    !Trad = (rho0*(etot-vxyzu(4,i))/a)**(1./4.)
    !print*, Tref, Trad, Tgas
 enddo

end subroutine setup_radiation_diffusion_problem_sinusoid
!---------------------------------------------------------
!+
!  function to return 'implicit' or 'explicit' string
!  based on logical flag
!+
!---------------------------------------------------------
function get_tag(implicit_flag) result(str)
 logical, intent(in) :: implicit_flag
 character(len=9) :: str

 if (implicit_flag) then
    str = ' implicit'
 else
    str = ' explicit'
 endif

end function get_tag

end module testradiation
