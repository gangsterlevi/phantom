!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testsetstar
!
! Unit tests of set_star and relax_star modules
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: checksetup, datafiles, dim, eos, io, mpidomain, options,
!   orbits, part, physcon, radiation_utils, readwrite_mesa, relaxstar, setorbit,
!   setstar, setstar_utils, sortutils, table_utils, testutils, units, vectorutils
!
 use testutils, only:checkval,update_test_scores
 implicit none
 public :: test_setstar

 private

contains
!-----------------------------------------------------------------------
!+
!   Unit tests of the set_star routine and associated utilities
!+
!-----------------------------------------------------------------------
subroutine test_setstar(ntests,npass)
 use checksetup, only:check_setup
 use io,         only:id,master
 use dim,        only:gravity,do_radiation
 integer, intent(inout) :: ntests,npass

 if (.not.gravity) then
    if (id==master) write(*,"(/,a)") '--> SKIPPING STAR SETUP TESTS'
    return
 else
    if (id==master) write(*,"(/,a,/)") '--> TESTING STAR SETUP'
 endif

 ! test the get_mass_coordinate routine in relax_star
 call test_get_mass_coord(ntests,npass)
 call test_radrelax_energy_change(ntests,npass)
 call test_constant_entropy_core(ntests,npass)
 call test_stellar_orbital_placement(ntests,npass)

 ! test polytrope
 if (.not. do_radiation) call test_polytrope(ntests,npass)

 ! test red supergiant
 call test_redsupergiant(ntests,npass)

 ! test white dwarf
 call test_whitedwarf(ntests,npass)

end subroutine test_setstar

! Check that gas-radiation exchange cancels locally, while energy moving
! between particles cannot hide evolution. Exclude other bodies and dead gas.
subroutine test_radrelax_energy_change(ntests,npass)
 use relaxstar, only:get_radrelax_energy_change
 use part, only:massoftype,igas,aprmassoftype,apr_level
 use dim, only:use_apr
 integer, intent(inout) :: ntests,npass
 real :: xyzh(4,5),u(5),xi(5),u0(5),xi0(5),change,energy,rmax,mass_prev,energy_ref
 real :: apr_mass_prev(size(aprmassoftype,2))
 integer :: nerr(1)
 integer(kind=1) :: levels_prev(5)

 mass_prev = massoftype(igas)
 massoftype(igas) = 2.
 if (use_apr) then
    apr_mass_prev = aprmassoftype(igas,:)
    levels_prev = apr_level(1:5)
    aprmassoftype(igas,:) = 2.
    apr_level(1:5) = 1
 endif
 xyzh = 0.
 xyzh(4,:) = 1.
 xyzh(1,:) = [1000.,1.,2.,3.,1000.]
 xyzh(4,5) = -1.
 u0 = 10.
 xi0 = 4.
 u = u0
 xi = xi0
 ! Changes on another body and a dead particle must not enter the diagnostic.
 u(1) = 1000.
 u(5) = 1000.
 call get_radrelax_energy_change(1,5,xyzh,u,xi,u0,xi0,change,energy_ref,rmax)
 call checkval(change,0.,1.e-14,nerr(1),'radiative error excludes other bodies and dead particles')
 call update_test_scores(ntests,nerr,npass)
 call checkval(energy_ref,84.,1.e-14,nerr(1),'radiative reference includes active gas and photons with mass weights')
 call update_test_scores(ntests,nerr,npass)
 call checkval(rmax,3.,1.e-14,nerr(1),'radiative radius excludes other bodies and dead particles')
 call update_test_scores(ntests,nerr,npass)

 u(2) = u0(2)-2.
 xi(2) = xi0(2)+2.
 call get_radrelax_energy_change(1,5,xyzh,u,xi,u0,xi0,change,energy,rmax)
 call checkval(change,0.,1.e-14,nerr(1),'local gas-photon exchange cancels in combined energy error')
 call update_test_scores(ntests,nerr,npass)

 u(3) = u0(3)+1.
 u(4) = u0(4)-1.
 call get_radrelax_energy_change(1,5,xyzh,u,xi,u0,xi0,change,energy,rmax)
 call checkval(change,4.,1.e-14,nerr(1),'energy redistribution cannot cancel between particles')
 call update_test_scores(ntests,nerr,npass)
 call checkval(energy,energy_ref,1.e-14,nerr(1),'redistribution preserves total thermal energy')
 call update_test_scores(ntests,nerr,npass)

 call get_radrelax_energy_change(1,5,xyzh,10.*u,10.*xi,10.*u0,10.*xi0,change,energy,rmax)
 call checkval(change/energy,4./energy_ref,1.e-14,nerr(1),'fractional radiative error independent of energy units')
 call update_test_scores(ntests,nerr,npass)
 massoftype(igas) = mass_prev
 if (use_apr) then
    aprmassoftype(igas,:) = apr_mass_prev
    apr_level(1:5) = levels_prev
 endif

end subroutine test_radrelax_energy_change

! Test the final stellar centres and velocities, not only the Keplerian
! elements: synchronous spin must preserve the body's orbital translation.
subroutine test_stellar_orbital_placement(ntests,npass)
 use setstar, only:star_t,set_defaults_star,centre_stars_for_orbit,shift_stars,ipoly
 use setorbit, only:orbit_t,set_roche_lobe_orbit,set_orbit
 use part, only:massoftype,igas,iphase,set_particle_type,aprmassoftype,apr_level
 use dim, only:use_apr
 use vectorutils, only:cross_product3D
 use orbits, only:get_eccentricity
 integer, intent(inout) :: ntests,npass
 type(star_t) :: stars(2)
 type(orbit_t) :: orbit
 real :: xyz(4,2),vel(4,2),sinks(6,2),vsinks(3,2),orbxyz(6,2),orbvel(3,2)
 real :: mass_prev,omega(3),vspin(3),xcm(3),vcm(3),dx(3),dv(3)
 real :: apr_mass_prev(size(aprmassoftype,2))
 integer(kind=kind(apr_level)) :: levels_prev(2)
 integer :: i,ierr,ns,norb,types(8),nerr(7)
 integer(kind=kind(iphase)) :: phase_prev(2)

 mass_prev = massoftype(igas)
 phase_prev = iphase(1:2)
 massoftype(igas) = 1.
 if (use_apr) then
    apr_mass_prev = aprmassoftype(igas,:)
    levels_prev = apr_level(1:2)
    aprmassoftype(igas,:) = 1.
    apr_level(1:2) = 1
 endif
 do i=1,2
    call set_defaults_star(stars(i))
    call set_particle_type(i,4) ! star 1's temporary tag
 enddo
 stars(1)%iprofile = ipoly
 stars(1)%isinkcore = .true.
 stars(2)%iprofile = 0
 stars(2)%m = '1.0'
 xyz = 0.; vel = 0.; sinks = 0.; vsinks = 0.
 xyz(1,:) = [9.,11.]
 xyz(4,:) = 1.
 vel(1,:) = 0.3; vel(2,:) = 0.4
 sinks(1,1) = 10.5
 sinks(4,1) = 2.
 vsinks(:,1) = [0.3,0.4,0.]
 call centre_stars_for_orbit(2,stars,2,xyz,vel,sinks,vsinks)
 call checkval(stars(1)%r_particles_code,1.25,1.e-13,nerr(1),'donor radius is measured about gas-plus-core centre')
 call set_roche_lobe_orbit(orbit,4.,1.,stars(1)%r_particles_code,10.,ierr)
 norb = 0
 call set_orbit(orbit,4.,1.,0.,0.,orbxyz,orbvel,norb,.false.,ierr)
 types = 0; types(4) = 2
 ns = 1
 call shift_stars(2,stars,orbxyz,orbvel,xyz,vel,sinks,vsinks,2,types,ns,corotate=.true.)
 xcm = (sum(xyz(1:3,:),dim=2)+2.*sinks(1:3,1))/4.
 vcm = (sum(vel(1:3,:),dim=2)+2.*vsinks(:,1))/4.
 call checkval(maxval(abs(xcm-orbxyz(1:3,1))),0.,1.e-12,nerr(2),'stellar centre follows requested orbit')
 call checkval(maxval(abs(vcm-orbvel(:,1))),0.,1.e-12,nerr(3),'synchronous spin preserves orbital velocity')
 call cross_product3D(orbxyz(1:3,1),orbvel(:,1),omega)
 omega = omega/dot_product(orbxyz(1:3,1),orbxyz(1:3,1))
 call cross_product3D(omega,xyz(1:3,1)-xcm,vspin)
 call checkval(maxval(abs(vel(1:3,1)-vcm-vspin)),0.,1.e-12,nerr(4),'stellar spin has prograde omega cross r sign')
 call checkval(maxval(abs(sinks(1:3,1)-orbxyz(1:3,1)-[0.25,0.,0.])),0.,1.e-12,nerr(5), &
               'orbital placement preserves internal core offset')
 dx = sinks(1:3,2)-xcm
 dv = vsinks(:,2)-vcm
 call checkval(get_eccentricity(5.,dx,dv),0.,1.e-11,nerr(6),'final stellar centres have a circular orbit')
 call checkval(ns,2,0,nerr(7),'point companion is created alongside the donor core')
 call update_test_scores(ntests,nerr,npass)
 massoftype(igas) = mass_prev
 iphase(1:2) = phase_prev
 if (use_apr) then
    aprmassoftype(igas,:) = apr_mass_prev
    apr_level(1:2) = levels_prev
 endif

end subroutine test_stellar_orbital_placement

! Check the shooting solve from different guesses, including zero (which
! multiplicative mass updates cannot escape), and closure of the output mass.
subroutine test_constant_entropy_core(ntests,npass)
 use datafiles, only:find_phantom_datafile
 use readwrite_mesa, only:read_mesa
 use setfixedentropycore, only:set_fixedS_softened_core
 use table_utils, only:interpolator
 use physcon, only:solarr,solarm
 use eos, only:entropy,get_mean_molecular_weight
 integer, intent(inout) :: ntests,npass
 real, allocatable :: rho(:),r(:),pres(:),m(:),en(:),temp(:),X(:),Y(:),mu(:)
 real :: mt,mc,rc,mh,mc_ref,Sc,muc,errS
 integer :: itest,icore,i,ierr,nerr(1)
 character(len=500) :: filename

 filename = find_phantom_datafile('P12_Phantom_Profile.data','star_data_files')
 mc_ref = 0.
 do itest=1,4
    call read_mesa(filename,rho,r,pres,m,en,temp,0.7,0.02,X,Y,mu,mt,ierr,cgsunits=.true.)
    call checkval(ierr,0,0,nerr(1),'read constant-entropy core test profile')
    call update_test_scores(ntests,nerr,npass)
    if (ierr/=0) return
    rc = 0.1*r(size(r))/solarr
    call interpolator(r,rc*solarr,icore)
    ! Strongly nonuniform central spacing exposes cancellation/alternating
    ! errors in a pressure recurrence based on neighbouring absolute values.
    if (itest==4) then
       do i=1,icore-1
          r(i) = r(icore)*(real(i)/real(icore))**2
       enddo
       call interpolator(r,rc*solarr,icore)
    endif
    mh = m(icore)
    mc = 0.5*mh/solarm
    if (itest==2) mc = 0.
    if (itest==3) mc = 0.9*mh/solarm
    muc = get_mean_molecular_weight(X(icore),1.-X(icore)-Y(icore))
    Sc = entropy(rho(icore),pres(icore),muc,2)
    call set_fixedS_softened_core(12,mc,rc,rho,r,pres,m,X(icore),Y(icore),ierr)
    call checkval(ierr,0,0,nerr(1),'constant-entropy core converges')
    call update_test_scores(ntests,nerr,npass)
    if (ierr/=0) return
    call checkval((m(icore)-mh)/(mh-mc*solarm),0.,1.e-10,nerr(1),'softened gas mass closes at boundary')
    call update_test_scores(ntests,nerr,npass)
    if (itest==1) mc_ref = mc
    if (itest<=3) then
       call checkval(mc,mc_ref,1.e-10,nerr(1),'core mass independent of initial guess')
       call update_test_scores(ntests,nerr,npass)
    endif
    call checkval(all(pres(1:icore-1)>=pres(2:icore)),.true.,nerr(1),'core pressure increases inward')
    call update_test_scores(ntests,nerr,npass)
    errS = 0.
    do i=1,icore
       errS = max(errS,abs((entropy(rho(i),pres(i),muc,2)-Sc)/Sc))
    enddo
    call checkval(errS,0.,1.e-8,nerr(1),'softened core preserves boundary entropy')
    call update_test_scores(ntests,nerr,npass)
    deallocate(rho,r,pres,m,en,temp,X,Y,mu)
 enddo
end subroutine test_constant_entropy_core

!-----------------------------------------------------------------------
!+
!   Check that the procedure to get the Lagrangian mass coordinate
!   from the radius works
!+
!-----------------------------------------------------------------------
subroutine test_get_mass_coord(ntests,npass)
 use sortutils,     only:r2func,find_rank
 use part,          only:xyzh,massoftype,igas,iorder=>ll
 use setstar_utils, only:get_mass_coord
 use table_utils,   only:linspace
 use testutils,     only:checkvalbuf,checkvalbuf_end
 integer, intent(inout) :: ntests,npass
 real,    allocatable   :: mass_enclosed_r(:)
 integer :: i,i1,itest,np,nfail(1),ncheck
 real    :: massri,errmax,x0(3)
 real, parameter :: tol = 1.e-14

 i1 = 10  ! start at non-zero offset to check this
 np = 100
 x0 = 0.

 ! place the particles in a line between x=[0,1]
 xyzh(:,1:np) = 0.
 call linspace(xyzh(1,i1+1:np),0.,1.)

 ! give them equal masses
 massoftype(igas) = 3.e-6

 ! call the routine we are trying to test
 call get_mass_coord(i1,np,xyzh,mass_enclosed_r,x0)

 ! check memory was allocated correctly
 call checkval(size(mass_enclosed_r),np-i1,0,nfail(1),'size of mass_enclosed_r array')
 call update_test_scores(ntests,nfail,npass)

 ! check that the mass ranking is correct
 ncheck = 0; nfail = 0; errmax = 0.
 do i=i1+1,np
    call checkvalbuf(mass_enclosed_r(i-i1),(i-i1-1)*massoftype(igas),&
                     tol,'m(x) for particles in x=[0,1]',nfail(1),ncheck,errmax)
 enddo
 call checkvalbuf_end('m(x) for particles in x=[0,1]',ncheck,nfail(1),errmax,tol)
 call update_test_scores(ntests,nfail,npass)

 do itest=1,2
    if (itest==2) then
       ! now check pathological case where half the particles are at y=0, and half are at y=1
       xyzh(:,1:np) = 0.
       xyzh(2,i1+np/2:np) = 1.
       call get_mass_coord(i1,np,xyzh,mass_enclosed_r,x0)
    endif
    !
    ! check that this agrees with our previous method that
    ! works only for equal mass particles
    !
    mass_enclosed_r = mass_enclosed_r / ((np-i1)*massoftype(igas))  ! make it a mass fraction
    call find_rank(np-i1,r2func,xyzh(:,i1+1:np),iorder)

    ncheck = 0; nfail = 0; errmax = 0.
    do i=i1+1,np
       massri = real(iorder(i-i1)-1) / real(np-i1)
       call checkvalbuf(mass_enclosed_r(i-i1),massri,tol,'m(<r) agrees with previous method',nfail(1),ncheck,errmax)
    enddo
    call checkvalbuf_end('m(<r) same as previous method',ncheck,nfail(1),errmax,tol)
    call update_test_scores(ntests,nfail,npass)
 enddo

end subroutine test_get_mass_coord

!-----------------------------------------------------------------------
!+
!   test that we can successfully relax a polytrope
!   and that the profile matches the analytic solution
!+
!-----------------------------------------------------------------------
subroutine test_polytrope(ntests,npass)
 use io,        only:id,master,iverbose
 use part,      only:npart,npartoftype,xyzh,vxyzu,eos_vars,rad,massoftype,hfact,&
                     xyzmh_ptmass,vxyz_ptmass,nptmass
 use mpidomain, only:i_belong
 use options,   only:ieos
 use physcon,   only:solarr,solarm,pi
 use eos,       only:gamma,X_in,Z_in,polyk
 use setstar,   only:star_t,set_star,set_defaults_star,ipoly
 use units,     only:set_units
 use checksetup, only:check_setup
 integer, intent(inout) :: ntests,npass
 type(star_t) :: star
 real :: rhozero,rmserr,ekin,x0(3)
 integer(kind=8) :: ntot
 integer :: ierr,nfail(1),i,nerror,nwarn

 npart = 0
 npartoftype = 0
 massoftype = 0.
 iverbose = 0
 call set_units(dist=10.*solarr,mass=10.*solarm,G=1.d0)
 ieos = 2
 gamma = 5./3.
 call set_defaults_star(star)
 star%m = '1*msun'
 star%r = '1*rsun'
 star%iprofile = ipoly  ! a polytrope
 star%np = 1000
 x0 = 0.
 ! do this test twice, to check the second star relaxes...
 do i=1,2
    if (i==2) x0 = [3.,0.,0.]

    call set_star(id,master,star,xyzh,vxyzu,eos_vars,rad,&
               npart,npartoftype,massoftype,hfact,&
               xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,gamma,X_in,Z_in,&
               relax=.true.,use_var_comp=.false.,write_rho_to_file=.false.,&
               rhozero=rhozero,npart_total=ntot,mask=i_belong,ierr=ierr,&
               write_files=.false.,density_error=rmserr,energy_error=ekin,x0=x0)

    call checkval(ierr,0,0,nfail(1),'set_star runs with ierr = 0')
    call update_test_scores(ntests,nfail,npass)

    call check_setup(nerror,nwarn,restart=.false.)
    call checkval(nerror+nwarn,0,0,nfail(1),'no errors or warnings')
    call update_test_scores(ntests,nfail,npass)

    call checkval(rhozero,100./(4./3.*pi),1e-6,nfail(1),'mean density')
    call update_test_scores(ntests,nfail,npass)

    call checkval(star%polyk,1.9694457e-2,1e-6,nfail(1),'polyk value for M=0.1,R=0.1')
    call update_test_scores(ntests,nfail,npass)

    call checkval(polyk,1.9694457e-2,1e-6,nfail(1),'polyk value for M=0.1,R=0.1')
    call update_test_scores(ntests,nfail,npass)

    call checkval(rmserr,0.0,0.04,nfail(1),'error in density profile')
    call update_test_scores(ntests,nfail,npass)

    call checkval(ekin,0.,1e-7,nfail(1),'ekin/epot < 1.e-7')
    call update_test_scores(ntests,nfail,npass)
 enddo

 call checkval(npart,2*star%np,0,nfail(1),'np = 2000')
 call update_test_scores(ntests,nfail,npass)

end subroutine test_polytrope

!-----------------------------------------------------------------------
!+
!   test that we can successfully relax a red supergiant star from
!   Lau et al. (2022a,b)
!+
!-----------------------------------------------------------------------
subroutine test_redsupergiant(ntests,npass)
 use io,        only:id,master,iverbose,warning
 use dim,       only:do_radiation
 use datafiles, only:find_phantom_datafile
 use part,      only:init_part,npart,npartoftype,xyzh,vxyzu,eos_vars,rad,massoftype,hfact,&
                     xyzmh_ptmass,vxyz_ptmass,nptmass,rho,igas,igasP,imu,iX,iZ,iradxi,&
                     eos_vars,itemp
 use mpidomain, only:i_belong
 use options,   only:ieos
 use physcon,   only:solarr,solarm
 use eos,       only:init_eos,gamma,gmw,X_in,Z_in,irecomb,eos_works_with_radiation
 use setstar,   only:star_t,set_star,set_defaults_star,imesa,ierr_radiation_conflict
 use units,     only:set_units
 use relaxstar, only:maxits
 use checksetup,  only:check_setup
 use testutils,   only:checkvalbuf,checkvalbuf_end
 use table_utils, only:yinterp
 use readwrite_mesa,  only:read_mesa
 use radiation_utils, only:Trad_from_radxi
 integer, intent(inout) :: ntests,npass
 type(star_t) :: star
 character(len=500) :: filepath
 real :: rhozero,rmserr,rmserr_mu,rmserr_X,rmserr_Z,ekin,x0(3),errmax(2)
 real :: Mstar,tolprof,tolcomposition,rhoj,rhoj_mesa,rj,Tgas,Trad
 integer(kind=8) :: ntot
 integer :: ierr,nfail(2),ncheck(2),i,j,nerror,nwarn,iunit,expected_error
 logical :: relax_star,var_comp
 real, allocatable :: r(:),den(:),pres(:),temp(:),en(:),mtab(:),Xfrac(:),Yfrac(:),Zfrac(:),mu(:)

 filepath = find_phantom_datafile('lau22.data','star_data_files')
 open(newunit=iunit,file=filepath,status="old",iostat=ierr)
 if (ierr /= 0) then
    call warning('test_setstar','lau22.data not found, skipping red supergiant test')
    return
 else
    close(unit=iunit)
 endif

 iverbose = 0
 call set_units(dist=solarr,mass=solarm,G=1.d0)
 gamma = 5./3.
 gmw = 0.61821
 X_in = 0.698
 Z_in = 1.-X_in-0.287
 irecomb = 0
 relax_star = .true.
 maxits = 100
 var_comp = .false.
 call set_defaults_star(star)
 star%isinkcore   = .true.
 star%hsoft          = '9.25'
 star%hacc           = '0.0'
 star%mcore          = '3.8405'
 star%np             = 10000
 star%input_profile  = trim(filepath)
 star%iprofile = imesa
 x0 = 0.

 nfail = 0; ncheck = 0; errmax = 0.
 tolprof = 0.08
 tolcomposition = 0.08

 call read_mesa(filepath,den,r,pres,mtab,en,temp,X_in,Z_in,Xfrac,Yfrac,mu,Mstar,ierr)
 Zfrac = 1.-Xfrac-Yfrac

 do i=1,4
    call init_part()
    if (i==1) then  ! ideal gas EoS
       ieos = 2
    elseif (i==2) then  ! ideal gas + radiation EoS
       ieos = 12
    elseif (i==3) then  ! MESA EoS
       ieos = 10
    elseif (i==4) then  ! Gas + radiation + recombination energy EoS
       ieos = 20
    elseif (i==5) then  ! ideal gas + radiation EoS, with variable composition
       ieos = 12
       var_comp = .true.
    endif

    if (id==master) write(*,"(/,a,i2/)") '--> Testing ieos = ',ieos
    call init_eos(ieos,ierr)

    rmserr = 0.
    call set_star(id,master,star,xyzh,vxyzu,eos_vars,rad,&
               npart,npartoftype,massoftype,hfact,&
               xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,gamma,X_in,Z_in,&
               relax=relax_star,use_var_comp=var_comp,write_rho_to_file=.false.,&
               rhozero=rhozero,npart_total=ntot,mask=i_belong,ierr=ierr,&
               write_files=.false.,density_error=rmserr,energy_error=ekin)

    if (do_radiation .and. (.not. eos_works_with_radiation(ieos))) then
       expected_error = ierr_radiation_conflict
    else
       expected_error = 0
    endif
    call checkval(ierr,expected_error,0,nfail(1),'set_star runs with expected ierr')
    call update_test_scores(ntests,nfail,npass)
    if (ierr == ierr_radiation_conflict) cycle  ! skip the rest of the checks if we get a radiation conflict error

    call check_setup(nerror,nwarn,restart=.false.)
    call checkval(nerror+nwarn,0,0,nfail(1),'no errors or warnings')
    call update_test_scores(ntests,nfail,npass)

    call checkval(rmserr,0.0,0.04,nfail(1),'error in density profile')
    call update_test_scores(ntests,nfail,npass)

    call checkval(ekin,0.,5e-6,nfail(1),'ekin/epot < 1.e-7')
    call update_test_scores(ntests,nfail,npass)

    rmserr = 0.
    rmserr_mu = 0.
    rmserr_X = 0.
    rmserr_Z = 0.
    do j=1,npart
       rj = sqrt(dot_product(xyzh(1:3,j),xyzh(1:3,j)))
       rhoj = rho(j)
       rhoj_mesa = yinterp(den,r,rj)
       rmserr = rmserr + (1-rhoj/rhoj_mesa)**2

       if (do_radiation) then
          Trad = Trad_from_radxi(rhoj, rad(iradxi,j))
          Tgas = eos_vars(itemp,j)
          call checkvalbuf(Trad/Tgas,1.,1e-15,'Trad/Tgas = 1 for radiative star',nfail(2),ncheck(2),errmax(2),use_rel_tol=.true.)
       endif

       if (i==5) then  ! calculate rms error in mu, X, Z, not implemented yet
          rmserr_mu = rmserr_mu + (1-abs(eos_vars(imu,j)/yinterp(mu,r,rj)))**2
          rmserr_X = rmserr_X + (1-abs(eos_vars(iX,j)/yinterp(Xfrac,r,rj)))**2
          rmserr_Z = rmserr_Z + (1-abs(eos_vars(iZ,j)/yinterp(Zfrac,r,rj)))**2
       endif
    enddo
    rmserr = sqrt(rmserr/npart)  ! root mean square relative error
    call checkval(rmserr,0.,tolprof,nfail(1),'density matches MESA profile')
    call update_test_scores(ntests,nfail,npass)

    call checkvalbuf_end('Trad/Tgas ~ 1 for radiative star',ncheck(2),nfail(2),errmax(2),1.e-15)
    call update_test_scores(ntests,nfail,npass)

    if (i==5) then  ! not implemented yet
       rmserr_mu = sqrt(rmserr_mu/npart)
       rmserr_X = sqrt(rmserr_X/npart)
       rmserr_Z = sqrt(rmserr_Z/npart)
       call checkval(rmserr_mu,0.,tolcomposition,nfail(1),'mu matches MESA profile')
       call update_test_scores(ntests,nfail,npass)
       call checkval(rmserr_X,0.,tolcomposition,nfail(1),'X matches MESA profile')
       call update_test_scores(ntests,nfail,npass)
       call checkval(rmserr_Z,0.,tolcomposition,nfail(1),'Z matches MESA profile')
       call update_test_scores(ntests,nfail,npass)
    endif

 enddo

end subroutine test_redsupergiant

!-----------------------------------------------------------------------
!+
!   test that we can successfully relax a white dwarf with the
!   zero temperature EoS and Helmholtz eos(Ali Pourmand, similar to test_redsupergiant)
!+
!-----------------------------------------------------------------------
subroutine test_whitedwarf(ntests,npass)
 use io,        only:id,master,iverbose,warning
 use dim,       only:do_radiation
 use datafiles, only:find_phantom_datafile
 use part,      only:init_part,npart,npartoftype,xyzh,vxyzu,eos_vars,rad,massoftype,hfact,&
                     xyzmh_ptmass,vxyz_ptmass,nptmass,rho,igas,&
                     eos_vars
 use mpidomain, only:i_belong
 use options,   only:ieos
 use physcon,   only:solarr,solarm
 use eos,       only:init_eos,gamma,X_in,Z_in
 use setstar,   only:star_t,set_star,set_defaults_star,imesa
 use units,     only:set_units
 use relaxstar, only:maxits
 use checksetup,  only:check_setup
 use testutils,   only:checkvalbuf,checkvalbuf_end
 use table_utils, only:yinterp
 use readwrite_mesa,  only:read_mesa
 integer, intent(inout) :: ntests,npass
 type(star_t) :: star
 character(len=500) :: filepath
 real :: rhozero,rmserr,ekin,x0(3),errmax(2)
 real :: Mstar,tolprof,rhoj,rhoj_mesa,rj
 integer(kind=8) :: ntot
 integer :: ierr,nfail(2),ncheck(2),i,j,nerror,nwarn,iunit,expected_error
 logical :: relax_star,var_comp
 real, allocatable :: r(:),den(:),pres(:),temp(:),en(:),mtab(:),Xfrac(:),Yfrac(:),Zfrac(:),mu(:)

 filepath = find_phantom_datafile('pourmandwd.data','star_data_files')
 open(newunit=iunit,file=filepath,status="old",iostat=ierr)
 if (ierr /= 0) then
    call warning('test_setstar','pourmandwd.data not found, skipping white dwarf test')
    return
 else
    close(unit=iunit)
 endif

 iverbose = 0
 call set_units(dist=solarr,mass=solarm,G=1.d0)
 gamma = 5./3.
 X_in = 0.
 Z_in = 1.
 relax_star = .true.
 maxits = 1500
 var_comp = .false.
 call set_defaults_star(star)
 star%np             = 10000
 star%input_profile  = trim(filepath)
 star%iprofile = imesa
 x0 = 0.

 nfail = 0; ncheck = 0; errmax = 0.
 tolprof = 0.1

 call read_mesa(filepath,den,r,pres,mtab,en,temp,X_in,Z_in,Xfrac,Yfrac,mu,Mstar,ierr)
 Zfrac = 1.-Xfrac-Yfrac

 do i=1,2
    call init_part()
    if (i==1) then  ! EOS Helmholtz
       ieos = 15
    elseif (i==2) then  ! zero Temperature EOS
       ieos = 25
    endif

    if (id==master) write(*,"(/,a,i2/)") '--> Testing ieos = ',ieos
    call init_eos(ieos,ierr)

    rmserr = 0.
    call set_star(id,master,star,xyzh,vxyzu,eos_vars,rad,&
               npart,npartoftype,massoftype,hfact,&
               xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,gamma,X_in,Z_in,&
               relax=relax_star,use_var_comp=var_comp,write_rho_to_file=.false.,&
               rhozero=rhozero,npart_total=ntot,mask=i_belong,ierr=ierr,&
               write_files=.false.,density_error=rmserr,energy_error=ekin)

    call checkval(ierr,expected_error,0,nfail(1),'set_star runs with expected ierr')
    call update_test_scores(ntests,nfail,npass)

    call check_setup(nerror,nwarn,restart=.false.)
    call checkval(nerror+nwarn,0,0,nfail(1),'no errors or warnings')
    call update_test_scores(ntests,nfail,npass)

    call checkval(rmserr,0.0,0.05,nfail(1),'error in density profile')
    call update_test_scores(ntests,nfail,npass)

    call checkval(ekin,0.,5e-6,nfail(1),'ekin/epot < 1.e-7')
    call update_test_scores(ntests,nfail,npass)

    rmserr = 0.
    do j=1,npart
       rj = sqrt(dot_product(xyzh(1:3,j),xyzh(1:3,j)))
       rhoj = rho(j)
       rhoj_mesa = yinterp(den,r,rj)
       rmserr = rmserr + (1-rhoj/rhoj_mesa)**2

    enddo
    rmserr = sqrt(rmserr/npart)  ! root mean square relative error
    call checkval(rmserr,0.,tolprof,nfail(1),'density matches MESA profile')
    call update_test_scores(ntests,nfail,npass)

 enddo

end subroutine test_whitedwarf

end module testsetstar
