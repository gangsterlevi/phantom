!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup of two stars or sink particles in a binary
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - corotate   : *set stars in corotation*
!   - deltat     : *output interval as fraction of binary period*
!   - gwinspiral : *set gravitational wave inspiral*
!   - norbits    : *maximum number of binary orbits*
!   - itype      : *orbital input type; 5 derives a circular two-body Roche-lobe orbit*
!   - roche_overflow1 : *percentage overflow, 100*(R1/RL1-1), using the relaxed particle radius*
!
! :Dependencies: centreofmass, dim, eos, extern_gwinspiral, externalforces,
!   infile_utils, io, kernel, mpidomain, options, part, physcon,
!   sethier_utils, sethierarchical, setorbit, setstar, setunits,
!   setup_params, timestep, units
!
 use setstar,       only:star_t
 use setorbit,      only:orbit_t,iorbit_roche
 use sethier_utils, only:max_hier_levels
 use dim,      only:gr
 use eos,      only:ieos
 implicit none
 public :: setpart

 logical :: relax,write_rho_to_file,corotate,gwinspiral
 integer, parameter :: max_stars = max_hier_levels
 type(star_t)  :: star(max_stars)
 type(orbit_t) :: orbit
 real :: norbits,deltat
 integer :: nstar

 private

contains

!----------------------------------------------------------------
!+
!  setup for binary star simulations (with or without gas)
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,&
                   polyk,gamma,hfact,time,fileprefix)
 use part,           only:gr,nptmass,xyzmh_ptmass,vxyz_ptmass,&
                          ihacc,ihsoft,eos_vars,rad,nsinkproperties,iJ2,iReff,ispinx,ispinz
 use setorbit,       only:set_defaults_orbit,set_orbit,set_roche_lobe_orbit
 use sethierarchical,only:set_hierarchical,set_hierarchical_default_options,hs,generate_hierarchy_string
 use options,        only:iexternalforce,alphau
 use externalforces, only:iext_corotate,iext_geopot,iext_star,iext_gwinspiral,omega_corotate,mass1,accradius1
 use extern_gwinspiral, only:Nstar_gw
 use io,             only:master,fatal
 use setstar,        only:set_defaults_stars,set_stars,shift_stars,centre_stars_for_orbit
 use eos,            only:X_in,Z_in,use_var_comp
 use setup_params,   only:rhozero,npart_total
 use mpidomain,      only:i_belong
 use centreofmass,   only:reset_centreofmass
 use setunits,       only:mass_unit,dist_unit
 use physcon,        only:deg_to_rad
 use kernel,         only:hfact_default
 use infile_utils,   only:get_options,infile_exists
 use timestep,       only:tmax,dtmax
 integer,          intent(in)    :: id
 integer,          intent(inout) :: npart
 integer,          intent(out)   :: npartoftype(:)
 real,             intent(out)   :: xyzh(:,:)
 real,             intent(out)   :: massoftype(:)
 real,             intent(out)   :: polyk,gamma,hfact
 real,             intent(inout) :: time
 character(len=*), intent(in)    :: fileprefix
 real,             intent(out)   :: vxyzu(:,:)
 integer :: ierr,nptmass_in,iextern_prev,i
 logical :: add_spin
 real :: xyzmh_ptmass_in(nsinkproperties,max_stars),vxyz_ptmass_in(3,max_stars),angle
 real :: m1,m2,hacc1,hacc2,overflow1
 logical, parameter :: set_oblateness = .false.
!
!--general parameters
!
 dist_unit = 'solarr'
 mass_unit = 'solarm'
 time = 0.
 polyk = 0.
 gamma = 5./3.
 ieos  = 2
 hfact = hfact_default
 norbits = 10.
 deltat = 0.1
!
!--space available for injected gas particles
!  in case only sink particles are used
!
 npart = 0
 npartoftype(:) = 0
 massoftype = 0.

 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 nptmass = 0
 nstar = 2
 call set_defaults_stars(star)
 call set_defaults_orbit(orbit)
 call generate_hierarchy_string(nstar)
 call set_hierarchical_default_options()
 relax = .true.
 corotate = .false.
 gwinspiral = .false.
 use_var_comp = .false.
 write_rho_to_file = .true.

 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",&
   ' Welcome to the Ultimate Binary & Multiple Setup'

 !
 ! read/write from .setup file
 !
 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 !
 !--setup and relax stars as needed
 !
 iextern_prev = iexternalforce
 iexternalforce = 0

 call set_stars(id,master,nstar,star,xyzh,vxyzu,eos_vars,rad,npart,npartoftype,&
                massoftype,hfact,xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,gamma,&
                X_in,Z_in,relax,use_var_comp,write_rho_to_file,&
                rhozero,npart_total,i_belong,ierr)

 add_spin = .false.
 xyzmh_ptmass_in(:,:) = 0.
 vxyz_ptmass_in(:,:) = 0.
 nptmass_in = 0
 if (nstar > 2) then
    ! setup hierarchical multiple system
    do i=1,nstar
       hs%sinks(i)%mass = star(i)%m_code
       hs%sinks(i)%accr = star(i)%hacc_code
    enddo
    call set_hierarchical(fileprefix,nptmass_in,xyzmh_ptmass_in,vxyz_ptmass_in,ierr)
    if (ierr /= 0) call fatal ('setup_binary','error in call to set_hierarchical')
 elseif (nstar == 2) then
    ! setup binary system
    m1 = star(1)%m_code
    m2 = star(2)%m_code
    hacc1 = star(1)%hacc_code
    hacc2 = star(2)%hacc_code

    if (orbit%input_type==iorbit_roche) then
       call centre_stars_for_orbit(nstar,star,npart,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass)
       m1 = star(1)%m_code
       m2 = star(2)%m_code
       overflow1 = orbit%roche_overflow
       call set_roche_lobe_orbit(orbit,m1,m2,star(1)%r_particles_code,overflow1,ierr)
       if (ierr /= 0) call fatal('setup_binary','Roche orbit requires positive masses and a live donor radius')
       if (id==master) then
          print "(/,a)",' CIRCULAR ROCHE-LOBE ORBIT (star 1 is the donor)'
          print "(a,es14.6)",'  relaxed particle radius R1 (code units) = ',star(1)%r_particles_code
          print "(a,es14.6)",'  requested overflow 100*(R1/RL1-1) (%)  = ',overflow1
          print "(a,es14.6)",'  target Roche-lobe radius RL1           = ',star(1)%r_particles_code/(1.+overflow1/100.)
       endif
    endif

    if (iextern_prev==iext_corotate) then
       call set_orbit(orbit,m1,m2,hacc1,hacc2,xyzmh_ptmass_in,vxyz_ptmass_in,&
                      nptmass_in,(id==master),ierr,omega_corotate)
    else
       call set_orbit(orbit,m1,m2,hacc1,hacc2,xyzmh_ptmass_in,vxyz_ptmass_in,&
                      nptmass_in,(id==master),ierr)
       add_spin = corotate
    endif
    if (ierr /= 0) call fatal ('setup_binary','error in call to set_orbit')
 endif
 !
 !--place stars into orbit, or add real sink particles if iprofile=0
 !
 call shift_stars(nstar,star,xyzmh_ptmass_in,vxyz_ptmass_in,xyzh,vxyzu,&
                  xyzmh_ptmass,vxyz_ptmass,npart,npartoftype,nptmass,corotate=add_spin)
 !
 !--restore options
 !
 iexternalforce = iextern_prev

 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 if (nstar >= 2 .and. (iexternalforce==iext_geopot .or. iexternalforce==iext_star) .and. star(1)%iprofile == 0) then
    ! delete first sink particle and copy its properties to the central potential
    nptmass = nptmass - 1
    mass1 = xyzmh_ptmass(4,nptmass+1)
    accradius1 = xyzmh_ptmass(ihacc,nptmass+1)
    xyzmh_ptmass(:,nptmass) = xyzmh_ptmass(:,nptmass+1)
    vxyz_ptmass(:,nptmass) = vxyz_ptmass(:,nptmass+1)
 elseif (nstar == 2 .and. (gwinspiral .or. iexternalforce==iext_gwinspiral)) then
    Nstar_gw(1) = star(1)%np
    Nstar_gw(2) = npart - star(1)%np
    iexternalforce = iext_gwinspiral
 elseif (nstar >= 1 .and. set_oblateness .and. star(1)%iprofile == 0) then
    ! set J2 for sink particle 1 to be equal to oblateness of Saturn
    xyzmh_ptmass(iJ2,1) = 0.01629
    angle = 30.*deg_to_rad
    xyzmh_ptmass(ispinx,1) = sin(angle)
    xyzmh_ptmass(ispinz,1) = cos(angle)
    xyzmh_ptmass(iReff,1) = xyzmh_ptmass(ihacc,1)
 endif

 alphau = min(0.1,alphau)
 if (nstar == 2) then
    ! Roche separation is derived from the newly constructed star, so refresh
    ! its period-based run times even when a previous .in file exists.
    if (orbit%input_type==iorbit_roche .or. .not.infile_exists(fileprefix)) then
       tmax = norbits * orbit%period
       dtmax = deltat * orbit%period
       if (id==master .and. orbit%input_type==iorbit_roche) then
          print "(/,a)",' ROCHE ORBIT SIMULATION TIMES (code units)'
          print "(a,es14.6)",'  orbital period P = ',orbit%period
          print "(a,es14.6,a,f10.4,a)",'  tmax             = ',tmax,' (',norbits,' orbits)'
          print "(a,es14.6,a,f10.4,a)",'  dtmax            = ',dtmax,' (',deltat,' orbits between dumps)'
       endif
    endif
 endif

end subroutine setpart

!----------------------------------------------------------------
!+
!  write options to .setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils,    only:write_inopt
 use setstar,         only:write_options_stars
 use setorbit,        only:write_options_orbit
 use setunits,        only:write_options_units
 use sethierarchical, only:write_hierarchical_setupfile
 character(len=*), intent(in) :: filename
 integer :: iunit

 print "(a)",' writing setup options file '//trim(filename)
 open(newunit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for binary and multiple star setup routines'

 call write_options_units(iunit,gr)
 call write_options_stars(star,relax,write_rho_to_file,ieos,iunit,nstar)
 if (nstar > 2) then
    call write_hierarchical_setupfile(iunit,nstar)
 elseif (nstar == 2) then
    call write_inopt(corotate,'corotate','set stars in corotation',iunit)
    call write_inopt(gwinspiral,'gwinspiral','set gravitational wave inspiral',iunit)
    call write_options_orbit(orbit,iunit,allow_roche=.true.)

    write(iunit,"(/,a)") '# timestepping'
    call write_inopt(norbits,'norbits','maximum number of binary orbits',iunit)
    call write_inopt(deltat,'deltat','output interval as fraction of binary period',iunit)
 endif
 close(iunit)

end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  read options from .setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils,    only:open_db_from_file,inopts,read_inopt,close_db
 use io,              only:error,fatal
 use setstar,         only:read_options_stars
 use setorbit,        only:read_options_orbit
 use setunits,        only:read_options_and_set_units
 use sethierarchical, only:read_hierarchical_setupfile
 use units,           only:in_code_units
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr,ierr1,input_type
 logical :: legacy_roche
 type(inopts), allocatable :: db(:)
 real :: m1,m2

 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_options_and_set_units(db,nerr,gr)
 call read_options_stars(star,ieos,relax,write_rho_to_file,db,nerr,nstar)
 call read_inopt(input_type,'itype',db,errcount=nerr,min=0,max=iorbit_roche,default=0)
 ! Migrate the previous switch to itype=5 on the normal setup-file rewrite.
 call read_inopt(legacy_roche,'roche_circular',db,errcount=nerr,default=.false.)
 if (legacy_roche) then
    input_type = iorbit_roche
    print "(a)",' Replacing roche_circular=T with itype=5 in the setup file; review and rerun.'
    nerr = nerr+1
 endif
 if (input_type==iorbit_roche) then
    if (nstar /= 2) call fatal('setup_binary','itype=5 requires exactly two bodies (nstars=2)')
    if (star(1)%iprofile <= 0) call fatal('setup_binary','itype=5 requires star 1 to contain gas particles')
 endif
 if (nstar > 2) then
    call read_hierarchical_setupfile(db,nerr,nstar)
 elseif (nstar == 2) then
    call read_inopt(corotate,'corotate',db,errcount=nerr)
    call read_inopt(gwinspiral,'gwinspiral',db,errcount=nerr,default=.false.)
    m1 = 1.
    m2 = 1.
    if (input_type /= iorbit_roche) then
       m1 = in_code_units(star(1)%m,ierr,unit_type='mass')
       m2 = in_code_units(star(2)%m,ierr1,unit_type='mass')
       if (ierr /= 0 .or. ierr1 /= 0) then
          print*,' ERROR: error converting masses to code units'
          nerr = nerr + 1
       endif
    endif
    call read_options_orbit(orbit,m1,m2,db,nerr,input_type=input_type,allow_roche=.true.)
    call read_inopt(norbits,'norbits',db,errcount=nerr)
    call read_inopt(deltat,'deltat',db,errcount=nerr,default=0.1)
 endif
 call close_db(db)

 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
