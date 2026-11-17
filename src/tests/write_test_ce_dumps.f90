!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
program write_test_ce_dumps
!
! Write small EOS 20 radiation fixtures for test_analysis_ce_rad.py.
! Gas/radiation temperatures differ and the last full dump has a dead particle.
!
! :References: None
!
! :Owner:
!
! :Usage: write_test_ce_dumps [no arguments]
!
! :Dependencies: dim, eos, eos_gasradrec, initial, io, ionization_mod,
!   kernel, memory, mpiutils, part, physcon, radiation_utils,
!   readwrite_dumps, units
!
 use dim, only:do_radiation
 use io, only:id,nprocs,set_io_unit_numbers
 use mpiutils, only:init_mpi,finalise_mpi
 use initial, only:initialise,finalise
 use memory, only:allocate_memory
 use part, only:init_part,npart,npartoftype,igas,iphase,isetphase,xyzh,vxyzu,massoftype,hfact,&
                rho,nptmass,xyzmh_ptmass,vxyz_ptmass,poten,eos_vars,iX,iZ,imu,itemp,igasP,rad,iradxi
 use kernel, only:hfact_default
 use physcon, only:solarm,solarr
 use units, only:set_units,unit_density,unit_ergg,unit_pressure,unit_energ
 use eos, only:ieos,init_eos,use_var_comp,gamma,polyk
 use eos_gasradrec, only:irecomb,calc_uP_from_rhoT_gasradrec
 use ionization_mod, only:get_erec_cveff
 use radiation_utils, only:radxi_from_Trad
 use readwrite_dumps, only:write_fulldump,write_smalldump
 implicit none
 integer :: i,ierr,iu
 real :: temp(2),trad(2),u,p,imui,erec,cveff,chemical,grav,internal,kinetic,potential,total

 if (.not.do_radiation) stop 'fixtures require RADIATION=yes'
 call init_mpi(id,nprocs)
 call set_io_unit_numbers
 call initialise()
 call allocate_memory(4_8)
 call init_part()
 call set_units(mass=solarm,dist=solarr,G=1.d0)
 ieos = 20
 gamma = 5./3.
 polyk = 0.
 irecomb = 0
 use_var_comp = .true.
 call init_eos(ieos,ierr)
 if (ierr /= 0) stop 1
 npart = 2
 npartoftype = 0
 npartoftype(igas) = npart
 massoftype(igas) = 1.e-3
 hfact = hfact_default
 iphase(1:3) = isetphase(igas,iactive=.true.)
 xyzh = 0.
 xyzh(1,1:2) = [1000.,-1500.]
 rho(1:2) = 1.e-9/unit_density
 xyzh(4,1:2) = hfact*(massoftype(igas)/rho(1:2))**(1./3.)
 vxyzu = 0.
 vxyzu(2,1:2) = [0.02,-0.01]
 poten(1:2) = real([-1.e-6,-2.e-6],kind=kind(poten))
 nptmass = 2
 xyzmh_ptmass = 0.
 xyzmh_ptmass(4,1:2) = [3.,1.]
 xyzmh_ptmass(1,2) = 4000.
 vxyz_ptmass = 0.
 vxyz_ptmass(2,1:2) = [0.01,-0.03]
 eos_vars(iX,1:2) = [0.37,0.74]
 eos_vars(iZ,1:2) = [0.02,0.01]
 temp = [4000.,10000.]
 trad = [8000.,3000.]
 chemical = 0.
 do i=1,npart
    call calc_uP_from_rhoT_gasradrec(rho(i)*unit_density,temp(i),eos_vars(iX,i),&
                                    1.-eos_vars(iX,i)-eos_vars(iZ,i),u,p,imui,do_radiation=.true.)
    vxyzu(4,i) = u/unit_ergg
    eos_vars(imu,i) = 1./imui
    eos_vars(itemp,i) = temp(i)
    eos_vars(igasP,i) = p/unit_pressure
    rad(iradxi,i) = radxi_from_Trad(rho(i),trad(i))
    call get_erec_cveff(log10(rho(i)*unit_density),temp(i),eos_vars(iX,i),&
                       1.-eos_vars(iX,i)-eos_vars(iZ,i),erec,cveff)
    chemical = chemical + massoftype(igas)*erec/unit_ergg
 enddo
 grav = sum(real(poten(1:2))) - massoftype(igas)*(3./1000.+3./1500.)
 internal = massoftype(igas)*sum(vxyzu(4,1:2)+rad(iradxi,1:2))
 kinetic = 0.5*massoftype(igas)*sum(vxyzu(1:3,1:2)**2) &
           +0.5*sum(xyzmh_ptmass(4,1:2)*vxyz_ptmass(2,1:2)**2)
 potential = grav - massoftype(igas)*(1./3000.+1./5500.) - 3./4000.
 total = kinetic + potential + internal
 open(newunit=iu,file='ce_reference.txt',status='replace')
 write(iu,*) grav,grav+internal-chemical,grav+internal,total,massoftype(igas)*sum(rad(iradxi,1:2))
 close(iu)
 open(newunit=iu,file='ce_unit.txt',status='replace')
 write(iu,*) unit_energ
 close(iu)
 call write_fulldump(0.,'test_00000')
 call write_smalldump(0.5,'test_00001')
 npart = 3
 npartoftype(igas) = npart
 xyzh(4,3) = -1.
 vxyzu(4,3) = 0.
 poten(3) = -100.
 call write_fulldump(1.,'test_00002')
 open(newunit=iu,file='test.in',status='replace')
 write(iu,'(a)') 'ieos = 20','irecomb = 0','idust_opacity = 0'
 close(iu)
 call finalise()
 call finalise_mpi()

end program write_test_ce_dumps
