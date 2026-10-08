!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup routine for uniform distribution
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - Bzero       : *magnetic field strength in code units*
!   - cs0         : *initial sound speed in code units*
!   - dust_to_gas : *dust-to-gas ratio*
!   - ilattice    : *lattice type (1=cubic, 2=closepacked)*
!   - nx          : *number of particles in x direction*
!   - rhozero     : *initial density in code units*
!   - xmax        : *xmax boundary*
!   - xmin        : *xmin boundary*
!   - ymax        : *ymax boundary*
!   - ymin        : *ymin boundary*
!   - zmax        : *zmax boundary*
!   - zmin        : *zmin boundary*
!
! :Dependencies: boundary, cooling, cooling_ism, dim, eos, infile_utils,
!   io, mpidomain, options, part, physcon, radiation_utils, set_dust,
!   setunits, setup_params, timestep, unifdis, units
!
 use dim,          only:use_dust,mhd,gr
 use options,      only:use_dustfrac
 use setup_params, only:rhozero
 implicit none
 public :: setpart

 integer           :: npartx,ilattice
 real              :: cs0,xmini,xmaxi,ymini,ymaxi,zmini,zmaxi,Bzero

 !--change default defaults to reproduce the test from Section 5.6.7 of Price+(2018)
 logical :: BalsaraKim = .false.

 !--dust
 real    :: dust_to_gas

 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:maxvxyzu,h2chemistry,use_dustgrowth,do_radiation
 use setup_params, only:npart_total,ihavesetupB
 use io,           only:master
 use unifdis,      only:set_unifdis,latticetype,get_xyzmin_xyzmax_exact
 use boundary,     only:xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound,set_boundary
 use part,         only:Bxyz,periodic,abundance,igas,iHI,dustfrac,ndustsmall,&
                        ndusttypes,grainsize,graindens,dustprop,rad
 use physcon,      only:pi,mass_proton_cgs,kboltz,years,pc,solarm,micron
 use set_dust,     only:set_dustfrac
 use setunits,     only:dist_unit,mass_unit
 use units,        only:unit_density,udist
 use mpidomain,    only:i_belong
 use eos,          only:gmw
 use options,      only:icooling,alpha,alphau
 use timestep,     only:dtmax,tmax,C_cour,C_force,C_cool,tolv
 use cooling,      only:Tfloor
 use cooling_ism,  only:abundc,abundo,abundsi,abunde,dust_to_gas_ratio,iphoto
 use radiation_utils, only:set_radiation_and_gas_temperature_equal
 use infile_utils, only:get_options,infile_exists
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma
 real,              intent(inout) :: hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 real    :: totmass,deltax
 integer :: i,ierr

 !--NEW VARIABLES FOR GAUSSIAN STRETCH
  real :: u_x, u_y, u_z, r_sq
  real :: sigma_gauss
  real :: center_x, center_y, center_z
 !
 !--general parameters
 !
 time = 0.
 if (maxvxyzu < 4) then
    gamma = 1.
 else
    gamma = 5./3.
 endif
 !
 ! default units
 !
 mass_unit = 'solarm'
 dist_unit = 'pc'
 !
 ! set boundaries to default values
 !
 xmini = xmin; xmaxi = xmax
 ymini = ymin; ymaxi = ymax
 zmini = zmin; zmaxi = zmax
 !
 ! set default values for input parameters
 !
 npartx   = 64
 ilattice = 1
 rhozero  = 1.
 if (gr) then
    cs0 = 1.e-4
 else
    cs0 = 1.
 endif
 if (use_dust) then
    use_dustfrac = .true.
    dust_to_gas  = 0.01
 endif
 if (BalsaraKim) then
    ! there is a typo in Price+ (2018) in stating the physical density;
    ! this mass_unit yields the correct value of 2.3e-24g/cm^3
    mass_unit = '33982786.25*solarm'
    dist_unit = 'kpc'
    xmini = -0.1; xmaxi = 0.1
    ymini = -0.1; ymaxi = 0.1
    zmini = -0.1; zmaxi = 0.1
    Bzero = 0.056117
    cs0   = sqrt(0.3*gamma)
    ilattice = 2
    if (.not.infile_exists(fileprefix)) then
       tmax    = 0.035688
       dtmax   = 1.250d-4
       C_cour  = 0.200
       C_force = 0.150
       C_cool  = 0.15
       tolv    = 1.0d3
       alpha   = 1
       alphau  = 0.1
       gmw     = 1.22
       Tfloor  = 3.
       if (h2chemistry) then
          ! flags controlling h2chemistry
          icooling = 1
          abundc  = 2.0d-4
          abundo  = 4.5d-4
          abundsi = 3.0d-5
          abunde  = 2.0d-4
          iphoto  = 0
          dust_to_gas_ratio = 0.010
       else
          icooling = 0
       endif
    endif
 endif
 !
 ! get disc setup parameters from file or interactive setup
 !
 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'
 !
 ! set dust grain sizes
 !
 if (use_dust) then
    ndustsmall   = 1
    ndusttypes   = 1
    grainsize(1) = 1.*micron/udist
    graindens(1) = 3./unit_density
 endif
 !
 ! set boundaries
 !
 call set_boundary(xmini,xmaxi,ymini,ymaxi,zmini,zmaxi)
 !
 ! setup particles
 !
 deltax = dxbound/npartx
 npart = 0
 npart_total = 0

 call set_unifdis(latticetype(ilattice),id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                  deltax,hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong)
 call get_xyzmin_xyzmax_exact(latticetype(ilattice),xmin,xmax,ymin,ymax,zmin,zmax,ierr,deltax,npartx)

! --- GAUSSIAN COORDINATE STRETCHING ---
  ! Define the width of your Gaussian (e.g., 15% of the box size)
  sigma_gauss = 0.15 * (xmax - xmin) 
  
  ! Define the center of the box
  center_x = xmin + (xmax - xmin) / 2.0
  center_y = ymin + (ymax - ymin) / 2.0
  center_z = zmin + (zmax - zmin) / 2.0

  do i=1, npart
     ! 1. Normalize current uniform coordinates to a 0.0 to 1.0 range
     u_x = (xyzh(1,i) - xmin) / (xmax - xmin)
     u_y = (xyzh(2,i) - ymin) / (ymax - ymin)
     u_z = (xyzh(3,i) - zmin) / (zmax - zmin)

     ! 2. Stretch coordinates using the inverse normal CDF
     ! (Calls a custom function we will add at the bottom)
     xyzh(1,i) = center_x + sigma_gauss * inv_ndf(u_x)
     xyzh(2,i) = center_y + sigma_gauss * inv_ndf(u_y)
     xyzh(3,i) = center_z + sigma_gauss * inv_ndf(u_z)

     ! 3. CRITICAL: Adjust the smoothing length (h)
     ! Because particles are now packed tighter in the center, 
     ! we must shrink 'h' so they don't have too many neighbors.
     ! In 3D SPH, h scales with density as h_new = h_old * (rho_old / rho_new)^(1/3)
     r_sq = (xyzh(1,i)-center_x)**2 + (xyzh(2,i)-center_y)**2 + (xyzh(3,i)-center_z)**2
     
     ! The density increases exponentially toward the center, so h must decrease.
     ! exp(r_sq / (2 * sigma_gauss^2)) ^ (1/3) = exp(r_sq / (6 * sigma_gauss^2))
     xyzh(4,i) = xyzh(4,i) * exp(r_sq / (6.0 * sigma_gauss**2))
  enddo
  ! --- END GAUSSIAN STRETCHING ---

 npartoftype(:) = 0
 npartoftype(igas) = npart
 print*,' npart = ',npart,npart_total

 totmass = rhozero*dxbound*dybound*dzbound
 massoftype = totmass/npart_total
 if (id==master) print*,' particle mass = ',massoftype(igas)
 if (id==master) print*,' initial sound speed = ',cs0,' pressure = ',cs0**2/gamma

 if (maxvxyzu < 4 .or. gamma <= 1.) then
    polyk = cs0**2
 else
    polyk = 0.
 endif
 do i=1,npart
    vxyzu(1:3,i) = 0.
    if (maxvxyzu >= 4 .and. gamma > 1.) vxyzu(4,i) = cs0**2/(gamma*(gamma-1.))
 enddo

 if (use_dustfrac) then
    do i=1,npart
       call set_dustfrac(dust_to_gas,dustfrac(:,i))
    enddo
 endif

 if (use_dustgrowth) then
    do i=1,npart
       dustprop(1,i) = grainsize(1)
       dustprop(2,i) = graindens(1)
    enddo
 endif

 if (h2chemistry) then
    do i=1,npart
       abundance(:,i)   = 0.
       abundance(iHI,i) = 1.  ! assume all atomic hydrogen initially
    enddo
 endif

 if (mhd .and. balsarakim) then
    Bxyz = 0.
    do i = 1,npart
       Bxyz(1,i) = Bzero
    enddo
    ihavesetupB = .true.
 endif

 if (do_radiation) call set_radiation_and_gas_temperature_equal(npart,xyzh,vxyzu,massoftype,rad)

end subroutine setpart

!------------------------------------------------------------------------
!+
!  write setup file
!+
!------------------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 use setunits,     only:write_options_units
 character(len=*), intent(in) :: filename
 integer :: iunit

 print "(/,a)",' writing setup options file '//trim(filename)
 open(newunit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for uniform setup routine'

 call write_options_units(iunit)
 !
 ! boundaries
 !
 write(iunit,"(/,a)") '# boundaries'
 call write_inopt(xmini,'xmin','xmin boundary',iunit)
 call write_inopt(xmaxi,'xmax','xmax boundary',iunit)
 call write_inopt(ymini,'ymin','ymin boundary',iunit)
 call write_inopt(ymaxi,'ymax','ymax boundary',iunit)
 call write_inopt(zmini,'zmin','zmin boundary',iunit)
 call write_inopt(zmaxi,'zmax','zmax boundary',iunit)
 !
 ! other parameters
 !
 write(iunit,"(/,a)") '# setup'
 call write_inopt(npartx,'nx','number of particles in x direction',iunit)
 call write_inopt(rhozero,'rhozero','initial density in code units',iunit)
 call write_inopt(cs0,'cs0','initial sound speed in code units',iunit)
 if (use_dustfrac) then
    call write_inopt(dust_to_gas,'dust_to_gas','dust-to-gas ratio',iunit)
 endif
 if (mhd .and. balsarakim) then
    call write_inopt(Bzero,'Bzero','magnetic field strength in code units',iunit)
 endif
 call write_inopt(ilattice,'ilattice','lattice type (1=cubic, 2=closepacked)',iunit)
 close(iunit)

end subroutine write_setupfile

!------------------------------------------------------------------------
!+
!  read setup file
!+
!------------------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use setunits,     only:read_options_and_set_units
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",' reading setup options from '//trim(filename)
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 !
 ! units
 !
 call read_options_and_set_units(db,nerr,gr)
 !
 ! boundaries
 !
 call read_inopt(xmini,'xmin',db,errcount=nerr)
 call read_inopt(xmaxi,'xmax',db,min=xmini,errcount=nerr)
 call read_inopt(ymini,'ymin',db,errcount=nerr)
 call read_inopt(ymaxi,'ymax',db,min=ymini,errcount=nerr)
 call read_inopt(zmini,'zmin',db,errcount=nerr)
 call read_inopt(zmaxi,'zmax',db,min=zmini,errcount=nerr)
 !
 ! other parameters
 !
 call read_inopt(npartx,'nx',db,min=8,errcount=nerr)
 call read_inopt(rhozero,'rhozero',db,min=0.,errcount=nerr)
 call read_inopt(cs0,'cs0',db,min=0.,errcount=nerr)
 if (use_dustfrac) then
    call read_inopt(dust_to_gas,'dust_to_gas',db,min=0.,errcount=nerr)
 endif
 if (mhd .and. balsarakim) then
    call read_inopt(Bzero,'Bzero',db,min=0.,errcount=nerr)
 endif
 call read_inopt(ilattice,'ilattice',db,min=1,max=2,errcount=nerr)
 call close_db(db)

 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

!------------------------------------------------------------------------
!+
!  Rational approximation for the Inverse Normal CDF
!  Maps a uniform probability (0 to 1) to a standard Gaussian (-inf to inf)
!+
!------------------------------------------------------------------------
  function inv_ndf(p) result(x)
      implicit none
      real, intent(in) :: p
      real :: x, p_clamped
      real, parameter :: a1 = -3.9696830, a2 =  2.2094609, a3 = -2.7592851
      real, parameter :: a4 =  1.3835775, a5 = -0.0680889, a6 =  0.0000000
      real, parameter :: b1 = -0.5447609, b2 =  1.6158583, b3 = -0.1556989
      real, parameter :: b4 =  0.0668013, b5 = -0.0013280
      real :: q, r

      ! Prevent infinities at the absolute edges of the box
      p_clamped = max(1.e-6, min(p, 0.999999))
      q = p_clamped - 0.5
      r = q * q
      
      ! Simple rational approximation (Acklam's algorithm)
      x = (((((a6*r+a5)*r+a4)*r+a3)*r+a2)*r+a1)*q / &
          (((((b5*r+b4)*r+b3)*r+b2)*r+b1)*r+1.0)
  end function inv_ndf

end module setup
