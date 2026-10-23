!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module relaxstar
!
! Automated relaxation of stellar density profile,
!   iterating towards hydrostatic equilibrium
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - damp_radrelax    : *damp velocities during radiative relaxation*
!   - maxits           : *maximum number of relaxation iterations*
!   - maxits_radrelax  : *maximum number of steps in radiative relaxation*
!   - relax_radiation  : *relax star with radiation transport after setting composition*
!   - tol_ekin         : *tolerance on ekin/epot to stop relaxation*
!   - tol_radrelax     : *tolerance on average absolute change in u+xi per dynamical time*
!
! :Dependencies: apr, checksetup, damping, deriv, dim, dump_utils,
!   energies, eos, externalforces, fileutils, infile_utils, initial, io,
!   io_summary, memory, neighkdtree, options, part, physcon, ptmass,
!   readwrite_dumps, setstar_utils, step_lf_global, table_utils, units
!
 implicit none
 public :: relax_star,relax_star_radiation,write_options_relax,read_options_relax

 real,    public  :: tol_ekin = 1.e-7 ! criteria for being converged
 integer, public  :: maxits = 1000

 logical, public  :: relax_radiation = .false. ! relax with radiation transport after the hydrodynamic relaxation
 logical, public  :: damp_radrelax = .false.   ! damp velocities during the radiative relaxation
 real,    public  :: tol_radrelax = 1.e-3      ! average absolute change in u+xi per dynamical time
 integer, public  :: maxits_radrelax = 10000

 real,    private :: gammaprev,hfactprev,mass1prev
 integer, private :: ieos_prev

 integer, public :: ierr_setup_errors = 1, &
                    ierr_no_pressure  = 2, &
                    ierr_unbound = 3, &
                    ierr_notconverged = 4, &
                    ierr_radrelax_unsupported = 5, &
                    ierr_radrelax_notconverged = 6

 private

contains

!----------------------------------------------------------------
!+
!  relax a star to hydrostatic equilibrium. We run the main
!  code but with a fake equation of state, low neighbour number
!  and fixing the entropy as a function of r
!
!  IN:
!    rhotab(nt) - tabulated density as function of r (in code units)
!    pr(nt)     - tabulated pressure as function of r (in code units)
!    r(nt)      - radius for each point in the table
!
!  IN/OUT:
!    xyzh(:,:) - positions and smoothing lengths of all particles
!+
!----------------------------------------------------------------
subroutine relax_star(nt,rhotab,pr,temp,r,npart,xyzh,use_var_comp,Xfrac,Yfrac,mu,&
                      iptmass_core,xyzmh_ptmass,ierr,npin,label,write_dumps,density_error,energy_error,mtab)
 use table_utils,     only:yinterp
 use deriv,           only:get_derivs_global
 use dim,             only:maxp,maxvxyzu,gr,gravity,use_apr
 use part,            only:vxyzu,rad,eos_vars,massoftype,igas,apr_level,fxyzu,init_rho_from_h
 use step_lf_global,  only:init_step,step
 use initial,         only:initialise
 use memory,          only:allocate_memory
 use energies,        only:compute_energies,ekin,epot,etherm
 use checksetup,      only:check_setup
 use io,              only:error,warning,fatal,id,master
 use fileutils,       only:getnextfilename
 use readwrite_dumps, only:write_fulldump
 use eos,             only:gamma,eos_outputs_mu,ieos,eos_has_pressure_without_u,polyk
 use physcon,         only:pi
 use options,         only:iexternalforce
 use io_summary,      only:summary_initialise
 use setstar_utils,   only:set_star_thermalenergy,set_star_composition
 use apr,             only:init_apr,update_apr
 use neighkdtree,     only:allocate_neigh
 integer,           intent(in)    :: nt,iptmass_core
 integer,           intent(inout) :: npart
 real,              intent(in)    :: rhotab(nt),pr(nt),temp(nt),r(nt)
 logical,           intent(in)    :: use_var_comp
 real, allocatable, intent(in)    :: Xfrac(:),Yfrac(:),mu(:)
 real,              intent(inout) :: xyzh(:,:),xyzmh_ptmass(:,:)
 integer,           intent(out)   :: ierr
 integer,           intent(in),  optional :: npin
 character(len=*),  intent(in),  optional :: label
 logical,           intent(in),  optional :: write_dumps
 real,              intent(out), optional :: density_error,energy_error
 real,              intent(in),  optional :: mtab(nt)
 integer :: nits,nerr,nwarn,iunit,i1
 real    :: t,dt,dtmax,rmserr,rstar,mstar,tdyn,x0(3),mtot
 real    :: entrop(nt),utherm(nt),mr(nt),rmax,dtext,dtnew
 logical :: converged,use_step,restart
 logical, parameter :: fix_entrop = .true. ! fix entropy instead of thermal energy
 logical :: write_files
 character(len=20) :: filename,mylabel

 i1 = 0
 if (present(npin)) i1 = npin  ! starting position in particle array
 !
 ! optional argument to not write files to disk
 !
 write_files = .true.
 if (present(write_dumps)) write_files = write_dumps
 !
 ! label for relax_star snapshots
 !
 mylabel = ''
 if (present(label)) mylabel = label
 !
 ! save settings and set a bunch of options
 !
 ierr = 0
 x0 = 0.
 if (gr) x0 = [1000.,0.,0.]   ! for GR need to shift star away from origin to avoid singularities
 rstar = maxval(r)

 !The mass profile is constructed from the density profile unless mtab is present.
 !If mtab is present, it is used as the mass profile.
 if (present(mtab)) then
    mr=mtab
 else
    mr = get_mr(rhotab,r)
 endif

 mstar = mr(nt)   ! mstar is the mass of the star excluding the core
 mtot = mstar
 if (iptmass_core > 0) mtot = mtot + xyzmh_ptmass(4,iptmass_core) ! mtot is the total mass of the star
 tdyn  = 2.*pi*sqrt(rstar**3/(32.*mtot))
 if (id==master) print*,'rstar = ',rstar,' mstar = ',mtot, ' tdyn = ',tdyn
 !
 ! see if we can restart or skip the relaxation process
 ! based on a previous run
 !
 filename = 'relax'//trim(mylabel)//'_00000'
 call check_for_existing_file(filename,npart,massoftype(igas),&
                              xyzh,vxyzu,restart,ierr)
 !
 ! quit with fatal error if non-matching file found, otherwise
 ! this will be overwritten
 !
 if (write_files .and. ierr /= 0) then
    if (id==master) print "(a)",' ERROR: pre-existing relaxed star dump(s) do not match'
    call fatal('relax_star','please delete relax'//trim(mylabel)//'_* and restart...')
 endif

 call set_options_for_relaxation(tdyn)
 call summary_initialise()
 if (gr) call shift_star_origin(x0,i1,npart,xyzh,iptmass_core,xyzmh_ptmass)
 !
 ! check particle setup is sensible
 !
 call check_setup(nerr,nwarn,restart=.true.) ! restart=T allows accreted/masked particles
 if (nerr > 0) then
    call error('relax_star','cannot relax star because particle setup contains errors')
    call restore_original_options(i1,npart)
    ierr = ierr_setup_errors
    return
 endif
 use_step = .false.

 if (iexternalforce > 0 .and. (.not. gr)) then
    call warning('relax_star','asynchronous shifting not implemented with external forces: evolving in time instead')
    use_step = .true.
 endif

 ! if using apr, options set in setup file but needs to be initialised here
 if (use_apr) then
    call allocate_neigh
    call init_apr(apr_level,ierr)
    call update_apr(npart,xyzh,vxyzu,fxyzu,apr_level)
 endif

 !
 ! define utherm(r) based on P(r) and rho(r)
 ! and use this to set the thermal energy of all particles
 !
 where (rhotab > epsilon(0.) .and. gamma > 1.)
    entrop = pr/rhotab**gamma
    utherm = pr/(rhotab*(gamma-1.))
 elsewhere
    entrop = 0.
    utherm = 0.
 end where

 if (any(utherm(1:nt-1) <= 0.) .and. .not.eos_has_pressure_without_u(ieos)) then
    call error('relax_star','relax-o-matic needs non-zero pressure array set in order to work')
    call restore_original_options(i1,npart)
    ierr = ierr_no_pressure
    return
 endif
 call init_rho_from_h(i1+1,npart)
 call reset_u_and_get_errors(i1,npart,xyzh,vxyzu,x0,rad,nt,mr,rhotab,&
                             utherm,entrop,fix_entrop,rmax,rmserr)
 !
 ! compute derivatives the first time around (needed if using actual step routine)
 !
 t = 0.
 call allocate_memory(int(min(2*npart,maxp),kind=8))
 call get_derivs_global()
 call reset_u_and_get_errors(i1,npart,xyzh,vxyzu,x0,rad,nt,mr,rhotab,&
                             utherm,entrop,fix_entrop,rmax,rmserr)
 call compute_energies(t)
 !
 ! perform sanity checks
 !
 if (etherm > abs(epot) .and. ieos /= 9) then
    call error('relax_star','cannot relax star because it is unbound (etherm > epot)')
    if (id==master) print*,' Etherm = ',etherm,' Epot = ',Epot
    if (maxvxyzu < 4) print "(/,a,/)",' *** Try compiling with ISOTHERMAL=no instead... ***'
    call restore_original_options(i1,npart)
    ierr = ierr_unbound
    return
 endif
 if (id==master) print "(/,3(a,1pg11.3),/,a,1pg11.3,a,i0)",&
   ' RELAX-A-STAR-O-MATIC: Etherm:',etherm,' Epot:',Epot, ' R*:',maxval(r), &
   '       WILL stop when Ekin/Epot < ',tol_ekin,' OR Iter=',maxits

 if (write_files) then
    if (gr) call shift_star_origin(-x0,i1,npart,xyzh,iptmass_core,xyzmh_ptmass)
    if (.not.restart) call write_fulldump(t,filename)
    ! move star back to 100,0,0
    if (gr) call shift_star_origin(x0,i1,npart,xyzh,iptmass_core,xyzmh_ptmass)
    open(newunit=iunit,file='relax'//trim(mylabel)//'.ev',status='replace')
    write(iunit,"(a)") '# nits,rmax,etherm,epot,ekin/epot,L2_{err}'
 endif
 converged = .false.
 dt = epsilon(0.) ! To avoid error in sink-gas substepping
 dtext = huge(dtext)
 if (use_step) then
    dtmax = tdyn
    call init_step(npart,t,dtmax)
 endif
 nits = 0
 rmserr = 0.
 ekin = 0.
 epot = 1. ! to avoid compiler warning
 do while (.not. converged .and. nits < maxits)
    nits = nits + 1
    !
    ! shift particles by one "timestep"
    !
    t = t + dt
    if (use_step) then
       call step(npart,npart,t,dt,dtext,dtnew)
       dt = dtnew
    else
       call shift_particles(i1,npart,xyzh,vxyzu,dt)
    endif
    ! if using apr, update here
    if (use_apr) call update_apr(npart,xyzh,vxyzu,fxyzu,apr_level)
    !
    ! reset thermal energy and calculate information
    !
    call reset_u_and_get_errors(i1,npart,xyzh,vxyzu,x0,rad,nt,mr,&
         rhotab,utherm,entrop,fix_entrop,rmax,rmserr)
    !
    ! compute energies and check for convergence
    !
    call compute_energies(t)
    converged = (ekin > 0. .and. ekin/abs(epot) < tol_ekin) !.and. rmserr < 0.01*tol_dens)
    !
    ! print information to screen
    !
    if (use_step) then
       if (id==master .and. mod(nits,10)==0 .or. nits==1) then
          print "(a,es10.3,a,2pf6.2,2(a,1pg11.3))",&
                ' Relaxing star: t/dyn:',t/tdyn,', dens error:',rmserr,'%, R*:',rmax, &
                ' Ekin/Epot:',ekin/abs(epot)
       endif
    else
       if (id==master .and. mod(nits,10)==0 .or. nits==1) then
          print "(a,i0,a,i0,a,2pf6.2,2(a,1pg11.3))",&
                ' Relaxing star: Iter ',nits,'/',maxits, &
                ', dens error:',rmserr,'%, R*:',rmax,' Ekin/Epot:',ekin/abs(epot)
       endif
    endif
    !
    ! additional diagnostic output, mainly for debugging/checking
    !
    if (write_files) then
       !
       ! write information to the relax.ev file
       !
       write(iunit,*) nits,rmax,etherm,epot,ekin/abs(epot),rmserr
       !
       ! write dump files
       !
       if (mod(nits,100)==0 .or. ((nits==maxits .or. converged).and.nits > 1)) then
          filename = getnextfilename(filename)

          !
          ! before writing a file, set the real thermal energy profile
          ! so the file is useable as a starting file for the main calculation
          !
          if (use_var_comp) call set_star_composition(eos_outputs_mu(ieos_prev),&
                                 npart,xyzh,Xfrac,Yfrac,mu,mr,eos_vars,npin=i1,x0=x0)

          ! Use the same composition treatment as the final conversion in set_star.
          if (maxvxyzu==4) call set_star_thermalenergy(ieos_prev,rhotab,pr,temp,&
                                r,nt,npart,xyzh,vxyzu,rad,eos_vars,.true.,&
                                use_var_comp=use_var_comp,initialtemp=1.e3,polyk_in=polyk,npin=i1,x0=x0)

          ! write relaxation snapshots
          if (write_files) then
             ! move star back to 0,0,0
             if (gr) call shift_star_origin(-x0,i1,npart,xyzh,iptmass_core,xyzmh_ptmass)
             ! write snapshot
             call write_fulldump(t,filename)
             ! move star back to 100,0,0
             if (gr) call shift_star_origin(x0,i1,npart,xyzh,iptmass_core,xyzmh_ptmass)
          endif

          ! flush the relax.ev file
          call flush(iunit)

          ! restore the fake thermal energy profile
          call reset_u_and_get_errors(i1,npart,xyzh,vxyzu,x0,rad,nt,mr,rhotab,&
               utherm,entrop,fix_entrop,rmax,rmserr)
       endif
    endif
 enddo
 if (write_files) close(iunit)
 !
 ! optional diagnostics
 !
 if (present(density_error)) density_error = rmserr
 if (present(energy_error))  energy_error = ekin/abs(epot)
 !
 ! warn if relaxation finished due to hitting nits=nitsmax
 !
 if (.not.converged) then
    call warning('relax_star','relaxation did not converge, just reached max iterations')
    ierr = ierr_notconverged
 else
    if (id==master) print "(5(a,/),/,a,/,/,a,/)",&
    "                             _                    _ ",&
    "                    _ __ ___| | __ ___  _____  __| |",&
    "                   | '__/ _ \ |/ _` \ \/ / _ \/ _` |",&
    "                   | | |  __/ | (_| |>  <  __/ (_| |",&
    "          o  o  o  |_|  \___|_|\__,_/_/\_\___|\__,_|  o  o  o",&
    "          ************  Check the profile using: ************", &
    "                  splash "//trim(filename)//" -y density -x r",&
    "          ***************************************************"
 endif
 !
 ! unfake some things
 !
 call restore_original_options(i1,npart)

 ! move star back to 0,0,0
 if (gr) call shift_star_origin(-x0,i1,npart,xyzh,iptmass_core,xyzmh_ptmass)

end subroutine relax_star

!----------------------------------------------------------------
!+
!  relax the star further with radiation transport, using the real
!  equation of state and the implicit radiation solver. Particles
!  evolve hydrodynamically (damped) together with u and xi until
!  the mass-weighted mean absolute change in u+xi per dynamical time is below
!  tol_radrelax. Must be called after the composition and thermal
!  energy / radiation energy have been set.
!
!  IN:
!    rstar,mtot - radius and total mass of the star (code units)
!    npin       - number of particles to ignore at start of array
!    label      - label for the rad_relax<label>_* snapshots
!    write_dumps- if false, do not write snapshots or the .ev file
!                 if previous rad_relax<label>_* snapshots exist, the
!                 relaxation is resumed from the latest matching dump
!+
!----------------------------------------------------------------
subroutine relax_star_radiation(npart,xyzh,rstar,mtot,npin,ierr,label,write_dumps)
 use deriv,          only:get_derivs_global
 use dim,            only:maxvxyzu,do_radiation,gr
 use eos_gasradrec,  only:report_newton_stalls
 use eos,            only:ieos
 use part,           only:vxyzu,rad,iradxi,init_rho_from_h,massoftype,igas,apr_level,aprmassoftype
 use step_lf_global, only:init_step,step
 use damping,        only:damp
 use options,        only:idamp,implicit_radiation
 use io,             only:warning,id,master
 use physcon,        only:pi
 use fileutils,      only:getnextfilename
 use readwrite_dumps, only:write_fulldump
 use apr,             only:use_apr
 use radiation_utils, only:write_options_radiation
 integer,  intent(in)    :: npin
 integer,  intent(inout) :: npart
 real,     intent(inout) :: xyzh(:,:)
 real,     intent(in)    :: rstar,mtot
 integer,  intent(out)   :: ierr
 character(len=*), intent(in), optional :: label
 logical,          intent(in), optional :: write_dumps
 integer :: i,nits,ncheck,iunit,nits_dumped
 real    :: t,dt,dtext,dtnew,tdyn,tlast,dtcheck,err,rmax,du_sum,dxi_sum,duxi_sum,mtot_sum,pmassi
 real, allocatable :: u0(:),xi0(:)
 logical :: converged,implicit_prev,write_files,restart
 character(len=20) :: filename,mylabel
 integer, parameter :: nstep_check = 10
 integer, parameter :: nstep_dump = 100

 ierr = 0
 if (.not.do_radiation .or. maxvxyzu < 4 .or. gr) then
    call warning('relax_star_radiation','radiative relaxation needs radiation, thermal energy and no GR')
    ierr = ierr_radrelax_unsupported
    return
 endif

 write_files = .true.
 if (present(write_dumps)) write_files = write_dumps
 mylabel = ''
 if (present(label)) mylabel = label

 tdyn = 2.*pi*sqrt(rstar**3/(32.*mtot))
 if (id==master) print "(/,a,es10.3,a,i0)",' RADIATIVE RELAXATION: average |max(du, dxi)|/tdyn tolerance = ',tol_radrelax,&
                                           ', max steps = ',maxits_radrelax

 filename = 'rad_relax'//trim(mylabel)//'_00000'
 call check_for_existing_file(filename,npart,massoftype(igas),xyzh,vxyzu,restart,ierr)
 if (ierr /= 0) return

 implicit_prev = implicit_radiation
 implicit_radiation = .true.
 if (id==master) call write_options_radiation(6)
 if (damp_radrelax) then
    idamp = 1
    damp = 1./tdyn  ! damping timescale of one dynamical time
 else
    idamp = 0
    damp = 0.
 endif
 vxyzu(1:3,npin+1:npart) = 0.

 call init_rho_from_h(npin+1,npart)
 t  = 0.
 dt = epsilon(0.)
 dtext = huge(dtext)
 call get_derivs_global()
 call init_step(npart,t,tdyn)

 allocate(u0(npart),xi0(npart))
 u0  = vxyzu(4,1:npart)
 xi0 = rad(iradxi,1:npart)
 tlast = t
 converged = .false.
 ncheck = 0
 nits = 0
 nits_dumped = 0
 if (write_files) then
    if (.not.restart) call write_fulldump(t,filename)
    if (restart) then
       open(newunit=iunit,file='rad_relax'//trim(mylabel)//'.ev',status='unknown',position='append')
    else
       open(newunit=iunit,file='rad_relax'//trim(mylabel)//'.ev',status='replace')
       write(iunit,"(a)") '# nits,t/tdyn,mass-averaged |d(max(u, xi))| per tdyn'
    endif
 endif
 do while (.not.converged .and. nits < maxits_radrelax)
    nits = nits + 1
    t = t + dt
    call step(npart,npart,t,dt,dtext,dtnew)
    dt = dtnew

    if (mod(nits,nstep_check)==0) then
       du_sum = 0.
       dxi_sum = 0.
       duxi_sum = 0.
       mtot_sum = 0.
       rmax = 0.
       do i = npin+1,npart
          if (use_apr) then
             pmassi = aprmassoftype(igas,apr_level(i))
          else
             pmassi = massoftype(igas)
          endif
          du_sum = du_sum + pmassi*abs(vxyzu(4,i)-u0(i))
          dxi_sum = dxi_sum + pmassi*abs(rad(iradxi,i)-xi0(i))
          duxi_sum = duxi_sum + pmassi*((vxyzu(4,i)-u0(i)) + (rad(iradxi,i)-xi0(i)))
          mtot_sum = mtot_sum + pmassi
          rmax = max(rmax,sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i))))
       enddo
       dtcheck = t - tlast
       err = max(du_sum,dxi_sum)/mtot_sum*tdyn/dtcheck
       u0  = vxyzu(4,1:npart)
       xi0 = rad(iradxi,1:npart)
       tlast = t
       converged = (err < tol_radrelax)
       ncheck = ncheck + 1
       if (write_files) write(iunit,*) nits,t/tdyn,err
       if (id==master) then
          write(*,'(a,i0,a,i0,a,1pg11.3,a,1pg11.3,a,1pg11.3,a,1pg11.3,a,1pg11.3)') &
             ' Relaxing star with radiation: Iter ',nits,'/',maxits_radrelax, &
             ', R*:',rmax,', t/tdyn:',t/tdyn, &
             ', average absolute change in u or xi per tdyn:',err, &
             ', radiation (du+dxi)/tdyn:',abs(duxi_sum)/mtot_sum*tdyn/dtcheck
       endif
       call report_newton_stalls(6,int(npart,kind=8),ieos)
    endif

    if (write_files .and. mod(nits,nstep_dump)==0) then
       filename = getnextfilename(filename)
       call write_fulldump(t,filename)
       call flush(iunit)
       nits_dumped = nits
    endif
 enddo

 if (write_files) then
    if (nits /= nits_dumped) then
       filename = getnextfilename(filename)
       call write_fulldump(t,filename)
    endif
    close(iunit)
 endif

 if (.not.converged) then
    call warning('relax_star_radiation','radiative relaxation did not converge, just reached max steps')
    ierr = ierr_radrelax_notconverged
 endif

 implicit_radiation = implicit_prev
 idamp = 0
 damp = 0.
 vxyzu(1:3,npin+1:npart) = 0.
 deallocate(u0,xi0)

end subroutine relax_star_radiation

!----------------------------------------------------------------
!+
!  shift particles: this is like timestepping but done
!  asynchronously. Each particle shifts by dx = 0.5*dt^2*a
!  where dt is the local courant timestep, i.e. h/c_s
!+
!----------------------------------------------------------------
subroutine shift_particles(i1,npart,xyzh,vxyzu,dtmin)
 use deriv, only:get_derivs_global
 use part,  only:fxyzu,fext,rho,xyzmh_ptmass,nptmass
 use ptmass,only:get_accel_sink_gas
 use eos,   only:get_spsound
 use options, only:ieos
 integer, intent(in)    :: i1,npart
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:)
 real,    intent(out)   :: dtmin
 real :: dx(3),dti,phi,rhoi,cs,hi
 integer :: i,nlargeshift
!
! shift particles asynchronously
!
 dtmin = huge(dtmin)
 nlargeshift = 0
 !$omp parallel do schedule(guided) default(none) &
!$omp shared(i1,npart,xyzh,vxyzu,fxyzu,fext,xyzmh_ptmass,nptmass,ieos,rho) &
 !$omp private(i,dx,dti,phi,cs,rhoi,hi) &
 !$omp reduction(min:dtmin) &
 !$omp reduction(+:nlargeshift)
 do i=i1+1,npart
    fext(1:3,i) = 0.
    if (nptmass > 0) then
       call get_accel_sink_gas(nptmass,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),&
                               xyzmh_ptmass,fext(1,i),fext(2,i),fext(3,i),phi)
    endif
    hi = xyzh(4,i)
    rhoi = rho(i)
    cs = get_spsound(ieos,xyzh(:,i),rhoi,vxyzu(:,i))
    dti = 0.3*hi/cs   ! h/cs
    dx  = 0.5*dti**2*(fxyzu(1:3,i) + fext(1:3,i))
    if (dot_product(dx,dx) > hi**2) then
       dx = dx / sqrt(dot_product(dx,dx)) * hi  ! Avoid large shift in particle position
       nlargeshift = nlargeshift + 1
    endif
    xyzh(1:3,i) = xyzh(1:3,i) + dx(:)
    vxyzu(1:3,i) = dx(:)/dti ! fake velocities, so we can measure kinetic energy
    dtmin = min(dtmin,dti)   ! used to print a "time" in the output (but it is fake)
 enddo
 !$omp end parallel do
 if (nlargeshift > 0) print*,'Warning: Restricted dx for ', nlargeshift, 'particles'
 !
 ! get forces on particles
 !
 call get_derivs_global()

end subroutine shift_particles

!----------------------------------------------------------------
!+
!  shift the origin of the star, mainly to avoid coordinate
!  singularities in GR code
!+
!----------------------------------------------------------------
subroutine shift_star_origin(x0,i1,npart,xyzh,iptmass_core,xyzmh_ptmass)
 integer, intent(in)    :: i1,npart,iptmass_core
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(in)    :: x0(3)
 real,    intent(inout) :: xyzmh_ptmass(:,:)
 integer :: i

 !$omp parallel do schedule(guided) default(none) &
 !$omp shared(xyzh,x0,i1,npart) &
 !$omp private(i)
 do i=i1+1,npart
    xyzh(1:3,i) = xyzh(1:3,i) + x0
 enddo
 !$omp end parallel do

 if (iptmass_core > 0 .and. iptmass_core <= size(xyzmh_ptmass,2)) then
    xyzmh_ptmass(1:3,iptmass_core) = xyzmh_ptmass(1:3,iptmass_core) + x0
 endif

end subroutine shift_star_origin

!----------------------------------------------------------------
!+
!  reset the thermal energy to be exactly p(r)/((gam-1)*rho(r))
!  according to the desired p(r) and rho(r)
!  also compute error between true rho(r) and desired rho(r)
!+
!----------------------------------------------------------------
subroutine reset_u_and_get_errors(i1,npart,xyzh,vxyzu,x0,rad,nt,mr,rhotab,&
                                  utherm,entrop,fix_entrop,rmax,rmserr)
 use table_utils,   only:yinterp
 use part,          only:rho,maxvxyzu
 use dim,           only:do_radiation
 use eos,           only:gamma
 use setstar_utils, only:get_mass_coord
 integer, intent(in)    :: i1,npart,nt
 real,    intent(in)    :: xyzh(:,:),x0(3),mr(nt),rhotab(nt),utherm(nt),entrop(nt)
 real,    intent(inout) :: vxyzu(:,:),rad(:,:)
 real,    intent(out)   :: rmax,rmserr
 logical, intent(in)    :: fix_entrop
 real :: ri,rhor,rhoi,rho1,mstar,massri
 real, allocatable :: mass_enclosed_r(:)
 integer :: i

 rho1 = yinterp(rhotab,mr,0.)
 rmax = 0.
 rmserr = 0.

 call get_mass_coord(i1,npart,xyzh,mass_enclosed_r,x0)
 mstar = mr(nt)

 !$omp parallel do schedule(guided) default(none) &
 !$omp shared(i1,npart,xyzh,vxyzu,x0,mass_enclosed_r,mr,rhotab,utherm,entrop) &
 !$omp shared(fix_entrop,gamma,rho) &
 !$omp private(i,ri,rhor,rhoi,massri) &
 !$omp reduction(+:rmserr) &
 !$omp reduction(max:rmax)
 do i = i1+1,npart
    ri = sqrt(dot_product(xyzh(1:3,i)-x0,xyzh(1:3,i)-x0))
    massri = mass_enclosed_r(i-i1)
    rhor = yinterp(rhotab,mr,massri) ! analytic rho(r)

    rhoi = rho(i) ! actual rho
    if (maxvxyzu >= 4) then
       if (fix_entrop .and. gamma > 1.) then
          vxyzu(4,i) = (yinterp(entrop,mr,massri)*rhoi**(gamma-1.))/(gamma-1.)
       else
          vxyzu(4,i) = yinterp(utherm,mr,massri)
       endif
    endif
    rmserr = rmserr + (rhor - rhoi)**2
    rmax   = max(rmax,ri)
 enddo
 !$omp end parallel do
 if (do_radiation) rad = 0.
 rmserr = sqrt(rmserr/npart)/rho1

 deallocate(mass_enclosed_r)

end subroutine reset_u_and_get_errors

!----------------------------------------------------------------
!+
!  set code options specific to relaxation calculations
!+
!----------------------------------------------------------------
subroutine set_options_for_relaxation(tdyn)
 use eos,  only:ieos,gamma,eos_has_pressure_without_u
 use part, only:hfact,maxvxyzu,gr
 use damping, only:damp,tdyn_s
 use options, only:idamp
 use units,          only:utime
 use externalforces, only:mass1
 real, intent(in) :: tdyn

 gammaprev = gamma
 hfactprev = hfact
 ieos_prev = ieos
 mass1prev = mass1
 !
 ! turn on settings appropriate to relaxation
 !
 if (maxvxyzu >= 4 .and. .not.eos_has_pressure_without_u(ieos)) ieos = 2
 if (tdyn > 0.) then
    idamp = 2
    tdyn_s = tdyn*utime
 else
    idamp = 1
    damp = 0.05
 endif
 if (gr) mass1 = 0. ! use Minkowski metric during relaxation

end subroutine set_options_for_relaxation
!----------------------------------------------------------------
!+
!  check if a previous snapshot exists of the relaxed star
!+
!----------------------------------------------------------------
subroutine check_for_existing_file(filename,npart,mgas,xyzh,vxyzu,restart,ierr)
 use dump_utils, only:open_dumpfile_r,read_header,dump_h,lenid,extract
 use fileutils,  only:getnextfilename
 use io,         only:idump,idisk1,id,nprocs,iprint,iverbose
 use readwrite_dumps, only:read_dump
 character(len=*), intent(inout) :: filename
 integer,          intent(in)    :: npart
 real,             intent(in)    :: mgas
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:)
 logical,          intent(out)   :: restart
 integer,          intent(out)   :: ierr
 logical :: iexist
 character(len=len(filename)) :: restart_file,filetmp
 character(len=lenid) :: fileid
 type(dump_h) :: hdr
 integer :: npartfile,iverboseprev
 real :: hfactfile,tfile,mfile
 !
 ! check for the last file in the list relax_00000, relax_00001 etc
 !
 ierr = 0
 iexist = .true.
 filetmp = filename
 restart_file = ''
 restart = .false.
 do while (iexist)
    inquire(file=filetmp,exist=iexist)
    if (iexist) then
       restart_file = filetmp
       filetmp = getnextfilename(filetmp)
    endif
 enddo
 if (len_trim(restart_file) <= 0) return

 print "(/,1x,a,/)",'>> RESTARTING relaxation from '//trim(restart_file)
 iverboseprev = iverbose
 call open_dumpfile_r(idump,restart_file,fileid,ierr)
 call read_header(idump,hdr,ierr)
 close(idump)
 if (ierr /= 0) then
    print "(a)",' ERROR: could not read file header'
    return
 else
    call extract('nparttot',npartfile,hdr,ierr,default=0)
    if (npartfile /= npart) then
       print "(a,i0,a,i0)",' ERROR: np=',npartfile,' in '//&
         trim(restart_file)//' differs from current np=',npart
       ierr = 2
       return
    else
       call extract('massoftype',mfile,hdr,ierr,default=0.)
       if (abs(mfile-mgas) > epsilon(0.)) then
          print "(a,es10.3,a,es10.3)",' ERROR: M=',npart*mfile,' in '//&
            trim(restart_file)//' differs from current M=',npart*mgas
          ierr = 3
          return
       else
          iverbose = -1
          restart = .true.
          call read_dump(restart_file,tfile,hfactfile,&
                         idisk1,iprint,id,nprocs,ierr)
          filename = restart_file
       endif
    endif
 endif
 iverbose = iverboseprev

end subroutine check_for_existing_file

!--------------------------------------------------
!+
!  get mass coordinate m(r) = \int 4.*pi*rho*r^2 dr
!+
!--------------------------------------------------
function get_mr(rho,r) result(mr)
 use physcon, only:pi
 real, intent(in) :: rho(:),r(:)
 real :: mr(size(r))
 integer :: i

 mr(1) = 0.
 do i=2,size(rho)
    mr(i) = mr(i-1) + 4./3.*pi*rho(i) * (r(i) - r(i-1)) * (r(i)**2 + r(i)*r(i-1) + r(i-1)**2)
 enddo

end function get_mr

!----------------------------------------------------------------
!+
!  restore previous settings
!+
!----------------------------------------------------------------
subroutine restore_original_options(i1,npart)
 use eos,     only:ieos,gamma
 use damping, only:damp
 use options, only:idamp
 use part,    only:hfact,vxyzu,gr
 use externalforces, only:mass1
 integer, intent(in) :: i1,npart

 gamma = gammaprev
 hfact = hfactprev
 ieos  = ieos_prev
 idamp = 0
 damp = 0.
 vxyzu(1:3,i1+1:npart) = 0.
 if (gr) mass1 = mass1prev

end subroutine restore_original_options

!----------------------------------------------------------------
!+
!  write relaxation options to .setup file
!+
!----------------------------------------------------------------
subroutine write_options_relax(iunit,relax,ieos)
 use infile_utils, only:write_inopt
 use dim,          only:do_radiation
 use radiation_utils, only:exchange_radiation_energy,limit_radiation_flux,tol_rad,itsmax_rad,cv_type,kappa_cgs
 use eos,          only:iopacity_type
 integer, intent(in) :: iunit,ieos
 logical, intent(in) :: relax

 call write_inopt(tol_ekin,'tol_ekin','tolerance on ekin/epot to stop relaxation',iunit)
 call write_inopt(maxits,'maxits','maximum number of relaxation iterations',iunit)
 if (do_radiation .and. relax .and. ieos==20) then
    call write_inopt(relax_radiation,'relax_radiation','relax with radiation transport after setting composition',iunit)
    call write_inopt(damp_radrelax,'damp_radrelax','damp velocities during radiative relaxation',iunit)
   call write_inopt(tol_radrelax,'tol_radrelax',&
                'tolerance on mass-weighted mean absolute change in u+xi per dynamical time (code specific energy)',iunit)
    call write_inopt(maxits_radrelax,'maxits_radrelax','maximum number of steps in radiative relaxation',iunit)
    write(iunit,"(/,a)") '# radiation options for relaxation (override .in values during setup)'
    call write_inopt(exchange_radiation_energy,'gas-rad_exchange','exchange energy between gas and radiation',iunit)
    call write_inopt(limit_radiation_flux,'flux_limiter','limit radiation flux',iunit)
    call write_inopt(iopacity_type,'iopacity_type','opacity method (0=inf,1=mesa,2=constant,-1=preserve)',iunit)
    if (iopacity_type==2) call write_inopt(kappa_cgs,'kappa_cgs','constant opacity value in cm2/g',iunit)
    call write_inopt(tol_rad,'tol_rad','tolerance on backwards Euler implicit solve of dxi/dt',iunit)
    call write_inopt(itsmax_rad,'itsmax_rad','max number of iterations for radiation implicit solve',iunit)
    call write_inopt(cv_type,'cv_type','temperature conversion (0=constant,1=MESA,20=gas+rec)',iunit)
 endif

end subroutine write_options_relax

!----------------------------------------------------------------
!+
!  read relaxation options from .setup file
!+
!----------------------------------------------------------------
subroutine read_options_relax(db,nerr,relax,ieos)
 use infile_utils, only:inopts,read_inopt
 use dim,          only:do_radiation
 use radiation_utils, only:exchange_radiation_energy,limit_radiation_flux,tol_rad,itsmax_rad,cv_type,kappa_cgs
 use eos,          only:iopacity_type
 type(inopts), allocatable, intent(inout) :: db(:)
 integer,                   intent(inout) :: nerr
 logical,                   intent(in)    :: relax
 integer,                   intent(in)    :: ieos

 call read_inopt(tol_ekin,'tol_ekin',db,errcount=nerr)
 call read_inopt(maxits,'maxits',db,errcount=nerr)
 if (do_radiation .and. relax .and. ieos==20) then
    call read_inopt(relax_radiation,'relax_radiation',db,errcount=nerr)
    call read_inopt(damp_radrelax,'damp_radrelax',db,errcount=nerr)
    call read_inopt(tol_radrelax,'tol_radrelax',db,errcount=nerr,min=0.)
    call read_inopt(maxits_radrelax,'maxits_radrelax',db,errcount=nerr,min=1)
    ! Require explicit radiation settings when enabled. Missing keys trigger
    ! the normal setup-file rewrite, making the effective defaults visible.
    if (relax_radiation) then
       call read_inopt(exchange_radiation_energy,'gas-rad_exchange',db,errcount=nerr)
       call read_inopt(limit_radiation_flux,'flux_limiter',db,errcount=nerr)
       call read_inopt(iopacity_type,'iopacity_type',db,errcount=nerr,min=-1,max=2)
       if (iopacity_type==2) call read_inopt(kappa_cgs,'kappa_cgs',db,errcount=nerr,min=0.)
       call read_inopt(tol_rad,'tol_rad',db,errcount=nerr,min=epsilon(tol_rad))
       call read_inopt(itsmax_rad,'itsmax_rad',db,errcount=nerr,min=1)
       call read_inopt(cv_type,'cv_type',db,errcount=nerr,min=0,max=20)
    endif
 endif

end subroutine read_options_relax

end module relaxstar
