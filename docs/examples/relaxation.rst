
Stellar relaxation
==================


When mapping 1D profile from MESA to a 3D code such as Phantom, a relaxation procedure is required to ensure that the particles are in hydrostatic equilibrium.

In phantom this is achieved by a process called ``relax-o-matic``
(`Lau et al., 2022 <https://doi.org/10.1093/mnras/stac049>`_ - **please see its appendix C for full details of how it works**).


To run ``relax-o-matic`` after setting up a star,
just follow the instructions during setup (from ``phantomsetup``) and type ``yes`` when prompted "Relax star automatically during setup?";

Then, modify the ``.setup`` file and fill in the relaxation options, which are pretty self-explanatory::

	relax_star        =           T    ! relax star automatically during setup
	tol_ekin          =   1.000E-07    ! tolerance on ekin/epot to stop relaxation
	tol_dens          =       1.000    ! % error in density to stop relaxation
	maxits            =        1000    ! maximum number of relaxation iterations
	write_rho_to_file =           F    ! write density profile to file

Run the ``phantomsetup`` again to start the auto relaxation process.
You are now good to go!

If interruptted during the relaxation, run ``phantomsetup`` again and it will pick up from where it left off automatically.

If in doubt, you can always run the star in isolation for a few years to see if relaxation procedure worked.

Radiative relaxation
--------------------

For a radiation build with ``ieos=20`` and ``relax=T``, the ``.setup`` file
also supports a second relaxation stage after composition and thermal energies
have been assigned::

    relax_radiation = T
    damp_radrelax = F
    tol_radrelax = 1.e-3
    maxits_radrelax = 10000
    gas-rad_exchange = T
    flux_limiter = T
    iopacity_type = 1
    tol_rad = 1.e-6
    itsmax_rad = 250
    cv_type = 0

This stage evolves hydrodynamics and radiation using the implicit radiation
solver. Its timestep initialisation also sets the sink/external-force
integrator's kick and drift coefficients, as in normal evolution. Without
these coefficients, particles can acquire velocities without moving, while
radiation pressure work removes energy from a spatially frozen model.
The radiation options in ``.setup`` override values loaded from
``.in`` during setup. When ``relax_radiation=T``, missing radiation options
trigger the normal setup-file rewrite, inserting the inherited values or
program defaults for review before relaxation starts. Run setup again after
reviewing the updated file. These are shared radiation settings, so an input file subsequently
written by setup also contains them. The effective radiation settings are
printed when radiative relaxation starts. For constant opacity
(``iopacity_type=2``), also specify ``kappa_cgs`` in cm2/g.

``gas-rad_exchange=F`` disables thermal energy exchange between gas and
radiation, while radiation diffusion, radiation pressure work and hydrodynamics
remain active. ``cv_type`` controls the conversion from gas specific energy to
temperature in the radiation solver: 0 uses constant heat capacity, 1 uses
the MESA EoS, and 20 uses the gas/recombination EoS. With ``cv_type=20`` and
variable composition, the implicit solver uses each particle's hydrogen
and metal fractions for both its energy-to-temperature conversion and the
temperature supplied to the opacity lookup. The opacity tables themselves
remain initialized for the global composition. A constant conversion
(``cv_type=0``) does not generally reproduce EoS 20 temperatures.

``tol_radrelax`` measures the larger of the mass-averaged absolute changes in
gas and radiation specific energies, scaled to one dynamical time (in code
specific-energy units). It is not a fractional energy tolerance.
Convergence is checked every 10 steps, and snapshots are normally written
every 100 steps. Thus ``rad_relax1_00001`` usually represents step 100.
Existing radiative relaxation snapshots are automatically reused on restart;
compare changed physics in a fresh run directory to avoid resuming an already
cooled model.

