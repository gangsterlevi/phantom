
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

Constant-entropy core replacement
---------------------------------

The constant-entropy core option adjusts the sink mass to match the enclosed
mass at the replacement boundary. A bounded bisection solve requires the
residual central gas mass to be less than ``1.e-10`` times the gas mass in
the softened region. This is a mass-closure tolerance, not a measure of the
radial profile's discretisation error. The reconstructed mass profile includes
the innermost gas shell. If the solve fails, setup reports an error rather
than using an unconverged profile.

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

``tol_radrelax`` limits the fractional absolute change in combined gas and
radiation energy per dynamical time. For each particle, the changes in gas
and radiation specific energies are added before taking the absolute value;
the mass-weighted sum is divided by the total gas plus radiation energy at
initialisation and multiplied by ``tdyn/dtcheck``. Thus gas--radiation exchange
cancels locally, while energy redistribution between particles does not.
The reference energy is fixed for each invocation, including a resumed run.
A tolerance of ``1.e-3`` corresponds to 0.1 per cent of that reference energy
per dynamical time; this replaces the previous tolerance in code
specific-energy units.

Convergence requires this energy error below ``tol_radrelax``. The ratio
``Ekin/abs(Epot)`` is monitored using the same global energy calculation as
hydrodynamic relaxation, but does not restrict radiative convergence;
``tol_ekin`` applies only to hydrodynamic relaxation. There is no minimum
elapsed-time requirement. Both diagnostics are printed to the console and
written to the radiative relaxation ``.ev`` file (columns 3 and 4). Checks occur every 10
steps, and snapshots are normally written every 100 steps. Thus
``rad_relax1_00001`` usually represents step 100.
Existing radiative relaxation snapshots are automatically reused on restart;
compare changed physics in a fresh run directory to avoid resuming an already
cooled model.

