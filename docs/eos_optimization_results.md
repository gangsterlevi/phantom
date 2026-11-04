# EOS optimization measurements

Baseline: `c7d68f465` (three-tier solver). `eoserr=1e-15` and the
finite-difference derivatives are unchanged.

## Method

GNU Fortran 15.2.0, production `SETUP=radbinary` flags, Intel Core i7-1360P,
one pinned CPU. Values are median CPU times from five interleaved batches.
The benchmark compiles each feature independently in a temporary directory;
it does not modify the production source to switch features on and off.

The ordinary workload samples 1,000 particles from `sim_00088`. Warm guesses
use each stored temperature deliberately offset by alternating -1% and +1%.
This is a controlled proxy for a recent temperature, not a measurement of a
complete simulation or a guarantee about temperature changes in every step.

## Ordinary radiation cv conversions

Times include the EOS call used to obtain `u/T`. Each single-feature row
compares against the same baseline; speedups must not be multiplied together.

| Variant | µs/call | Speedup | Tier-2 calls / 1,000 |
|---|---:|---:|---:|
| Three-tier baseline | 29.049 | 1.00× | 223 |
| Recent temperature only | 2.445 | 11.88× | 18 |
| Temperature-only solve | 25.050 | 1.16× | 223 |
| Early escalation only | 29.476 | 0.99× | 225 |
| Recent temperature + temperature-only output | 1.796 | 16.17× | 18 |
| All three changes | 1.736 | 16.73× | 20 |

None of these snapshot calls used tier 3. Temperature-only evaluation with the
same guess returned exactly the baseline temperature. Warm guesses changed
returned temperatures by at most `3.0e-15` relative, within the existing EOS
round-trip test tolerance (`2e-14`). Every solver retains its `1e-15` stopping
criterion.

## Avoiding a long tier-1 plateau

Tier 1 now switches to tier 2 when the absolute energy or pressure residual
stops decreasing near the solution. The check is active only when
`abs(correction) <= sqrt(epsilon(T))*T` (about `1.5e-8` relative in double
precision). This selects when to retry at higher precision; it does not
relax the `1e-15` convergence test, which is checked first. It needs no
additional ionization evaluation, fixed iteration cutoff, or repeat counter.

One snapshot state previously cycled among three nearby temperatures for all
500 tier-1 iterations, then needed only one tier-2 iteration. A two-cycle
detector misses this pattern. With the residual check, this state uses 22
tier-1 iterations and one tier-2 iteration. With identical analytic guesses
and the complete EOS output, its cv conversion falls from 82.078 to 23.905
microseconds (3.43× faster); the returned temperature differs by `4.44e-16`
relative.

Across the 1,000 snapshot states, the largest tier-1 count falls from 500 to
27 with early escalation alone, or 21 with all changes enabled. These are
observed counts, not iteration limits: improving solves can continue.
Early escalation sends two additional states to tier 2. It has little effect
on the ordinary workload's average cost, where quad ionization already
dominates, but removes the long tier-1 stall in the affected state.

## Actual stalled solves

A separate workload searched 10,500 states: floating-point neighbours of an
observed failure and a broad rho-T grid, with `irecomb=1`. Six states exhausted
tier 2 in the baseline. These six were timed using the complete EOS interface
with identical guesses, so temperature reuse/output omission do not confound
the early-escalation comparison.

- Baseline: 3.490 ms/call; up to 500 tier-2 iterations.
- Early escalation: 1.119 ms/call; up to 105 tier-2 iterations; 3.12× faster.
- All six still used the full quad fallback, which checks the rounded output.
- Maximum relative temperature difference: `1.11e-16`.

Tier-2 escalation detects an exact fixed point or two-cycle of the complete iteration
state `(T, Tdot, damping)` after the final scheduled damping change. It switches
to tier 3 immediately; there is no `nstagnant` counter or requirement to wait
for eight repetitions. Comparing the auxiliary state prevents a repeated
temperature during a changing update from being mistaken for a true cycle.
It never accepts a stalled state as converged.

## Compatibility and validation

- Optimized radiation build: 31/31 focused EOS and radiation groups passed.
- Bounds/FPE-checked changed modules: the same 31/31 groups passed.
- Full EOS suite: 67/69 groups passed, versus 62/64 in the baseline. The same
  two pre-existing failure counts come from the truncated Lombardi table;
  radiation-incompatible stratified tests are skipped. All five additional
  groups passed.
- A separate build with radiation disabled (`SETUP=test2`) passed 13/13
  focused EOS groups. Its general EOS suite passed 77/79 groups, with the
  same Lombardi-table failure counts.
- Constant-cv and MESA conversion tests show exactly unchanged results when
  the optional temperature hint is supplied, including a NaN hint.
- Invalid hints for `cv_type=20` fall back to the original analytic guess.
- Warm guesses are also checked at factors of 0.001 and 1,000 from the true
  temperature to cover large energy changes between calls.
- Temperature-only/full-output equivalence is tested with radiation included
  and excluded from the gas energy, for all four recombination settings.

The benchmark measures local EOS cost, not Heian wall time. The frequency of
large temperature changes and quad retries, and compiler/thread scheduling,
will determine the benefit in a full run.

## Reproduce

The optional `get_cv` argument changes its Fortran interface. Rebuild all
callers rather than reusing objects compiled against the old interface.
For an existing build, clean the compiled objects before rebuilding:

```bash
make -C build SYSTEM=gfortran SETUP=radbinary clean
make -C build SYSTEM=gfortran SETUP=radbinary phantom phantomtest
OMP_NUM_THREADS=2 bin/phantomtest eosgasradrec radiation
python3 scripts/bench_eos_optimizations.py \
  --baseline-ref c7d68f465 \
  --snapshot /path/to/rad_heian/sim_00088 \
  --output /tmp/eos_bench
```

The snapshot directory must also contain its `1d_phantom_profiles.py` binary
reader. The runner uses NumPy, records compiler flags, individual timings,
retry and iteration counts and result differences in `results.json`, and saves each isolated
source variant for inspection. The baseline revision is required explicitly.

The checked-in `eos_optimization_results.json` contains the measurements above.
The snapshot and tier-1 plateau cohorts use `irecomb=0`; the tier-2 stalled
cohort uses `irecomb=1`. No production source or input flags are altered by
the runner.
