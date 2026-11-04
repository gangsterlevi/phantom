#!/usr/bin/env python3
"""Measure EOS changes independently against a committed three-tier baseline.

Example (after building radbinary):
  python3 scripts/bench_eos_optimizations.py --baseline-ref c7d68f465 \
      --snapshot /home/levi/rad_heian/sim_00088 --output /tmp/eos_bench

Compiles isolated variants with the production flags, changes no working-tree
sources, pins timing runs to one CPU, and reports median CPU time from five
interleaved runs. Snapshot warm guesses are deliberately offset by +/-1%.
The generated stall cohort tests early escalation separately from normal calls.
Requires numpy and the dump reader alongside the supplied snapshot.
"""
import argparse
import ast
from collections import Counter
import json
import math
import os
from pathlib import Path
import shlex
import statistics
import subprocess

import numpy as np

DRIVER = r'''
program benchmark
 use eos, only:get_cv
 use eos_gasradrec, only:equationofstate_gasradrec,calc_uP_from_rhoT_gasradrec,init_eos_gasradrec,&
                        nstall_newton,nfallback_quad,bench_iter1,bench_iter2,irecomb
 use units, only:set_units,unit_density,unit_ergg
 use physcon, only:solarm,solarr
 implicit none
 real,allocatable::state(:,:)
 real::T,p,cs,gamma,imu,cv,t0,t1,checksum,rhocode,ucode,u
 integer::n,i,j,reps,ierr,unit,n2,n3,warm
 character(len=256)::mode,input,arg,operation
 call get_command_argument(1,mode)
 call get_command_argument(2,input)
 call get_command_argument(3,arg);read(arg,*)warm
 call get_command_argument(4,arg);read(arg,*)reps
 call get_command_argument(5,operation)
 if(trim(operation)=='solver')irecomb=1
 call set_units(mass=solarm,dist=solarr,G=1.d0)
 call init_eos_gasradrec(ierr)
 if(trim(mode)=='generate') then
    open(newunit=unit,file=trim(input),status='replace')
    write(unit,*)10500
    ! Floating-point neighbours of a previously observed tier-2 failure.
    do i=1,500
       u=1.4174741629268076e-15+real(mod(i,7)-3)*spacing(1.4174741629268076e-15)
       cv=10545316870190.395+real(i-250)*spacing(10545316870190.395)
       write(unit,'(5(es25.17,1x))')u,u*cv,1.,0.69843,0.28731
    enddo
    do i=1,100
       do j=1,100
          T=10.**(2.+6.*real(j-1)/99.)
          u=10.**(-18.+24.*real(i-1)/99.)
          call calc_uP_from_rhoT_gasradrec(u,T,0.69843,0.28731,cv,p,imu,do_radiation=.true.)
          write(unit,'(5(es25.17,1x))')u,u*cv,1.,0.69843,0.28731
       enddo
    enddo
    close(unit)
    stop
 endif
 open(newunit=unit,file=trim(input),status='old')
 read(unit,*)n
 allocate(state(5,n))
 do i=1,n
    read(unit,*)state(:,i)
 enddo
 close(unit)
 checksum=0.
 if(trim(mode)=='classify')then
    do i=1,n
       n2=nstall_newton;n3=nfallback_quad;bench_iter1=0;bench_iter2=0
       call evaluate(i,T)
       write(*,'(5(i0,1x),es25.17)') i,nstall_newton-n2,nfallback_quad-n3,bench_iter1,bench_iter2,T
    enddo
 else
    do i=1,n
       call evaluate(i,T)
    enddo
    call cpu_time(t0)
    do j=1,reps
       do i=1,n
          call evaluate(i,T)
          checksum=checksum+T
       enddo
    enddo
    call cpu_time(t1)
    write(*,'(a,es25.17,1x,a,i0,1x,a,es25.17)') 'seconds=',t1-t0,'calls=',n*reps,'checksum=',checksum
 endif
contains
subroutine evaluate(k,temp)
 integer,intent(in)::k
 real,intent(out)::temp
 if(trim(operation)=='solver')then
    temp=state(3,k)
    call equationofstate_gasradrec(state(1,k),state(2,k),temp,imu,state(4,k),state(5,k),&
                                  p,cs,gamma,do_radiation=.true.)
 else
    rhocode=state(1,k)/unit_density
    ucode=(state(2,k)/state(1,k))/unit_ergg
    if(warm==1)then
       cv=get_cv(20,rhocode,ucode,X_local=state(4,k),Z_local=1.-state(4,k)-state(5,k),&
                 temperature_guess=state(3,k))
    else
       cv=get_cv(20,rhocode,ucode,X_local=state(4,k),Z_local=1.-state(4,k)-state(5,k))
    endif
    temp=ucode/cv
 endif
end subroutine
end program
'''


def checked(args, **kwargs):
    result = subprocess.run(args, text=True, capture_output=True, **kwargs)
    if result.returncode:
        raise RuntimeError(f"Command failed: {shlex.join(map(str, args))}\n{result.stdout}\n{result.stderr}")
    return result.stdout


def without_stagnation(source):
    start_marker = '       ! An exact fixed point or two-cycle of the complete iteration state'
    while start_marker in source:
        start = source.index(start_marker)
        end = source.index('       state_prev2 = state_prev', start)
        source = source[:start] + source[end:]
    progress_marker = '       ! Double precision has stopped improving in the local Newton'
    while progress_marker in source:
        start = source.index(progress_marker)
        end = source.index('       endif\n', start) + len('       endif\n')
        source = source[:start] + source[end:]
    return source


def full_output_cv(source):
    a = source.index('real function get_cv(')
    b = source.index('end function get_cv', a)
    body = source[a:b].replace('only:get_temperature_gasradrec', 'only:equationofstate_gasradrec')
    body = body.replace('temp,X,Z,mu,', 'temp,imu,pres_cgs,cs_cgs,gamma_eff,X,Z,mu,')
    old = '    call get_temperature_gasradrec(rho_cgs,rho_cgs*u_cgs,temp,X,1.-X-Z,do_radiation=do_radiation)'
    assert old in body
    body = body.replace(old, '    call equationofstate_gasradrec(rho_cgs,rho_cgs*u_cgs,temp,imu,X,1.-X-Z,pres_cgs,&\n'
                             '                                   cs_cgs,gamma_eff,cveff,do_radiation=do_radiation)')
    return source[:a] + body + source[b:]


def baseline_optional_argument(source):
    # ABI adapter only: the baseline continues to ignore the temperature hint.
    return source.replace('gamma_local) result(cv)', 'gamma_local,temperature_guess) result(cv)', 1).replace(
        ' integer, intent(in) :: cv_type\n',
        ' integer, intent(in) :: cv_type\n real, intent(in), optional :: temperature_guess\n', 1)


def instrument(source):
    source = source.replace('\ncontains\n', '\n public :: nstall_newton,nfallback_quad\n'
                            ' integer, public :: bench_iter1=0,bench_iter2=0\ncontains\n', 1)
    return source.replace('    do n = 1,nmax\n', '    do n = 1,nmax\n'
                          '       if (use_quad) then\n'
                          '          bench_iter2=n\n'
                          '       else\n'
                          '          bench_iter1=n\n'
                          '       endif\n')


def read_snapshot(path, count):
    # Load only the existing binary-reader functions, without plotting imports.
    namespace = {'np': np, 'Counter': Counter}
    reader = path.parent / '1d_phantom_profiles.py'
    names = {'read_binary', 'read_fortran_record', 'read_string', 'read_names',
             'read_record_with_names', 'read_phantom_dump'}
    for node in ast.parse(reader.read_text()).body:
        if isinstance(node, ast.FunctionDef) and node.name in names:
            exec(compile(ast.Module(body=[node], type_ignores=[]), str(reader), 'exec'), namespace)
    dump = namespace['read_phantom_dump'](path)
    data, units, quantities = dump['blocks'][0]['data'], dump['units'], dump['quantities']
    indices = np.linspace(0, len(data['u']) - 1, count, dtype=int)
    rho = float(np.asarray(quantities['massoftype']).flat[0]) * float(units['umass']) * (
        float(quantities['hfact']) / (np.asarray(data['h'], dtype=float)[indices] * float(units['udist'])))**3
    u = np.asarray(data['u'])[indices] * (float(units['udist']) / float(units['utime']))**2
    x = np.asarray(data['H fraction'])[indices]
    y = 1. - x - np.asarray(data['metallicity'])[indices]
    guess = np.asarray(data['temperature'])[indices] * np.where(np.arange(count) % 2 == 0, .99, 1.01)
    return np.column_stack([rho, rho * u, guess, x, y])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-ref', required=True)
    parser.add_argument('--snapshot', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--samples', type=int, default=1000)
    parser.add_argument('--repetitions', type=int, default=5)
    parser.add_argument('--seconds', type=float, default=.25, help='Target CPU seconds per timed batch')
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    flags = shlex.split((repo / 'build/.make_lastfflags').read_text())
    driver = root / 'benchmark.f90'
    driver.write_text(DRIVER)
    base_gas = checked(['git', 'show', args.baseline_ref + ':src/main/eos_gasradrec.f90'], cwd=repo)
    base_eos = checked(['git', 'show', args.baseline_ref + ':src/main/eos.f90'], cwd=repo)
    gas = (repo / 'src/main/eos_gasradrec.f90').read_text()
    eos = (repo / 'src/main/eos.f90').read_text()
    variants = {
        'baseline': (base_gas, baseline_optional_argument(base_eos), 0),
        'warm_only': (base_gas, full_output_cv(eos), 1),
        'temperature_only': (without_stagnation(gas), eos, 0),
        'early_only': (gas, baseline_optional_argument(base_eos), 0),
        'warm_temperature': (without_stagnation(gas), eos, 1),
        'all': (gas, eos, 1),
    }
    executables = {}
    for name, (gas_src, eos_src, warm) in variants.items():
        folder = root / name
        folder.mkdir(exist_ok=True)
        (folder / 'eos_gasradrec.f90').write_text(instrument(gas_src))
        (folder / 'eos.f90').write_text(eos_src)
        for source in ['eos_gasradrec.f90', 'eos.f90']:
            checked(['gfortran', *flags, '-I' + str(repo / 'build'), '-c', source], cwd=folder)
        exe = folder / 'benchmark'
        checked(['gfortran', *flags, '-I' + str(repo / 'build'), str(driver),
                 'eos_gasradrec.o', 'eos.o', str(repo / 'bin/libphantom.a'), '-o', str(exe)], cwd=folder)
        executables[name] = (exe, warm)
    cpu = min(os.sched_getaffinity(0))

    def invoke(name, mode, path, reps=1, operation='cv'):
        exe, warm = executables[name]
        return checked(['taskset', '-c', str(cpu), str(exe), mode, str(path), str(warm), str(reps), operation],
                       cwd=exe.parent, env={**os.environ, 'OMP_NUM_THREADS': '1'})

    def save(name, states):
        path = root / (name + '.dat')
        np.savetxt(path, states, fmt='%.17e', header=str(len(states)), comments='')
        return path

    def classify(name, path, operation='cv'):
        out = invoke(name, 'classify', path, operation=operation)
        rows = [line.split() for line in out.splitlines() if len(line.split()) == 6 and line.split()[0].isdigit()]
        (root / (path.stem + '-' + name + '.txt')).write_text(out)
        return np.array([[int(r[1]), int(r[2]), int(r[3]), int(r[4]), float(r[5])] for r in rows])

    def timing(name, path, reps, operation):
        line = next(line for line in invoke(name, 'time', path, reps, operation).splitlines()
                    if line.startswith('seconds='))
        fields = line.replace('=', ' ').split()
        return 1e6 * float(fields[1]) / int(fields[3])

    report = {'baseline_ref': args.baseline_ref, 'snapshot': str(args.snapshot), 'flags': flags,
              'compiler': checked(['gfortran', '--version']).splitlines()[0], 'cpu': cpu, 'workloads': {}}
    normal = save('snapshot', read_snapshot(args.snapshot, args.samples))
    normal_reference = classify('baseline', normal)
    slow_tier1 = normal_reference[:, 2] >= 400
    plateau_path = None
    if np.any(slow_tier1):
        plateau_path = save('tier1_plateau', np.loadtxt(normal, skiprows=1)[slow_tier1])
    candidates = root / 'stall_candidates.dat'
    invoke('baseline', 'generate', candidates, operation='solver')
    records = classify('baseline', candidates, operation='solver')
    stalled = records[:, 1] > 0
    stall_path = None
    if np.any(stalled):
        stall_path = save('stalled', np.loadtxt(candidates, skiprows=1)[stalled][:100])
    print(f'Found {int(sum(stalled))} genuine tier-2 failures in {len(records)} states.', flush=True)
    for workload, path, names, operation in [
        ('snapshot', normal, list(variants), 'cv'),
        ('tier1_plateau', plateau_path, ['baseline', 'early_only'], 'cv'),
        ('stalled', stall_path, ['baseline', 'early_only'], 'solver'),
    ]:
        if path is None:
            continue
        reference = classify('baseline', path, operation)
        entries = {}
        repeats = {}
        for name in names:
            data = reference if name == 'baseline' else classify(name, path, operation)
            delta = float(np.max(np.abs(data[:, 4] / reference[:, 4] - 1.)))
            if not np.all(np.isfinite(data[:, 4])) or delta > 5e-14:
                raise RuntimeError(f'{workload}/{name}: temperature disagreement {delta}')
            entries[name] = {'tier2': int(sum(data[:, 0])), 'tier3': int(sum(data[:, 1])),
                             'max_tier1_iterations': int(max(data[:, 2])),
                             'total_tier1_iterations': int(sum(data[:, 2])),
                             'max_tier2_iterations': int(max(data[:, 3])), 'max_relative_T_difference': delta,
                             'samples_us': []}
            first = timing(name, path, 1, operation)
            repeats[name] = max(1, math.ceil(args.seconds * 1e6 / (first * len(data))))
        for trial in range(args.repetitions):
            order = names[trial % len(names):] + names[:trial % len(names)]
            for name in order:
                entries[name]['samples_us'].append(timing(name, path, repeats[name], operation))
        for name, entry in entries.items():
            entry['median_us'] = statistics.median(entry['samples_us'])
            entry['speedup'] = statistics.median(entries['baseline']['samples_us']) / entry['median_us']
            print(workload, name, json.dumps(entry), flush=True)
        report['workloads'][workload] = {'n': len(reference), 'operation': operation, 'variants': entries}
        (root / 'results.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
