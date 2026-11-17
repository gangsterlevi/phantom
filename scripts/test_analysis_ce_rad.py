#!/usr/bin/env python3
"""Regress EOS 20 CE analysis with synthetic full/small radiation dumps.

Build first with ``make SETUP=radbinary NOWARN=yes phantomanalysis``. This
test uses that build's compiler flags and libraries to write its fixtures.
All generated files are kept in a temporary directory.
"""
import argparse
import math
import os
from pathlib import Path
import re
import resource
import shlex
import subprocess
import tempfile


def run_checks(repo, build, binary, compiler):
    env = {**os.environ, 'OMP_NUM_THREADS': '1', 'PHANTOM_DIR': str(repo)}
    _, hard = resource.getrlimit(resource.RLIMIT_STACK)
    resource.setrlimit(resource.RLIMIT_STACK, (hard, hard))
    flags = shlex.split((build / '.make_lastfflags').read_text())
    with tempfile.TemporaryDirectory(prefix='phantom-ce-rad-') as tmp:
        work = Path(tmp)
        fixture = work / 'write_test_ce_dumps'
        command = [compiler, *flags, '-I', str(build),
                   str(repo / 'src/tests/write_test_ce_dumps.f90'),
                   str(build.parent / 'bin/libphantom.a'), '-o', str(fixture)]
        for invocation in (command, [str(fixture)]):
            p = subprocess.run(invocation, cwd=work, env=env, capture_output=True, text=True)
            assert p.returncode == 0, p.stdout + p.stderr
        reference = [float(x) for x in (work / 'ce_reference.txt').read_text().split()]

        # First file is deliberately small: it must not consume prompts or state.
        dumps = [work / 'test_00001', work / 'test_00000', work / 'test_00001', work / 'test_00002']
        for choice, eos in ((21, 20), (3, 20), (26, 20), (21, 2), (3, 2)):
            p = subprocess.run([str(binary), *map(str, dumps)], cwd=work, env=env,
                               input=f'{choice}\nno\n{eos}\n\n', text=True,
                               capture_output=True, timeout=60)
            output = p.stdout + p.stderr
            assert p.returncode == 0, output
            assert 'FATAL ERROR' not in output, output
            assert output.count('common-envelope analysis requires a full dump') == 2, output
            assert output.count('Choose analysis type') == 1, output
            if choice == 21:
                matches = re.findall(r'including radiation\)\s*\n([^\n]+)', output)
                assert len(matches) == 2, output
                rows = [[float(x) for x in row.split()] for row in matches]
                # Fixture references are in code units; printed values are erg.
                # The writer provides the unit conversion alongside the energies.
                energy_unit = float((work / 'ce_unit.txt').read_text())
                expected_binding = reference[:3] if eos == 20 else [reference[0], reference[2], reference[2]]
                for row in rows:
                    assert len(row) == 3, row
                    for actual, expected in zip(row, expected_binding):
                        assert math.isfinite(actual)
                        assert math.isclose(actual / energy_unit, expected, rel_tol=2e-6), (row, reference)
            elif choice == 3:
                rows = [[float(x) for x in line.split()]
                        for line in (work / 'energy.ev').read_text().splitlines()
                        if line.strip() and not line.lstrip().startswith('#')]
                assert len(rows) == 2, rows
                for row in rows:
                    assert all(math.isfinite(x) for x in row), row
                    assert math.isclose(row[1], reference[3], rel_tol=2e-6), (row, reference)
                    assert math.isclose(row[-1], reference[4], rel_tol=2e-6), (row, reference)
            else:
                rows = [[float(x) for x in line.split()]
                        for line in (work / 'erec_vs_t.ev').read_text().splitlines()
                        if line.strip() and not line.lstrip().startswith('#')]
                assert len(rows) == 2, rows
                for row in rows:
                    assert math.isclose(sum(row[1:]), reference[2] - reference[1], rel_tol=2e-6), (row, reference)
            print(f'PASS: CE option {choice}: EOS {eos}, separate radiation, composition, dead particles, small dumps')


def main():
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, default=repo / 'build')
    parser.add_argument('--binary', type=Path, default=repo / 'bin/phantomanalysis')
    parser.add_argument('--compiler', default='gfortran')
    args = parser.parse_args()
    run_checks(repo, args.build.resolve(), args.binary.resolve(), args.compiler)


if __name__ == '__main__':
    main()
