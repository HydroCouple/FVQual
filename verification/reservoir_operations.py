#!/usr/bin/env python3
"""Qualify scheduled operations and compare FVQual CLI with Composer outputs.

Run after importing examples/swmm_reservoir/reservoir.2dm as reservoir.mesh.json.
Artifacts stay in tests/output; no third-party Python packages are needed.
"""
import argparse
import copy
import csv
import json
import math
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def read_rows(path):
    rows = {}
    for row in csv.DictReader(path.open()):
        key = (float(row['time']), row['item'], int(row['index']))
        assert key not in rows, key
        rows[key] = float(row['value'])
    return rows


def close(actual, expected, tolerance):
    assert math.isfinite(actual) and abs(actual - expected) <= tolerance, (actual, expected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--composer', type=Path, required=True)
    parser.add_argument('--fvqual', type=Path, default=ROOT/'build/apps/fvqual')
    parser.add_argument('--components', type=Path, default=ROOT/'build/components')
    parser.add_argument('--output', type=Path, default=ROOT/'tests/output/reservoir_operations_phase4')
    args = parser.parse_args()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    example = ROOT/'examples/reservoir_operations/composition.json'
    spec = json.loads(example.read_text())
    spec['components'][0]['arguments']['mesh']['@uri'] = str((example.parent/spec['components'][0]['arguments']['mesh']['@uri']).resolve())
    for host in ('cli', 'composer'):
        document = copy.deepcopy(spec)
        document['writers'][0]['path'] = host+'/fields.csv'
        document['run']['manifest'] = host+'/run.json'
        path = args.output/(host+'.composition.json')
        path.write_text(json.dumps(document, indent=2)+'\n')
        command = ([str(args.fvqual.resolve()), 'run', str(path)] if host == 'cli' else
                   [str(args.composer.resolve()), '--components', str(args.components.resolve()), '--run', str(path)])
        result = subprocess.run(command, env=dict(os.environ, QT_QPA_PLATFORM='offscreen'), capture_output=True, text=True, timeout=120)
        (args.output/(host+'.log')).write_text(result.stdout+result.stderr)
        assert result.returncode == 0, result.stdout+result.stderr
        manifest = json.loads((args.output/host/'run.json').read_text())
        assert all(c['final_status'] == 'Finished' for c in manifest['components'])
    cli = read_rows(args.output/'cli/fields.csv')
    gui_host = read_rows(args.output/'composer/fields.csv')
    # The CLI uses elapsed seconds; Composer records absolute Julian days.
    # Infer the origin for this numerical parity check; the separate clock
    # tests qualify epoch semantics. Both hosts use 7-second exchanges.
    times = sorted({t for t, _, _ in gui_host})
    cli_times = sorted({t for t, _, _ in cli})
    assert len(times) == len(cli_times)
    epoch = times[0]-cli_times[0]/86400
    for t, elapsed in zip(times, cli_times):
        close((t-epoch)*86400, elapsed, 5e-5)
    mapped = {(cli_times[times.index(t)], item, index): value for (t, item, index), value in gui_host.items()}
    assert cli.keys() == mapped.keys()
    difference = max(abs(value-mapped[key]) for key, value in cli.items())
    assert difference == 0, difference
    max_residual = {}
    for time in cli_times:
        v = {item: value for (t, item, index), value in cli.items() if t == time and index == 0}
        # Independently integrate rates and the deliberately irregular knots.
        source = .01*min(time,33.5)+.02*max(time-33.5,0)
        release = .0005*min(time,40.5)+.001*max(time-40.5,0)
        boundary = .0001*10*time  # known 10 m exterior edge
        close(v['budget_water_source_in'], source, 1e-10)
        close(v['budget_water_structure_out'], release, 1e-10)
        close(v['budget_water_boundary_in'], boundary, 1e-9)
        close(v['water_volume'], 1000+source-release+boundary, 1e-8)
        close(v['budget_tracer_source_in'], source*.01, 1e-12)
        close(v['budget_heat_source_in'], source*18*1000*4186, .01)
        close(v['budget_heat_shortwave'], 200*.94*(100*min(time,33.5)+200*max(time-33.5,0)), .001)
        close(v['budget_particles_settling']+v['budget_sediment_settling'], 0, 1e-10)
        for quantity in ('water','heat','tracer','reactive','particles','sediment'):
            residual = abs(v['budget_'+quantity+'_residual'])
            close(residual,0,.1 if quantity == 'heat' else 1e-8)
            max_residual[quantity] = max(max_residual.get(quantity,0),residual)
        for quantity in ('water','heat','tracer','reactive','particles'):
            close(v['budget_'+quantity+'_transfer_in'],v['budget_'+quantity+'_transfer_out'],.01 if quantity == 'heat' else 1e-12)
    report = {'status':'passed','compared_values':len(cli),'maximum_host_difference':difference,
              'maximum_residuals':max_residual,'irregular_control_knots_s':[33.5,40.5],
              'cli_manifest':str(args.output/'cli/run.json'),'composer_manifest':str(args.output/'composer/run.json')}
    (args.output/'qualification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
