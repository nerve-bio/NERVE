#!/usr/bin/env python3
import argparse
import os
import subprocess
import time


parser = argparse.ArgumentParser()
parser.add_argument('--root', required=True)
parser.add_argument('--image', required=True)
parser.add_argument('--output', required=True)
args = parser.parse_args()

name = f'nerve-checkpoint-validation-{os.getpid()}'
checkpoint = os.path.join(args.output, 'checkpoint.pkl')
container_output = '/workdir/' + os.path.relpath(args.output, args.root)
command = [
    'docker', 'run', '--rm', '--name', name, '--network', 'nerve-network',
    '-v', f'{args.root}:/workdir', args.image,
    '--proteome1', '/workdir/tests/input_data/proteome1.fasta', '--gram', 'n',
    '--annotation', 'False', '--razor', 'False', '--mouse', 'False',
    '--virulent', 'False', '--select', 'True', '--epitopes', 'False',
    '--ep_plots', 'False', '--threads', '1', '--epitope_workers', '1',
    '--keep_snapshots', 'True', '--checkpoint', 'True',
    '--working_dir', container_output,
]

process = subprocess.Popen(command)
deadline = time.monotonic() + 120
try:
    while time.monotonic() < deadline:
        if process.poll() is not None:
            raise SystemExit('Pipeline exited before a checkpoint could be interrupted')
        if os.path.isfile(checkpoint):
            subprocess.run(['docker', 'stop', name], check=True,
                           stdout=subprocess.DEVNULL)
            break
        time.sleep(0.05)
    else:
        raise SystemExit('Timed out waiting for a completed checkpoint')
    try:
        process.wait(timeout=20)
    except subprocess.TimeoutExpired:
        subprocess.run(['docker', 'stop', name], check=False,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        process.wait(timeout=20)
finally:
    if process.poll() is None:
        subprocess.run(['docker', 'stop', name], check=False,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        process.wait()

if process.returncode == 0:
    raise SystemExit('Pipeline completed instead of being interrupted')
print('Interrupted after a completed checkpoint was written.')
