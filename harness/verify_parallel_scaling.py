#!/usr/bin/env python3
"""Smoke-test ProcessPool autoscaling under the container's actual CPU limits."""
import argparse, hashlib, os, time
from concurrent.futures import ProcessPoolExecutor
from code.Epitope import available_cpu_count, resolve_epitope_workers

def work(rounds):
    value = b'NERVE parallel scaling probe'
    for _ in range(rounds): value = hashlib.sha256(value).digest()
    return os.getpid(), value.hex()

def run(workers, tasks, rounds):
    started = time.monotonic()
    with ProcessPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(work, [rounds] * tasks))
    return time.monotonic() - started, len(set(pid for pid, _ in results))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tasks', type=int, default=8)
    parser.add_argument('--rounds', type=int, default=300000)
    args = parser.parse_args()
    workers = resolve_epitope_workers(0, args.tasks)
    serial, serial_pids = run(1, args.tasks, args.rounds)
    parallel, parallel_pids = run(workers, args.tasks, args.rounds)
    print(f'available_cpus={available_cpu_count()} auto_workers={workers} tasks={args.tasks}')
    print(f'serial={serial:.3f}s ({serial_pids} process) parallel={parallel:.3f}s '
          f'({parallel_pids} processes) speedup={serial / parallel:.2f}x')
    if workers > 1 and parallel_pids < 2:
        raise SystemExit('FAIL: auto mode selected multiple workers but work stayed in one process')

if __name__ == '__main__': main()
