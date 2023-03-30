#! /bin/env sage-python
#
#SBATCH --partition m
#SBATCH --tasks=1
#SBATCH --mem-per-cpu=3000
#SBATCH --nice=10000
#SBATCH --time=7-00:00
#SBATCH --output=slurm_out/%j
#SBATCH --error=slurm_error/%j

import snappy
import taskdb2.worker
import quickdisorder
import sys, os, json
sys.path.append('../check_proof')
import check_proof

dir = 'proofs/'
if not os.path.exists(dir):
    dir = '/pkgs/tmp/proofs/'
    assert os.path.exists(dir)

def check_one(task):
    name = task['name']
    proof = json.loads(open(dir + name, 'r').read())
    success, bits = check_proof.check_proof_harder(proof)
    if success:
        task['nonord_pf_prec'] = bits
    task['done'] = True

task = {'name':'m003(-3, 1)'}
taskdb2.worker.run_function('QHSpheres', 'task_check_proofs', check_one)
