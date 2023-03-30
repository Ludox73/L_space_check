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


def save_proof(task):
    name = task['name']
    manifold = snappy.Manifold(name)
    radius = task['cayley_radius']
    if manifold.solution_type() == 'all tetrahedra positively oriented':
        success, proof = quickdisorder.has_non_orderable_group(manifold,
                                                               silent=True,
                                                               ball_radius=radius, return_proof = True)
        assert success
        file = open('proofs/' + name, 'w')
        file.write(proof + '\n')
        file.close()
    task['done'] = True

task = {'name':'m003(-3, 1)', 'cayley_radius':4}
taskdb2.worker.run_function('QHSpheres', 'task_save_proofs', save_proof)
