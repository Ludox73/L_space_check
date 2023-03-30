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

def printer(depth, string):
    pass

quickdisorder.printer = printer

def try_to_nonorder(task):
    manifold = snappy.Manifold(task['name'])
    radius = task['cayley_radius']
    if radius is None:
        radius = 4
    else:
        radius += 1
    ans = quickdisorder.has_non_orderable_group(manifold, ball_radius = radius, 
                                           fundamental_group_args = [True, True, False])
    if ans:
        task['orderable'] = -1
        task['done'] = True
    task['cayley_radius'] = radius

task = {'name':'v0123(2,3)', 'cayley_radius':3}
#taskdb2.worker.run_function('QHSpheres', 'task_disorder', try_to_nonorder)
