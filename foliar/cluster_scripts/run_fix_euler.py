#! /bin/env sage-python
#
#SBATCH --partition m
#SBATCH --tasks=1
#SBATCH --mem-per-cpu=3000
#SBATCH --nice=10000
#SBATCH --time=7-00:00
#SBATCH --output=slurm_out/%j
#SBATCH --error=slurm_error/%j

import taskdb2.worker
import snappy
import snappy.snap.t3mlite as t3m
from edge_orient import EdgeOrientation

def compute_euler(task):
    N = t3m.Mcomplex(str(task['foliar_tri']))
    signs = eval(task['foliar_orients'])
    orients =  [EdgeOrientation(N, s) for s in signs]
    assert all(eo.gives_foliation() for eo in orients)
    task['vertices_foliar'] = len(N.Vertices)
    if len(N.Vertices) == 1:
        task['new_taut_euler_0'] = repr([1 if eo.euler_class_vanishes() else 0
                                     for eo in orients]).replace(' ', '')
    task['done'] = True


task1 = {'name':'m376(-2, 3)',
         'foliar_tri':'oLLvLMLPQQcacgikmkimjlnnnmnkjaaagnnnwkwkw',
         'foliar_orients':'[[1,1,-1,1,1,-1,-1,1,1,1,1,-1,1,1,-1]]'}

taskdb2.worker.run_function('QHSpheres', 'task_double_check_euler', compute_euler)
