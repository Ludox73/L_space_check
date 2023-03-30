#! /bin/env sage-python
#
#SBATCH --partition m
#SBATCH --tasks=1
#SBATCH --mem-per-cpu=3000
#SBATCH --nice=10000
#SBATCH --time=7-00:00
#SBATCH --output=slurm_out/%j
#SBATCH --error=slurm_error/%j

import taskdb2, snappy, foliar

def search_for_taut(task):
    M = snappy.Manifold(task['isosig'])
    fol = foliar.first_foliation(M, 40)
    if fol is not None:
        task['taut'] = True
        task['laminar_tri'] = fol.mcomplex.name
        task['done'] = True
        
exdb = taskdb2.ExampleDatabase('foliate_double')
exdb.run_function('task_fol', search_for_taut, columns=['isosig'])
