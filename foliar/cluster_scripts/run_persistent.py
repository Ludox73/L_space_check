#! /bin/env sage-python
#
#SBATCH --partition m
#SBATCH --tasks=1
#SBATCH --mem-per-cpu=3000
#SBATCH --nice=10000
#SBATCH --time=7-00:00
#SBATCH --output=slurm_out/%j
#SBATCH --error=slurm_error/%j


import taskdb2, snappy, edge_orient

def create_database():
    examples = [M.name() for M in snappy.OrientableCuspedCensus(cusps=1)]
    columns = [('degeneracy_slopes', 'text'), ('persistent_triangulations', 'text')]
    db = taskdb2.ExampleDatabase('persistent_cusped', examples, columns)
    db.new_task_table('task_find_laminar')
    return db

def search_for_persistent(task):
    M = snappy.Triangulation(task['name'])
    slopes, tris = edge_orient.degeneracy_slopes_with_search(M)
    task['degeneracy_slopes'] = repr(slopes)
    task['persistent_triangulations'] = repr(tris)
    task['done'] = True

exdb = taskdb2.ExampleDatabase('persistent_cusped')
exdb.run_function('task_find_laminar', search_for_persistent)




