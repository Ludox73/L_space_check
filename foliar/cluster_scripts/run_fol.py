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
import snappy, foliar
import snappy.snap.t3mlite as t3m
import edge_orient
import search

def search_for_taut(task):
    for D in eval(task['descriptions']):
        M = snappy.Manifold(D)
        fol = foliar.first_foliation(M, 1000, 40)
        if fol is not None:
            task['taut'] = True
            task['laminar_tri'] = fol.mcomplex.name
            task['done'] = True
            return

def search_for_taut_2_vertex(task):
    M = snappy.Manifold(task['name'])
    fol = search.examine_two_vertex(M, 1000, 40)
    if fol is not None:
        task['taut'] = True
        task['laminar_tri'] = fol.mcomplex.name
        task['done'] = True
        return 

def count_taut(task):
    N = t3m.Mcomplex(str(task['laminar_tri']))
    laminar_orients = [eo for eo in edge_orient.edge_orientations(N) if eo.gives_foliation()]
    task['laminar_orients'] = repr([eo.signs for eo in laminar_orients]).replace(' ', '')
    task['taut_euler_0'] = repr([1 if eo.euler_class_vanishes() else 0 for eo in laminar_orients]).replace(' ', '')
    task['done'] = True

task1 = {'name':'m003(-1, 3)', 'laminar_tri':'jLLvMQQcdfigihghihsafroggnw',
         'descriptions':"['m003(-1, 3)', 'm003(-2, 3)', 'm004(-5, 1)', 'm004(5, 1)', 'm011(1, 2)', 'm015(3, 1)', 'm019(2, 1)', 's000(3, 1)', 's467(1, 1)', 's586(1, 1)', 's885(1, 1)', 's932(-1, 1)', 'v1315(0, 1)', 'v2302(1, 1)', 'v2334(-1, 1)', 'v2344(1, 1)', 'v2542(-1, 1)', 'v2848(-1, 1)', 'v3052(1, 1)', 'v3266(1, 1)', 'v3274(-1, 1)', 't00448(-2, 1)', 't06645(1, 1)', 't10297(1, 1)', 't10708(-1, 1)', 't10885(0, 1)', 't11312(0, 1)', 't11456(1, 1)', 't11688(0, 1)', 't11764(1, 0)', 't11897(0, 1)', 't12059(0, 1)', 't12615(-1, 1)', 't12786(1, 1)', 'o9_03815(0, 1)', 'o9_18326(1, 1)', 'o9_19336(-1, 1)', 'o9_20865(1, 1)', 'o9_22656(-1, 1)', 'o9_24484(1, 1)', 'o9_24890(-1, 1)', 'o9_28719(0, 1)', 'o9_28915(-1, 1)', 'o9_29719(1, 1)', 'o9_31461(1, 1)', 'o9_31677(1, 0)', 'o9_32886(0, 1)', 'o9_33515(1, 1)', 'o9_34039(0, 1)', 'o9_35072(-1, 1)', 'o9_36122(0, 1)', 'o9_36194(0, 1)', 'o9_36640(0, 1)', 'o9_37670(1, 1)', 'o9_37757(0, 1)', 'o9_37994(-1, 1)', 'o9_38895(1, 0)', 'o9_39767(0, 1)', 'o9_39821(0, 1)', 'o9_40290(1, 0)', 'o9_41249(0, 1)', 'o9_41302(1, 0)', 'o9_42355(1, 0)', 'o9_42457(0, 1)', 'o9_42715(1, 0)', 'o9_43725(0, 1)', 'o9_43880(0, 1)', 'o9_43943(1, 1)', 'o9_44080(-1, 1)']"}
task2 = {'name':'o9_34819(5, 1)', 'laminar_tri':'uLALvvLPMQvAQQccbbeilkjpknmqoprrtsstqqnnbmxeonkvtngpfrkkk'}
    
#taskdb2.worker.run_function('QHSpheres', 'task_fol', search_for_taut_2_vertex)
