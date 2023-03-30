#! /bin/env sage-python
#
#SBATCH --partition m
#SBATCH --tasks=1
#SBATCH --mem-per-cpu=3000
#SBATCH --nice=10000
#SBATCH --time=7-00:00
#SBATCH --output=slurm_out/%j
#SBATCH --error=slurm_error/%j

import taskdb2, taskdb2.worker
import snappy, spherogram, edge_orient, util, peripheral

def search_for_persist(task):
    for isosig in util.cusped_isosigs(task['name']):
        M = peripheral.Triangulation(isosig)
        try:
            for eo in edge_orient.edge_orientations(M):
                if eo.gives_foliation():
                    task['foliar_tri'] = isosig
                    task['foliar_orient'] = repr(eo.signs).replace(' ', '')
                    task['done'] = True
                    return
        except AssertionError:
            return 

def add_alex(task):
    M = snappy.Triangulation(task['name'])
    alex = M.alexander_polynomial()
    task['alex'] = repr(alex).replace(' ', '')
    coeffs = alex.coefficients()
    if (coeffs != (len(coeffs)//2)*[1, -1] + [1] or coeffs == [1] or
        (coeffs == [1, -1, 1] and len(M.link()) != 3)):
        task['floer_simple'] = -1
    task['done'] = True

def nonalt_data_for_paper():
    """
    (1210608,
    492,
    12,
    ['K8n3',
    'K10n21',
    'K12n242',
    'K13n4587',
    'K13n4639',
    'K14n6022',
    'K14n21881',
    '15n40211',
    '15n41185',
    '15n124802',
    '16n184868',
    '16n783154'])
    """
    import sage.all
    import pandas
    R = sage.all.PolynomialRing(sage.all.ZZ, 'a')
    db = taskdb2.ExampleDatabase('KnotsInS3')
    df = db.dataframe()
    df = df[df.name.apply(lambda x:'n' in x)]
    for db_name in ['KnotsInS3NonAlt15','KnotsInS3NonAlt16']:
        db = taskdb2.ExampleDatabase(db_name)
        df = pandas.concat([df, db.dataframe()], ignore_index=True)

    alex_compat = 0
    for i, row in df.iterrows():
        alex = R(row.alex)
        coeffs = alex.coefficients()
        if (coeffs != (len(coeffs)//2)*[1, -1] + [1] or coeffs == [1] or
            (coeffs == [1, -1, 1] and row['name'] != 'K3a1')):
            pass
        else:
            alex_compat += 1
    exceptions = df[df.foliar_orient.isnull()]
    return len(df), alex_compat, len(exceptions), list(exceptions.name)


def cable_trefoil(n):
    """
    Follows Figure 1 of Hedden "On knot Floer homology and cabling"
    """
    
    braid = 3*[2, 1, 3, 2] + 6*[-1] + n*[1]
    K = spherogram.Link(braid_closure=braid)
    K.simplify()
    return K

def brute_check():
    to_check = [('K8n3', 'T(3, 4)'),
                ('K10n21', 'T(3, 5)'),
                ('K14n21881', 'T(3, 7)'),
                ('15n41185', 'T(4, 5)'),
                ('16n783154', 'T(3, 8)'),
                ('13n4587', 'cable_trefoil(7)'),
                ('13n4639', 'cable_trefoil(5)'),
                ('15n40211', 'cable_trefoil(9)'),
                ('15n124802', 'cable_trefoil(3)'),
                
    ]
    for knot, desc in to_check:
        M = snappy.Triangulation(knot)
        if desc[0] == 'T':
            K = spherogram.Link(desc)
        else:
            K = eval(desc)
        D = K.exterior().without_hyperbolic_structure()
        while M != D:  # checks for isomorphic triangulations
            D.randomize()
            if D == M:
                break
            M.randomize()
        print('Found isomorphism between %s and %s' % (knot, desc))
            
        
task1 = {'name':'K10a12'}
task2 = {'name':'K7a7'}
#search_for_persist(task1)


#taskdb2.worker.run_function('KnotsInS3', 'task_alex', add_alex)
#taskdb2.worker.run_function('KnotsInS3NonAlt16', 'task_persist', search_for_persist)
