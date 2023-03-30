import snappy
import snappy.snap.t3mlite as t3m
from sage.all import ZZ, QQ, matrix
import pandas as pd
import edge_orient, link

df = pd.read_csv('laminar_tris.csv.bz2')
df['taut_euler_0'] = df['taut_euler_0'].apply(eval)
df['laminar_orients'] = df['laminar_orients'].apply(eval)

manifolds = list(df.laminar_tri.apply(t3m.Mcomplex))
T = manifolds[0]

def test_homology(T):
    M = T.snappy_manifold()
    divs = M.homology().elementary_divisors()
    D = T.boundary_maps()[1]
    D = matrix(ZZ, D.nrows(), D.ncols(), D.list())
    other_divs = [e for e in D.elementary_divisors() if e != 1]
    return divs == other_divs

def test_cohomology(row):
    N = t3m.Mcomplex(row['laminar_tri'])
    vlink = link.LinkSphere(N)
    for signs in row['laminar_orients']:
        eo = edge_orient.EdgeOrientation(N, vlink, signs)
        print (eo.gives_foliation(), eo.euler_class_vanishes())
