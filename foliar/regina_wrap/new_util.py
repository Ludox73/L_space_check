"""
This module uses the "pyregina" interface from the "regina_wrap"
subdirectory to try to get a broader collection of triangulations to
examine for foliar orientations.
"""

import snappy
import snappy.snap.t3mlite as t3m
import pyregina, regina
import foliar

def to_iso(snappy_manifold):
    M = snappy_manifold
    return (M.num_tetrahedra(), M.triangulation_isosig(decorated=False))
        
def fancy_closed_isosigs(snappy_manifold, height, print_starts=True):
    """
    >>> M = snappy.Manifold('m004(1,2)')
    >>> len(fancy_closed_isosigs(M,1, print_starts=False)) > 0
    True
    """
    M = snappy_manifold.copy()
    starts = set()

    for i in range(10):
        starts.add(to_iso(M.filled_triangulation()))
        for curve in M.dual_curves():
            N = M.drill(curve)
            N.dehn_fill((1,0), 1)
            starts.add(to_iso(N.filled_triangulation()))
        M.randomize()

    min_tets = min(starts)[0]
    max_tets = min_tets + height
    starts = sorted(N for N in starts if N[0] <= max_tets)
    isosigs = set()
    
    for n, iso in starts:
        if (n, iso) not in isosigs:
            isosigs.add((n, iso))
            if print_starts:
                print('Starting massive search from %s' % iso)
            T = pyregina.Triangulation(iso)
            isosigs.update(T.retriangulate(max_tets - n))
        else:
            print('Already found %s' % iso)
    return sorted(isosigs)

def two_vertex_tris(snappy_manifold, height=0):
    M = snappy_manifold.filled_triangulation()
    T = regina.NTriangulation(M._to_string())
    t = T.tetrahedra()[0]
    T.oneFourMove(t)
    isosig = T.isoSig()
    P = pyregina.Triangulation(isosig)
    isosigs = P.retriangulate(height)
    return isosigs

def first_foliation(snappy_manifold, height=0, print_starts=True):
    """
    >>> M = snappy.Manifold('m004(1,2)')
    >>> eo = first_foliation(M, height=0, print_starts=False)
    >>> eo.gives_foliation()
    True
    """
    sigs = fancy_closed_isosigs(snappy_manifold, height, print_starts)
    if print_starts:
        print('Found %d distinct triangulations' % len(sigs))
    for n, iso in sigs:
        T = t3m.Mcomplex(iso)
        T.name = iso
        if len(T.Vertices) == 1 and T.Vertices[0].link_genus() == 0:
            orient = list(foliar.edge_orientations(T))
            for eo in orient:
                if eo.gives_foliation():
                    return eo

def first_persisent(snappy_manifold, height=0):
    sigs = fancy_closed_isosigs(snappy_manifold, height)
    print('Found %d distinct triangulations' % len(sigs))
    for n, iso in sigs:
        ans = edge_orient.is_persistent_tri(iso)
        if ans is not None:
            return ans

                
if __name__ == '__main__':
    import doctest
    print(doctest.testmod())
