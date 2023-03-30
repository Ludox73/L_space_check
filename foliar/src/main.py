import snappy
import snappy.snap.t3mlite as t3m
from . import util, edge_orient

def has_compatible_foliation(snappy_manifold):
    ans = first_foliation(snappy_manifold)
    return ans is not None

def first_foliation_mcomplex(mcomplex):
    for eo in edge_orient.edge_orientations(mcomplex):
        if eo.gives_foliation():
            return eo
    
def first_foliation(snappy_manifold, rand_max, max_size):
    """
    Given a SnapPy Manifold which is closed, searches for a taut
    foliation as certified as a by a foliar orientation. The
    parameters are:

    * rand_max: controls the amount of randomization used to search
      for different triangulations.

    * max_size: bounds the number of tetrahedra of any triangulation
      that will be examined in detail.

    >>> M = snappy.Manifold('m004(1, 2)')
    >>> eo = first_foliation(M, 5, 25)
    >>> eo.gives_foliation()
    True

    >>> M = snappy.Manifold('m003(-3, 1)')
    >>> eo = first_foliation(M, 5, 25)
    >>> eo is None
    True
    """
    t = 0
    for iso in util.closed_isosigs(snappy_manifold, rand_max, max_size):
        T = t3m.Mcomplex(iso)
        if len(T.Vertices) == 1 and T.Vertices[0].link_genus() == 0:
            orient = list(edge_orient.edge_orientations(T))
            for eo in orient:
                if eo.gives_foliation():
                    return eo

def nonorderable(snappy_manifold, max_triangulations=10):
    """
    A quick test for nonorderability which makes the *assumption* that
    all edges of the triangualtions are homotopy esssential.
    """
    for iso in util.closed_isosigs(snappy_manifold)[:max_triangulations]:
        T = t3m.Mcomplex(iso)
        T.name = iso
        if len(T.Vertices) == 1 and T.Vertices[0].link_genus() == 0:
            orients = edge_orient.edge_orientations(T)
            found = False
            for eo in orients:
                found = True
                break
            if not found:
                return iso
                    
            

