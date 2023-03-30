import snappy
import snappy.snap.t3mlite as t3m
import regina

def closed_from_isosig(isosig):
    """
    Function not used, but kept around for doctests

    >>> isosig = 'jLLvQPQcdfhghigiihshhgfifme'
    >>> N = t3m.Mcomplex(isosig)
    >>> len(N)
    9
    >>> N.snappy_manifold().identify()
    [m003(-2,3)]
    >>> N.isosig() == isosig
    True
    """
    return t3m.Mcomplex(isosig)
        
def closed_isosigs(snappy_manifold, trys=20, max_tets=50):
    """
    >>> M = snappy.Manifold('m004(1,2)')
    >>> len(closed_isosigs(M, trys=5)) > 0
    True
    """
    M = snappy_manifold.copy()
    assert M.cusp_info('complete?') == [False]
    surgery_descriptions = [M]

    for curve in M.dual_curves():
        N = M.drill(curve)
        N.dehn_fill((1,0), 1)
        surgery_descriptions.append(N.filled_triangulation([0]))

    ans = set()
    for N in surgery_descriptions:
        for i in range(trys):
            T = N.filled_triangulation()
            if T._num_fake_cusps() == 1:
                n = T.num_tetrahedra()
                if n <= max_tets:
                    ans.add((n, T.triangulation_isosig(decorated=False)))
            N.randomize()

    return [iso for n, iso in sorted(ans)]

def cusped_triangulations(snappy_manifold, trys=1000):
    """
    >>> M = snappy.Manifold('m004')
    >>> len(cusped_triangulations(M, trys=100))
    1
    """
    M = snappy.Triangulation(snappy_manifold)
    ans = [M.copy()]
    for i in range(trys):
        M.randomize()
        if {len(M.isomorphisms_to(N)) > 0 for N in ans} == {False}:
            ans.append(M.copy())
    return ans

def cusped_isosigs(snappy_manifold, trys=1000):
    """
    >>> M = snappy.Manifold('m004')
    >>> len(list(cusped_isosigs(M, trys=100)))
    1
    """
    M = snappy.Triangulation(snappy_manifold)
    seen = set()
    for i in range(trys):
        isosig = M.triangulation_isosig()
        isobase = isosig.split('_')[0]
        if isobase not in seen:
            seen.add(isobase)
            yield(isosig)
        M.randomize()



# Code for trying to create interesting 2 vertex triangulations.  
            
def highest_degree_vertex(regina_tri):
    return max(regina_tri.vertices(), key=lambda v:v.degree())

def triangles_on_high(regina_tri):
    T = regina_tri
    v = highest_degree_vertex(T)
    for F in T.triangles():
        if F.vertex(0) == F.vertex(1) == F.vertex(2) == v:
            return F

def obvious_simplify(regina_tri):
    T = regina_tri
    any_progress = False
    progress = True
    while progress:
        progress = False
        edges = sorted(T.edges(), key=lambda e:e.degree())
        for e in edges:
            d = e.degree()
            if d == 1:
                progress = T.twoOneMove(e, 0)
            elif d == 2:
                progress = T.twoZeroMove(e)
            elif d == 3:
                progress = T.threeTwoMove(e)
            else:
                break

            if progress:
                any_progress = True
                break
    return any_progress

def simplify_via_randomization(regina_tri, max_failed_attempts=100):
    T = regina_tri
    any_progress = obvious_simplify(T)
    failures = 0
    while failures < max_failed_attempts:
        edges = [e for e in T.edges() if e.degree() == 4]
        if len(edges) == 0:
            return 
        edge = random.choice(edges)
        T.fourFourMove(edge, random.choice([0, 1]))
        progress = obvious_simplify(T)
        if progress:
            any_progress = True
        else:
            failures += 1

def interesting_two_vertex_triangulation(regina_tri):
    T = regina_tri
    tets = T.tetrahedra()
    for tet in tets:
        T.oneFourMove(tet)

    face = triangles_on_high(T)
    while not face is None:
        T.twoThreeMove(face)
        face = triangles_on_high(T)

    progress = True
    while progress and T.countVertices() > 2:
        progress = False
        edges = sorted(T.edges(), key=lambda e:e.vertex(0).degree() + e.vertex(1).degree())
        for edge in edges:
            if T.collapseEdge(edge):
                progress = True
                break

    obvious_simplify(T)
    
    
if __name__ == '__main__':
    import doctest
    doctest.testmod()
