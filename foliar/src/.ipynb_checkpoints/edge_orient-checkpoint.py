from . import find_orient, link, dual_cellulation, surface, peripheral, util
import snappy
import snappy.snap.t3mlite as t3m
from snappy.snap.t3mlite.simplex import (Head, Tail,
                                         ZeroSubsimplices, OneSubsimplices,
                                         RightFace, LeftFace)
import networkx as nx

class EdgeOrientation(object):
    """
    An orientation on the edges of a triangulation of a closed
    3-manifold where no face is a directed cycle.

    >>> N = t3m.Mcomplex('jLLvQPQcdfhghigiihshhgfifme')
    >>> orients = edge_orientations(N)
    >>> good = [eo for eo in orients if eo.num_sutures() == 1]
    >>> len(good)
    8
    >>> [eo(e) for e in N.Edges] == eo.signs
    True
    >>> tet = N.Tetrahedra[3]
    >>> [eo.is_very_long(tet, e) for e in OneSubsimplices]
    [True, False, False, False, False, False]
    >>> eo.has_sink_edge()
    False
    >>> eo.gives_foliation()
    True
    >>> {eo.euler_class_vanishes() for eo in good}
    set([True])
    >>> T = t3m.Mcomplex('tLLLLMLLwPMQPkacfihjinmlpmoqrpsrssjkgqqthqkwtvxofsqcaa')
    >>> orients = list(edge_orientations(T))
    >>> foliations = [eo for eo in orients if eo.gives_foliation()]
    >>> len(orients), len(foliations)
    (55, 30)
    >>> {eo.euler_class_vanishes() for eo in foliations}
    set([False, True])

    Sanity checks:

    >>> EdgeOrientation(T, len(T.Edges)*[1])
    Traceback (most recent call last):
       ...
    ValueError: EdgeOrientation is not acyclic
    >>> EdgeOrientation(t3m.Mcomplex('m004'), len(T.Edges)*[1])
    Traceback (most recent call last):
       ...
    ValueError: Some vertex link is not a sphere

    Examples where there are more than one vetex.

    >>> R = t3m.Mcomplex('nLLLwAPLQkcdefhhihklmlmmhsdarkdjselaxj')
    >>> {eo.strongly_connected() for eo in edge_orientations(R)}
    set([False])

    >>> Y = t3m.Mcomplex('uLLLMvPzvwPQQQcacfgghimlrnonspsqtttsjkwwjmbxwagoronokvkwr')
    >>> orients = list(edge_orientations(Y))
    >>> foliations = [eo for eo in orients if eo.gives_foliation()]
    >>> len(orients), len(foliations)
    (115, 0)
    >>> Z = t3m.Mcomplex('sLLLvLLLQAPQQcdghmljnpmlrqoqoprrhshvxuulhrrptftvgpk')
    >>> orients = list(edge_orientations(Z))
    >>> foliations = [eo for eo in orients if eo.gives_foliation()]
    >>> len(orients), len(foliations)
    (202, 116)
    """
    def __init__(self, mcomplex, signs, check=True):
        self.mcomplex, self.signs = mcomplex, signs
        if not all(v.link_genus()==0 for v in mcomplex.Vertices):
            raise ValueError('Some vertex link is not a sphere')
        if check:
            self._check_acyclic()
        self._setup_local_structure()
        

    def _check_acyclic(self):
        for face in find_orient.oriented_edges_around_faces(self.mcomplex):
            face_signs = []
            for e in face:
                s = 1 if e > 0 else -1
                face_signs.append(s*self.signs[abs(e) - 1])
            if face_signs in ([1, 1, 1], [-1, -1, -1]):
                raise ValueError('EdgeOrientation is not acyclic')

    def _setup_local_structure(self):
        """
        The "local_structure[tet, vertex]" attribute records the number of
        "out" and "in" arrows of the vertex in the given tetrahedra.
        """

        self.local_structure = local_structure = dict()
        for tet in self.mcomplex.Tetrahedra:
            for a in ZeroSubsimplices:
                out_arrows, in_arrows = 0, 0
                for (edge, sign) in tet.edge_info[a]:
                    if sign * self.signs[edge] > 0:
                        out_arrows += 1
                    else:
                        in_arrows += 1
                local_structure[tet, a] = (out_arrows, in_arrows)

    def __call__(self, edge):
        """
        Return the sign self assigns to the given edge.
        """
        return self.signs[edge.Index]


    def local_structure_edge(self, tet, edge):
        loc_str = self.local_structure
        return {loc_str[tet, Head[edge]], loc_str[tet, Tail[edge]]}

    def is_very_long(self, tet_or_corner, edge=None):
        if edge is None:  # Called with a corner
            tet = tet_or_corner.Tetrahedron
            edge = tet_or_corner.Subsimplex
        else:
            tet = tet_or_corner

        for a in [Head[edge], Tail[edge]]:
            if 0 not in self.local_structure[tet, a]:
                return False
        return True
                        
    def has_sink_edge(self):
        """
        In an acyclically oriented tetrahedron, the unique edge that runs
        from the "source vertex" to the "sink vertex" is called very
        long.  A sink edge is an edge which is very long in every
        adjacent tetrahedron.  Such sink edges lead to sink discs in
        the associated branched surface.
        """
        return self.num_sink_edges() > 0

    def num_sink_edges(self):
        ans = 0
        for edge in self.mcomplex.Edges:
            if all(self.is_very_long(c) for c in edge.Corners):
                ans += 1
        return ans

    def one_skeleton_digraph(self):
        """
        Returns the 1-skeleton of the ambient 3-manifold triangulation
        oriented via self.
        """
        T = self.mcomplex
        G = nx.DiGraph()
        G.add_nodes_from(T.Vertices)
        for edge in T.Edges:
            a, b = edge.Vertices
            if self(edge) == -1:
                a, b = b, a
            G.add_edge(a, b)
        return G
            
    def strongly_connected(self):
        """
        An EdgeOrientation is strongly connected if the 1-skeleton is
        strong connected as a digraph (that is, every vertex can be
        reached from any other via a directed path) *and* every vertex
        has a loop.
        """
        G = self.one_skeleton_digraph()
        return (nx.number_strongly_connected_components(G) == 1 and 
                len(list(nx.nodes_with_selfloops(G))) == G.number_of_nodes())

    def euler_cocycle(self):
        """
        Assuming self gives a foliation, returns the values of the
        standard cocycle representing the Euler class evaluated on the
        dual face to the given edge.
        """
        cocycle = []
        for edge in self.mcomplex.Edges:
            mixed_count = 0
            for c in edge.Corners:
                tet = c.Tetrahedron
                e = c.Subsimplex
                data = {self.local_structure[tet, v] for v in [Head[e], Tail[e]]}
                if data in [{(2, 1), (0, 3)}, {(3, 0), (1, 2)}]:
                    mixed_count += 1

            assert mixed_count % 2 == 0
            val = -mixed_count/2 + 1
            cocycle.append(val * self.signs[edge.Index])
        return cocycle

    def euler_class_vanishes(self):
        # Let T be self.mcomplex and D be the dual cellulation.  Then
        # the boundary map C_2(D) -> C_1(D) is the transpose of the
        # boundary map C_2(T) -> C_1(T).  Which mean the coboundary map
        # C^1(D) -> C^2(D) is C_2(T) -> C_1(T) on the nose.
        assert len(self.mcomplex.Vertices) == 1
        d = self.mcomplex.boundary_maps()[1]
        cohomology_elem_div = d.pari.matsnf(flag=4)
        euler = t3m.linalg.Vector(self.euler_cocycle()).pari.Col()
        elem_div_after_quot_by_euler = d.pari.concat(euler).matsnf(flag=4)
        return cohomology_elem_div == elem_div_after_quot_by_euler

    def num_sutures(self):
        M = self.mcomplex
        G = nx.Graph()
        for tet in M.Tetrahedra:
            for e in OneSubsimplices:
                if self.local_structure_edge(tet, e) in [{(0,3), (1,2)}, {(2,1), (3,0)}]:
                    i = tet.Class[RightFace[e]].Index
                    j = tet.Class[LeftFace[e]].Index
                    G.add_edge(i,j)
        assert G.number_of_nodes() == len(M.Faces)
        return nx.number_connected_components(G)

    def gives_foliation(self):
        return (not self.has_sink_edge()
                and self.num_sutures() == len(self.mcomplex.Vertices)
                and self.strongly_connected())

class IdealEdgeOrientation(EdgeOrientation):
    """
    An orientation on the edges of an ideal triangulation of a
    1-cusped 3-manifold where no face is a directed cycle.

    >>> N = snappy.Triangulation('v1234')
    >>> orients = edge_orientations(N)
    >>> [eo.num_sutures() for eo in orients]
    [2, 2]
    >>> N = snappy.Triangulation('m016')
    >>> list(edge_orientations(N))
    []

    In favorable circumstances, the induced branched surface is
    laminar.  In this case, every Dehn filling except along the
    degeneracy slope gives a manifold with a co-orientable taut
    foliation.
    """
    def __init__(self, mcomplex, signs):
        self.triangulation = T = mcomplex.original_triangulation
        self.signs = signs
        self.mcomplex = mcomplex
        self.vertex_link = mcomplex.cusp_triangulation
        self.link_dual_cellulation = mcomplex.cusp_dual_cellulation
        assert len(self.mcomplex.Vertices) == 1
        assert self.mcomplex.Vertices[0].link_genus() == 1
        self._setup_local_structure()
        self._add_link_vertex_signs()

    def _add_link_vertex_signs(self):
        self.link_vertex_signs = dict()
        for vert in self.vertex_link.vertices:
            i = vert.index
            edge_sign = self.signs[abs(i) - 1]
            vert_sign = 1 if edge_sign*i > 0 else -1
            self.link_vertex_signs[vert] = vert_sign

    def link_subgraphs(self):
        pos, neg = [], []
        for vert in self.vertex_link.vertices:
            if self.link_vertex_signs[vert] > 0:
                pos.append(vert.index)
            else:
                neg.append(vert.index)
        assert {abs(p) for p in pos} == {abs(n) for n in neg}
        G = self.vertex_link.edge_graph()
        return G.subgraph(vertices=pos), G.subgraph(vertices=neg)

    def sutures(self):
        """
        Returns a list of cycles on the 1-skeleton of the dual cellulation
        of the vertex link, each representing one of the sutures
        induced by the corresponding branched surface.

        >>> M = snappy.Triangulation('m015')
        >>> orients = edge_orientations(M)
        >>> len([eo.sutures() for eo in orients])
        2
        """
        T = self.vertex_link
        D = self.link_dual_cellulation
        weights = len(T.edges) * [0]
        for tri in T.triangles:
            signs = [self.link_vertex_signs[v] for v in tri.vertices]
            if len(set(signs)) > 1:
                if sorted(signs) == [-1, 1, 1]:
                    i = signs.index(-1)
                    j = (i + 1) % 3
                    side = surface.Side(tri, (i, j))
                else:
                    assert sorted(signs) == [-1, -1, 1]
                    i = signs.index(1)
                    j = (i - 1) % 3
                    side = surface.Side(tri, (j, i))
                edge = side.edge()
                sign = edge.orientation_with_respect_to(side)
                # Edge and its dual have the same index
                weights[edge.index] = sign
        return dual_cellulation.OneCycle(D, weights).components()
                
    def link_compatible_with_foliation(self):
        sutures = self.sutures()
        assert len(sutures) % 2 == 0
        D = self.link_dual_cellulation
        return all(D.slope(suture) != 0 for suture in sutures)

    def gives_foliation(self):
        if self.has_sink_edge():
            return False

        link_ok = self.link_compatible_with_foliation()
        return link_ok

    def degeneracy_slope(self):
        """
        >>> M = snappy.Triangulation('m004')
        >>> orients = edge_orientations(M)
        >>> [eo.degeneracy_slope() for eo in orients]
        [(1, 0), (1, 0)]
        """
        assert self.gives_foliation()
        suture = self.sutures()[0]
        a, b = self.link_dual_cellulation.slope(suture)
        if a*b == 0:
            a, b = abs(a), abs(b)
        elif b < 0:
            a, b = -a, -b
        return (a, b) 

def edge_orientations(manifold):
    if isinstance(manifold, t3m.Mcomplex):  # closed manifold
        N = manifold
        for signs in find_orient.cycle_free_orientations(N):
            yield EdgeOrientation(N, signs, check=False)
    else: # 1-cusped manifold
        N = peripheral.peripheral_curve_package(manifold)[0]
        assert len(N.Vertices) == 1 and N.Vertices[0].link_genus() == 1
        for signs in find_orient.cycle_free_orientations(N):
            yield IdealEdgeOrientation(N, signs)

def degeneracy_slopes(manifold):
    """
    Finds all persistently foliar orientations on the given
    triangulation and returns a list of their degeneracy slopes.

    >>> degeneracy_slopes('m003')
    []
    >>> degeneracy_slopes('m004')
    [(1, 0)]
    """
    M = snappy.Triangulation(manifold)
    degeneracy_slopes = []
    for eo in edge_orientations(M):
        if eo.gives_foliation():
            degeneracy_slopes.append(eo.degeneracy_slope())
    return sorted(set(degeneracy_slopes))

def degeneracy_slopes_with_search(manifold, tries=1000):
    """
    >>> degeneracy_slopes_with_search('m004', tries=0)
    ([(1, 0)], ['cPcbbbiht_BaCB'])
    """
    manifold = snappy.Triangulation(manifold)
    slopes, triangulations = list(), list()
    for M in util.cusped_triangulations(manifold, tries):
        for eo in edge_orientations(M):
            if eo.gives_foliation():
                slope = eo.degeneracy_slope()
                if slope not in slopes:
                    slopes.append(slope)
                    triangulations.append(M.triangulation_isosig())
    return slopes, triangulations
    
def test_cusped(n=100, tries=1000, progress=True):
    """
    >>> test_cusped(5, 10, False)
    """
    census = snappy.OrientableCuspedCensus(cusps=1)
    for i in range(n):
        M = census.random()
        ans =  degeneracy_slopes_with_search(M, tries)
        if progress:
            print(M.name() + ' ' + repr(ans))

def has_taut_fol_with_euler_0(spec):
    N = t3m.Mcomplex(spec)
    orients = edge_orientations(N)
    good = [eo for eo in orients if eo.gives_foliation()]
    return any(eo.euler_class_vanishes() for eo in good)


def search_for_persistent(manifold, tries=10):
    """
    >>> degeneracy_slopes_with_search('m004', tries=0)
    ([(1, 0)], ['cPcbbbiht_BaCB'])
    """
    manifold = snappy.Triangulation(manifold)
    slopes, triangulations = list(), list()
    for M in util.cusped_triangulations(manifold, tries):
        for eo in edge_orientations(M):
            if eo.gives_foliation():
                return eo

def is_persistent_tri(manifold):
    """
    >>> M = snappy.Manifold('m016')
    >>> is_persistent_tri(M) is None
    True
    """
    M = snappy.Triangulation(manifold)
    for eo in edge_orientations(M):
        if eo.gives_foliation():
            return eo

if __name__ == '__main__':
    import doctest
    doctest.testmod()
    #N = t3m.Mcomplex('jLLvQPQcdfhghigiihshhgfifme')
    #N = t3m.Mcomplex('kLLLLMQkccfhijhhjijlnacshncljt')
    #orients = list(edge_orientations(N))
    #fol = [eo for eo in orients if eo.gives_foliation()]
    #eo = fol[0]
    #[eo.num_sutures() for eo in orients]
    #orients = edge_orientations(N)
    #eo = next(orients)
    #C = eo.vertex_link
    #s = eo.sutures()
    #ans = eo.link_compatible_with_foliation()

    def quick_test(N):
        for _ in range(N):
            M = snappy.OrientableClosedCensus.random()
            for filled in util.closed_isosigs(M)[:10]:
                N = t3m.Mcomplex(filled)
                [eo.num_sutures() for eo in edge_orientations(N)]
                
                
