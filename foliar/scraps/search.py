import edge_orient, find_orient
import snappy
import snappy.snap.t3mlite as t3m
import util
from collections import Counter
import regina
import random

def edge_orientation_stats(manifold):
    if isinstance(manifold, str):
        T = t3m.Mcomplex(manifold)
    else:
        T = manifold
    ans = dict()
    if len(T.Vertices) == 1 and T.Vertices[0].link_genus() == 0:
        orient = list(edge_orient.edge_orientations(T))
        ans['num_orient'] = len(orient)
        ans['sutures'] = Counter(eo.num_sutures() for eo in orient)
        ans['num_sinks'] = Counter(eo.num_super_long_edges() for eo in orient)
        ans['fol'] = Counter(eo.gives_foliation() for eo in orient)
    return ans


def first_foliation(snappy_manifold, max_triangulations=10):
    for iso in util.closed_isosigs(snappy_manifold)[:max_triangulations]:
        T = t3m.Mcomplex(iso)
        T.name = iso
        if len(T.Vertices) == 1 and T.Vertices[0].link_genus() == 0:
            orient = list(edge_orient.edge_orientations(T))
            for eo in orient:
                if eo.gives_foliation():
                    return eo

def try_persistent(snappy_manifold):
    Y = snappy_manifold
    n = len(Y.dual_curves())
    for i in range(n):
        M = Y.drill(i).filled_triangulation()
        if M.solution_type() == 'all tetrahedra positively oriented': 
            slopes = edge_orient.degeneracy_slopes_with_search(M)
            print slopes

def tri_supports_foliation(iso):
    T = t3m.Mcomplex(iso)
    for eo in edge_orient.edge_orientations(T):
        if eo.gives_foliation():
            return eo
        


        
def examine_two_vertex(snappy_manifold, rand_max, max_size):
    M = snappy_manifold
    for iso in util.closed_isosigs(M, rand_max, max_size):
        for i in range(1):
            R = regina.NTriangulation(iso)
            interesting_two_vertex_triangulation(R)
            new_iso = R.isoSig()
            T = t3m.Mcomplex(new_iso)
            T.name = new_iso
            print(len(T))
            for eo in edge_orient.edge_orientations(T):
                if eo.gives_foliation():
                    print(M, new_iso, eo.signs)
                    return eo

def basic_examine_two_vertex(iso):
    R = regina.NTriangulation(iso)
    interesting_two_vertex_triangulation(R)
    new_iso = R.isoSig()
    T = t3m.Mcomplex(new_iso)
    print new_iso
    for eo in edge_orient.edge_orientations(T):
        if eo.gives_foliation():
            print(M, new_iso, eo.signs)
            return 

                

if __name__ == '__main__':
    #iso = 'oLLLALAPzQcbedgfjihlkmnlnnxxnxxqaxaxqqnhx'
    #R = regina.NTriangulation(iso)
    #M = snappy.Manifold('m146(6,1)')
    #bad = 'pLALvAzAQzQabcegiklijkmnooobrwbbclhjjnjhxoj_cDbC'
    M = snappy.Manifold('s137(5, 4)')
    #for ex in zhs:
    #    M = snappy.Manifold(ex)
    #    examine_two_vertex(M)

    onevert = 'zLALLzvLLQALzMPQMkcbbeghlpqkpnmtruwsvxxxvyyyqqjnajwahjhcnqrehsqqaqwaet'
    twovert = 'zLLLvwMwMwzPMMLQQkbcgikkllolspttqvrwuyuyxxxytsmjxpiiagpftixavwrakpaxat'

    #Winner for 'm146(6, 1)' sLLLzvQzPPwQQacfglhimkoplnprqrqrnkwkqexxbuavubupipj
    
