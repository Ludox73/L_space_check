import nsagetools, foliar, snappy, edge_orient
T = snappy.RationalTangle

def pretzel_link(a, b, c):
    tangle = T(1, a) + T(1, b) + T(1, c)
    link = tangle.numerator_closure()
    return link.exterior()

def some_dehn_fillings(M, num_fillings):
    lattice = nsagetools.NormalizedCuspLattice(M)
    L = 5
    while True:
        slopes = lattice.primitive_elements(L)
        if len(slopes) > num_fillings:
            break
        L = 1.5*L
    return slopes[:num_fillings]

def other_dehn_fillings(M, max_denom):
    d = M.alexander_polynomial().degree()
    slopes = []
    for i in range(d, 3*d):
        slopes += [(i, 1), (-i, 1)]
    return slopes

def have_taut_foliation(M, fillings):
    ans = []
    for slope in fillings:
        N = M.copy()
        N.dehn_fill(slope)
        F = foliar.first_foliation(N)
        if F is not None:
            print(slope)
            ans.append(slope)
    return ans

pretzels = [(-2, 5, 5), (-2, 5, 7), (-2, 7, 7)]

for p in pretzels:
    M = pretzel_link(*p)
    ans = edge_orient.degeneracy_slopes_with_search(M)
    print(p, ans)




