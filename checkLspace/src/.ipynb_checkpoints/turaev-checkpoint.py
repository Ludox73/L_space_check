""" This code was downloaded by https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/LCYXPO and written by Nathan Dunfield. """

""" 
Computing the possible endpoints of L-space intervals of Floer simple
manifolds using: 

[RR] Jacob Rasmussen and Sarah Dean Rasmussen,
     Floer Simple Manifolds and L-Space Intervals,
     http://arxiv.org/abs/1508.05900v2

Notation matches [RR] closely and I assume familiarity with the
introduction to [RR] in the comments and code below.

Note: This should get merged into SnapPy proper, though there is the
issue of needing Heegaard to check realizability.
"""

import snappy
import heegaard
from sage.all import (ZZ, QQ, AbelianGroup, GroupAlgebra,
                      LaurentPolynomialRing, cartesian_product_iterator,
                      gcd, xgcd, lcm, floor, vector, matrix, infinity)
from snappy.snap.nsagetools import (MapToAbelianization,
                                    MapToGroupRingOfAbelianization,
                                    fox_derivative)
from collections import defaultdict
import slopes
import util

# Need a realizable presentation, that is, one comming from a Heegaard
# splitting, to apply [RR].

def realizable_presentation(manifold):
    """
    >>> G = realizable_presentation(snappy.Manifold('m004'))
    """
    G = manifold.fundamental_group()
    if heegaard.is_realizable(G.relators()):
        return G
    for T in util.cusped_triangulations(manifold):
        G = T.fundamental_group()
        if heegaard.is_realizable(G.relators()):
            return G
    G = manifold.fundamental_group(True, True, False)
    if heegaard.is_realizable(G.relators()):
        return G 
    raise ValueError('Could not find realizable presentation')

def inverse_word(word):
    return "".join(reversed(word.swapcase()))

def free_rank(abelian_group):
    return abelian_group.elementary_divisors().count(0)

def torsion_subgroup_elements(abelian_group):
    """
    >>> H = AbelianGroup([2,3,0], names='uvt')
    >>> torsion_subgroup_elements(H)
    [1, v, v^2, u, u*v, u*v^2]
    """
    A = abelian_group
    r = free_rank(A)
    gens = [g for g in A.gens() if g.order() != infinity]
    if len(gens) == 0:
        return [A.identity()]
    exponents = cartesian_product_iterator([range(g.order()) for g in gens])
    return [A(e + r*(0,)) for e in exponents]

def sum_of_torsion_elements(abelian_group):
    """
    >>> H = AbelianGroup([3,0], names='ut')
    >>> sum_of_torsion_elements(H)
    1 + u + u^2
    """
    R = GroupAlgebra(abelian_group)
    return sum(R(g) for g in torsion_subgroup_elements(abelian_group))

# Helper functions for elements of the GroupAlgebra of an abelian
# group of the form (torsion) + Z.

def t_deg(gp_ring_monomial):
    if not hasattr(gp_ring_monomial, 'exponents'):
        support =  gp_ring_monomial.support()
        assert len(support) == 1
        gp_ring_monomial = support[0]
    return gp_ring_monomial.exponents()[-1]

def t_degs(gp_ring_elt):
    return sorted([t_deg(s) for s in gp_ring_elt.terms()])

def t_deg_min(gp_ring_elt):
    return min(t_degs(gp_ring_elt))

def t_deg_max(gp_ring_elt):
    return max(t_degs(gp_ring_elt))

def collect(gp_ring_elt):
    """
    Return a dictionary whose keys are t_degrees and whose values are
    coefficents.
    """
    ans = defaultdict(lambda:0)
    for term in gp_ring_elt.terms():
        ans[t_deg(term)] += term
    return ans

def stable_t_degs(elt):
    """
    The powers of t which appear with coefficient equal to the sum of
    all the torsion elements.

    >>> H = AbelianGroup([2,0], names='ut')
    >>> A = GroupAlgebra(H)
    >>> u, t = A.gens()
    >>> stable_t_degs(1 + u + t + t**2 + u*t**2 + 3*t**3)
    set([0, 2])
    """
    R = elt.parent()
    H = R.group()
    t = H.gens_dict()['t']
    ans = set()
    f = sum_of_torsion_elements(H)
    for k, part in collect(elt).items():
        if part == f * t**k:
            ans.add(k)
    return ans

def lop_off_stable_range(elt):
    """
    >>> H = AbelianGroup([2,0], names='ut')
    >>> A = GroupAlgebra(H)
    >>> u, t = A.gens()
    >>> lop_off_stable_range(t + t**3 + (1 + u)*(1 + t**2 + t**4 + t**5))
    1 + t + t^2 + t^3 + t^4 + u + u*t^2 + u*t^4
    """
    d = t_deg_max(elt)
    assert t_deg_min(elt) >= 0
    if d == 0:
        first_stable = 0
    else:
        degs = stable_t_degs(elt)
        assert d in degs
        first_stable = max(set(range(d)) - degs) + 1
    return sum(term for term in elt.terms() if t_deg(term) <= first_stable)


def multiply_by_geometric_series(elt, x):
    """
    Multiply the GroupAlgebra element elt by (part of) the formal
    power series corresponding to 1/(1 - x) which is to say sum(x^k
    for k in range(infinity)) where x is a monomial in the
    GroupAlgebra.  The answer is cut off at t_deg equal two times that
    of elt.
    """
    n = t_deg_max(elt)
    d = t_deg(x)
    assert t_deg_min(elt) >= 0 and d > 0
    # If the following fails, then the product can't stabilize.  
    assert d <= n + 1
    k = (3*n)//d + 2
    ans = sum(x**i for i in range(k)) * elt
    return sum(term for term in ans.terms() if t_deg(term) <= 2*n)

def homological_framing(manifold, group, psi):
    alpha, beta = group.peripheral_curves()[0]

    def to_word(x, y):
        ans = ''
        for (u, word) in [(x, alpha), (y, beta)]:
            if u < 0:
                word = inverse_word(word)
            ans += abs(u)*word
        return ans
        
    l = c, d = manifold.homological_longitude()
    g, a, b = xgcd(c, d)
    m = (b, -a)
    assert m[0]*l[1] - m[1]*l[0] == 1
    u, v = [psi(g).support()[0] for g in [alpha, beta]]
    m_ab = u**m[0] * v**m[1]
    if t_deg(m_ab) < 0:
        m = (-m[0], -m[1])
        l = (-l[0], -l[1])
    m_ab = u**m[0] * v**m[1]
    l_ab = u**l[0] * v**l[1]
    m_word = to_word(*m)
    l_word = to_word(*l)
    assert psi(m_word).support()[0] == m_ab
    assert psi(l_word).support()[0] == l_ab
    assert m_ab.order() == infinity and l_ab.order() < infinity
    return m, l, m_word, l_word, m_ab, l_ab

class PartitionedInterval(object):
    """
    A partition of a closed interval [a, b], that is:

        a = x_0 < x_1 < ... < x_n = b

    >>> I = PartitionedInterval([0, 1, 3, 5, 9, 10])
    >>> -1 in I
    False
    >>> 4 in I
    True
    >>> I.subinterval(4)
    [3, 5]
    >>> I.subinterval(5)
    [5]
    >>> I.common_subinterval([3.5, 4, 4.5])
    [3, 5]
    >>> I.common_subinterval([0, 1])
    [0, 1]
    >>> I.common_subinterval([1, 1.5])
    [1, 3]
    >>> I.common_subinterval([3, 3, 3])
    [3]
    >>> I.common_subinterval([3, 9])
    Traceback (most recent call last):
    ...
    ValueError: Points are not all in a common subinterval
    >>> I.common_subinterval([0.5, 1.5])
    Traceback (most recent call last):
    ...
    ValueError: Points are not all in a common subinterval
    """
    def __init__(self, points):
        assert sorted(points) == points
        self.points = points
        self.a, self.b = points[0], points[-1]

    def __contains__(self, point):
        return self.a <= point <= self.b

    def subinterval(self, y):
        """
        Returns the subinterval of self containing y; if y is actually in
        the partition, just returns [y].
        """
        assert y in self
        if y in self.points:
            return [y]
        u = max(x for x in self.points if x <= y)
        v = min(x for x in self.points if x >= y)
        return [u, v]

    def common_subinterval(self, other_points):
        """
        Assuming all of other_points lie in a common subinterval of self,
        return said subinterval.  Raises an exception if this is not
        the case.
        """
        intervals = {tuple(self.subinterval(p)) for p in other_points}
        if len(intervals) == 1:
            return list(intervals.pop())
        proper_intervals = [I for I in intervals if len(I) == 2]
        singletons = [I[0] for I in intervals if len(I) == 1]
        if len(proper_intervals) == 1:
            I = proper_intervals[0]
            if all(s in I for s in singletons):
                return list(I)
        elif len(proper_intervals) == 0:
            if len(singletons) == 1:
                return [singletons[0]]
            if len(singletons) == 2:
                a, b = sorted(singletons)
                if {a < x < b for x in self.points} == {False}:
                    return [a, b]
        raise ValueError('Points are not all in a common subinterval')
                

class IotaInverseDtau(object):
    """
    Computing the preimage P of D_tau_plus under iota.  Since P does
    not include the longitude (0, 1), it can be regarded as a subset
    of the affine chart (1, b) in the homological framing.  While P is
    infinite, it is periodic under a Z-action corresponding to adding
    elements in the kernel of H_1(bdry M; Z) -> H_1(M; Z).

    >>> M = snappy.Manifold('m016')
    >>> tau = TuraevTorsion(M)
    >>> D = IotaInverseDtau(tau)
    >>> D.points
    [0, 1/9, 1/6, 2/9, 1/4, 1/3, 4/9, 1/2, 5/9, 2/3, 3/4, 7/9, 5/6, 8/9, 1]
    >>> D.to_chart((1, 0))
    0
    >>> D.to_chart((18, 1))
    Traceback (most recent call last):
    ...
    ValueError: Longitude not in chart
    >>> D.to_chart((0, 1))
    -1/18
    >>> D.from_chart(QQ('-1/18'))
    Slope(0, 1)
    >>> D.minimal_interval((-13, 7), (15, -8))
    [-5/9, -1/2]
    >>> QQ('-3/4') in D
    True
    >>> QQ('31/10') in D
    False
    >>> D.possible_non_L_space_cones((1, 0))
    [SlopeCone((1, 0), (9, 1)), SlopeCone((27, 1), (1, 0))]
    >>> D.possible_non_L_space_cones((0, 1))
    [SlopeCone((1, 0), (9, 1))]
    >>> D.non_L_space_cone([(1, 0), (0, 1), (1, -1)])
    SlopeCone((1, 0), (9, 1))
    >>> D
    IotaInverseDtau(L=1,m=(-1,0),l=(-18,-1),values=[(1,0),(2,0),(3,0),(4,0),(6,0),(9,0)])
    >>> E = eval(repr(D))
    >>> repr(D) == repr(E)
    True
    >>> E.non_L_space_cone([(1, 0), (0, 1), (1, -1)])
    SlopeCone((1, 0), (9, 1))
    """
    def __init__(self, tau=None, L=None, m=None, l=None, values=None):
        if tau is not None:
            L = tau.longitude_ab.order()
            m, l = tau.meridian, tau.longitude
            values = sorted(tau.D_tau_plus().values())
            
        self.L, self.meridian, self.longitude, self.values = L, m, l, values
        self._from_homological = matrix(ZZ, [self.meridian, self.longitude]).transpose()
        self._to_homological = self._from_homological**-1
        points = []
        for x, y in values:
            x = QQ(x)
            p = QQ(y)/x
            step = L/x
            while p <= L:
                points.append(p)
                p += step
        self.points = P = sorted(set(points))

        # It is convenient to extend self.points slightly to deal with
        # the periodicity of self.  Specifically, in this version we
        # can trivially look up the subinterval containing any point
        # (1, b) with b in [0, L).

        if len(P) > 0:
            if P[0] == 0:
                P = [P[-2] - L] + P
            else:
                P = [P[-1] - L] + P + [P[0] + L]
            self.partition = PartitionedInterval(P)

    def to_chart(self, slope):
        """
        For a slope in the default framing, given as a tuple or a Slope,
        return the rational number b so that slope is projectively
        equivalent to (1, b).
        """
        hom_slope = self._to_homological * vector(slope)
        if hom_slope[0] == 0:
            raise ValueError('Longitude not in chart')
        return QQ(hom_slope[1])/QQ(hom_slope[0])

    def from_chart(self, b):
        """
        For a rational number b, return the slope in the default framing
        associated to the point (1, b) in the homological framing. 
        """
        ans = self._from_homological * vector((1, b))
        ans = (ans.denominator() * ans).change_ring(ZZ)
        assert gcd(ans) == 1
        return slopes.Slope(ans)

    def shift_to_fund_domain(self, b):
        L = self.L
        return -L*floor(b/L)

    def __contains__(self, x):
        x += self.shift_to_fund_domain(x)
        return x in self.points
        
    def minimal_interval(self, x, y):
        """
        Given a pair of slopes (in the homological framing), returns a
        closed interval which

        a. contains both slopes, 

        b. has endpoints in self,

        c. contains no point of self in its interior.  

        If there is no such interval, then an exception is raised.
        """
        L, P = self.L, self.partition
        # normalize to chart (1, b) and put in natural order
        x = QQ(x[1])/QQ(x[0])
        y = QQ(y[1])/QQ(y[0])
        x, y = sorted([x, y])
        if y - x > L:
            raise ValueError('Slopes too far apart')
        shift = self.shift_to_fund_domain(x)
        x += shift
        y += shift
        assert x in P and y in P
        return [u - shift for u in P.common_subinterval([x, y])]

    def is_empty(self):
        return len(self.points) == 0

    def possible_non_L_space_cones(self, slope, non_L_slopes=None):
        """
        Given an L-space slope in the default (non-homological) framing,
        return the possible non-L-space cones as SlopeCones. Depending
        on whether slope is in self, there will be either one or two
        such cones.

        If the argument non_L_slopes is provided, will only return
        those cones which contain all of the non_L_slopes
        """
        if len(self.points) == 0:
            return [slopes.SingleSlope(self.longitude)]
        x = self.to_chart(slope)
        s = self.shift_to_fund_domain(x)
        x += s
        P = self.partition
        if x not in P.points:
            a, b = P.subinterval(x)
            u = self.from_chart(a - s)
            v = self.from_chart(b - s)
            C = slopes.SlopeCone(v, u)
            assert slope not in C
            assert self.longitude in C
            ans = [C]
        else:
            i = P.points.index(x)
            a, b = P.points[i-1], P.points[i+1]
            u = self.from_chart(a - s)
            v = self.from_chart(b - s)
            C0 = slopes.SlopeCone(slope, u)
            C1 = slopes.SlopeCone(v, slope)
            assert self.longitude in C0
            assert self.longitude in C1
            ans = [C0, C1]

        if non_L_slopes is not None:
            ans = [C for C in ans if all(slope in C for slope in non_L_slopes)]
        return ans
    
    def non_L_space_cone(self, L_space_slopes):
        """
        Given at least two L-space slopes in the default (non-homological)
        framing, return the non-L-space cones a SlopeCone.
        """
        L_space_slopes = {tuple(slope) for slope in L_space_slopes}
        if len(self.points) == 0 and len(L_space_slopes) > 0:
            return slopes.SingleSlope(self.longitude)
        
        if len(L_space_slopes) == 1:
            raise ValueError('Need two distinct slopes to determine cone')
        slope = L_space_slopes.pop()
        ans = set(self.possible_non_L_space_cones(slope))
        for slope in L_space_slopes:
            cones = set(self.possible_non_L_space_cones(slope))
            ans.intersection_update(cones)
        assert len(ans) == 1
        return ans.pop()

    def __repr__(self):
        spec = (self.L, self.meridian, self.longitude, self.values)
        ans = 'IotaInverseDtau(L=%d, m=%s, l=%s, values=%s)' % spec
        return ans.replace(' ', '')
    

class MapToPolynomialRingOfAbelianization(MapToAbelianization):
    """
    Arithmetic with group rings of (multiplicative) abelian groups is
    *extremely* slow in Sage. When computing the TuraevTorsion, there
    is a 10x speedup by working in the corresponding Laurent
    polynomial ring and then converting at the end.

    >>> M = snappy.Manifold('m003')
    >>> G = M.fundamental_group()
    >>> psi_fake = MapToPolynomialRingOfAbelianization(G)
    >>> elt = psi_fake('ab') + psi_fake('AAAAB'); elt
    u*t + u^4*t^-4
    >>> ans = psi_fake.convert_to_group_ring(elt**3); ans
    3*u*t^-2 + u^2*t^-12 + u^3*t^3 + 3*u^4*t^-7
    >>> ans.parent() == psi_fake.group_ring
    True
    """
    def __init__(self, fund_group, base_ring=ZZ):
        MapToAbelianization.__init__(self, fund_group)
        self.H = H = self._range  # The abelian group
        self.R = LaurentPolynomialRing(base_ring, H._names)
        self.group_ring = GroupAlgebra(H, base_ring)

    def __call__(self, word):
        ans = self.R(1)
        for x, e in zip(self.R.gens(), self._exponents_of_word(word)):
            ans = ans * x**e
        return ans

    def range(self):
        return self.R

    def convert_to_group_ring(self, p):
        H, R, RH = self.H, self.R, self.group_ring
        D = self.elementary_divisors
        ans = RH(0)
        for c, e in zip(p.coefficients(), p.exponents()):
            if R.ngens() == 1:
                e = [e]
            e = tuple(self._normalize_exponents(e))
            ans += c*RH.monomial(H(e))
        return ans
    
class TuraevTorsion(object):
    """
    As normalized in [RR], the Turaev torsion tau for a rational
    homology solid torus M lives in the following ring.  First, write
    H = H_1(M; Z) as Z + T where T is the torsion subgroup.  If R is
    the group ring Z[T], we will regard tau as element of R[[t]], that
    is, formal power series in non-negative powers of a variable t
    which corresponds to the (multiplicative) generator of free part
    of H.

    If we set f in R to be the sum of all elements of T, then tau has
    the property that all but finitely many powers of t^n have
    coefficient exactly f.  Internally, we use this to store tau as an
    element of the group ring Z[H]; specifically, the highest power of
    t^n appearing in the internal representation will have coefficient
    f, and all higher terms are implicitly f*t^m.

    >>> M = snappy.Manifold('m003')
    >>> tau = TuraevTorsion(M)
    """
    def __init__(self, manifold):
        self.manifold = M = manifold
        self.group = G = realizable_presentation(manifold)

        # Check the assumptions of this algorithm.
        assert manifold.num_cusps() == 1
        assert manifold.homology().elementary_divisors().count(0) == 1
        
        self.psi = psi = MapToGroupRingOfAbelianization(G)
        def phi(word):
            return t_deg(psi(word))
        self.phi = phi
        self.QH = QH = psi.range()
        self.H = H = QH.group()
        self.t = t = H.gens()[-1]
        
        # Fix a peripheral framing which is homologically natural.
        m, l, m_word, l_word, m_ab, l_ab = homological_framing(M, G, psi)
        self.meridian, self.longitude = m, l
        self.meridian_word, self.longitude_word = m_word, l_word
        self.meridian_ab, self.longitude_ab = m_ab, l_ab

        # Compute the torsion and normalize it
        rels = G.relators() + [m_word]
        gens = G.generators()
        psi_fake = MapToPolynomialRingOfAbelianization(G)
        A = matrix(psi_fake.R,
                   [[fox_derivative(R, psi_fake, g) for R in rels] for g in gens])
        d = psi_fake.convert_to_group_ring(A.det())
        if sum(d.coefficients()) < 0:
            d = -d
        d = t**(-t_deg_min(d))*d
        assert t_deg_min(d) == 0
        lowest_terms = [x for x in d.terms() if t_deg(x) == 0]
        lowest_supports = [x.support()[0] for x in lowest_terms]
        if not any(x.is_one() for x in lowest_terms):
            x = min(lowest_supports, key=lambda x:x.exponents())
            d = x**-1 * d
        tau = multiply_by_geometric_series(d, QH(m_ab))
        stable = stable_t_degs(tau)
        assert stable.issuperset(range(t_deg_max(d), t_deg_max(tau) + 1))
        self.tau = lop_off_stable_range(tau)
        assert any(x.is_one() for x in self.support())

    def could_be_floer_simple(self):
        """
        For M to be Floer simple, it is necessarily (but not sufficient)
        that all the coefficients of the torsion are 0 or 1.

        >>> M = snappy.Manifold('m015')
        >>> tau = TuraevTorsion(M)
        >>> tau.could_be_floer_simple()
        False

        >>> M = snappy.Manifold('m016')
        >>> tau = TuraevTorsion(M)
        >>> tau.could_be_floer_simple()
        True
        """
        return set(self.tau.coefficients()).issubset(set([0,1]))

    def support(self):
        return self.tau.support()
    
    def complement_of_support(self):
        """
        >>> M = snappy.Manifold('m007')
        >>> tau = TuraevTorsion(M)
        >>> sorted(tau.complement_of_support())
        [u*t, u^2, u^2*t]
        """
        tau, t, H = self.tau, self.QH(self.t), self.H
        n = t_deg_max(tau)
        f = sum_of_torsion_elements(H)
        X = f*sum(t**i for i in range(n+1))
        return set(X.support()) - set(tau.support())

    def lower_bound_on_thurston(self):
        """
        Proposition 2.2 of [RR]
        """
        return t_deg_max(self.tau) - 1

    def iota_image(self):
        """
        Part of the image of H_1(bdry; Z) in H_1(M; Z), namely a chunk
        large enough to contain all the interesting terms in self, up
        to the action of the integral longitude.

        >>> M = snappy.Manifold('m043')
        >>> M.homology()
        Z/5 + Z
        >>> tau = TuraevTorsion(M)
        >>> iota = tau.iota_image()
        >>> len(iota) == 5*(tau.lower_bound_on_thurston() + 1)
        True
        """
        m_ab, l_ab = self.meridian_ab, self.longitude_ab
        n = t_deg_max(self.tau)
        ans = dict()
        for i in range(n):
            for j in range(l_ab.order()):
                ans[m_ab**i * l_ab**j] = (i, j)
        return ans

    def D_tau_plus_pre(self):
        """
        The elements of H_1(Y) of the form x - y where x is not in the
        support of tau, y, is in the support of tau, and phi(x) >
        phi(y).

        >>> M = snappy.Manifold('K3_1')
        >>> sorted(TuraevTorsion(M).D_tau_plus_pre())
        [t, t^2, t^3, t^4, t^6, t^9]
        """
        S, S_comp = self.support(), self.complement_of_support()
        ans = set()
        for x in S_comp:
            for y in S:
                if t_deg(x) > t_deg(y):
                    ans.add(x*(y**-1))
        return ans

    def D_tau_plus(self):
        """
        The key set from Definition 1.5 of [RR].  Returns both the
        elements of the set and a linear combination of the
        *homological* peripheral basis that maps to it.

        >>> M = snappy.Manifold('m179')
        >>> tau = TuraevTorsion(M)
        >>> sorted(tau.D_tau_plus().items())
        [(t^2, (1, 0)), (t^4, (2, 0)), (t^8, (4, 0)), (u*t^2, (1, 1)), (u*t^4, (2, 1)), (u*t^6, (3, 1)), (u*t^12, (6, 1))]
        """
        D = self.D_tau_plus_pre()
        iota = self.iota_image()
        return {h:v for h, v in iota.items() if h in D}

        
def time_test():
    census = snappy.OrientableCuspedCensus(betti=1, cusps=1)
    manifolds = [census.random() for i in range(100)]
    for M in manifolds:
        TuraevTorsion(M)

        

if __name__ == '__main__':
    import doctest
    print(doctest.testmod())
