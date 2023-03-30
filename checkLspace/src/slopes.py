""" This code was downloaded by https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/LCYXPO and written by Nathan Dunfield. """

from sage.all import (gcd, xgcd, vector, matrix, Polyhedron, 
                      Graphics, point2d, RR, cos, sin)

class Slope(object):
    """
    An unoriented multicurve on a torus T, with respect to an implicit
    positively oriented basis for H_1(T; Z).

    >>> w, x, y, z = Slope(3, -2), Slope(-3, 2), Slope(-1, 0), Slope(0, -2)
    >>> w, x, y, z
    (Slope(-3, 2), Slope(-3, 2), Slope(1, 0), Slope(0, 2))
    """
    def __init__(self, a, b=None):
        if b is None:
            a, b = a
        if a*b == 0:
            a, b = abs(a), abs(b)
        elif b < 0:
            a, b = -a, -b
        self.tuple = (a,b)

    def __hash__(self):
        return hash(self.tuple)

    def __add__(self, other):
        """
        Orient each slope and take the algebraic sum; there are four ways to
        do this with two possible results.
        >>> Slope(-1, -2) + Slope(3, -2)
        (Slope(-2, 4), Slope(4, 0))
        """
        a, b = self.tuple
        c, d = other.tuple
        return Slope(a+c, b+d), Slope(a-c, b-d)

    def num_components(self):
        """
        >>> Slope(6, 4).num_components()
        2
        """
        return abs(gcd(*self.tuple))

    def primitive_slope(self):
        """
        >>> Slope(6, 4).primitive_slope()
        Slope(3, 2)
        """
        a, b = self.tuple
        g = gcd(a, b)
        return Slope(a/g, b/g)

    def complement(self):
        """
        >>> Slope(2, 3).complement()
        Slope(1, 1)
        """
        g, s, t = xgcd(*self.tuple)
        ans = Slope(-t, s)
        assert ans*self == 1
        return ans

    def __getitem__(self, i):
        return self.tuple[i]
    
    def __mul__(self, other):
        """
        Geometric intersection number if another slope,
        otherwise tries multiplication as a vector.

        >>> x, y = Slope(3, -2), Slope(-1, 1)
        >>> x*y
        1
        """
        if isinstance(other, Slope):
            c, d = other.tuple
            a, b = self.tuple
            return abs(a*d - b*c)

    def __rmul__(self, other):
        """
        >>> x = Slope(1, 0)
        >>> A = matrix([[2, 1], [-3, -1]])
        >>> A*x
        Slope(-2, 3)
        >>> -5*x
        Slope(5, 0)
        """
        return Slope(other*vector(self.tuple))

    def __eq__(self, other):
        """
        >>> Slope(3, -5) == Slope(-3, 5)
        True
        >>> Slope(2, 0) == [2, 0]
        True
        >>> Slope(3, 5) != Slope(-3, -5)
        False
        """
        return self.tuple == Slope(other).tuple

    def __ne__(self, other):
        return self.tuple != Slope(other).tuple

    def __lt__(self, other):
        self.tuple < Slope(other).tuple

    def __le__(self, other):
        self.tuple <= Slope(other).tuple

    def __gt__(self, other):
        self.tuple > Slope(other).tuple

    def __ge__(self, other):
        self.tuple >= Slope(other).tuple
        
    def __repr__(self):
        return 'Slope(%d, %d)' % self.tuple


class SlopeCone(object):
    """
    An open connected arc in P^1(Q).  Note this circle inherits a
    counter-clockwise orientation from the orientation of Q^2.
    """
    def __init__(self, u=None, v=None, contains=None, avoids=None):
        
        """
        Specify the SlopeCone by the two endpoints in anticlockwise
        orientation or via the "avoids" and "contains" arguments.

        >>> C0 = SlopeCone( (1, 0), (0, 1) )
        >>> (1, 1) in C0
        True
        >>> (1, 0) in C0
        False
        >>> C1 = SlopeCone((-1, 0), (0, 1), avoids=(1, -1))
        >>> C2 = SlopeCone((-1, 0), (0, 1), contains=(1, 1))
        >>> C1 == C2
        True
        >>> (1, 1) in SlopeCone((2, 1), (1, 2), contains=(1,1))
        True
        >>> (1, 1) in SlopeCone((1, 2), (2, 1), avoids=(1,1))
        False

        The case where u and v are equal, or where v is None corresponds
        to the complementary interval to u.

        >>> C3 = SlopeCone((1, 0))
        >>> (1,0) in C3
        False
        >>> (0, 1) in C3
        True

        Can also create from a collection of vectors in Z^2 which span
        a proper cone, or equivalently lie in some half-space through
        the origin.

        >>> SlopeCone(contains=[(0,1), (1, 1), (-1, 1), (2, 1), (-5, 3)])
        SlopeCone((2, 1), (-5, 3))
        >>> SlopeCone(contains=[(1,0), (1, 1), (-2, -1)])
        Traceback (most recent call last):
        ...
        ValueError: Postive cone of vectors is all of R^2
        >>> SlopeCone(contains=[(-1,0), (1, 0), (0, 1)])
        SlopeCone((1, 0))
        """
        if u is None: 
            assert v is None and avoids is None
            P = Polyhedron(vertices=contains + [(0,0)])
            if P.interior_contains((0,0)):
                raise ValueError('Postive cone of vectors is all of R^2')
            if P.dimension() < 2:
                raise ValueError('Postive cone is contained in a line')
            vert_at_zero = [vert for vert in P.vertices() if vert.vector() == 0]
            if vert_at_zero:
                z = vert_at_zero[0]
                u, v = [n.vector() for n in z.neighbors()]
                contains = u + v
            else: # half-plane
                for face in P.faces(face_dimension=1):
                    if face.as_polyhedron().contains((0,0)):
                        u = v = tuple(face.vertices()[0])
        elif v is None:
            v = u
        self.u, self.v = Slope(u), Slope(v)
        if self.u != self.v:
            if contains is not None:
                assert avoids is None
                if contains not in self:
                    self.u, self.v = self.v, self.u
            elif avoids is not None:
                assert contains is None
                if avoids in self:
                    self.u, self.v = self.v, self.u
            
    def  __contains__(self, x):
        """
        Test whether the *ordered* tripple of points (u, x, v) has the same 
        the orientation on P^1(Q) as the triple ( (-1, 1), (1, 0), (1, 1) ), which
        is to say "anticlockwise".
        """
        u0, u1 = self.u
        x0, x1 = x
        v0, v1 = self.v
        if  self.u != self.v:
            return (u0*x1 - u1*x0)*(x0*v1-x1*v0)*(v0*u1-v1*u0) < 0
        else:
            return u0*x1 - u1*x0 != 0

    def __repr__(self):
        if self.u == self.v:
            return 'SlopeCone((%d, %d))' % tuple(self.u)
        return 'SlopeCone((%d, %d), (%d, %d))' % (tuple(self.u) + tuple(self.v))

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self, other):
        if isinstance(other, SlopeCone):
            return self.u == other.u and self.v == other.v
        
    def draw(self, N=100):
        G = Graphics()
        s = 2*RR.pi()/N
        points = []
        for k in range(N):
            x = (cos(s*k), sin(s*k))
            if x in self:
                points.append(x)
        ans = point2d(points, color='black')
        ans.set_aspect_ratio(1)
        u, v = vector(self.u), vector(self.v)
        u, v = u/u.norm(), v/v.norm()
        ans += point2d([u, -u], color='red')
        ans += point2d([v, -v], color='blue')
        return ans
        
class AllSlopes(object):
    """
    All of P^1(Q).  

    >>> A = AllSlopes()
    >>> (1, 0) in A
    True
    >>> B = AllSlopes()
    >>> A == B
    True
    >>> C = SlopeCone((1, 0), (0, 1))
    >>> A == C
    False
    """
    def __contains__(self, x):
        return True

    def __repr__(self):
        return 'AllSlopes()'

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self, other):
        if isinstance(other, AllSlopes):
            return True
        if isinstance(other, SlopeCone):
            return False


class SingleSlope(object):
    """
    A single point in P^1(Q).  

    >>> S = SingleSlope(-8, -4)
    >>> S
    SingleSlope(2, 1)
    >>> (-8, -4) in S
    True
    >>> (3, -1) in S
    False
    >>> C = SlopeCone((1, 0), (0, 1))
    >>> S == C
    False
    >>> S == AllSlopes()
    False
    >>> S == SingleSlope(8, 4)
    True
    >>> S == SingleSlope(1, 0)
    False
    """
    def __init__(self, a, b=None):
        self.slope = Slope(a, b).primitive_slope()
        
    def __contains__(self, x):
        return self.slope == Slope(x).primitive_slope()
    
    def __repr__(self):
        return 'SingleSlope(%d, %d)' % self.slope.tuple

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self, other):
        if isinstance(other, AllSlopes):
            return False
        if isinstance(other, SlopeCone):
            return False
        if isinstance(other, SingleSlope):
            return self.slope == other.slope
                
                
    

if __name__ == '__main__':
    import doctest
    doctest.testmod()
    
