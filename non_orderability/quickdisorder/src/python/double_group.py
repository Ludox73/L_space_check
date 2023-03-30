import snappy
from . import sl2matrix, ball

class Double3ManifoldGroup(object):
    """
    >>> M = snappy.Manifold('m004(1,2)')
    >>> G = Double3ManifoldGroup(M)
    >>> G('abcABC')
    <NGE: abcABC; (-0.25857712442742253+0.48081251833270267j) (-0.783996070719998+2.150934815529558j) (-0.2774615146731634+0.23771079574212495j) (-1.8761415167194206-0.4598487432176286j)>
    >>> B = G.ball(2)
    >>> len(B)
    37
    """
    def __init__(self, manifold, min_bits_accuracy=15, fundamental_group_args = [True, True, False]):
        self.manifold = manifold
        self.rho = snappy.snap.polished_holonomy(manifold, 100,
                                                 lift_to_SL2=True,
                                                 fundamental_group_args=fundamental_group_args)
        self.min_bits_accuracy = min_bits_accuracy

    def __call__(self, word):
        A = self.rho(word)
        return sl2matrix.DoubleGroupElement(A, self.min_bits_accuracy, word)

    def ball(self, radius):
        return ball.CayleyBall(self, radius)



