import json
from . import double_group

class MonoidInGroup(object):
    """
    If track is True, remember how each element in self can be
    expressed in terms of the given monoid generators.
    """
    def __init__(self, elements, ball, saturate=True, track=False):
        self.ball, self.track = ball, track
        if track:
            self.expressed_in_gens = dict()
            self._word_rep_one = None
        if saturate:
            self.elements = set()
            self._has_one = self.saturate(elements)
        else:
            self.elements = elements.copy()
        
    def saturate(self, new_elements):
        """
        Returns whether 1 is in self after saturation
        """
        active = set(new_elements)
        ans = self.elements | set(active)

        track = self.track
        if track:
            in_gens = self.expressed_in_gens
            for x in active:
                in_gens[x] = [x.word]

        while len(active) > 0:
            new_elts = set()
            for x in ans:
                for y in active:
                    for a, b in [(x, y), (y, x)]:
                        z = a*b
                        if z in self.ball.elements and z not in ans:
                            # The next application of the "identity" map is
                            # to try to prevent accumulation of numerical error.
                            z = self.ball.element_dict[z]
                            new_elts.add(z)
                            if track:
                                in_gens[z] = in_gens[a] + in_gens[b]                            
                            if z.is_one():
                                self.elements = ans | new_elts
                                self._has_one = True
                                if track:
                                    self._word_rep_one = in_gens[z]
                                return True

            ans = ans | new_elts
            active = new_elts

        self.elements = ans
        return False

    def has_one(self):
        return self._has_one

    def copy(self):
        M = MonoidInGroup(self.elements, self.ball, False, self.track)
        M._has_one = self._has_one
        if self.track:
            M._word_rep_one = self._word_rep_one
            M.expressed_in_gens = self.expressed_in_gens.copy()
        return M
        
    def words(self):
        ans = [a.word for a in self.elements]
        ans.sort(key=lambda x: (len(x), x))
        return ans

    def __contains__(self, x):
        return x in self.elements

    def __len__(self):
        return len(self.elements)

class Printer(object):
    def __init__(self, silent=False):
        self.silent = silent

    def size_of_ball(self, B, depth):
        self.write('Ball has %d elements' % len(B.elements), depth)
        
    def add_monoid_gen(self, element, depth):
        self.write('Adding %s' % element.word, depth)

    def size_of_monoid(self, P, depth):
        self.write('Size of P is %d' % len(P.elements), depth)

    def contradiction(self, word, depth):
        if word:
            word = '*'.join(word)
            self.write('Contradiction: ' + word + ' = 1 in P', depth)
        else:
            self.write('Contradiction: 1 in P', depth)
            
    def write(self, string, depth):
        if not self.silent:
            print('  '*depth + '%d: ' % depth + string)

    def proof_string(self, manifold, args):
        G = manifold.fundamental_group(*args)
        ans = {'name':repr(manifold),
               'group_args':[1 if x else 0 for x in args],
               'gens':'.'.join(G.generators()),
               'rels':G.relators(),
               'proof':self.value,
               }
        
        return json.dumps(ans).replace(' ', '')

    
    
class ProofPrinter(Printer):
    """
    Stores a proof that a group G is nonorderable.  The proof itself
    is a rooted binary tree, oriented away from the root, with the
    following structure:

    1. Each edge is labelled by a word that represents a *nontrivial*
    element of G; the labels on the two edges leaving a vertex are
    inverses *words* of one another.

    2. Each leaf is labelled by a word in the edge labels that appear
    along the unique path from the leaf back to the root.

    3. The word associated to each leaf is trivial in G.

    We output the proof as a sequence of pairs

         [(edge path from root to leaf), (word at leaf)]

    Both items in the pair are strings consisting of words in [a-zA-Z]
    separated by periods.

    To check a such a proof, one must verify:

    a) All edge labels give non-trivial elements of G.

    b) The collection of edge paths forms a binary tree satisfying (2)
    
    c) Each leaf word gives 1 in G and is in fact a product of the
    preceding edge labels.
    """
    def __init__(self, silent):
        self.silent = silent
        self.edges_back_to_root = []
        self.value = []

    def add_monoid_gen(self, element, depth):
        word = element.word
        self.edges_back_to_root = self.edges_back_to_root[:depth] + [word]
        Printer.add_monoid_gen(self, element, depth)

    def contradiction(self, word, depth):
        if word is not None:
            self.value.append(['.'.join(self.edges_back_to_root), '.'.join(word)])
        self.edges_back_to_root = self.edges_back_to_root[:-1]
        Printer.contradiction(self, word, depth)
    
def ball_has_order(B, P, printer, recur_depth):
    printer.size_of_monoid(P, recur_depth)
    if P.has_one():
        if P.track:
            word = P._word_rep_one
        else:
            word = None
        printer.contradiction(word, recur_depth)
        return False, P

    # If we get close to building a full P, we almost always get
    # there.  It is better to stop now and later try again with an
    # increased radius, because adding those last few elements is very
    # expensive.
    if len(P) > 0.9 * 0.5 * len(B.elements):
        return True, P

    for x, y in B.non_id_element_pairs:
        if (x not in P) and (y not in P):
            printer.add_monoid_gen(x, recur_depth)
            newP = P.copy()
            newP.saturate([x])
            ans = ball_has_order(B, newP, printer, recur_depth+1)
            if ans[0]:
                return ans
            else:
                printer.add_monoid_gen(y, recur_depth)
                newP = P.copy()
                newP.saturate([y])
                return ball_has_order(B, newP, printer, recur_depth+1)

    return True, P

def has_non_orderable_group(manifold, ball_radius=3,
                            silent=False, track=False, return_proof=False,
                            min_bits_accuracy=15,
                            fundamental_group_args = [True, True, False]):
    """
    >>> import snappy
    >>> M = snappy.Manifold('m003(-3,1)')
    >>> has_non_orderable_group(M, track=True)  # doctest: +SKIP
    0: Ball has 151 elements
    0: Adding a
      1: Size of P is 3
      1: Adding b
        2: Size of P is 104
        2: Contradiction: b*a*b*a*b*b*a*a*a*b = 1 in P
      1: Adding B
        2: Size of P is 14
        2: Adding c
          3: Size of P is 97
          3: Contradiction: a*B*c*a*a*c*c = 1 in P
        2: Adding C
          3: Size of P is 54
          3: Contradiction: C*a*C*a*a*B = 1 in P
    True
    >>> ans, proof = has_non_orderable_group(M, silent=True, return_proof=True)
    >>> ans
    True
    >>> json.loads(proof)['proof']   # doctest: +SKIP
    [[u'a.b', u'b.a.b.a.b.b.a.a.a.b'], [u'a.B.c', u'a.B.c.a.a.c.c'], [u'a.B.C', u'C.a.C.a.a.B']]
    >>> N = snappy.Manifold('m004(1, 2)')
    >>> has_non_orderable_group(N)
    0: Ball has 159 elements
    0: Adding a
      1: Size of P is 3
      1: Adding b
        2: Size of P is 14
        2: Adding c
          3: Size of P is 36
          3: Adding aB
            4: Size of P is 49
            4: Adding aC
              5: Size of P is 67
              5: Adding bC
                6: Size of P is 70
                6: Adding aBB
                  7: Size of P is 72
    False


    >>> names = ['m003(-1, 3)', 'm003(-1, 5)', 'm003(-3, 1)', 'm003(-4, 3)', 'm003(-5, 2)', 'm003(-5, 4)', 'm004(-3, 2)', 'm004(-6, 1)']
    >>> manifolds = [snappy.Manifold(name) for name in names]
    >>> [has_non_orderable_group(M, silent=True) for M in manifolds]
    [False, False, True, True, True, True, False, False]
    """
    if return_proof:
        track = True
        printer = ProofPrinter(silent)
    else:
        printer = Printer(silent)
    G = double_group.Double3ManifoldGroup(
              manifold, min_bits_accuracy, fundamental_group_args)
    a = G('a')
    B = G.ball(ball_radius)
    printer.size_of_ball(B, 0)
    printer.add_monoid_gen(a, 0)
    P = MonoidInGroup([a], B, track=track)
    ans = not ball_has_order(B, P, printer, 1)[0]
    if return_proof:
        if ans:
            return ans, printer.proof_string(manifold, fundamental_group_args)
        else:
            return ans, None
    return ans

