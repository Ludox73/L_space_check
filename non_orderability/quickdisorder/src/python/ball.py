
def inverse_word(word):
    return word.swapcase()[::-1]

def next_gen_dict(gens):
    all_gens = gens + [g.swapcase() for g in gens]
    next_gen = { g:[h for h in all_gens if h.swapcase() != g] for g in all_gens }
    return next_gen

def ball_in_free_group(gens,  length):
    next_gen = next_gen_dict(gens)
    curr = gens + [g.swapcase() for g in gens]
    ans = [''] + curr
    for i in range(length - 1):
        new_words = []
        for word in curr:
            new_words += [word + g for g in next_gen[word[-1]]]
        ans += new_words
        curr = new_words
    return ans

class CayleyBall(object):
    """
    Ball about 1 in the Cayley graph of a group.
    """
    def __init__(self, group, radius):
        all_words = ball_in_free_group(group.rho.generators(), radius)
        ordered_elements = [group(w) for w in all_words]
        elements = set(ordered_elements)
        # Create the list of {g, g^-1} pairs
        seen_one = dict()
        non_id_element_pairs = []
        for g in ordered_elements:
            h = g.inverse()
            if h in seen_one:
                h = seen_one[h]
                # Want the words associated to g and h be inverses in
                # the free group.
                if h.word != inverse_word(g.word):
                    g._set_word(inverse_word(h.word))
                non_id_element_pairs.append( (seen_one[h], g) )
            else:
                seen_one[g] = g

        self.ordered_elements = ordered_elements
        self.elements = elements
        self.non_id_element_pairs = non_id_element_pairs
        self.element_dict = {e:e for e in ordered_elements}

    def __len__(self):
        return len(self.elements)
