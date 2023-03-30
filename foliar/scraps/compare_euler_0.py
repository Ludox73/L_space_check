import edge_orient, link
import snappy.snap.t3mlite as t3m

def compare_new_with_old(datum):
    M = t3m.Mcomplex(datum['laminar_tri'])
    L = link.LinkSphere(M)
    old_euler_0 = eval(datum['taut_euler_0'])
    new_euler_0 = []
    for signs in eval(datum['laminar_orients']):
        eo = edge_orient.EdgeOrientation(M, L, signs)
        new_euler_0.append(1 if eo.euler_class_vanishes() else 0)
    return sum(old_euler_0), sum(new_euler_0)
