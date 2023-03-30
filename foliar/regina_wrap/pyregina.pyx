from libcpp.string cimport string
from libcpp.functional cimport function
import os, sys, re, tempfile

cdef extern from "engine.h" namespace "regina":
    cdef int versionMajor()
    cdef int versionMinor()

cdef extern from "maths/integer.h" namespace "regina":
    cdef cppclass LargeInteger:
        long longValue()

cdef extern from "triangulation/dim3.h" namespace "regina":
    cdef cppclass NTriangulation:
        NTriangulation()
        NTriangulation(string)
        unsigned long countTetrahedra()
        unsigned long countVertices()
        bint intelligentSimplify()
        bint simplifyToLocalMinimum(bint perform)
        bint retriangulate(int, unsigned int, void*, function[bint(NTriangulation&)])
        string isoSig()

_global_action_isosigs = set()
_global_action_config = {'print':True}

cdef bint action(NTriangulation& triangulation):
    if triangulation.simplifyToLocalMinimum(False):
        isosig = triangulation.isoSig()
        n = int(triangulation.countTetrahedra())
        if _global_action_config['print']:
            print(isosig)
        _global_action_isosigs.add((n, isosig))
        
def version():
    return "%d.%d" % (versionMajor(), versionMinor())

cdef class Triangulation:
    """
    A Regina triangulation of a 3-manifold, which can be specified
    either by name of a triangulation file in SnapPea format or by giving
    a SnapPy triangulation.
    """
    cdef NTriangulation* triangulation

    def __cinit__(self, spec):
        if not isinstance(spec, str):
            spec = spec._to_string()
        self.triangulation = new NTriangulation(spec)

    def __dealloc__(self):
        del self.triangulation
            
    def num_tetrahedra(self):
        return int(self.triangulation.countTetrahedra())

    def num_vertices(self):
        return int(self.triangulation.countVertices())

    def simplify(self):
        cdef bint did_simplify = self.ntriangulation.intelligentSimplify()
        return did_simplify

    def isosig(self):
        return self.triangulation.isoSig()

    def retriangulate(self, height, print_each_isosig=False):
        _global_action_isosigs.clear()
        _global_action_config['print'] = print_each_isosig
        cdef function[bint(NTriangulation&)] *callback
        callback = new function[bint(NTriangulation&)](action)
        self.triangulation.retriangulate(height, 1, NULL, callback[0])
        return _global_action_isosigs
