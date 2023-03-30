This package was downloaded by https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/LCYXPO and written by Nathan Dunfield.

It is contained here just for convenience. We deleted the database of foliations present in Dunfield's version, because we do not need it.


===================
Foliar orientations
===================

To install the "foliar" module into SageMath+SnapPy do the following
in this directory::

  sage -pip install .

You can test to make sure all is well by doing::

  sage -python -m foliar.test

which should conclude with a message such as::

  All doctests:
    0 failures out of 73 tests.

Example of searching for a foliar orientation on various
triangulations of a closed manifold::

  sage: import snappy, foliar
  sage: M = snappy.Manifold('m004(1,2)')
  sage: eo = foliar.first_foliation(M, 5, 25)
  sage: eo.gives_foliation()
  True
  sage: eo.euler_class_vanishes()
  True

Example of searching for a persistently foliar orientation on a
1-cusped manifold::

  sage: m004 = snappy.Manifold('m004')
  sage: foliar.degeneracy_slopes(m004)
  [(1, 0)]
  sage: m003 = snappy.Manifold('m003')
  sage: foliar.degeneracy_slopes_with_search(m003)
  ([], [])

Note here that m003 is Floer simple and hence nothing is found.

For many more useful examples, see the docstrings the files
"src/main.py" and "src/edge_orient.py".


Data
====

The file "foliar_tris/persistent_knots.csv.bz2" contains the
persistent foliar orientations for some 1.2 million knot exteriors.


Todo
====

It would not be hard to merge the code here into SnapPy proper, and so
this should be done.
