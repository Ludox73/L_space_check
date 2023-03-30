import doctest
from . import double_group, disorder

for module in [double_group, disorder]:
    print(module.__name__ + ': ' + repr(doctest.testmod(module)))
