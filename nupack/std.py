import operator, numpy
from .rebind import forward

################################################################################

@forward
class Tuple:
    '''Class modelling a (fixed-length) C++ std::tuple, std::pair, or std::array'''

    def __getitem__(self, i, _fun_=None):
        if isinstance(i, slice):
            return tuple(map(self.__getitem__, tuple(range(len(self)))[i]))
        return _fun_(self, i)._set_ward(self)

    def __setitem__(self, i, value):
        self[i].copy_from(value)

    def __len__(self) -> int:
        pass

    def __bool__(self):
        return bool(len(self))

    def __iter__(self):
        return map(self.__getitem__, range(len(self)))

################################################################################

@forward
class Vector(Tuple):
    '''Class modelling a (variable-length) C++ std::vector, std::deque, etc.'''

    def value_type(self) -> 'TypeIndex':
        pass

    def append(self, item) -> None:
        '''Append an item'''

    def __iadd__(self, other):
        for x in other:
            self.append(x)
        return self

    def view(self, copy=True):
        out = self.cast(numpy.ndarray)
        return out.copy() if copy else out

################################################################################

@forward
class Iterator:
    '''Class modelling an iterator'''

    def good(self) -> bool:
        pass

    def __next__(self, _fun_=None):
        _fun_(self)
        if self.good():
            return self.get()
        raise StopIteration

    def __iter__(self):
        return self

    def get(self, _fun_=None):
        out = _fun_(self)
        if out.qualifier() != 0:
            out._set_ward(self)
        return out

################################################################################

@forward
class Map:
    '''Class modelling a C++ std::map'''

    def __setitem__(self, key, value) -> None:
        '''Emplace an item'''

    def __getitem__(self, key, _fun_=None):
        return _fun_(self, key)._set_ward(self)

    def __len__(self) -> int:
        pass

################################################################################

@forward
class String:
    '''Class modelling a C++ std::string'''

    def __pos__(self):
        return self.cast(str)

    def __str__(self):
        return self.cast(str)

    def __repr__(self):
        return 'String({})'.format(repr(self.cast(str)))

    def __len__(self) -> int:
        pass

    def __iter__(self):
        return iter(self.cast(str))

    def __hash__(self):
        return self.cast(str).__hash__()

for k in 'add mul mod eq ne lt gt le ge'.split():
    op = getattr(operator, k)
    setattr(String, '__%s__' % k, lambda self, other, op=op: op(+self, other))
    setattr(String, '__r%s__' % k, lambda self, other, op=op: op(other, +self))
    setattr(String, '__i%s__' % k, lambda self, other, op=op: self.move_from(op(+self, other)))

################################################################################

@forward
class _Scalar:

    def __float__(self):
        return self.cast(float)

    def __int__(self):
        return self.cast(int)

    def __bool__(self):
        return self.cast(bool)

    def __getattr__(self, key):
        return getattr(+self, key)

    def __str__(self):
        return str(+self)

    def __neg__(self):
        return -(+self)

    def __abs__(self):
        return abs(+self)

    def __invert__(self):
        return ~(+self)

for k in 'add sub mul floordiv truediv mod pow lshift rshift and xor or pow eq ne lt gt le ge'.split():
    op = getattr(operator, {'and': 'and_', 'or': 'or_'}.get(k, k))
    setattr(_Scalar, '__%s__' % k, lambda self, other, op=op: op(+self, other))
    setattr(_Scalar, '__r%s__' % k, lambda self, other, op=op: op(other, +self))
    setattr(_Scalar, '__i%s__' % k, lambda self, other, op=op: self.move_from(op(+self, other)))

################################################################################

@forward
class Float(_Scalar):
    '''Class modelling a C++ floating point object'''

    def __pos__(self):
        return self.cast(float)

    def __float__(self):
        return self.cast(float)

    def __repr__(self):
        return 'Float({})'.format(self.cast(float))

################################################################################

@forward
class Integer(_Scalar):
    '''Class modelling a C++ integer object'''

    def __pos__(self):
        return self.cast(int)

    def __int__(self):
        return self.cast(int)

    def __repr__(self):
        return 'Integer({})'.format(self.cast(int))

################################################################################

@forward
class Bool(_Scalar):
    '''Class modelling a C++ boolean object'''

    def __pos__(self):
        return self.cast(bool)

    def __bool__(self):
        return self.cast(bool)

    def __repr__(self):
        return 'Bool({})'.format(self.cast(bool))

################################################################################

@forward
class Optional:
    def __bool__(self) -> bool:
        '''Return whether self contains a value'''

    def value(self):
        '''Return held object or throw exception if there is none held'''

    def get(self):
        '''Return held object or None'''
        return self.value() if self else None

################################################################################

@forward
class Variant:
    def get(self):
        '''Return held object as its specific type'''
