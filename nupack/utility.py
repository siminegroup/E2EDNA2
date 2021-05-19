'''
General purpose utility functions
'''
import re, json, uuid, datetime, pickle, copy
import numpy as np, pandas
from random import choice
from multiprocessing import cpu_count, pool
from typing import NamedTuple
import collections, weakref

################################################################################

class Rules:
    def __init__(self, rules, strict=True):
        self.rules = rules.rules() if hasattr(rules, 'rules') else dict(rules)
        self.strict = bool(strict)

    def get(self, key, default=None):
        '''Direct lookup as a dict'''
        return self.rules.get(key, default)

    def map(self, collection):
        '''Map a collection like list, tuple, set'''
        return type(collection)(map(self, collection))

    def __call__(self, item):
        return item.substitute(self)

################################################################################

def substitute(self, rules):
    '''Method to return a new version of self but with rules applied to each contained variable'''
    out = copy.copy(self)
    out.replace(rules if isinstance(rules, Rules) else Rules(rules))
    return out

def mappable(cls):
    '''Mixin decorator for substitute()'''
    assert cls.replace is not None
    cls.substitute = substitute
    return cls

################################################################################

def complement_name(string):
    if not isinstance(string, str):
        raise TypeError('Expected instance of str')
    if string and string[-1] == '*':
        return string[:-1]
    return string + '*'

################################################################################

def create_name(name):
    if isinstance(name, str):
        return name
    if name is None:
        raise TypeError('name must not be None')
    try:
        first, *last = name
    except TypeError:
        raise TypeError("If not a str, name should be a list like ['my-name', 0]")
    return str(name[0]) + ''.join(map('[{}]'.format, last))

################################################################################

class BaseResult(collections.abc.Mapping):
    def as_dict(self):
        return dict(zip(self.names, self.fields))

    def __getitem__(self, key):
        if isinstance(key, str):
            for field in self.fields:
                try:
                    return next(v for k, v in field.items() if getattr(k, 'name', '') == key)
                except StopIteration:
                    pass
        else:
            for field in self.fields:
                result = field.get(key)
                if result is not None:
                    return result

        raise KeyError(key)

    def __iter__(self):
        return (k for f in self.fields for k in f)

    def __len__(self):
        return sum(len(f) for f in self.fields)

    def __repr__(self):
        return '%s(%s)' % (type(self).__name__, ', '.join('%s=%r' % (k, v)
            for k, v in zip(self.names, self.fields) if v))

    def __str__(self):
        return '{%s}' % ', '.join('%s: %s' % (k,
            '{%s}' % ', '.join('%s: %s' % i for i in v.items()))
            for k, v in zip(self.names, self.fields) if v)

##########################################################################

def from_argument(cls, options):
    if options is None:
        return cls()
    if isinstance(options, cls):
        return options
    if isinstance(options, collections.abc.Mapping):
        return cls(**dict(options))
    raise TypeError('Expected None, mapping, or instance of %s; got %r' % (cls.__name__, options))

##########################################################################

class Loadable:
    def save(self, file):
        '''Save this object to a given binary file path'''
        with open(str(file), 'wb') as f:
            pickle.dump(self, f)

    @classmethod
    def load(cls, file):
        '''Load this type from a given binary'''
        with open(str(file), 'rb') as f:
            out = pickle.load(f)
        if cls is not Loadable and not isinstance(out, cls):
            raise TypeError('Loaded object is of type %r (expected %r)' % (type(out), cls))
        return out

##########################################################################

def long_output():
    return pandas.option_context('display.max_rows', 1000, 'display.max_colwidth', None)

class Printable:
    def save_text(self, file):
        with open(str(file), 'w') as f:
            with long_output():
                print(self, file=f)

def printable(cls):
    cls.save_text = Printable.save_text
    return cls

##########################################################################

def check_instance(item, expected_type):
    if not isinstance(item, expected_type):
        raise TypeError('Expected instance of %s (got %r)' % (expected_type.__name__, item))
    return item

def check_instances(items, expected_type):
    try:
        items = tuple(items)
    except TypeError:
        raise TypeError('Expected iterable (got %r)' % items)
    for i in items:
        if not isinstance(i, expected_type):
            raise TypeError('Expected only instances of %s (got %r)' % (expected_type.__name__, i))
    return items

##########################################################################

def str_repr(obj):
    return '{}.{}({})'.format(obj.__module__, type(obj).__name__, str(obj))

##########################################################################

def dtype_index(dtype):
    from . import rendered_document
    return rendered_document['scalars'][getattr(dtype, '__name__', str(dtype))]

##########################################################################

def match(gen):
    for x in gen:
        return x
    raise ValueError('Match not found')

##########################################################################

def nbits(x):
    try:
        x = np.dtype(x).itemsize * 8
    except TypeError:
        pass
    return int(x)

##########################################################################

def is_iterable(obj):
    try: return (True, iter(obj))[0]
    except TypeError: return False

##########################################################################

def thread_map(Iterable, function, threads=None, chunksize=None, index=False):
    if not is_iterable(Iterable):
        Iterable = range(Iterable)
    if threads is None:
        threads = cpu_count()
    if isinstance(threads, int):
        threads = pool.ThreadPool(threads)
    stuff = list(Iterable)
    if index:
        def f(x, i=[0]):
            function(i[0], stuff.pop(0))
            i[0] += 1
        return threads.map(f, tuple(range(len(stuff))), chunksize)
    else:
        def f(x): function(stuff.pop(0))
        return threads.map(f, tuple(range(len(stuff))), chunksize)

##########################################################################

def time_string(time=None, delta=None):
    if time is None: time = datetime.datetime.now()
    if delta is not None: time += datetime.timedelta(seconds=delta)
    return time.strftime('%H:%M')

##########################################################################

def can_convert(x, type):
    try: return (True, type(x))[0]
    except Exception: return False

##########################################################################

def specified_args(*args):
    return (i for i in args if i is not None)

##########################################################################

# see https://stackoverflow.com/questions/12397279/custom-json-encoder-in-python-with-precomputed-literal-json
class JSONEncoder(json.JSONEncoder):
    """
    JSON encoder which will work for objects containing NUPACK C++ types
    Works by using the default encoder, then carrying out replacements in the result
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._replacements = {}

    def default(self, o):
        if isinstance(o, Value):
            key = uuid.uuid4().hex
            self._replacements[key] = o.dump(-1 if self.indent is None else self.indent, True)
            return key
        else:
            return super().default(o)

    def encode(self, o):
        out = super().encode(o)
        for k, v in self._replacements.items():
             out = out.replace('"{}"'.format(k), v)
        return out

##########################################################################

def pair_matshow(ax, matrix, strands, period=4, line_options={'color': 'white'}, **kws):
    ax.matshow(matrix, **kws)
    s = list(''.join(strands))
    ax.set_xticks(range(len(s)))
    ax.set_yticks(range(len(s)))
    ax.set_xticklabels(['%d\n%s' % (i, s) if (i % period == 0) else s for i, s in enumerate(s)])
    ax.set_yticklabels(['%d %s' % (i, s) if (i % period == 0) else s for i, s in enumerate(s)])
    for i in np.cumsum(tuple(map(len, strands)))[:-1]:
        ax.plot([i+0.5,i+0.5], [-0.5, len(s)-0.5], **line_options)
        ax.plot([-0.5, len(s)-0.5], [i+0.5,i+0.5], **line_options)

##########################################################################

def DUp_recursion(toks):
    def match_left_paren(toks):
        toks = toks[2:]
        stack = 1
        for i, z in enumerate(toks):
            if isinstance(z, str):
                if z == "(":
                    stack += 1
                elif z == ")":
                    stack -= 1

            if stack == 0:
                return i + 2

    if not toks:
        return "", []
    tok = toks[0]

    if len(tok) == 1:
        assert tok == "+", "misplaced {}".format(tok)
        return tok, toks[1:]

    typ, n = tok

    if typ == "U":
        return "." * n, toks[1:]

    elif typ == "D":
        assert len(toks) > 1, "Duplex must have either strand break or internal structure"
        nxt = toks[1]
        if len(nxt) == 2 and nxt[0] == "D":
            mid, rest = DUp_recursion(toks[1:])
            intermediate = "("*n + mid + ")"*n
            return intermediate, rest

        if nxt == "(":
            matching = match_left_paren(toks)
            middle = toks[2:matching]
        else:
            matching = 1
            middle = [nxt]


        intermediate = ""
        while middle:
            string, middle = DUp_recursion(middle)
            intermediate += string
        return "(" * n + intermediate + ")" * n, toks[matching+1:]

##########################################################################

def DUp_to_dpp(string):
    def tokenize(string):
        token_types = r"\(|\)|D[1-9][0-9]*|U[1-9][0-9]*|[+]"
        components = r"D|U|[0-9]+"
        temp = re.findall(token_types, string)
        ret = list()
        for x in temp:
            j = re.findall(components, x)
            if j:
                ret.append((j[0], int(j[1])))
            else:
                ret.append(x)
        return ret

    toks = tokenize(string)
    ret = ""
    while toks:
        string, toks = DUp_recursion(toks)
        ret += string
    return ret

