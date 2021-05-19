import itertools

################################################################################

def partition(collection):
    """Returns a generator of tuples"""
    collection = tuple(collection)

    if len(collection) == 1:
        yield (collection,)
        return

    first = (collection[0],)
    for smaller in partition(collection[1:]):
        # insert `first` in each of the subpartition's subsets
        for n, subset in enumerate(smaller):
            yield smaller[:n] + (first + subset,)  + smaller[n+1:]
        # put `first` in its own subset
        yield (first,) + smaller

################################################################################

def all_permutations(stuff):
    '''All possible assignments of any permutations (for a list of lists)'''
    if len(stuff) == 1:
        yield from ((p,) for p in itertools.permutations(stuff[0]))
    else:
        for p in itertools.permutations(stuff[0]):
            yield from ((p,) + i for i in all_permutations(stuff[1:]))

################################################################################

def lowest_rotation(perm):
    '''Lowest rotation of a list'''
    i = perm.index(min(perm))
    return perm[i:] + perm[:i]

################################################################################

def unique_permutations(stuff):
    '''All rotationally unique assignments of any permutations (for a list of lists)'''
    ret = set()
    for perm in all_permutations(stuff):
        ret.add(tuple(map(lowest_rotation, perm)))
    return tuple(ret)

################################################################################

def unique_partitions(elements):
    '''All partitionings of rotationally unique assignments of any permutations (for a list of lists)'''
    return list(j for i in partition(tuple(elements)) for j in unique_permutations(i))

################################################################################

def sample_state(n, engine, seq, model):
    return [State(seq, s.dp([len(seq)]), model=model, kind='slow') for s in engine.sample(seq, n)]

################################################################################

def disconnected_pf(engine, scale, strands, perm):
    ret = 1.
    for i in perm:
        seq = '+'.join(strands[j] for j in i)
        ret *= engine(seq)
    return ret * scale ** (len(perm) - len(strands))

################################################################################

def partitioned_pfs(engine, strands):
    partitions = unique_partitions(range(len(strands)))
    scale = (mod.unimolecular_scaling / mod.bimolecular_scaling / mod.molarity)
    pfs = np.asarray([calculate_pf(engine, scale, strands,) for p in partitions])
    pfs /= pfs.sum()

################################################################################

def sample_complexes():
    indices = np.random.choice(len(partitions), p=bias_pfs, size=n_launch)
    launch_perms = [partitions[i] for i in launch_indices]
    weights = [pfs[i] / bias_pfs[i] for i in launch_indices]

################################################################################

def perm_to_state(perm):
    # where each strand starts in standard representation
    nicks = [1]
    for s in seqs:
        nicks.append(nicks[-1] + len(s) + 2) # where each strand starts, but 2 null bases between each strand

    # where each strand starts in this representation

    ret = np.arange(sum(len(s) for s in seqs) + 2 * len(seqs)) # standard pair vec

    for cx in perm:
        pairs = np.asarray(eng.sample('+'.join(seqs[i] for i in cx))[0]) # no null bases included

        transfer = [] # where bases are supposed to go
        for s in cx:
            for i in range(len(seqs[s])):
                transfer.append(nicks[s] + i)
        transfer = np.asarray(transfer, dtype=int)

        start = 0
        for s in cx:
            ret[nicks[s]:nicks[s]+len(seqs[s])] = transfer[pairs[start:start+len(seqs[s])]]
            start += len(seqs[s])

    pairs = State('+'.join(seqs)).pairs
    for i, j in enumerate(ret):
        pairs[i] = j

    return State('+'.join(seqs), model=mod, kind='slow').with_structure(pairs)

################################################################################

#print(perm_to_state(launch_perms[0]))

