# currently unused or old functions

def extractPeptide(file):
    '''
    if ssDNA and peptide are concatenated, unlabelled in a PDB, extract the peptide, if it's second
    also extracts the DNA - given it's a lightdock complex file
    '''
    proteinResidues = ['ALA', 'ASN', 'CYS', 'GLU', 'HIS', 'LEU', 'MET', 'PRO', 'THR', 'TYR', 'ARG', 'ASP', 'GLN', 'GLY', 'ILE', 'LYS', 'PHE', 'SER', 'TRP', 'VAL']

    f = open(file, 'r')
    text = f.read()
    f.close()
    text = text.split('\n')

    proteins = []
    lineInd = 1
    for line in text:
        for string in proteinResidues:
            if string in line:
                proteins.append(lineInd)  # index from 1
                break
        lineInd += 1

    out = [text[i-1] for i in proteins]
    out = "\n".join(out)
    f = open('extractedPeptide.pdb','w')
    f.write(out)
    f.close()

    out2 = text[:proteins[0]-1]
    out2 = "\n".join(out2)
    f = open('extractedAptamer.pdb','w')
    f.write(out2)
    f.close()


def dnaBasePairDist(u,sequence): # nuclinfo is slow
    pairDists = np.zeros((len(u.trajectory),len(sequence),len(sequence))) # matrix of watson-crick distances
    tt = 0
    for ts in u.trajectory:
        for i in range(len(sequence)):
            for j in range(len(sequence)):
                if i > j:
                    pairDists[tt,i,j] = nuclinfo.wc_pair(u,i+1,j+1,seg1='A',seg2='A') # also do WC base-pair distances (can set to only follow secondary structure prediction)
                    pairDists[tt,j,i] = pairDists[tt,i,j]

        tt += 1

    return pairDists


def dnaBaseDihedrals(u,sequence): # nuclinfo is slow
    angles = np.zeros((len(u.trajectory),len(sequence)-3,7))
    tt = 0
    for ts in u.trajectory:
        for j in range(2,len(sequence)-1):
            angles[tt,j-2,:] = nuclinfo.tors(u,seg="A",i=j)

        tt += 1

    return angles


def recenterPDB(structure):
    pdb = app.PDBFile(structure)
    topology = pdb.topology
    positions = pdb.positions
    modeller = app.Modeller(topology, positions)
    app.PDBFile.writeFile(modeller.topology, modeller.positions, open(structure.split('.')[0] + '_recentered.pdb','w'))

