from utils import *

# we want to collate the relevant trajectories and all the relevant analyses

os.chdir('../../../../Desktop/mmruns/cluster')


run_num = 998
dir = 'run{}'.format(run_num) # go into the directory

dirList = os.listdir(dir)
os.chdir(dir)

outputs = np.load('opendnaOutput.npy',allow_pickle=True).item() # load outputs

runPairs = []

# recenter trajectory files
for file in dirList: # recenter .dcd files for better visualization
    if 'clean' in file:
        if file[-4:] == '.dcd':
            trajectory = file
            for otherFile in dirList: # find pairing .pdb file
                fileNoExt = otherFile.split('.')[0]
                if fileNoExt.split('_')[:4] == file.split('_')[:4]:
                    if otherFile[-4:] == '.pdb':
                        topology = otherFile

                        if 'recentered' in trajectory:
                            pass
                        else:
                            try:
                                recenterDCD(topology, aaa) # never actually recenter, we don't need it and sometimes it breaks everything
                                runPairs.append([topology,trajectory.split('.')[0] + '_recentered.dcd'])
                            except:
                                copyfile(trajectory, trajectory.split('.')[0] + '_recentered.dcd')

                        break

# combine folding trajectories
mmbDirs = []
for dirs in dirList:
    if 'mmbFiles' in dirs:
        mmbDirs.append(dirs)

for i in range(len(mmbDirs)):
    mmbdir = mmbDirs[i]
    copyfile(mmbdir + '/last.1.pdb','foldFrame_{}.pdb'.format(i))
    mmbdirList = os.listdir(mmbdir)

    filenames = []
    for file in mmbdirList:
        if 'trajectory' in file:
            filenames.append(mmbdir + '/' +file)

    with open('foldingTraj_%d'%i + '.pdb', 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

    replaceText('foldingTraj_%d'%i + '.pdb', '*', "'")  # due to a bug in this version of MMB - structures are encoded improperly - this fixes it

    u = mda.Universe('foldingTraj_%d'%i + '.pdb')
    with mda.Writer('foldingTraj_%d'%i +'.dcd', u.atoms.n_atoms) as W:
        for ts in u.trajectory:
            W.write(u)

#====================

nFolds = len(mmbDirs)
nDocks = sum(['top_' in dir for dir in dirList])

for i in range(nFolds):
    trajectories = []
    trajectories.append('foldingTraj_{}.dcd'.format(i)) # folding trajectory
    trajectories.append('clean_foldedSequence_{}_processed_trajectory_recentered.dcd'.format(i))
    trajectories.append('clean_relaxedSequence_{}_processed_complete_trajectory_recentered.dcd'.format(i))

    u = mda.Universe('foldFrame_{}.pdb'.format(i), trajectories)

    with mda.Writer('fullPipeTraj_{}.dcd'.format(i), u.atoms.n_atoms) as W:
        for ts in u.trajectory:
            W.write(u)


complexes = []
for i in range(nFolds):
    for j in range(nDocks):
        # combine all the complex trajectories
        complexes.append('clean_complex_{}_{}_processed_complete_trajectory_recentered.dcd'.format(i,j))


os.system('cp clean_complex_0_0_processed.pdb ./complexTopology.pdb')
u = mda.Universe('complexTopology.pdb', complexes)

with mda.Writer('combinedComplexTraj.dcd', u.atoms.n_atoms) as W:
    for ts in u.trajectory:
        W.write(u)


os.mkdir('fragments')
os.system('mv clean*  fragments')
