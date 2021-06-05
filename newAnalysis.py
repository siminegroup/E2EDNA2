from utils import *
import os
os.chdir('C:/Users/mikem/Desktop/mmruns/run47')
u = mda.Universe('cleanaptamer.pdb','cleanaptamer.dcd')

class atomDistances(AnalysisBase):
    '''
    calculate watson-crick bond lengths
    '''
    def __init__(self, atomgroups, **kwargs):
        '''
        :param atomgroups: a list of atomgroups for which the wc lengths are calculated
        '''
        super(atomDistances,self).__init__(atomgroups[0].universe.trajectory,**kwargs)
        self.atomgroups = atomgroups
        self.ag1 = atomgroups[0]
        self.ag2 = atomgroups[1]

    def _prepare(self):
        self.results = []

    def _single_frame(self):
        distance = calc_bonds(self.ag1.positions,self.ag2.positions,box=self.ag1.dimensions)
        self.results.append(distance)

    def _conclude(self):
        self.results = np.asarray(self.results)


aind = []
bind = []
for i in range(50):
    for j in range(50):
        aind.append(i)
        bind.append(j)

a = mda.AtomGroup(aind,u)
b = mda.AtomGroup(bind,u)

na = atomDistances([a,b]).run()
traj = na.results