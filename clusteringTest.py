from utils import *
import sklearn.cluster as cluster

#sequence = numbers2letters(np.random.randint(0,4,size=(1,80)))[0]
sequence = 'AATTCCACAGAAGGGTTGAGACCGGTTCGTTCGGTTCGAGGAACCCCTGACTCGGCCAATATAAAAGTTGGTTAGAAAGT' # sequence with tons of suboptimal structures
dict, subopts = nupackAnalysis(sequence,310,0.163)

lists = []
configs = []
strings = []
for i in range(len(subopts[0])):
    strings.append(str(subopts[0][i].structure))
    lists.append(ssToList(strings[-1]))
    configs.append(pairListToConfig(lists[-1],sequence))


dists = getSecondaryStructureDistance(configs)

# ok but I'd like to have a threshold
#kmeans = cluster.KMeans(n_clusters = 4, random_state=0).fit(dists) # kmeans cluster
#labels = kmeans.labels_

# nice, nice - getting somewhere - this one uses a minimum threshold rather than a maximum
# then algorithm depends on quality of distance metric
converged = False
threshold = 0.1 # threshold (on normalized distance metric
while not converged:
    agglomerate = cluster.AgglomerativeClustering(n_clusters=None,affinity='precomputed',linkage='average',compute_full_tree=True,distance_threshold=threshold).fit(dists)
    labels = agglomerate.labels_
    nClusters = agglomerate.n_clusters_

    clusters = []
    clusterProbs = []
    clusterDists = []
    probSums = np.zeros(len(np.unique(labels)))
    for i in range(len(np.unique(labels))):
        inds = np.where(labels == i)[0].astype(int)
        clusters.append([strings[j] for j in inds])
        clusterProbs.append([subopts[1][j] for j in inds])
        probSums[i] = np.sum(clusterProbs[-1]) # unnormalized
        clusterDists.append(getSecondaryStructureDistance([configs[j] for j in inds]))

    normedProbSums = probSums / np.sum(probSums)

    if np.average(normedProbSums) < 0.1: # if the average normed probability of our clusters is less than 10%, we're not converged
        threshold *= 1.1 # if it isn't converged, boost the threshold
    else:
        converged = True

# find a representative from each cluster
weightedAvgDistance = []
clusterReps = []
for i in range(nClusters):
    weightedAvgDistance.append((clusterDists[i] + np.eye(len(clusterDists[i]))) @ (1/np.asarray(clusterProbs[i]))) # weight the distance to each config against its probability of being seen
    clusterReps.append(clusters[i][np.argmin(weightedAvgDistance[i])])
