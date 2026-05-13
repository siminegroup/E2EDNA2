# The active learning codebase is adapted from: https://github.com/AntoninDuval/Active-Learning-For-Drug-Discovery

import os, joblib
import numpy as np
# import itertools
import pandas as pd

from acquirer.acquirer import *
from model.model import *
from molecule_pool.molecule_pool import DofMaskPool

# import seaborn as sns
# import matplotlib.pyplot as plt

# integer internal coordinates (multiplied by 10 to integerize)
data_bound = np.load(os.path.join('<path_to_file>', '<integerized_internal_coords_dna_backbone_bound_state>.npy'), allow_pickle=False)
data_free = np.load(os.path.join('<path_to_file>', '<integerized_internal_coords_dna_backbone_free_state>.npy'), allow_pickle=False)

init_data = np.load('init_100_data.npy', allow_pickle=False)
init_target = np.load('init_100_target.npy', allow_pickle=False)
molecule_pool = DofMaskPool(init_data, init_target)
num_dof = 633  # number of internal coordinates of the dna backbone atoms

BATCH_SIZE = 10000 # randomly sample data points from which find pormising candidates
MAX_ITERATIONS = 2000
# NB_EXPERIMENTS = 1
k = 100 # selective, only choose the top k candidates based on surrogate model prediction

lig_index = 0 # The first ligand: the strongest binder to the DNA, as the test system
experiment = f'2000iter_k100_minimize_deltaCompression_ligand{lig_index}'        

# models = [RandomForest(n_estimators=100, max_depth=8)]
models = [RandomForest()]

acquisition_functions = [Greedy(batch_size=k, require_var=True)]

objective_function = DeltaCompression(data_bound, data_free)

df_results = pd.DataFrame(columns=["experiment", "model", "batch_size", "acquisition_function"])


# record the train/test loss/variance of the RF model given each batch of 50 data points
train_loss = np.zeros(MAX_ITERATIONS)
train_var = np.zeros(MAX_ITERATIONS)
test_loss = np.zeros(MAX_ITERATIONS)
test_var = np.zeros(MAX_ITERATIONS)

verbose = True

# for i, model in enumerate(models):
model = models[0]
# for j, acquisition_function in enumerate(acquisition_functions):
acquisition_function = acquisition_functions[0]

print('Model used  :', model.name)
print('Acquisition function used  :', acquisition_function.name)
print('Batch Size  :', BATCH_SIZE)

# for experiment in range(NB_EXPERIMENTS):
# experiment = ''
print('EXPERIMENT N°', experiment)
# ratio_found = []

print('Initialize training set...')
# train_set, test_set = molecule_pool.initialize_batch(batch_size=BATCH_SIZE) # randomly select BATCH_SIZE data points
# For now, use all the evaluated data that we have
train_set = molecule_pool

# # Get the index of top k selections according to the evaluated compression size: in the whole dataset or in the train set
# top_k_selection_idx = molecule_pool.sort_idx_by_true_score()[:k]
# # top_k_mol = set(molecule_pool.df[idx_best, 0])
# top_k_found_in_train_set = train_set.get_top_k(k, top_k_selection_idx)

if verbose:
    # print(f"% of top {k} molecules found in the train set: {(len(top_k_found_in_train_set)/k)*100}%")
    print('='*50)

# if model.name == 'NN' and acquisition_function.name == 'UCB':
#     print('Variance not yet implemented for NN')
#     continue

iteration = 0

result = {'experiment':f'Experiment_{experiment}',
            'model': model.name,
            'batch_size': BATCH_SIZE,
            'acquisition_function' : acquisition_function.name}#,
        #  len(train_set.arr):(len(top_k_found_in_train_set)/k)*100}

# save the trainging set: data and targets, as they are the ground truth.
BO_params = {'BATCH_SIZE': BATCH_SIZE, 'MAX_ITERATIONS': MAX_ITERATIONS, 'k': 50, 'experiment_index': 1, 'lig_index': lig_index}
try:
    os.mkdir(f'Exp_N_{experiment}')
    np.save(f'./Exp_N_{experiment}/train_dof_mask.npy', train_set.data, allow_pickle=False)
    np.save(f'./Exp_N_{experiment}/train_target.npy', train_set.target, allow_pickle=False)

    np.save(f'./Exp_N_{experiment}/train_loss.npy', train_loss, allow_pickle=False)
    np.save(f'./Exp_N_{experiment}/train_var.npy', train_var, allow_pickle=False)
    np.save(f'./Exp_N_{experiment}/test_loss.npy', test_loss, allow_pickle=False)
    np.save(f'./Exp_N_{experiment}/test_var.npy', test_var, allow_pickle=False)

    metadata = {"model": model.model,
                "params": model.model.get_params(),
                "BO_params": BO_params
                }
    # Save everything in one file
    joblib.dump(metadata, f"./Exp_N_{experiment}/random_forest_with_metadata.pkl")
except:
    print('Error: cannot create a folder. Folder may have existed.')


while iteration < MAX_ITERATIONS:
    if verbose:
        print(f'ITERATION: {iteration} out of {MAX_ITERATIONS}')
        print('Train set shape : ', train_set.data.shape)
        print('Training the model...')
    
    print('train_data:', train_set.data.shape, train_set.data.dtype)
    # print('train_target:', train_set.target)
    model.train(train_set)

    print('Model trained. Predicting the score on the train set...')
    # record the train loss using MSE: the loss of the model on the train set
    model.predict(train_set, require_var=True, verbose=False)
    train_loss[iteration] = np.mean((train_set.score - train_set.target)**2)
    train_var[iteration]  = np.mean(train_set.variance)

    # find the batch of data from test set to be added to train set:
    # try a bunch test data
    # randomly generate other permutation of the selection mask. Assign -1 as the target
    # may add a seed here
    print(f'Generating random candidate data ({BATCH_SIZE}x{num_dof})...')
    candidate_set_data = np.random.randint(0, 2, (BATCH_SIZE, num_dof), dtype=np.int8) # [0,2) = {0,1}
    # check for unique data points
    candidate_set_data = np.unique(candidate_set_data, axis=0)

    # # Keep only the rows in candidate_set that are NOT in train_set.
    # # candidate_set_data = np.setdiff1d(candidate_set_data, train_set.data, assume_unique=True, axis=0)
    # print('Removing the data points that are already in the train set...')
    # mask = np.any(np.all(candidate_set_data[:, None] == train_set.data, axis=2), axis=1)
    # candidate_set_data = candidate_set_data[~mask]
    # # OOM WARNING: `np.all()` step will temporarily create an integer array of (10000, train_set_rows, 633). 5000 rows in train_set will cause 29 GB => OOM 
    # # Solution: use `np.all()` to compare the promising candidates (much less) with training set. Then deal with array of (50, train_set_rows, 633).
    # # Better solution: use concatenate two datasets and np.unique():
    combined = np.vstack((train_set.data, candidate_set_data)) # must put candidate_set_data at the bottom
    _, indices_unique = np.unique(combined, axis=0, return_index=True) # if excluding the last column, use: combined[:,:-1]
    mask = indices_unique[indices_unique >= train_set.data.shape[0]] - train_set.data.shape[0] # get the indices of the unique rows in candidate_set_data
    candidate_set_data = candidate_set_data[mask]
    print('Size of pre-selected candidate data (novel from training set):', candidate_set_data.shape)
    # assign -1 as the target. can evaluate this if want to monitor how well the model selecting the candidate data points
    candidate_set_target = np.zeros(candidate_set_data.shape[0], dtype=np.int8) - 1 # just dummy values, no use
    # assemble a candidate_set object
    # candidate_set = DofMaskPool(np.hstack((candidate_set_data, candidate_set_target.reshape(BATCH_SIZE,1)))
    candidate_set = DofMaskPool(data_arr=candidate_set_data, target_arr=candidate_set_target)

    print('Predicting the score...')
    # print(acquisition_function.require_var)
    score = model.predict(candidate_set, acquisition_function.require_var, True)
    # `model.predict()` also assigns candidate_set.score with `score` under the hood
    print(f'Max, min, mean, std of the predicted score: {score.max():.4f}, {score.min():.4f}, {score.mean():.4f}, {score.std():.4f}')

    print(f'Selecting the {k} most promising DOF_masks (more negative predicted deltaCompression), evaluating their target values (deltaCompression) & adding to train set...')
    most_promising_mol_set = acquisition_function.select_train_set(candidate_set)

    '''# I don't need to split the evaluated data into train & test, unless I want to evaluate the loss.
    I can just use all the evaluated data to train the model, then sample candidate_set data for surrogate model,
    and find the promising data to be evaluated by oracle and added to the dataset
    '''
    # evaluate the objective function for the newly selected data (ie, query oracle, ie, run compression)
    # and use them to update the overall dataset's target (ground truth)
    new_targets = objective_function.evaluate(most_promising_mol_set, lig_index)
    # print('Newly evaluated target: size and values:', new_targets.shape, new_targets)

    # print('Adding the evaluated target values')
    most_promising_mol_set.update_target(new_targets)
    # save the score for the candidate_set or the selected top-k candidates? Purpose: compare with the true target values to evaluate the test loss.

    # record the test loss using MSE: the loss of the model on the next batch of data points
    test_loss[iteration] = np.mean((most_promising_mol_set.score - new_targets)**2)
    test_var[iteration]  = np.mean(most_promising_mol_set.variance)

    # merge with current trian_set
    # train_set = DofMaskPool(np.vstack((train_set.arr, most_promising_mol_set.arr)))
    new_target_arr = np.hstack((train_set.target, most_promising_mol_set.target))
    train_set = DofMaskPool(data_arr=np.vstack((train_set.data, most_promising_mol_set.data)),
                            target_arr=new_target_arr)
    n = 10
    # nth_largest = np.sort(np.partition(new_target_arr, -n)[-n:])[::-1]
    # print(f'>> Latest 10 highest target value: {nth_largest}')
    nth_largest = np.sort(np.partition(new_target_arr, n)[:n])                                
    print(f'>> Latest {n} lowest target value: {nth_largest}')  

    # # # merge the index of the new train data with the previous one
    # # new_train_mol_indx = np.concatenate((train_set.arr[:, 0], most_promising_mol_set.arr[:, 0]))
    # # train_set, test_set = molecule_pool.create_batch(new_train_mol_indx)

    iteration += 1

    # # top_k_found = train_set.get_top_k(k, top_k_mol)
    # top_k_found_in_train_set = train_set.get_top_k(k, top_k_selection_idx)

    if verbose:
        # print('Get top k molecules...')
        # print(f"% of top {k} molecules found in the train set: {(len(top_k_found_in_train_set)/k)*100}%")
        print('='*50)
    
    # result[len(train_set.arr)] =  (len(top_k_found_in_train_set)/k)*100

    # save the trainging set: data and targets, as they are the ground truth.
    if iteration % 3 == 0 or iteration == MAX_ITERATIONS:
        try:
            np.save(f'./Exp_N_{experiment}/train_dof_mask.npy', train_set.data, allow_pickle=False)
            np.save(f'./Exp_N_{experiment}/train_target.npy', train_set.target, allow_pickle=False)

            np.save(f'./Exp_N_{experiment}/train_loss.npy', train_loss, allow_pickle=False)
            np.save(f'./Exp_N_{experiment}/train_var.npy', train_var, allow_pickle=False)
            np.save(f'./Exp_N_{experiment}/test_loss.npy', test_loss, allow_pickle=False)
            np.save(f'./Exp_N_{experiment}/test_var.npy', test_var, allow_pickle=False)

            metadata = {"model": model.model,
                        "params": model.model.get_params(),
                        "BO_params": BO_params
                        }
            # Save everything in one file
            joblib.dump(metadata, f"./Exp_N_{experiment}/random_forest_with_metadata.pkl")
        except:
            print('>> Error when saving checkpoint files.')
            # experiment += 1
            # save again here?

# # add the results of the experiment
# df_results = df_results.append(result, ignore_index=True)
# # df_results.to_csv('./result.csv', index=False)