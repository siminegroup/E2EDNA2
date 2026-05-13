from __future__ import annotations

from typing import List

import numpy as np
from sklearn.preprocessing import StandardScaler


class DofMaskPool(object):

    def __init__(self, data_arr, target_arr):
        # self.arr = arr
        # self.target = self.arr[:, -1].copy()
        # # self.data   = self.arr[:, 1:-1].copy() # the first column is the index
        # self.data   = self.arr[:, :-1].copy()
        self.data = data_arr.copy()
        self.target = target_arr.copy()
        # use a panda df will be better: store different data types

        self.score = None
        self.variance = None

    # Need to update the following methods to work with the new data structure
    # def initialize_batch(self, batch_size: int) -> (DofMaskPool, DofMaskPool):
    #     """
    #     Split the dataset into a random train set of size batch_size and a test set.
    #     :param batch_size:
    #     :return:
    #     """
    #     init_train_index = np.random.choice(len(self.arr), size=batch_size, replace=False)
    #     train = DofMaskPool(self.arr[init_train_index, :])
    #     test = DofMaskPool(np.delete(self.arr, init_train_index, axis=0))
    #     return train, test

    # def create_batch(self, train_index) -> (DofMaskPool, DofMaskPool):
    #     """
    #     Create a new training and testing set from a list of indexes
    #     :param train_index:
    #     :return:
    #     """
    #     train = DofMaskPool(self.arr[train_index,:])
    #     test = DofMaskPool(np.delete(self.arr, train_index, axis=0))

    #     return train, test

    # def get_top_k(self, k, top_k_indx) -> set:
    #     """
    #     Return the top k selections found in the current dataset matching the top k molecules in the whole dataset
    #     :param k:
    #     :param top_k:
    #     :return:
    #     """
    #     # the first column is the index
    #     top_k_current_dataset  = self.sort_idx_by_true_score()[:k]
    #     current_dataset_top_k_indx = self.arr[top_k_current_dataset, 0] # this returns the original index in the whole dataset

    #     overall_top_k_in_current_dataset = set(current_dataset_top_k_indx).intersection(top_k_indx)
    #     return overall_top_k_in_current_dataset

    def preprocess_data(self):
        """
        Process the data by scaling the features and removing categorical variables
        :return:
        """
        scaler = StandardScaler() # Standardize features by removing the mean and scaling to unit variance.
        preprocessed_data = scaler.fit_transform(self.data)
        return preprocessed_data

    def add_score(self, score):
        """
        Add the predicted score
        :param score:
        :return:
        """
        self.score = score
    
    def update_target(self, new_targets):
        """
        Add the evaluated target values
        :param new_targets:
        :return:
        """
        self.target = new_targets
        # self.arr[:,-1] = new_targets

    def add_variance(self, variance):
        """
        Add the variance
        :param variance:
        :return:
        """
        self.variance = variance

    def sort_idx_by_true_score(self) -> List:
        """
        Return the sorted index by true docking score
        :return:
        """
        return self.target.argsort()

    # def sort_idx_best_preds(self) -> List:
    #     """
    #     Return the sorted index by predicted score from max to min (descending order)
    #     :return:
    #     """
    #     return self.score.argsort()[::-1] # top k molecules with the highest predicted score
    
    def sort_idx_best_preds(self) -> List:
        """
        Return the sorted index by predicted score from max to min (ascending order)
        :return:
        """
        return self.score.argsort()  # top k molecules with the lowest predicted score: push the deltaCompression to negative infinity
