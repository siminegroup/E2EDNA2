# import sys
# sys.path.append("../src")

from abc import ABC, abstractmethod
from molecule_pool.molecule_pool import DofMaskPool
import numpy as np


class Acquirer(ABC):

    def __init__(self, name, batch_size):
        self.name = name
        self.batch_size = batch_size
        self.require_var = False

    @abstractmethod
    def select_train_set(self, moleculepool):
        pass


# class RandomSearch(Acquirer):

#     def __init__(self, batch_size):
#         super().__init__("RandomSearch", batch_size)
#         self.require_var = False

#     def select_train_set(self, moleculepool: MoleculePool) -> MoleculePool:
#         """
#         Return a random subset of the molecule dataset of size batch size
#         :param moleculepool:
#         :return:
#         """
#         train_idx = np.random.choice(len(moleculepool.df), size=self.batch_size, replace=False)
#         train_set = MoleculePool(moleculepool.df[train_idx])
#         return train_set


class Greedy(Acquirer):
    def __init__(self, batch_size, require_var):
        super().__init__("Greedy", batch_size)
        self.require_var = require_var

    def select_train_set(self, input_set: DofMaskPool) -> DofMaskPool:
        """
        Select the top molecules that have the highest predicted score.
        :param input_set:
        :param batch_size:
        :return:
        """
        idx_best_preds = input_set.sort_idx_best_preds()[:self.batch_size]
        output_set = DofMaskPool(input_set.data[idx_best_preds,:], input_set.target[idx_best_preds]) # target values are dummy values
        output_set.score = input_set.score[idx_best_preds].copy()
        if input_set.variance is not None:
            output_set.variance = input_set.variance[idx_best_preds].copy()
        return output_set


# class UCB(Acquirer):
#     def __init__(self, batch_size, beta=2):
#         super().__init__("UCB", batch_size)
#         self.beta = beta
#         self.dict_ = {}
#         self.require_var = True

#     def select_train_set(self, moleculepool: MoleculePool) -> MoleculePool:
#         """
#         Select the top molecules that have the highest UCB score, which is a combination of exploration and exploitation
#         :param moleculepool:
#         :param batch_size:
#         :return:
#         """
#         ucb_score = moleculepool.score + self.beta*np.sqrt(moleculepool.variance)
#         index_sorted = np.argsort(ucb_score)[:self.batch_size]
#         return MoleculePool(moleculepool.df[index_sorted])