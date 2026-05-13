from abc import ABC
from sklearn.ensemble import RandomForestRegressor
import gzip
# import torch
# import torch.nn as nn
# from torch.optim import Adam
# from torch.utils.data import DataLoader
import tqdm
from molecule_pool.molecule_pool import *
# from network.mlp import MLP

class Model(ABC):
    def __init__(self, name):
        self.name = name

    # def predict(self, moleculepool: MoleculePool, require_var: bool):
    #     pass

    # def train(self, moleculepool: MoleculePool):
        # pass

# objective function
class DeltaCompression(Model):
    def __init__(self, data_bound, data_free, name="Lossless_Compression"):
        # data to be compressed from bound & free states.
        # it's integer internal coordinate data in this case.
        self.data_bound = data_bound
        self.data_free = data_free
        self.name = name
        
    def calc_deltaS(self, dof_mask, lig_index):
        # compression size in bit
        C_each_frame_bound = 8*np.array([len(gzip.compress(np.ndarray.flatten(self.data_bound[lig_index, j, dof_mask]), compresslevel=9)) for j in range(self.data_bound.shape[1])])
        C_each_frame_free  = 8*np.array([len(gzip.compress(np.ndarray.flatten(self.data_free[j, dof_mask]), compresslevel=9)) for j in range(self.data_free.shape[0])])

        obj_fun = np.mean(C_each_frame_bound) - np.mean(C_each_frame_free)
        return obj_fun

    def evaluate(self, data_set, lig_index):

        batch_size = data_set.data.shape[0]
        delta_compression_size = np.zeros(batch_size)

        for i in tqdm.tqdm(range(batch_size)):
            dof_mask = data_set.data[i,:].astype(bool) # covnert {0,1} array to boolean array!
            
            delta_compression_size[i] = self.calc_deltaS(dof_mask, lig_index)

        return delta_compression_size


# surrogate model
class RandomForest(Model):
    # def __init__(self, **kwargs):
    def __init__(self, n_estimators=100, max_depth=8):
        super().__init__("RandomForestRegressor") # call the superclass's constructor __init__()
        # print(kwargs)
        self.model = RandomForestRegressor(n_estimators=n_estimators, max_depth=max_depth)
        # print(self.model)

    def train(self, moleculepool, verbose: bool = False):
        """

        :param moleculepool:
        :param verbose:
        :return:
        """
        # does it make sense to normalize the data from {0,1}?
        # data = moleculepool.preprocess_data() # masks in the training set. The compression of those mask options are evaluated by compression    
        data   = moleculepool.data
        target = moleculepool.target # compression size of the training masks

        # print('start training...')
        self.model.fit(data, target)
        # print('Done training')
        if verbose: print('R2 score on train: ', self.model.score(data, target))
        

    def predict(self, test_set, require_var, verbose=False):
        """
        Predict the score using the RF model on the test set. Compute the variance by getting
        the prediction of each tree.
        :param test_set:
        :param require_var:
        :param verbose:
        :return:
        """
        preds = np.zeros((len(test_set.data), len(self.model.estimators_)))

        # the mask options that have not been compressed
        # test_prepro = test_set.preprocess_data() # does it make sense to normalize the data from {0,1}?
        test_prepro = test_set.data.copy()

        if verbose: print('Model R2 score: ', self.model.score(test_prepro, test_set.target))

        score = self.model.predict(test_prepro)

        if require_var:
            # print('Computing the variance...')
            for j, submodel in enumerate(self.model.estimators_):
                preds[:, j] = submodel.predict(test_prepro)
            test_set.add_variance(np.var(preds, axis=1))

        test_set.add_score(score)
        return score


# class NN(Model):

#     def __init__(self, param, epoch):
#         super().__init__("NN")
#         self.model = MLP(**param).double()
#         self.optimiser = Adam(self.model.parameters(), lr=0.01, weight_decay=0.01)
#         self.epoch = epoch
#         self.criterion = nn.MSELoss()

#     def train(self, moleculepool, verbose: bool = False):
#         """

#         :param moleculepool:
#         :param verbose:
#         :return:
#         """
#         data = moleculepool.preprocess_data().astype(float)
#         inputs = torch.tensor(np.concatenate([data, moleculepool.target.reshape(-1, 1).astype(float)], axis=1))

#         dataloader = DataLoader(
#             inputs,
#             batch_size=4096,
#             shuffle=False,
#         )

#         for _ in range(self.epoch):
#             for i, batch in enumerate(dataloader):
#                 self.optimiser.zero_grad()
#                 x = batch[:, :-1]
#                 y = batch[:, -1]
#                 y = self._normalize(y)
#                 preds = self.model(x.double())
#                 loss = self.criterion(preds.squeeze(), y)
#                 loss.backward()
#                 self.optimiser.step()
#         if verbose:
#             print(loss)

#     def predict(self, test_set: MoleculePool, require_var=False) -> np.array:
#         """

#         :param test_set:
#         :param require_var:
#         :return:
#         """
#         test_prepro = torch.tensor(test_set.preprocess_data().astype(float))
#         with torch.no_grad():
#             score = self.model(test_prepro.detach())

#         score = score * self.std + self.mean
#         test_set.add_score(score)
#         return score

#     def reset_params(self):
#         for layer in self.model.children():
#             if hasattr(layer, 'reset_parameters'):
#                 layer.reset_parameters()

#     def _normalize(self, target):
#         """
#         Normalize the target to make the training of the nn easier
#         :param target:
#         :return:
#         """
#         self.mean = np.nanmean(target)
#         self.std = np.nanstd(target)
#         return (target - self.mean) / self.std


