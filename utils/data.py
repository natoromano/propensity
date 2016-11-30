"""
Defines a few utility functions for data
load and processing.
"""

import numpy as np

DEFAULT_PATH = "/scratch/users/naromano/OSIM2"


def load_data(path=None):
    """
    Loads datasets at provided path, in a dataset object.
    """
    if path is None:
        path = DEFAULT_PATH
        
    dtrain = Dataset(path + "/dtrain.txt")
    dval1 = Dataset(path + "/dval1.txt")
    # Don't touch!
    # dval2 = Dataset(path + "/dval2.txt")
    # dtest = Dataset(path + "/dtest.txt")

    return Datasets(train=dtrain, val1=dval1)



class Datasets(object):
    """
    Simple container for datasets.
    """
    
    def __init__(self, **kwargs):
        for argname, dataset in kwargs.iteritems():
            setattr(self, argname, dataset)



class Dataset(object):
    """
    Simple dataset with patients and treatments.
    Offers a next_batch() method.
    """

    def __init__(self, path):
        data = np.loadtxt(
            path, 
            skiprows=1, 
            usecols=range(704)[:-4]
        )
        self._patients = data[:, :-1]
        self._labels = data[:, -1]
        self._num_patients = data.shape[0]
        self._index_in_epoch = 0
        self._epochs_completed = 0


    def next_batch(self, batch_size):
        start = self._index_in_epoch
        self._index_in_epoch += batch_size
 
        if self._index_in_epoch > self._num_patients:
            # Finished epoch
            self._epochs_completed += 1
            # Shuffle the data
            perm = numpy.arange(self._num_patients)
            numpy.random.shuffle(perm)
            self._patients = self._patients[perm]
            self._labels = self._labels[perm]
            # Start next epoch
            start = 0
            self._index_in_epoch = batch_size
            assert batch_size <= self._num_examples

        end = self._index_in_epoch
        return self._patients[start:end], self._labels[start:end]
