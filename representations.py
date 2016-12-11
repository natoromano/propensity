"""
Learns representation on EHR data using TensforFlow, and saves them.

Does this for many different dimensions if necessary.

Nathanael Romano
"""

import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt

from utils.data import load_data
from utils.nn import (
    add_fully_connected, 
    precision_recall,
    add_criteria,
    extract_activations
)

if __name__ == "__main__":
    # Load data
    datasets = load_data()
    dim = datasets.dimension
    
    # Instantiate session and create base variables
    sess = tf.InteractiveSession()
    x = tf.placeholder(tf.float32, shape=[None, dim])
    y_ = tf.placeholder(tf.float32, shape=[None, 2])
    
    # Hidden layers dimensions
    hidden1 = 500
    hidden2 = 50
    hidden3 = outDim

    # Actual network creation
    h1 = add_fully_connected(x, dim, hidden1, relu=False)
    h2 = add_fully_connected(h1, hidden1, hidden2, tanh=True)
    h3 = add_fully_connected(h1, hidden1, hidden3, tanh=False)
    y_scores = add_fully_connected(h3, hidden3, 2, tanh=False)

    cross_entropy, train_step, \ 
    correct_prediction, accuracy = add_criteria(y_, y_scores, 7.2e-7)
    
    # Actually train data
    sess.run(tf.initialize_all_variables())

    for i in range(150000):
        batch = datasets.train.next_batch(512)
        if i % 1000 == 0:
            val_accuracy = accuracy.eval(
                feed_dict={x: datasets.val1._patients, 
                           y_: datasets.val1._labels}
            )
            if verbose >= 2:
                print "Step %d:" % i, val_accuracy

        train_step.run(feed_dict={x: batch[0], y_: batch[1]})

    if verbose >= 1:
        print "Validation accuracy for dim %d is %g" % (
            accuracy.eval(feed_dict={x: datasets.val1._patients,
                                     y_: datasets.val1._labels}),
            outDim
        )
    
    # Extract activations
    extract_activations(
        x, y_, y_scores, datasets.val2, "val2_{}_{}".format(outDim)
    )
    extract_activations(
        x, y_, y_scores, datasets.test, "test_{}_{}".format(outDim)
    )
