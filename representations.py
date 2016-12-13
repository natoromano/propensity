"""
Learns representation on EHR data using TensforFlow, and saves them.
The representations are learned training a simple fully
connected neural network on the data, using the treatment
assignment as target, and by extracting the its last
layer, which ought to be low-dimensional.

Does this for many different dimensions if necessary.

Nathanael Romano
"""

import argparse

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
    parser = argparse.ArgumentParser(description='Parameters')
    parser.add_argument('-d', '--dims', nargs='+', type=int, default=range(2, 11))
    parser.add_argument('-n', '--numIter', type=int, default=150000)
    parser.add_argument('-m', '--mu', type=float, default=7.2e-7)
    parser.add_argument('-b', '--batchSize', type=int, default=512)
    parser.add_argument('--h1', type=int, default=500)
    parser.add_argument('--h2', type=int, default=50)
    parser.add_argument('-v', '--verbose', type=int, default=0)
    args = parser.parse_args()
    
    # Load data
    if args.verbose >= 1:
        print "Loading data..."
        print
    datasets = load_data()
    dim = datasets.dimension
    
    for outDim in args.dims:
        if args.verbose >= 1:
            print "Starting to learn representation for d = %d..." % (outDim)
            
        # Instantiate session and create base variables
        tf.reset_default_graph()
        sess = tf.InteractiveSession()
        x = tf.placeholder(tf.float32, shape=[None, dim])
        y_ = tf.placeholder(tf.float32, shape=[None, 2])

        # Hidden layers dimensions
        hidden1 = args.h1
        hidden2 = args.h2
        hidden3 = outDim

        # Actual network creation
        h1 = add_fully_connected(x, dim, hidden1, relu=False)
        h2 = add_fully_connected(h1, hidden1, hidden2, tanh=True)
        h3 = add_fully_connected(h1, hidden1, hidden3, tanh=False)
        y_scores = add_fully_connected(h3, hidden3, 2, tanh=False)

        cross_entropy, train_step, correct_prediction, accuracy = add_criteria(y_, y_scores, args.mu)

        # Actually train data
        sess.run(tf.initialize_all_variables())

        for i in range(args.numIter):
            batch = datasets.train.next_batch(args.batchSize)
            if i % 1000 == 0:
                val_accuracy = accuracy.eval(
                    feed_dict={x: datasets.val1._patients, 
                               y_: datasets.val1._labels}
                )
                if args.verbose >= 2:
                    print "Step %d:" % i, val_accuracy

            train_step.run(feed_dict={x: batch[0], y_: batch[1]})

        if args.verbose >= 1:
            print "Accuracyy(d=%d): %g" % (
                outDim,
                accuracy.eval(feed_dict={x: datasets.val1._patients,
                                         y_: datasets.val1._labels})
            )
            print

        # Extract activations
        extract_activations(
            x, y_, h3, datasets.val2, "val2_{}".format(outDim)
        )
        extract_activations(
            x, y_, h3, datasets.test, "test_{}".format(outDim)
        )
        
        sess.close()

