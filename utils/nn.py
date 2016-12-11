"""
Defines a few utility functions for neural networks,
using TensorFlow.
"""

import tensorflow as tf
import numpy as np


def add_fully_connected(x, inshape, outshape, relu=False, tanh=False):
    """
    Creates a fully connected layer with specified in and 
    out shape.
    """
    W = create_weights([inshape, outshape])
    b = create_bias([outshape])
    if relu:
        return tf.nn.relu(tf.matmul(x, W) + b)
    elif tanh:
        return tf.nn.tanh(tf.matmul(x, W) + b)
    else:
        return tf.matmul(x, W) + b


def create_weights(shape, stdDev=0.25):
    """
    Creates a TF weight matrix.
    """
    init = tf.truncated_normal(shape, stddev=stdDev)
    return tf.Variable(init)


def create_bias(shape, defaultVal=0.1):
    """
    Creates a TF bias vector.
    """
    init = tf.constant(defaultVal, shape=shape)
    return tf.Variable(init)


def conv2d(x, W):
    """
    Creates a simple convolution with specified weights.
    """
    return tf.nn.conv2d(
        x, W, strides=[1, 1, 1, 1], padding='SAME'
    )


def max_pool_2x2(x):
    """
    Creates a simple max pooling layer.
    """
    return tf.nn.max_pool(
        x, ksize=[1, 2, 2, 1], strides=[1, 2, 2, 1], padding='SAME'
    )


def add_criteria(y_, y_scores, step=1e-5):
    """
    Adds the last layer of a neural network, using a simple
    gradient descent and cross-entropy.
    """
    cross_entropy = tf.reduce_mean(
        tf.nn.softmax_cross_entropy_with_logits(y_scores, y_)
    )
    train_step = tf.train.AdamOptimizer(step).minimize(cross_entropy)

    correct_prediction = tf.equal(tf.argmax(y_scores, 1), tf.argmax(y_, 1))
    accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
    return cross_entropy, train_step, correct_prediction, accuracy


def precision_recall(x, y, scores, dataset):
    """
    Computes precision recall for given model, on the
    provided dataset
    """
    predictions = tf.argmax(scores, 1).eval(
        feed_dict={x: dataset._patients,
                   y: dataset._labels}
    )
    actual = np.argmax(dataset._labels, axis=1)
    fn = float(sum((1 - np.equal(predictions, actual)) * np.equal(predictions, 0)))
    fp = float(sum((1 - np.equal(predictions, actual)) * np.equal(predictions, 1)))
    tp = float(sum(np.equal(predictions, actual) * np.equal(predictions, 1)))
    return round(tp / (fp + tp), 2), round(tp / (fn + tp), 2)


def extract_activations(x, y, layer, dataset, name):
    """
    Extracts layer's activations from a trained network, on
    the provided dataset.
    """
    activations = layer.eval(
        feed_dict={x: dataset._patients,
                   y: dataset._labels}  
    )
    fileName = "/scratch/users/naromano/OSIM2/layer" + name + ".txt"
    np.savetxt(fileName, activations)
