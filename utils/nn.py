"""
Defines a few utility functions for neural networks,
using TensorFlow.
"""

import tensorflow as tf


def add_fully_connected(x, inshape, outshape, relu=False):
    """
    Creates a fully connected layer with specified in and 
    out shape.
    """
    W = create_weights([inshape, outshape])
    b = create_bias([outshape])
    if relu:
        return tf.nn.relu(tf.matmul(x, W) + b)
    else:
        return tf.matmul(x, W) + b


def create_weights(shape, stdDev=0.1):
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
