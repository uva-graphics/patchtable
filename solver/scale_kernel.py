
"""
Functions to find optimal scale and error between an 'output' (approximated) kernel and a target (ground truth) kernel.
"""

import numpy

def optimal_scale(out, target):
    out = out.flatten()
    target = target.flatten()
    assert len(out) == len(target)
    num = 0.0
    denom = 0.0
    for i in range(len(target)):
        num += target[i] * out[i]
        denom += out[i] * out[i]
    return num / max(denom, 1e-12)

def kernel_error(out, target, scale=None):
    assert out.shape == target.shape
    tscale = 1.0/numpy.linalg.norm(target.flatten(), 2)    # Normalize 'target' kernel so Euclidean norm is 1
    #print 'tscale', tscale
    if scale is None:
        out = out * optimal_scale(out, target) * tscale           # Scale 'out' kernel by optimal scale factor
    else:
        out = out * scale

    diff = out.flatten() - target.flatten()
    return numpy.sqrt(numpy.dot(diff, diff))                  # Return L2 distance (with sqrt)

