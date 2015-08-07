
import sys
import os
import numpy, numpy.linalg
from parse_args import parse_args
from box_filters import write_pareto

use_1d = False   # Set to True to validate Van Vliet's numbers (need to divide Pareto results by 2.66267, the sum of the filter)
                 # It seems also the L2 error with 0.91e-6 should be 0.91e-3.

def usage():
    print >> sys.stderr, 'vanvliet [-pareto filename] [-pareto_full filename]'
    print >> sys.stderr, 'Generates Pareto that compares with van Vliet et al 1998, Recursive Gaussian Derivative Filters:'
    print >> sys.stderr, 'http://homepage.tudelft.nl/e3q6n/publications/1998/ICPR98LVTYPV/ICPR98LVTYPV.pdf'
    sys.exit(1)

def zpoles_to_b(d):
    b = numpy.product(1.0/d)
    N = len(d)
    ans = [1.0]
    for i in range(1, N+1):
        m = N-i
        lead = b * (-1)**(N-m)
        if i == N:
            ans.append(lead)
        elif i == N-1:
            ans.append(lead * numpy.sum(d))
        elif i == N-2:
            ans.append(lead * sum([sum([d[i]*d[j] for j in range(0, i)]) for i in range(1, N)]))
        elif i == N-3:
            ans.append(lead * sum([sum([sum([d[i]*d[j]*d[k] for k in range(0, j)]) for j in range(1, i)]) for i in range(2, N)]))
        elif i == N-4:
            ans.append(lead * sum([sum([sum([sum([d[i]*d[j]*d[k]*d[l] for l in range(0, k)]) for k in range(1, j)]) for j in range(2, i)]) for i in range(3, N)]))
        else:
            raise ValueError((N, i, m))
    return numpy.array(ans)

def to_real(v):
    assert (numpy.linalg.norm(v.imag) < 1e-8)
    return v.real
    
def van_vliet_program(b):
    K = [1 + numpy.sum(b[1:])]
    F = list(-b[1:])
    if use_1d:
        return {'Input': ['ImageParam', 0],
                'Filter0': ['IIRFilter', 'Input',   K, F,  1,  0],
                'OUT':     ['IIRFilter', 'Filter0', K, F, -1,  0]}
    else:
        return {'Input': ['ImageParam', 0],
                'Filter0': ['IIRFilter', 'Input',   K, F,  1,  0],
                'Filter1': ['IIRFilter', 'Filter0', K, F, -1,  0],
                'Filter2': ['IIRFilter', 'Filter1', K, F,  0,  1],
                'OUT':     ['IIRFilter', 'Filter2', K, F,  0, -1]}

def main():
    (args, kw) = parse_args('pareto pareto_full'.split(), usage)
        
    dL = [numpy.array([1.41650+1.00829j, 1.41650-1.00829j, 1.86543]),
          numpy.array([1.13228+1.28114j, 1.13228-1.28114j, 1.78534+0.46763j, 1.78534-0.46763j]),
          numpy.array([0.86430+1.45389j, 0.86430-1.45389j, 1.61433+0.83134j, 1.61433-0.83134j, 1.87504])]
    bL = [to_real(zpoles_to_b(d)) for d in dL]

    print 'van Vliet first program:'
    print van_vliet_program(bL[0])
    #print bL

    ans = []
    program_count = 0
    for b in bL:
        ans.append({'name': 'program%03d'%program_count,
            'time': 0.0,
            'error': 0.0,
            'program': van_vliet_program(b)})
        program_count += 1
    
    if use_1d:
        K = numpy.loadtxt('gaussian2_large.txt')
        w = K.shape[0]
        K[:w/2,:] = 0.0
        K[w/2+1:,:] = 0.0
        K = K / numpy.linalg.norm(K)
        numpy.savetxt('gaussian2_large1d.txt', K)
    kernel = 'gaussian2_large1d.txt' if use_1d else 'kernels/gaussian2.txt' #'gaussian2_large.txt'
#    kernel = 'gaussian2_large.txt'
    write_pareto(ans, kernel, kw, '-recalc_solve 1 -quantize 0')
    
if __name__ == '__main__':
    main()
    
