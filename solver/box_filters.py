
import sys
import numpy
import json
import scipy.signal
import hashlib
import os
from parse_args import parse_args
from scale_kernel import kernel_error
import itertools
numpy.set_printoptions(linewidth=150)

maxw = 15

def system(s):
    #print s
    assert os.system(s) == 0, s
    
def error(box, K):
    #print box
    #K = K/numpy.linalg.norm(K.flatten(), 2)
    #boxf = numpy.zeros(K.shape)
    hw = K.shape[0]/2
    #w2 = w/2.0
    #w2i = int(w2)
    #w2f = w2-w2i
    #horiz[hw-w2i:hw+w2i+1] = 1.0
    #horiz[hw-w2i-1] = w2f
    #horiz[hw+w2i+1] = w2f
    #print horiz
#    boxf[hw-w2i:hw+w2i+1, hw-w2i:hw+w2i+1] = 1.0
    #print boxf

#    ans = numpy.array(boxf)
    for i in range(len(box)):
        w = box[i][0]
        dx = box[i][1]
        horiz = numpy.zeros(K.shape[0])
        horiz[hw-w/2+dx:hw-w/2+w+1+dx] = 1.0
        boxf = numpy.outer(horiz, horiz)

        if i == 0:
            ans = numpy.array(boxf)
        else:
            ans = scipy.signal.convolve2d(ans, boxf, 'same')
    return kernel_error(ans, K)
    #scale = optimal_scale(ans, K)
    #ans *= scale
#    print '-'*80
#    print 'w=%d, nbox=%d'%(w, nbox)
#    print '-'*80
#    print 'ans:'
#    print ans
#    print 'K:'
#    print K
#    print
    

def box_program(box):
    filter_w = maxw
    ans = {'Input': ['ImageParam', 0]}
    for i in range(len(box)):
        idx = i*2
        idx2 = i*2+1
        key_prev = 'Input' if i == 0 else 'Filter%d' % (i*2-1)
        key = 'Filter%d' % idx
        key2 = 'Filter%d' % idx2
        if i == len(box)-1:
            key2 = 'OUT'
        w = box[i][0]
        dx = box[i][1]
        
        K = numpy.zeros(filter_w)
        K[filter_w/2-w/2+dx]     = -1.0
        K[filter_w/2-w/2+w+1+dx] =  1.0
        
        F = numpy.zeros(filter_w)
        F[0] = 1.0
        
        ans[key]  = ['IIRFilter', key_prev, list(K), list(F), 1, 0]
        ans[key2] = ['IIRFilter', key,      list(K), list(F), 0, 1]
    return ans

def usage():
    print >> sys.stderr, 'box_filters kernel.txt [-nmin 1] [-nmax 1] [-pareto filename] [-pareto_full filename]'
    print >> sys.stderr, '  Approximate target kernel using from 1...nmax (square) box filters. Writes Pareto json.'
    sys.exit(1)

def write_pareto(ans, kernel, kw, extra_args=''):
    pareto_filename = kw.get('pareto', 'pareto_summary.txt')
    pareto_full_filename = kw.get('pareto_full', 'pareto_full.txt')
    s = json.dumps(ans)
    filename = hashlib.md5(s).hexdigest()
    temp_filename = '.' + filename + '.txt'
    with open(temp_filename, 'wt') as f:
        f.write(s)
    system('./solve_kernel recalc %s %s -pareto_full %s -pareto %s%s' % (kernel, temp_filename, pareto_full_filename, pareto_filename, ' ' + extra_args if extra_args else ''))
    os.remove(temp_filename)

def box_list(nbox):
    base_list = []
    for w in range(1, maxw+1):
        if w % 2 == 1:
            base_list.append((w, 0))
        else:
            base_list.append((w, 0))
            base_list.append((w, 1))
    L = list(itertools.product(*(base_list,)*nbox))
    
    # Filter Cartesian product for speed
    ans = []
    for box in L:
        sizes = [x[0] for x in box]
        if max(sizes) - min(sizes) <= 2:
            ans.append(box)
    #print nbox, len(L)
    return ans
    
def main():
    (args, kw) = parse_args('nmin nmax pareto pareto_full verbose'.split(), usage)
    if len(args) == 0:
        usage()
        
    nmin = int(kw.get('nmin', '1'))
    nmax = int(kw.get('nmax', '1'))
    verbose = int(kw.get('verbose', '0'))
    
    K = numpy.loadtxt(args[0])
    if len(args) >= 2:
        nmax = int(args[1])

    #maxw = 15
    ans = []
    program_count = 0
    for n in range(nmin, nmax+1):
        boxL = box_list(n)
        #for boxw in range(1, maxw+1, 2):
        #boxw_best = 
        #for boxw in numpy.arange(1, maxw+1, 0.5):
        #    err = error(boxw, n, K)
#            err = error(boxw, n, K)
#            print '%03d %03d %f' % (n, boxw, err)
        #res = scipy.optimize.minimize(lambda boxw: error(boxw, n, K), [5.0])
        err_best = 1e100
        box_best = None
        for box in boxL:
            err = error(box, K)
            if err < err_best:
                err_best = err
                box_best = box
        if verbose:
            print >> sys.stderr, 'n=%d box=%r err=%f' % (n, box_best, err_best)
        ans.append({'name': 'program%03d'%program_count,
                    'time': 0.0,
                    'error': err_best,
                    'program': box_program(box_best)})
        program_count += 1
    write_pareto(ans, args[0], kw, '-recalc_solve 1')
    
if __name__ == '__main__':
    main()
    