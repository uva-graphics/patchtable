
import sys
import numpy
import os
from parse_args import parse_args
from scale_kernel import *
import json
numpy.set_printoptions(linewidth=150, precision=3)
import pylab

DEFAULT_LEVELS = 1
DEFAULT_OUTSPEC = 'filter_out%02d.txt' 
DEFAULT_FACTOR = 1
DEFAULT_W = 30
move_target = False

def usage():
    print >> sys.stderr, 'python check_alias.py kernel.txt pareto_full.txt i [-target_w w] [-levels n=1] [-outspec "filter_out%02d.txt"]'
    print >> sys.stderr, '  [-factor 1|2]'
    print >> sys.stderr, '  Retrieves filter i from the Pareto Frontier. Can pass - to use a kernel produced from shift (0, 0).'
    print >> sys.stderr, '  Computes 4**levels responses to a shifted unit impulse function, writing with numpy.savetxt to outspec.'
    print >> sys.stderr, '  Prints approximation errors of each shift against kernel.'
    sys.exit(1)

def pad_center(A, (h, w)):
    ans = numpy.zeros((h, w))
    for y in range(h):
        for x in range(w):
            ysrc = y-h/2+A.shape[0]/2
            xsrc = x-w/2+A.shape[1]/2
            if 0 <= ysrc < A.shape[0] and 0 <= xsrc < A.shape[1]:
                ans[y,x] = A[ysrc, xsrc]
    return ans
    
def system(s):
    #print s
    #print os.system(s)
    os.system(s)
    
def check_alias(kernel, pareto, index, w=DEFAULT_W, levels=DEFAULT_LEVELS, factor=DEFAULT_FACTOR, outspec=DEFAULT_OUTSPEC, verbose=True):
    max_shift = 2**levels
    current = 0
    all_error = []
    
    scaleL = []
    outL = []
    targetL = []
    deltaL = range(max_shift*factor) #range(-max_shift+1, 1)
    for dy in deltaL:#range(max_shift*factor):
        for dx in deltaL: #range(max_shift*factor):
            unit = numpy.zeros((w, w))
            unit[w/2+dy,w/2+dx] = 1
            unit_filename = '.unit.txt'
            numpy.savetxt(unit_filename, unit)
            system('./solve_kernel apply %s %s %d %s -target_w %d' % (unit_filename, pareto, index, outspec%current, w))
            
            out = numpy.loadtxt(outspec%current)
            #print w, unit.shape, out.shape
            if kernel is None:
                kernel = numpy.array(out)
            if move_target:
                target = numpy.zeros((w, w))
                for y in range(w):
                    for x in range(w):
                        ysrc = y-dy
                        xsrc = x-dx
                        if 0 <= ysrc < w and 0 <= xsrc < w:
                            target[y,x] = kernel[ysrc,xsrc]
            else:
                target = numpy.array(kernel)
                for y in range(w):
                    for x in range(w):
                        ysrc = y+dy
                        xsrc = x+dx
                        if 0 <= ysrc < w and 0 <= xsrc < w:
                            out[y,x] = out[ysrc,xsrc]
                        else:
                            out[y,x] = 0.0
                #print out
            outL.append(out)
            targetL.append(target/numpy.linalg.norm(target.flatten(), 2))
            #scaleL.append(optimal_scale(out, target))
            current += 1
    outL_flatten = numpy.concatenate([x.flatten() for x in outL[:1]])
    targetL_flatten = numpy.concatenate([x.flatten() for x in targetL[:1]])
    #tscale = 1.0 / numpy.linalg.norm(targetL_flatten, 2)
    #print 'tscale0', tscale
    #targetL_flatten = targetL_flatten * tscale# * len(deltaL) * len(deltaL)
    #for i in range(len(targetL)):
    #    targetL[i] *= tscale
    
    scale = optimal_scale(outL_flatten, targetL_flatten)
#    print scale
    #print kernel_error(outL_flatten, targetL_flatten, scale)
    #scale = numpy.mean(scaleL)
#    print scaleL
    #print 'scale:', scale
    #print outL[0]
    #print targetL[0]
    #print scaleL[0]
    #pylab.imshow(targetL[0])
    #pylab.show()
    #pylab.imshow(outL[0])
    #pylab.show()
    #print 'norm:', numpy.linalg.norm(outL[0]-targetL[0])
    #print 'kernel_error:', kernel_error(outL[0], targetL[0], scale)
    
    current = 0
    for dy in deltaL: #range(max_shift*factor):
        for dx in deltaL: #range(max_shift*factor):
            current_error = kernel_error(outL[current], targetL[current], scale)
            all_error.append(current_error)
            if verbose:
                print '%-30s %.5f' % (outspec%current, current_error)
            current += 1
        if verbose:
            print
    ans = numpy.mean(numpy.array(all_error)**2)**0.5
    if verbose:
        print ans
    return ans

def main():
    (args, kw) = parse_args('target_w levels outspec factor'.split(), usage)
    if len(args) < 3:
        usage()

    kernel_filename = args[0]
    if kernel_filename == '-':
        kernel = None
    else:
        kernel = numpy.loadtxt(kernel_filename)
    pareto = args[1]
    paretoL = json.loads(open(pareto, 'rt').read())
    index = int(args[2])
    
    if index < 0:
        index = len(paretoL) + index
    
    w = int(kw.get('target_w', kernel.shape[0] if kernel is not None else DEFAULT_W))
    levels = int(kw.get('levels', DEFAULT_LEVELS))
    factor = int(kw.get('factor', DEFAULT_FACTOR))
    outspec = kw.get('outspec', DEFAULT_OUTSPEC)
    if kernel is not None and w != kernel.shape[0]:
        kernel = pad_center(kernel, (w, w))
        numpy.savetxt('filter_kernel.txt', kernel)
    check_alias(kernel, pareto, index, w, levels, factor, outspec)
    
if __name__ == '__main__':
    main()

