
"""
Finds Gaussian sigmas for shift-invariant downsampling/upsampling.

Results:

xmax=1
0.9103125
0.0320246446732
max:   0.0688866249307
mean:  0.0514008558052

xmax=2
1.1625
0.00826367983468
max:   0.0100686787493
mean:  0.00741959259365

xmax=3
1.335
0.00193814926287
max:   0.00427119443061
mean:  0.00190132326245

"""
import json
import numpy
import sys
import check_alias
import scipy.optimize

xmax=2
is_gaussian = True
use_uniform = False

def gaussian(sigma=1.0, sz=None):
    if sz is None:
        sz = xmax
    X = numpy.arange(-float(sz), float(sz)+1)
    Gh1 = numpy.exp(-X**2/(2*sigma**2))
    K1 = numpy.outer(Gh1, Gh1)
    K1 = K1 / numpy.sum(K1.flatten())
#    print >> sys.stderr, Gh1
    return K1

def program(sigma_vec, seed=0):
    if is_gaussian:
        sigma1 = sigma_vec[0]
        sigma2 = sigma_vec[0]
        K1 = gaussian(sigma=sigma1)
        K2 = gaussian(sigma=sigma2)
    else:
        K1 = numpy.concatenate((sigma_vec, [1.0], sigma_vec[::-1]))
        K1 = numpy.outer(K1, K1)
        K1 = K1 / numpy.sum(K1.flatten())
        K2 = K1
    numpy.random.seed(seed)
    if use_uniform:
        arb_K1 = 1*(numpy.random.random((3, 3))*2-1)
        arb_K2 = 1*(numpy.random.random((5, 5))*2-1)
    else:
        arb_K1 = gaussian(numpy.random.uniform(0.5, 1.0), 2) + numpy.random.random((5, 5))*0.1
        arb_K2 = gaussian(numpy.random.uniform(0.5, 1.0), 2) + numpy.random.random((5, 5))*0.1
    return {'Input': ['ImageParam', 0],
            'Filter0': ['FIRFilter', 'Input', K1.tolist()],
            'Filter1': ['FilterDownsample', 'Filter0'],
            'Filter2': ['FIRFilter', 'Filter1', arb_K1.tolist()],
            'Filter3': ['FilterUpsample', 'Filter2'],
            'Filter4': ['FIRFilter', 'Filter3', K2.tolist()],
            'Filter5': ['FIRFilter', 'Input', arb_K2.tolist()],
            'OUT': ['FilterAdd', 'Filter4', 'Filter5']}

def gaussian_band_error(sigma_vec, seed_start=0, seed_count=50, is_max=True, print_result=False):
    ans = []
    for seed in range(seed_start, seed_start+seed_count):
        filename = 'gaussian_band.txt'
        L = [{'name': 'program000', 'time': 0.0, 'error': 0.0, 'scale': 0.0, 'program': program(sigma_vec, seed)}]
        s = json.dumps(L)
        with open(filename, 'wt') as f:
            f.write(s)
        ans.append(check_alias.check_alias(None, filename, 0, verbose=False))
    if print_result:
        print 'max:  ', numpy.max(ans)
        print 'mean: ', numpy.mean(ans)
    if is_max:
        return numpy.max(ans)
    else:
        return numpy.sum(numpy.array(ans)**2)
    
def main():
    #args = sys.argv[1:]
    global xmax
    #if len(args) < 1:
    #    print >> sys.stderr, 'gaussian_band xmax'
    #    sys.exit(1)
    #xmax = int(args[0])
    #x0 = [ 1.21 ]   # xmax=1, err=0.11656882282
    #x0 = [0.24966785127548971, 0.70687179946533907]
    x0 = [1.2]
    #gaussian_band_error(x0, seed_start=13040, seed_count=20, is_max=True, print_result=True)
    #x0 = [ 1.113, 1.479]   # xmax=2, err=0.0287925412356
    #x0 = [ 1.251, 1.657]   # xmax=3, err=0.00582298367526, max=0.00694371600185
#    x0 = [ 1.422, 1.743]   # xmax=4, err=0.00164909720092
#    x0 = [ 1.568, 1.87 ]   # xmax=5, err=0.000527602634343
    
    # Gaussian coefficients:
    #[ 0.057  0.279  0.727  1.     0.727  0.279  0.057]
    #[ 0.195  0.484  0.834  1.     0.834  0.484  0.195]

    for xmax in range(1, 6):
        #errL = []
        #print gaussian_band_error([ 1.421, 1.748])
        #return
        res = scipy.optimize.fmin(gaussian_band_error, x0) #, seed_count=10)
        print 'xmax=%d'%xmax
        print res[0]
        print gaussian_band_error(res)
        
        #for seed in range(100):
        #    errL.append(gaussian_band_error(x0, seed=seed))
        #    #print errL[-1]
        gaussian_band_error(x0, seed_start=13040, seed_count=100, is_max=True, print_result=True)#'max:', numpy.max(errL)

if __name__ == '__main__':
    main()
