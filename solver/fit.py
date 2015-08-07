
import numpy
import scipy.optimize
import pylab
import sys

l2 = True #False
l1_regularize = 0.1
bounds = True #True #True
include_bias = False

def train_test(L_train, L_test, plot=False):
    Amat = numpy.array([[float(x) for x in row[1:-1]] for row in L_train], 'float')
    time = numpy.array([float(x[-1]) for x in L_train], 'float')
    
    if include_bias:
        Amat = numpy.concatenate((Amat, numpy.ones((len(time), 1))), 1)
    nvars = Amat.shape[1]
    Bmat = time
#    print Amat.shape
    #print Bmat.shape
    assert Amat.shape[0] == Bmat.shape[0]
    if l2:
        if bounds:
            def f(x, *args):
                d = numpy.dot(Amat, x) - Bmat
                #if numpy.any(x<0):
                #    return 1e100
                return numpy.dot(d, d)
            #print 'Start min'
#            min_res = scipy.optimize.minimize(f, numpy.array([1]*nvars), method='Nelder-Mead', bounds=[(0.0, None) for i in range(nvars)])
            min_res = scipy.optimize.minimize(f, numpy.array([1.0]*nvars), method='L-BFGS-B', bounds=[(0.0, None) for i in range(nvars)])
            #print 'End min'
            sol = min_res.x
            
        else:
            (sol,resid,rank,s) = numpy.linalg.lstsq(Amat, Bmat)
    else:
        def f(x):
            d = numpy.dot(Amat, x) - Bmat
            #ans = numpy.mean(numpy.absolute(d))
            ans = numpy.dot(d, d)
            if l1_regularize:
                ans += numpy.mean(numpy.abs(x)) * l1_regularize
            return ans
        min_res = scipy.optimize.minimize(f, numpy.array([1]*nvars)) #, method='L-BFGS-B')
        sol = min_res.x
    
    if plot:
        print 'Solution:'
        print sol
#    print
    #print 'Model:', model
#    print 'Solution:', sol
#    if len(sol) == 2:
#        print 'sol[1]/sol[0]:', sol[1]/sol[0]

    Amat_test = numpy.array([[float(x) for x in row[1:-1]] for row in L_test], 'float')
    time_test = numpy.array([float(x[-1]) for x in L_test], 'float')

    if include_bias:
        Amat_test = numpy.concatenate((Amat_test, numpy.ones((len(time_test), 1))), 1)

    time_model_test = numpy.inner(Amat_test, sol)
    #if model == [TAPS]:
    #    time_model = sol[0] * taps
    #elif model == [PASSES]:
    #    time_model = sol[0] * passes
    #elif model == [TAPS, PASSES]:
    #    time_model = sol[0] * taps + sol[1] * passes
    #elif model == [TAPS, PASSES, SURFACEX, SURFACEY]:
    #    time_model = sol[0] * taps + sol[1] * passes + sol[2] * surfacex + sol[3] * surfacey
    #else:
    #    raise ValueError
        
    delta = (time_test - time_model_test)
    
#    print '%-20s %s %s %s' % ('name', 'measured', 'model', 'delta')
#    for i in range(len(time)):
#        print '%-20s %f %f %f' % (name[i], time[i], time_model[i], delta[i])
#    print 'avg delta (Mean):', numpy.mean(numpy.absolute(delta))
#    print 'avg delta (RMS): ', numpy.sqrt(numpy.mean((delta*delta)))
#    print
    if plot:
        print 'Average T_actual in test set:', numpy.mean(time_test)
        print 'Correlation coefficient:', numpy.corrcoef(time_test, time_model_test)[0, 1]
    
    if plot:
        pylab.scatter(time_model_test, time_test)
        pylab.xlabel('time_model')
        pylab.ylabel('time_actual')
        pylab.show()

    return numpy.mean(numpy.absolute(delta))

def main():
    args = sys.argv[1:]
    if len(args) == 0:
        print >> sys.stderr, 'fit.py time.txt'
        sys.exit(1)
    
    L = [x.split() for x in open(args[0], 'rt').read().strip().split('\n')]

    print 'Leave one out:'
    ans = []
    for idx_leave_out in range(len(L)):
        print idx_leave_out, len(L)
        ans.append(train_test([L[i] for i in range(len(L)) if i != idx_leave_out],
                              [L[i] for i in range(len(L)) if i == idx_leave_out]))
    print 'avg delta (Mean):', numpy.mean(ans)
    print
    print 'No cross-validation:'
    print 'avg delta (Mean):', train_test(L, L, True)

if __name__ == '__main__':
    main()
    