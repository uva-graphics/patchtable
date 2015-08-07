
import numpy
import scipy.optimize

l2 = True

K = numpy.array([1.0, 2.0, 3.0, 2.0, 1.0])
numpy.random.seed(0)
A = numpy.random.random((3,)) #numpy.array([1.0, 2.0, 1.0])
B = numpy.random.random((3,)) #numpy.array([1.0, 1.0, 1.0])

def f(x, disp=False):
#    res = numpy.convolve(x, B)
#    res = numpy.convolve(A, x)
    res = numpy.convolve(numpy.convolve(A, x), B)
    if disp:
        print
        print 'Target Kernel:', K
        print 'Current Kernel:', res
#    print res.shape
#    print res
    ndiff = (len(res)-len(K))/2
    Kpad = numpy.pad(K, (ndiff,), 'constant')
    resid = Kpad-res
    if l2:
        err = numpy.sum(resid*resid)
    else:
        err = numpy.mean(numpy.absolute(resid))
    if disp:
        print 'Error:', err
    return err
    #
    #return numpy.mean(numpy.absolute(d))

x0 = numpy.ones(15)
result = scipy.optimize.minimize(f, x0, options={'maxiter': 1000})#, method='TNC')#, method='SLSQP')
sol = result.x
print 'Solution:', sol
print 'Objective:', result.fun
f(sol, True)
#sol[0] = 0.0
#f(sol, True)

#print f(x0)