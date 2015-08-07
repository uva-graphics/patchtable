
"""
Preconditioned conjugate gradient for Toeplitz systems
"""

import scipy, scipy.linalg, scipy.optimize
import numpy

def chan_precond_vector(a):
    """
    Get Tony Chan 1998 preconditioner -- An Optimal Circulant Preconditioner for Toeplitz Systems. Vector form
    """
    i0 = len(a)/2
    n = i0+1
    print 'i0', i0, 'n', n
    c = [(i*a[i0-(n-i)] + (n-i)*a[i0+i])*1.0/n for i in range(n)]  # range(-(n-1),n)
    print 'c', c
    return c

def chan_precond_matrix(A):
    """
    Get Chan 1998 preconditioner in matrix form
    """
    (h, w) = A.shape
    assert h == w
    diags = [A[h-1-i,0] for i in range(h-1)] + [A[0,i] for i in range(w)]
    c = chan_precond_vector(diags)
    return scipy.linalg.toeplitz(numpy.roll(c[::-1], 1), c)

def precond_func(A):
    Ainv = numpy.linalg.inv(A)
    def f(x):
        return numpy.dot(Ainv, x)
    return f
    
def pcg(A, b, x, precond=lambda q: q, iters=10, is_fr=True, verbose=True):
    """
    Preconditioned conjugate gradient.
    """
    def lin_f(xv):
        d = numpy.dot(A, xv) - b
        return numpy.dot(d, d)
    if verbose:
        print 'cond(A):', numpy.linalg.cond(A, 'fro')
        print 'cond(precond(A)):', numpy.linalg.cond(precond(A), 'fro')
    r = b - numpy.dot(A, x)
    z = precond(r)
    p = numpy.array(z)
    for k in range(iters):
        beta_denom = numpy.dot(z, r)
        Ap = numpy.dot(A, p)
        alpha = numpy.dot(r, z) / numpy.dot(p, Ap)
        x = x + alpha * p
        rprev = r
        r = r - alpha * Ap
        # if rk is sufficiently small then exit loop
        z = precond(r)
        if is_fr:        
            beta = numpy.dot(z, r) / beta_denom
        else:
            beta = numpy.dot(z, (r-rprev)) / beta_denom
        p = z + beta * p
        if verbose:
            print 'x(%d) ='%k, x, 'r = ', r, 'f = ', lin_f(x)
    return x

def test_chan():
    a0  = numpy.random.random()
    a1  = numpy.random.random()
    a2  = numpy.random.random()
    an1 = numpy.random.random()
    an2 = numpy.random.random()
    n = 3
    
    c0 = ((n)*a0) * 1.0/n
    c1 = (an2 + (n-1)*a1) * 1.0/n
    c2 = (2*an1 + (n-2)*a2) * 1.0/n
    
    A = numpy.array([[a0, a1, a2], [an1, a0, a1], [an2, an1, a0]])
    C = numpy.array([[c0, c1, c2], [c2, c0, c1], [c1, c2, c0]])    
    
    c0 = a0
    c1 = (2*a1 + an2) / 3.0
    c2 = (a2 + 2*an1) / 3.0
    
    C2 = numpy.array([[c0, c1, c2], [c2, c0, c1], [c1, c2, c0]])    
    
    print 'A:'
    print A
    print
    print 'C:'
    print C
    print
    print 'C(chan):'
    print chan_precond_matrix(C)    
    print
    print 'C2:'
    print C2

def solve_chan(A):
    def f(c):
        C = scipy.linalg.circulant(c)
        return numpy.linalg.norm(A-C, 'fro')
    n = A.shape[0]
    res = scipy.optimize.minimize(f, numpy.zeros(n))
    return scipy.linalg.circulant(res.x)

def pcg_toeplitz(A, b, x):
    C = chan_precond_matrix(A)
    return pcg(A, b, x, precond_func(C))
    
def main():
    numpy.random.seed(0)
    n = 4
    A = scipy.linalg.toeplitz(numpy.random.random(n), numpy.random.random(n))
    print 'A:'
    print A
    print
    print 'eigenvalues of A:'
    print numpy.linalg.eigvals(A)
    print
    #print 'C:'
    #print chan_preconditioner_matrix(A)
    b = numpy.random.random(n)
    print 'b:'
    print b
    print
    print 'Solution (linalg):'
    print numpy.linalg.solve(A, b)
    print
    print 'Solution (cg):'
    print pcg(A, b, numpy.zeros(n))
    print
    C = chan_precond_matrix(A)
    print 'chan:'
    print C
    print
    spec_A = numpy.dot(numpy.linalg.inv(C), (A))
    print 'eigenvalues:'
    print numpy.abs(numpy.linalg.eigvals(spec_A))
    print
    precond = precond_func(C)
    print
    print 'norm_fro: %f' % (numpy.linalg.norm(A - C, 'fro'))
    print 'solve_chan:'
    print solve_chan(A)
    print
    print 'Solution (pcg):'
    print pcg(A, b, numpy.zeros(n), precond)
    
if __name__ == '__main__':
    main()
    #test_chan()
    