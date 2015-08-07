
import numpy
import pylab
#import scipy.signal
from scipy.ndimage.filters import convolve as convolve_orig
import scipy.optimize, scipy.signal
import cg
# scipy.signal.convolve2d

def convolve(A, B):
    return convolve_orig(A, B, mode='nearest')

def convolve_full(A, B):
    return scipy.signal.convolve2d(A, B)

def promote(A, B):
    npad = ((B.shape[0]-A.shape[0])/2, (B.shape[1]-A.shape[1])/2)
    return numpy.pad(A, npad, 'constant')
    
numpy.set_printoptions(threshold=numpy.nan, linewidth=1000, suppress=True)

def truncate_kernel(K, maxsize):
    assert maxsize % 2 == 1
    
    h = K.shape[0]
    w = K.shape[1]
    hsel = min(h, maxsize)
    wsel = min(w, maxsize)
    top = h/2-hsel/2
    left = w/2-wsel/2
    return K[top:top+hsel, left:left+wsel]    

def print_matrix(desc, A):
    A = numpy.array(A)
    A[numpy.abs(A)<1e-10] = 0
    Af = A.flatten()
    num_nonzero = sum(Af!=0)
    print desc
    print '%dx%d' % (A.shape[0], A.shape[1])
    print A
    print '  (Nonzero: %.2f%%)' % (100.0*num_nonzero/len(Af))
    print

def solve_kernel(A, K, A_size=None, B_size=None, Apad='nearest', plot=False, use_pcg=False, verbose=True):
    """
    Solve for optimal kernel B s.t. A[0] * A[1] * ... * B = K (* = convolution, in least squares sense).
    
    Apad options are 'nearest', 'constant', methods for looking up missing values in A.
    """
    Aconv = A[0]
    for Asub in A[1:]:
        Aconv = convolve(Aconv, Asub)
    if A_size is not None:
        Aconv = truncate_kernel(Aconv, A_size)
        
    # Aconv * B = K (convolution),   convert to matrix form    Amat   *   B.flatten()  =   K.flatten()
    #                                                        nK * nB        nB x 1           nK x 1
    if B_size is None:
        B_size = (K.shape[0], K.shape[1])
        
    Kh = K.shape[0]
    Kw = K.shape[1]
    Ah = Aconv.shape[0]
    Aw = Aconv.shape[1]

    nK = Kw*Kh
    nB = B_size[0]*B_size[1]
    #Amat = scipy.sparse.lil_matrix((nK, nB))
    Amat = numpy.zeros((nK, nB))
    #Kf = K.flatten()
    Bshift = numpy.array(B_size)/2#-7
    
    for iK in range(Kh):          # K row (y)
        for jK in range(Kw):      # K column (x)
            idxK = iK*Kw + jK
#            assert K.flatten()[idxK] == K[iK, jK]
            
            # Find all coords in A, B that add up to (iK, jK)
            if Apad == 'nearest':
                iAlist = range(-Ah, 2*Ah)
                jAlist = range(-Aw, 2*Aw)
            elif Apad == 'constant':
                iAlist = range(Ah)
                jAlist = range(Aw)
            else:
                raise ValueError(Apad)
                
            for iA in iAlist:
                iB = iK - iA + Bshift[0]
                if 0 <= iB < B_size[0]:
                    for jA in jAlist:
                        jB = jK - jA + Bshift[1]
                        if 0 <= jB < B_size[1]:
                            idxB = iB*B_size[1] + jB
                            #idxA = iA*Aw + jA
                            Amat[idxK, idxB] = Aconv[numpy.clip(iA, 0, Ah-1), numpy.clip(jA, 0, Aw-1)]
    if verbose:
        print_matrix('Amat:', Amat)
    if plot:
        pylab.imshow(Amat)
        pylab.show()

    Kf = K.flatten()
    
    if use_pcg:
        B = cg.pcg_toeplitz(Amat, Kf, numpy.zeros(len(Kf)))
    else:
        B = numpy.linalg.lstsq(Amat, Kf)[0]
    return (B.reshape(B_size), Amat)
    
def test_outer():
    sub = [0.1, 0.3, 1.0, 0.4, 2.0]
    nA = len(sub)
    A = numpy.zeros((nA, nA))
    xA = numpy.arange(nA)
    sigma = 1
    #A[nA/2,xA] = numpy.exp(-(xA-nA/2)**2.0/(2*sigma**2))
    A[nA/2,:] = sub
    print_matrix('A:', A)
    #A = numpy.random.random((nA,nA))
    
    nK = nA
    yA = numpy.arange(nA)
    [x, y] = numpy.meshgrid(xA, yA)

    #K = numpy.zeros((nK, nK))
    K = numpy.outer(sub, sub)
    #K = numpy.exp(-((x-nA/2)**2.0+(y-nA/2)**2.0)/(2.0*sigma**2))
    print_matrix('K:', K)

    (B, Amat) = solve_kernel([A], K)
    print_matrix('B (solution):', B)

def A_box(nA=15):
    A = numpy.zeros((nA, nA))
    A[:nA/2+1,:nA/2+1] = 1
    print_matrix('A:', A)
    return A

def test_box(boxw=3):
    hw = boxw/2
    nA = 17

    A = A_box(nA)
    
    K = numpy.zeros((nA, nA))
    K[nA/2-hw:nA/2+hw+1,nA/2-hw:nA/2+hw+1] = 1
    print_matrix('K:', K)
    
    (B, Amat) = solve_kernel([A], K)#, B_size=3)
    print_matrix('B (solution):', B)

def test_box1d(boxw=5):
    hw = boxw/2
    nA = 17
    A = numpy.zeros((1,nA))
    A[:,:nA/2+1] = 1
    print_matrix('A:', A)
    
    K = numpy.zeros((1, nA))
    K[:,nA/2-hw:nA/2+hw+1] = 1
    #K[:,nA/2+1] = 2
    print_matrix('K:', K)
    
    (B, Amat) = solve_kernel([A], K) #, Apad='nearest')#, B_size=3)
    print_matrix('B (solution):', B)
    
    prod = convolve(A, B)
    print_matrix('A * B:', prod)
    
    #B0 = numpy.array(B)
    #B0[:,:3] = 0
    #print_matrix('B0:', B0)
    #prod = convolve(A, B0)
    #print_matrix('A * B0:', prod)
    
    return (A, K, B)

def test_rand():
    numpy.random.seed(0)
    nA = 5 #25 #17

    A = numpy.random.random((nA, nA))
    
    K = numpy.random.random((nA, nA))
    print_matrix('K:', K)
    
    (B, Amat) = solve_kernel([A], K, plot=True)#, B_size=3)
    print_matrix('B (solution):', B)

convolve_func = convolve
def resid((B0, B1), K, convolve=False, full=True):
    if full:
        B = convolve_full(B0, B1)
    else:
        B = convolve_func(B0, B1)

    Bf = B.flatten()
    Kf = K.flatten()
    beta = numpy.sum(Kf) / numpy.sum(Bf)
    if convolve:
        return beta*B
    d = beta * B.flatten()-K.flatten()
    return numpy.dot(d, d)

def split_arg(x, B0_orig):
    return (numpy.array(x[:len(x)/2]).reshape(B0_orig.shape),
            numpy.array(x[len(x)/2:]).reshape(B0_orig.shape))

def optimize_func(B0_orig, B1_orig, K, f, x0=None, Bmul=None):
    if Bmul is None:
        B0_mul = numpy.ones_like(B0_orig)
        B1_mul = numpy.ones_like(B1_orig)
    else:
        (B0_mul, B1_mul) = split_arg(Bmul, B0_orig)
        
    if x0 is None:
        x0 = numpy.array(list(B0_orig.flatten()) + list(B1_orig.flatten()))
    res = scipy.optimize.minimize(f, x0)
    (B0, B1) = split_arg(res.x, B0_orig)
    print_matrix('B0:', B0*B0_mul)
    print_matrix('B1:', B1*B1_mul)
    print_matrix('B0 * B1:', resid((B0*B0_mul, B1*B1_mul), K, convolve=True))
    print_matrix('K:', K)
    print 'resid:', resid((B0*B0_mul, B1*B1_mul), K)
    return res.x

def solve_simul((B0_orig, B1_orig), K):
    def f(x):
        (B0, B1) = split_arg(x, B0_orig)
        return resid((B0, B1), K)
    return optimize_func(B0_orig, B1_orig, K, f)

def test_split_l1((B0_orig, B1_orig), K):
    C0 = convolve_full(B0_orig, B1_orig)
    K = promote(K, C0)
    
    gamma = 0.0
    def f(x):
        (B0, B1) = split_arg(x, B0_orig)
        B0_ss = numpy.sqrt(numpy.sum(B0*B0))
        B1_ss = numpy.sqrt(numpy.sum(B1*B1))
        B0_norm = B0 / B0_ss
        B1_norm = B1 / B1_ss
        C = convolve_full(B0_norm, B1_norm)
        Cf = C.flatten()
        Kf = K.flatten()
        beta = numpy.sum(Kf) / numpy.sum(Cf)

        resid1 = beta*Cf - Kf
        err1 = numpy.dot(resid1, resid1)

        B0_sparsity = numpy.sum(numpy.abs(B0)) / B0_ss
        B1_sparsity = numpy.sum(numpy.abs(B1)) / B1_ss
        err2 = gamma * (B0_sparsity + B1_sparsity)
        
        ans = err1 + err2
        #print 'ans', ans
        #print 'err1', err1
        #print 'err2', err2
        return ans
    x0 = None
    gamma = 0.0
    while gamma <= 0.1+1e-8:
        x0 = optimize_func(B0_orig, B1_orig, K, f, x0)
        gamma += 0.01

def test_split_l0((B0_orig, B1_orig), K, proper=False):
    C0 = convolve_full(B0_orig, B1_orig)
    K = promote(K, C0)
    Bmul = numpy.ones(len(B0_orig.flatten()) + len(B1_orig.flatten()))

    def f(x):
        (B0, B1) = split_arg(x, B0_orig)
        (B0_mul, B1_mul) = split_arg(Bmul, B0_orig)
        C = convolve_full(B0*B0_mul, B1*B1_mul)

        resid = (C - K).flatten()
        err1 = numpy.dot(resid, resid)

        return err1
    x0 = None
    
    i_done = set()
    while sum(Bmul.flatten()) > 0:
        x0 = optimize_func(B0_orig, B1_orig, K, f, x0, Bmul)
        (B0, B1) = split_arg(x0, B0_orig)
        ichosen = -1
        if not proper:
            B0_scale = 1.0/numpy.sum(B0.flatten())
            B1_scale = 1.0/numpy.sum(B1.flatten())
            x0 = list(B0.flatten()*B0_scale) + list(B1.flatten()*B1_scale)
            for i in range(len(x0)):
                if i not in i_done:
                    if ichosen == -1 or abs(x0[i]) < abs(x0[ichosen]):
                        ichosen = i
        else:
            #f0 = f(x0)
            fmin = 1e100
            for i in range(len(x0)):
                if i not in i_done:
                    orig = Bmul[i]
                    assert orig > 0
                    Bmul[i] = 0
                    fp = f(x0)
                    if fp < fmin:
                        fmin = fp
                        ichosen = i
                    Bmul[i] = orig
        i = ichosen
        i_done.add(i)
        print 'ichosen', ichosen
        print Bmul
        assert Bmul[ichosen] > 0
        Bmul[ichosen] = 0

def random_sym((m, n)):
    ans = numpy.random.random((m, n))
    ans = (ans + numpy.flipud(ans)) / 2
    ans = (ans + numpy.fliplr(ans)) / 2
    return ans
    
def test_split(K, iters=30, simul=True):
    numpy.random.seed(0)
    #B = [numpy.random.random((K.shape[0], K.shape[1])) for i in range(2)]
    #B = [numpy.ones((K.shape[0], K.shape[1])) for i in range(2)]
    B = [random_sym((K.shape[0], K.shape[1])) for i in range(2)]
    #B = [numpy.array(K) for i in range(2)]
    #B[1][:,:2] = 0
    #B[1][:,3:] = 0
    if simul:
        #solve_simul((B[0], B[1]), K)
        test_split_l0((B[0], B[1]), K)        
        return
    
    for j in range(iters):
        for k in range(len(B)):
            print_matrix('B[%d] [iter %d]:'%(k,j), B[k])
            (Bp, Amat) = solve_kernel([B[i] for i in range(len(B)) if i != k], K, verbose=False, Apad='constant')
            B[k] = Bp
    print_matrix('B0 * B1:', convolve(B[0], B[1]))
    print_matrix('K:', K)
    print 'resid:', resid((B[0], B[1]), K, full=False)
    
def test_split_gaussian():
    sigma = 1.3
    n = 2
    x = numpy.arange(-n,n+1)
    y = numpy.arange(-n,n+1)
    X, Y = numpy.meshgrid(x, y)
    K = numpy.exp(-(X**2 + Y**2) / (2.0*sigma**2))
    K = K / numpy.sum(K.flatten())
    print_matrix('K', K)
    test_split(K)    
    
if __name__ == '__main__':
    test_split_gaussian()
    #test_rand()
    #test_box()
    #test_box1d()
    