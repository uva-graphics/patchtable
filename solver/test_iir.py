
import numpy
import random
import scipy.signal
numpy.set_printoptions(linewidth=150, formatter={'all':lambda x: '%.5f'%x})

def height(I):
    return I.shape[0]

def width(I):
    return I.shape[1]

def filter_iir(in_, K, F, xdir=1, ydir=0):
    out = numpy.array(in_, 'double')
    w = width(in_)
    h = height(in_)
    ystart = 0; yend = h; ystep = 1
    xstart = 0; xend = w; xstep = 1
    if ydir < 0:
        ystart = h-1; yend = -1; ystep = -1
    if xdir < 0:
        xstart = w-1; xend = -1; xstep = -1
    hw = len(K)/2
    
    y = ystart
    while y != yend:
        x = xstart
        while x != xend:
            ans = 0.0
#            for dx in range(len(K)):
#                yp = y-dx*ydir; xp = x-dx*xdir
#                if 0 <= yp < h and 0 <= xp < w:
#                    ans += in_[yp, xp] * K[dx]
            for dx in range(-hw, hw+1):
                yp = y + (dx if ydir else 0)
                xp = x + (dx if xdir else 0)
                if 0 <= yp < h and 0 <= xp < w:
                    ans += in_[yp, xp] * K[hw-dx]
            for i in range(len(F)):
                yp = y-(i+1)*ydir; xp = x-(i+1)*xdir
                if 0 <= yp < h and 0 <= xp < w:
                    ans += out[y-(i+1)*ydir,x-(i+1)*xdir] * F[i]
            out[y, x] = ans
            x += xstep
        y += ystep
    return out

def main():
    nK = 1
    nF = 1
    gamma = 3.0
    
    hw = 12
    target = numpy.array([[(1.0/abs(i)**gamma if i != 0 else 1.0) for i in range(-hw, hw+1)]])
    weight = 1.0/target
    weight[0,hw] = 0.0
    
    w = hw*2+1
    I = numpy.zeros((1, w))
    I[0,hw] = 1.0
    
    def fit(v, do_print=False):
        K = v[:nK]
        F = v[nK:]
        out1 = filter_iir(I, K, F)
        out  = filter_iir(out1, K, F, -1)
            
        out_weight = weight * out
        num = numpy.sum(weight * target * out_weight)
        denom = numpy.sum(out_weight * out_weight)
        scale = num/denom

        delta = target - out*scale
        ans = numpy.sum(delta * delta)

        if do_print:
            print 'err:', ans
            print 'target', target
            out_s = out*scale
            print 'out', out_s
            for i in range(len(out_s[0])):
                print 'weight[%d] = %f' % (i, out_s[0, i])
        
        return ans
    
    res = scipy.optimize.minimize(fit, numpy.random.random(nK+nF))
    print 'Error:', res.fun
    print 'Coefficients:'
    print res.x
    #res.x[0] = 1.0
    #fit(res.x, do_print=True)
    
def filter_box(I, boxw):
    box = numpy.ones(boxw)
    unit = numpy.zeros(boxw)
    unit[boxw/2] = 1
    F = numpy.outer(unit, box)
    return scipy.signal.convolve2d(I, F, 'same')
    
def main2():
    K = numpy.array([  0.000000,   0.000000,  1,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   -1,   0.000000])
    F = numpy.array([  1.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000])
    w = 15
    #I = numpy.random.random((w, w))  
    I = numpy.zeros((w, w))
    #I[w/2,w/2] = 1.0
    I[w/2-2:w/2+3,w/2-2:w/2+3] = numpy.random.random((5,5))
    out1 = filter_iir(I, K, F)
    out2 = filter_box(I, 11)
    
    print out1
    print out2
    print numpy.mean(out1.flatten())
    print numpy.mean(out2.flatten())
    d = out1.flatten()-out2.flatten()
    print numpy.linalg.norm(d)
    
if __name__ == '__main__':
    main()
    