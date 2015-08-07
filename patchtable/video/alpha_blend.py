
"""
python alpha_blend.py indir1/file%05d.png indir2/file%05d.png outdir/file%05d.png alpha2
"""

import sys
import skimage
import skimage.io
import os, os.path
import threadmap

def main():
    args = sys.argv[1:]
    if len(args) != 4:
        print >> sys.stderr, __doc__
        sys.exit(1)
    
    (in1, in2, out, alpha) = args
    alpha = float(alpha)
    
    L = []
    for i in range(1, 1000**2):
        if not os.path.exists(in1%i):
            break
        L.append(i)
    
    def process(i):
        print 'Processing frame', i
        I1 = skimage.img_as_float(skimage.io.imread(in1%i))
        I2 = skimage.img_as_float(skimage.io.imread(in2%i))
        I = I1 + (I2-I1) * alpha
        skimage.io.imsave(out%i, I)
    
    threadmap.map(process, L[::-1], dynamic=True)
    
if __name__ == '__main__':
    main()
    