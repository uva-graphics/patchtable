
import os
import sys
import glob
import numpy
import shutil
import random
import skimage, skimage.io

def cat_images(L):
#    print len(L)
    if len(L) == 100:
        return numpy.vstack((cat_images(L[:50]), cat_images(L[50:])))
    return numpy.hstack(tuple(L))

def imread(filename):
    return skimage.io.imread(filename)

def main():
    args = sys.argv[1:]
    if len(args) < 2:
        print >> sys.stderr, 'python one_vs_many.py indir outdir [ncat]'
        sys.exit(1)
      
    ncat = 100
    if len(args) >= 3:
        ncat = int(args[2])
    
    (indir, outdir) = args[:2]
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    indir_a = os.path.join(indir, 'a')
    indir_b = os.path.join(indir, 'b')
    outdir_a = os.path.join(outdir, 'a')
    outdir_b = os.path.join(outdir, 'b')
    
    if not os.path.exists(outdir_a):
        os.mkdir(outdir_a)
    if not os.path.exists(outdir_b):
        os.mkdir(outdir_b)
    
    aL = sorted(glob.glob(os.path.join(indir_a, '*.png')))#[:5]
    bL = sorted(glob.glob(os.path.join(indir_b, '*.png')))#[:5]
    #indir_aL = [imread(x) for x in aL]
    indir_bL = [imread(x) for x in bL]
    for i in range(len(indir_bL)):
        print i, indir_bL[i].shape
    h = min([x.shape[0] for x in indir_bL])
    w = min([x.shape[1] for x in indir_bL])
    #print h
    #print w
    #print x.shape
    for i in range(len(indir_bL)):
        #print i, aL[i], indir_aL[i].shape, bL[i], indir_bL[i].shape
        #indir_aL[i] = indir_aL[i][0:h,0:w,:]
        indir_bL[i] = indir_bL[i][0:h,0:w,:]
    
    for i in range(len(aL)):
        filename = os.path.split(aL[i])[1]
        print i, filename
        #shutil.copyfile(aL[i], os.path.join(outdir_a, filename))
        A = imread(aL[i])
        skimage.io.imsave(os.path.join(outdir_a, filename), A)
        remaining_indices = [j for j in range(len(bL)) if j != i]
        b = cat_images([indir_bL[i]] + [indir_bL[j] for j in random.sample(remaining_indices, ncat-1)])
#        remaining_indices = [j for j in range(len(bL)) if j > i]
#        b = cat_images([indir_bL[i]] + [indir_bL[j] for j in remaining_indices[:ncat-1]])
        print bL[i]
        skimage.io.imsave(os.path.join(outdir_b, filename), b)

if __name__ == '__main__':
    main()

