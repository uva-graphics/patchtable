
import sys
import skimage, skimage.io
import sys; sys.path += ['../patchmatch']
import pfm
import numpy
import scipy.io
import os
import os.path
import subprocess

patch_w = 8
use_c = True

def read_mat(filename):
    nnf = scipy.io.loadmat(filename).values()[0]
    nnf = numpy.array(nnf, 'float')
    return nnf

def main():
    args = sys.argv[1:]
    if len(args) < 3:
        print >> sys.stderr, 'check_dist.py a.png b.png out.pfm'
        sys.exit(1)
    
    if use_c:
        pfm_arg = args[2]
        if os.path.splitext(pfm_arg)[1] == '.mat':
            nnf = read_mat(pfm_arg)
            if nnf.shape[2] != 3:
                nnf = numpy.dstack((nnf, numpy.zeros((nnf.shape[0], nnf.shape[1]))))
            pfm_arg = '_temp_dist.pfm'
            pfm.writepfm(nnf, pfm_arg) 
        ans = subprocess.check_output('./patchtable check_dist %s %s %s' % (args[0], args[1], pfm_arg), shell=True)
        sys.stdout.write(ans)
        return
    a = skimage.img_as_float(skimage.io.imread(args[0]))
    b = skimage.img_as_float(skimage.io.imread(args[1]))
    if os.path.splitext(args[2])[1] == '.pfm':
        nnf = pfm.readpfm(args[2])
    else:
        nnf = read_mat(args[2])
        
    beh = b.shape[0]-patch_w+1
    bew = b.shape[1]-patch_w+1
    
    nnf_h = a.shape[0]-patch_w+1
    nnf_w = a.shape[1]-patch_w+1
    dimg = numpy.zeros((nnf_h, nnf_w))
    for ay in range(nnf_h):
        for ax in range(nnf_w):
            bx = nnf[ay,ax,0]
            by = nnf[ay,ax,1]
            assert 0 <= bx < bew, (bx, bew)
            assert 0 <= by < beh, (by, beh)
            apatch = a[ay:ay+patch_w, ax:ax+patch_w, :].flatten()
            bpatch = b[by:by+patch_w, bx:bx+patch_w, :].flatten()
            d = (apatch-bpatch)
            d = numpy.sqrt(numpy.sum(d*d))
            dimg[ay,ax] = d
            
            #d_nnf = numpy.sqrt(nnf[ay,ax,2])
            #delta = abs(d-d_nnf)
            #if delta > 1e-4:
            #    print 'error:', ay, ax, d, dcorrect
    #print '0,0 nnf:', nnf[0,0,:]
    #print '0,0 correct distance:', dimg[0,0]
    #print '%d,0 nnf:'%(nnf_w-1), nnf[0,nnf_w-1,:]
    #print '%d,0 correct distance:'%(nnf_w-1), dimg[0,nnf_w-1]
    print numpy.mean(dimg.flatten())
    
if __name__ == '__main__':
    main()
    
