
# PatchTable match stress test.
#
# Usage:
# $ make clean; make debug; python stress.py

import Image
import numpy
import random
import os
import sys; sys.path += ['../patchmatch']
import pfm
import subprocess
from bench import float_last_line

imgcollection_bug = False

min_speed = 0
max_speed = 10

def gen_image():
    h = random.randrange(20, 1100)
    w = random.randrange(20, 1100)
    if random.randrange(2) == 0:
        a = numpy.random.random((h, w, 3))
    else:
        i = random.randrange(133)
        a = Image.open('../patchmatch/vidpairs/a/%03d.png' % i)
        a = a.resize((h, w), Image.ANTIALIAS)
        a = numpy.asarray(a, 'double')/255
    return a
    
def imwrite(a, filename):
    Image.fromarray(numpy.asarray(a*255, 'uint8')).save(filename)

def system(s):
    print s
    os.system(s)

def gen_allowed(b):
    a = numpy.random.random(b[:,:,0].shape) >= 0.5
    return numpy.dstack((a, a, a))

def main():
    def make_prev_nnf():
        prev_nnf_x = numpy.asarray(numpy.random.randint(bew, size=(aeh, aew)), float)
        prev_nnf_y = numpy.asarray(numpy.random.randint(beh, size=(aeh, aew)), float)
        if use_allowed:
            for y in range(aeh):
                for x in range(aew):
                    while not allowed[prev_nnf_y[y,x], prev_nnf_x[y,x], 0]:
                        prev_nnf_y[y,x] = random.randrange(beh)
                        prev_nnf_x[y,x] = random.randrange(bew)
        prev_nnf_dist = numpy.zeros((aeh, aew))
        prev_nnf = numpy.dstack((prev_nnf_x, prev_nnf_y, prev_nnf_dist))
        prev_nnf_filename = 'prev_nnf.pfm'
        pfm.writepfm(prev_nnf, prev_nnf_filename)

        return prev_nnf_filename


    ntest = 100
    patch_w = 8

    if not imgcollection_bug:
        print '=' * 80
        print 'Stress Test Descriptor'
        print '=' * 80
        s = subprocess.check_output('./patchtable test_descriptor', shell=True)
        last_line = s.strip().split('\n')[-1].strip()
        if last_line != 'test_descriptor: OK':
            raise ValueError('descriptor test failed')
        else:
            print 'Passed'

    simple_test_count = 10

    for i in range(ntest):
        random.seed(i)
        numpy.random.seed(i)
        print '=' * 80
        print 'Stress Test %d/%d' % (i, ntest)
        print '=' * 80
        print
        a = gen_image()
        b = gen_image()
        allowed = gen_allowed(b)
        
        aeh = a.shape[0]-patch_w+1
        aew = a.shape[1]-patch_w+1
        beh = b.shape[0]-patch_w+1
        bew = b.shape[1]-patch_w+1
        
        imwrite(a, 'atest.png')
        imwrite(b, 'btest.png')
        imwrite(allowed, 'allowed.png')
        if os.path.exists('out.pfm'):
            os.remove('out.pfm')
        suffix = ''

        use_allowed = False
        r2 = random.randrange(2) == 0
        if i <= simple_test_count:
            r2 = 1
        if r2:
            use_allowed = True

        r = random.randrange(5)
        if i <= simple_test_count:
            r = 3
        speed = random.randrange(min_speed,max_speed)
        if r == 0:
            suffix = ' -speed %d' % speed
        elif r == 1:
            suffix = ' -speed %d -coherence_spatial %f' % (speed, random.uniform(0, 100))
        elif r == 2:
            suffix = ' -speed %d -is_descriptor 1' % speed
        elif r == 3:
            prev_nnf_filename = make_prev_nnf()
            suffix = ' -speed %d -prev_nnf %s -recalc_dist_temporal 1 -coherence_spatial %f -coherence_temporal %f' % (speed, prev_nnf_filename, random.uniform(0, 100), random.uniform(0, 100))

        if use_allowed:
            suffix += ' -allowed_patches allowed.png'

        if imgcollection_bug:
            prev_nnf_filename = make_prev_nnf()
            coherence_spatial = coherence_temporal = 10.0
            suffix = ' -prev_nnf %s -coherence_spatial %f -coherence_temporal %f -limit 1000000 -allowed_patches allowed.png' % (prev_nnf_filename, coherence_spatial, coherence_temporal)

        cmd = './patchtable match atest.png btest.png out.pfm%s' % suffix
        print cmd
        s = subprocess.check_output(cmd, shell=True)
        print s
        mean_dist = float_last_line(s, 'mean_dist:')
        mean_dist_recomputed = float_last_line(s, 'mean_dist_recomputed:')
        if abs(mean_dist-mean_dist_recomputed) > 1e-4:
            raise ValueError('mean_dist (%f) and mean_dist_recomputed (%f) differ' % (mean_dist, mean_dist_recomputed))

        ann = pfm.readpfm('out.pfm')
        xmin = numpy.min(ann[:,:,0])
        xmax = numpy.max(ann[:,:,0])
        ymin = numpy.min(ann[:,:,1])
        ymax = numpy.max(ann[:,:,1])
        assert xmin >= 0 and ymin >= 0, (xmin, ymin)
        assert xmax < bew and ymax < beh, (xmax, ymax, bew, beh)
        
        if use_allowed:
            for y in range(ann.shape[0]):
                for x in range(ann.shape[1]):
                    bx = ann[y,x,0]
                    by = ann[y,x,1]
                    allowed_val = allowed[by,bx,0]
                    if not allowed_val:
                        raise ValueError('accessed disallowed patch: %d, %d => %d, %d, %d' % (x, y, bx, by, allowed_val))
    print
    print 'Unit tests passed.'
    
if __name__ == '__main__':
    main()
    