
import sys
import glob
import subprocess
import numpy
import os, os.path
from numpy.distutils.lib2def import DEFAULT_NM

include_precompute = False #True     # Include precomputation time in result
different_categories = True # True
subtract_dt = False             # Remove the dt time from the result 

DEFAULT_N = 4
#rootdirL = ['../patchmatch/vidpairs', '../patchmatch/vidpairs_small']
rootdirL = ['../patchmatch/vidpairs']

def usage():
    print >> sys.stderr, './bench ours|pm|csh|treecann nimages|. "switches to patchdb" [options]'
    print >> sys.stderr, 'Options: [-one_vs_many 0|1] [-lores 0|1] [-different_categories 0|1] [-skip m]'
    print >> sys.stderr, 'e.g. ./bench ours 5 "-ndims 6 -limit 340 -do_rs 1 -randomize_dt 1"'
    print >> sys.stderr, 'Prints mean_time, mean_error, switches (if pm or csh, writes bench_pm.m|bench_csh.m which does this)'
    print >> sys.stderr, 'Here . uses default n'
    sys.exit(1)

def image_pairs(n, skip=0):
    def add_pairs(aL, bL):
        assert len(aL) == len(bL)
        if n < len(aL):
            aL = aL[:n]
            bL = bL[:n]
        for i in range(len(aL)):
            ans.append((aL[i], bL[i]))

    ans = []
    for rootdir in rootdirL:
        aL = sorted(glob.glob(os.path.join(rootdir, 'a', '*.png')))
        bL = sorted(glob.glob(os.path.join(rootdir, 'b', '*.png')))
        add_pairs(aL, bL)
        
        if different_categories:
            mid = len(bL)/2
            bL_rot = bL[mid:] + bL[:mid]
            add_pairs(aL, bL_rot)
    return ans[skip:]
        
def float_last_line(s, msg):
    found = False
    ans = 0.0
    L = s.split('\n')
    for line in L:
        if line.startswith(msg):
            ans = float(line[len(msg):].strip().split()[0])
            found = True
    assert found
    return ans
    
def main():
    args = sys.argv[1:]
    if len(args) < 3:
        usage()
    method = args[0]
    switches = args[2]
    
    verbose = True
    timeL = []
    distL = []
    if method in ['pm', 'csh', 'treecann']:
        mfile = open('bench_%s.m'%method, 'wt')
        print >> mfile, 'maxNumCompThreads(1);'
        print >> mfile, 'patch_w = 8;'
        print >> mfile, 'anames = cell(0,0); bnames = cell(0, 0);'
    
    global rootdirL
    global different_categories
    skip = 0
    if '-different_categories' in args:
        different_categories = int(args[args.index('-different_categories')+1])
    if '-lores' in args and int(args[args.index('-lores')+1]) == 1:
        rootdirL = ['../patchmatch/vidpairs_small']
    if '-one_vs_many' in args and int(args[args.index('-one_vs_many')+1]) == 1:
        global DEFAULT_N
        DEFAULT_N = 1
        rootdirL = ['../patchmatch/one_vs_many']
        different_categories = False
    if '-skip' in args:
        skip = int(args[args.index('-skip')+1])

    if args[1] == '.':
        n = DEFAULT_N
    else:
        n = int(args[1])
    for (a, b) in image_pairs(n, skip):
        if method == 'ours':
            cmd = './patchtable match %s %s out.pfm %s' % (a, b, switches)
            if verbose:
                print >> sys.stderr, cmd
            s = subprocess.check_output(cmd, shell=True)
            if verbose:
                print >> sys.stderr, s
    #        timeL.append(float_last_line(s, 'k-NN overall lookup time:'))
#            lookup_T = float_last_line(s, '    1-NN overall lookup time (excluding exact patch dist computation):')
            lookup_T = float_last_line(s, '    1-NN overall lookup time:')
            precompute_T = float_last_line(s, 'table total precomputation time:')
            T = (lookup_T + precompute_T) if include_precompute else lookup_T
            if subtract_dt:
                dt_T = float_last_line(s, 'table dt time (prop):')
                T -= dt_T
            timeL.append(T)
            distL.append(float_last_line(s, 'mean_dist:'))
        elif method in ['pm', 'csh', 'treecann']:
            print >> mfile, "anames{length(anames)+1} = '%s';" % a
            print >> mfile, "bnames{length(bnames)+1} = '%s';" % b
        else:
            raise ValueError
    if method in ['pm', 'csh', 'treecann']:
        print >> mfile, 'for i=1:length(anames)'
        prefix = suffix = ''
        if method == 'pm':
            prefix = 'im2double('
            suffix = ')'
        
        print >> mfile, '    aimg{i} = %(prefix)simread(anames{i})%(suffix)s;' % locals()
        print >> mfile, '    bimg{i} = %(prefix)simread(bnames{i})%(suffix)s;' % locals()
        print >> mfile, 'end'
        if method == 'pm':
            print >> mfile, "addpath 'patchmatch-2.1';"
        elif method == 'csh':
            print >> mfile, "addpath 'CSH_code_v2';"
            print >> mfile, 'cd CSH_code_v2; AddPaths; cd ..;'
        elif method == 'treecann':
            print >> mfile, "addpath 'TreeCANN';"
            print >> mfile, "addpath 'TreeCANN/C_code';"
            print >> mfile, "addpath 'TreeCANN/matlab_tools/ann_wrapper/';"
        else:
            raise ValueError
        iters = '1 2 3 4 5 6 7 8 9 10'
        threadL = [1]
        if method == 'pm':
            threadL = [1] #[8, 1]
            iters += ' 50'
        elif method == 'treecann':
            iters = '10 9 8 7 6 5 4 3 2 1'
            
        print >> mfile, 'for threads=%(threadL)r' % locals()
        print >> mfile, '    for iters=[%(iters)s]' % locals()
        if method == 'treecann':
            print >> mfile, '    for bgrid=[%(iters)s]' % locals()
        print >> mfile, '        dL = []; tL = [];'
        print >> mfile, '        for i=1:length(anames)'
        print >> mfile, '            a = aimg{i};'
        print >> mfile, '            b = bimg{i};'
        py = subprocess.check_output('which python', shell=True).strip()
        if method == 'pm':
            print >> mfile, "            tic; nnf=nnmex(a,b,'cputiled',8,iters,[],[],[],[],threads); T = toc;"
            print >> mfile, '            %d=nnf(:,:,3); dvalid=d(1:end-8,1:end-8,:); davg = mean(sqrt(double(dvalid(:))/255^2));'
            print >> mfile, '            nnf_d = double(nnf);'
            print >> mfile, '            nnf_d = nnf_d(1:end-7,1:end-7,:);'
            print >> mfile, "            save('pm_nnf.mat', 'nnf_d');"
            print >> mfile, "            [status, out] = system(sprintf('%(py)s check_dist.py %%s %%s pm_nnf.mat', anames{i}, bnames{i}));" % locals()
            print >> mfile, "            davg = str2num(out);"
        elif method == 'csh':
            print >> mfile, "            T0=clock; nnf=CSH_nn(a,b,8,iters); T = etime(clock, T0);"
            print >> mfile, "            nnf = nnf(1:end-7,1:end-7,:);"
            print >> mfile, "            nnf1 = nnf-1;"
            print >> mfile, "            delete 'csh_nnf.mat';"
            print >> mfile, "            save('csh_nnf.mat', 'nnf1');"
            print >> mfile, "            [status, out] = system(sprintf('%(py)s check_dist.py %%s %%s csh_nnf.mat', anames{i}, bnames{i}));" % locals()
            print >> mfile, "            davg = str2num(out);"
        elif method == 'treecann':
            print >> mfile, "            A=a; B=b;"
            print >> mfile, "            T0=clock; [nnf_dist, nnf_X , nnf_Y, runtime] =run_TreeCANN(uint8(A),uint8(B),patch_w,iters,bgrid, [], [], [], [], [], [], [], 1); T = etime(clock, T0);"
#            run_TreeCANN (A, B, patch_w ,A_grid, B_grid, num_of_train_patches, num_PCA_dims, eps, num_of_ann_matches, A_win, B_win, second_phase, precompute_b)
            if not include_precompute:
                print >> mfile, "        T = runtime(7);"
            print >> mfile, "            bew=size(B,2)-patch_w+1;"
            print >> mfile, "            beh=size(B,1)-patch_w+1;"
            print >> mfile, "            aew=size(A,2)-patch_w+1;"
            print >> mfile, "            aeh=size(A,1)-patch_w+1;"
            print >> mfile, "            nnf=zeros(aeh, aew, 3);"
            print >> mfile, "            nnf(:,:,1)=min(max(nnf_X(1:aeh, 1:aew), 1), bew)-1;"
            print >> mfile, "            nnf(:,:,2)=min(max(nnf_Y(1:aeh, 1:aew), 1), beh)-1;"
            print >> mfile, "            delete 'treecann_nnf.mat';"
            print >> mfile, "            save('treecann_nnf.mat', 'nnf');"
            print >> mfile, "            [status, out] = system(sprintf('%(py)s check_dist.py %%s %%s treecann_nnf.mat', anames{i}, bnames{i}));" % locals()
            print >> mfile, "            davg = str2num(out);"
        print >> mfile, '            dL = [dL davg];'
        print >> mfile, '            tL = [tL T];'
        print >> mfile, '        end'
        method_cap = method.upper()
        if method == 'treecann':
            print >> mfile, """        fprintf('%%f %%f "%(method_cap)s(agrid=%%d, bgrid=%%d, threads=%%d)"\\n', mean(tL), mean(dL), iters, bgrid, threads);""" % locals()
        else:
            print >> mfile, """        fprintf('%%f %%f "%(method_cap)s(iters=%%d, threads=%%d)"\\n', mean(tL), mean(dL), iters, threads);""" % locals()

        if method == 'treecann':
            print >> mfile, '    end'
        print >> mfile, '    end'
        print >> mfile, 'end'
        return
        
    if verbose:
        print >> sys.stderr, 'timeL:', timeL
        print >> sys.stderr, 'distL:', distL
    
    print '%f %f "PatchTable(%s)"' % (numpy.mean(timeL), numpy.mean(distL), switches)
    
if __name__ == '__main__':
    main()
    
