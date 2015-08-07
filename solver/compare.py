
import glob
import os
import sys
import forkmap
from distutils.spawn import find_executable

skip_existing = True
print_only = False #False

def system(s):
    print s
    if not print_only:
        os.system(s)

def qsub_list(L, qsub_args='', name_list=None):
    for (i, x) in enumerate(L):
        name = 'filter-approx'
        if name_list is not None:
            name = name_list[i]
        system('''echo "cd %s; %s" | qsub -N %s -o /dev/null -e /dev/null -l select=1:ncpus=16:mem=16gb %s''' % (os.getcwd(), x, name, qsub_args))

qsub_exists = find_executable('qsub') is not None

def main():
    global print_only
    global skip_existing
    args = sys.argv[1:]
    if len(args) == 0:
        print >> sys.stderr, 'compare out_dir [--print] [--print-all] [options for solver]'
        sys.exit(1)
    out_dir = args[0]
    if len(args) >= 2 and args[1] == '--print':
        print_only = True
        args = args[1:]
    if len(args) >= 2 and args[1] == '--print-all':
        print_only = True
        skip_existing = False
        args = args[1:]
    #print print_only
    #print skip_existing
    #print args
    #return
    rest = ' '.join(args[1:])
    try:
        os.mkdir(out_dir)
    except:
        pass
    L = []
    kernels = glob.glob('kernels/*.txt')
    for topology in ['fir_fir', 'fir_iir', 'firh_firv', 'iir_iir', 'iir_iir2', 'iir_iir3', 'iir4']:
        for kernel in kernels:
            if topology == 'iir_iir' and 'box' not in kernel:
                continue
            if topology in ['iir_iir3', 'iir_iir2'] and not kernel.endswith('gaussian4.txt'):
                continue
            if topology == 'iir4' and not kernel.endswith('gaussian2.txt'):
                continue
            prefix = os.path.splitext(os.path.split(kernel)[1])[0]
            out_filename = os.path.join(out_dir, '%s_%s' % (prefix, topology))
            pareto = out_filename + '_pareto_summary.txt'
            filter_w = 7
            if 'gabor19' in kernel or 'gaussian4' in kernel:
                filter_w = 19
            elif 'gaborasym' in kernel:
                filter_w = 17
            elif 'box3' in kernel or 'box4' in kernel or 'box5' in kernel or 'gaussian2' in kernel or 'gaussian4' in kernel or topology == 'firh_firv' or 'gabor' in kernel:
                filter_w = 15
            if 'gaussian8' in kernel or 'gabor19' in kernel:
                continue
                #filter_w = 19                
            if (not skip_existing) or (not os.path.exists(pareto)): # or print_only:
                L.append('./solve_kernel solve %s -topology %s -filter_w %d -out %s %s > %s' % (kernel, topology, filter_w, out_filename, rest, out_filename + '_out.txt'))
#    threadmap.map(system, L, n=1)
    #forkmap.map(system, L)
    if qsub_exists:
        qsub_list(L)
    elif print_only:
        map(system, L)
    else:
        forkmap.map(system, L, n=8)
    
if __name__ == '__main__':
    main()

