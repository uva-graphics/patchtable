
import os
import sys

def usage():
    print >> sys.stderr, 'plot_box outdir'
    sys.exit(1)

def system(s):
    print s
    os.system(s)

def compare(outdir, kernel, nbox, topology):
    system('python box_filters.py kernels/%s.txt -nmin %d -nmax %d > /dev/null' % (kernel, nbox, nbox))
    system('python plot_pareto.py pareto_summary.txt %s_%s_pareto_summary.txt > /dev/null'%(os.path.join(outdir, kernel), topology))

def main():
    args = sys.argv[1:]
    if len(args) == 0:
        usage()
    outdir = args[0]
    
    for kernel in ['box2', 'box3', 'box4', 'box5']:
        compare(outdir, kernel, 1, 'iir_iir')
    
    for kernel in ['gaussian4']:
        for (nbox, topology) in [(2, 'iir_iir2'), (3, 'iir_iir3')]:
            compare(outdir, kernel, nbox, topology)
            
if __name__ == '__main__':
    main()
    