
import sys
import os

def system(s):
    print s
    os.system(s)
    
def usage():
    print >> sys.stderr, 'plot_compare.py outdir  -- Plot comparisons with box filters and separable'
    sys.exit(1)

def compare_van_vliet(outdir, kernel, topology):
    system('python vanvliet.py > /dev/null')
    system('python plot_pareto.py pareto_summary.txt %s_%s_pareto_summary.txt > /dev/null'%(os.path.join(outdir, kernel), topology))

def main():
    args = sys.argv[1:]
    if len(args) == 0:
        usage()
    outdir = args[0]
    
    compare_van_vliet(outdir, 'gaussian2', 'iir4')
#    return
    os.system('python plot_separable.py %s 1' % outdir)
    os.system('python plot_box.py %s' % outdir)

if __name__ == '__main__':
    main()
    
