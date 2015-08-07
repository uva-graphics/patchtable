
import sys
import glob
import os

def system(s):
    print s
    os.system(s)
    
def main():
    args = sys.argv[1:]
    if len(args) == 0:
        print >> sys.stderr, 'plot_separable out_dir [compare_all=0[|1]]'
        sys.exit(1)
    out_dir = args[0]
    compare_all = False
    if len(args) >= 2:
        compare_all = bool(int(args[1]))
    
    suffix = 'firh_firv_pareto_summary.txt'
    suffix2 = 'fir_fir_pareto_summary.txt'
    filenames = glob.glob(os.path.join(out_dir, '*' + suffix))
    for filename in filenames:
        pre = os.path.split(filename)[1]
        pre = pre.split('_')[0]
        #pre = pre.replace('gabor_asym', 'gaborasym').split('_')[0].replace('gaborasym', 'gabor_asym')
        filename2 = os.path.join(out_dir, pre + '_' + suffix2)
        ok = os.path.exists(filename2)
        if not ok:
            print filename, filename2, 'missing'
        #print filename, filename2, 'OK' if ok else 'missing'
        if ok:
            if compare_all:
                all_L = [x for x in glob.glob(os.path.join(out_dir, pre + '_' + '*.txt')) if not '_full.txt' in x and not '_out.txt' in x]
                system('python plot_pareto.py %s' % (' '.join(all_L)))
            else:
                system('python plot_pareto.py %s %s' % (filename, filename2))
    
if __name__ == '__main__':
    main()
    
