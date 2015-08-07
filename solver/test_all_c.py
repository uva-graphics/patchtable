
import sys
import os
import json

def main():
    args = sys.argv[1:]
    if len(args) < 1:
        print >> sys.stderr, 'test_all_c.py pareto.json [outdir]  -- Tests compiler on all programs in Pareto frontier'
        sys.exit(1)
    
    filename = args[0]
    outdir = '.temp_out'
    if len(args) > 1:
        outdir = args[1]
    print 'Writing test output to %s' % outdir
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    L = json.loads(open(filename, 'rt').read())
    outL = [os.path.join(outdir, 'f%05d.cpp'%i) for i in range(len(L))]
    compileL = ['python compile_to_c.py "%s" %d "%s" -compile 1' % (filename, i, outL[i]) for i in range(len(L))]

    compile_script = os.path.join(outdir, 'compile_all.sh')
    with open(compile_script, 'wt') as f:
        f.write('\n'.join(compileL))
    os.system('python parallel.py %s' % compile_script)

    for i in range(len(L)):
        if os.system('python compile_to_c.py "%s" %d "%s" -compile 0 -test 1' % (filename, i, outL[i])):
            print 'Failed on %d' % i
            sys.exit(1)

if __name__ == '__main__':
    main()

