
import sys
import json
import glob
import string
import os, os.path

def usage():
    print >> sys.stderr, 'python resolve.py dirname [additional args for solver]'
    print >> sys.stderr, '  Re-runs the solver for each of the top topologies (determined by collect_pareto.py)'
    print >> sys.stderr, '  Output is assumed to be from the G.A. which stores the commandline for solver (in solver_cmd)'
    print >> sys.stderr, '  Supplies same arguments as original solver call by default.'
    sys.exit(1)

def system(s):
    print s
    os.system(s)
    
def main():
    args = sys.argv[1:]
    if len(args) == 0:
        usage()
        
    dirname = args[0]
    solver_extra_args = (' ' + ' '.join(args[1:])) if len(args)>1 else '' 

    rank_filename = os.path.join(dirname, 'rank_top.txt')
    system('python collect_pareto.py %s -rank_top %s' % (dirname, rank_filename))
    
    prefix = 'g'
    while glob.glob(os.path.join(dirname, prefix + '*pareto_summary.txt')):
        prefix = string.ascii_lowercase[string.ascii_lowercase.index(prefix)+1]
        
    cmdL = []
    for (i, filename) in enumerate(open(rank_filename, 'rt').read().strip().split('\n')):
        pareto_full = json.loads(open(filename, 'rt').read())
        cmd = pareto_full[0]['solver_cmd']
        outfile = os.path.join(dirname, prefix + ('%05d'%i))
        cmd += ' -out %s -resume %s%s' % (outfile, filename, solver_extra_args)
        cmdL.append(cmd)

    scriptfile = os.path.join(dirname, 'resolve.sh')
    with open(scriptfile, 'wt') as f:
        print >> f, '\n'.join(cmdL)
    
    print '\n'.join(cmdL)
    
    system('python parallel.py %s' % scriptfile)
    
    system('python collect_pareto.py %s %s' % (dirname, os.path.join(dirname, 'collect_top.json')))
    system('python collect_pareto.py %s %s -subsample 2' % (dirname, os.path.join(dirname, 'collect_top_sub2.json')))
    
if __name__ == '__main__':
    main()
    