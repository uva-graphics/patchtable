
import sys
import string
import subprocess
import parse_args
import os

def usage():
    print >> sys.stderr, 'python ga.py kernel.txt -out [outdir] -filter_w [maxw] [more args to solver]'
    print >> sys.stderr, '  Runs GA on FIR cluster. Then re-solves for the top topologies.'
    print >> sys.stderr, '  Use -out /bigtmp/your_userid/outdir which has sufficient space.'
    sys.exit(1)

def main():
    def add_default_param(p0, add_str):
        (prefix, suffix) = add_str.split()
        if prefix in p0:
            return p0
        return p0 + ' ' + add_str
        
    args = sys.argv[1:]
    if len(args) < 1:
        usage()
        
    params = ' '.join(args)
    
    params = add_default_param(params, '-iters 100')
    params = add_default_param(params, '-ga_population 100')
    params = add_default_param(params, '-ga_generations 30')
    params = add_default_param(params, '-verbose 2')
    
    cmd = """echo "cd %s; ./solve_kernel ga %s" | qsub -N filter-ga -o ga_out.txt -e ga_err.txt -l mem=40gb""" % (os.getcwd(), params)
    print cmd
    jobid = subprocess.check_output(cmd, shell=True)
    
    outdir = parse_args.parse_args(allow_all=True)[1]['out']
    
    if len(jobid) and jobid[0] in string.digits:
        cmd_p = """echo "cd %s; python resolve.py %s" | qsub -N filter-ga-res -o resolve_out.txt -e resolve_err.txt -l mem=40gb -W  depend=afterany:%s""" % (os.getcwd(), outdir, jobid)
        print cmd_p
        os.system(cmd_p)
    else:
        print >> sys.stderr, 'qsub failed: %s' % jobid
        sys.exit(1)
    
if __name__ == '__main__':
    main()

