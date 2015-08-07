
import sys
import json
import os, os.path
import subprocess
import md5
import time
import traceback
from parse_args import parse_args

def usage():
    print >> sys.stderr, 'time_all_c.py pareto_in.json pareto_out.json [outdir] [-image img.png] [-parallel 0|1] [-compile 1[|0]]'
    print >> sys.stderr, '  Times all members of pareto_in.json, writing pareto_out.json, using compile_to_c.'
    print >> sys.stderr, '  Writes programs to outdir which is a temp dir if not specified.'
    sys.exit(1)

def system(s):
    print s
    os.system(s)
    
def main():
    (args, kw) = parse_args('image parallel compile'.split(), usage)
    if len(args) < 2:
        usage()

    image = kw.get('image', 'bird.png')
    parallel = int(kw.get('parallel', '1'))
    compile = int(kw.get('compile', '1'))
    in_filename = args[0]
    out_filename = args[1]
    
    is_temp = True
    outdir = '.temp%s' % md5.md5(str(time.time())).hexdigest()
    if len(args) > 2:
        outdir = args[2]
        is_temp = False

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    L = json.loads(open(in_filename, 'rt').read())        
    outL = [os.path.join(outdir, 'f%05d.cpp'%i) for i in range(len(L))]
    compileL = ['python compile_to_c.py "%s" %d "%s" -compile 1' % (in_filename, i, outL[i]) for i in range(len(L))]
    
    #for i in range(len(L)):
    #    system(compileL[i])         # Do not compile with C compiler for first pass
    
    if parallel and compile:
        compile_script = os.path.join(outdir, 'compile_all.sh')
        with open(compile_script, 'wt') as f:
            f.write('\n'.join(compileL))
        os.system('python parallel.py %s' % compile_script)
    
    for i in range(len(L)):
        try:
            if (not parallel) and compile:
                os.system(compileL[i])
            binary = os.path.abspath(os.path.splitext(outL[i])[0])
            out_image = os.path.join(outdir, 'f%05d.png'%i)
            T = subprocess.check_output('%s %s %s' % (binary, image, out_image), shell=True)
            try:
                T = float(T)
            except:
                print 'Error on program %d' % i
                continue
            print 'Program %d: %f sec' % (i, T)
            L[i]['tune'] = {'time': T}
            L[i]['time'] = T
        except:
            traceback.print_exc()
    
    with open(out_filename, 'wt') as f:
        f.write(json.dumps(L))
        
if __name__ == '__main__':
    main()
    