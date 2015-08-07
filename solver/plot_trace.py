
import pylab
import numpy
import sys
import json
import os, os.path
import itertools
import time
import plot_pareto

RETRIES = 4

def read_trace(filename):
#    done = True
#    if not s.endswith(']]}]'):
#        done = False
#        s += ']'
    for it in range(RETRIES):
        s = open(filename, 'rt').read()
        try:
            try:
                L = json.loads(s)
                done = True
            except:
                L = json.loads(s + ']')
                done = False
            break
        except:
            if it < RETRIES-1:
                time.sleep(0.5)
            else:
                raise
            
    for x in L:
        x['pareto'] = numpy.array(x['pareto'])
    return (L, done)

def main():
    args = sys.argv[1:]
    if len(args) < 2:
        print >> sys.stderr, 'plot_trace trace.json plotdir [step]'
        sys.exit(1)
        
    ihandled = set()
    for it in itertools.count(0):
        #L = [numpy.array(x) for x in L]
        (L, done) = read_trace(args[0])
        outdir = args[1]
        step = 1
        if len(args) >= 3:
            step = int(args[2])
        T_max = max(numpy.max(x['pareto'][:,0]) for x in L)
        E_max = max(numpy.max(x['pareto'][:,1]) for x in L)
        
        if it == 0:
            if not os.path.exists(outdir):
                os.mkdir(outdir)
            os.system('open %s' % outdir)
            
        for i in range(0,len(L),step):
            if i in ihandled:
                continue
            ihandled.add(i)
            pylab.clf()
            A = L[i]['pareto']
            #pylab.plot(A[:,1], A[:,0])
            #pylab.scatter(A[:,1], A[:,0], facecolors='none')
            #pylab.xlabel('Error')
            #pylab.ylabel('Time')
            plot_pareto.plot_array(A)
            pylab.xlim(0, E_max)
            pylab.ylim(0, T_max)
            filename = os.path.join(outdir, 'plot%04d.png'%i)
            pylab.savefig(filename, dpi=200)
            print filename
    #    pylab.legend()
    #    pylab.show()
        time.sleep(1.0)
        if done:
            break

if __name__ == '__main__':
    main()


