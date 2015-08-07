
import sys
import numpy
import re

speed_count = 10 #30

def main():
    args = sys.argv[1:]
    if len(args) == 0:
        print >> sys.stderr, 'Generates code for set_speed() in patchtable.cpp by sampling a Pareto frontier at equally spaced (log) time.'
        print >> sys.stderr, ' $ python ../solver/pareto.py bench_output.txt > bench_pareto.txt'
        print >> sys.stderr, ' $ python ../solver/gen_speed.py bench_pareto.txt'
        print >> sys.stderr, 'Prints to stdout the selected lines from the input text file.'
        sys.exit(1)
    
    L = open(args[0], 'rU').read().strip().split('\n')
    L = [x.strip() for x in L]
    timeL = [float(x.split(' ', 2)[0]) for x in L]
    log_timeL = numpy.log(timeL)

    print >> sys.stderr, L #'timeL length:', len(timeL)#timeL
    
    targetL = numpy.linspace(log_timeL[0], log_timeL[-1], speed_count)
    idxL = []
    min_idx = 0
    for i in range(len(targetL)):
        selected = min_idx
        #print i, min_idx, len(targetL), len(log_timeL)
        dselected = (targetL[i]-log_timeL[min_idx])**2
        for j in range(min_idx, len(log_timeL)):
            dcurrent = (targetL[i]-log_timeL[j])**2
            if dcurrent < dselected:
                dselected = dcurrent
                selected = j
        idxL.append(selected)
        min_idx = selected+1
    
    for (i, idx) in enumerate(idxL):
        print >> sys.stderr, L[idx]
    print >> sys.stderr
    
    print 'void PatchTableParams::set_speed(int i) {'
    for (i, idx) in enumerate(idxL):
        current = L[idx]
        current = current.replace('"PatchTable(', 'PatchTable(')
        current = current.replace(')"', ')')
        if i == 0:
            sys.stdout.write('    if      (i == 0) { ')
        else:
            sys.stdout.write('    else if (i == %d) { '%i)
        switches = re.findall(re.escape('PatchTable(') + '(.*?)' + re.escape(')'), current)
        assert len(switches) == 1, switches
        switches = switches[0]
        switches = switches.replace('-threads 1', '')
        while '  ' in switches:
            switches = switches.replace('  ', ' ')  
        sys.stdout.write('set_from_switches(string("%s")); }\n' % switches)
    print r'    else { fprintf(stderr, "Error: set_speed(i), i=%d, expected integer in 0...9\n"); exit(1); }'
    print '}'
    
if __name__ == '__main__':
    main()
    
