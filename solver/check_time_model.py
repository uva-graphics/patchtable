
import json
import sys
import pylab

def usage():
    print >> sys.stderr, 'check_time_model.py tuned.json'
    sys.exit(1)
    
def main():
    args = sys.argv[1:]
    if len(args) < 1:
        usage()
    L = json.loads(open(args[0], 'rt').read())
    Tmodel = []
    Tactual = []
    E = []
    
    for i in range(len(L)):
        try:
            Tactual.append(L[i]['tune']['time'])
            Tmodel.append(L[i]['time'])
            E.append(L[i]['error'])
        except KeyError:
            pass
    print 'Got %d/%d data points' % (len(Tmodel), len(L))
    
    print 'Writing to time.txt'
    with open('time.txt', 'wt') as f:
        for i in range(len(L)):
            try:
                name = L[i]['name']
                featureL = L[i]['features']
                try:
                    time = L[i]['tune']['time']
                except KeyError:
                    time = float(L[i]['baseline'])
                print >>f, name, ' '.join(str(x) for x in featureL), time
            except KeyError:
                pass
        
    print 'Done'
    
    pylab.xlabel('Time (model)')
    pylab.ylabel('Time (actual)')
    pylab.scatter(Tmodel, Tactual)
    pylab.show()
    
    pylab.xlabel('Error')
    pylab.ylabel('Time (actual)')
    pylab.scatter(E, Tactual)
    pylab.show()

    pylab.xlabel('Error')
    pylab.ylabel('Time (model)')
    pylab.scatter(E, Tmodel)
    pylab.show()
    
    
if __name__ == '__main__':
    main()
    