
import pylab
import numpy
import sys
import logplot
import json
from parse_args import parse_args
import matplotlib.cm

T_max = 0.0
E_max = 0.0

marker_str = '.+xsD8^p*v<>h_|d' * 100
color_list = ['#3bac5a', '#0d90c2', '#8d69b3', '#ff4d1a', '#b3591f', '#f0e15e']*100
colors = matplotlib.cm.rainbow(numpy.linspace(0, 1, 7))

def plot_array(A, label='', autolim=False, markeridx=0, scatter=False):
    global T_max, E_max
    #print markeridx, marker_str[markeridx]
    
    if not scatter:
        pylab.plot(A[:,1], A[:,0], marker_str[markeridx] + '-', label=label) #, c=color_list[markeridx])
    else:
        pylab.scatter(A[:,1], A[:,0], label=label, color=colors[markeridx%len(colors)]) #, marker_str[markeridx]) #, label=label)
#    kw = {}
#    if marker_str[markeridx] == 'o':
#        kw = dict(facecolors='none')
    #pylab.scatter(A[:,1], A[:,0], marker=marker_str[markeridx], c=color_list[markeridx], **kw)
    pylab.xlabel('Error')
    pylab.ylabel('Time')
#    if autolim:
#    if markeridx > 1:
    if 'firh_firv' in label:
        T_max = max(T_max, numpy.max(A[:,0]))
    E_max = max(E_max, numpy.max(A[:,1]))


def usage():
    print >> sys.stderr, 'plot_pareto.py trace1.json|pareto1.txt [trace2.json|pareto2.txt] [-legend b]'
    sys.exit(1)
    
def main():
    import plot_trace
    import plot_converge
#    args = sys.argv[1:]
#    if len(args) == 0:
#        args = ['pareto.txt']
    (args, kw) = parse_args('legend scatter'.split(), usage)
    legend = int(kw.get('legend', '1'))
    scatter = int(kw.get('scatter', '0'))
    
    for (i, filename) in enumerate(args):
        success = False
        if filename.endswith('.json'):
            try:
                A = plot_trace.read_trace(filename)[0][-1]['pareto']
                success = True
            except KeyError:
                pass
        if not success:
            try:
                if filename.endswith('.txt'):
                    A = numpy.loadtxt(filename)
                else:
                    raise ValueError
            except ValueError:
                J = json.loads(open(filename, 'rt').read())
                A = numpy.array([[J[j]['time'], J[j]['error']] for j in range(len(J))]) 
            if len(A.shape) == 1:
                A = numpy.array([A])
        print filename, plot_converge.avg_time_func(A)
        plot_array(A, filename, autolim=True, markeridx=i, scatter=scatter)
    pylab.xlim(-0.005, E_max+0.005)
    if T_max:
        pylab.ylim(0, T_max+0.2)
    logplot.logplot()
    if legend:
        pylab.legend()
    pylab.show()

if __name__ == '__main__':
    main()

