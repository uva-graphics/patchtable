
import os
import pylab
import glob

def hide_legend():
    pylab.gca().legend_ = None #().set_visible(False)
    pylab.draw()

def logplot():
    # Log plot to 'plots' directory
    plotdir = 'plots'
    def filename(idx, legend=True):
        return os.path.join(plotdir, 'plot%05d%s.png'%(idx, ('_legend' if legend else '_nolegend')))

    if not os.path.exists(plotdir):
        os.mkdir(plotdir)
    i = 0
    while os.path.exists(filename(i)):
        i += 1
    #print i
    pylab.legend()
    pylab.savefig(filename(i), dpi=200)
    hide_legend()
    pylab.savefig(filename(i, False), dpi=200)
    pylab.legend()

    L = sorted(glob.glob(os.path.join(plotdir, '*.png')))
    L = [os.path.split(x)[1] for x in L]
    with open(os.path.join(plotdir, 'index.html'), 'wt') as f:
        print >> f, '<html><body><table border=0>'
        print >> f, '<tr><td align=center>With Legend</td><td align=center>Without Legend</td></tr>'
        for i in range(0,len(L),2):
            print >> f, '<tr>'
            print >> f, '<td><img src=%s width=400></td>' % L[i]
            if i+1 < len(L):
                print >> f, '<td><img src=%s width=400></td>' % L[i+1]
            print >> f, '</tr>'
        print >> f, '</table></body></html>'