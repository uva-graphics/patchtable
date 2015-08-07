
"""
analogy.py -- Video Analogies without Temporal Coherence

python analogy.py a.png aprime.png b_video_filespec.png out_video_filespec.png

Example:
python analogy.py a.png ap.png video/in/bike%05d.png video/out/bike%05d.png
"""

import sys
import os
import os.path

def system(s):
    print s
    os.system(s)
    
def main():
    args = sys.argv[1:]
    if len(args) < 4:
        print >> sys.stderr, __doc__
        sys.exit(1)
        
    (a, ap) = args[:2]
    
    for i in range(1, 1000**2):
        b = args[2] % i
        bp = args[3] % i
        if not os.path.exists(b):
            break
        print '-'*72
        print 'Frame %d' % i
        print '-'*72
        print
        system('./patchtable analogy %s %s %s %s' % (a, ap, b, bp))        
        
if __name__ == '__main__':
    main()
    