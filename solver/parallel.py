
import sys
import os
from compare import qsub_exists, qsub_list
import hashlib
import subprocess
import time
import getpass
import traceback

def main():
    args = sys.argv[1:]
    if len(args) == 0:
        print >> sys.stderr, 'python parallel.py scriptfile'
        print >> sys.stderr, '  Runs lines of scriptfile in parallel, using GNU parallel or qsub (if available)'
        sys.exit(1)
        
    filename = args[0]
    s = open(filename, 'rt').read()
    #print s
    
    if not qsub_exists:
        os.system('parallel < %s' % filename)
    else:
        L = s.strip().split('\n')
        sig = hashlib.md5(s+str(time.time())).hexdigest()[:15]
        name_list = [sig for x in L]
        qsub_list(L, '', name_list)
        username = getpass.getuser()
        while True:
            try:
                s = subprocess.check_output('qstat', shell=True)
                L = s.strip().split('\n')
                L = [x for x in L if (username in x and sig in x)]
                n = len(L)
                if n == 0:
                    break
            except:
                traceback.print_exc()
            time.sleep(5.0)
            
if __name__ == '__main__':
    main()
