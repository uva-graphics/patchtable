
import os, os.path
import sys
import time
from pipes import quote
import random
#from ForkedWatchdog import Watchdog

#TIMEOUT = 10.0

# From logconsolidate.py
arch_d = {'E5649': ['lc5-compute-3-24', 'lc5-compute-3-26', 'lc5-compute-3-29', 'lc5-compute-3-31', 'lc5-compute-3-33', 'lc5-compute-3-34', 'lc5-compute-3-35', 'lc5-compute-3-37', 'lc5-compute-3-39', 'lc5-compute-3-40', 'lc5-compute-3-41', 'lc5-compute-3-42', 'lc5-compute-3-46'], 'X5550': ['lc5-compute-4-0', 'lc5-compute-4-1', 'lc5-compute-4-10', 'lc5-compute-4-11', 'lc5-compute-4-12', 'lc5-compute-4-13', 'lc5-compute-4-14', 'lc5-compute-4-15', 'lc5-compute-4-16', 'lc5-compute-4-17', 'lc5-compute-4-18', 'lc5-compute-4-2', 'lc5-compute-4-20', 'lc5-compute-4-21', 'lc5-compute-4-22', 'lc5-compute-4-23', 'lc5-compute-4-24', 'lc5-compute-4-25', 'lc5-compute-4-26', 'lc5-compute-4-27', 'lc5-compute-4-28', 'lc5-compute-4-29', 'lc5-compute-4-3', 'lc5-compute-4-31', 'lc5-compute-4-32', 'lc5-compute-4-34', 'lc5-compute-4-35', 'lc5-compute-4-36', 'lc5-compute-4-37', 'lc5-compute-4-38', 'lc5-compute-4-39', 'lc5-compute-4-4', 'lc5-compute-4-40', 'lc5-compute-4-41', 'lc5-compute-4-42', 'lc5-compute-4-43', 'lc5-compute-4-44', 'lc5-compute-4-45', 'lc5-compute-4-46', 'lc5-compute-4-47', 'lc5-compute-4-48', 'lc5-compute-4-49', 'lc5-compute-4-5', 'lc5-compute-4-50', 'lc5-compute-4-51', 'lc5-compute-4-52', 'lc5-compute-4-54', 'lc5-compute-4-55', 'lc5-compute-4-6', 'lc5-compute-4-8', 'lc5-compute-4-9']}

def system(s):
    print >> sys.stderr, s
    os.system(s)
    
def main():
    args = sys.argv[1:]
    if len(args) < 2:
        print >> sys.stderr, 'qsub_arch.py ARCH "cmdline" [qsub args]'
        print >> sys.stderr, 'Equivalent to echo "cd $PWD; command line" | qsub [qsub args], but only allow certain architectures.'
        print >> sys.stderr, '  Use any for ARCH to allow any architecture.'
        print >> sys.stderr
        print >> sys.stderr, 'Architectures: X5550'
        sys.exit(1)
    
    arch = args[0].upper()
    cmdline = args[1]
    qsub_args = ' '.join(args[2:])
    
    pyfile = os.path.abspath(__file__)
    
    suffix = ''
    if arch != 'any':
        suffix = ' -l host=%s' % random.choice(arch_d[arch])

    system('echo "cd %s; %s" | qsub %s%s' % (quote(os.getcwd()), cmdline, qsub_args, suffix))

if __name__ == '__main__':
    main()
    
