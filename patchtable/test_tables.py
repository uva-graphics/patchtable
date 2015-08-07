
import os
import time
import subprocess
from bench import float_last_line

cmd = './patchtable match vidpair0/a.png vidpair0/b.png %s -speed 3 -verbose 1 -dt_algo downsample -partition_step 16 -kcoherence 0 -limit %f -ndims %d'
#early_suffix = ''
#early_suffix = ' -do_prop 0'
early_suffix = ' -do_prop 0 -spatial 0 -do_rs 0'
#last_suffix = ''
#last_suffix = ' -kcoherence 1'
last_suffix = ' -kcoherence 10 -kcoherence_step 6'
#early_suffix = ' -do_prop 0 -query_step 2'

#nL = [1.0, 0.5, 2.0, 4.0]
#nL = [8.0]*5
#nL = [8.0]
#dimL = [6]
nL = [1.0, 0.5, 2.0, 4.0, 8.0]
dimL = [6, 7, 6, 7, 6, 7]
#dimL = [6, 7, 8, 6, 7, 8]
#nL = [0.1, 0.2, 0.3, 0.4, 8.0]
#nL = [10.0, 0.5, 0.25, 0.125, 0.125*0.5]
#nL = [1.0, 0.5, 0.25, 0.125, 0.125*0.5]

os.system('rm out_prev*.pfm')

def system(s):
    print s
    return subprocess.check_output(s, shell=True)

T0 = time.time()
tL = []
lookup_tL = []

for (i, n) in enumerate(nL):
    T0_i = time.time()
    suffix = (' -prev_nnf out_prev%d.pfm' % (i-1)) if i > 0 else ''
    cmd_full = cmd % ('out_prev%d.pfm' % i, n, dimL[i]) + suffix
    if i < len(nL)-1:
        cmd_full += early_suffix
    else:
        cmd_full += last_suffix
    s = system(cmd_full)
    tL.append(float_last_line(s, 'Combined precompute and lookup time:'))
    lookup_tL.append(float_last_line(s, '    1-NN overall lookup time:'))
    print s

print 'Individual times:', tL
print 'Individual lookup times:', lookup_tL
print 'Total lookup time:', sum(lookup_tL)
print 'Total time: %f' % sum(tL)
