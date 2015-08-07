
import sys
import os

param_L = """
-do_prop 1 -ndims 7 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 1 -threads 1 -prop_iters 3 -partition_step 1 -limit 100.000000 -kcoherence 20 -ntables 6 -nchroma 1 -dt_iters -1 -triangle_factor 2.000000 -spatial 1
-do_prop 1 -ndims 7 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 1 -threads 1 -prop_iters 3 -partition_step 1 -limit 100.000000 -kcoherence 20 -ntables 4 -nchroma 1 -dt_iters -1 -triangle_factor 1.500000 -spatial 1
-do_prop 1 -ndims 7 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 1 -threads 1 -prop_iters 2 -partition_step 1 -limit 100.000000 -kcoherence 20 -ntables 6 -nchroma 1 -dt_iters -1 -triangle_factor 2.000000 -spatial 1
-do_prop 1 -ndims 7 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 1 -threads 1 -prop_iters 2 -partition_step 1 -limit 100.000000 -kcoherence 20 -ntables 4 -nchroma 1 -dt_iters -1 -triangle_factor 1.500000 -spatial 1
-do_prop 1 -ndims 7 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 1 -threads 1 -prop_iters 2 -partition_step 1 -limit 100.000000 -kcoherence 20 -ntables 3 -nchroma 1 -dt_iters -1 -triangle_factor 1.000000 -spatial 1
-do_prop 1 -ndims 7 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 1 -threads 1 -prop_iters 2 -partition_step 1 -limit 100.000000 -kcoherence 10 -ntables 3 -nchroma 1 -dt_iters -1 -triangle_factor 1.000000 -spatial 1
-do_prop 1 -ndims 7 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 1 -threads 1 -prop_iters 2 -partition_step 1 -limit 100.000000 -kcoherence 4 -ntables 3 -nchroma 1 -dt_iters -1 -triangle_factor 1.500000 -spatial 1
-do_prop 1 -ndims 6 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 1 -threads 1 -prop_iters 1 -partition_step 1 -limit 100.000000 -kcoherence 20 -ntables 3 -nchroma 1 -dt_iters -1 -triangle_factor 1.500000 -spatial 1
-do_prop 1 -ndims 7 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 1 -threads 1 -prop_iters 1 -partition_step 1 -limit 100.000000 -kcoherence 10 -ntables 3 -nchroma 1 -dt_iters -1 -triangle_factor 1.000000 -spatial 1
-do_prop 1 -ndims 7 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 1 -threads 1 -prop_iters 1 -partition_step 1 -limit 100.000000 -kcoherence 4 -ntables 3 -nchroma 1 -dt_iters -1 -triangle_factor 1.000000 -spatial 1
-do_prop 1 -ndims 6 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 2 -threads 1 -prop_iters 2 -partition_step 1 -limit 100.000000 -kcoherence 20 -ntables 4 -nchroma 1 -dt_iters -1 -triangle_factor 2.000000 -spatial 1
-do_prop 1 -ndims 6 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 2 -threads 1 -prop_iters 2 -partition_step 1 -limit 100.000000 -kcoherence 20 -ntables 2 -nchroma 1 -dt_iters -1 -triangle_factor 1.500000 -spatial 1
-do_prop 1 -ndims 6 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 2 -threads 1 -prop_iters 2 -partition_step 1 -limit 100.000000 -kcoherence 10 -ntables 3 -nchroma 1 -dt_iters -1 -triangle_factor 1.500000 -spatial 1
-do_prop 1 -ndims 6 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 2 -threads 1 -prop_iters 2 -partition_step 1 -limit 100.000000 -kcoherence 4 -ntables 3 -nchroma 1 -dt_iters -1 -triangle_factor 2.000000 -spatial 1
-do_prop 1 -ndims 6 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 3 -threads 1 -prop_iters 2 -partition_step 1 -limit 100.000000 -kcoherence 20 -ntables 4 -nchroma 1 -dt_iters -1 -triangle_factor 2.000000 -spatial 1
-do_prop 1 -ndims 6 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 3 -threads 1 -prop_iters 2 -partition_step 1 -limit 100.000000 -kcoherence 20 -ntables 4 -nchroma 1 -dt_iters -1 -triangle_factor 1.000000 -spatial 1
-do_prop 1 -ndims 6 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 3 -threads 1 -prop_iters 2 -partition_step 1 -limit 100.000000 -kcoherence 10 -ntables 3 -nchroma 1 -dt_iters -1 -triangle_factor 1.000000 -spatial 1
-ndims 6 -nchroma 1 -limit 100.000000 -do_rs 0 -threads 1 -run_dt 1 -partition_step 1 -dt_iters -1 -do_prop 1 -query_step 3 -kcoherence 4 -kcoherence_iter -1 -prop_iters 2 -spatial 1 -ntables 2 -triangle_factor 1.500000
-do_prop 1 -ndims 6 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 3 -threads 1 -prop_iters 1 -partition_step 1 -limit 100.000000 -kcoherence 20 -ntables 4 -nchroma 1 -dt_iters -1 -triangle_factor 1.000000 -spatial 1
-do_prop 1 -ndims 6 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 3 -threads 1 -prop_iters 1 -partition_step 1 -limit 100.000000 -kcoherence 4 -ntables 4 -nchroma 1 -dt_iters -1 -triangle_factor 2.000000 -spatial 1
-do_prop 1 -ndims 6 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 3 -threads 1 -prop_iters 1 -partition_step 1 -limit 100.000000 -kcoherence 3 -ntables 3 -nchroma 1 -dt_iters -1 -triangle_factor 2.000000 -spatial 1
-do_prop 1 -ndims 6 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 3 -threads 1 -prop_iters 1 -partition_step 1 -limit 100.000000 -kcoherence 1 -ntables 2 -nchroma 1 -dt_iters -1 -triangle_factor 4.000000 -spatial 1
-do_prop 1 -ndims 6 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 3 -threads 1 -prop_iters 1 -partition_step 1 -limit 100.000000 -kcoherence 1 -ntables 1 -nchroma 1 -dt_iters -1 -triangle_factor 2.000000 -spatial 1
-do_prop 1 -ndims 6 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 4 -threads 1 -prop_iters 1 -partition_step 1 -limit 100.000000 -kcoherence 3 -ntables 2 -nchroma 1 -dt_iters -1 -triangle_factor 1.000000 -spatial 1
-do_prop 1 -ndims 6 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 4 -threads 1 -prop_iters 1 -partition_step 1 -limit 100.000000 -kcoherence 1 -ntables 1 -nchroma 1 -dt_iters -1 -triangle_factor 4.000000 -spatial 1
-do_prop 1 -ndims 6 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 6 -threads 1 -prop_iters 1 -partition_step 1 -limit 100.000000 -kcoherence 3 -ntables 3 -nchroma 1 -dt_iters -1 -triangle_factor 2.000000 -spatial 1
-do_prop 1 -ndims 6 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 6 -threads 1 -prop_iters 1 -partition_step 1 -limit 100.000000 -kcoherence 3 -ntables 1 -nchroma 1 -dt_iters -1 -triangle_factor 2.000000 -spatial 1
-do_prop 1 -ndims 6 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 6 -threads 1 -prop_iters 1 -partition_step 1 -limit 100.000000 -kcoherence 2 -ntables 1 -nchroma 1 -dt_iters -1 -triangle_factor 1.000000 -spatial 0
-do_prop 1 -ndims 6 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 10 -threads 1 -prop_iters 1 -partition_step 1 -limit 100.000000 -kcoherence 1 -ntables 3 -nchroma 1 -dt_iters -1 -triangle_factor 1.000000 -spatial 1
-do_prop 1 -ndims 6 -kcoherence_iter -1 -run_dt 1 -do_rs 0 -query_step 10 -threads 1 -prop_iters 1 -partition_step 1 -limit 100.000000 -kcoherence 0 -ntables 1 -nchroma 1 -dt_iters -1 -triangle_factor 1.000000 -spatial 0
""".strip().split('\n')

print len(param_L)
sys.exit(1)

full = True

if full:
    cmd_spec = 'python bench.py ours 10 "%s" -one_vs_many 1 -skip 1'
else:
    cmd_spec = 'python bench.py ours 2 "%s" -one_vs_many 1 -skip 1'

def system(s):
    #print >> sys.stderr, s
    os.system(s)

for i in [10]: #[6,2,3,4,0]:
    if i == 0:
        args = ''
        name = 'Full Algorithm'
    elif i == 1:
        args = '-kcoherence 0'
        name = 'No k-coherence'
    elif i == 2:
        args = '-spatial 0'
        name = 'No spatial search'
    elif i == 3:
        args = '-query_step 1'
        name = 'No sparse grid'
    elif i == 4:
        args = '-init_random 1 -do_table_lookup 0'
        name = 'No table lookup'
    elif i == 5:
        args = '-kcoherence 0 -spatial 0'
        name = 'No k-coherence, no spatial search'
    elif i == 6:
        args = '-ntables 1'
        name = 'One table only'
    elif i == 7:
        args = '-ntables 1 -kcoherence 0'
        name = 'One table only, no k-coherence'
    elif i == 8:
        args = '-triangle_factor 10000000'
        name = 'No early termination'
    elif i == 9:
        args = '-prop_iters 1'
        name = 'One query iteration'
    elif i == 10:
        args = '-run_dt 0'
        name = 'No distance transform'
    if len(args):
        args = ' ' + args
    print '=' * 80
    print name
    print '=' * 80
    print
    sys.stdout.flush()
    n = len(param_L)
    if not full:
        n = 1
    j_min = 0
    if i == 10:
        j_min = 13 + 4
    for j in range(j_min, n):
        param = param_L[j] + args
        cmd = cmd_spec % param
        system(cmd + ' 2> /dev/null')
        sys.stdout.flush()
