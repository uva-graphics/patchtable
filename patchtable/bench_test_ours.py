
import os
L = """
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


os.system('rm bench_test_ours.txt')

for n in [10, 30, 60]:
    os.system('echo "Test set %d" >> bench_test_ours.txt' % n)
    os.system('rm -rf ../patchmatch/one_vs_many')
    os.system('cp -r ../patchmatch/one_vs_many%d ../patchmatch/one_vs_many' % n)
    for switches in L:
        os.system('python bench.py ours 10 "%s" -one_vs_many 1 -skip 1 >> bench_test_ours.txt' % switches)
    os.system('echo "" >> bench_test_ours.txt')