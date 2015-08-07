
import os
L = """
-kcoherence 7 -limit 0.01 -init_random 1 -prop_iters 2 -query_step 1 -do_table_lookup 0 -verbose 1 -spatial 0 -rs_max 0 -triangle_factor 1000000000
-kcoherence 6 -limit 0.01 -init_random 1 -prop_iters 2 -query_step 1 -do_table_lookup 0 -verbose 1 -spatial 0 -rs_max 0 -triangle_factor 1000000000
-kcoherence 5 -limit 0.01 -init_random 1 -prop_iters 2 -query_step 1 -do_table_lookup 0 -verbose 1 -spatial 0 -rs_max 0 -triangle_factor 1000000000
-kcoherence 4 -limit 0.01 -init_random 1 -prop_iters 2 -query_step 1 -do_table_lookup 0 -verbose 1 -spatial 0 -rs_max 0 -triangle_factor 1000000000
-kcoherence 3 -limit 0.01 -init_random 1 -prop_iters 2 -query_step 1 -do_table_lookup 0 -verbose 1 -spatial 0 -rs_max 0 -triangle_factor 1000000000
-kcoherence 2 -limit 0.01 -init_random 1 -prop_iters 2 -query_step 1 -do_table_lookup 0 -verbose 1 -spatial 0 -rs_max 0 -triangle_factor 1000000000
-kcoherence 12 -limit 0.01 -init_random 1 -prop_iters 1 -query_step 1 -do_table_lookup 0 -verbose 1 -spatial 0 -rs_max 0 -triangle_factor 1000000000
-kcoherence 9 -limit 0.01 -init_random 1 -prop_iters 1 -query_step 1 -do_table_lookup 0 -verbose 1 -spatial 0 -rs_max 0 -triangle_factor 1000000000
-kcoherence 8 -limit 0.01 -init_random 1 -prop_iters 1 -query_step 1 -do_table_lookup 0 -verbose 1 -spatial 0 -rs_max 0 -triangle_factor 1000000000
-kcoherence 6 -limit 0.01 -init_random 1 -prop_iters 1 -query_step 1 -do_table_lookup 0 -verbose 1 -spatial 0 -rs_max 0 -triangle_factor 1000000000
-kcoherence 4 -limit 0.01 -init_random 1 -prop_iters 1 -query_step 1 -do_table_lookup 0 -verbose 1 -spatial 0 -rs_max 0 -triangle_factor 1000000000
-kcoherence 3 -limit 0.01 -init_random 1 -prop_iters 1 -query_step 1 -do_table_lookup 0 -verbose 1 -spatial 0 -rs_max 0 -triangle_factor 1000000000
-kcoherence 2 -limit 0.01 -init_random 1 -prop_iters 1 -query_step 1 -do_table_lookup 0 -verbose 1 -spatial 0 -rs_max 0 -triangle_factor 1000000000
-kcoherence 1 -limit 0.01 -init_random 1 -prop_iters 1 -query_step 1 -do_table_lookup 0 -verbose 1 -spatial 0 -rs_max 0 -triangle_factor 1000000000
""".strip().split('\n')

#print len(L)
#sys.exit(1)

os.system('rm bench_test_kcoherence.txt')

for n in [60]: #[10, 30, 60]:
    os.system('echo "Test set %d" >> bench_test_kcoherence.txt' % n)
    os.system('rm -rf ../patchmatch/one_vs_many')
    os.system('cp -r ../patchmatch/one_vs_many%d ../patchmatch/one_vs_many' % n)
    for switches in L:
        os.system('python bench.py ours 10 "%s" -one_vs_many 1 -skip 1 >> bench_test_kcoherence.txt' % switches)
    os.system('echo "" >> bench_test_kcoherence.txt')
