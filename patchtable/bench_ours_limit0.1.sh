python bench.py ours 6 "-ndims 6 -limit 0.1 -do_rs 1 -randomize_dt 1 -threads 1 -populate_nearest 1" > bench_ours_limit0.1.txt
python bench.py ours 6 "-ndims 6 -limit 0.1 -do_rs 1 -randomize_dt 0 -threads 1 -populate_nearest 1" >> bench_ours_limit0.1.txt
python bench.py ours 6 "-ndims 6 -limit 0.1 -do_rs 0 -randomize_dt 0 -threads 1 -populate_nearest 1" >> bench_ours_limit0.1.txt
python bench.py ours 6 "-ndims 7 -limit 0.1 -do_rs 1 -randomize_dt 1 -threads 1 -populate_nearest 1" >> bench_ours_limit0.1.txt
python bench.py ours 6 "-ndims 7 -limit 0.1 -do_rs 1 -randomize_dt 0 -threads 1 -populate_nearest 1" >> bench_ours_limit0.1.txt
python bench.py ours 6 "-ndims 7 -limit 0.1 -do_rs 0 -randomize_dt 0 -threads 1 -populate_nearest 1" >> bench_ours_limit0.1.txt
python bench.py ours 6 "-ndims 8 -limit 0.1 -do_rs 1 -randomize_dt 1 -threads 1 -populate_nearest 1" >> bench_ours_limit0.1.txt
python bench.py ours 6 "-ndims 8 -limit 0.1 -do_rs 1 -randomize_dt 0 -threads 1 -populate_nearest 1" >> bench_ours_limit0.1.txt
python bench.py ours 6 "-ndims 8 -limit 0.1 -do_rs 0 -randomize_dt 0 -threads 1 -populate_nearest 1" >> bench_ours_limit0.1.txt

python bench.py ours 6 "-ndims 6 -limit 0.1 -do_rs 0 -randomize_dt 0 -do_prop 0 -threads 1 -populate_nearest 1" >> bench_ours_limit0.1.txt
python bench.py ours 6 "-ndims 7 -limit 0.1 -do_rs 0 -randomize_dt 0 -do_prop 0 -threads 1 -populate_nearest 1" >> bench_ours_limit0.1.txt
python bench.py ours 6 "-ndims 8 -limit 0.1 -do_rs 0 -randomize_dt 0 -do_prop 0 -threads 1 -populate_nearest 1" >> bench_ours_limit0.1.txt
