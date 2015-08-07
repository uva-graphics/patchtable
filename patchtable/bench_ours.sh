python bench.py ours 6 "-ndims 6 -limit 340 -do_rs 1 -randomize_dt 1 -threads 1" > bench_ours.txt
python bench.py ours 6 "-ndims 6 -limit 340 -do_rs 1 -randomize_dt 0 -threads 1" >> bench_ours.txt
python bench.py ours 6 "-ndims 6 -limit 340 -do_rs 0 -randomize_dt 0 -threads 1" >> bench_ours.txt
python bench.py ours 6 "-ndims 7 -limit 340 -do_rs 1 -randomize_dt 1 -threads 1" >> bench_ours.txt
python bench.py ours 6 "-ndims 7 -limit 340 -do_rs 1 -randomize_dt 0 -threads 1" >> bench_ours.txt
python bench.py ours 6 "-ndims 7 -limit 340 -do_rs 0 -randomize_dt 0 -threads 1" >> bench_ours.txt
python bench.py ours 6 "-ndims 8 -limit 340 -do_rs 1 -randomize_dt 1 -threads 1" >> bench_ours.txt
python bench.py ours 6 "-ndims 8 -limit 340 -do_rs 1 -randomize_dt 0 -threads 1" >> bench_ours.txt
python bench.py ours 6 "-ndims 8 -limit 340 -do_rs 0 -randomize_dt 0 -threads 1" >> bench_ours.txt

python bench.py ours 6 "-ndims 6 -limit 340 -do_rs 1 -randomize_dt 1 -threads 8" >> bench_ours.txt
python bench.py ours 6 "-ndims 6 -limit 340 -do_rs 1 -randomize_dt 0 -threads 8" >> bench_ours.txt
python bench.py ours 6 "-ndims 6 -limit 340 -do_rs 0 -randomize_dt 0 -threads 8" >> bench_ours.txt
python bench.py ours 6 "-ndims 7 -limit 340 -do_rs 1 -randomize_dt 1 -threads 8" >> bench_ours.txt
python bench.py ours 6 "-ndims 7 -limit 340 -do_rs 1 -randomize_dt 0 -threads 8" >> bench_ours.txt
python bench.py ours 6 "-ndims 7 -limit 340 -do_rs 0 -randomize_dt 0 -threads 8" >> bench_ours.txt
python bench.py ours 6 "-ndims 8 -limit 340 -do_rs 1 -randomize_dt 1 -threads 8" >> bench_ours.txt
python bench.py ours 6 "-ndims 8 -limit 340 -do_rs 1 -randomize_dt 0 -threads 8" >> bench_ours.txt
python bench.py ours 6 "-ndims 8 -limit 340 -do_rs 0 -randomize_dt 0 -threads 8" >> bench_ours.txt


