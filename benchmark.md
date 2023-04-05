commit 30b53a3780cfe6a6be16611968a4210523fa3b5a
argon_2916.inp in my notebook: i5 2 cores with -g -O0
Starting simulation with 2916 atoms for 1000 steps.
Serial Simulation Done. Run time:    129.456s
1 thread Simulation Done. Run time:    133.827s
2 threads Simulation Done. Run time:    104.084s
4 threads Simulation Done. Run time:     75.539s
8 threads Simulation Done. Run time:     61.835s

-O3
1 thread Simulation Done. Run time:     47.017s
2 threads Simulation Done. Run time:     38.781s
4 threads Simulation Done. Run time:     32.115s
8 threads Simulation Done. Run time:     25.753s

-O3 -ffast-math -fomit-frame-pointer
1 thread Simulation Done. Run time:     54.320s
4 threads Simulation Done. Run time:     31.142s
8 threads Simulation Done. Run time:     27.023s


commit 341f79b44ba760cec110e45ddebcc0e91b5c840d
-O3
Parallel region alla Bussi
1 thread Simulation Done. Run time:     45.972s
2 threads Simulation Done. Run time:     36.478s
4 threads Simulation Done. Run time:     27.247s
8 threads Simulation Done. Run time:     20.165s

1000          60.24508534         523.47323926       -4425.11992013       -3901.64668087
