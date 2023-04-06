# Results

## Initial tests on my notebook

### commit 30b53a3780cfe6a6be16611968a4210523fa3b5a (reduction/critical)

argon_2916.inp in my notebook: i5 2 cores with -g **-O0**

```
Starting simulation with 2916 atoms for 1000 steps.
Serial Simulation Done. Run time:    129.456s
1 thread Simulation Done. Run time:    133.827s
2 threads Simulation Done. Run time:    104.084s
4 threads Simulation Done. Run time:     75.539s
8 threads Simulation Done. Run time:     61.835s
```

### Same commit compiled with **-O3**

```
1 thread Simulation Done. Run time:     47.017s
2 threads Simulation Done. Run time:     38.781s
4 threads Simulation Done. Run time:     32.115s
8 threads Simulation Done. Run time:     25.753s
```

### Same commit compiled with **-O3 -ffast-math -fomit-frame-pointer**

```
1 thread Simulation Done. Run time:     54.320s
4 threads Simulation Done. Run time:     31.142s
8 threads Simulation Done. Run time:     27.023s
```

### commit 341f79b44ba760cec110e45ddebcc0e91b5c840d (parallel region) compiled with **-O3**

```
1 thread Simulation Done. Run time:     45.972s
2 threads Simulation Done. Run time:     36.478s
4 threads Simulation Done. Run time:     27.247s
8 threads Simulation Done. Run time:     20.165s
```

## Running on Ulysses

At `/home/fsantana/ljmd`

## Steps to compile

```
module load cmake
module load intel/2021.1
cmake -D CMAKE_C_COMPILER=gcc ..
make -j2 VERBOSE=1
```

### commit e35a0b05b2ed2324aaf9fbf2615c28f6fc9512ae (Parallel region) (with -O3)

File omp_md.o9147871

Allocated 10 cpus

| Num threads     | Time (s) 108atoms | 2916atoms |
| --------------- | ----------------- | --------- |
|20               | 0.129             | 9.595
|10               | 0.062             | 11.936
|5                | 0.083             | 18.094
|2                | 0.155             | 32.965
|1                | 0.194             | 42.235

### commit c3cba54507dec408ff548cd6ff2f2235a935b164 (OpenMP version compiled with -O3) (reduction/critical)

File omp_md.o9147876

Allocated 10 cpus

| Num threads     | Time (s) 108atoms | 2916atoms |
| --------------- | ----------------- | --------- |
| 20              | 0.879             | 61.665    |
| 10              | 1.399             | 50.857    |
| 5               | 0.744             | 36.887    |
| 2               | 0.382             | 41.398    |
| 1               | 0.174             | 50.889    |

### The serial version from the same commit, but compiled with -O3

File omp_md.o9147913

| Num threads     | Time (s) 108atoms | 2916atoms |
| --------------- | ----------------- | --------- |
| 1               | 0.137             | 43.643

### The serial version from the same commit, but compiled with -O0

File omp_md.o9147922

| Num threads     | Time (s) 108atoms | 2916atoms |
| --------------- | ----------------- | --------- |
| 1               | 0.331             | 118.622   |

