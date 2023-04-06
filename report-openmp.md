# Report: OpenMP Paralellization of Molecular Dynamics Code

## Fernando Santana Pacheco

## Serial version (baseline code), compiled with -O3

The original code has no paralellization nor optimization. Compiled with gcc 8.3.0  and `-O3`, the reference timings in Ulysses were

Commit 8b4306d392cd5bd9c9261c7395286a5b51869de9. Raw timings on file omp_md.o9148952.

| Num threads     | Time (s) 108atoms | 2916atoms |
| --------------- | ----------------- | --------- |
| 1               | 3.475             | 121.887   |

## Serial version with optimized forces, compiled with -O3

After some initial tests, I decided to implement directly the OpenMP paralellization over a first force optimized version from my colleague. The timings are the following:

Commit c3cba54507dec408ff548cd6ff2f2235a935b164. Raw timings on file omp_md.o9147913.

Includes Newton's 3rd law and avoids power and sqrt inside the loop

| Num threads     | Time (s) 108atoms | 2916atoms |
| --------------- | ----------------- | --------- |
| 1               | 0.137             | 43.643    |

## OpemMP with reduction/critical

A first attempt with OpenMP was the use of a parallel for with reduction for potential energy followed by a critical region, as in the lecture notes (page 23). As stated in the notes, the results were not really faster when increasing the number of cores.

Commit c3cba54507dec408ff548cd6ff2f2235a935b164 (OpenMP version compiled with -O3) (reduction/critical). Raw timings on file omp_md.o9147876

Allocated 10 cpus on Ulysses

| Num threads     | Time (s) 108atoms | 2916atoms |
| --------------- | ----------------- | --------- |
| 1               | 0.174             | 50.889    |
| 2               | 0.382             | 41.398    |
| 5               | 0.744             | 36.887    |
| 10              | 1.399             | 50.857    |
| 20              | 0.879             | 61.665    |

## OpenMP with private array for each thread

Then we moved to another strategy: each thread with a private array, that is, without any data race while writing the forces of the inner loop (`j`), followed by a critical region to reduce each thread computation to the full array.

Commit e35a0b05b2ed2324aaf9fbf2615c28f6fc9512ae (Parallel region) (with -O3). Raw timings on file omp_md.o9147871.

Allocated 10 cpus on Ulysses

| Num threads     | Time (s) 108atoms | 2916atoms |
| --------------- | ----------------- | --------- |
| 1               | 0.194             | 42.235    |
| 2               | 0.155             | 32.965    |
| 5               | 0.083             | 18.094    |
| 10              | 0.062             | 11.936    |
| 20              | 0.129             | 9.595     |

One thing that can be improved is the memory allocation. In the current implementation it is done at each call of the forces function.


