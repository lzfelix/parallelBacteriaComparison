# Assignment

## Initial timming


| Tidy     | Improved | Improved w/o IO |
|----------|----------|-----------------|
| 36 mins? | ~48 s    | ~45 s           |


These results are single run times using time() functionality and **-O3** flag.

## Code versions

- **v1.0:** is the same code as improved, but with the bacteria class extracted;
- **v1.1:** the creation of bacterias and calculation of its vectors were
parallelised on the same thread. This version also replaced the timer by high
resolution timers

## v1.1 Runs

This version improved a little the amount of time needed, at cost of using heaps of RAM. It divdes the workload evenly among $k$ threads, responsible for creating the stochastic vectors for the bacterias. Some run results are shown below:


|   #Threads   |   RAM (GB)   |   Processing   |   Time(ms)   |
|--------------|--------------|----------------|--------------|
|    10        |    ~ 8.5     |     310        |  43,728.262  |
|    08        |      7.0     |     357        |  34,988.207  |
|    06        |      5.9     |     380        |  27,946.863  |
|    04        |      3.5     |     350        |**25,488.096**|
|    03        |      2.5     |     270        |  27,628.929  |
|    02        |      1.2     |     195        |  29,800.941  |
|    01        |      1.7     |     099        |  44,862.351  |

**OBS:**
- The RAM and Processing are peak values informed from xCode's built in profiler
- Processing (%) is informed by xCode's built in profiler and varies between 0
(no cores in use) and 400 (4 cores being completely used). Since the processor is
shared between other processes, even with most of programs closed. The OS schedules
other background processes to work, that's why it's virtually impossible to reach
400% of CPU usage. **Absolute values may be obtained by using Activity profiler**.

## v1.2 Runs

This version interleaves Bacteria creation and stochastic processing, since as displayed on (one of the profiler figures) this activity also takes considerable time.

Such approach wasn't very successful, as it increased the time to compute the similarity
between the 41 bacterias. This is probably due to the fact that the threads responsible for reading and doing the calculation must be synchronised.

|   #Threads   |   RAM (GB)   |   Processing   |   Time(ms)   |
|--------------|--------------|----------------|--------------|
|    10        |     11,0     |     300        |  53,838.934  |
|    08        |     11.0     |     300        |  41,007.902  |
|    06        |      8.8     |     365        |  38,021.304  |
|    04        |      5.0     |     350        |  27,534.663  |
|    03        |      4.0     |     327        |**26,914.542**|
|    02        |      2.5     |     247        |  28,731.686  |
|    01        |      1.9     |     125        |  40,923.022  |

Notice that #$k$ threads mean $k$ threads for reading files and $k$ threads for computation, interleaved.

A probable cause of slowdown is the significant increase on amount of used RAM, as the Bacteria constructor allocate some big arrays (values here), which are dealocated only after the stochastic precomputation is done. This probably increased the amount of page faults and virtual memory pagination on the computer, which leads to more disk accesses thus more time spent.

Previously to this strategy, the temporary arrays were allocated, the pre-computation was done and these arrays were removed from memory, decreasing the amount of actually used RAM.

The approach was the following: the reader thread is initialised and starts parsing each file on
its work block. Meanwhile the analysing thread is started and waits for a signal from its coworker. When the first file is read, the working thread wakes up the analyser and moves to the next bacteria.

When awaken, the analyser calculates the vectors and moves to the next bacteria of its block. If this other file was already parsed by worker, then it continues, otherwise it waits a signal meaning that this new bacteria is ready.

An alternative approach was decreasing the amount of syncrhonisations needed by sending signals only at set intervals of time. Such interval was determined by dividing the worker's thread workload into stages. Uppon completion of each stage, a signal is sent to the analyser. Even though this approach wasn't successful.

***Runs varying delta on 4 threads***

|   $\delta$   |   Time(ms)   |
|--------------|--------------|
| $\frac{1}{workload}$     |  29,534.261  |
|    1/8       |  32,969.799  |
|    1/4       |  32,027.620  |
|    1/2       |  42,774.848  |
|    1/1       |  47,706.509  |

Notice that $\delta = 1/2$ consumed 12GB of ram and $\delta = 1/1$ consumed 20 GB of RAM due the previously discussed aspects.

*Delta indeed reduced the idle time on the CPUs, but probably due to extra amount of RAM usage, the system faces extra overhead, slowing down the process due pagination.*

**Idea:** instead of keeping $k$ threads, keep just one or two, that read the first bacteria from each thread, then the second, and so on...

**This was done by using *pthread's* *signals*.**

## v1.3 Parallelsing the comparison

As it is possible to see on the instrumentation's CPU use graphs, the comparison part of code consists in a long flat valley, which could be parallelised, in order to increse the CPU use and to decrease the time spent on this processing.

This part consists in comparing each bacteria $b_i$ with all ther other $b_j$, where $i \neq j$. Since the $C(b_i,b_j) = C(b_j,b_i)$, there's no need to perform two comparisons. If such values are stored in a bidimensional matrix, only the elements above the main diagonal have to be compared, leading to the following iteration pattern:

**<center>Matrix figure</center>**

Parallelising the iterations over this matrix by distributing one row or column to each thread leads to an imbalancing, since the thread responsible for processing the first row does much more processing than the thread responsible for the last one.

In order to prevent such scenario, this workload was rearranged in a one-dimensional array and then split evenly across $k$ threads.

In order to do so, first it's needed to find out how many elements need to be effectively compared, which are only the elements above the main diagonal. To find this value, consider the matrix as a continuos square. It is possible to find the square area by

$$A = d \times d $$

where $d$ is the dimension of the (square) matrix.

Since just the upper elements are computated, these can be estimated as half of the square area, hence

$$A' = \frac{1}{2}A $$

Unfortunately the amount of elements in the matrix is discrete, instead of continuous. To adjust the obtained value to this fact, it's possible to consider that since the matrix is squared, it's diagonal intercepts each of the diagonal elements in half, making the area of two small squares equals the area of a cell square. It happens, however, that the area below the trace only on the main diagonal is at most the matrix dimension $d$.

Putting these elements togheter, it's possible then to determine the amount of elements below such trace, as follows:

$$\bar A = \frac{1}{2}d^2 +  d = \frac{d}{2}(d + 1)$$

Knowing this information, it's easy to find how many elements are above the main diagonal (thus need to be computed) by subtracting the "matrix area" from $\bar A$. Namely:

$$\bar B = d^2 - \Big[ \frac{d}{2}(d + 1) \Big] = d^2 - \frac{d^2}{2} + d$$
$$ \therefore \bar B(d) = \frac{d}{2}(d - 1)$$

With this information, it's just a matter of dividing the workload evenly among the threads. The problem, however, is that the sequence that needs to be computed is not evenly distributed on the matrix, as Figure 2 shows

**<center>Teatris piece</center>**

By analysing the previous figure it's possible to notice that the last element on each row increases by $d - t$, where $t$ is the number of the row. It is possible to verify that by generating the corresponding sequence, on this case, $(6, 12, 17, 21, 24, 26, 27)$.

Such sequence can be seen as the *upper bound* of each row. The *lower bound* can be found then just by subtracting 1 from the previous upper limit, except for the first row, which value is 0.

By keeping an array with these information in memory, it's possible to easily find the next bacterias to be compared and where the resulting value has to be stored.

*talk about pseudo-binary-seach*

|   #Threads   |   RAM (GB)   |   Processing   |   Time(ms)   |
|--------------|--------------|----------------|--------------|
|    10        |     9.2      |     355        |  37,932.374  |
|    08        |     7.8      |     358        |  33,442.493  |
|    06        |     5.5      |     350        |  31,046.430  |
|    05        |     4.0      |     360        |**20,257.896**|
|    04        |     4.0      |     360        |  20,356.844  |
|    03        |     3.0      |     276        |  23,069.363  |
|    02        |     1.6      |     196        |  25,454.075  |
|    01        |     1.8      |     99         |  45,824.223  |
