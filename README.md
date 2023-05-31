# suchrandom

## To Compile

First, you need to <pre><code> make </code></pre> the executable file.

Then run <pre><code> ./main <i>run_settings</i> </code></pre>

where *run_settings* has the following structure ```{flag} {arg}```

| Flag  | Argument range                | Default         | Description                                  |
| ----- | ----------------------------- | --------------- | -------------------------------              |
| -alg  | {bfs, dfs, bibfs, sv_*k*}     | sv_1            | Choose your favorite algorithm.              |
| -qp   | [0, 100]                      | 50              | Query Percentage                             |
| -qt   | Any Int                       | 0               | Query Timestamp                              |
| -trc  | Any Int                       | 10              | Test Run Count (#times to run the experiment)|
| -ts   | Any Int                       | 1800            | Set time out in seconds                      |
| -os   | Any Int                       | 1223            | Set Operation Seed                           |
| -qs   | Any Int                       | 2334            | Set Query Seed                               |
| -inp  | Any string                    | sample.txt      | Input file name                              |
| -meta | Any string                    | meta-sample.txt | Meta file name                               |
| -out  | Any string                    | output.txt      | Output file name                             |
| -log  | Any string                    | log.txt         | Log file name                                |

### Example

```./main -alg sv_2 -qp 40 -log log.txt```

> **Note**
> You may need to set stack size to unlimited for dfs; use
```ulimit -s unlimited``` before running the executable.

## Samples

For now, all samples are available [here](https://dyreach.taa.univie.ac.at/).
You can download and extract them and give the path as an input (with
 flag *-inp*) to the executable.

> **Note**
> Since the operations for our project are insertions-only (in contrast to
 the fully dynamic case), we have to generate our own input. This new input is stored
 in a *build* folder in the same directory as the original input.

## Inputs

Our inputs consist of the following line(s) ```u v t``` where *u* and *v* are
vertices for our edge insertions (sequentially), and *t* is the timestamp.
The timestamps are non-decreasing and by setting *-qt*, you can specify the timestamp
for our initial graph (i.e. edge insertions before *qt* make up the initial graph).
After we have generated our initial graph, by *-qp* chance, there will be a
query operation after each insertion. The query operations are generated randomly
over all the vertices.

<!-- add log description later-->