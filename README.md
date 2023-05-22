# suchrandom

## to compile:

first you need to <pre><code> make </code></pre> the executable file.

Then run <pre><code> ./main <i>run_settings</i> </code></pre>

where *run_settings* has the following structure ```{flag} {arg}```

| Flag  | Argument range                | Default         | Description                                  |
| ----- | ----------------------------- | --------------- | -------------------------------------------- |
| -alg  | {bfs, dfs, bibfs, sv_1, sv_2} | sv_1            | Choose your favorite algorithm               |
| -qp   | [0, 100]                      | 50              | Query Percentage                             |
| -trc  | Any Int                       | 10              | Test Run Count (#times to run the experiment)|
| -ts   | Any Int                       | 1800            | Set time out in seconds                      |
| -os   | Any Int                       | 1223            | Set Operation Seed                           |
| -qs   | Any Int                       | 2334            | Set Query Seed                               |
| -inp  | Any string                    | sample.txt      | Input file name                              |
| -meta | Any string                    | meta-sample.txt | Meta file name                               |
| -out  | Any string                    | output.txt      | Output file name                             |
| -log  | Any string                    | log.txt         | Log file name                                |


### example:
```./main -alg sv_2 -qp 40 -log log.txt```

> **Note**
> you may need to set stack size to unlimited for dfs; use ```ulimit -s unlimited``` before running the executable.
