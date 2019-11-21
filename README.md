# DSEP
This repo is the implementation and measurement of "DSEP-Dynamic sparsity and expansion packets for Fulcrum network coding" paper.
We have implemented three variations of sparse Fulcrum codes using the Kodo library (kodo-fulcrum version 7.0) and have improved the measuremnt of the encoding and decoding throughput as well as the decoding probability with the standard benchmarks available in the library.

How to build and run this benchmark

1) download the kodo-fulcrum version 7.0 at here: https://github.com/steinwurf
and build kodo-fulcrum at here http://docs.steinwurf.com/kodo-fulcrum/master/index.html
(make sure that you configure and build successfully kodo-fulcrum library in your pc, and then carry out next steps)
2) copy both folder src/kodo_fulcrum_sparse and src/tunable_fulcrum to src folder which is downloaded in step 1
3) copy benchmark/throughput_fulcrum to the same folder name which is downloaded in step 1
4) add a command: bld.recurse('benchmark/throughput_fulcrum/') to the file wscript in step 1.
5) rebuild the kodo-fulcrum in your pc. After building successfully the library, you can run the benchmark.



