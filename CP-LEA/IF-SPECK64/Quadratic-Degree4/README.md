# CP-LEA against IF-SPECK64-K128

This is a script for CP-LEA against an implicit white-box SPECK implementation. The IF-SPECK is constructed with block size 64 and key size 128. The encodings are affine-quadratic self-equivalences. The degree of the implicit function is 4.

Due to size constraints, the source file "white_box_backend_64K128_degree4_non-trivialqua." has not been uploaded. This file can be generated using the [implicit white-box arx framework](https://github.com/ranea/whiteboxarx) with the option "--irf-degree 4".


# Experiment Environment
Ubuntu 20.04.6 LTS

Intel(R) Core(TM) i7-8700K CPU @ 3.70GHz   3.70 GHz

16.0 GB RAM

cmake version 3.16.3

gcc version 9.4.0 (Ubuntu 9.4.0-1ubuntu1~20.04.2)

# Build

```
$ mkdir build
$ cd build
$ cmake ..
$ make
```

## Run

```
$ ./CPLEA
```