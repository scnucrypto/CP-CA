# CP-DCA against SE-SPECK64-K128

This is a script for CP-DCA against a self-equivalence white-box SPECK implementation. The SE-SPECK is constructed with block size 64 and key size 128. The encodings are affine self-equivalences.

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
$ ./CPDCA
```