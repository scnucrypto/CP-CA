# CP-DCA against IF-SPECK32-K64

This is a script for CP-DCA against an implicit white-box SPECK implementation. The IF-SPECK is constructed with block size 32 and key size 64. The encodings are affine self-equivalences.

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