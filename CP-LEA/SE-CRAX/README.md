# CP-LEA against SE-CRAX64-K128

This is a script for CP-LEA against a self-equivalence white-box CRAX implementation. The SE-CRAX is constructed with block size 64 and key size 128. The encodings are affine self-equivalences.

# Experiment Environment
Ubuntu 20.04.6 LTS

Intel(R) Core(TM) i7-8700K CPU @ 3.70GHz   3.70 GHz

16.0 GB RAM

cmake version 3.16.3

gcc version 9.4.0 (Ubuntu 9.4.0-1ubuntu1~20.04.2)

# Build
```
$ gcc SE_default_white_box_crax.c
```

## Run

```
$ ./a.out
```