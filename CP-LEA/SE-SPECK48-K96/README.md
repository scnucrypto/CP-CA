# CP-LEA against SE-SPECK48-K96

This is a script for CP-LEA against a self-equivalence white-box SPECK implementation. The SE-SPECK is constructed with block size 48 and key size 96. The encodings are affine self-equivalences.

# Experiment Environment
Ubuntu 20.04.6 LTS

Intel(R) Core(TM) i7-8700K CPU @ 3.70GHz   3.70 GHz

16.0 GB RAM

cmake version 3.16.3

gcc version 9.4.0 (Ubuntu 9.4.0-1ubuntu1~20.04.2)

# Build

```
$ gcc default_white_box_speck_48_K96.c
```

## Run

```
$ ./a.out
```