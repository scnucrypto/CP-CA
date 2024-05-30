# CP-LEA against IF-SPECK64-K128

This is a script for CP-LEA against an implicit white-box SPECK implementation. The IF-SPECK is constructed with block size 64 and key size 128. The encodings are affine-quadratic self-equivalences. The degree of the implicit function is 4.

Due to size constraints, the source file "white_box_backend_64K128_degree4_non-trivialqua.c" (691 MB) of the target white-box implementation has not been uploaded. This file can be generated using the [implicit white-box arx framework](https://github.com/ranea/whiteboxarx) with the option "--irf-degree 4".

Detailed generation:

```
$ ./sage -python speck.py --key 0f0e0d0c 0b0a0908 07060504 03020100 --block-size 64 --output-file speck64_128_affine_layers.sobj
$ ./sage -python generate_wb.py --input-file speck64_128_affine_layers.sobj --irf-degree 4 --trivial-external-encodings  --disable-redundant-perturbations --output-file speck64_128_irf.sobj 
$ ./sage -python export_wb.py --input-file speck64_128_irf.sobj --irf-degree 4 --output-file white_box_backend_64K128_degree4_non-trivialqua.c --first-explicit-round "x = ((x >> 7) | (x << (WORD_SIZE - 7))); x = (x + y) & WORD_MASK;" --disabled-redundant-perturbations --cancel-external-encodings
```

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