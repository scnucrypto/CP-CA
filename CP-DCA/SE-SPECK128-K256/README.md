# CP-DCA against SE-SPECK128-K256

This is a script for CP-DCA against a self-equivalence white-box SPECK implementation. The SE-SPECK is constructed with block size 128 and key size 256. The encodings are affine self-equivalences.

# Experiment Environment
Windows 10

Intel(R) Core(TM) i7-8700K CPU @ 3.70GHz   3.70 GHz

16.0 GB RAM

gcc version 8.2.0 (MinGW.org GCC-8.2.0-3)

# Build

```
$ mkdir build
$ cd build
$ cmake -DBUILD_TARGET=win32 -G "MinGW Makefiles" ..
$ mingw32-make
```

## Run

```
$ .\CPDCA.exe
```