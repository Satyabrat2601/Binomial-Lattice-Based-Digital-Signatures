# BIRDS: BInomial Rejection-based Digital Signature in Lattices

This repository contains the reference implementation of the BIRDS signature scheme and a security estimator. The scheme provides three parameter sets targeting ~120-bit, ~180-bit, and ~260-bit security.

## Build Instructions

The implementation contains correctness tests and benchmarking programs with a Makefile to facilitate compilation.

### Prerequisites

The implementation has no external library dependencies. Only a C11 compiler and GNU Make are required.

```sh
# Ubuntu/Debian
sudo apt install gcc make

# macOS
xcode-select --install
```

### Test Programs

To compile and run the test programs, go to the `BIRDS_SCHEME/` directory and run

```sh
make
```

This builds and tests all three parameter sets, producing the executables

```sh
build/mode1/test_mode1
build/mode3/test_mode3
build/mode5/test_mode5
```

where `mode1`, `mode3`, and `mode5` correspond to BIRDS-I (~120-bit), BIRDS-II (~180-bit), and BIRDS-III (~260-bit) respectively.

To build and test a single mode:

```sh
make test MODE=1
make test MODE=3
make test MODE=5
```


### Key and Signature Sizes

|                  | BIRDS-I | BIRDS-II| BIRDS-III |
|------------------|--------|---------|----------|
| Public key       | 544 B  | 800 B   | 1056 B   |
| Signature        | 773 B  | 1001 B  | 1408 B   |

## Security Estimator

The `SECURITY_ESTIMATOR/` directory contains the Python security estimator. It requires `numpy`:

```sh
pip install numpy
```

To run:

```sh
cd SECURITY_ESTIMATOR/SECURITY_ESTIMATOR
python3 BIRDS_main.py
```

This prints the full MSIS and MLWE security analysis (classical and quantum bits), rejection rates, and signature sizes for all three parameter sets.

## Project Structure

```
BIRDS_SCHEME/BIRDS_SCHEME/
├── src/            # Core implementation (sign, poly, NTT, encoding, packing, …....)
├── include/        # Headers — (params.h , api.h , config.h .....)
├── test/           # Correctness test program (main.c)
├── benchmark/      # Benchmarking program (speed.c, speed.h .....)
├── build/          # Compiled binaries (created by make)
└── Makefile

SECURITY_ESTIMATOR/SECURITY_ESTIMATOR/
├── BIRDS_main.py        # Main entry point
├── MSIS_security.py    # MSIS hardness estimator
├── MLWE_security.py    # MLWE hardness estimator
├── model_BKZ.py        # BKZ cost model
└── proba_util.py       # Probability utilities for CBD analysis
```
