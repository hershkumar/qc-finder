# qc-finder
CMSC/PHYS457 Group Project: Quantum Circuit Finding and Optimization

## Roadmap

### Circuit Finder
Use the SAT encodings in [this paper](https://www.informatik.uni-bremen.de/agra/doc/konf/13_africon_sat_encoding_quantum_circuits.pdf) and then pass that to a SAT solver, possibly [Z3 from Microsoft](https://github.com/Z3Prover/z3), or we can make our own basic one if we want.

#### TODO:
- [] Figure out what language we want to write it in (c++/python?)
- [] Define the pipeline (matrix -> SAT encoding -> sat solver -> output -> circuit representation)?
- [] Define matrix file encoding (csv? plaintext?)
- [] Implement SAT encoding
- [] Pick SAT solver (z3 from microsoft has c++ and python bindings)
- [] Parse output and get circuit representation and any other important properties

### Circuit Optimization

