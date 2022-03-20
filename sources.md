Here's a list of sources that I think might be worth taking a look at when writing a SAT solver encoding.

### Towards a SAT Encoding for Quantum Circuits: A Journey From Classical Circuits to Clifford Circuits and Beyond
[Link to paper.](https://arxiv.org/abs/2203.00698)
They define an encoding that can be used for arbitrary quantum circuits, but then talk about how this is practically infeasible in general. However, they do show that it is doable for Clifford gates.


### SAT-based {CNOT, T} quantum circuit synthesis
[Link to paper.](http://www.downloadmaghaleh.com/wp-content/uploads/edd/maghaleh/1398/13398.pdf)
This looks at the restricted gate set of just CNOTs and T gates, and define a SAT encoding that attempts to minimize the number of CNOT gates in a CNOT, T circuit. This might be a good starter SAT encoding to work with, because it also attempts to do circuit optmization as part of the SAT process.
