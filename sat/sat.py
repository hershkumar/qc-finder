# SAT solver for quantum circuits
# based on http://www.downloadmaghaleh.com/wp-content/uploads/edd/maghaleh/1398/13398.pdf
# made to work with CNOT and T circuits

import qiskit

# we first have a function that converts a given circuit to a phase polynomial representation.
# We can iterate through the circuit gate by gate, and keep track of a couple things
# we have a matrix that represents the action of the circuit so far, and we can update that gate by gate
# to keep track of the coefficients of the phase polynomial, we use a dict
# the dict is keyed by the function that the phase is being applied to, which is a list
# the values will just be a coefficient, which we mod by 8.

def circ_to_pr(circuit):
