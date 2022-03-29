# SAT solver for quantum circuits
# based on http://www.downloadmaghaleh.com/wp-content/uploads/edd/maghaleh/1398/13398.pdf
# made to work with CNOT and T circuits

from qiskit import *
import numpy as np

# we first have a function that converts a given circuit to a phase polynomial representation.
# We can iterate through the circuit gate by gate, and keep track of a couple things
# we have a matrix that represents the action of the circuit so far, and we can update that gate by gate
# to keep track of the coefficients of the phase polynomial, we use a dict
# the dict is keyed by the function that the phase is being applied to, which is a list
# the values will just be a coefficient, which we mod by 8.

def circ_to_pr(circuit):
	# we first begin by getting the number of qubits
	num_qubits = circuit.num_qubits
	# access the circuits instruction data
	qc = circuit.data
	# 
	# loop through every instruction in the circuit
	for i in range(len(qc)):
		# to get the matrix, we only care about cnot gates, we ignore the T gates
		# first we check what the instruction is
		name = qc[i][0].name
		if name == "cx":
			# get the target qubits
		print(qc[i][1])
	


# make a testing circuit to test the phase rep on
qc = QuantumCircuit(3)
qc.cnot(0,1)
qc.cnot(1,2)
circ_to_pr(qc)
print(qc.draw())