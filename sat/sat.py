# SAT solver for quantum circuits
# based on http://www.downloadmaghaleh.com/wp-content/uploads/edd/maghaleh/1398/13398.pdf
# another useful source is https://uwspace.uwaterloo.ca/bitstream/handle/10012/14480/Amy_Matthew.pdf?isAllowed=y&sequence=5
# made to work with CNOT and T circuits

from qiskit import *
import numpy as np
import z3
import warnings
# gets rid of some qiskit deprecation warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

"""
This function takes in a circuit composed of CNOTS, T, T^2, T^3,... T^7, and returns a circuit with
just CNOTS and Ts. This essentially expands the exponentials of the T gate into just repeated T gates.
This can then be passed into the function that finds the phase polynomial representation of a circuit.
"""
def expand_t(circuit):
	num_qubits = circuit.num_qubits
	# access the circuits instruction data
	qc = circuit.data
	
	# loop through all the instructions in the circuit
	for i in range(len(qc)):
		# get the instruction
		instruction = qc[i]
		# get the name of the current instruction
		name = instruction[0].name
		if name == "s": # T^2
			# change the instruction name to "t"
			instruction[0].name = "t"
			# get the instruction
			copy = instruction
			qc.insert(i, copy)
		if name == "z": # T^4
			# change the instruction name to "t"
			instruction[0].name = "t"
			# get the instruction
			copy = instruction
			for k in range(3):
				qc.insert(i, copy)
		if name == "tdg": # T^7
			# change the instruction name to "t"
			instruction[0].name = "t"
			# get the instruction
			copy = instruction
			for k in range(6):
				qc.insert(i, copy)
	return circuit


"""
we have a function that converts a given circuit to a phase polynomial representation.
We can iterate through the circuit gate by gate, and keep track of a couple things
we have a matrix that represents the action of the circuit so far, and we can update that gate by gate
to keep track of the coefficients of the phase polynomial, we use a dict
the dict is keyed by the function that the phase is being applied to, which is a list
the values will just be a coefficient, which we mod by 8.
"""

"""
Function that returns the phase polynomial representation of a passed in qiskit circuit.
"""
def circ_to_pr(circuit):
	# first expand all T gates in the circuit:
	circuit = expand_t(circuit)
	# we begin by getting the number of qubits
	num_qubits = circuit.num_qubits
	# access the circuits instruction data
	qc = circuit.data
	# we start the overall matrix at the identity
	# note that we store it as an int array, to make things easier
	# nothing in this matrix will every be anything other than a 1 or a 0
	mat = np.eye(num_qubits, dtype=np.int8)
	# and we need a dictionary to store the phases
	phases = {}
	# loop through every instruction in the circuit
	for i in range(len(qc)):
		# get the instruction
		instruction = qc[i]
		# get the name of the current instruction
		name = instruction[0].name

		if name == "cx": # if the instruction is a CNOT
			# get the indices of the qubits that the cnot acts on
			control = instruction[1][0].index
			target = instruction[1][1].index
			# now we modify the matrix accordingly
			# we take the elements of the control'th row, and if its 1, we flip the corresponding
			# matrix element in the target row
			for j in range(num_qubits):
				if mat[control][j] == 1: # if the matrix element is 1
					mat[target][j] = 1 - mat[target][j] # flip the bit in the target row
		elif name == "t": # if its a T gate
			# get the qubit that its acting on
			qbit = instruction[1][0].index
			# get the current state the qubit is in, via the matrix
			f = mat[qbit]
			# convert the array to a string to index by in the dictionary
			fstring = "".join(map(str, f))
			# either insert or increment the coefficient in the dictionary
			if fstring not in phases:
				phases[fstring] = 1
			else:
				phases[fstring] += 1
				# mod the phase by 8
				phases[fstring] = phases[fstring] % 8
	# return the phase polynomial
	return (mat, phases)
"""
Takes in the phase polynomial representation of a CNOT,T circuit and returns a circuit with a 
minimal number of CNOT gates.
"""
def find(phase_rep):
	# unpack the phase polynomial representation
	mat, phases = phase_rep
	# the number of qubits is given by the size of the array
	num_qubits = len(mat)
	# make a quantum circuit object
	circ = QuantumCircuit(num_qubits)
	#TODO: Write the SAT algorithm from the paper here

	return circ

# make a testing circuit to test the phase rep on
qc = QuantumCircuit(3)
qc.cnot(0,1)
qc.t(1)
qc.cnot(1,2)

qc.tdg(2)

print(qc.draw())

print(circ_to_pr(qc))

