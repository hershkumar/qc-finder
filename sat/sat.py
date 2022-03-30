# SAT solver for quantum circuits
# based on http://www.downloadmaghaleh.com/wp-content/uploads/edd/maghaleh/1398/13398.pdf
# another useful source is https://uwspace.uwaterloo.ca/bitstream/handle/10012/14480/Amy_Matthew.pdf?isAllowed=y&sequence=5
# made to work with CNOT and T circuits

from qiskit import *
import numpy as np
from z3 import *
import re
import warnings
# gets rid of some qiskit deprecation warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

def conv(lst):
    return lst.index(1)


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
#TODO: Write the SAT algorithm from the paper here
def find(phase_rep):
	# unpack the phase polynomial representation
	mat, phases = phase_rep
	print(mat)
	print(phases)
	# the number of qubits is given by the size of the array
	num_qubits = len(mat)
	# make a quantum circuit object
	circ = QuantumCircuit(num_qubits)
	K = 1
	while K <= 10:
		# initialize the z3 solver here
		s = Solver()
		
		# now we generate constraints
		matrices = []
		cnots = []
		hs = []

		# generate the variables that we will be using
		for k in range(K):
			# we generate k matrices, and k CNOT operations
			# each matrix is num_qubits by num_qubits
			Ak = []
			for i in range(num_qubits):
				Ak.append([])
				for j in range(num_qubits):
					Ak[i].append(Bool("A^"+str(k)+"_"+str(i)+"_"+str(j))) # add the matrix element variables
			matrices.append(Ak)
			# each cnot has two row vectors, one for the control and one for the target
			qk = []
			tk = []
			for i in range(num_qubits):
				qk.append(Bool("q^"+str(k) + "_"+str(i)))
				tk.append(Bool("t^"+str(k) + "_"+str(i)))
			cnots.append((qk,tk))
			# now we have the auxiliary variables h
			hk = []
			for i in range(num_qubits):
				hk.append(Bool("h^"+str(k)+"_"+str(i)))
			hs.append(hk)
		# now we generate the actual constraints
		# the first constraint is that the initial matrix is the identity matrix
		for i in range(num_qubits):
			s.add(matrices[0][i][i] == True)
			for j in range(num_qubits):
				if i != j:
					s.add(matrices[0][i][j] == False)
		# Then we have the final matrix clause, the final matrix is equal to the matrix in the phase
		# representation
		for i in range(num_qubits):
			for j in range(num_qubits):
				s.add(matrices[-1][i][j] == bool(mat[i][j]))

		# now we need to encode that every function that we care about will at some point show up in the circuit
		for j in phases:
			# convert the function back into a list
			F = list(j)
			for index in range(len(F)):
				F[index] = bool(int(F[index]))
			# we have the value of the row, now we need to ensure that this row shows up somewhere
			# in the matrices
			clause = False
			# loop over all matrices
			for matrix in matrices:
				# loop over all rows
				for row_ind in range(num_qubits):
					# the check we want is whether the row is equal to the function
					row_eql_check = True
					for index in range(num_qubits):
						row_eql_check = And(row_eql_check, matrix[row_ind][index] == F[index])

					# Or all the clauses together, so that at least 1 must be true
					clause = Or(clause, row_eql_check)
			# add the clause to the SAT problem
			s.add(clause == True)

		# now we need the cnot gates to be valid gates
		# for each cnot, in both q and t, only of the terms can be 1
		for k in range(1, K):
			q = cnots[k][0]
			t = cnots[k][1]
			h = hs[k]

			tmpclause = True
			for i in range(num_qubits):
				tmpclause = Or(tmpclause, q[i])
			# now we AND this against 
			tmpclause2 = True
			for i in range(1, num_qubits):
				for j in range(1, num_qubits):
					if (i < j):
						tmpclause2 = And(tmpclause2, Or(Not(q[i]), Not(q[j])))
			s.add(And(tmpclause, tmpclause2))
			
			# now we do the same thing for the t vectors
			tmpclauset = True
			for i in range(num_qubits):
				tmpclauset = Or(tmpclauset, t[i])
			# now we AND this against 
			tmpclauset2 = True
			for i in range(1, num_qubits):
				for j in range(1, num_qubits):
					if (i < j):
						tmpclauset2 = And(tmpclauset2, Or(Not(t[i]), Not(t[j])))
			s.add(And(tmpclauset, tmpclauset2))

			# the cnots need to have a control and target on a different qubit
			targ_clause = True
			for i in range(1, num_qubits):
				targ_clause = And(targ_clause, q[i] != t[i])
			s.add(targ_clause)

			A = matrices[k]
			A_prev = matrices[k-1]

			for i in range(num_qubits):
				for j in range(num_qubits):
					s.add(A[i][j] == Xor(A_prev[i][j], And(t[i], h[j])))
		# now check if the model is satisfiable with K cnots
		if s.check() == sat:
			print("Constraints satisfied, with K = " + str(K))
			break
		else:
			# if its not satisfiable, we add another CNOT and try again
			K += 1
	# now that the constraints are satisfied, we use the solution to generate a circuit
	# get the results of the model
	model = s.model()
	q_k = []
	t_k = []
	A_k = []

	for k in range(K):
		A_k.append([])
		q_k.append([])
		t_k.append([])
		for i in range(num_qubits):
			q_k[k].append(int(bool(model.eval(cnots[k][0][i], model_completion=True))))
			t_k[k].append(int(bool(model.eval(cnots[k][1][i], model_completion=True))))
		# get the matrices as well
		for i in range(num_qubits):
			A_k[k].append([])
			for j in range(num_qubits):
				A_k[k][i].append(int(bool(model.eval(matrices[k][i][j], model_completion=True))))
	print(q_k)
	print(t_k)
	
	# now we actually make the circuit in qiskit
	for k in range(1, K):
		ctrl_index = conv(q_k[k])
		targ_index = conv(t_k[k])
		circ.cx(ctrl_index, targ_index)
	return circ

# make a testing circuit to test the phase rep on
qc = QuantumCircuit(3)
qc.cnot(0,1)
#qc.t(0)
#qc.t(1)
qc.cnot(1,2)
qc.cnot(0,2)

#qc.tdg(2)

print(qc.draw())

pr = circ_to_pr(qc)
print(find(pr))
