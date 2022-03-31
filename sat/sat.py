# SAT solver for quantum circuits
# based on http://www.downloadmaghaleh.com/wp-content/uploads/edd/maghaleh/1398/13398.pdf
# another useful source is https://uwspace.uwaterloo.ca/bitstream/handle/10012/14480/Amy_Matthew.pdf?isAllowed=y&sequence=5
# also useful: https://www.borealisai.com/en/blog/tutorial-9-sat-solvers-i-introduction-and-applications/
# made to work with CNOT and T circuits

from qiskit import *
from qiskit.quantum_info import Statevector
import numpy as np
from z3 import *
import re
import warnings



# gets rid of some qiskit deprecation warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

# function that converts a binary array representation of a qubit to an integer
def conv(lst):
    return lst.index(1)

def insert_gates(circ, coeff, qubit):
	if coeff == 1:
		circ.t(qubit)
	if coeff == 2:
		circ.s(qubit)
	if coeff == 3:
		circ.s(qubit)
		circ.t(qubit)
	if coeff == 4:
		circ.z(qubit)
	if coeff == 5:
		circ.z(qubit)
		circ.t(qubit)
	if coeff == 6:
		circ.s(qubit)
		circ.z(qubit)
	if coeff == 7:
		circ.tdg(qubit)
	return circ

KMAX = 10

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
Uses the Z3 theorem prover from Microsoft to solve SAT clauses that were encoded in the paper linked at the top.
"""
def find(pr):
	# unpack the phase representation into the matrix
	# and the phases
	G, phases = pr
	# the number of qubits is the length of the matrix, since its square
	num_qubits = len(G)
	# starting number of gates
	K = 1
	# while we haven't hit the limit of the number of gates
	while K <= KMAX:
		# create the solver instance
		s = Solver()
		# first we want to set up all the variables and structures that we'll be using
		# we have A^0, A^1, A^2 ... A^K
		A = []
		# A^0 is the identity
		ide = []
		for i in range(num_qubits):
			ide.append([])
			for j in range(num_qubits):
				ide[i].append(False)
		for i in range(num_qubits):
			ide[i][i] = True
		A.append(ide)

		# set up the rest of the matrices
		for k in range(1, K+1): # A^k, where k in [1,K]
			Ak = []
			for i in range(num_qubits):
				Ak.append([])
				for j in range(num_qubits):
					Ak[i].append(Bool("A^"+str(k)+"_"+str(i)+"_"+str(j)))
			A.append(Ak)
		# set up the cnot variables
		# the first index is a dummy index, we should never be accessing it
		q = [[-1]]
		for k in range(1, K + 1):
			qk = [[-1]]
			for i in range(1, num_qubits + 1):
				qk.append(Bool("q^"+str(k)+"_"+str(i)))
			q.append(qk)
		
		t = [[-1]]
		# the first index is a dummy index, we should never be accessing it
		for k in range(1, K + 1):
			tk = [[-1]]
			for i in range(1, num_qubits + 1):
				tk.append(Bool("t^"+str(k)+"_"+str(i)))
			t.append(tk)
		# set up the auxiliary variables
		# same as the rest, the first index is a dummy variable
		h = [[-1]]
		for k in range(1, K+1):
			hk = [[-1]]
			# fill the rest, from 1 to num_qubits
			for i in range(1, num_qubits + 1):
				hk.append(Bool("h^"+str(k)+"_"+str(i)))
			h.append(hk)

		# Now we've set up all the variables, we can start encoding the actual clauses
		# We encode that A^K is G
		for i in range(num_qubits):
			for j in range(num_qubits):
				s.add(A[K][i][j] == bool(G[i][j]))
		# now we need to encode that every function that we care about will at some point show up in the circuit
		for j in phases:
			# convert the function back into a list
			F = list(j)
			for index in range(len(F)):
				F[index] = bool(int(F[index]))
			# we have the value of the row, now we need to ensure that this row shows up somewhere
			# in the matrices
			phase_clause = False
			# loop over all matrices
			for matrix in A:
				# loop over all rows
				for row_ind in range(num_qubits):
					# the check we want is whether the row is equal to the function
					row_eql_check = True
					for index in range(num_qubits):
						row_eql_check = And(row_eql_check, matrix[row_ind][index] == F[index])

					# Or all the phase_clauses together, so that at least 1 must be true
					phase_clause = Or(phase_clause, row_eql_check)
			# add the phase_clause to the SAT problem
			s.add(phase_clause == True)

		# And then we get to the CNOT clauses
		# these make sure that any CNOT we make is a valid cnot
		for k in range(1, K + 1):
			# set up the variables for this iteration
			q_curr  = q[k]
			t_curr = t[k]
			A_curr = A[k]
			h_curr = h[k]
			# this clause ensures that the q array has one and only one 1 among its elements
			q_clause = False
			for i in range(1, num_qubits + 1):
				q_clause = Or(q_clause, q_curr[i])
			for i in range(1, num_qubits + 1):
				for j in range(1, num_qubits + 1):
					if i < j:
						q_clause = And(q_clause, Not(And(q_curr[i], q_curr[j])))
			s.add(q_clause == True)

			# now we do the same clause for the t arrays
			t_clause = False
			for i in range(1, num_qubits + 1):
				t_clause = Or(t_clause, t_curr[i])
			for i in range(1, num_qubits + 1):
				for j in range(1, num_qubits + 1):
					if i < j:
						t_clause = And(t_clause, Not(And( t_curr[i],t_curr[j])))
			s.add(t_clause == True)

			# now make sure that the control and target are different qubits
			cnot_clause = True
			for i in range(1, num_qubits + 1):
				cnot_clause = And(cnot_clause, q_curr[i] == t_curr[i])
			s.add(cnot_clause == False)

			# Now define the auxiliary variables
			for j in range(1, num_qubits + 1):
				h_clause = False
				for i in range(1, num_qubits + 1):
					h_clause = Xor(h_clause, And(A[k - 1][i - 1][j - 1], q_curr[i]))
				h_curr[j] = h_clause
		
			# and then how the matrices relate to the CNOTS and the auxiliary variables
			# this clause governs how we get from A^k-1 to A^k, its a sort of matrix evolution clause
			for i in range(1, num_qubits + 1):
				for j in range(1, num_qubits + 1):
					s.add(A[k][i - 1][j - 1] == Xor(A[k-1][i - 1][j - 1], And(t_curr[i], h_curr[j])))
		# now we check if the model is satisfiable with K CNOTS
		if s.check() == sat:
			print("sat for K=" + str(K))
			break
		else:
			# if its not, we try it with 1 more CNOT
			print("unsat for K=" + str(k))
			K += 1
	# now we make the circuit
	if K <= KMAX:
		# now that the constraints are satisfied, we use the solution to generate a circuit
		# get the results of the model
		model = s.model()
		# make a quantum circuit object
		circ = QuantumCircuit(num_qubits)
		
		model_A = []
		model_q = []
		model_t = []
		# fill the lists with dummy variables to start
		ide = []
		for i in range(num_qubits):
			ide.append([])
			for j in range(num_qubits):
				ide[i].append(0)
		for i in range(num_qubits):
			ide[i][i] = 1
		model_A.append(ide)
		# fill the matrices
		# this is actually not really necessary for building the circuit, so I've commented it out for now

		# for k in range(1, K+1):
		# 	Ak = []
		# 	for i in range(num_qubits):
		# 		Ak.append([])
		# 		for j in range(num_qubits):
		# 			Ak[i].append(int(bool(model.eval(A[k][i][j], model_completion=True))))
		# 	model_A.append(Ak)

		# the first index is a dummy index, we should never be accessing it
		
		for k in range(1, K + 1):
			qk = []
			for i in range(1, num_qubits + 1):
				qk.append(int(bool(model.eval(q[k][i], model_completion=True))))
			model_q.append(qk)

		# the first index is a dummy index, we should never be accessing it
		
		for k in range(1, K + 1):
			tk = []
			for i in range(1, num_qubits + 1):
				tk.append(int(bool(model.eval(t[k][i], model_completion=True))))
			model_t.append(tk)

		# start state of the circuit
		start = np.eye(num_qubits, dtype=np.int8)
		# now we actually make the circuit in qiskit
		# get the circuit data
		qc = circ.data

		# store the T instructions for later
		t_inst = []
		for i in range(num_qubits):
			circ.t(i)
			t_inst.append(qc[0])
			del qc[0]
		# loop through all the CNOTs we have
		for k in range(0, K):
			ctrl_index = conv(model_q[k])
			targ_index = conv(model_t[k])
			# attach the cnot gate to the circuit
			circ.cx(ctrl_index, targ_index)
			
			# now update the matrix
			for j in range(num_qubits):
				if start[ctrl_index][j] == 1: # if the matrix element is 1
					start[targ_index][j] = 1 - start[targ_index][j] # flip the bit in the target row
			# now we check if this matches any of the phases

			for key in phases:
				# for every row in the matrix
				for row_ind in range(num_qubits):
					equal = True
					row = start[row_ind] # get the row of the matrix
					# check element by element and see if they are equal
					for col in range(num_qubits):
						if row[col] != int(key[col]):
							equal = False
					# if they are equal, we add the phase gates
					if equal:
						# number of times we have to add a T gate
						coeff = phases[key]
						# get the qubit number to add the gate to
						qubit = row_ind
						# depending on the coefficient, we insert the correct gate(s)
						circ = insert_gates(circ, coeff, qubit)
						# after we finish adding it, we decrement the pair down to 0
						# this way we don't repeat adding the phase gates
						phases[key] = 0
		# return the circuit object
		return circ
"""
This function compares the unitary matrices of the two circuits, and takes into account global phase.
e^iphi A != A
"""
def comp_global(c1,c2):
	backend = Aer.get_backend('unitary_simulator')

	job = execute(c1, backend)
	result = job.result()
	orig = result.get_unitary(c1)

	fjob = execute(c2, backend)
	fresult = fjob.result()
	new = fresult.get_unitary(c2)

	return orig.equiv(new)
"""
This function compares the statevectors of the two circuits, and does not take into account the global phase.
e^iphi A == A
"""
def comp(c1,c2):
	return Statevector.from_instruction(c1).equiv(Statevector.from_instruction(c2))

# make a testing circuit to test the phase rep on
qc = QuantumCircuit(3)
# qc.cnot(0,1)
# qc.t(0)
# qc.cnot(1,2)
# qc.t(1)
# qc.s(2)
# qc.cnot(2,1)
# qc.s(2)
# qc.cnot(0,2)
# qc.cnot(0,1)
# qc.t(1)
qc.cnot(0,1)
qc.t(0)
qc.cnot(1,0)
qc.cnot(0,1)
qc.t(1)
qc.cnot(1,2)
qc.tdg(2)
print(qc.draw())

pr = circ_to_pr(qc)

found = find(pr)

print(found.draw())
print(comp(qc, found))