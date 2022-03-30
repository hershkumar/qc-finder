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
#TODO: Debug this :(
def find(phase_rep):
	# unpack the phase polynomial representation
	mat, phases = phase_rep
	print(mat)
	print(phases)
	# the number of qubits is given by the size of the array
	num_qubits = len(mat)
	# make a quantum circuit object
	circ = QuantumCircuit(num_qubits)
	K = 2
	while K <= KMAX:
		# initialize the z3 solver here
		s = Solver()

		# generate the variables that we will be using
		matrices = []
		cnots = []
		hs = []

		
		for k in range(0, K + 1):
			# we generate K+1 matrices, and K CNOT operations
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
			for j in range(num_qubits):
				if i != j:
					s.add(matrices[0][i][j] == False)
				else:
					s.add(matrices[0][i][j] == True)
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
			phase_clause = False
			# loop over all matrices
			for matrix in matrices:
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

		# now we need the cnot gates to be valid gates
		# for each cnot, in both q and t, only of the terms can be 1
		for k in range(K):
			q = cnots[k][0]
			t = cnots[k][1]
			h = hs[k]
			A = matrices[k]

			tmpclause = False
			for i in range(len(q)):
				tmpclause = Or(tmpclause, q[i])

			# now we AND this against 
			for i in range(num_qubits):
				for j in range(num_qubits):
					if (i < j):
						tmpclause = And(tmpclause, Or(Not(q[i]), Not(q[j])))
			s.add(tmpclause == True)


			# now we do the same thing for the t vectors
			tmpclauset = False
			for i in range(num_qubits):
				tmpclauset = Or(tmpclauset, t[i])
			# now we AND this against 
			for i in range(num_qubits):
				for j in range(num_qubits):
					if (i < j):
						tmpclauset = And(tmpclauset, Or(Not(t[i]), Not(t[j])))
			s.add(tmpclauset == True)

			# the cnots need to have a control and target on a different qubit
			targ_clause = True
			for i in range(num_qubits):
				targ_clause = And(targ_clause, Xor(q[i], t[i]))
			s.add(targ_clause == True)
			# these encoding map how each matrix turns into the next matrix
			if k != 0:
				A_prev = matrices[k - 1]

				for j in range(num_qubits):
					aux_clause = False
					for i in range(num_qubits):
						aux_clause = Xor(aux_clause, And(A_prev[i][j], q[i]))
					s.add(h[j] == aux_clause)
					
				for i in range(num_qubits):
					for j in range(num_qubits):
						s.add(A[i][j] == Xor(A_prev[i][j], And(t[i], h[j])))
		
		# now check if the model is satisfiable with K cnots
		if s.check() != unsat:
			print("Constraints satisfied, with K = " + str(K))
			break
		else:
			# if its not satisfiable, we add another CNOT and try again
			print("unsat for K=" + str(K))
			K += 1
	if K <= KMAX:
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
		for k in range(K):
			ctrl_index = conv(q_k[k])
			targ_index = conv(t_k[k])
			circ.cx(ctrl_index, targ_index)
	return circ

# make a testing circuit to test the phase rep on
qc = QuantumCircuit(2)
qc.cnot(0,1)
qc.cnot(1,0)
#qc.t(0)
#qc.t(1)
#qc.cnot(1,2)
#qc.tdg(2)

#qc.tdg(2)

print(qc.draw())

pr = circ_to_pr(qc)
found = find(pr)
print(found.draw())
# run it to get the unitary
backend = Aer.get_backend('unitary_simulator')

job = execute(qc, backend)
result = job.result()
orig = result.get_unitary(qc)

fjob = execute(found, backend)
fresult = fjob.result()
new = fresult.get_unitary(found)

print(orig.equiv(new))