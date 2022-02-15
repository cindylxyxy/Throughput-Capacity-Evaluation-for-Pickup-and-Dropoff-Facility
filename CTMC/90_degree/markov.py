########################################################################################################################################################
# Last checked: Feb 06, 2022 by Xinyu Liu
# This file runs the Markov Chain model, generates all states, and calculates the stationary distribution.
 
########################################################################################################################################################
import sys
import time
import pprint
import random
import json
import pickle
from copy import deepcopy

from scipy.sparse import csr_matrix, diags, save_npz, load_npz
from scipy.linalg import solve, svd
from scipy.sparse.linalg import spsolve, cg, inv, eigsh, bicg, bicgstab
import numpy as np
from numpy.linalg import eig, eigh, eigvals, eigvalsh

from inputSetDef import allspot, spot2blk
from transitionDef import stateObject
from SLinkedList import SLinkedList
from params import *

###############################################################################
############################# generate all states #############################
###############################################################################
def generate_all(N, rmu, rlambda):

	# create initial state
	start = stateObject(N)
	state_space = dict()
	state_space[start] = 0 
	state_count = 1
	toVisit = SLinkedList()
	toVisit.addHead(start)

	matData = np.array([]) 
	rowIdx = np.array([])
	colIdx = np.array([])

	curr = toVisit.head
	while curr:
		total = 0.0
		curr_idx = state_space[curr.data]
		if curr_idx in rowIdx:
			import pdb; pdb.set_trace()

		for j in range(1, curr.data.N + 1):
			if curr.data.x[j-1] == 1:
				nex = deepcopy(curr.data)
				nex.mu_transition(j)
				matData = np.append(matData, rmu)
				total -= rmu
				rowIdx = np.append(rowIdx, curr_idx)
				if nex not in state_space:
					state_space[nex] = state_count
					toVisit.addEnd(nex)
					state_count += 1
				nex_idx = state_space[nex]
				colIdx = np.append(colIdx, nex_idx)

		if (np.sum(curr.data.y) > 0) | ( (curr.data.x != np.ones(curr.data.N)).any() ):
			nex = deepcopy(curr.data)
			if control == 'full':
				nex.lambda_transition()
			else:
				assert control == 'partial'
				nex.lambda_transition_partial()
			matData = np.append(matData, rlambda)
			total -= rlambda
			rowIdx = np.append(rowIdx, curr_idx)
			if nex not in state_space:
				state_space[nex] = state_count
				toVisit.addEnd(nex)
				state_count += 1
			nex_idx = state_space[nex]
			colIdx = np.append(colIdx, nex_idx)
			assert curr_idx != nex_idx

		matData = np.append(matData, total)
		rowIdx = np.append(rowIdx, curr_idx)
		colIdx = np.append(colIdx, curr_idx)

		curr = curr.next

	mask = (colIdx != (state_count - 1))
	matData = matData[mask]
	rowIdx = rowIdx[mask]
	colIdx = colIdx[mask]
	assert toVisit.length == state_count
	assert rowIdx.size == colIdx.size == matData.size
	return (matData, rowIdx, colIdx, state_space, state_count)

##############################################################################
##############################################################################
##############################################################################

print ('Average service time (sec):', meanSERV)
print ('Average driving speed (ft/sec):', avgSPEED)

for LOT_COUNT in range(1, LOT_TOTAL + 1):

	print ('started_%s_%s_%s_%s_%s_%s'%(angle, mode, control, m_out, m_in, LOT_COUNT))	
	np.savetxt(dirname + 'started_%s_%s_%s.csv'%(control, filename, LOT_COUNT), np.array([]), delimiter=",")
	print ('N:', LOT_COUNT)
	sys.stdout.flush()

	##############################################################################
	##############################################################################
	##############################################################################
	start = time.time()
	matData, rowIdx, colIdx, state_space, state_count = generate_all(LOT_COUNT, rmu, rlambda)
	print ('Time to generate all states,', time.time() - start)
	sys.stdout.flush()

	start = time.time()
	colIdx = np.append(colIdx, np.ones(state_count) * (state_count - 1) )
	rowIdx = np.append(rowIdx, np.arange(state_count))
	matData = np.append(matData, np.ones(state_count))
	rowIdx = rowIdx.astype(int)
	colIdx = colIdx.astype(int)
	Q = csr_matrix( (matData, (rowIdx, colIdx)) )
	assert state_count == Q.shape[0]
	save_npz(dirname + 'mc_%s_%s_%s.npz'%(control, filename, LOT_COUNT), Q)
	pickle.dump(state_space, open(dirname + 'mc_%s_%s_%s.p'%(control, filename, LOT_COUNT), 'wb'))

	''' to load saved data '''
	# Q = load_npz( dirname + 'mc_%s_%s_%s.npz'%(control, filename, LOT_COUNT) )
	# state_count = Q.shape[0]
	# state_space = pickle.load(open(dirname + 'mc_%s_%s_%s.p'%(control, filename, LOT_COUNT), 'rb'))
	# assert state_count == len(state_space) 

	print ('Total number of states:', state_count)
	print ('Total number of non-zero entries:', Q.getnnz())
	sys.stdout.flush()

	###############################################################################
	################# BiCG with Diagonal Preconditioner ###########################
	###############################################################################
	start = time.time()
	b = np.append(np.zeros(state_count - 1), 1)

	def calc_diag(rowIdx):
		return 1 / Q[rowIdx,:].dot(Q[rowIdx,:].transpose()).A[0][0]

	if True:

		###############################################################################
		######################## Sparse Linear System #################################
		###############################################################################	
		start = time.time()
		x0 = spsolve( Q.transpose(), b )
		# x0 = np.genfromtxt( '%s_%s_%s_%s.csv'%(m_out, m_in, nUnit, LOT_COUNT), delimiter = ',' )
		print ( 'Time to solve sparse matrix:', time.time() - start )
		print ('Sup norm of residual:', np.max( np.abs( Q.transpose().dot(x0) - b ) ) )
		sys.stdout.flush()
		np.savetxt( dirname + 'mc_%s_%s_%s.csv'%(control, filename, LOT_COUNT), x0, delimiter = ',' )

		###############################################################################
		######################## Throughput Evaluation ################################
		###############################################################################
		effec = 0.0
		for state in state_space:
			if state.y[spot2blk(state.N) - 1] == state.N + 1:
				effec += x0[state_space[state]]
			else:
				for j in allspot(spot2blk(state.N)):
					if j <= state.N and state.x[j-1] == m_out + 1:
						effec += x0[state_space[state]]
						break
		print ('Effective rate:', effec)
		thruput = rateDRIV * effec
		print ('Total throughput:', thruput)
		sys.stdout.flush()