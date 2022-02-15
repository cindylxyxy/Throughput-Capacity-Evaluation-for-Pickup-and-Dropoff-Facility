########################################################################################################################################################
# Last updated: Feb 06, 2022 by Xinyu Liu
# This file runs the Markov Decision Process model and generates all state-action pairs.
 
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

from transitionDef import stateObject
from SLinkedList import SLinkedList
from params import *


find_state = False
find_x = [5,1,0,5,0]
find_y = [0,7,4,0,0]


###############################################################################
############################# generate all states #############################
###############################################################################
def generate_all(N):

	start_time = time.time()

	# create the placeholders for output
	state_space = []
	state_dict = dict()
	lambda_func = dict()
	state_count = 0
	action_count = []

	# create initial state
	start_state = stateObject(N)
	start_state.get_actions()
	state_dict[start_state] = state_count
	state_space.append(start_state)
	assert state_space[state_count] == start_state
	assert state_count == start_state.idx
	state_count += 1

	# initialize the while-loop to generate all states
	toVisit = SLinkedList()
	toVisit.addHead(start_state)
	curr = toVisit.head

	while curr:
		assert curr.data.idx not in lambda_func
		curr_idx = curr.data.idx
		action_count.append(curr.data.action_count)

		if curr_idx % 10000 == 0:
			print (curr_idx)
			sys.stdout.flush()

		# first evaluate and record all mu transitions
		for j in range(1, curr.data.N + 1):
			if curr.data.x[j-1] == 1:
				nex = deepcopy(curr.data)
				nex.mu_transition(j)
				if find_state and nex.x == find_x and nex.y == find_y:
					import pdb; pdb.set_trace()
				if nex not in state_dict:
					nex.mu_next = []
					nex.lambda_next = []
					nex.idx = state_count
					nex.prevs = []
					nex.get_actions()
					state_dict[nex] = state_count
					state_space.append(nex)
					assert state_space[state_count] == nex
					state_count += 1
					toVisit.addEnd(nex)
				else:
					nex = state_space[state_dict[nex]]
				curr.data.mu_next.append(nex.idx)
				nex.prevs.append(curr.data.idx)

		# now evaluate the lambda transitions for all feasible actions
		lambda_func[curr.data.idx] = []
		if (np.sum(curr.data.y) > 0) | ( (curr.data.x != np.ones(curr.data.N)).any() ):
			for (idx, action) in enumerate(curr.data.actions):
				nex = deepcopy(curr.data)
				nex.lambda_transition(action)
				if find_state and nex.x == find_x and nex.y == find_y:
					import pdb; pdb.set_trace()
				if nex not in state_dict:
					nex.mu_next = []
					nex.lambda_next = []
					nex.idx = state_count
					nex.prevs = []
					nex.get_actions()
					state_dict[nex] = state_count
					state_space.append(nex)
					assert state_space[state_count] == nex
					state_count += 1
					toVisit.addEnd(nex)
				else:
					nex = state_space[state_dict[nex]]
				curr.data.lambda_next.append(nex.idx)
				lambda_func[curr.data.idx].append(nex.idx)
				nex.prevs.append(curr.data.idx)

		if (0 < len(curr.data.lambda_next) != len(curr.data.actions)):
			print ('Encountered a state with unequal number of lambda transitions and actions!')
			print (curr.data.x, curr.data.y, curr.data.actions)
			for idx in curr.data.lambda_next:
				print (state_space[idx])
			import pdb; pdb.set_trace()
				
		# go to the next state in queue
		curr = curr.next

	assert toVisit.length == state_count
	assert state_count == len(state_space)
	assert state_count == len(state_dict)
	assert state_count == len(lambda_func)

	print ('Time to generate all states and actions:', round( time.time() - start_time, 4) )
	print ('Number of states:', state_count)
	print ('Average number of actions per state:', round( sum(action_count)/state_count, 4) )
	sys.stdout.flush()

	return (state_space, state_dict, lambda_func, state_count, toVisit, action_count)


###############################################################################
############################### value iteration ###############################
###############################################################################
def val_iter(state_space, state_count, alpha = DISC_RATE, eps = 1e-7, v0 = None, max_iter = 1e07, gauss_seidel = True):

	start_time = time.time()

	# keep track of the number of iterations
	iter_count = 0

	# initialize the value function
	if v0 is None:
		v0 = np.zeros( state_count )

	# initialize thee placeholder for sup norm of v_{n+1} - v_n
	residuals = []

	while True:
		
		iter_count += 1
		if iter_count % 1000 == 0:
			print (iter_count)

		v1 = np.zeros(state_count) * np.nan

		for curr in state_space:
			
			discount = alpha

			v1[curr.idx] = 0
			
			if len(curr.lambda_next) > 0:
				discount += rateDRIV

			for (idx, nex) in enumerate(curr.lambda_next):
				assert curr.prob_lambda == rateDRIV
				if nex < curr.idx and gauss_seidel:
					assert not np.isnan( v1[nex] )
					v1[curr.idx] = max( v1[curr.idx], curr.prob_lambda * v1[nex] )
				else:
					v1[curr.idx] = max( v1[curr.idx], curr.prob_lambda * v0[nex] )
			
			for (idx, nex) in enumerate(curr.mu_next):
				assert (state_space[nex].y == curr.y)
				j = np.where(np.array(state_space[nex].x) != np.array(curr.x))[0]
				assert len(j) == 1
				j = j[0]
				assert state_space[nex].x[j] == curr.x[j] + 1
				assert curr.prob_mu[j] == rateSERV
				if nex < curr.idx and gauss_seidel:
					assert not np.isnan( v1[nex] )
					v1[curr.idx] += (curr.prob_mu[j] * v1[nex])
				else:
					v1[curr.idx] += (curr.prob_mu[j] * v0[nex])
				discount += rateSERV

			if np.abs(alpha + curr.prob_lambda + sum(curr.prob_mu) - discount) > 1e-10:
				import pdb; pdb.set_trace()

			if curr.reward > 0 and curr.reward != rateDRIV:
				if np.abs( curr.reward * (rateDRIV + rateSERV * curr.x.count(1)) - rateDRIV ) > 1e-10:
					if not curr.x.count(2) > 0:
						import pdb; pdb.set_trace()
				curr.reward = rateDRIV

			assert curr.reward == 0.0 or curr.prob_lambda == rateDRIV

			v1[curr.idx] = (curr.reward + v1[curr.idx]) / (alpha + curr.prob_lambda + sum(curr.prob_mu))

		residuals.append( np.max( np.abs( v1 - v0 ) ) )
		# residuals.append( np.max(v1 - v0) - np.min(v1 - v0) )
		v0 = v1

		if residuals[-1] < eps:
			break

		if iter_count > max_iter:
			print ('Reaching the maximum iterations, but not obtaining the required tolerance ...')
			break
			
	actions = [ [] for _ in range(state_count) ]
	for curr in state_space:
		if len(curr.lambda_next) == 0:
			if len(curr.actions) == 0:
				print ('Encountered a state with no lambda transition but zero action!')
				print (curr.x, curr.y)
				import pdb; pdb.set_trace()
			assert len(curr.actions) == 1
			actions[curr.idx] = curr.actions[0]
		elif len(curr.lambda_next) != len(curr.actions):
			print ('Encountered a state with unequal number of lambda transitions and actions!')
			print (curr.x, curr.y)
			import pdb; pdb.set_trace()			
		else:
			val = 0
			actions[curr.idx] = curr.actions[0]
			for (idx, nex) in enumerate(curr.lambda_next): 
				if val < curr.prob_lambda * v0[nex]:
					actions[curr.idx] = curr.actions[idx]
					val = curr.prob_lambda * v0[nex]

	print ('Time taken by the value iteration:', round( time.time() - start_time, 4) )
	sys.stdout.flush()

	return (v0, actions, residuals)

##############################################################################
##############################################################################
##############################################################################

if __name__ == "__main__":

	print ('Service rate:', rateSERV)
	print ('Average driving speed:', rateDRIV)

	for LOT_COUNT in range(1, LOT_TOTAL + 1):

		# print ('started_mdp_%s_%s'%(filename, LOT_COUNT))	
		# np.savetxt(dirname + 'started_mdp_%s_%s.csv'%(filename, LOT_COUNT), np.array([]), delimiter=",")
		print ('N:', LOT_COUNT)
		sys.stdout.flush()

		# to generate the (state, action) space 
		(state_space, state_dict, lambda_func, state_count, state_list, action_count) = generate_all(LOT_COUNT)
		pickle.dump(state_space, open(dirname + 'mdp_state_list_%s_%s.p'%(filename, LOT_COUNT), 'wb'))
		pickle.dump(state_dict, open(dirname + 'mdp_state_dict_%s_%s.p'%(filename, LOT_COUNT), 'wb'))
		pickle.dump(lambda_func, open(dirname + 'mdp_lambda_%s_%s.p'%(filename, LOT_COUNT), 'wb'))
		pickle.dump(action_count, open(dirname + 'mdp_actions_%s_%s.p'%(filename, LOT_COUNT), 'wb'))
		# import pdb; pdb.set_trace()

		# # to load the generated (state, action) space 
		# state_space = pickle.load(open(dirname + 'mdp_state_list_%s_%s.p'%(filename, LOT_COUNT), 'rb'))
		# state_dict = pickle.load(open(dirname + 'mdp_state_dict_%s_%s.p'%(filename, LOT_COUNT), 'rb'))
		# lambda_func = pickle.load(open(dirname + 'mdp_lambda_%s_%s.p'%(filename, LOT_COUNT), 'rb'))
		# action_count = pickle.load(open(dirname + 'mdp_actions_%s_%s.p'%(filename, LOT_COUNT), 'rb'))
		# state_count = len(state_space)
		# print (state_count)
		
		# to run the value iteration algorithm and save the output
		v0, actions, residuals = val_iter(state_space, state_count)
		pickle.dump(tuple((v0, actions, residuals)), open(dirname + 'mdp_val_iter_%s_%s.p'%(filename, LOT_COUNT), 'wb'))		