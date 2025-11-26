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
# import matplotlib.pyplot as plt

from transitionDef import stateObject
from SLinkedList import SLinkedList
from params import *
from inputSetDef import mu, spot2ass

###############################################################################
############################# generate all states #############################
###############################################################################
def generate_mdp_states(N, print_progress = False):

	start_time = time.time()

	''' create placeholders for output '''
	state_space = []
	state_dict = dict()
	lambd_func = dict()
	state_count = 0
	action_count = []

	''' create initial state '''
	s0 = stateObject(N)
	s0.get_mdp_actions()
	state_space.append(s0)
	state_dict[s0] = state_count
	assert state_space[state_count] == s0
	assert state_count == s0.idx
	state_count += 1

	''' generate all states '''
	toVisit = SLinkedList()
	toVisit.addHead(s0)
	curr = toVisit.head

	while curr: # list of unexplored states is non-emptry
		
		assert curr.data.idx not in lambd_func
		action_count.append(curr.data.action_count)

		if print_progress and curr.data.idx % 10000 == 0:
			print (curr.data.idx)
			sys.stdout.flush()
			
		''' first evaluate all service completion transitions '''
		for j in range(1, curr.data.N + 1):
			if curr.data.x[j-1] == 1:
			
				nex = deepcopy(curr.data)
				nex.mu_transition(j)
				if nex not in state_dict: # has not been transitioned into
					nex.mu_next = dict()
					nex.lambda_next = []
					nex.idx = state_count
					nex.prevs = []
					nex.get_mdp_actions()
					state_dict[nex] = state_count
					state_space.append(nex)
					assert state_space[state_count] == nex
					state_count += 1
					toVisit.addEnd(nex)
				else:
					nex = state_space[state_dict[nex]]
				curr.data.mu_next[j] = nex.idx
				nex.prevs.append(curr.data.idx)

		''' for each feasible action, evaluate the vehicle movement and maneuver transition '''
		lambd_func[curr.data.idx] = []
		if (np.sum(curr.data.y) > 0) | ( (curr.data.x != np.ones(curr.data.N)).any() ):
			for (idx, action) in enumerate(curr.data.actions):
				nex = deepcopy(curr.data)
				nex.lambda_transition(action)
				if nex not in state_dict:
					nex.mu_next = dict()
					nex.lambda_next = []
					nex.idx = state_count
					nex.prevs = []
					nex.get_mdp_actions()
					state_dict[nex] = state_count
					state_space.append(nex)
					assert state_space[state_count] == nex
					state_count += 1
					toVisit.addEnd(nex)
				else:
					nex = state_space[state_dict[nex]]

				curr.data.lambda_next.append(nex.idx)
				lambd_func[curr.data.idx].append(nex.idx)
				nex.prevs.append(curr.data.idx)


		if len(curr.data.lambda_next) != len(curr.data.actions) and len(curr.data.lambda_next) > 0:
			import pdb; pdb.set_trace()
				
		''' move onto the next state in queue '''
		# pickle.dump(curr.data, open(outdir + 'state_in_process_%s_%s_%s.p'%(filename, curr.data.N, curr.data.idx), 'wb'))
		curr = curr.next

	''' sanity check '''
	assert toVisit.length == state_count
	assert state_count == len(state_space)
	assert state_count == len(state_dict)
	assert state_count == len(lambd_func)

	print ('Time to generate all states and actions:', round( time.time() - start_time, 4) )
	print ('Number of states:', state_count)
	print ('Average number of actions per state:', round( sum(action_count)/state_count, 4) )
	sys.stdout.flush()

	return (state_space, state_dict, lambd_func, action_count)

###############################################################################
################# generate all states under a specific policy #################
###############################################################################
def generate_subspace_optimal(state_space):

	''' create placeholder for the state subspace '''
	state_subspace = [0]
	assert 0 == state_space[0].idx
	state_count = 1

	''' create placeholders for the generator matrix '''
	matData = np.array([]) 
	rowIdx = np.array([])
	colIdx = np.array([])

	''' start enumerating the states '''
	toVisit = SLinkedList()
	toVisit.addHead(0)
	curr = toVisit.head

	while curr: # while the unexplored list is non-empty
		total = 0.0
		state = state_space[curr.data]

		''' generate the set of next states after service completion transitions '''
		for (j, nex_idx) in state.mu_next.items():
			assert (state_space[nex_idx].y == state.y)
			assert len( np.where(np.array(state_space[nex_idx].x) != np.array(state.x))[0] ) == 1
			assert state_space[nex_idx].x[j-1] == state.x[j-1] + 1

			if nex_idx not in state_subspace:
				state_subspace.append(nex_idx)
				toVisit.addEnd(nex_idx)
				state_count += 1

			matData = np.append(matData, mu(j))
			total -= mu(j)
			rowIdx = np.append(rowIdx, state.idx)
			colIdx = np.append(colIdx, nex_idx)

		''' generate the **single** next states after movement and maneuver transition under the optimal control '''
		if (np.sum(state.y) > 0) | ( (state.x != np.ones(state.N)).any() ):
			assert state.opt_lambda_next is not None
			assert state.prob_lambda > 0.

			nex = deepcopy(state)
			nex.lambda_transition(state.opt_access_control)
			assert nex in state_space
			nex_idx = state_space.index(nex)
			assert nex_idx == state.opt_lambda_next
			assert state.idx != nex_idx

			if nex_idx not in state_subspace:
				state_subspace.append(nex_idx)
				toVisit.addEnd(nex_idx)
				state_count += 1

			matData = np.append(matData, rateDRIV)
			total -= rateDRIV
			rowIdx = np.append(rowIdx, state.idx)
			colIdx = np.append(colIdx, nex_idx)
		
		# add a transition from self to self	
		matData = np.append(matData, total)
		rowIdx = np.append(rowIdx, state.idx)
		colIdx = np.append(colIdx, state.idx)

		curr = curr.next

	mask = (colIdx != (state_count - 1))
	matData = matData[mask]
	rowIdx = rowIdx[mask]
	colIdx = colIdx[mask]
	assert toVisit.length == state_count
	assert len(state_subspace) == state_count
	assert rowIdx.size == colIdx.size == matData.size

	return (matData, rowIdx, colIdx, state_subspace)

def generate_policy_states(N, control, print_progress = False):

	start_time = time.time()

	''' create placeholders for output '''
	state_space = []
	state_dict = dict()
	state_count = 0

	''' create initial state '''
	s0 = stateObject(N)
	s0.lambda_next = None
	if control == 'full':
		s0.get_full_access_control()
	elif control == 'partial':
		s0.get_partial_access_control()
	state_space.append(s0)
	state_dict[s0] = state_count
	assert state_space[state_count] == s0
	assert state_count == s0.idx
	state_count += 1

	''' generate all states '''
	toVisit = SLinkedList()
	toVisit.addHead(s0)
	curr = toVisit.head

	while curr: # list of unexplored states is non-emptry

		if print_progress and curr.data.idx % 10000 == 0:
			print (curr.data.idx)
			sys.stdout.flush()

		''' first evaluate all service completion transitions '''
		for j in range(1, curr.data.N + 1):
			if curr.data.x[j-1] == 1:
			
				nex = deepcopy(curr.data)
				nex.mu_transition(j)
				if nex not in state_dict: # has not been transitioned into
					nex.mu_next = dict()
					nex.lambda_next = None
					nex.idx = state_count
					nex.prevs = []
					if control == 'full':
						nex.get_full_access_control()
					elif control == 'partial':
						nex.get_partial_access_control()
					state_dict[nex] = state_count
					state_space.append(nex)
					assert state_space[state_count] == nex
					state_count += 1
					toVisit.addEnd(nex)
				else:
					nex = state_space[state_dict[nex]]
				curr.data.mu_next[j] = nex.idx
				nex.prevs.append(curr.data.idx)

		''' for each feasible action, evaluate the vehicle movement and maneuver transition '''
		if (np.sum(curr.data.y) > 0) | ( (curr.data.x != np.ones(curr.data.N)).any() ):

			if control == 'full':
				action = curr.data.full_access_control
			elif control == 'partial':
				action = curr.data.partial_access_control

			nex = deepcopy(curr.data)
			nex.lambda_transition(action, control)
			if nex not in state_dict:
				nex.mu_next = dict()
				nex.lambda_next = None
				nex.idx = state_count
				nex.prevs = []
				if control == 'full':
					nex.get_full_access_control()
				elif control == 'partial':
					nex.get_partial_access_control()
				state_dict[nex] = state_count
				state_space.append(nex)
				assert state_space[state_count] == nex
				state_count += 1
				toVisit.addEnd(nex)
			else:
				nex = state_space[state_dict[nex]]

			curr.data.lambda_next = nex.idx
			nex.prevs.append(curr.data.idx)
				
		''' move onto the next state in queue '''
		curr = curr.next

	''' sanity check '''
	assert toVisit.length == state_count == len(state_space) == len(state_dict)

	print ('Time to generate all states and actions:', round( time.time() - start_time, 4) )
	print ('Number of states:', state_count)
	sys.stdout.flush()

	return (state_space, state_dict)

###############################################################################
############################### value iteration ###############################
###############################################################################
def val_iter(state_space, alpha = DISC_RATE, eps = 1e-7, v0 = None, max_iter = 1e07, gauss_seidel = True):

	state_count = len(state_space)
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
			
			for (j, nex_idx) in curr.mu_next.items():
				assert (state_space[nex_idx].y == curr.y)
				assert len( np.where(np.array(state_space[nex_idx].x) != np.array(curr.x))[0] ) == 1
				assert state_space[nex_idx].x[j-1] == curr.x[j-1]+1
				assert curr.prob_mu[j-1] == rateSERV
				if nex_idx < curr.idx and gauss_seidel:
					assert not np.isnan( v1[nex_idx] )
					v1[curr.idx] += (curr.prob_mu[j-1] * v1[nex_idx])
				else:
					v1[curr.idx] += (curr.prob_mu[j-1] * v0[nex_idx])
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
			assert curr.prob_lambda == 0
			actions[curr.idx] = None
			curr.opt_access_control = None
			curr.opt_lambda_next = None

		elif len(curr.lambda_next) != len(curr.actions):
			print ('Encountered a state with unequal number of lambda transitions and actions!')
			print (curr.x, curr.y)
			import pdb; pdb.set_trace()		

		else:
			val = 0
			actions[curr.idx] = curr.actions[0]
			nex_idx = curr.lambda_next[0]
			for (idx, nex) in enumerate(curr.lambda_next): 
				if val < curr.prob_lambda * v0[nex]:
					actions[curr.idx] = curr.actions[idx]
					nex_idx = curr.lambda_next[idx]
					val = curr.prob_lambda * v0[nex]
			curr.opt_access_control = deepcopy(actions[curr.idx])
			curr.opt_lambda_next = nex_idx

	print ('Time taken by the value iteration:', round( time.time() - start_time, 4) )
	sys.stdout.flush()

	return (v0, actions, residuals)


def val_iter_policy(state_space, alpha = DISC_RATE, eps = 1e-7, v0 = None, max_iter = 1e07, gauss_seidel = True):

	state_count = len(state_space)
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
			
			if curr.lambda_next is not None: # there is a lambda transition under the specified control policy
				assert curr.prob_lambda == rateDRIV
				discount += rateDRIV
				nex = curr.lambda_next
				if nex < curr.idx and gauss_seidel:
					assert not np.isnan( v1[nex] )
					v1[curr.idx] = max( v1[curr.idx], curr.prob_lambda * v1[nex] )
				else:
					v1[curr.idx] = max( v1[curr.idx], curr.prob_lambda * v0[nex] )
			
			for (j, nex_idx) in curr.mu_next.items():
				assert (state_space[nex_idx].y == curr.y)
				assert len( np.where(np.array(state_space[nex_idx].x) != np.array(curr.x))[0] ) == 1
				assert state_space[nex_idx].x[j-1] == curr.x[j-1]+1
				assert curr.prob_mu[j-1] == rateSERV
				if nex_idx < curr.idx and gauss_seidel:
					assert not np.isnan( v1[nex_idx] )
					v1[curr.idx] += (curr.prob_mu[j-1] * v1[nex_idx])
				else:
					v1[curr.idx] += (curr.prob_mu[j-1] * v0[nex_idx])
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
		v0 = v1

		if residuals[-1] < eps:
			break

		if iter_count > max_iter:
			print ('Reaching the maximum iterations, but not obtaining the required tolerance ...')
			break

	print ('Time taken by the value iteration:', round( time.time() - start_time, 4) )
	sys.stdout.flush()

	return (v0, residuals)

###############################################################################
######################### access controls (simulation) ########################
###############################################################################
def simulate_markov_chain(N, control):

	state_list = pickle.load(open('outdir/%s_state_list_%s_%s.p'%(control, filename, N), 'rb'))
	state_dict = pickle.load(open('outdir/%s_state_dict_%s_%s.p'%(control, filename, N), 'rb'))

	valiter_rewards, actions, residuals = pickle.load(open('outdir/mdp_val_iter_%s_%s.p'%(filename, N), 'rb'))
	import pdb; pdb.set_trace()

	results = dict()

	for start_state in state_list:

		rewards = [[0.0 for _ in range(SIM_ITER)] for _ in TIME_SERIES]
		outCounts = [[0 for _ in range(SIM_ITER)] for _ in TIME_SERIES]
		inCounts = [[0 for _ in range(SIM_ITER)] for _ in TIME_SERIES]

		sampleDRIV = random.Random(SEED_DRIV[LOT_COUNT-1]).sample(range( int(1e12) ), SIM_ITER)
		sampleSERV = random.Random(SEED_SERV[LOT_COUNT-1]).sample(range( int(1e12) ), SIM_ITER) 

		for n_iter in range(SIM_ITER):

			next_rec = 0
			curr_time = 0.0
			
			# initialize the random number generators
			if simType == 'mc':
				timeDRIV = ParamGen(Expo(rateDRIV, seed = sampleDRIV[n_iter]))
				timeSERV = ParamGen(Expo(rateSERV, seed = sampleSERV[n_iter]))
			if simType == 'cav':
				timeDRIV = ParamGen(Cons(meanDRIV, seed = sampleDRIV[n_iter]))
				timeSERV = ParamGen(Expo(rateSERV, seed = sampleSERV[n_iter]))

			# initialize the state 
			state = start_state
			assert (len(state.mu_next) == state.N or len(state.lambda_next) > 0)
			inCount = len(state.x) - state.x.count(0) + len(state.y) - state.y.count(0)
			outCount = 0
			reward = 0
			eventheap = []

			# add the mu-transition events
			for j in range(1, state.N + 1):
				if state.x[j-1] == 1:
					assert curr_time == 0.
					next_serv = timeSERV.next()
					heappush(eventheap, event(next_serv, j, 'mu'))

			# add the lambda-transition event (if any)
			if (np.sum(state.y) > 0) | ( (state.x != np.ones(state.N)).any() ):	
				heappush(eventheap, event(timeDRIV.next(), None, 'lambd'))

			# print ('starting properly...')

			while curr_time <= SIM_UNIT:
					
				if curr_time >= TIME_SERIES[next_rec]:
					rewards[next_rec][n_iter] = reward
					outCounts[next_rec][n_iter] = outCount
					inCounts[next_rec][n_iter] = inCount
					next_rec += 1
					# import pdb;pdb.set_trace()

				curr_event = heappop(eventheap)
				curr_time = float(curr_event.time)
				curr_vehicle = curr_event.vehicle
				curr_typ = curr_event.typ

				################################### update system ###################################				
				prev = state 
				if curr_typ == 'mu':
					# import pdb; pdb.set_trace()
					state = state_list[ state.mu_next[ np.where(np.where(np.array(state.x) == 1)[0] == curr_vehicle - 1)[0][0] ] ]
					# import pdb; pdb.set_trace()
					assert (np.sum(state.y) > 0) | ( (state.x != np.ones(state.N)).any() )
					if (np.sum(prev.y) == 0) and (prev.x == np.ones(state.N)).all():
						heappush(eventheap, event(curr_time + timeDRIV.next(), None, 'lambd'))
				
				else:
					assert curr_typ == 'lambd'
					assert len(state.actions) >= 1
					# import pdb; pdb.set_trace()
					curr_action = actions[state_dict[state]]
					action_idx = state._actions.index(curr_action)
					state = state_list[state.lambda_next[action_idx]]

					# import pdb; pdb.set_trace()
					if curr_action[0] > 0:
						inCount += 1
					if prev.y[-1] == state.N + 1:
						reward += exp(- DISC_RATE * curr_time)
						outCount += 1
					else:
						for j in allspot(spot2blk(state.N)):
							if prev.x[j-1] == m_out + 1:
								reward += exp(- DISC_RATE * curr_time)
								outCount += 1
								break

					# import pdb; pdb.set_trace()
					for j in range(1, state.N + 1):
						if state.x[j-1] == 1 and prev.x[j-1] == m_in + m_out:
							next_serv = timeSERV.next()
							heappush(eventheap, event(curr_time + next_serv, j, 'mu') )

					if (np.sum(state.y) > 0) | ( (state.x != np.ones(state.N)).any() ):
						heappush(eventheap, event(curr_time + timeDRIV.next(), None, 'lambd'))
					
				state.prev = prev

			assert curr_time >= SIM_UNIT
			if not TIME_SERIES[next_rec] == SIM_UNIT:
				import pdb; pdb.set_trace()
			rewards[next_rec][n_iter] = reward
			outCounts[next_rec][n_iter] = outCount
			inCounts[next_rec][n_iter] = inCount

		results[start_state] = {'reward': np.mean(rewards[-1]), 'outCount': np.mean(outCounts[-1]), 'inCount': np.mean(inCounts[-1])}

		plt.figure(1)
		plt.plot(TIME_SERIES, [np.mean(reward) for reward in rewards], label = str(start_state.x) + ', ' + str(start_state.y) )
		# plt.figure(2)
		# plt.plot(TIME_SERIES, [np.mean(outCount) for outCount in outCounts], label = str(start_state.x) + ', ' + str(start_state.y) )
		# plt.figure(3)
		# plt.plot(TIME_SERIES, [np.mean(inCount) for inCount in inCounts], label = str(start_state.x) + ', ' + str(start_state.y) )
	
	plt.figure(1)
	if LOT_COUNT == 1:
		plt.legend(loc = 'lower right', bbox_to_anchor = (1.01, .0))
	plt.xlabel('Simulation Run Time (sec)')
	plt.ylabel('Cumulative Discounted Rewards')
	plt.savefig('mdp_sim_%s_%s.jpg'%(filename, LOT_COUNT))

	# plt.figure(2)
	# plt.legend()
	# plt.xlabel('Simulation Run Time (sec)')
	# plt.ylabel('Cumulative System Outflow')
	# plt.savefig('mdp_sim_outCounts_%s_%s.jpg'%(filename, LOT_COUNT))

	# plt.figure(1)
	# plt.legend()
	# plt.xlabel('Simulation Run Time (sec)')
	# plt.ylabel('Cumulative System Inflow')
	# plt.savefig('mdp_sim_inCounts_%s_%s.jpg'%(filename, LOT_COUNT))
	# import pdb; pdb.set_trace()
	return results

