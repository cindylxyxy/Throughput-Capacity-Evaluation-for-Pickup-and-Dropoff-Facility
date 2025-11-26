import io
import sys
import csv
import time
import pprint
import random
import json
import pickle
import contextlib
from copy import deepcopy
from math import ceil, floor, exp
from heapq import *
from collections import *

from scipy.sparse import csr_matrix, diags, save_npz, load_npz
from scipy.linalg import solve, svd
from scipy.sparse.linalg import spsolve, cg, inv, eigsh, bicg, bicgstab
import numpy as np
from numpy.linalg import eig, eigh, eigvals, eigvalsh
# import matplotlib.pyplot as plt

from inputSetDef import *
from transitionDef import stateObject
from eventDef import simulationObject
from SLinkedList import SLinkedList
from utils import *
from params import *

find_state = False
find_x = [1, 1, 0, 1, 0]
find_y = [0, 5, 4, 0, 6]

###############################################################################
############################# generate all states #############################
###############################################################################
def generate_all(N):

	start_time = time.time()

	# create the placeholders for output
	state_space = []
	state_dict = dict()
	state_count = 0
	
	# create initial state
	start_state = stateObject(N)
	state_dict[start_state] = state_count
	state_space.append(start_state)
	assert state_space[state_count] == start_state
	assert start_state.idx == state_count
	state_count += 1

	# initialize the while-loop to generate all states
	toVisit = SLinkedList()
	toVisit.addHead(start_state)
	curr = toVisit.head

	while curr:
		curr_idx = state_dict[curr.data]

		# first evaluate and record all mu transitions
		for j in range(1, curr.data.N + 1):
			if curr.data.x[j-1] == 1:
				nex = deepcopy(curr.data)
				nex.mu_transition(j)
				if find_state and nex.x == find_x and nex.y == find_y:
					import pdb; pdb.set_trace()	  
				if nex not in state_dict:
					nex.mu_next = []
					nex.lambda_next = None
					nex.idx = state_count
					state_dict[nex] = state_count
					state_space.append(nex)
					assert state_space[state_count] == nex
					state_count += 1
					toVisit.addEnd(nex)
				else:
					nex = state_space[state_dict[nex]]
				curr.data.mu_next.append(state_dict[nex])

		# second if a lambda transition makes a change to the system state
		if (np.sum(curr.data.y) > 0) | ( (curr.data.x != np.ones(curr.data.N)).any() ):
			nex = deepcopy(curr.data)
			nex.lambda_transition()
			if find_state and nex.x == find_x and nex.y == find_y:
				import pdb; pdb.set_trace()
			if nex not in state_dict:
				nex.mu_next = []
				nex.lambda_next = None
				nex.idx = state_count
				state_dict[nex] = state_count
				state_space.append(nex)
				assert state_space[state_count] == nex
				state_count += 1
				toVisit.addEnd(nex)
			else:
				nex = state_space[state_dict[nex]]
			curr.data.lambda_next = state_dict[nex]
			if state_dict[curr.data] == state_dict[nex]:
				import pdb; pdb.set_trace()

		# go to the next state in queue 
		curr = curr.next

	assert toVisit.length == state_count
	assert state_count == len(state_space)
	assert state_count == len(state_dict)

	print ('Time to generate all states for %s access control:'%control, round(time.time() - start_time, 4))
	print ('Number of states:', state_count)
	sys.stdout.flush()

	return (state_space, state_count, toVisit)


###############################################################################
###################### access controls (value iteration) ######################
###############################################################################
def val_iter_control(state_space, state_count, alpha = DISC_RATE, eps = 1e-7, v0 = None, max_iter = 1e07, gauss_seidel = True):

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

			# first consider the lambda transition (if any)
			nex = curr.lambda_next			
			if nex is not None:
				if nex < curr.idx and gauss_seidel:
					assert not np.isnan( v1[nex] )
					v1[curr.idx] = rateDRIV * v1[nex]
				else:
					v1[curr.idx] = rateDRIV * v0[nex]
				discount += rateDRIV
			elif (np.sum(curr.y) > 0) | ( (curr.x != np.ones(curr.N)).any() ):
				import pdb; pdb.set_trace()
			else:
				v1[curr.idx] = 0.

			# second consider the mu transition(s)
			for (idx, nex) in enumerate(curr.mu_next):
				assert (state_space[nex].y == curr.y)
				j = np.where(np.array(state_space[nex].x) != np.array(curr.x))[0]
				assert len(j) == 1
				j = j[0]
				assert state_space[nex].x[j] == curr.x[j] + 1
				if nex < curr.idx and gauss_seidel:
					assert not np.isnan( v1[nex] )
					v1[curr.idx] += (rateSERV * v1[nex])
				else:
					v1[curr.idx] += (rateSERV * v0[nex])
				discount += rateSERV

			v1[curr.idx] = (curr.reward + v1[curr.idx]) / discount

		residuals.append( np.max( np.abs( v1 - v0 ) ) )
		# residuals.append( np.max(v1 - v0) - np.min(v1 - v0) )
		v0 = v1

		if residuals[-1] < eps:
		# if residuals[-1] < eps * (1 - alpha) / (2 * alpha):
		# if (np.max(v1 - v0) - np.min(v1 - v0)) < eps * ((alpha + rateDRIV + curr.N * rateSERV) / (rateDRIV + curr.N * rateSERV) - 1):
			break

		if iter_count > max_iter:
			print ('Reaching the maximum iterations, but not obtaining the required tolerance ...')
			break

	print ('Time to generate the value function for %s access control:'%control, round( time.time() - start_time, 4) )
	sys.stdout.flush()

	return (v0, residuals)


###############################################################################
################### check states of controls visited in MDP ###################
###############################################################################
def check_states(LOT_COUNT, control_states, control_vals, print_progress = False):

	# mdp_statelist = pickle.load(open('../v2 - in manuscript/mdp_state_list_%s_%s.p'%(filename, LOT_COUNT), 'rb'))
	# mdp_statedict = pickle.load(open('../v2 - in manuscript/mdp_state_dict_%s_%s.p'%(filename, LOT_COUNT), 'rb'))
	# mdp_rewards, mdp_actions, mdp_residuals = pickle.load(open('../v2 - in manuscript/mdp_val_iter_%s_%s.p'%(filename, LOT_COUNT), 'rb'))

	# mdp_statelist = pickle.load(open('2025-11-17p/mdp_state_list_%s_%s.p'%(filename, LOT_COUNT), 'rb'))
	mdp_statedict = pickle.load(open('2025-11-21p/mdp_state_dict_%s_%s.p'%(filename, LOT_COUNT), 'rb'))
	mdp_rewards, mdp_actions, mdp_residuals = pickle.load(open('2025-11-21p/mdp_val_iter_%s_%s.p'%(filename, LOT_COUNT), 'rb'))
	opt_vals = []
	del_idx = []

	for (idx, state) in enumerate(control_states):

		if print_progress and idx % 10000 == 0:	
			print (idx)
			sys.stdout.flush()

		state_ = deepcopy(state)
		y_ = deepcopy(state.y)
		for i in range(1, len(state.y) + 1):
			if 1 <= state.y[i-1] <= state.N and spot2ass(state.y[i-1]) >= i:
				y_[i-1] = state.N + 2
		
		if state.y != y_:
			state_.y = y_
	
		if state_ not in mdp_statedict:
			print ('State in %s access control but not in MDP!' %control)
			print (state.x, state.y, state_.y)
			del_idx.append(idx)
			# import pdb; pdb.set_trace()
		else:			
			mdp_idx = mdp_statedict[state_]
			# assert mdp_statelist[mdp_idx].x == state_.x
			# assert mdp_statelist[mdp_idx].y == state_.y
			opt_vals.append(mdp_rewards[mdp_idx])
			
	control_vals = np.delete(control_vals, del_idx)

	print ('Diff. between Opt. rewards and %s access control rewards (Max, Min, Mean):'%control)
	print (np.max(opt_vals - control_vals), np.min(opt_vals - control_vals), np.mean(opt_vals - control_vals))
	print ('Pct. Diff. between Opt. rewards and %s access control rewards (Max, Min, Mean):'%control)
	print ('{:.6%}'.format(np.max((opt_vals - control_vals) / opt_vals)), 
		'{:.6%}'.format(np.min((opt_vals - control_vals) / opt_vals)), 
		'{:.6%}'.format(np.mean((opt_vals - control_vals) / opt_vals))
		)
	# import pdb; pdb.set_trace()

	# with open('mdp_actions.csv', 'w', newline = '') as csvfile:
	# 	writer = csv.writer(csvfile, delimiter = ',')
	# 	for (idx, state) in enumerate(mdp_statelist):
	# 		if len(state._actions) > 1:
	# 			writer.writerow([1])
	# 		else:
	# 			writer.writerow([0])

	print ('Finished checking %s states in %s access control!' %(len(control_states), control))
	
	return


###############################################################################
######################### access controls (simulation) ########################
###############################################################################
def simulation_control(start_state, N, alpha = DISC_RATE, eps = 1e-7, max_iter = 1e07):

	reward = []
	outCount = []
	inCount = []
	departure_times = []
	sampleDRIV = random.Random(SEED_DRIV[N-1]).sample(range( int(1e12) ), SIM_ITER)
	sampleSERV = random.Random(SEED_SERV[N-1]).sample(range( int(1e12) ), SIM_ITER) 

	for n_iter in range(SIM_ITER):
			
		# initialize the system and state objects
		test = system(N, seedDRIV = sampleDRIV[n_iter], seedSERV = sampleSERV[n_iter])
		test.x = deepcopy(start_state.x)
		test.y = deepcopy(start_state.y)
		test.waiting = deepcopy(start_state.waiting)
		test.inCount = len(start_state.x) - start_state.x.count(0) + len(start_state.y) - start_state.y.count(0)
		# import pdb; pdb.set_trace()
		test.idx = 0
		test.mu_next = []
		test.lambda_next = None
		test.reward = 0.

		# add the mu-transition events
		for j in range(1, test.N + 1):
			if start_state.x[j-1] == 1:
				assert test.curr == 0.
				next_serv = test.timeSERV.next()
				test.add_event( event(next_serv, j, 'mu') )
				test.service_times[j].append(next_serv)

		# add the lambda-transition event (if any)
		if (np.sum(start_state.y) > 0) | ( (start_state.x != np.ones(start_state.N)).any() ):	
			test.add_event( event(test.timeDRIV.next(), None, 'lambd') )
		 
		test.state = test.var
		with contextlib.redirect_stdout(io.StringIO()):
			val = test.run()

		if first_departure:
			assert val == 1
		
		# print (test.reward, test.outCount, test.inCount)
		reward.append(test.reward)
		outCount.append(test.outCount)
		inCount.append(test.inCount)
		# import pdb; pdb.set_trace()

	return (reward, outCount, inCount)


###############################################################################
######################### access controls (simulation) ########################
###############################################################################
def simulation_mdp(LOT_COUNT):

	state_list = pickle.load(open('../mdp_state_list_%s_%s.p'%(filename, LOT_COUNT), 'rb'))
	state_dict = pickle.load(open('../mdp_state_dict_%s_%s.p'%(filename, LOT_COUNT), 'rb'))
	valiter_rewards, actions, residuals = pickle.load(open('../mdp_val_iter_%s_%s.p'%(filename, LOT_COUNT), 'rb'))
	# import pdb; pdb.set_trace()

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
					assert len(state._actions) >= 1
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


##############################################################################
##############################################################################
##############################################################################

control = 'full'
# control = 'partial'

if __name__ == "__main__":

	print ('printing params for %s sided layout ...'%side)
	print ('angle:', angle)
	if angle == 0:
		print ('mode:', mode)
	print ('loading spot configuration:', 'm_{out}_%s_m_{in}_%s_nUnit_%s'%(MOUT, MIIN, nUnit))
	print ('input to distributions:', 'meanSERV_%s_meanDRIV_%s'%(meanSERV, meanDRIV))
	sys.stdout.flush()

	for LOT_COUNT in range(4, 9):

		print ('N:', LOT_COUNT)
		sys.stdout.flush()

		# result = simulation_mdp(LOT_COUNT)
		# pickle.dump(result, open('mdp_sim_%s_%s.p'%(filename, LOT_COUNT), 'wb')

		# generate all states and record the transitions
		# (state_space, state_count, state_list) = generate_all(LOT_COUNT)
		# pickle.dump(state_space, open('%s_states_%s_%s.p'%(control, filename, LOT_COUNT), 'wb'))
		# import pdb; pdb.set_trace()

		# # if needed, read the state space for the specified control and total lot count
		# state_space = pickle.load(open('mc_results/mc_%s_%s_%s.p'%(control, filename, LOT_COUNT), 'rb'))
		state_space = pickle.load(open('2025-11-21p/%s_state_list_%s_%s.p'%(control, filename, LOT_COUNT), 'rb'))
		state_count = len(state_space)
		# import pdb; pdb.set_trace()

		# v0, residuals = val_iter_control(state_space, state_count)
		# pickle.dump(tuple((v0, residuals)), open('%s_val_iter_%s_%s.p'%(control, filename, LOT_COUNT), 'wb'))
		# v0, residuals = pickle.load(open('mdp_related_results/%s_val_iter_%s_%s.p'%(control, filename, LOT_COUNT), 'rb'))
		v0, residuals = pickle.load(open('2025-11-21p/%s_val_iter_%s_%s.p'%(control, filename, LOT_COUNT), 'rb'))
		# import pdb; pdb.set_trace()

		# # # make sure that all states are also in the MDP
		check_states(LOT_COUNT, state_space, v0, print_progress = True)

		# # # run the simulation to validate the value iteration algorithm
		# # rewards = []
		# # rewards_stddev = []
		# # thrupt = []
		# # if first_departure:
		# # 	first_dept_states = []
		# # 	first_dept_expected = []

		# # for (idx, state) in enumerate(state_space):

		# # 	state_departure = False
		# # 	if state.y[-1] == state.N + 1:
		# # 		state_departure = True

		# # 	for j in allspot(spot2blk(state.N)):
		# # 		if state.x[j-1] == m_out + 1:
		# # 			state_departure = True

		# # 	if first_departure and not state_departure:
		# # 		continue

		# # 	# print (idx)
		# # 	reward, outCount, inCount = simulation_control(state, LOT_COUNT) 
		# # 	rewards.append(np.mean(reward))
		# # 	rewards_stddev.append(np.std(reward) / sqrt(len(reward)))
		# # 	thrupt.append(np.mean(outCount))
		# # 	if first_departure:
		# # 		first_dept_states.append(idx)
		# # 		first_dept_expected.append( rateDRIV / (0.005 + rateSERV * state.x.count(1) + rateDRIV) )

		# # import pdb; pdb.set_trace()

		# # pickle.dump(tuple((rewards, rewards_stddev, thrupt)), open(dirname + '%s_sim_%s_%s.p'%(control, filename, LOT_COUNT), 'wb'))