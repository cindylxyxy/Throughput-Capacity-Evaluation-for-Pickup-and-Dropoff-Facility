########################################################################################################################################################
# Last checked: Feb 06, 2022 by Xinyu Liu
# This file initiates and runs the Markov Chain simulation for the specified configuration and the specified access control. 
# Changes added: [Aug 31, 2021] time-varied simulation
# 				 [Dec 22, 2021] separated simulation by sample path from the mc simulation; added density control
########################################################################################################################################################

import sys
import csv
import json
from math import ceil, floor
import numpy as np
from heapq import *
from collections import *

from eventDef import *
from inputSetDef import spot2blk, spot2in
from utils import *
from params import *


np.savetxt(dirname + 'started_mc_sim_%s.csv'%filename, np.array([]), delimiter=",")
print ('printing params for %s sided layout ...'%side)
print ('angle:', angle)
if angle == 0:
	print ('mode:', mode)
print ('access control:', pcontrol)
print ('type of probability distributions:', simType)
print ('loading spot configuration:', 'm_{out}_%s_m_{in}_%s_nUnit_%s'%(m_out, m_in, nUnit))
print ('number of hours in simulation:', SIM_HOUR)
print ('number of simulation iterations:', SIM_ITER) 
print ('input to distributions:', 'meanSERV_%s_meanDRIV_%s'%(meanSERV, meanDRIV))
print ('sample path match?', spmatch)
sys.stdout.flush()


def simulation():

	if (not debug):

		outfile = open(dirname + 'res_mc_sim_%s_%s.csv'%(control, filename), 'w')
		writer = csv.writer(outfile, delimiter = ',')

		if side == 'single':
			writer.writerow( ['N'] + list(range(1, 1 + SIM_ITER)) )
		else:
			writer.writerow( ['half_N', 'N'] + list(range(1, 1 + SIM_ITER)) )

	for N in range(1, 1 + LOT_TOTAL):

		print ('number of boarding spots on each lane:', N)
		sys.stdout.flush()

		test = system(N, seedDRIV = SEED_DRIV[N-1], seedSERV = SEED_SERV[N-1])

		for j in range(1, test.N + 1):

			assert test.curr == 0.0
			next_serv = test.timeSERV.next()
			test.add_event( event(next_serv, j, 'mu') )
			test.service_times[j].append(next_serv)

		test.state = test.var
		count = []
			
		for i in range(SIM_ITER):
			count.append( test.run() )

		print (test.inCount, test.outCount)
		print (test.x)
		print (test.y)
		count = [count[0]] + [count[i+1] - count[i] for i in range(SIM_ITER-1)]

		if side == 'single':
			writer.writerow( [N] + count )
		else:
			assert side == 'double'
			writer.writerow( [N, 2 * N] + count )

	return 


def simulation_by_density():

	if (not debug):

		outfile = open(dirname + 'outflow_mc_sim_%s_%s.csv'%(control, filename), 'w')
		distfile = open(dirname + 'sysdist_mc_sim_%s_%s.csv'%(control, filename), 'w')
		timefile = open(dirname + 'systime_mc_sim_%s_%s.csv'%(control, filename), 'w')
		writer = csv.writer(outfile, delimiter = ',')
		d_writer = csv.writer(distfile, delimiter = ',')
		t_writer = csv.writer(timefile, delimiter = ',')

		if side == 'single':
			writer.writerow( ['N', 'maxEntry', # 'max.in.sys', 'avg.in.sys', 
							  'sysFlow', 'sysDens'] + list(range(1, 1 + SIM_ITER)) )
			d_writer.writerow(['N', 'maxEntry'] + list(range(1, 1 + SIM_ITER)))
			t_writer.writerow(['N', 'maxEntry'] + list(range(1, 1 + SIM_ITER)))
		else:
			writer.writerow( ['half_N', 'N', 'maxEntry', # 'max.in.sys', 'avg.in.sys', 
							  'sysFlow', 'sysDens'] + list(range(1, 1 + SIM_ITER)) )
			d_writer.writerow(['half_N', 'N', 'maxEntry'] + list(range(1, 1 + SIM_ITER)))
			t_writer.writerow(['half_N', 'N', 'maxEntry'] + list(range(1, 1 + SIM_ITER)))

	resDict = dict()

	for N in range(38, 39):

		print ('number of boarding spots on each lane:', N)
		sys.stdout.flush()
		resDict[N] = dict()

		if side == 'single':
			Nmax = N
		else:
			assert side == 'double'
			Nmax = 2 * N

		if angle == 0:
			Nmax += (N + 1)
		else:
			assert angle == 90
			Nmax += (spot2blk(N) + 1)

		# for probEntry in np.arange(0.1, 1.1, 0.1):
		for maxEntry in np.arange(1, Nmax + 2, 2):
		# for maxEntry in [43, 44, 45, 46, 47, 63, 64, 65]:

			print ('max. number of vehicles allowed in the facility:', maxEntry)
			sys.stdout.flush()
			resDict[N][maxEntry] = dict()

			test = system(N, seedDRIV = SEED_DRIV[N-1], seedSERV = SEED_SERV[N-1])

			# test.probEntry = prob
			# test.randomEntry = ParamGen(Unif(2., seed = SEED_ENTR[N-1]))

			test.maxEntry = maxEntry

			for j in range(1, test.N + 1):
				assert test.curr == 0.0
				next_serv = test.timeSERV.next()
				test.add_event( event(next_serv, j, 'mu') )
				test.service_times[j].append(next_serv)

			test.state = test.var

			outCount = []
			inCount = []
			sysDist = []
			sysTime = []

			curr_out = test.run(warm_up = True)
			
			for i in range(SIM_ITER + 1):

				print (i)
				sys.stdout.flush()

				curr_in = curr_out
				curr_dist = curr_out * (test.n * CAR_LENGTH)
				curr_time = test.totalTime
				outCount.append( curr_out )
				
				assert test.probEntry is not None or test.maxEntry is not None
				for j in range(len(test.x) + 1):
					if test.x[j-1] > 0:
						curr_in += 1
						curr_dist += (j * LOT_LENGTH)
				for j in range(len(test.y) + 1):
					if test.y[j-1] > 0:
						curr_in += 1
						curr_dist += (j * CAR_LENGTH)

				sysDist.append( curr_dist )
				sysTime.append( curr_time )
				inCount.append( curr_in )
				if i < SIM_ITER:
					curr_out = test.run()

			total_out = outCount[-1]
			total_in = inCount[-1]
			iter_outCounts = [outCount[i+1] - outCount[i] for i in range(SIM_ITER)]
			iter_inCounts  = [inCount[i+1] - inCount[i] for i in range(SIM_ITER)]
			iter_sysDists  = [sysDist[i+1] - sysDist[i] for i in range(SIM_ITER)]
			iter_sysTimes  = [sysTime[i+1] - sysTime[i] for i in range(SIM_ITER)]

			resDict[N][maxEntry]['inCount'] = {'val': np.mean(iter_inCounts) / SIM_HOUR, 'stdev': np.std(iter_inCounts) / sqrt(SIM_ITER) / SIM_HOUR}
			resDict[N][maxEntry]['outCount'] = {'val': np.mean(iter_outCounts) / SIM_HOUR, 'stdev': np.std(iter_outCounts) / sqrt(SIM_ITER) / SIM_HOUR}
			resDict[N][maxEntry]['sysFlow'] = {'val': np.mean(iter_sysDists) / SIM_HOUR / (test.n * CAR_LENGTH), 'stdev': np.std(iter_sysDists) / sqrt(SIM_ITER) / SIM_HOUR / (test.n * CAR_LENGTH)}
			resDict[N][maxEntry]['sysDens'] = {'val': np.mean(iter_sysTimes) / 3600. / SIM_HOUR / (test.n * CAR_LENGTH), 'stdev': np.std(iter_sysTimes) / sqrt(SIM_ITER) / 3600. / SIM_HOUR / (test.n * CAR_LENGTH)}

			if (not debug):

				if side == 'single':
					writer.writerow([N, maxEntry,
									 resDict[N][maxEntry]['sysFlow']['val'], resDict[N][maxEntry]['sysDens']['val']] + iter_outCounts )
					d_writer.writerow([N, maxEntry] + iter_sysDists)
					t_writer.writerow([N, maxEntry] + iter_sysTimes)
				else:
					writer.writerow([N, 2 * N, maxEntry,
									 resDict[N][maxEntry]['sysFlow']['val'], resDict[N][maxEntry]['sysDens']['val']] + iter_outCounts )
					d_writer.writerow([N, 2 * N, maxEntry] + iter_sysDists)
					t_writer.writerow([N, 2 * N, maxEntry] + iter_sysTimes)

			print (total_out, total_in)
			sys.stdout.flush()

	with open(dirname + 'plotdata_%s_%s.pkl'%(control, filename), 'wb') as f:
		pickle.dump(resDict, f)

	return 


def simulation_by_time(LOT_COUNT):

	state_list = pickle.load(open('%s_states_%s_%s.p'%(control, filename, LOT_COUNT), 'rb'))
	valiter_rewards, residuals = pickle.load(open('%s_val_iter_%s_%s.p'%(control, filename, LOT_COUNT), 'rb'))

	results = dict()

	for start_state in state_list:
		
		rewards = [[0.0 for _ in range(SIM_ITER)] for _ in TIME_SERIES]
		outCounts = [[0 for _ in range(SIM_ITER)] for _ in TIME_SERIES]
		inCounts = [[0 for _ in range(SIM_ITER)] for _ in TIME_SERIES]

		sampleDRIV = random.Random(SEED_DRIV[LOT_COUNT-1]).sample(range( int(1e12) ), SIM_ITER)
		sampleSERV = random.Random(SEED_SERV[LOT_COUNT-1]).sample(range( int(1e12) ), SIM_ITER) 
		
		for n_iter in range(SIM_ITER):

			next_rec = 0

			# initialize the system and state objects
			test = system(LOT_COUNT, seedDRIV = sampleDRIV[n_iter], seedSERV = sampleSERV[n_iter])
			test.x = deepcopy(start_state.x)
			test.y = deepcopy(start_state.y)
			test.waiting = deepcopy(start_state.waiting)
			test.inCount = len(start_state.x) - start_state.x.count(0) + len(start_state.y) - start_state.y.count(0)
			test.outCount = 0
			test.reward = 0.

			# add the mu-transition events
			for j in range(1, test.N + 1):
				if test.x[j-1] == 1:
					assert test.curr == 0.
					next_serv = test.timeSERV.next()
					test.add_event( event(next_serv, j, 'mu'))
					test.service_times[j].append(next_serv)

			# add the lambda-transition event (if any)
			if (np.sum(test.y) > 0) | ( (test.x != np.ones(test.N)).any() ):
				test.add_event( event(test.timeDRIV.next(), None, 'lambd') )

			while test.curr <= SIM_UNIT:

				if test.curr >= TIME_SERIES[next_rec]:
					rewards[next_rec][n_iter] = test.reward
					outCounts[next_rec][n_iter] = test.outCount
					inCounts[next_rec][n_iter] = test.inCount
					next_rec += 1

				curr_event = heappop(test.eventheap)
				test.curr = float(curr_event.time)
				curr_vehicle = curr_event.vehicle
				curr_typ = curr_event.typ

				################################### update system ###################################				
				x = deepcopy(test.x)
				y = deepcopy(test.y)
				if curr_typ == 'mu':
					test.mu_transition(curr_vehicle)
					assert (np.sum(test.y) > 0) | ( (test.x != np.ones(test.N)).any() )
					if (np.sum(y) == 0) and (x == np.ones(test.N)).all():
						test.add_event( event(test.curr + test.timeDRIV.next(), None, 'lambd'))
				
				else:
					assert curr_typ == 'lambd'
					test.lambda_transition()
					if (np.sum(test.y) > 0) | ( (test.x != np.ones(test.N)).any() ):
						test.add_event( event(test.curr + test.timeDRIV.next(), None, 'lambd'))

			assert test.curr >= SIM_UNIT
			if not TIME_SERIES[next_rec] == SIM_UNIT:
				import pdb; pdb.set_trace()
			rewards[next_rec][n_iter] = test.reward
			outCounts[next_rec][n_iter] = test.outCount
			inCounts[next_rec][n_iter] = test.inCount

		results[start_state] = {'reward': np.mean(rewards[-1]), 'outCount': np.mean(outCounts[-1]), 'inCount': np.mean(inCounts[-1])}

	return results


def simulation_by_sample_path():

	for N in range(38, 39):

		print ('number of boarding spots on each lane:', N)
		sys.stdout.flush()

		test = system(N, seedDRIV = SEED_DRIV[N-1], seedSERV = SEED_SERV[N-1])

		###### sample path example starts here #####
		''' Case 1: SS0S Partial Congestion in Loop '''
		test.x = [2, 1, 1, 0, 1, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		test.y = [0, 9, 4, 8, 6, 7, 0, 0, 0, 0, 39, 39, 39, 39, 39, 39, 39, 39, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 0, 29, 0, 0, 0, 0, 0, 0, 0, 0, 0]

		''' Case 2: Small perturbation to the initial state (decreasing sequence of spot indices from exit to entrance) '''
		# test.x = [2, 1, 1, 0, 1, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		# test.y = [0, 4, 6, 7, 8, 9, 0, 0, 0, 0, 39, 39, 39, 39, 39, 39, 39, 39, 39, 29, 30, 31, 32, 33, 34, 35, 36, 37, 0, 38, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		
		test.maxEntry = 45

		test.waiting = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 1, 25]
		test.inSystem = 47
		test.onlane = 24
		test.totalTime = None

		test.outfile = open(dirname + 'sample_path_mc_sim_%s_%s_%s.csv'%(control, filename, N), 'w')
		test.writer = csv.writer(test.outfile, delimiter = ',')
		
		if side == 'single':
			if len(test.x) == len(test.y): 	
				test.writer.writerow( ['', 'x'] + test.x )
			else:
				assert len(test.x) + 1 == len(test.y)
				test.writer.writerow( ['', 'x'] + test.x + [0])
			test.writer.writerow( ['', 'y'] + test.y )
		else:
			assert side == 'double'
			if np.abs( 0.5 * len(test.x) - len(test.y) ) <= 1e-10: 	
				test.writer.writerow( ['', 'x'] + list(np.array(test.x)[list(range(0, len(test.x), 2))]) )
			else:
				assert np.abs( 0.5 * len(test.x) + 1 - len(test.y) ) <= 1e-10 
				test.writer.writerow( ['', 'x'] + list(np.array(test.x)[list(range(0, len(test.x), 2))] + [0]) )
			test.writer.writerow( ['', 'y'] + test.y )
			if np.abs( 0.5 * len(test.x) - len(test.y) ) <= 1e-10: 	
				test.writer.writerow( ['', 'x'] + list(np.array(test.x)[list(range(1, len(test.x), 2))]) )
			else:
				assert np.abs( 0.5 * len(test.x) + 1 - len(test.y) ) <= 1e-10 
				test.writer.writerow( ['', 'x'] + list(np.array(test.x)[list(range(1, len(test.x), 2))] + [0]) )

		for j in range(1, test.N + 1):
			assert test.curr == 0.0
			if test.x[j - 1] == 1:
				next_serv = test.timeSERV.next()
				test.add_event( event(next_serv, j, 'mu') )
				test.service_times[j].append(next_serv)
		test.add_event( event(test.curr + test.timeDRIV.next(), None, 'lambd') )
		test.run()
	
	return


'''
for LOT_COUNT in range(1, LOT_TOTAL + 1):
	print ('N:', LOT_COUNT)
	sys.stdout.flush()

	result = simulation_by_time(LOT_COUNT)
	pickle.dump(result, open('%s_sim_%s_%s.p'%(control, filename, LOT_COUNT), 'wb'))
'''

# simulation_by_density()

simulation()

# simulation_by_sample_path()