########################################################################################################################################################
# Last checked: Feb 06, 2022 by Xinyu Liu
# This file initiates and runs the trajectory-based simulation for the specified configuration and the specified access control. 

########################################################################################################################################################
import sys
import csv
import pickle
from copy import deepcopy
import numpy as np
from params import *
from event import *
from eventShort import *
from event90 import *
from eventSingleNew import *
from utils import *

np.savetxt(dirname + 'started_traj_%s.csv'%filename, np.array([]), delimiter=",")
print ('printing params for %s sided layout ...'%side)
print ('angle:', angle)
if angle == 0:
	print ('mode:', mode)
print ('access control:', pcontrol)
print ('type of simulation:', simType)
print ('loading spot configuration:', 'Lot length is %s, lot width is %s, and lane width is %s'%(LOT_LENGTH, LOT_WIDTH, LANE_WIDTH) )
print ('number of hours in simulation:', SIM_HOUR)
print ('number of simulation iterations:', SIM_ITER) 
print ('input to distributions:', 'meanSERV is %s, speed is %s, meanPOUT is %s, and meanPLIN is %s'%(meanSERV, rateDRIV, meanPOUT, meanPLIN) ) 
print ('sample path match?', spmatch)
sys.stdout.flush()

if (not debug):
	outfile = open(dirname + 'res_traj_%s.csv'%filename, 'w')
	writer = csv.writer(outfile, delimiter = ',')
	if side == 'single':
		writer.writerow( ['N', 'max_idle_time', 'max_idle_veh', 'max_idle_pct', 'max_idle_pct_veh', 'avg_prod', 'avg_idle', 'avg_pct_idle'] + list(range(1, 1 + SIM_ITER)) )
	else:
		writer.writerow( ['half_N', 'N', 'max_idle_time', 'max_idle_veh', 'max_idle_pct', 'max_idle_pct_veh', 'avg_prod', 'avg_idle', 'avg_pct_idle'] + list(range(1, 1 + SIM_ITER)) )

for N in range(1, 1 + LOT_TOTAL):

	print (N)
	sys.stdout.flush()

	if angle == 0 and mode == 'long':
		test = system(N, SEED_SERV[N-1], SEED_POUT[N-1], SEED_PLIN[N-1], SEED_DRIV[N-1])
	elif angle == 0 and mode == 'short':
		test = systemShort(N, SEED_SERV[N-1], SEED_POUT[N-1], SEED_PLIN[N-1], SEED_DRIV[N-1])
	elif angle == 0 and mode == 'singlelane':
		test = SingleLaneSystem(N)
	elif angle == 90:
		test = system90(N, SEED_SERV[N-1], SEED_POUT[N-1], SEED_PLIN[N-1], SEED_DRIV[N-1])
	else:
		print ('Configuration not specified!')
		import pdb; pdb.set_trace()

	if mode == 'singlelane':
		test.inCount += 1
		car = SLvehicle(test, getin = True, stop_idx = test.N)
		test.add_event( event(car.serv_time, car, 'prepare_pulling_out') )
		for j in range(test.N - 1, 0, -1):
			test.inCount += 1
			car.nex = SLvehicle(test, getin = True, stop_idx = j, prev = car)
			car = car.nex
			test.add_event( event(car.serv_time, car, 'prepare_pulling_out') )
	else:
		for j in range(1, test.N + 1):
			test.inCount += 1
			car = vehicle(test, getin = True, stop_idx = j)
			test.inservice[j-1] = car
			test.add_event( event(car.serv_time, car, 'prepare_pulling_out') )

	count = []
	for i in range(SIM_ITER):
		count.append( test.run() )

	if delay:
		with open(dirname + 'delay_traj_%s_%s.csv'%(filename, N), 'w') as delayfile:
			delay_writer = csv.writer(delayfile, delimiter = ',')
			delay_writer.writerow( ['N'] + list(sorted(test.wait_out[0].keys())) )
			for i in range(1, test.N + 1):
				wait_out = deepcopy(test.wait_out[i-1])
				delay_writer.writerow( [i] + [wait_out[key] for key in list(sorted(test.wait_out[0].keys()))] )

	total_out = count[-1]
	total_in = total_out
	count = [count[0]] + [count[i+1] - count[i] for i in range(SIM_ITER-1)]

	car = test.head
	while car is not None:
		total_in += 1
		car.update_loc()
		if car.status == 1:
			if not car.prod_time == 0.0:
				import pdb; pdb.set_trace()
			if car.curr_loc == 0.0:
				total_in -= 1
				car = car.nex
				continue
			else:
				car.prod_time += (car.curr_loc / car.driv)

		elif car.status == 2:
			if not car.plin_end >= test.curr:
				import pdb; pdb.set_trace()
			car.prod_time -= (car.plin_end - test.curr)

		elif car.status == 3:
			assert mode == 'singlelane'
			total_in += 1
			if not car.serv_end >= test.curr:
				import pdb; pdb.set_trace()
			car.prod_time -= (car.serv_end - test.curr)

		elif car.status == 5:
			if not car.pout_end >= test.curr:
				import pdb; pdb.set_trace()
			car.prod_time -= (car.pout_end - test.curr)

		else:
			assert car.status == 6
			if spmatch and angle == 90:
				start_x = car.block_idx * CAR_LENGTH
			elif angle == 90:
				start_x = (car.stop + dgap) * LOT_LENGTH
			else:
				start_x = car.stop * LOT_LENGTH
			if not car.curr_loc >= start_x:
				import pdb; pdb.set_trace()
			car.prod_time += ((car.curr_loc - start_x) / car.driv)

		test.idle_time_calc(car)
		car = car.nex

	if not mode == 'singlelane':
		for stop_idx in range(1, test.N + 1):
			if test.inservice[stop_idx - 1] is not None and test.inservice[stop_idx - 1].status in [3, 4]:
				car = test.inservice[stop_idx - 1]
				total_in += 1
				if car.status == 3:
					if not car.serv_end >= test.curr:
						import pdb; pdb.set_trace()
					car.prod_time -= (car.serv_end - test.curr)
				test.idle_time_calc(car)

	total = total_in

	if (not debug):

		if side == 'single':
			writer.writerow([N, 
							 test.max_idle, test.max_idle_idx, test.max_pct_idle, test.max_pct_idle_idx,
							 np.sum(test.prod_time) / total, 
							 np.sum(test.idle_time) / total,
							 np.sum(test.pct_idle) / total ] + count )
		else:
			writer.writerow([N, 2 * N, 
							 test.max_idle, test.max_idle_idx, test.max_pct_idle, test.max_pct_idle_idx,
							 np.sum(test.prod_time) / total, 
							 np.sum(test.idle_time) / total,
							 np.sum(test.pct_idle) / total ] + count )

	print (count, np.sum(test.prod_time) / total, np.sum(test.idle_time) / total, np.sum(test.pct_idle) / total)
