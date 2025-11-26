########################################################################################################################################################
# Last updated: Feb 22, 2021 by Xinyu Liu
# This file initiates and runs the Markov Chain simulation for the specified configuration and the specified access control. 

########################################################################################################################################################
import sys
import csv
import json
import numpy as np
import pandas as pd
from math import ceil, floor

from heapq import *
from collections import *

from utils import *
from params import *
from eventDef import *


print ('printing params for %s sided layout ...'%side)
print ('angle:', angle)
if angle == 0:
	print ('mode:', mode)
# print ('access control:', control)
print ('type of model/simulation:', simType)
print ('loading spot configuration:', 'm_{out}_%s_m_{in}_%s_nUnit_%s'%(MOUT, MIIN, nUnit)) 
print ('input to distributions:', 'meanSERV_%s_meanDRIV_%s'%(meanSERV, meanDRIV)) 
print ('number of hours in simulation:', SIM_HOUR)
print ('number of simulation iterations:', SIM_ITER)
sys.stdout.flush()


outfile = open(outdir + 'mc_sim_%s.csv'%filename, 'w')
writer = csv.writer(outfile, delimiter = ',')
if side == 'single':
	writer.writerow( ['N'] + list(range(1, 1 + SIM_ITER)) )
else:
	assert side == 'double'
	writer.writerow( ['half_N', 'N'] + list(range(1, 1 + SIM_ITER)) )


for N in range(1, 1 + LOT_TOTAL):

	print (N)
	sys.stdout.flush()
	state_dict = pickle.load(open(outdir + 'mdp_state_dict_%s_%s.p'%(filename, N), 'rb'))
	opt_actions = pickle.load(open(outdir + 'mdp_val_iter_%s_%s.p'%(filename, N), 'rb'))[1]

	test = simulationObject(N, seedDRIV = SEED_DRIV[N-1], seedSERV = SEED_SERV[N-1])

	count = []
	times = []
	disc_rewards = []
	for i in range(SIM_ITER):
		count.append( test.run(state_dict, opt_actions) )
		times.append( test.curr )
		disc_rewards.append( test.simReward)
	
	thruput = [count[0] / times[0]] + [ (count[i+1] - count[i]) / (times[i+1] - times[i]) for i in range(SIM_ITER-1)]
	rewards = [disc_rewards[0] / times[0]] + [ (disc_rewards[i+1] - disc_rewards[i]) / (times[i+1] - times[i]) for i in range(SIM_ITER-1)]
	cycle_times = [times[0]] + [times[i+1] - times[i] for i in range(SIM_ITER-1)]
	
	if side == 'single':
		writer.writerow( [N] + thruput )
	else:
		assert side == 'double'
		writer.writerow( [N, 2 * N] + thruput )

	print ('Throughput/hr:', np.mean(thruput) * 3600, np.std(thruput) / sqrt(len(thruput)) )
	print ('Disc. Thrupt/hr:', np.mean(rewards) * 3600, np.std(rewards) / sqrt(len(rewards)) )

	with open(outdir + 'mdp_service_times_%s_%s.json'%(N, filename), 'w') as write_file:
		json.dump(test.service_times, write_file)
