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
from markov import *

''' Set hyperparameters '''
MIN_LOT_COUNT = 4
MAX_LOT_COUNT = 10
policy_type = 'mdp'
# policy_type = 'full'
# policy_type = 'partial'

''' Print model setup '''
print ('printing params for %s sided layout ...'%side)
print ('angle:', angle)
if angle == 0:
	print ('mode:', mode)
print ('type of policy being evaluated:', policy_type)
print ('loading spot configuration:', 'm_{out}_%s_m_{in}_%s_nUnit_%s'%(MOUT, MIIN, nUnit))
print ('input to distributions:', 'meanSERV_%s_meanDRIV_%s'%(meanSERV, meanDRIV))
sys.stdout.flush()

for LOT_COUNT in range(MIN_LOT_COUNT, MAX_LOT_COUNT+1):
	print ('N:', LOT_COUNT)
	sys.stdout.flush()
	
	if policy_type == 'mdp': # generate MDP and MC under the optimal policy
		(state_list, state_dict, lambda_func, action_count) = generate_mdp_state(LOT_COUNT, print_progress = True)
		v0, opt_actions, residuals = val_iter(state_list)
		# (matData, rowIdx, colIdx, state_subspace) = generate_subspace_optimal(state_list)
		# print ('sizes of state spaces of MDP and optimal MC:', len(state_list), len(state_subspace))

		pickle.dump(state_list, open(outdir + 'mdp_state_list_%s_%s.p'%(filename, LOT_COUNT), 'wb'))
		pickle.dump(state_dict, open(outdir + 'mdp_state_dict_%s_%s.p'%(filename, LOT_COUNT), 'wb'))
		pickle.dump(lambda_func, open(outdir + 'mdp_lambda_%s_%s.p'%(filename, LOT_COUNT), 'wb'))
		pickle.dump(action_count, open(outdir + 'mdp_actions_%s_%s.p'%(filename, LOT_COUNT), 'wb'))
		# import pdb; pdb.set_trace()
		pickle.dump(tuple((v0, opt_actions, residuals)), open(outdir + 'mdp_val_iter_%s_%s.p'%(filename, LOT_COUNT), 'wb'))	
		# import pdb; pdb.set_trace()

	else: # generate MC under an approximate policy
		(state_list, state_dict) = generate_policy_states(LOT_COUNT, control = policy_type)
		v0, residuals = val_iter_policy(state_list)

		pickle.dump(state_list, open(outdir + '%s_state_list_%s_%s.p'%(policy_type, filename, LOT_COUNT), 'wb'))
		pickle.dump(state_dict, open(outdir + '%s_state_dict_%s_%s.p'%(policy_type, filename, LOT_COUNT), 'wb'))
		# import pdb; pdb.set_trace()
		pickle.dump(tuple((v0, residuals)), open(outdir + '%s_val_iter_%s_%s.p'%(policy_type, filename, LOT_COUNT), 'wb'))	
		# import pdb; pdb.set_trace()

	# state_space = pickle.load(open(outdir + 'mdp_state_list_%s_%s.p'%(filename, LOT_COUNT), 'rb'))
	# state_dict = pickle.load(open(outdir + 'mdp_state_dict_%s_%s.p'%(filename, LOT_COUNT), 'rb'))
	# lambda_func = pickle.load(open(outdir + 'mdp_lambda_%s_%s.p'%(filename, LOT_COUNT), 'rb'))
	# action_count = pickle.load(open(outdir + 'mdp_actions_%s_%s.p'%(filename, LOT_COUNT), 'rb'))
	# state_count = len(state_space)
	# print (state_count)

	continue

	for curr in state_space:
		curr.opt_access_control = deepcopy(opt_actions[curr.idx])
		curr.get_full_access_control()
		curr.get_partial_access_control()
		# if curr.full_access_control not in curr.actions:
		# 	full = deepcopy(curr.full_access_control)
		# 	if full[0] == 0 or spot2ass( full[0], control = 'mdp' ) <= 0:
		# 		import pdb; pdb.set_trace()
		# if curr.partial_access_control not in curr.actions:
		# 	partial = deepcopy(curr.partial_access_control)
		# 	if partial[0] == 0 or spot2ass( partial[0], control = 'mdp' ) <= 0:
		# 		import pdb; pdb.set_trace()

	fullMatData, fullRow, fullCol, full_subspace = generate_subspace_control(state_space, 'full')
	partialMatData, partialRow, partialCol, partial_subspace = generate_subspace_control(state_space, 'partial')
	print (state_count, len(full_subspace), len(partial_subspace))

	import pdb; pdb.set_trace()

	for curr in state_space:
		print (curr.idx, curr.x, curr.y, curr.actions)
		print (curr.opt_access_control, curr.full_access_control, curr.partial_access_control)
		if curr.full_access_control not in curr.actions and curr.idx in full_subspace:
			full = deepcopy(curr.full_access_control)
			import pdb; pdb.set_trace()
			if full[0] > 0 and spot2ass( full[0], control = 'mdp' ) <= 0:
				import pdb; pdb.set_trace()
		if curr.partial_access_control not in curr.actions and curr.idx in partial_subspace:
			partial = deepcopy(curr.partial_access_control)
			import pdb; pdb.set_trace()
			if partial[0] > 0 and spot2ass( partial[0], control = 'mdp' ) <= 0:
				import pdb; pdb.set_trace()

	import pdb; pdb.set_trace()

	pickle.dump(state_space, open(outdir + 'mdp_state_list_%s_%s.p'%(filename, LOT_COUNT), 'wb'))



''' log files for full access control
printing params for double sided layout ...
angle: 0
mode: long
type of policy being evaluated: full
loading spot configuration: m_{out}_4_m_{in}_2_nUnit_1
input to distributions: meanSERV_60.0_meanDRIV_2.215909090909091
N: 1
Time to generate all states and actions: 0.0034
Number of states: 20
N: 2
Time to generate all states and actions: 0.0731
Number of states: 272
N: 3
Time to generate all states and actions: 1.4192
Number of states: 3072
N: 4
Time to generate all states and actions: 32.2996
Number of states: 24320
N: 5
^CTraceback (most recent call last):
  File "/Users/cindy/Downloads/2025 PUDO/0 - MDP code/runfile.py", line 58, in <module>
    (state_list, state_dict) = generate_policy_states(LOT_COUNT, control = policy_type)
  File "/Users/cindy/Downloads/2025 PUDO/0 - MDP code/markov.py", line 288, in generate_policy_states
    curr.data.lambda_next = nex.idx
  File "/Users/cindy/Downloads/2025 PUDO/0 - MDP code/SLinkedList.py", line 91, in addEnd
    curr = curr.next
KeyboardInterrupt


########################## SS0L - FULL ###########################
loading spot configuration: m_{out}_4_m_{in}_2_nUnit_1
input to distributions: meanSERV_60.0_meanDRIV_2.215909090909091

N: 1
Time to generate all states and actions: 0.0009
Number of states: 6
Time taken by the value iteration: 0.0056

N: 2
Time to generate all states and actions: 0.0042
Number of states: 28
Time taken by the value iteration: 0.0448

N: 3
Time to generate all states and actions: 0.027
Number of states: 120
Time taken by the value iteration: 0.361

N: 4
Time to generate all states and actions: 0.1245
Number of states: 422
Time taken by the value iteration: 2.104

N: 5
Time to generate all states and actions: 0.7476
Number of states: 1800
Time taken by the value iteration: 13.2836

N: 6
Time to generate all states and actions: 5.7628
Number of states: 8390
1000
Time taken by the value iteration: 87.0787

N: 7
Time to generate all states and actions: 71.2817
Number of states: 40068
1000
Time taken by the value iteration: 558.0241


########################## SS0S - FULL ###########################
loading spot configuration: m_{out}_4_m_{in}_10_nUnit_1
input to distributions: meanSERV_60.0_meanDRIV_1.7045454545454546

N: 1
Time to generate all states and actions: 0.0019
Number of states: 16
Time taken by the value iteration: 0.0262

N: 2
Time to generate all states and actions: 0.0121
Number of states: 82
Time taken by the value iteration: 0.255

N: 3
Time to generate all states and actions: 0.1022
Number of states: 494
1000
Time taken by the value iteration: 2.5988

N: 4
Time to generate all states and actions: 0.9067
Number of states: 2792
1000
Time taken by the value iteration: 21.4928

N: 5
Time to generate all states and actions: 11.6384
Number of states: 15217
1000
Time taken by the value iteration: 171.186

########################## SS90 - FULL ###########################
loading spot configuration: m_{out}_3_m_{in}_5_nUnit_3
input to distributions: meanSERV_60.0_meanDRIV_2.0454545454545454

N: 1
Time to generate all states and actions: 0.0014
Number of states: 10
Time taken by the value iteration: 0.0127

N: 2
Time to generate all states and actions: 0.0055
Number of states: 36
Time taken by the value iteration: 0.1072

N: 3
Time to generate all states and actions: 0.0225
Number of states: 104
Time taken by the value iteration: 0.5628

N: 4
Time to generate all states and actions: 0.1118
Number of states: 400
Time taken by the value iteration: 3.0251

N: 5
Time to generate all states and actions: 0.5808
Number of states: 1488
Time taken by the value iteration: 16.0646

N: 6
Time to generate all states and actions: 2.894
Number of states: 4992
Time taken by the value iteration: 71.2709

N: 7
Time to generate all states and actions: 24.3071
Number of states: 20112

########################## DS0L - FULL ###########################
loading spot configuration: m_{out}_4_m_{in}_2_nUnit_1
input to distributions: meanSERV_60.0_meanDRIV_2.215909090909091

N: 1
Time to generate all states and actions: 0.0035
Number of states: 20
Time taken by the value iteration: 0.0361

N: 2
Time to generate all states and actions: 0.0726
Number of states: 272
Time taken by the value iteration: 1.219

N: 3
Time to generate all states and actions: 1.4251
Number of states: 3072
Time taken by the value iteration: 27.9119

N: 4
Time to generate all states and actions: 30.8243
Number of states: 24320
Time taken by the value iteration: 375.1169

########################## DS0S - FULL ###########################
loading spot configuration: m_{out}_4_m_{in}_10_nUnit_1
input to distributions: meanSERV_60.0_meanDRIV_1.7045454545454546

N: 1
Time to generate all states and actions: 0.0102
Number of states: 60
Time taken by the value iteration: 0.2848

N: 2
Time to generate all states and actions: 0.3115
Number of states: 968
Time taken by the value iteration: 9.3748

N: 3
Time to generate all states and actions: 20.8856
Number of states: 20046
Time taken by the value iteration: 303.7198


########################## DS90 - FULL ###########################
loading spot configuration: m_{out}_3_m_{in}_5_nUnit_3
input to distributions: meanSERV_60.0_meanDRIV_2.0454545454545454

N: 1
Time to generate all states and actions: 0.0056
Number of states: 36
Time taken by the value iteration: 0.1034

N: 2
Time to generate all states and actions: 0.071
Number of states: 272
Time taken by the value iteration: 2.1628

N: 3
Time to generate all states and actions: 0.6735
Number of states: 1600
Time taken by the value iteration: 22.8675

N: 4
Time to generate all states and actions: 17.2279
Number of states: 16768
Time taken by the value iteration: 355.2056
'''
