########################################################################################################################################################
# Last updated: Oct 10, 2025 by Xinyu Liu
# This file defines the state objects in the Markov Decision Process model.
# The state object has two main methods for the mu- and the lambda-transitions, and a few helper functions.
# This file also contains two test functions for the state object: One for random sample path and one for a specifed sample path.
 
########################################################################################################################################################
import sys
import random
import pickle
import itertools
import collections 
import numpy as np
import pandas as pd

from heapq import *
from collections import *
from copy import deepcopy
from itertools import product
from math import ceil, floor, sqrt, exp

from utils import *
from params import *
from inputSetDef import *

class stateObject():

	def __init__(self, N, seedDRIV = None, seedSERV = None):

		if side == 'double':
			self.half_N = N
			self.N = 2 * N
		else:
			assert side == 'single'
			self.N = N

		self.n = spot2blk(self.N)
		if (angle == 0 and mode == 'short') or angle == 90:
			self.n += 1
		self.x = [1 for _ in range(self.N)]
		self.y = [0 for _ in range(self.n)]
		self._var = None
		# self._suffix = None
		
		self.actions = []
		self.spot_actions = None
		self.access_actions = None
		self.action_count = 0
		self.prob_mu = [mu(j) for j in range(1, self.N+1)]
		self.prob_lambda = 0

		self.idx = 0
		self.mu_next = dict()
		self.lambda_next = []
		self.prevs = []
		self.reward = 0.

		self.C_out = [None for _ in range(self.N)]
		self.C_in  = [None for _ in range(self.N)]
		self.C_f   = [None for _ in range(self.n)]
		self.J_in  = None
		self.I_b   = None

	@property
	def var(self):
		self.evaluate()
		# return ''.join(map(str, self._var)) + self._suffix
		return ''.join(map(str, self._var))

	def evaluate(self):

		self._var = []
		self._var.append( str(self.x.count(0)).zfill(2) )
		self._var.append( str(self.x.count(1)).zfill(2) )
		self._var.append( str(self.x.count(2)).zfill(2) )
		self._var.append( str(self.y.count(0)).zfill(2) )
		self._var.append( str(self.y.count(self.N + 2)).zfill(2) )
		self._var.append( str(self.y.count(self.N + 1)).zfill(2) )

		type_2 = []
		type_other = []
		for (j, val) in enumerate(self.x):
			if val == 1:
				self._var.append( str(j).zfill(2) )
			elif val == 2:
				type_2.append( str(j).zfill(2) )
			elif val > 0:
				type_other.append( str(j).zfill(2) )
				type_other.append( str(val).zfill(2) )
		self._var += type_2
		self._var += type_other

		type_2 = []
		type_other = []
		for (i, val) in enumerate(self.y):
			if val == self.N + 2:
				self._var.append( str(i).zfill(2) )
			elif val == self.N + 1:
				type_2.append( str(i).zfill(2) )
			elif val > 0:
				type_other.append( str(i).zfill(2) )
				type_other.append( str(val).zfill(2) )
		self._var += type_2
		self._var += type_other

	def evaluate2(self):

		self._var = []
		self._suffix = ''
		for j in range(self.N):
			if self.x[j] <= 2:
				self._var.append(self.x[j])
			else:
				self._var.append(3)
				self._suffix += str(self.x[j]).zfill(2)
		for i in range(self.n):
			if self.y[i] == 0:
				self._var.append(0)
			elif self.y[i] == self.N + 1:
				self._var.append(1)
			elif self.y[i] == self.N + 2:
				self._var.append(2)
			else:
				self._var.append(3)
				self._suffix += str(self.y[i]).zfill(2)

	def __eq__(self, other):
		return self.var == other.var

	def __hash__(self):
		return int(self.var)

	def mu_transition(self, j):

		assert j <= self.N
		assert self.x[j-1] == 1
		self.x[j-1] = 2
		self.prob_mu[j-1] = 0.
		if self.prob_lambda == .0:
			self.prob_lambda = rateDRIV
		return

	# Algorithm 6 Computation of (x',y') = f_{\lambda}((x,y),a) in the MDP
	def lambda_transition(self, a, control = 'mdp'):
		
		# Step 0: Initialize x and z
		x = deepcopy(self.x)
		y = deepcopy(self.y)
		self.y = [0 for _ in range(self.n)]
		
		a0 = a[0]
		z = deepcopy(y)
		# update state vector y with the selected spot assignment decisions
		for (i, ai) in enumerate(a[1:]):
			if ai > 0:
				assert z[i] == self.N+2
				z[i] = ai

		# Pre-calculations and simulation recorders
		if control == 'full' or control == 'partial':
			self.J_in, self.I_b = self.schedule_enter_and_forward(x, z)
		if control == 'full' and len(self.I_b) > 0:
			import pdb; pdb.set_trace()
			
		self.effec_arrival = False
		self.departure = False
		self.start_service = []

		for j in range(self.N, 0, -1):

			# Step 1: Vehicles that completed some but not all steps of their exit or enter maneeuvers
			#         complete one more step, i.e.,
			if 1 <= x[j-1]-2 <= m_out(j)-2:
				# i.e. if the vehicle in spot j has completed the first but not all steps of exit maneuver
				assert self.x[j-1] == x[j-1]
				self.x[j-1] = x[j-1] + 1

			elif x[j-1]-2 == m_out(j)-1:
				# i.e. if the vehicle in spot j has completed all but the last step of exit maneuver
				assert self.x[j-1] == x[j-1]
				self.x[j-1] = 0
				
				if spot2out(j) <= self.n:
					# i.e. vehicle from spot j will be in lane block i(j)+ upon completing the last step
					assert self.y[spot2out(j)-1] == 0
					self.y[spot2out(j)-1] = self.N + 1 
				else: 	
					# i.e. vehicle in spot j leaves the facility with its last step of exit maneuver
					self.departure = True
					assert self.reward > 0.
					self.reward = 0.

			elif 1 <= x[j-1] - (m_out(j)+1) <= m_in(j)-2:
				# i.e. if the vehicle in spot j has completed the first but not all steps of enter maneuver
				assert self.x[j-1] == x[j-1]
				self.x[j-1] = x[j-1] + 1

			elif x[j-1] - (m_out(j)+1) == m_in(j)-1:
				# i.e. if the vehicle in spot j has completed all but the last step of enter maneuver
				assert self.x[j-1] == x[j-1]
				self.x[j-1] = 1
				self.prob_mu[j-1] = rateSERV
				self.start_service.append(j)

		# Step 2: A vehicle ready to leave the facility may proceed to leave
		#         Other vehicles in the driving lane either start an enter maneuver,
		#         or move forward unless blocked by the downstream traffic, i.e.,
		for i in range(self.n, 0, -1):

			if z[i-1] in self.J_in:
				assert self.x[z[i-1]-1] == 0
				if m_in(j) == 1:
					self.x[z[i-1]-1] = 1
				else:
					self.x[z[i-1]-1] = m_out(z[i-1]) + 2
			
			elif 1 <= z[i-1] <= self.N and i == spot2in(z[i-1]):
				assert self.y[i-1] == 0
				self.y[i-1] = z[i-1]
			
			elif i in self.I_b:
				if not self.y[i-1] == 0:
					print (x, y, self.x, self.y, self.idx)
				self.y[i-1] = z[i-1]

			elif (1 <= z[i-1] <= self.N and i < spot2in(z[i-1])) or (z[i-1] == self.N+1 and blk2blk(i) <= self.n) or (z[i-1] == self.N+2):
				assert self.y[blk2blk(i)-1] == 0
				self.y[blk2blk(i)-1] = z[i-1]
		
			elif (z[i-1] == self.N+1 and blk2blk(i) > self.n):
				self.departure = True
				assert self.reward > 0.
				self.reward = 0.

		# Step 3: A replacement vehicle may enter the facility
		if 1 <= a0 <= self.N:	
			if spot2in(a0) == 0:
				assert self.x[a0-1] == 0

				if m_in(a0) == 1:
					self.x[a0-1] = 1
					self.prob_mu[a0-1] = rateSERV
					self.start_service.append(a0)
				else:
					self.x[a0-1] = m_out(a0)+2
			else:
				assert self.y[blk2blk(0)-1] == 0
				self.y[blk2blk(0)-1] = a0 

			self.effec_arrival = True

		elif a0 == self.N+1:
			assert self.y[blk2blk(0)-1] == 0
			self.y[blk2blk(0)-1] = self.N+2

			self.effec_arrival = True

		# Step 4: Vehicles that completed service but have not yet started to exit will, 
		#         if possible, complete the first step of their exit maneuvers
		if control == 'full':
			self.schedule_exit_full(x, z, a0)
		else:
			scheduled_exit = set([])
			if 1 <= a0 <= self.N and x[a0-1] == 2:
				if self.get_Cout(x, z, a0) > 0:
					if z[0] > 0 and 0 in self.I_b:
						import pdb; pdb.set_trace()
				else:
					scheduled_exit.add(a0)
					assert self.x[a0-1] == 2
					if m_out(a0) == 1:
						self.x[a0-1] = 0
						if spot2out(a0) <= self.n:
							# i.e. vehicle from spot a0 will be in lane block i_out(a0) upon completing the last step
							assert self.y[spot2out(a0)-1] == 0
							self.y[spot2out(a0)-1] = self.N+1 	
						else: # i.e. vehicle in spot a0 leaves the facility with its last step of exit maneuver
							assert self.reward > 0.
							self.reward = 0.
							self.departure = True
					else:
						self.x[a0-1] = 3

			for j in range(self.N, 0, -1):
				if j in scheduled_exit:
					continue

				if x[j-1] == 2 and self.get_Cout(x, z, j) == 0 and np.sum([in_Koo(j, j_, 0) for j_ in scheduled_exit]) == 0:
					if (spot2wo(j, 0) > 0 or self.y[blk2blk(0)-1] == 0) and (a0 in [0, j, self.N+1] or spot2in(a0) > 0 or not in_Koi(a0, j, 0)):
						scheduled_exit.add(j)
						assert self.x[j-1] == 2
						if m_out(j) == 1:
							self.x[j-1] = 0
							if spot2out(j) <= self.n:
								# i.e. vehicle from spot j will be in lane block i_out(j) upon completing the last step
								assert self.y[spot2out(j)-1] == 0
								self.y[spot2out(j)-1] = self.N+1 				
							else: # i.e. vehicle in spot j leaves the facility with its last step of exit maneuver
								assert self.reward > 0.
								self.reward = 0.
								self.departure = True
						else:
							self.x[j-1] = 3

		# Last Step: Update the reward for the new state 
		assert self.reward == 0.
		if np.abs( np.sum( self.prob_mu ) - np.sum([mu(j) for j in self.x if j == 1]) ) > 1e-12:
			import pdb; pdb.set_trace()
		self.prob_lambda = int( ( sum(self.y) > 0) | (self.x.count(1) < self.N) ) * rateDRIV
		
		if self.y[-1] == self.N + 1:
			self.reward = rateDRIV
			return

		for j in range(1, self.N+1):
			if self.x[j-1] == m_out(j)+1 and spot2out(j) > self.n:
				self.reward = rateDRIV
				break
		return

	def schedule_exit_full(self, x, z, a0):

		scheduled_exit = set()
		if a0 > 0 and x[a0-1] == 2:
			scheduled_exit.add(a0)

		for j in range(self.N, 0, -1):
			if j not in scheduled_exit and x[j-1] == 2 and self.step_schedule_exit_full(x, z, j, a0, scheduled_exit):
				scheduled_exit.add(j)

		for j in scheduled_exit:
			assert self.x[j-1] == 2
	
			if m_out(j) == 1:
				self.x[j-1] = 0
				if spot2out(j) <= self.n:
					# i.e. vehicle from spot j will be in lane block i_out(j) upon completing the last step
					assert self.y[spot2out(j)-1] == 0
					self.y[spot2out(j)-1] = self.N+1 				
				else: # i.e. vehicle in spot j leaves the facility with its last step of exit maneuver
					assert self.reward > 0.
					self.reward = 0.
					self.departure = True

			else:
				self.x[j-1] = 3

		return

	def step_schedule_exit_full(self, x, z, j, a0, scheduled_exit):
		
		for j_ in range(1, self.N+1):
			if j_ != j and 1 <= x[j_-1]-(m_out(j_)+1) <= m_in(j_)+1 and in_Ein(j, j_, x[j_-1]-(m_out(j_)+1)):
				return False
			if j_ != j and 1 <= x[j_-1]-2 <= m_out(j_)+1 and in_Eout(j, j_, x[j_-1]-2):
				return False
			if j_ != a0 and np.sum([z[i-1] == j_ for i in range_Iout(j, j_, self.n) if i > 0]):
				return False
		if np.sum([z[i-1] == self.N+1 for i in range_Iout(j, self.N+1, self.n) if i > 0]) > 0:
			return False
		if np.sum([in_Eout(j, j_, 0) for j_ in scheduled_exit]) > 0:
			return False
		if a0 > 0 and in_Iout(j, a0, 0, self.n):
			return False

		return True

	# Algorithm 1 + the overall action generation
	def get_mdp_actions(self):

		# housekeeping
		x = deepcopy(self.x)
		y = deepcopy(self.y)
		spot_actions = []

		# Algorithm 1 Determining the sets A_i(s) of feasible spot assignment decisions.
		
		for i in range(1, self.n + 1):
			
			# if y_i = 0 or 1 <= y_i <= N or y_i = N+1, 
			# then no assignment can be done 
			if y[i-1] <= self.N + 1:
				spot_actions.append( [0] )

			# else y_i has a vehicle without a spot assignment yet
			else:
				assert y[i-1] == self.N + 2

				if len(assignable(i, self.N, 'mdp')) == 0:
					spot_actions.append( [0] )
					continue

				j_star = max( assignable(i, self.N, 'mdp') )  # relative to block i
				
				# Count the available spots downstream of spot j_star
				S_count = np.sum([ ((x[j-1] == 0) or (2 <= x[j-1] <= m_out(j) + 1)) for j in range(j_star+1, self.N+1)])
				# Count the vehicles on the driving lane that are or have to be assgined to the available spots downstream of j_star
				V_count = np.sum([ ((j_star < y[i_-1] <= self.N) or (i_ >= blk2blk(i) and y[i_-1] == self.N+2)) for i_ in range(1, self.n+1)])				
				if S_count < V_count:
					print ('Error Warning: More replacement vehicles downstream than required!')
					import pdb; pdb.set_trace()

				# Identify the available spots to be assigned, i.e., A_{i}(s)
				avail_spots = []
				for j in assignable(i, self.N, 'mdp'):
					check_assigned = False
					for i_ in range(1, self.n+1):
						if y[i_-1] == j:
							check_assigned = True

					if ((x[j-1] == 0) or (2 <= x[j-1] <= m_out(j) + 1)) and not check_assigned:
						avail_spots.append(j)

				if S_count == V_count:
					for j in avail_spots:
						if spot2ass(j, 'mdp') != i:
							avail_spots.remove(j)
				else:
					assert S_count > V_count
					avail_spots.append(0)

				spot_actions.append( avail_spots )

		# generate the set of all feasible decisions 
		self.actions = []
		self.access_actions = []
		self.action_count = 0
		self.spot_actions = [a1 for a1 in product(*spot_actions)]

		for a1 in product(*spot_actions):

			# check for single spot assigned to multiple vehicles in the driving lane
			VALID = True
			assigned = set([])
			for ai in a1:
				if ai in assigned:
					VALID = False
					break
				if ai > 0:
					assigned.add(ai)

			if not VALID:
				self.spot_actions.remove(a1)

			else:
				# update state vector y with the selected spot assignment decisions
				z = deepcopy(y)
				for (i, ai) in enumerate(a1):
					if ai > 0:
						assert z[i] == self.N+2
						z[i] = ai

				# generate access control decisions for each feasible combination of spot assignment decisions
				self.access_actions.append( self.get_access_control_actions(x, z, 'mdp') )
				for a0 in self.access_actions[-1]:
					if np.sum(self.y) + np.sum( np.array(self.x) >= 2) + a0 > 0:
						self.actions.append( [a0] + list(a1) )
						self.action_count += 1

		assert len(self.spot_actions) == len(self.access_actions)
		return
	
	# Algorithm 2 Determining scheduled enter maneuvers J^in and blocked forward movements I^b
	def schedule_enter_and_forward(self, x, z):

		# reset placeholders
		self.C_out = [None for _ in range(self.N)]
		self.C_in  = [None for _ in range(self.N)]
		self.C_f   = [None for _ in range(self.n)]

		# Step 0: Initialization
		J_in = set([])
		I_b  = set([])

		for i in range(self.n, 0, -1):

			# Step 1: Computation of C^{in}(j) if z_i == j
			if 1 <= z[i-1] <= self.N and i == spot2in(z[i-1]):

				assert self.C_in[z[i-1]-1] is None
				if self.step_schedule_enter(x, z, J_in, I_b, z[i-1], i):
					self.C_in[z[i-1]-1] = 0
					J_in.add(z[i-1])
				else:
					self.C_in[z[i-1]-1] = 1

			# Step 2: Computation of C^f(i)
			elif z[i-1] > 0:
				assert self.C_f[i-1] is None
				if blk2blk(i) in I_b or self.step_schedule_blocked(x, z, J_in, i):
					self.C_f[i-1] = 1
					I_b.add(i)
				else:
					self.C_f[i-1] = 0
		
		return J_in, I_b
			
	# helper function for Algorithm 2 --- Computation of C^{in}(j)
	def step_schedule_enter(self, x, z, J_in, I_b, j, i):
		
		if i == 0:
			assert spot2in(j) == 0
		else:
			assert j == z[i-1] and 1 <= j <= self.N and i == spot2in(j)

		for j_ in range(1, self.N+1):

			# Step 1(a): Count the number of ongoing exit maneuvers that block the enter maneuver into spot j
			if x[j_-1] != 2 and in_Kio(j, j_, x[j_-1]-2):
				return False
			if x[j-1] == 2 and in_Kio(j, j, 0):
				return False

			# Step 1(b): Count the number of ongoing enter maneuvers that block the enter maneuver into spot j
			if j_ != j and x[j_-1] != m_out(j_)+1 and in_Kii(j, j_, x[j_-1]-(m_out(j_)+1)):
				return False

		# Step 1(c): Count the number of scheduled downstream enter maneuvers that block the enter maneuver into spot j
		for i_ in range(blk2blk(i), self.n+1):
			if j in J_in and in_Kii(j, z[i_-1], 0):
				return False

		# Step 1(d): Count the number of vehicles in blocks i_ \in {i^+, \cdots, i(j)} staying in 
		#            blocks that are required by the enter maneuver into spot j by Rule A.1.
		for i_ in range(blk2blk(i), spot2blk(j)):
			if z[i_-1] > 0 and z[i_-1] not in J_in:
				return False
		if 1 <= z[spot2blk(j)-1] <= self.N and spot2blk(j) == spot2in(z[spot2blk(j)-1]) and z[spot2blk(j)-1] not in J_in:
			return False
		if spot2blk(j) in I_b:
			return False

		return True

	# helper function for Algorithm 2 --- Computation of C^f(i)
	def step_schedule_blocked(self, x, z, J_in, i):

		assert i == 0 or z[i-1] > 0
		assert i == 0 or z[i-1] > self.N or i < spot2in(z[i-1])

		# needs to be revisited
		try:
			if i > 0 and z[i-1] <= self.N and x[z[i-1]-1] == 2 and i == spot2wo(z[i-1], 0):
				return True
		except IndexError:
			import pdb; pdb.set_trace()

		for j in range(1, self.N+1):

			# Step 2(a): Count the number of ongoing exit maneuvers that require the forward-moving vehicle to stay in block i
			if 3 <= x[j-1] <= m_out(j)+1 and i == spot2wo(j, x[j-1]-2):
				return True
			# Step 2(b): Count the number of ongoing enter maneuvers that require the forward-moving vehicle to stay in block i
			if m_out(j)+2 <= x[j-1] <= m_in(j) + m_out(j) and i == spot2wi(j, x[j-1]-(m_out(j) + 1)):
				return True
		# Step 2(c): Count the number of vehicles that have been assigned spots (including the scheduled enter maneuvers) 
		#                                     and that require the forward-moving vehicle to stay in block i
		for i_ in range(i+1, self.n+1):
			if z[i_-1] in J_in and i == spot2wi(z[i_-1], 0):
				return True
			if 1 <= z[i_-1] <= self.N and z[i_-1] not in J_in and i == spot2wi(z[i_-1], -1):
				return True
		
		return False

	# Algorithm 4 Computation of \tau_{s,a_{-0}}(j)
	def compute_tau(self, x, z, j):

		# Initialization
		I_wi = spot2wo(j,0)
		T_wi = 0

		for i in range(spot2wo(j,0), 0, -1):
			# Scenario 1: if block i has a vehicle that will go pass i_{wo}(j,0)
			if (1 <= z[i-1] <= self.N and spot2in(z[i-1]) > spot2wo(j,0) ) or z[i-1] >= self.N+1:
				return max(I_wi - i, T_wi) + spot2wo(j,0) - I_wi

			# in all other scenarios, block i either has no vehicle
			# or it has a vehicle that will start an enter maneuver before reaching i_{wo}(j,0)

			# Scenario 2: otherwise block i has a vehicle that will go pass I_wi
			if 1 <= z[i-1] <= self.N and spot2in(z[i-1]) > I_wi:
				T_wi = max(I_wi - i, T_wi) + spot2in(z[i-1]) - I_wi + m_in(z[i-1])
				I_wi = spot2wi( z[i-1], m_in(z[i-1])-1 )

			# Scenario 3: otherwise block i has a vehicle that will not go pass I_wi	
			elif 1 <= z[i-1] <= self.N and spot2in(z[i-1]) <= I_wi:
				# if it will impose additional delay
				if spot2in(z[i-1]) - i + m_in(z[i-1]) + I_wi - spot2wi(z[i-1], m_in(z[i-1])-1) > T_wi:
					T_wi = spot2in(z[i-1]) - i + m_in(z[i-1])
					I_wi = spot2wi(z[i-1], m_in(z[i-1])-1)
			
			# Scenario 4: block i is empty but it is close to an ongoing maneuver such that 
			#             any vehicle delayed by this ongoing maneuver will move to i when blocking is lifted
			else:
				for j_ in range(1, self.N+1):
					if 1 <= x[j_-1]-2 <= m_out(j_)-1 and i == blk2blk(spot2wo(j_, x[j_-1]-2)):
						return max( m_out(j_) - (x[j_-1]+2) + I_wi - spot2out(j_), T_wi ) + spot2wo(j,0) - I_wi

					if 1 <= x[j_-1]-(m_out(j_)+1) <= m_in(j_)-1 and i == blk2blk(spot2wi(j_, x[j_-1]-(m_out(j_)+1))):
						# if it will impose additional delay
						if m_in(j_) - (x[j_-1]-(m_out(j_)+1)) + I_wi - spot2wi(j_, m_in(j_)-1) > T_wi:
							T_wi = m_in(j_) - (x[j_-1]-(m_out(j_)+1))
							I_wi = spot2wi(j_, m_in(j_)-1)
						break
	
	# Algorithm 3 Determining C^{out}_{s,a_{-0}}(j)
	def compute_Cout(self, x, z, j):
		
		# Steps 1 - 4
		if self.step_compute_Cout(x, z, j) == 1:
			return 1

		# Step 5: If no vehicle in blocks i <= i_{wo}(j,0) needs to enter the conflict region
		ENTER_CONFICT = False
		for i in range(1, spot2wo(j,0)+1):
			if 1 <= z[i-1] <= self.N and spot2in(z[i-1]) > spot2wo(j,0):
				ENTER_CONFICT = True
			if z[i-1] >= self.N + 1:
				ENTER_CONFICT = True
		if not ENTER_CONFICT:
			return 0

		# Step 6: If the first vehicle in blocks i <= i_{wo}(j,0) that may conflict with the exit maneuver from spot j
		#                              is assigned to the same spot j
		for i in range(1, min(spot2in(j), spot2wo(j,0))+1):
			if z[i-1] == j:
				ENTER_CONFICT = False
				for i_ in range(i+1, spot2wo(j,0)+1):
					if 1 <= z[i_-1] <= self.N and spot2in(z[i_-1]) > spot2wo(j,0):
						ENTER_CONFICT = True
					if z[i_-1] >= self.N+1:
						ENTER_CONFICT = True
				if not ENTER_CONFICT:
					return 0

		# Step 7: calculate tau
		if self.compute_tau(x, z, j) >= g_out(j):
			return 0
		return 1

	# helper function for Algorithm 3
	def step_compute_Cout(self, x, z, j):

		for j_ in range(1, self.N+1):
			if j_ != j:
				
				# Step 1: Count the ongoing exit maneuvers that block the exit maneuver from spot j
				if x[j_-1] != 2 and in_Koo(j, j_, x[j_-1]-2):
					return 1

				# Step 2: Count the ongoing and scheduled enter maneuvers that block the exit maneuver from spot j
				if x[j_-1] != m_out(j_)+1 and in_Koi(j, j_, x[j_-1]-(m_out(j_)+1)):
					return 1
				if j_ in self.J_in and in_Koi(j, j_, 0): 
					return 1

		# Step 3: Count the downstream vehicles that have been assigned spots but not scheduled to complete 
		#               the first step of their enter maneuvers, which conflict with the exit maneuver from spot j 
		for i in range(blk2blk(spot2wo(j,0)), self.n+1):
			if 1 <= z[i-1] <= self.N and z[i-1] not in self.J_in and z[i-1] != j and spot2out(j) > spot2wi(z[i-1],-1):
				return 1

		# Step 4: Count the number of vehicles already in blocks as below, 
		#         staying in the range of blocks that are required by the exit maneuver from spot j
		bar_I = min(self.n, max(spot2blk(j), spot2out(j)))

		for i in range(blk2blk(spot2wo(j,0)), bar_I):
			if z[i-1] > 0 and z[i-1] not in self.J_in:
				return 1
		if 1 <= z[bar_I-1] <= self.N and z[bar_I-1] not in self.J_in and bar_I == spot2in(z[bar_I-1]):
			return 1
		if bar_I in self.I_b:
			return 1

		return 0

	# helper function to call C^{out} in Alg. 5 without repeated computation
	def get_Cout(self, x, z, j):
		
		if self.C_out[j-1] is None:
			self.C_out[j-1] = self.compute_Cout(x, z, j)
		return self.C_out[j-1]

	# Algorithm 5 Determining feasible access control decisions A_0(s,a_{-0})
	def get_access_control_actions(self, x, z, control):

		# Step 0: Initialize the pre-calculation
		self.J_in, self.I_b = self.schedule_enter_and_forward(x, z)

		# Step 1: Determine any available spots for a replacement vehicle.
		S_count = np.sum([ (x[j-1] == 0 or 2 <= x[j-1] <= m_out(j)+1) for j in range(1, self.N+1)])
		V_count = np.sum([ (1 <= z[i-1] <= self.N or z[i-1] == self.N+2) for i in range(1, self.n+1)])
		E_count = np.sum([ (spot2ass(j, control) == 0 and (x[j-1] == 0 or 2 <= x[j-1] <= m_out(j)+1) and 0 == np.sum(np.array(z) == j)) for j in range(1, self.N+1) ])

		if S_count == V_count:
			return [0]

		# Step 2: Determin any waiting at the facility entrance due to blocking in the driving lane.
		if blk2blk(0) in self.I_b:
			return [0]

		# Initialize some pre-calculated quantities 
		from_entrance = []
		from_general = []
		for j in assignable(0, self.N, control):
			if spot2in(j) == 0 and ((x[j-1] == 2 and self.get_Cout(x, z, j) == 0) or x[j-1] == 0 or 3 <= x[j-1] <= m_out(j)+1): 
				if self.step_schedule_enter(x, z, self.J_in, self.I_b, j, 0) and not in_Kio(j, j, x[j-1]-2):
					from_general.append(j)
					from_entrance.append(j)

			if spot2in(j) > 0 and (x[j-1] == 0 or (2 <= x[j-1] <= m_out(j)+1 and spot2wo(j,x[j-1]-2) > 0)) and 0 == np.sum(np.array(z) == j):
				from_general.append(j)

		# Step 3: Determine any waiting at the facility entrance due to other vehicles' maneuvers.
		#         such that a replacement vehicle must start an enter maneuver from entrance
		if self.step_schedule_blocked(x, z, self.J_in, 0):
			return [0] + from_entrance

		# Step 4: Determine the set of assignable spots in all other occasions.
		if S_count - E_count != V_count:
			return [0, self.N+1] + from_general

		return [0] + [j for j in from_general if spot2ass(j, control) == 0]

	# Helper function for full access control -- eq. (26)
	def step_full_access_control(self, j):

		if self.x[j-1] == 2 and not in_Kout(j,j,0):
			for j_ in range(1, self.N+2):
				if j_ != j and np.sum( [self.y[i-1] == j_ for i in range_Iout(j, j_, self.n) if i > 0]) > 0:
					return False
				if j_ != j and j_ <= self.N:
					if np.sum( [ self.y[i-1] == j_ for i in range_Iin(j, j_, self.n) ] ) > 0:
						return False
					if in_Ein(j, j_, self.x[j_-1]-(m_out(j_)+1)):
						return False
					if 3 <= self.x[j_-1] <= m_out(j_)+1 and in_Eout(j, j_, self.x[j_-1]-2):
						return False
					if in_Kin(j, j_, self.x[j_-1]-(m_out(j_)+1)):
						return False
					if in_Kout(j, j_, self.x[j_-1]-2):
						return False
			return True

		elif (self.x[j-1] == 0 or 1 <= self.x[j-1]-2 <= m_out(j)-1):
			if np.sum( np.array(self.y) == j ) > 0:
				return False
			for j_ in range(1, self.N+1):
				if j_ != j and np.sum( [self.y[i-1] == j_ for i in range_Iin(j, j_, self.n)] ) > 0:
					return False
				if j_ != j and in_Kin(j, j_, self.x[j_-1]-(m_out(j_)+1)):
					return False
				if in_Kout(j, j_, self.x[j_-1]-2):
					return False
			return True

		return False

	# Determining full access control
	def get_full_access_control(self):

		feas_access_actions = self.get_access_control_actions(self.x, self.y, 'full')
		
		j_star = None
		for j in sorted(feas_access_actions, reverse = True):
			if j > 0 and self.step_full_access_control(j):
				j_star  = j
				break
		if j_star is None:
			self.full_access_control = [0] + [0] * self.n
		else:
			self.full_access_control = [j_star] + [0] * self.n

		# j_star = None
		# for j in range(self.N, 0, -1):
		# 	if self.step_full_access_control(j):
		# 		j_star  = j
		# 		break
		# if j_star is None:
		# 	self.full_access_control = [0] + [0] * self.n
		# else:
		# 	self.full_access_control = [j_star] + [0] * self.n

		return 

	# Determining partial access control
	def get_partial_access_control(self):

		feas_access_actions = self.get_access_control_actions(self.x, self.y, 'partial')

		if len(feas_access_actions) >= 1:
			self.partial_access_control = [max(feas_access_actions)] + [0] * self.n
		else:
			self.partial_access_control = [0] + [0] * self.n

		return 

		# J_set = set([]) # self.get_access_control(self.x, self.y, 'partial')
		# if len(J_set) >= 1:
		# 	j_star = max(J_set)
		# else:
		# 	j_star = 0
		# self.partial_access_control = [j_star] + [0] * self.n
		# return

def simulationMC():

	def __init__():

		self.probRNG = ParamGen(Unif(2., seed = seedDRIV))
		self.timeRNG = ParamGen(Unif(2., seed = seedDRIV))

		self.curr = 0.0
		self.start_time = 0.0
		self.cycle_times = []
		self.eventheap = []
		self.debug = debug
		self.debug_unit = SIM_UNIT

		self.service_times = dict()
		for j in range(1, self.N + 1):
			self.service_times[j] = []
		self.out_times = []
		self.in_times  = []

	def add_event(self, event):
		heappush(self.eventheap, event)

	def run(self, state_space, control):

		test = state_space[0]
		assert self.curr == 0.0
		assert len(self.eventheap) == 0
		print ('starting properly ...')

		# initialize the single lambda transition under a specified control
		if test.prob_lambda > 0.:

			if control == 'mdp':
				nex_idx = test.opt_lambda_next
			elif control == 'full':
				nex_idx = test.full_lambda_next
			elif control == 'partial':
				nex_idx = test.partial_lambda_next
			else:
				print ('Unknown control policy!')
				import pdb; pdb.set_trace()
		
			nex = state_space[nex_idx]
			self.add_event( event(self.curr + self.timeDRIV.next(), nex_idx, 'lambd') )

		# initialize all possible mu transitions common to all control policies
		for (j, nex_idx) in test.mu_next.items():
			next_serv = self.timeSERV.next()
			self.add_event( event(self.curr + next_serv, (j, nex_idx), 'mu') )
			self.service_times[j].append(next_serv)

		self.run_iter()

	def run_iter(self, SIM_UNTIL):

		while self.curr - self.start_time <= SIM_UNTIL:
				
			curr_event = heappop(self.eventheap)
			# curr_time = curr_event.time
			self.curr = float(curr_time)
			curr_vehicle = curr_event.vehicle
			curr_typ = curr_event.typ

			################################### update system ###################################
			if curr_typ == 'mu':
				self.mu_transition(curr_vehicle)
				
			else:
				assert curr_typ == 'lambd'
				assert self in state_space
				self.lambda_transition(actions[state_space[self]])
				if self.effec_arrival:
					self.inCount += 1
					self.in_times.append(self.curr)
				if self.departure:
					self.outCount += 1
					self.out_times.append(self.curr)
					self.reward += (exp(- DISC_RATE * self.curr))
				if len(self.start_service) >= 1:
					for j in self.start_service:
						next_serv = self.timeSERV.next()
						self.add_event( event(self.curr + next_serv, j, 'mu') )
						self.service_times[j].append(next_serv)

				self.add_event( event(self.curr + self.timeDRIV.next(), None, 'lambd') )

		self.start_time = self.curr
		return (self.outCount)

###################################################

def TestRandomFull():
	# N = random.randint(8, 28)
	N = 8

	print ('\nTest Case (%s spaces)' %N)
		
	state = stateObject(N)
	ITE = 100

	print ('\n1')
	print (state.x)
	print (state.y)
	state.get_full_access_control()
	print ('full access control action:', state.full_access_control)

	# v = random.sample(list(range(1, N+1)), 1)[0]
	v = 2
	state.mu_transition(v)

	print ('\n2')
	print ('vehicle at %s finishing service ...' %v)
	print (state.x)
	print (state.y)
	state.get_full_access_control()
	print ('full access control action:', state.full_access_control)

	for ite in range(3, ITE+1):
		print ('\n%s' %ite)

		keys = []
		for j in range(len(state.x)):
			if state.x[j] == 1:
				keys.append(j + 1)

		if ((np.sum(state.y) > 0) | (state.x != np.ones(state.N)).any()) and (random.random() < rateDRIV / (state.x.count(1) * rateSERV + rateDRIV)):
			state.lambda_transition(state.full_access_control, 'full')
			print ('driving ...')
	
		elif keys != []:
			v = random.sample(keys, 1)[0]
			print ('vehicle at %s finishing service ...' %v)
			state.mu_transition(v)

		else:
			print ('The system is empty now!')
			import pdb; pdb.set_trace()

		print (state.x)
		print (state.y)
		state.get_full_access_control()
		print ('full access control action:', state.full_access_control)
		# import pdb; pdb.set_trace()

def TestRandomPartial():
	# N = random.randint(8, 28)
	N = 8

	print ('\nTest Case (%s spaces)' %N)
		
	state = stateObject(N)
	ITE = 100

	print ('\n1')
	print (state.x)
	print (state.y)
	state.get_partial_access_control()
	print ('partial access control action:', state.partial_access_control)

	# v = random.sample(list(range(1, N+1)), 1)[0]
	v = 2
	state.mu_transition(v)

	print ('\n2')
	print ('vehicle at %s finishing service ...' %v)
	print (state.x)
	print (state.y)
	state.get_partial_access_control()
	print ('partial access control action:', state.partial_access_control)

	for ite in range(3, ITE+1):
		print ('\n%s' %ite)

		keys = []
		for j in range(len(state.x)):
			if state.x[j] == 1:
				keys.append(j + 1)

		if ((np.sum(state.y) > 0) | (state.x != np.ones(state.N)).any()) and (random.random() < rateDRIV / (state.x.count(1) * rateSERV + rateDRIV)):
			state.lambda_transition(state.partial_access_control, 'partial')
			print ('driving ...')
	
		elif keys != []:
			v = random.sample(keys, 1)[0]
			print ('vehicle at %s finishing service ...' %v)
			state.mu_transition(v)

		else:
			print ('The system is empty now!')
			import pdb; pdb.set_trace()

		print (state.x)
		print (state.y)
		state.get_partial_access_control()
		print ('partial access control action:', state.partial_access_control)
		# import pdb; pdb.set_trace()

def TestRandom(policy_type):

	if policy_type == 'full':
		TestRandomFull()
		return
	elif policy_type == 'partial':
		TestRandomPartial()
		return

	# N = random.randint(8, 28)
	N = 8

	print ('\nTest Case (%s spaces)' %N)
		
	state = stateObject(N)
	ITE = 100

	print ('\n1')
	print (state.x)
	print (state.y)
	state.get_mdp_actions()
	print ('all feasible actions:', state.actions)

	# v = random.sample(list(range(1, N+1)), 1)[0]
	v = 2
	state.mu_transition(v)

	print ('\n2')
	print ('vehicle at %s finishing service ...' %v)
	print (state.x)
	print (state.y)
	state.get_mdp_actions()
	print ('all feasible actions:', state.actions)

	for ite in range(3, ITE+1):
		print ('\n%s' %ite)

		keys = []
		for j in range(len(state.x)):
			if state.x[j] == 1:
				keys.append(j + 1)

		if ((np.sum(state.y) > 0) | (state.x != np.ones(state.N)).any()) and (random.random() < rateDRIV / (state.x.count(1) * rateSERV + rateDRIV)):
			x = deepcopy(state.x)
			y = deepcopy(state.y)
			v = random.sample(list(range(len(state.actions))), 1)[0]
			print ('driving + selected action:', state.actions[v])
			state.lambda_transition(state.actions[v])

		elif keys != []:
			v = random.sample(keys, 1)[0]
			print ('vehicle at %s finishing service ...' %v)
			state.mu_transition(v)

		else:
			print ('The system is empty now!')
			import pdb; pdb.set_trace()

		print (state.x)
		print (state.y)
		state.get_mdp_actions()
		print ('all feasible actions:', state.actions)
		# import pdb; pdb.set_trace()

###################################################
if __name__ == "__main__":
	
	if len(sys.argv) > 1:
		policy_type = sys.argv[1]
	else:
		policy_type = 'mdp'

	TestRandom(policy_type)
