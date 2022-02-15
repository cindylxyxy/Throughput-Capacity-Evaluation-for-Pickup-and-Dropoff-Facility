########################################################################################################################################################
# Last updated: Feb 06, 2022 by Xinyu Liu
# This file defines the state objects in the Markov Decision Process model.
# The state object has two main methods for the mu- and the lambda-transitions, and a few helper functions.
# This file also contains two test functions for the state object: One for random sample path and one for a specifed sample path.
 
########################################################################################################################################################
from copy import deepcopy
from itertools import product
import random
import numpy as np
from params import *
from inputSetDef import *


class stateObject():

	def __init__(self, N):

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
		self.waiting = []
		# self.inCount = self.N
		# self.outCount = 0
		self._var = None
		# self._suffix = None
		self.idx = 0
		self.iter = 0
		self.vehi_count = 0
		self.spot_count = 0
		self.reward = 0.
		self.prob_mu = [rateSERV for _ in range(self.N)]
		self.prob_lambda = 0.
		self._actions = [ [0 for _ in range(self.n + 1)] ]
		self.mu_next = []
		self.lambda_next = []
		self.action_count = 0
		self.prevs = []

	@property
	def var(self):
		self.evaluate()
		# return ''.join(map(str, self._var)) + self._suffix
		return ''.join(map(str, self._var))

	@property
	def actions(self):
		return self._actions

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

	def get_actions(self):

		x = deepcopy(self.x)
		y = deepcopy(self.y)

		lane_actions = []
		for i in range(1, self.n + 1):
			if y[i-1] <= self.N + 1:
				lane_actions.append( [0] )
			else:
				assert y[i-1] == self.N + 2

				j = max( enterspot(i) )
				vehi_count = 0
				spot_count = 0

				for j_ in range(j + 1, self.N + 1):
					if (x[j_-1] == 0) or (2 <= x[j_-1] <= m_out + 1):
						spot_count += 1
				for i_ in range(i + 1, self.n + 1):
					if (y[i_-1] > j) and (y[i_-1] != self.N + 1):
						vehi_count += 1
				
				if spot_count < vehi_count:
					print ('More replacement vehicles downstream than required!')
					import pdb; pdb.set_trace()
				
				avail_spots = []
				for j in enterspot(i):
					if (j <= self.N) and ((x[j-1] == 0) or (2 <= x[j-1] <= m_out + 1)) and (j in self.waiting):
						avail_spots.append(j)
				
				if spot_count > vehi_count:
					avail_spots.append(0)
				lane_actions.append( avail_spots )

		self._actions = []
		self.action_count = 0
		for a_1 in product(*lane_actions):
			for a_0 in self.get_action_helper(x, y, a_1):
				self._actions.append( [a_0] + list(a_1) )
				self.action_count += 1

		return

	def get_action_helper(self, x, y, a):
		
		if self.spot_count < self.vehi_count:
			print ('More replacement vehicles in the facility than required!')
			import pdb; pdb.set_trace()

		if self.spot_count == self.vehi_count:
			return [0]

		for j in allspot(blk2blk(0)):
			if (j > self.N):
				assert angle == 90
			else:
				if (2 <= x[j-1] - m_out <= m_in) or (3 <= x[j-1] <= m_out):
					return [0]
				if (mode == 'short') or (angle == 90):
					if y[spot2in(j)-1] == j:
						return [0]

		if (mode == 'short') or (angle == 90):
			if g_in == 1 and y[blk2blk(0)-1] == self.N + 2 and a[blk2blk(0)] > 0:
				return [0]
			if spot2blk(y[blk2blk(0)-1]) == blk2blk(0):
				return [0]

		for j in range(max(allspot(blk2blk(blk2blk(0)))) + 1, self.N + 1):
			if 0 in y and y.index(0) >= spot2blk(j):
				if (2 <= x[j-1] - m_out <= m_in):
					return [0]
				if (3 <= x[j-1] <= m_out):
					return [0]
				if (g_in == 0) and (y[spot2in(j)-1] == self.N + 2) and (a[spot2in(j)] == j):
					return [0]
				if (y[spot2in(j)-1] == j) and (angle != 0 or mode != 'long' or (spot2blk(y[spot2in(j)-2]) != spot2blk(j) - 1)):
					return [0]
			else:
				break

		if (angle == 0) and (mode == 'long'):
			feas_actions = [0]
			if (x[1-1] == 0) and (1 in self.waiting):
				feas_actions.append(1)
			if (side == 'double') and (x[2-1] == 0) and (2 in self.waiting):
				feas_actions.append(2)

			if g_in == 0 and y[blk2blk(0)-1] == self.N + 2 and a[blk2blk(0)] > 0:
				return feas_actions
			if spot2in(y[blk2blk(0)-1]) == blk2blk(0):
				return feas_actions
			if side == 'single' and self.N >= 2 and 2 <= x[2-1] - m_out <= m_in - 1:
				return feas_actions
			if side == 'double' and self.N >= 3 and 2 <= x[3-1] - m_out <= m_in - 1:
				return feas_actions
			if side == 'double' and self.N >= 3 and 2 <= x[4-1] - m_out <= m_in - 1:
				return feas_actions

		feas_actions = [0]
		for j in range(1, min( [self.N + 1] + enterspot(blk2blk(0)) )):
			if j in self.waiting:
				if (x[j-1] == 0) or (3 <= x[j-1] <= m_out + 1): 
					feas_actions.append(j)
				if angle == 0 and mode == 'long' and spot2in(j) > 0 and x[j-1] == 2:
					feas_actions.append(j)
				if (angle == 90 or mode == 'short') and 1 < spot2blk(j) and x[j-1] == 2: 
					feas_actions.append(j)

		vehi_count = 0
		spot_count = 0

		for j in range(min( enterspot(blk2blk(0)) ), self.N + 1):
			if (x[j-1] == 0) or (2 <= x[j-1] <= m_out + 1):
				spot_count += 1

		for i in range(1, self.n + 1):
			if y[i-1] >= min(enterspot(blk2blk(0))) and y[i-1] != self.N + 1:
				vehi_count += 1
				
		if spot_count < vehi_count:
			print ('More replacement vehicles downstream than required!')
			import pdb; pdb.set_trace()

		if spot_count > vehi_count:
			feas_actions.append(self.N + 1) 
		return feas_actions	

	def mu_transition(self, j):

		assert j <= self.N
		assert self.x[j - 1] == 1
		self.x[j - 1] = 2
		self.spot_count += 1
		self.waiting.append(j)
		self.prob_mu[j - 1] = 0.
		if self.prob_lambda == .0:
			self.prob_lambda = rateDRIV
		return

	def lambda_transition(self, a):

		x = deepcopy(self.x)
		y = deepcopy(self.y)
		blockedlane = set()
		visitedlane = set()

		for j in range(self.N, 0, -1):

			# Step 1: vehicles that have completed some but not all steps of exit maneeuver
			#         complete their next step of exit maneuver
			# Step 2: vehicles that have completed some but not all steps of enter maneuver 
			#         complete their next stop of enter maneuver
			if (3 <= x[j - 1] <= m_out):
				# i.e. if the vehicle in spot j has completed the first but not all steps of exit maneuver
				self.x[j - 1] = x[j - 1] + 1
				blockedlane.add(spot2blk(j))
				blockedlane.add(spot2blk(j) + 1)

			elif (x[j - 1] == m_out + 1):
				# i.e. if the vehicle in spot j has completed (m_out - 1) steps of exit maneuver
				self.x[j - 1] = 0
				if spot2blk(j) == spot2blk(self.N):
					# i.e. vehicle in spot j leaves the facility with its last step of exit maneuver
					# self.outCount += 1
					assert self.reward > 0.
					self.reward = 0.
				else:
					# i.e. vehicle from spot j will be in lane block i(j) + 1 upon completing the last step
					# assert self.y[spot2blk(j)] == 0
					if spot2blk(j) + 1 in visitedlane:
						import pdb; pdb.set_trace()
					self.y[spot2blk(j)] = self.N + 1 	
					blockedlane.add(spot2blk(j) + 1)
					visitedlane.add(spot2blk(j) + 1)
		
			elif (2 + m_out <= x[j - 1] <= m_out + m_in - 1):
				# i.e. if the vehicle in spot j has completed the first but not all steps of enter maneuver
				self.x[j - 1] = x[j - 1] + 1
				blockedlane.add(spot2blk(j))
				blockedlane.add(spot2in(j))

			elif (x[j - 1] == m_out + m_in):
				# i.e. if the vehicle in spot j has completed (m_in - 1) steps of enter maneuver
				self.x[j - 1] = 1
				self.prob_mu[j - 1] = rateSERV
				blockedlane.add(spot2blk(j))

		# Step 3: vehicles move forward when possible
		for i in range(self.n, 0, -1):

			if y[i - 1] == 0:
				# i.e. if there is no vehicle in lane block i
				continue

			if self.y[i - 1] == y[i - 1] and i not in visitedlane:
				self.y[i - 1] = 0
			else:
				assert self.y[i - 1] == self.N + 1

			if y[i-1] == self.N + 2:

				if (blk2blk(i) in blockedlane):
					if not self.y[i - 1] == 0:
						import pdb; pdb.set_trace()
					self.y[i - 1] = y[i - 1]
					blockedlane.add(i)

				elif a[i] == 0:
					# i.e. the vehicle has not yet arrived at a lane block where it can start an enter maneuver
					assert self.y[blk2blk(i) - 1] == 0
					self.y[blk2blk(i) - 1] = y[i - 1]
					blockedlane.add(blk2blk(i))
				
				elif g_in >= 1:
					assert spot2in(a[i]) == i + g_in
					assert self.y[blk2blk(i) - 1] == 0
					self.y[blk2blk(i) - 1] = a[i]
					blockedlane.add(blk2blk(i))
					self.waiting.remove(a[i])

				else:
					assert spot2in(a[i]) == i 
					assert g_in == 0
					assert angle == 0 and mode == 'long'
					if self.x[a[i] - 1] == 0:
						assert self.x[a[i] - 1] == 0
						assert m_in > 1
						self.x[a[i] - 1] = m_out + 2
						blockedlane.add(i)
						blockedlane.add(spot2blk(y[i - 1]))
						self.spot_count -= 1
						self.vehi_count -= 1
						self.waiting.remove(a[i])
					elif self.x[a[i] - 1] == 2:
						assert self.y[i - 1] == 0
						self.y[i - 1] = a[i]
						blockedlane.add(i)
						self.waiting.remove(a[i])
					else:
						assert 3 <= x[a[i] - 1] <= m_out + 1
						assert self.y[i - 1] == 0
						self.y[i - 1] = y[i - 1]
						blockedlane.add(i)											

			elif (1 <= y[i - 1] <= self.N):

				if angle == 0 and mode == 'long':
					# i.e. if a vehicle enters or exits the boarding spot both by going forward 

					if (spot2in(y[i - 1]) == i) and (self.x[y[i - 1] - 1] == 0):
						# i.e. if a vehicle has arrived at the lane block i_in(j) = i(j) - 1, and is ready to start an enter maneuver,
						#      and there is no vehicle exiting from spot j 
						if (blk2blk(i) in blockedlane) and (1 <= y[blk2blk(i) - 1] <= self.N) and (self.x[y[blk2blk(i) - 1] - 1] != m_out + 2):
							assert self.y[i - 1] == 0
							self.y[i - 1] = y[i - 1]
							blockedlane.add(i)
							continue
						assert self.x[y[i - 1] - 1] == 0
						assert m_in > 1
						self.x[y[i - 1] - 1] = m_out + 2
						blockedlane.add(i)
						blockedlane.add(spot2blk(y[i - 1]))
						self.vehi_count -= 1
						self.spot_count -= 1

					elif spot2in(y[i - 1]) == i:
						# i.e. continuing from the previous case, if another vehicle is still exiting from spot j
						assert 2 <= x[y[i - 1] - 1] <= m_out + 1
						assert self.y[i - 1] == 0
						self.y[i - 1] = y[i - 1]
						blockedlane.add(i)

					else:
						# i.e. the vehicle has not yet arrived at a lane block where it can start an enter maneuver
						assert spot2in(y[i - 1]) > i
						if (blk2blk(i) in blockedlane):
							assert self.y[i - 1] == 0
							self.y[i - 1] = y[i - 1]
							blockedlane.add(i)
							continue
						assert self.y[blk2blk(i) - 1] == 0
						self.y[blk2blk(i) - 1] = y[i - 1]
						blockedlane.add(blk2blk(i))
				
				else:
					# i.e. if a vehicle enters a boarding spot by going backward and exits by going forward
					assert (angle == 0 and mode == 'short') or (angle == 90) 

					if (spot2blk(y[i - 1]) == i):
						# if a vehicle arrives at i(j) = i_in(j) - 1, it is not yet ready to start an enter maneuver
						# however, as it proceeds to (i(j))+, another vehicle should not proceed to i(j)
						if (blk2blk(i) in blockedlane):
							if angle == 90 and y[i - 1] + dgap + 1 <= self.N and 3 <= self.x[y[i - 1] + dgap] <= m_out + 1:
								pass
							elif angle == 90 and y[i - 1] % nUnit != 0 and y[i - 1] + dgap + 2 <= self.N and 3 <= self.x[y[i - 1] + dgap + 1] <= m_out + 1:
								pass
							elif angle == 90 and y[i - 1] % nUnit == 1 and y[i - 1] + dgap + 3 <= self.N and 3 <= self.x[y[i - 1] + dgap + 2] <= m_out + 1:
								pass
							else:
								assert self.y[i - 1] == 0
								self.y[i - 1] = y[i - 1]
								blockedlane.add(i)
								continue
						assert self.y[blk2blk(i) - 1] == 0
						self.y[blk2blk(i) - 1] = y[i - 1]
						blockedlane.add(i)
						blockedlane.add(blk2blk(i))					

					elif (spot2in(y[i - 1]) == i):
						# i.e. if a vehicle has arrived at the lane block i_in(j) = i(j) + 1, and is ready to start an enter maneuver
						# note that if the replacement vehicle has arrived at i(j) + 1, there should not be another vehicle exiting from j
						assert x[y[i - 1] - 1] == self.x[y[i - 1] - 1] == 0
						assert m_in > 1
						self.x[y[i - 1] - 1] = m_out + 2
						blockedlane.add(spot2blk(y[i - 1]))
						blockedlane.add(i)
						self.spot_count -= 1
						self.vehi_count -= 1

					elif (spot2blk(y[i - 1]) == blk2blk(i)):
						# i.e. if a vehicle has arrived at a lane block s.t. the block ahead is the one to start an enter maneuver 
						# if the previous vehicle in the same spot has not yet started its exit maneuver 
						# the replacement vheicle has to wait
						if (blk2blk(i) in blockedlane) or (x[y[i - 1] - 1] == 2):
							assert self.y[i - 1] == 0
							self.y[i - 1] = y[i - 1]
							blockedlane.add(i)
							continue
						assert self.y[blk2blk(i) - 1] == 0
						assert self.x[y[i - 1] - 1] == 0
						self.y[blk2blk(i) - 1] = y[i - 1]
						blockedlane.add(blk2blk(i)) 

					else:
						# i.e. the vehicle has not yet arrived at a lane block beside its assigned boarding spot
						assert spot2blk(y[i - 1]) > blk2blk(i)
						if (blk2blk(i) in blockedlane):
							assert self.y[i - 1] == 0
							self.y[i - 1] = y[i - 1]
							blockedlane.add(i)
							continue
						assert self.y[blk2blk(i) - 1] == 0
						self.y[blk2blk(i) - 1] = y[i - 1]
						blockedlane.add(blk2blk(i))					

			else:
				# i.e. if there is a vehicle in lane block i heading to the facility exit
				assert (y[i - 1] == self.N + 1)
				# self.y[i - 1] = 0
				# note that n = i(N) for 0-degree long spots
				#  and that n = i(N) + 1 for 0-degree short spots and 90-degree
				if (blk2blk(i) != spot2blk(self.N) + 1):
					if (blk2blk(i) in blockedlane):
						assert self.y[i - 1] == 0
						self.y[i - 1] = y[i - 1]
						blockedlane.add(i)
						continue
					assert self.y[blk2blk(i) - 1] == 0
					self.y[blk2blk(i) - 1] = self.N + 1
					blockedlane.add(blk2blk(i))
				else:
					assert (blk2blk(i) == spot2blk(self.N) + 1)
					# self.outCount += 1
					assert self.reward > 0.
					self.reward = 0.

		# Step 4: a replacement vehicle enter when a_0 > 0
		if a[0] > 0 and blk2blk(0) in blockedlane:
			if angle == 0 and mode == 'long' and spot2blk(y[blk2blk(0)-1]) == 2 and self.y[blk2blk(0)-1] == 0:
				pass
			else:
				a[0] = 0
	
		if a[0] == self.N + 1:
			assert self.y[blk2blk(0)-1] == 0
			assert blk2blk(0) not in blockedlane 
			self.y[blk2blk(0)-1] = self.N + 2
			# self.inCount += 1
			self.vehi_count += 1

		elif spot2in(a[0]) == 0:
			assert self.x[a[0]-1] == 0
			if m_in == 1:
				self.x[a[0]-1] = 1
				self.prob_mu[a[0]-1] = rateSERV
			else:
				self.x[a[0]-1] = 2 + m_out
			# self.inCount += 1
			self.spot_count -= 1
			self.waiting.remove(a[0])

		elif a[0] > 0:
			assert a[0] < min( enterspot(blk2blk(0)) )
			assert self.y[blk2blk(0)-1] == 0
			assert blk2blk(0) not in blockedlane 
			self.y[blk2blk(0)-1] = a[0]
			# self.inCount += 1
			self.vehi_count += 1
			self.waiting.remove(a[0])

		# Step 5: other vehicles that finished service but not the first step of exit maneuver
		# 		  complete the first step
		for j in range(self.N, 0, -1):
			if (x[j - 1] == 2) and (self.compl_service_helper(x, y, a, j)):
				assert (m_out > 1)
				self.x[j - 1] = 3
				blockedlane.add(spot2blk(j))
				blockedlane.add(spot2blk(j) + 1)

		# Last step: update the reward for the new state
		assert self.reward == 0.
		if np.abs( sum( self.prob_mu ) - self.x.count(1) * rateSERV ) > 1e-12:
			import pdb; pdb.set_trace()
		self.prob_lambda = int( ( sum(self.y) > 0) | (self.x.count(1) < self.N) ) * rateDRIV
		
		if angle == 0 and mode == 'long' and self.y[-1] == self.N + 1:
			self.reward = rateDRIV
			return

		if (angle == 90 or mode == 'short') and self.y[-2] == self.N + 1:
			self.reward = rateDRIV
			return
			
		for j in allspot(spot2blk(self.N)):
			if j <= self.N and self.x[j-1] == m_out + 1:
				self.reward = rateDRIV
				break
		return

	def compl_service_helper(self, x, y, a, j):

		for i in range_Hout(j, self.N + 1, self.n):
			assert (i <= self.n)	
			if (i <= spot2blk(j)) and (y[i-1] == self.N + 1) and tau(x, y, i, j):
				return False
			if (i <= spot2blk(j)) and (y[i-1] == self.N + 2) and tau(x, y, i, j):
				if a[i] > 0 and i in range_Hout(j, a[i], self.n):
					return False
				if a[i] == 0 and i + g_in > spot2in(j):
					return False

		for j_ in range(1, self.N + 1):
			if not self.compl_service_helper_step(x, y, a, j, j_):
				return False
				
		for i in range(spot2blk(j) - 1, max(1, spot2in(j) - g_in) - 1, -1 ):
			if y[i-1] == j:
				break
			if y[i-1] == self.N + 2 and a[i] > j:
				import pdb; pdb.set_trace()
				return False 

		if self.y[spot2blk(j) - 1] > 0:
			return False

		for i in range(spot2blk(j) + 1, self.n + 1):
			if self.y[i-1] == 0:
				break
			if self.y[i-1] <= self.N and spot2in(self.y[i-1]) == i and self.x[self.y[i-1] - 1] == 2:
				return False

		if not gamma(x, y, j):
				return False

		return True

	def compl_service_helper_step(self, x, y, a, j, j_):

		if j == j_:
			return True

		if in_Ein(j, j_, x[j_-1] - 1 - m_out):
			return False

		if (x[j_-1] != 2) and (in_Eout(j, j_, x[j_-1] - 2)):
			return False

		if (x[j_-1] == 2) and (j_ > j) and (in_Eout(j, j_, 0)): 
			if (m_out == 1) and (self.y[spot2blk(j_)] == self.N + 1):
				return False
			if (m_out > 1) and (self.x[j_-1] == 3):
				return False

		for i in range_Hout(j, j_, self.n):
			assert (i <= self.n)
			if (y[i-1] == j_ or a[i] == j_) and tau(x, y, i, j):
				return False

		return True