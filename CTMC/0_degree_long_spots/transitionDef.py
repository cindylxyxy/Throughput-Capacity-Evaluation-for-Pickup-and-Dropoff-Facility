########################################################################################################################################################
# Last checked: Feb 06, 2022 by Xinyu Liu
# This file defines the state objects in the Markov Chain model.
# The state object has two main methods for the mu- and the lambda-transitions, and a few helper functions.
# This file also contains two test functions for the state object: One for random sample path and one for a specifed sample path.
 
########################################################################################################################################################
from copy import deepcopy
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
		self.inCount = self.N
		self.outCount = 0
		self._var = None
		self._suffix = None

		self.idx = 0
		self.mu_next = []
		self.lambda_next = None
		self.reward = 0.
		self.prev = None

	@property
	def var(self):
		self.evaluate()
		return ''.join(map(str, self._var)) + self._suffix
	
	def evaluate(self):
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
			else:
				self._var.append(2)
				self._suffix += str(self.y[i]).zfill(2)

	def __eq__(self, other):
		return self.var == other.var

	def __hash__(self):
		return int(self.var)

	def compl_service_helper(self, x, y, j, j_new):

		if j == j_new:
			return True

		for i in range_Iout(j, self.N + 1, self.n):
			assert (i <= self.n)
			if (control == 'full') and (i > 0) and (i <= spot2blk(j)) and (y[i-1] == self.N + 1):
				return False	
			if (control == 'partial') and (i <= spot2blk(j)) and (y[i-1] == self.N + 1) and tau(x, y, i, j):
				return False

		for j_ in range(1, self.N + 1):
			if not self.compl_service_helper_step(x, y, j, j_new, j_):
				return False

		if j_new is not None:
			if in_Iout(j, j_new, 0, self.n):
				return False
			if (x[j_new-1] == 2) and (in_Eout(j, j_new, 0)):
				return False

		if control == 'partial' and m_in > m_out:
			if not gamma(x, y, j):
				return False

		return True

	def compl_service_helper_step(self, x, y, j, j_new, j_):

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

		if (j != j_new):
			for i in range_Iout(j, j_, self.n):
				assert (i <= self.n)
				if (control == 'full') and (i > 0) and (y[i-1] == j_):
					return False
				if (control == 'partial') and (y[i-1] == j_) and tau(x, y, i, j):
					return False

		return True

	def get_new_helper(self, x, y, j, j_):

		if control == 'partial':
			return self.get_new_partial_helper(x, y, j, j_)
		assert control == 'full'

		if (in_Kout(j, j_, x[j_-1] - 2)) | (in_Kin(j, j_, x[j_-1] - 1 - m_out)):
			return False

		for i in range_Iin(j, j_, self.n):
			if (i <= self.n) and (y[i-1] == j_):
				return False

		if x[j-1] == 2:
			for i in range_Iout(j, j_, self.n):
				assert (i <= self.n)
				if (i != 0) and (y[i-1] == j_):
					return False
			if (j != j_) and in_Ein(j, j_, x[j_-1] - 1 - m_out):
				return False
			if (j != j_) and (x[j_-1] != 2) and (in_Eout(j, j_, x[j_-1] - 2)):
				return False 
			return True

		return True

	def get_new_partial_helper(self, x, y, j, j_):

		if x[j-1] == 2:
			for i in range_Hout(j, j_, self.n):
				if (y[i-1] == j_) and tau(x, y, i, j):
					return False
			if (j != j_) and in_Ein(j, j_, x[j_-1] - 1 - m_out):
				return False
			if (j != j_) and (x[j_-1] != 2) and (in_Eout(j, j_, x[j_-1] - 2)):
				return False 
			return True

		assert (x[j-1] == 0) | (3 <= x[j-1] <= m_out + 1)
		for i in range(1, self.n + 1):
			if y[i-1] == j:
				return False

		return True

	def get_new_largest(self, J_star):

		if J_star != set():
			return sorted(J_star)[-1]
		else:
			return None

	def get_new(self, x, y):
		
		J_star = set()

		for j in range(1, self.N + 1):

			if (x[j-1] == 1) | (x[j-1] >= m_out + 2):
				continue

			new = True			
			if (x[j-1] == 2):
				for i in range_Iout(j, self.N + 1, self.n):
					assert (i <= self.n)
					if (control == 'full') and (i != 0) and (i <= spot2blk(j)) and (y[i-1] == self.N + 1):
						new = False
						break
					if (control == 'partial') and (i <= spot2blk(j)) and (y[i-1] == self.N + 1) and tau(x, y, i, j):
						new = False
						break

			for j_ in range(1, self.N + 1):
				if not self.get_new_helper(x, y, j, j_):
					new = False
					break

			if new and control == 'full' and ((x[j-1] == 0) or (3 <= x[j-1] <= m_out + 1)):
				for i in range(1, self.n + 1):
					if y[i-1] == j:
						new = False
						break

			if new and spmatch and control == 'full':
				if self.x[j-1] == 2:
					new = False
				for k in range(1, self.N + 1):
					if self.x[k - 1] == 3 and in_Iout(k, j, 0, self.n):
						new = False
						break

			if new and control == 'partial' and m_in > m_out and x[j-1] == 2:
				if not gamma(x, y, j):
					new = False

			if new:
				J_star.add(j)
	
		return self.get_new_largest(J_star)

	def mu_transition(self, j):

		assert j <= self.N
		assert self.x[j - 1] == 1
		self.x[j - 1] = 2
		self.waiting.append(j)
		return

	def lambda_transition(self):

		if control == 'partial':
			self.lambda_transition_partial()
			return

		assert control == 'full'

		x = deepcopy(self.x)
		y = deepcopy(self.y)
		
		# Step 1: at most 1 replacement vehicle enters
		if spmatch:
			j_new = None
		else:
			j_new = self.get_new(x, y)
			if j_new is not None and j_new not in self.waiting:
				import pdb; pdb.set_trace()
		
		for j in range(self.N, 0, -1):

			# Step 2: vehicles that finished service but not the first step of exit maneuver
			# 		  will complete the first step; if a replacement vehicle has been assigned,
			#		  the vehicle has precedence to complete the first step of exit maneuver.
			if (x[j - 1] == 2) and (self.compl_service_helper(x, y, j, j_new)):
				if (m_out == 1):
					self.x[j - 1] = 0
					if spot2blk(j) == spot2blk(self.N):
						self.outCount += 1
						assert self.reward > 0.
						self.reward = 0. 
					else:
						assert self.y[spot2blk(j)] == 0
						self.y[spot2blk(j)] = self.N + 1					
				else:
					self.x[j - 1] = 3
			
			# Step 3: other vehicles complete their next step of exit maneuver
			# and Step 4: other vehicles complete their next stop of enter maneuver
			elif (3 <= x[j - 1] <= m_out):
				# i.e. if the vehicle in spot j has completed the first but not all steps of exit maneuver
				self.x[j - 1] = x[j - 1] + 1

			elif (x[j - 1] == m_out + 1):
				# i.e. if the vehicle in spot j has completed (m_out - 1) steps of exit maneuver
				self.x[j - 1] = 0
				if spot2blk(j) == spot2blk(self.N):
					# i.e. vehicle in spot j leaves the facility with its last step of exit maneuver
					self.outCount += 1
					assert self.reward > 0.
					self.reward = 0. 
				else:
					# i.e. vehicle from spot j will be in lane block i(j) + 1 upon completing the last step
					assert self.y[spot2blk(j)] == 0
					self.y[spot2blk(j)] = self.N + 1
		
			elif (2 + m_out <= x[j - 1] <= m_out + m_in - 1):
				# i.e. if the vehicle in spot j has completed the first but not all steps of enter maneuver
				self.x[j - 1] = x[j - 1] + 1

			elif (x[j - 1] == m_out + m_in):
				# i.e. if the vehicle in spot j has completed (m_in - 1) steps of enter maneuver
				self.x[j - 1] = 1

		# Step 5: vehicles move forward when possible
		for i in range(self.n, 0, -1):

			if y[i - 1] == 0:
				# i.e. if there is no vehicle in lane block i
				continue

			if self.y[i - 1] == y[i - 1]:
				self.y[i - 1] = 0
			else:
				assert self.y[i - 1] == self.N + 1

			if (1 <= y[i - 1] <= self.N):
				# i.e. if there is a vehicle in lane block i which has not started an enter maneuver or service but is assigned to spot y[i - 1]

				if (i != spot2in(y[i - 1])):
					self.y[blk2blk(i) - 1] = y[i - 1]

				else:
					assert self.x[y[i - 1] - 1] == 0		
					if (m_in == 1):
						self.x[y[i - 1] - 1] = 1
					else:
						self.x[y[i - 1] - 1] = m_out + 2			
		
			else:
				# i.e. if there is a vehicle in lane block i heading to the facility exit
				assert (y[i - 1] == self.N + 1)
				# note that n = i(N) for 0-degree long spots
				#  and that n = i(N) + 1 for 0-degree short spots and 90-degree
				if (blk2blk(i) != spot2blk(self.N) + 1):
					self.y[blk2blk(i) - 1] = self.N + 1
				else:
					assert (blk2blk(i) == spot2blk(self.N) + 1)
					self.outCount += 1
					assert self.reward > 0.
					self.reward = 0.

		if spmatch:
			j_new = self.get_new(x, y)
			if j_new is not None and j_new not in self.waiting:
				import pdb; pdb.set_trace()

		if j_new is not None:
			assert self.y[blk2blk(0) - 1] == 0
			if spot2in(j_new) == 0:
				assert self.x[j_new - 1] == 0
				if (m_in == 1):
					self.x[j_new - 1] = 1
				else:
					self.x[j_new - 1] = m_out + 2
			else:		
				self.y[blk2blk(0) - 1] = j_new
			self.inCount += 1
			self.waiting.remove(j_new)

		# Last step: update the reward for the new state
		assert self.reward == 0.
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

	def lambda_transition_partial(self):

		x = deepcopy(self.x)
		y = deepcopy(self.y)
		blockedlane = set()
		visitedlane = set()
		
		# Step 1: at most 1 replacement vehicle enters
		j_new = self.get_new(x, y)
		
		if j_new is not None and j_new not in self.waiting:
			import pdb; pdb.set_trace()
		
		if spmatch:
			j_new = None

		for j in range(self.N, 0, -1):

			# Step 2: vehicles that finished service but not the first step of exit maneuver
			# 		  will complete the first step; if a replacement vehicle has been assigned,
			#		  the vehicle has precedence to complete the first step of exit maneuver.
			if (x[j - 1] == 2) and (self.compl_service_helper(x, y, j, j_new)):
				assert (m_out > 1)
				self.x[j - 1] = 3
				blockedlane.add(spot2blk(j))
				blockedlane.add(spot2blk(j) + 1)
			
			# Step 3: other vehicles complete their next step of exit maneuver
			# and Step 4: other vehicles complete their next stop of enter maneuver
			elif (3 <= x[j - 1] <= m_out):
				# i.e. if the vehicle in spot j has completed the first but not all steps of exit maneuver
				self.x[j - 1] = x[j - 1] + 1
				blockedlane.add(spot2blk(j))
				blockedlane.add(spot2blk(j) + 1)

			elif (x[j - 1] == m_out + 1):
				# i.e. if the vehicle in spot j has completed (m_out - 1) steps of exit maneuver
				self.x[j - 1] = 0
				if spot2blk(j) == spot2blk(self.N):
					# i.e. vehicle in spot j leaves the facility with its last step of exit maneuver
					self.outCount += 1
					assert self.reward > 0.
					self.reward = 0. 
				else:
					# i.e. vehicle from spot j will be in lane block i(j) + 1 upon completing the last step
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
				blockedlane.add(spot2blk(j))

		# Step 5: vehicles move forward when possible
		for i in range(self.n, 0, -1):

			if y[i - 1] == 0:
				# i.e. if there is no vehicle in lane block i
				continue

			if self.y[i - 1] == y[i - 1] and i not in visitedlane:
				self.y[i - 1] = 0
			else:
				assert self.y[i - 1] == self.N + 1

			if (1 <= y[i - 1] <= self.N):
				# i.e. if there is a vehicle in lane block i which has not started an enter maneuver or service but is assigned to spot y[i - 1]

				if angle == 0 and mode == 'long':
					# i.e. if a vehicle enters or exits the boarding spot both by going forward 

					if (spot2in(y[i - 1]) == i) and (self.x[y[i - 1] - 1] == 0):
					# if spot2in(y[i - 1]) == i and (x[y[i - 1] - 1] == 0):
						# i.e. if a vehicle has arrived at the lane block i_in(j) = i(j) - 1, and is ready to start an enter maneuver,
						#      and there is no vehicle exiting from spot j 
						if (blk2blk(i) in blockedlane) and (1 <= y[blk2blk(i) - 1] <= self.N) and (self.x[y[blk2blk(i) - 1] - 1] != m_out + 2):
							assert self.y[i - 1] == 0
							self.y[i - 1] = y[i - 1]
							blockedlane.add(i)
							continue
						assert m_in > 1
						self.x[y[i - 1] - 1] = m_out + 2
						blockedlane.add(i)
						blockedlane.add(spot2blk(y[i - 1]))

					elif spot2in(y[i - 1]) == i:
						# i.e. continuing from the previous case, if another vehicle is still exiting from spot j
						assert 3 <= x[y[i - 1] - 1] <= m_out + 1
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

					else:
						# i.e. the vehicle has not yet arrived at a lane block beside its assigned boarding spot
						assert spot2blk(y[i - 1]) > i
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
					self.outCount += 1
					assert self.reward > 0.
					self.reward = 0. 

		if spmatch:
			j_new = self.get_new(x, y)
			if j_new is not None and j_new not in self.waiting:
				import pdb; pdb.set_trace()

		if j_new is not None:
			new = (blk2blk(0) not in blockedlane)
			if (angle == 0 and mode == 'long') and (idx2spot(j_new) == 1):
				if (1 <= y[blk2blk(0) - 1] <= self.N and self.x[y[blk2blk(0) - 1] - 1] == m_out + 2):
					new = True
				elif side == 'double' and y[blk2blk(0) - 1] == 0: 
					if (1 in self.waiting and self.x[1] == 3 and self.x[0] != 2):
						assert self.x[0] == 0
						new = True
						j_new = 1
					elif (2 in self.waiting and self.x[0] == 3 and self.x[1] != 2):
						assert self.x[1] == 0
						new = True
						j_new = 2

			if new:
				assert self.y[blk2blk(0) - 1] == 0
				if spot2in(j_new) == 0:
					assert self.x[j_new - 1] == 0
					assert m_in > 1
					self.x[j_new - 1] = m_out + 2
				else:
					self.y[blk2blk(0) - 1] = j_new
				self.waiting.remove(j_new)
				self.enter = True
				self.inCount += 1
		
		# Last step: update the reward for the new state
		assert self.reward == 0.
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


def TestRandom():
	N = random.randint(8, 28)
	print ('\nTest Case (%s spaces)' %N)
		
	a = stateObject(N)
	ITE = 100

	print ('\n1')
	print (a.x)
	print (a.y)

	v = random.sample(list(range(1, N+1)), 1)[0]
	a.mu_transition(v)

	print ('\n2')
	print ('vehicle at %s finishing service ...' %v)
	print (a.x)
	print (a.y)

	for ite in range(3, ITE+1):
		print ('\n%s' %ite)

		keys = [] 
		for j in range(1, a.N + 1):
			if a.x[j-1] == 1:
				keys.append(j)

		V = len(keys)
		if ((np.sum(a.y) > 0) | (a.x != np.ones(a.N)).any()) and (random.random() < rlambda/(V*rmu+rlambda)):
			if control == 'full':
				a.lambda_transition()
			else:
				assert control == 'partial'
				a.lambda_transition_partial()		
			print ('driving ...')
			print (a.x)
			print (a.y)
		elif V > 0:
			v = random.sample(keys, 1)[0]
			print ('vehicle at %s finishing service ...' %v)
			a.mu_transition(v)
			print (a.x)
			print (a.y)
		else:
			print ('The system is empty now!')
			import pdb; pdb.set_trace()

		import pdb; pdb.set_trace()


def TestSamplePath():
	N = 9
	print ('\nTest a fixed sample path (%s spaces)' %N)
		
	a = stateObject(N)
	print ('\n1')
	print (a.x)
	print (a.y)

	a.mu_transition(2)
	print ('\n2')
	print (a.x)
	print (a.y)

	a.lambda_transition()
	print ('\n3')
	print (a.x)
	print (a.y)

	a.lambda_transition()
	print ('\n4')
	print (a.x)
	print (a.y)

	a.mu_transition(8)
	print ('\n5')
	print (a.x)
	print (a.y)

	a.lambda_transition()
	print ('\n6')
	print (a.x)
	print (a.y)

	a.lambda_transition()
	print ('\n7')
	print (a.x)
	print (a.y)

	a.lambda_transition()
	print ('\n8')
	print (a.x)
	print (a.y)

	a.lambda_transition()
	print ('\n9')
	print (a.x)
	print (a.y)

	a.lambda_transition()
	print ('\n10')
	print (a.x)
	print (a.y)

	a.lambda_transition()
	print ('\n11')
	print (a.x)
	print (a.y)

	a.lambda_transition()
	print ('\n12')
	print (a.x)
	print (a.y)

	a.lambda_transition()
	print ('\n13')
	print (a.x)
	print (a.y)

	a.lambda_transition()
	print ('\n14')
	print (a.x)
	print (a.y)

	a.lambda_transition()
	print ('\n15')
	print (a.x)
	print (a.y)

	a.lambda_transition()
	print ('\n16')
	print (a.x)
	print (a.y)

###################################################
if __name__ == "__main__":
	TestRandom()
