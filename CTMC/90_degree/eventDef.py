########################################################################################################################################################
# Last checked: Feb 06, 2022 by Xinyu Liu
# This file defines the objects which stores and updates the system state in a Markov Chain simulation.
# The system object is a child object of stateObject as defined in transitionDef.py
# which inherits the two main methods describing the mu- and the lambda-transitions.
# The system object has a slight modification on both lambda transitions (for full and partial access control), i.e. to generate random service times.
# The system object also has an additional "run" method to facilitate the simulation runs.

########################################################################################################################################################
from math import ceil, floor, sqrt, exp
from copy import deepcopy
import sys
import numpy as np
import pickle
from heapq import *
from collections import *
from utils import *
from params import *
from inputSetDef import *
from transitionDef import stateObject


class system(stateObject):

	def __init__(self, N, seedDRIV = None, seedSERV = None):

		super(system, self).__init__(N)
		self.randomEntry = None
		# self.randomEntry = ParamGen(Unif(2., seed = seedENTR))
		self.probEntry = None
		self.maxEntry = None
		self.inSystem = self.inCount
		self.onlane = 0
		self.totalTime = 0.
		self.last_lambda = 0.
		self.last_mu = [None for _ in self.x]

		# self.notRightSet = set()
		# self.ref_state = pickle.load(open('%s/%s_%s.p'%(MCdirname, MCfilename, N), 'rb'))
		# print ('%s/%s_%s.p is found and read!'%(MCdirname, MCfilename, N))
		# try:
		# 	self.ref_state = pickle.load(open('%s/%s_%s.p'%(MCdirname, MCfilename, N), 'rb'))
		# except:
		# 	print ('%s/%s_%s.p does not exist!'%(MCdirname, MCfilename, N))

		if simType == 'mc':
			self.timeDRIV = ParamGen(Expo(rateDRIV, seed = seedDRIV))
			self.timeSERV = ParamGen(Expo(rateSERV, seed = seedSERV))
		if simType == 'cav':
			self.timeDRIV = ParamGen(Cons(meanDRIV, seed = seedDRIV))
			self.timeSERV = ParamGen(Expo(rateSERV, seed = seedSERV))
		if simType == 'det':
			self.timeDRIV = ParamGen(Cons(meanDRIV, seed = seedDRIV))
			self.timeSERV = ParamGen(Cons(meanSERV, seed = seedSERV))			

		self.curr = 0.0
		self.start_time = 0.0
		self.eventheap = []
		self.debug = debug
		self.debug_unit = SIM_UNIT
		self.enter = False

		self.service_times = dict()
		for j in range(1, self.N + 1):
			self.service_times[j] = []
		self.out_times = []
		self.outfile = None
		self.writer = None

	def state_to_y(self):

		new_y = deepcopy(self.y)

		for i in range(1, len(self.y) + 1):

			if self.y[i-1] == 0:
				continue
			new_y[i-1] = 1

		for j in range(1, len(self.x) + 1):
			if self.x[j-1] in [0, 1, 2]:
				continue
			if self.x[j-1] <= m_out + 1:
				assert self.y[spot2blk(j)-1] == 0
				new_y[spot2blk(j)-1] = 3
				curr = spot2blk(j) - 1
			elif self.x[j-1] <= m_in + m_out:
				assert self.y[spot2blk(j)-1] == 0
				assert self.y[spot2in(j)-1] == 0
				new_y[min(spot2blk(j), spot2in(j))-1] = 4
				curr = min(spot2blk(j), spot2in(j)) - 1
			else:
				print (self.x[j-1])

			while curr >= 1 and new_y[curr - 1] == 1:
				new_y[curr - 1] = 2
				curr -= 1
		return new_y

	def add_event(self, event):
		heappush( self.eventheap, event) 

	def get_new_earliest(self, J_star):

		if J_star != set():
			for j in self.waiting:
				if j in J_star:
					return j
			import pdb; pdb.set_trace()
		else:
			return None

	
	'''
	# only uncomment the method below if you want to change the assignment policy of spots for the replacement vehicles
	# the default choice is the spot with the largest index
	# uncomment this section to change to the spot whose last vehicle finished the earliest 
	
	def get_new(self, x, y):

		J_star = set()
		for j in range(1, self.N + 1):
			if (x[j-1] == 1) | (x[j-1] >= m_out + 2):
				continue

			new = True			
			if (x[j-1] == 2):
				for i in range_Iout(j, self.N + 1, self.n):
					assert (i <= self.n)
					if (i != 0) and (i <= spot2blk(j)) and (y[i-1] == self.N + 1):
						new = False
						break

			for j_ in range(1, self.N + 1):
				if not self.get_new_helper(x, y, j, j_):
					new = False
					break
			
			if new:
				J_star.add(j)

		return self.get_new_earliest(J_star)

	'''

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
				self.onlane += 1
				if self.totalTime is not None:
					assert self.last_mu[j - 1] is not None
					self.totalTime += (self.curr - self.last_mu[j - 1])
					self.last_mu[j - 1] = None
				if (m_out == 1):
					self.x[j - 1] = 0
					if spot2blk(j) == spot2blk(self.N):
						self.outCount += 1
						assert self.inSystem >= 1
						self.inSystem -= 1
						if not self.onlane >= 1:
							import pdb; pdb.set_trace()
						self.onlane -= 1
						self.out_times.append(self.curr)
						self.reward += (exp(- DISC_RATE * self.curr))
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
					assert self.inSystem >= 1
					self.inSystem -= 1
					if not self.onlane >= 1:
						import pdb; pdb.set_trace()
					self.onlane -= 1
					self.out_times.append(self.curr)
					self.reward += (exp(- DISC_RATE * self.curr))
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
				next_serv = self.timeSERV.next()
				self.add_event( event(self.curr + next_serv, j, 'mu') )
				self.service_times[j].append(next_serv)
				if not self.onlane >= 1:
					import pdb; pdb.set_trace()
				self.onlane -= 1 
		
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
						next_serv = self.timeSERV.next()
						self.add_event( event(self.curr + next_serv, y[i - 1], 'mu') )
						self.service_times[y[i - 1]].append(next_serv)
						if not self.onlane >= 1:
							import pdb; pdb.set_trace()
						self.onlane -= 1
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
					assert self.inSystem >= 1
					self.inSystem -= 1
					if not self.onlane >= 1:
						import pdb; pdb.set_trace()
					self.onlane -= 1 
					self.out_times.append(self.curr)
					self.reward += (exp(- DISC_RATE * self.curr))

		if spmatch:
			j_new = self.get_new(x, y)
			if j_new is not None and j_new not in self.waiting:
				import pdb; pdb.set_trace()

		if j_new is not None and (self.maxEntry is None or self.inSystem < self.maxEntry):
			assert self.y[blk2blk(0) - 1] == 0
			self.onlane += 1
			if self.totalTime is not None:
				self.totalTime += (self.curr - self.last_lambda)
			if spot2in(j_new) == 0:
				assert self.x[j_new - 1] == 0
				if (m_in == 1):
					self.x[j_new - 1] = 1
					next_serv = self.timeSERV.next()
					self.add_event( event(self.curr + next_serv, j_new, 'mu') )
					self.service_times[j_new].append(next_serv)
					if not self.onlane >= 1:
						import pdb; pdb.set_trace()
					self.onlane -= 1
				else:
					self.x[j_new - 1] = m_out + 2
			else:		
				self.y[blk2blk(0) - 1] = j_new
			self.inCount += 1
			self.inSystem += 1
			self.waiting.remove(j_new)
			self.enter = True

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
				self.onlane += 1
				if self.totalTime is not None:
					self.last_mu[j - 1] is not None
					self.totalTime += (self.curr - self.last_mu[j - 1])
					self.last_mu[j - 1] = None

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
					assert self.inSystem >= 1
					self.inSystem -= 1
					if not self.onlane >= 1:
						import pdb; pdb.set_trace()
					self.onlane -= 1
					self.out_times.append(self.curr)
					self.reward += (exp(- DISC_RATE * self.curr))
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
				next_serv = self.timeSERV.next()
				self.add_event( event(self.curr + next_serv, j, 'mu') )
				self.service_times[j].append(next_serv)
				if not self.onlane >= 1:
					import pdb; pdb.set_trace()
				self.onlane -= 1

		# Step 5: vehicles move forward when possible
		for i in range(self.n, 0, -1):

			if y[i - 1] == 0:
				# i.e. if there is no vehicle in lane block i
				continue

			if self.y[i - 1] == y[i - 1] and i not in visitedlane:
				self.y[i - 1] = 0
			else:
				assert self.y[i - 1] == self.N + 1

			if (1 <= y[i-1] <= self.N):
				# i.e. if there is a vehicle in lane block i which has not started an enter maneuver or service but is assigned to spot y[i - 1]
	
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
						assert m_in > 1
						self.x[y[i - 1] - 1] = m_out + 2
						blockedlane.add(i)
						blockedlane.add(spot2blk(y[i - 1]))

					elif spot2in(y[i - 1]) == i:
						# i.e. continuing from the previous case, if another vehicle is still exiting from spot j
						assert 3 <= self.x[y[i - 1] - 1] <= m_out + 1
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
					assert self.inSystem >= 1
					self.inSystem -= 1
					if not self.onlane >= 1:
						import pdb; pdb.set_trace()
					self.onlane -= 1
					self.out_times.append(self.curr)
					self.reward += (exp(- DISC_RATE * self.curr))

		if spmatch:
			j_new = self.get_new(x, y)
			if j_new is not None and j_new not in self.waiting:
				import pdb; pdb.set_trace()

		if j_new is not None and (self.maxEntry is None or self.inSystem < self.maxEntry):
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
				self.inSystem += 1
				self.onlane += 1
				if self.totalTime is not None:
					self.totalTime += (self.curr - self.last_lambda)

		return

	def run(self, warm_up = False):

		STOP_TIME = SIM_UNIT
		if warm_up:
			STOP_TIME = WARMUP_UNIT

		if self.curr == 0.0:

			print ('starting properly...')

			first_event = heappop(self.eventheap)
			self.curr = float(first_event.time)
			curr_time = first_event.time
			curr_vehicle = first_event.vehicle
			curr_typ = first_event.typ
			
			if self.debug:
				print (curr_time, curr_vehicle, curr_typ)

			if curr_typ == 'mu':
				self.mu_transition(curr_vehicle)
				if self.totalTime is not None:
					assert self.last_mu[curr_vehicle - 1] is None
					self.last_mu[curr_vehicle - 1] = self.curr
					assert (np.sum(self.y) > 0) | ( (self.x != np.ones(self.N)).any() )
			else:
				assert curr_typ == 'lambd'
				if self.totalTime is not None:
					self.totalTime += (self.curr - self.last_lambda) * self.onlane
				self.lambda_transition()
				self.last_lambda = self.curr
			
			self.add_event( event(self.curr + self.timeDRIV.next(), None, 'lambd') )

			if self.debug:
				print (self.x)
				print (self.y)
				import pdb;pdb.set_trace()

		while self.curr - self.start_time <= STOP_TIME:

			if self.curr - self.start_time > self.debug_unit:
				import pdb; pdb.set_trace()
				
			curr_event = heappop(self.eventheap)
			curr_time = curr_event.time
			self.curr = float(curr_time)
			curr_vehicle = curr_event.vehicle
			curr_typ = curr_event.typ

			################################### update system ###################################
			if curr_typ == 'mu':
				self.mu_transition(curr_vehicle)
				if self.totalTime is not None:
					assert self.last_mu[curr_vehicle - 1] is None
					self.last_mu[curr_vehicle - 1] = self.curr

			else:
				assert curr_typ == 'lambd'
				if self.totalTime is not None:
					self.totalTime += (self.curr - self.last_lambda) * self.onlane				
				self.lambda_transition()
				self.last_lambda = self.curr
				self.add_event( event(self.curr + self.timeDRIV.next(), None, 'lambd') )
				if self.outfile is not None:
					if side == 'single':
						if len(self.x) == len(self.y): 	
							self.writer.writerow( ['lambd', 'x'] + self.x )
						else:
							assert len(self.x) + 1 == len(self.y)
							self.writer.writerow( ['lambd', 'x'] + self.x + [0])
						self.writer.writerow( ['lambd', 'y'] + self.y )
					else:
						assert side == 'double'
						if np.abs( 0.5 * len(self.x) - len(self.y) ) <= 1e-10: 	
							self.writer.writerow( ['lambd', 'x'] + list(np.array(self.x)[list(range(0, len(self.x), 2))]) )
						else:
							assert np.abs( 0.5 * len(self.x) + 1 - len(self.y) ) <= 1e-10 
							self.writer.writerow( ['lambd', 'x'] + list(np.array(self.x)[list(range(0, len(self.x), 2))] + [0]) )
						self.writer.writerow( ['lambd', 'y'] + self.y )

						if np.abs( 0.5 * len(self.x) - len(self.y) ) <= 1e-10: 	
							self.writer.writerow( ['lambd', 'x'] + list(np.array(self.x)[list(range(1, len(self.x), 2))]) )
						else:
							assert np.abs( 0.5 * len(self.x) + 1 - len(self.y) ) <= 1e-10 
							self.writer.writerow( ['lambd', 'x'] + list(np.array(self.x)[list(range(1, len(self.x), 2))] + [0]) )

			# try:
			# 	if self not in self.ref_state:
			# 		self.notRightSet.add(self.var)
			# 		import pdb; pdb.set_trace()
			# except:
			# 	import pdb; pdb.set_trace()
			# 	pass

			if self.debug:
				if spmatch:
					print (curr_time, curr_vehicle, curr_typ)				
					print (self.x)
					print (self.y)
					if self.N <= 8:
						service_dict = {}
						for an_event in self.eventheap:
							if an_event.typ == 'mu':
								service_dict[an_event.vehicle] = an_event.time
						for idx in range(1, self.N + 1):
							if idx in service_dict:
								print( idx, service_dict[idx] )
					if curr_typ == 'lambd' and self.enter:
						import pdb; pdb.set_trace()
			self.enter = False

		self.start_time = self.curr
		return (self.outCount)
