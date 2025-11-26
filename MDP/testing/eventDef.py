########################################################################################################################################################
# Last updated: Oct 10, 2025 by Xinyu Liu
# This file defines the objects which stores and updates the system state in a Markov Decision Process simulation.
# The system object is a child object of stateObject as defined in transitionDef.py
# which inherits the two main methods describing the mu- and the lambda-transitions.
# The system object has a slight modification on both lambda transitions (for full and partial access control), i.e. to generate random service times.
# The system object also has an additional "run" method to facilitate the simulation runs.

########################################################################################################################################################
from math import ceil, floor, sqrt, exp
from copy import deepcopy
import sys
import numpy as np
import pandas as pd
import pickle
from heapq import *
from collections import *
from utils import *
from params import *
from inputSetDef import *
from transitionDef import stateObject


class simulationObject(stateObject):

	def __init__(self, N, seedDRIV = None, seedSERV = None):

		super(simulationObject, self).__init__(N)

		if simType == 'mc':
			self.timeDRIV = ParamGen(Expo(rateDRIV, seed = seedDRIV))
			self.timeSERV = ParamGen(Expo(rateSERV, seed = seedSERV))
		if simType == 'cav':
			self.timeDRIV = ParamGen(Cons(meanDRIV, seed = seedDRIV))
			self.timeSERV = ParamGen(Expo(rateSERV, seed = seedSERV))

		self.curr = 0.0
		self.start_time = 0.0
		self.cycle_times = []
		self.eventheap = []

		self.service_times = dict()
		for j in range(1, self.N + 1):
			self.service_times[j] = []
		self.out_times = []
		self.in_times  = []
		self.inCount = 0
		self.outCount = 0
		self.simReward = 0.0

	def add_event(self, event):
		heappush( self.eventheap, event) 

	# Algorithm 6 Computation of (x',y') = f_{\lambda}((x,y),a) in the MDP
	def lambda_transition(self, a):

		self.get_actions()

		# Step 0: Initialize
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
					self.outCount += 1
					self.out_times.append(self.curr)
					self.simReward += (exp(- DISC_RATE * (self.curr - self.start_time)))

			elif 1 <= x[j-1] - (m_out(j)+1) <= m_in(j)-2:
				# i.e. if the vehicle in spot j has completed the first but not all steps of enter maneuver
				assert self.x[j-1] == x[j-1]
				self.x[j-1] = x[j-1] + 1

			elif x[j-1] - (m_out(j)+1) == m_in(j)-1:
				# i.e. if the vehicle in spot j has completed all but the last step of enter maneuver
				assert self.x[j-1] == x[j-1]
				self.x[j-1] = 1
				self.prob_mu[j-1] = rateSERV
				
				next_serv = self.timeSERV.next()
				self.add_event( event(self.curr + next_serv, j, 'mu') )
				self.service_times[j].append(next_serv)

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
				assert self.y[i-1] == 0
				self.y[i-1] = z[i-1]

			elif (1 <= z[i-1] <= self.N and i < spot2in(z[i-1])) or (z[i-1] == self.N+1 and blk2blk(i) <= self.n) or (z[i-1] == self.N+2):
				assert self.y[blk2blk(i)-1] == 0
				self.y[blk2blk(i)-1] = z[i-1]
		
			elif (z[i-1] == self.N+1 and blk2blk(i) > self.n):
				self.outCount += 1
				self.out_times.append(self.curr)
				self.simReward += (exp(- DISC_RATE * (self.curr - self.start_time)))

		# Step 3: A replacement vehicle may enter the facility
		if 1 <= a0 <= self.N:	
			if spot2in(a0) == 0:
				assert self.x[a0-1] == 0

				if m_in(a0) == 1:
					self.x[a0-1] = 1
					self.prob_mu[a0-1] = rateSERV

					next_serv = self.timeSERV.next()
					self.add_event( event(self.curr + next_serv, a0, 'mu') )
					self.service_times[a0].append(next_serv)
				else:
					self.x[a0-1] = m_out(a0)+2
			else:
				assert self.y[blk2blk(0)-1] == 0
				self.y[blk2blk(0)-1] = a0 

			self.inCount += 1
			self.in_times.append(self.curr)

		elif a0 == self.N+1:
			assert self.y[blk2blk(0)-1] == 0
			self.y[blk2blk(0)-1] = self.N+2

			self.inCount += 1
			self.in_times.append(self.curr)

		# Step 4: Vehicles that completed service but have not yet started to exit will, 
		#         if possible, complete the first step of their exit maneuvers
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
						self.outCount += 1
						self.out_times.append(self.curr)
						self.simReward += (exp(- DISC_RATE * (self.curr - self.start_time)))
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
							self.outCount += 1
							self.out_times.append(self.curr)
							self.simReward += (exp(- DISC_RATE * (self.curr - self.start_time)))
					else:
						self.x[j-1] = 3

		if np.abs( np.sum( self.prob_mu ) - np.sum([mu(j) for j in self.x if j == 1]) ) > 1e-12:
			import pdb; pdb.set_trace()
		self.prob_lambda = int( ( sum(self.y) > 0) | (self.x.count(1) < self.N) ) * rateDRIV
		
		return

	def run(self, state_dict, opt_actions):

		if self.curr == 0.0:

			assert len(self.eventheap) == 0
			assert self.inCount == 0
			assert self.outCount == 0
			print ('starting properly...')
			
			# initialize all possible mu transitions common to all control policies
			for j in range(1, self.N + 1):
				next_serv = self.timeSERV.next()
				self.add_event( event(self.curr + next_serv, j, 'mu') )
				self.service_times[j].append(next_serv)

			# run the first mu transition
			first_event = heappop(self.eventheap)
			self.curr = float(first_event.time)

			assert first_event.typ == 'mu'
			assert first_event.vehicle is not None
			self.mu_transition(first_event.vehicle)

			# add the first lambda transition
			assert (np.sum(self.y) > 0) | ( (self.x != np.ones(self.N)).any() )
			self.add_event( event(self.curr + self.timeDRIV.next(), None, 'lambd') )

		while self.curr - self.start_time <= SIM_UNIT:

			curr_event = heappop(self.eventheap)
			self.curr = float(curr_event.time)
			
			# update system 
			if curr_event.typ == 'mu':
				assert curr_event.vehicle is not None 
				self.mu_transition(curr_event.vehicle)				
			else:
				assert curr_event.typ == 'lambd'
				assert self in state_dict
				a = opt_actions[ state_dict[self] ]
				if a is not None:
					self.lambda_transition(a)
				self.add_event( event(self.curr + self.timeDRIV.next(), None, 'lambd') )

		self.start_time = self.curr
		return (self.outCount)