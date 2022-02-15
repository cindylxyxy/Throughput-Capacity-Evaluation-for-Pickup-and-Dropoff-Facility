########################################################################################################################################################
# Last checked: Feb 06, 2022 by Xinyu Liu
# This file defines 2 objects: the system object which stores and updates the system state in a trajectory-based simulation,
# the vehicle object which stores and updates the information related to each vehicle in the system;
# and several helper functions.
# This file assumes for 0-degree confogurations with long spots. The other configurations also uses the system object in this file as the parent class. 

########################################################################################################################################################
import sys
import json
from math import ceil, floor, sqrt
import numpy as np
from heapq import *
from utils import *
from params import *


class SingleLaneSystem():

	def __init__(self, N, seedSERV = None, seedDRIV = None):

		if simType == 'cav':
			self.timeSERV = ParamGen(Expo(rateSERV, seed = seedSERV))
			self.timeDRIV = ParamGen(Cons(rateDRIV, seed = seedDRIV))

		elif simType == 'exp':
			self.timeSERV = ParamGen(Expo(rateSERV, seed = seedSERV))
			self.timeDRIV = ParamGen(Unif2(rateDRIV, seed = seedDRIV))

		elif simType == 'exp1':
			self.timeSERV = ParamGen(Cons(meanSERV, seed = seedSERV))
			self.timeDRIV = ParamGen(Unif2(rateDRIV, seed = seedDRIV))

		elif simType == 'unif':
			self.timeSERV = ParamGen(Unif2(meanSERV, seed = seedSERV))
			self.timeDRIV = ParamGen(Unif2(rateDRIV, seed = seedDRIV))

		elif simType == 'unif1':
			self.timeSERV = ParamGen(Cons(meanSERV, seed = seedSERV))
			self.timeDRIV = ParamGen(Unif2(rateDRIV, seed = seedDRIV))

		elif simType == 'unif2':
			self.timeSERV = ParamGen(Unif2(meanSERV, seed = seedSERV))
			self.timeDRIV = ParamGen(Cons(rateDRIV, seed = seedDRIV))

		else:
			self.timeDRIV = ParamGen(Disc(rateDRIV, seed = seedDRIV))

			if simType == 'triag1':
				self.timeSERV = ParamGen(Tri1(rateSERV, seed = seedSERV))

			if simType == 'triag0':
				self.timeSERV = ParamGen(Tri0(rateSERV, seed = seedSERV))

			if simType == 'triag':
				self.timeSERV = ParamGen(Tria(rateSERV, seed = seedSERV))

		assert side == 'single'
		assert angle == 0
		assert mode == 'singlelane'
		assert LOT_LENGTH == 23.
		assert control == 'partial'
		self.N = N
		self.n = N * LOT_LENGTH

		self.curr = 0.0
		self.start_time = 0.0
		self.eventheap = []
		self.waiting = []
		self.head = None
		self.inCount = 0
		self.outCount = 0	
		self.entry_blocked = self.curr
		self.entry_cleared = self.curr
		self.debug = debug
		self.debug_unit = SIM_UNIT

		self.max_idle = 0.0
		self.max_idle_idx = None
		self.max_pct_idle = 0.0
		self.max_pct_idle_idx = None
		self.prod_time = []
		self.idle_time = []
		self.pct_idle = []

		self.speed = [None for _ in range(self.N)]
		self.service_times = [None for _ in range(self.N)]
		self.service_end_times = [None for _ in range(self.N)]
		self.min_moving_idx = self.N + 1

	def add_event(self, event):
		heappush( self.eventheap, event) 

	def debug_print(self):

		print (self.curr, self.curr_vehicle.stop_idx, self.curr_typ)
		car = self.head
		while car != None:
			car.update_loc()
			print ('Vehicle %s, Spot %s, Status %s, Service end %s, Location %s.' %(car.idx, car.stop_idx, car.status, car.serv_end, car.curr_loc) )
			car = car.nex			

	def idle_time_calc(self, curr_vehicle):

		total = self.curr - curr_vehicle.enter_time 
		if (total < curr_vehicle.prod_time) and np.abs(total - curr_vehicle.prod_time) > 21 * SMALL_INTERVAL:
			import pdb; pdb.set_trace()
		idle = max(0.0, total - curr_vehicle.prod_time)
		if idle > self.max_idle:
			self.max_idle = idle
			self.max_idle_idx = curr_vehicle.idx
		if idle + curr_vehicle.prod_time == 0.0:
			import pdb; pdb.set_trace()
		if idle / (idle + curr_vehicle.prod_time) > self.max_pct_idle:
			self.max_pct_idle = idle / (idle + curr_vehicle.prod_time)
			self.max_pct_idle_idx = curr_vehicle.idx
		self.prod_time.append(curr_vehicle.prod_time)
		self.idle_time.append(idle)
		self.pct_idle.append( idle / (idle + curr_vehicle.prod_time) )

	def run(self):

		while self.curr - self.start_time <= SIM_UNIT: 

			if self.curr - self.start_time >= self.debug_unit:
				import pdb; pdb.set_trace()
			
			curr_event = heappop(self.eventheap)
			curr_vehicle = curr_event.vehicle
			curr_typ  = curr_event.typ
			self.curr_vehicle = curr_vehicle
			self.curr = float(curr_event.time)
			self.curr_typ = curr_typ			
			try:
				stop = curr_vehicle.stop
			except AttributeError:
				assert curr_typ == 'enter_system'

			################################### update system ###################################
			if curr_typ == 'leave_system':
				self.leave_system(curr_vehicle)

			elif curr_typ == 'start_service':
				self.start_service(curr_vehicle)

			elif curr_typ == 'prepare_pulling_out':
				self.prepare_pulling_out(curr_vehicle)

			else:
				assert curr_typ == 'enter_system'
				self.enter_system(curr_vehicle, debug_idx = VEHICLE_IDX)

			if self.debug:
				if (not curr_typ == 'enter_system') or curr_vehicle.idx is not None:
					self.debug_print()
					import pdb; pdb.set_trace()			

		self.start_time = self.curr
		return (self.outCount)

	def leave_system(self, curr_vehicle):

		curr_vehicle.update_loc()

		if np.abs( curr_vehicle.curr_loc - curr_vehicle.dest_to_stop ) > 1e-5:
			if curr_vehicle.end_time <= self.curr + SMALL_INTERVAL:
				import pdb; pdb.set_trace() 
			self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'leave_system') )
			return

		assert self.head == curr_vehicle
		assert curr_vehicle.prev == None					
		self.head = curr_vehicle.nex
		if curr_vehicle.nex != None:
			curr_vehicle.nex.prev = None
			curr_vehicle.nex = None

		assert self.service_times[curr_vehicle.stop - 1] == curr_vehicle.serv_time
		self.service_times[curr_vehicle.stop - 1] = None
		self.service_end_times[curr_vehicle.stop - 1] = None

		self.outCount += 1
		curr_vehicle.prod_time += ( (curr_vehicle.dest_to_stop - curr_vehicle.stop * LOT_LENGTH) / curr_vehicle.driv )
		self.idle_time_calc(curr_vehicle)
		return

	def start_service(self, curr_vehicle):

		curr_vehicle.update_loc()

		if np.abs( curr_vehicle.curr_loc - curr_vehicle.dest_to_stop ) > 1e-5:
			if curr_vehicle.end_time <= self.curr + SMALL_INTERVAL:
				import pdb; pdb.set_trace() 
			self.add_event( event(curr_vehicle.end_time, curr_vehicle, 'start_service') )
			return

		if self.service_times[curr_vehicle.stop - 1] is not None:
			if curr_vehicle.stop != 1:
				import pdb;pdb.set_trace()
			if self.head.stop != 1:
				import pdb; pdb.set_trace()
			if self.head.end_time > self.curr + SMALL_INTERVAL:
				import pdb; pdb.set_trace()
			self.add_event( event(self.head.end_time + SMALL_INTERVAL, curr_vehicle, 'start_service') )
			return	

		curr_vehicle.curr_loc = curr_vehicle.dest_to_stop
		curr_vehicle.start_service()

		if curr_vehicle.nex != None:
			car = curr_vehicle.nex
			while car != None:
				car.update_loc()
				car.update_traj()
				car = car.nex

		assert curr_vehicle.serv_end > self.curr
		self.add_event( event(curr_vehicle.serv_end, curr_vehicle, 'prepare_pulling_out') )
		return

	def prepare_pulling_out(self, curr_vehicle):

		assert curr_vehicle.status == 3
		curr_vehicle.status = 6

		if curr_vehicle.stop == self.N:
			self.min_moving_idx = self.N
		elif curr_vehicle.stop == self.min_moving_idx - 1:
			self.min_moving_idx = curr_vehicle.stop

		curr_vehicle.update_traj()
		
		if curr_vehicle.nex != None:
			car = curr_vehicle.nex
			while car != None and car.status == 6:
				if self.min_moving_idx == car.prev.stop:
					self.min_moving_idx = car.stop
				car.update_loc()
				car.update_traj()
				car = car.nex	

		if self.min_moving_idx == 1:
			self.min_moving_idx = self.N + 1

		assert curr_vehicle.end_time > self.curr
		self.add_event( event(curr_vehicle.end_time, curr_vehicle, 'leave_system') )
		new_vehicle = SLvehicle(self)
		self.add_event( event(self.curr, new_vehicle, 'enter_system'))
		heappush( self.waiting, -curr_vehicle.stop )
		return

	def enter_system(self, curr_vehicle, debug_idx = None):

		assert curr_vehicle.status == 0
		if self.entry_blocked == self.curr:
			self.add_event( event( self.entry_cleared, curr_vehicle, 'enter_system' ) )
			return

		if self.head is None:
			self.head = curr_vehicle
			self.inCount += 1
			heapify(self.waiting)
			curr_vehicle.assign_spot( -heappop(self.waiting) )
			curr_vehicle.update_traj()
			assert curr_vehicle.end_time != None and curr_vehicle.end_time >= self.curr
			self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'start_service') )
			return
		
		car = self.head
		while car.nex is not None:
			car = car.nex	
		car.update_loc()

		req_time = self.curr

		if car.curr_loc == CAR_LENGTH and car.status == 3:
			assert car.stop == 1
			req_time = car.serv_end
		
		elif car.curr_loc <= CAR_LENGTH + SMALL_INTERVAL:
			assert car.status in [1, 6]
			car_time = car.calc_time(CAR_LENGTH + SMALL_INTERVAL)
			if car_time == self.curr:
				cp = car.traj.head
				assert cp.data.t <= self.curr
				while cp.nex is not None and cp.nex.data.t <= self.curr:
					cp = cp.nex
				if cp.nex is None:
					import pdb; pdb.set_trace()
				else:
					assert cp.nex.data.t > self.curr
					if cp.data.v > 0.0:
						pass
					else:
						import pdb; pdb.set_trace()
			else:
				assert car_time > self.curr
				req_time = max( req_time, car_time)

		if req_time > self.curr:
			self.entry_blocked = self.curr
			self.entry_cleared = req_time
			self.add_event( event( self.entry_cleared, curr_vehicle, 'enter_system') )
			return

		self.inCount += 1
		curr_vehicle.prev = car
		car.nex = curr_vehicle
		heapify(self.waiting)
		stop = -heappop(self.waiting)
		if car.stop != 1 and car.stop != stop + 1:
			if - (car.stop - 1) in self.waiting:
				heappush(self.waiting, -stop)
				self.waiting.remove( -(car.stop - 1) )
				heapify(self.waiting)
				stop = car.stop - 1
			else:
				self.debug_print()
				import pdb; pdb.set_trace() 
		curr_vehicle.assign_spot( stop )
		curr_vehicle.update_traj()
		if not (curr_vehicle.end_time != None and curr_vehicle.end_time >= self.curr):
			import pdb; pdb.set_trace()
		self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'start_service') )
		return

class SLvehicle():
	
	def __init__(self, sys, getin = False, stop_idx = None, prev = None):

		self.sys = sys
		self.driv = self.sys.timeDRIV.next()
		self.status = 0

		self.idx = None
		self.stop = None
		self.stop_idx = None
		self.curr_loc = None
		self.dest_to_stop = None
		self.enter_time = None
		self.end_time = None
		self.prod_time = 0.0	

		self.prev = None
		self.nex = None
		self.traj = None

		self.serv_start = None
		self.serv_time = None
		self.serv_end = None

		if getin:
			assert stop_idx is not None
			self.assign_spot(stop_idx)
			self.status = 1
			self.curr_loc = stop_idx * LOT_LENGTH
			if stop_idx == self.sys.N:
				assert prev is None
				self.sys.head = self
			else:
				assert prev is not None
				self.prev = prev
			self.start_service()
			self.prod_time -= (self.dest_to_stop / self.driv)
		else:
			assert stop_idx is None
			self.curr_loc = 0.0

	def assign_spot(self, stop_idx):
		assert self.status == 0
		self.status = 1
		self.idx = self.sys.inCount
		self.stop = stop_idx
		self.stop_idx = stop_idx
		self.enter_time = self.sys.curr
		self.dest_to_stop = self.stop * LOT_LENGTH
		self.sys.speed[self.stop - 1] = self.driv

	def start_service(self):
		assert (self.status == 1)
		self.status = 3
		self.serv_time = self.sys.timeSERV.next()
		self.serv_end = self.serv_time + self.sys.curr
		self.serv_start = self.sys.curr
		if not self.sys.service_times[self.stop - 1] is None:
			import pdb; pdb.set_trace()
		self.sys.service_times[self.stop - 1] = self.serv_time
		self.sys.service_end_times[self.stop - 1] = self.serv_end
		self.prod_time += ( self.serv_time + self.dest_to_stop / self.driv )
		self.dest_to_stop = self.sys.n + CAR_LENGTH
		self.update_traj()

	def calc_time(self, loc):

		if self.curr_loc > loc:
			import pdb; pdb.set_trace()

		cp = self.traj.head

		if self.curr_loc == loc:

			if cp.data.t > self.sys.curr:
				if np.abs(cp.data.t - self.sys.curr) > SMALL_INTERVAL:
					import pdb; pdb.set_trace()
				elif cp.nex is not None:
					import pdb; pdb.set_trace()
				else:
					return cp.data. t

			while cp.nex is not None and cp.nex.data.t <= self.sys.curr:
				cp = cp.nex
			assert cp.data.t <= self.sys.curr

			if cp.nex is None:
				return self.sys.curr
			else:
				assert cp.nex.data.t > self.sys.curr
				if cp.data.v > 0.0:
					return self.sys.curr
				while cp.data.v == 0.0:
					cp = cp.nex
				return cp.data.t

		while cp.nex != None:
			if cp.nex.data.x <= loc:
				cp = cp.nex
			else:
				break

		if cp.data.v == 0.0:
			import pdb;pdb.set_trace()

		if cp.nex == None:
			assert cp.data.v == 'D'
			return cp.data.t

		return cp.data.t + (loc - cp.data.x) / cp.data.v
	
	def update_loc(self):

		if (self.status == 3):
			assert self.curr_loc == self.stop * LOT_LENGTH
			return

		cp = self.traj.head
		if cp.nex is None:
			if not np.abs(cp.data.t - self.sys.curr) <= 1e-05:
				import pdb; pdb.set_trace()
			if not cp.data.x == self.curr_loc:
				if np.abs(cp.data.x - self.curr_loc) <= SMALL_INTERVAL:
					self.curr_loc = cp.data.x
				else:
					import pdb; pdb.set_trace()
			return 

		if cp.data.t > self.sys.curr:
			print ('Line 474: current time not accounted for in the trajectory!!!')
			import pdb; pdb.set_trace()

		if self.end_time < self.sys.curr - 2 * SMALL_INTERVAL:
			print ('Line 478: self.end_time, self.sys.curr', self.end_time, self.sys.curr)
			import pdb; pdb.set_trace()

		while cp.nex.data.t <= self.sys.curr:
			cp = cp.nex
			if cp.nex is None:
				if np.abs(self.end_time - self.sys.curr) > 3 * SMALL_INTERVAL:
					print ('Line 486: it should be small', np.abs(self.end_time - self.sys.curr))
					import pdb; pdb.set_trace()
				self.curr_loc = cp.data.x
				return

		assert cp.data.t <= self.sys.curr < cp.nex.data.t
		assert cp.data.v != 'D'
		self.curr_loc = cp.data.x + cp.data.v * (self.sys.curr - cp.data.t)
		if self.curr_loc > self.dest_to_stop:
			if not (self.curr_loc - self.dest_to_stop < 4.1e-05):
				import pdb; pdb.set_trace()
			self.curr_loc = self.dest_to_stop
		return

	def update_traj(self):

		traj = self.traj
		end_time = self.end_time

		self.traj = DLinkedList()

		if self.curr_loc == self.dest_to_stop and self.status == 6:
			assert self.end_time is not None
			assert (np.abs(self.end_time - self.sys.curr) < SMALL_INTERVAL)
			if (end_time is not None) and (self.end_time < end_time):
				if self.end_time < end_time - 2 * SMALL_INTERVAL:
					import pdb; pdb.set_trace()
				self.end_time = end_time
			assert self.end_time >= self.sys.curr
			self.traj.addEnd( changePoint(self.curr_loc, self.sys.curr, 0.0))
			self.traj.addEnd( changePoint(self.curr_loc, self.end_time, 'D') )
			return

		if self.curr_loc == self.dest_to_stop and self.status == 1:
			assert self.end_time is not None
			if (end_time is not None) and (self.end_time < end_time):
				if self.end_time < end_time - SMALL_INTERVAL:
					import pdb; pdb.set_trace()
				self.end_time = end_time
			if not self.end_time >= self.sys.curr:
				if not np.abs(self.end_time - self.sys.curr) < 2 * SMALL_INTERVAL:
					import pdb; pdb.set_trace()
				self.end_time = self.sys.curr
			self.traj.addEnd( changePoint(self.curr_loc, self.sys.curr, 0.0))
			self.traj.addEnd( changePoint(self.curr_loc, self.end_time, 'D') )
			return

		start_t = self.sys.curr

		if self.status == 3:
			if not (self.curr_loc == self.stop * LOT_LENGTH):
				import pdb; pdb.set_trace()
			assert (self.dest_to_stop == self.sys.n + CAR_LENGTH)

			if self.serv_end > self.sys.curr:
				self.traj.addEnd( changePoint(self.curr_loc, self.sys.curr, 0.0) )
				start_t = self.serv_end
			
		if self.prev is None:
			try:
				assert self.sys.head == self
			except:
				print ('Line 541: this car does not have a prev while self.sys.head is not itself!!!')
				import pdb; pdb.set_trace()

			self.traj.addEnd( changePoint(self.curr_loc, start_t, self.driv) )
			self.end_time = start_t + (self.dest_to_stop - self.curr_loc) / self.driv
			if (end_time is not None) and (self.end_time < end_time):
				if self.end_time < end_time - 2 * SMALL_INTERVAL:
					import pdb; pdb.set_trace()
				self.end_time = end_time
			self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D') )
			return

		if self.prev.end_time < start_t + 10 * SMALL_INTERVAL:
			if self.prev.status == 1 and self.prev.end_time > start_t:
				if self.curr_loc >= self.prev.dest_to_stop - CAR_LENGTH:
					self.curr_loc = self.prev.dest_to_stop - CAR_LENGTH
					self.traj.addEnd(  changePoint(self.curr_loc, start_t, 0.0) )
					start_t = self.prev.end_time 
			self.traj.addEnd( changePoint(self.curr_loc, start_t, self.driv) )
			self.end_time = start_t + (self.dest_to_stop - self.curr_loc) / self.driv
			if (end_time is not None) and (self.end_time < end_time):
				if self.end_time < end_time - 7 * SMALL_INTERVAL:
					import pdb; pdb.set_trace()
				self.end_time = end_time
			self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D') )
			return

		if self.prev.status == 3:
			if self.dest_to_stop <= (self.prev.stop - 1) * LOT_LENGTH:
				self.traj.addEnd( changePoint(self.curr_loc, start_t, self.driv) )
				self.end_time = start_t + (self.dest_to_stop - self.curr_loc) / self.driv
				if (end_time is not None) and (self.end_time < end_time):
					if self.end_time < end_time - 2 * SMALL_INTERVAL:
						import pdb; pdb.set_trace()
					self.end_time = end_time
				self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D') )
				return
		
		#####################################################################################
		cp = self.prev.traj.head
		if cp.nex is None:
			import pdb; pdb.set_trace()
		assert (cp.data.t <= start_t)
		while (cp.nex is not None) and (cp.nex.data.t <= start_t):
			assert cp.data.v != 'D'
			cp = cp.nex

		if (cp.nex is None):
			print ('line 575: testing if an option is possible')
			import pdb; pdb.set_trace()
			assert (self.prev.end_time == cp.nex.data.t > start_t)
			assert cp.data.v == 'D'
			self.traj.addEnd( changePoint(self.curr_loc, start_t, self.driv) )
			self.end_time = start_t + (self.dest_to_stop - self.curr_loc) / self.driv
			if (end_time is not None) and (self.end_time < end_time):
				if self.end_time < end_time - 2 * SMALL_INTERVAL:
					import pdb; pdb.set_trace()
				self.end_time = end_time
			self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D') )
			return
		
		# some sanity check before proceeding
		assert (cp.data.v != 'D') and (cp.nex != None)
		assert (cp.data.t <= start_t < cp.nex.data.t)

		# let us define some counters to keep track of time and location (in the projected traj)
		iter_t = start_t
		if start_t == self.sys.curr:
			self.prev.update_loc()
			prev_loc = self.prev.curr_loc
		else:
			prev_loc = cp.data.x + (iter_t - cp.data.t) * cp.data.v

		# the RHS below is the location of prev at start_t
		# check that the distance between the prev and the curr at start_t is at least CAR_LENGTH
		if 0 <= self.curr_loc + CAR_LENGTH - prev_loc < 12e-5:
			self.curr_loc = max( 0.0, prev_loc - CAR_LENGTH )
			if (prev_loc < CAR_LENGTH):
				self.traj.addEnd( changePoint(0.0, start_t, 0.0) )
				start_t = self.prev.calc_time(CAR_LENGTH + SMALL_INTERVAL)
				assert start_t > self.sys.curr
		start_x = self.curr_loc
		iter_x = self.curr_loc

		# check that the distance between the prev and the curr at start_t is at least CAR_LENGTH
		if self.curr_loc + CAR_LENGTH > prev_loc:
			import pdb; pdb.set_trace()

		if self.curr_loc == self.dest_to_stop and self.status == 1:
			assert self.end_time is not None
			if (end_time is not None) and (self.end_time < end_time):
				if self.end_time < end_time - 2 * SMALL_INTERVAL:
					import pdb; pdb.set_trace()
				self.end_time = end_time
			if not self.end_time >= self.sys.curr:
				if not np.abs(self.end_time - self.sys.curr) < 2 * SMALL_INTERVAL:
					import pdb; pdb.set_trace()
				self.end_time = self.sys.curr
			self.traj.addEnd( changePoint(self.curr_loc, self.sys.curr, 0.0))
			self.traj.addEnd( changePoint(self.curr_loc, self.end_time, 'D') )
			return

		# if the distance == CAR_LENGTH AND cp.data.v == self.driv
		# the curr has to travel at cp.data.v == self.driv for some time
		# i.e. the curr already catches up with the its prev
		#      and the curr has a free-flow rate == cp.data.v s.t the headway distance is maintained 
		if (start_x + CAR_LENGTH == prev_loc) and (cp.data.v == self.driv):
			
			while (cp.nex is not None):
				cp = cp.nex
				if (cp.data.v == 'D') or (cp.data.x >= self.dest_to_stop + CAR_LENGTH):
					self.traj.addEnd( changePoint(start_x, start_t, self.driv) )
					self.end_time = start_t + (self.dest_to_stop - start_x) / self.driv
					if (end_time is not None) and (self.end_time < end_time):
						if self.end_time < end_time - 2 * SMALL_INTERVAL:
							import pdb; pdb.set_trace()
						self.end_time = end_time
					self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D'))
					return
				elif (cp.data.v < self.driv):
					self.traj.addEnd( changePoint(start_x, start_t, self.driv) )
					break
				elif (cp.data.v > self.driv):
					break
				else:
					assert (cp.data.v == self.driv)

			assert (cp is not None) and (cp.data.v != self.driv)
			if not (start_x + (cp.data.t - start_t) * self.driv + CAR_LENGTH == cp.data.x):
				if np.abs(start_x + (cp.data.t - start_t) * self.driv + CAR_LENGTH - cp.data.x) > 2 * SMALL_INTERVAL:
					import pdb; pdb.set_trace()
			iter_t = cp.data.t
			prev_loc = cp.data.x
			iter_x = prev_loc - CAR_LENGTH
		
		while True:

			# if the distance == CAR_LENGTH AND cp.data.v < self.driv
			# the curr has to travel at cp.data.v <= self.driv for some time
			# i.e. the curr already catches up with the its prev
			#      and the curr has a free-flow rate geq cp.data.v s.t the headway distance is maintained 
			# the curr starts to travel at cp.data.v (not constant, but always the same as the prev)

			if (iter_x + CAR_LENGTH == prev_loc) and (cp.data.v < self.driv):

				while (cp.nex is not None):

					assert cp.data.v != 'D' 
					if not (iter_x + CAR_LENGTH == cp.data.x + (iter_t - cp.data.t) * cp.data.v):
						import pdb; pdb.set_trace()

					if cp.data.v <= self.driv:

						if cp.data.v == 0 and np.abs(iter_x - self.dest_to_stop) < SMALL_INTERVAL:
							self.end_time = iter_t
							if (end_time is not None) and (self.end_time < end_time):
								if self.end_time < end_time - SMALL_INTERVAL:
									import pdb; pdb.set_trace()
								self.end_time = end_time						
							self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D') )
							return
						else: 
							self.traj.addEnd( changePoint(iter_x, iter_t, cp.data.v) )

						if cp.nex.data.x - CAR_LENGTH >= self.dest_to_stop:
							if not (cp.data.v > 0.0):
								if ( np.abs(cp.nex.data.x - cp.data.x) < SMALL_INTERVAL ):
									self.end_time = cp.nex.data.t
									if (end_time is not None) and (self.end_time < end_time):
										if self.end_time < end_time - 2 * SMALL_INTERVAL:
											import pdb; pdb.set_trace()
										self.end_time = end_time
									self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D'))
								else:
									import pdb; pdb.set_trace()
							else:
								self.end_time = iter_t + (self.dest_to_stop - iter_x) / cp.data.v
								if (end_time is not None) and (self.end_time < end_time):
									if self.end_time < end_time - 2 * SMALL_INTERVAL:
										import pdb; pdb.set_trace()
									self.end_time = end_time
								self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D'))
							return

						if np.abs(iter_x + (cp.nex.data.t - iter_t) * cp.data.v + CAR_LENGTH - cp.nex.data.x) > 1.2e-04:
							import pdb; pdb.set_trace()
						cp = cp.nex 
						iter_x = cp.data.x - CAR_LENGTH
						iter_t = cp.data.t
					
					else:
						assert cp.data.v > self.driv
						break	

				if cp.nex is None:
					assert cp.data.v == 'D'
					assert (iter_x + CAR_LENGTH == cp.data.x)
					self.traj.addEnd( changePoint(iter_x, iter_t, self.driv) )
					self.end_time = iter_t + (self.dest_to_stop - iter_x) / self.driv
					if (end_time is not None) and (self.end_time < end_time):
						if self.end_time < end_time - 2 * SMALL_INTERVAL:
							import pdb; pdb.set_trace()
						self.end_time = end_time
					self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D'))
					return

				assert cp.nex is not None
				assert (cp.data.v != 'D') and (cp.data.v > self.driv)
				assert (iter_x + CAR_LENGTH == cp.data.x)
				assert (iter_x < self.dest_to_stop)
				start_t = iter_t
				start_x = iter_x
				prev_loc = cp.data.x

			# if the distance is > CAR_LENGTH 
			# OR if the distance == CAR_LENGTH while self.driv < cp.data.v
			# the curr can drive at its own free-flow rate for some time
			if not ( (iter_x + CAR_LENGTH < prev_loc) | (cp.data.v > self.driv) ):
				import pdb; pdb.set_trace()
			assert (iter_x + CAR_LENGTH <= prev_loc)

			# the curr starts to drive at self.driv
			self.traj.addEnd( changePoint(start_x, start_t, self.driv) )

			while (cp.nex is not None):

				if (iter_x >= self.dest_to_stop - SMALL_INTERVAL):
					break

				# as long as self.driv <= cp.data.v,
				# the curr continues to travel at self.driv
				# the inequality is not strict 
				# because the distance between two vehicles (> CAR_LENGTH) will be maintained if both travels at the same rate
				if (self.driv <= cp.data.v):

					if (self.driv == cp.data.v) and (iter_x + CAR_LENGTH == prev_loc):
						cp = cp.nex
						prev_loc = cp.data.x
						iter_x = prev_loc - CAR_LENGTH
						iter_t = cp.data.t
						if (cp.data.v != 'D') and (cp.data.v < self.driv):
							break
					else:
						cp = cp.nex
						prev_loc = cp.data.x
						iter_x = start_x + (cp.data.t - start_t) * self.driv
						iter_t = cp.data.t
						if not (iter_x + CAR_LENGTH <= prev_loc):
							if np.abs(iter_x + CAR_LENGTH - prev_loc) < 1e-05:
								iter_x = prev_loc - CAR_LENGTH
							else:
								import pdb; pdb.set_trace()
					continue

				assert (cp.data.v < self.driv)
				if (iter_x + CAR_LENGTH == prev_loc):
					break

				if np.abs(iter_x + (cp.nex.data.t - iter_t) * self.driv + CAR_LENGTH - cp.nex.data.x) < SMALL_INTERVAL:
					cp = cp.nex
					iter_x = cp.data.x - CAR_LENGTH
					iter_t = cp.data.t
					prev_loc = cp.data.x
					if (cp.data.v == 'D') or (cp.data.v < self.driv):
						break

				elif (iter_x + (cp.nex.data.t - iter_t) * self.driv + CAR_LENGTH < cp.nex.data.x):
					cp = cp.nex
					iter_x = start_x + (cp.data.t - start_t) * self.driv
					iter_t = cp.data.t
					prev_loc = cp.data.x
					assert (iter_x + CAR_LENGTH < cp.data.x)

				else:
					change_t = (prev_loc - CAR_LENGTH - iter_x) / (self.driv - cp.data.v)
					if not (0.0 < change_t < cp.nex.data.t - cp.data.t):
						import pdb; pdb.set_trace()
					if cp.data.v == 0.0:
						iter_x = prev_loc - CAR_LENGTH
						iter_t += change_t
					else:
						iter_t += change_t
						prev_loc = cp.data.x + (iter_t - cp.data.t) * cp.data.v
						iter_x  = prev_loc - CAR_LENGTH
						
					break

			# if the prev ends before the curr catches up with the prev
			# i.e. the curr can travel at self.driv until its own destination (from the info available up to now)
			# OR if the curr catches up with the prev after the destination of the curr
			if (cp.nex is None) | (iter_x >= self.dest_to_stop - SMALL_INTERVAL): 
				self.end_time = start_t + (self.dest_to_stop - start_x) / self.driv
				if (end_time is not None) and (self.end_time < end_time):
					if self.end_time < end_time - 2 * SMALL_INTERVAL:
						import pdb; pdb.set_trace()
					self.end_time = end_time
				self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D'))
				return

			assert (cp.data.v != 'D')
			assert (cp.data.v < self.driv)
			start_x = iter_x
			start_t = iter_t
			if not (prev_loc == start_x + CAR_LENGTH):
				import pdb; pdb.set_trace()

		return
