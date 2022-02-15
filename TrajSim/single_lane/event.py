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


def idx2spot(idx):
	if side == 'double':
		return ceil(idx / 2)
	else:
		assert side == 'single'
		return idx

def spot2blk(stop):
	# returns the lane block next to boarding stop 
	# for 1 <= stop <= system.half_N when double_sided
	if nUnit == 1:
		return stop
	else:
		return ceil(stop / nUnit)

def in_range(stop_idx, N):

	if (side == 'double'):
		
		if (angle == 0) and (mode == 'long'):
			if (stop_idx % 2 == 1):
				return range( min(N, stop_idx + 2), min(N, stop_idx + 3) + 1 )
			return range( min(N, stop_idx + 1), min(N, stop_idx + 2) + 1 )
		
		if (angle == 0) and (mode == 'short'):
			if (stop_idx % 2 == 1):
				return range( max(1, stop_idx - 2), stop_idx + 2 )
			return range( max(1, stop_idx - 3), stop_idx )

		if (angle == 90):
			if (stop_idx % 2 == 1):
				return range( max(1, stop_idx - 4), min(N, stop_idx + 5) + 1 )
			return range( max(1, stop_idx - 5), min(N, stop_idx + 4) + 1 )

	assert (side == 'single')
	if (angle == 0) and (mode == 'long'):
		return range( min(N, stop_idx + 1), min(N, stop_idx + 1) + 1)
	if (angle == 0) and (mode == 'short'):
		return range( max(1, stop_idx - 1), stop_idx )
	assert (angle == 90)
	return range( max(1, stop_idx - 2), min(N, stop_idx + 2) + 1 )

def out_range(stop_idx, N):
	
	if (side == 'double'):

		if (angle == 0) and (mode == 'long'):
			if (stop_idx % 2 == 1):
				return range( max(1, stop_idx - 2), stop_idx + 2 )
			return range( max(1, stop_idx - 3), stop_idx )

		if (angle == 0) and (mode == 'short'):
			if (stop_idx % 2 == 1):
				return range( max(1, stop_idx - 2), stop_idx + 2 )
			return range( max(1, stop_idx - 3), stop_idx )

		if (angle == 90):
			if (stop_idx % 2 == 1):
				return range( max(1, stop_idx - 4), min(N, stop_idx + 5) + 1 )
			return range( max(1, stop_idx - 5), min(N, stop_idx + 4) + 1 )

	assert side == 'single'

	if (angle == 0):
		return range( max(1, stop_idx - 1), stop_idx )

	assert (angle == 90)
	return range( max(1, stop_idx - 2), min(N, stop_idx + 2) + 1 )

class system():

	def __init__(self, N, seedSERV = None, seedPOUT = None, seedPLIN = None, seedDRIV = None):

		if simType == 'det':
			self.timeSERV = ParamGen(Cons(meanSERV, seed = seedSERV))
			self.timeDRIV = ParamGen(Cons(rateDRIV, seed = seedDRIV))
			self.timePOUT = ParamGen(Cons(meanPOUT, seed = seedPOUT))
			self.timePLIN = ParamGen(Cons(meanPLIN, seed = seedPLIN))

		elif simType == 'cav':
			self.timeSERV = ParamGen(Expo(rateSERV, seed = seedSERV))
			self.timeDRIV = ParamGen(Cons(rateDRIV, seed = seedDRIV))
			self.timePOUT = ParamGen(Cons(meanPOUT, seed = seedPOUT))
			self.timePLIN = ParamGen(Cons(meanPLIN, seed = seedPLIN))

		elif simType == 'exp':
			self.timeSERV = ParamGen(Expo(rateSERV, seed = seedSERV))
			self.timeDRIV = ParamGen(Unif2(rateDRIV, seed = seedDRIV))
			self.timePOUT = ParamGen(Erla(meanDRIV, m_out, seed = seedPOUT))
			self.timePLIN = ParamGen(Erla(meanDRIV, m_in, seed = seedPLIN))
		
		elif simType == 'exp1':
			self.timeSERV = ParamGen(Cons(meanSERV, seed = seedSERV))
			self.timeDRIV = ParamGen(Unif2(rateDRIV, seed = seedDRIV))
			self.timePOUT = ParamGen(Erla(meanDRIV, m_out, seed = seedPOUT))
			self.timePLIN = ParamGen(Erla(meanDRIV, m_in, seed = seedPLIN))

		elif simType == 'unif':
			self.timeSERV = ParamGen(Unif2(meanSERV, seed = seedSERV))
			self.timeDRIV = ParamGen(Unif2(rateDRIV, seed = seedDRIV))
			self.timePOUT = ParamGen(Unif2(meanPOUT, seed = seedPOUT))
			self.timePLIN = ParamGen(Unif2(meanPLIN, seed = seedPLIN))

		elif simType == 'unif1':
			self.timeSERV = ParamGen(Cons(meanSERV, seed = seedSERV))
			self.timeDRIV = ParamGen(Unif2(rateDRIV, seed = seedDRIV))
			self.timePOUT = ParamGen(Unif2(meanPOUT, seed = seedPOUT))
			self.timePLIN = ParamGen(Unif2(meanPLIN, seed = seedPLIN))

		elif simType == 'unif2':
			self.timeSERV = ParamGen(Unif2(meanSERV, seed = seedSERV))
			self.timeDRIV = ParamGen(Cons(rateDRIV, seed = seedDRIV))
			self.timePOUT = ParamGen(Cons(meanPOUT, seed = seedPOUT))
			self.timePLIN = ParamGen(Cons(meanPLIN, seed = seedPLIN))

		else:
			self.timeDRIV = ParamGen(Disc(rateDRIV, seed = seedDRIV))
			self.timePOUT = ParamGen(Tria(ratePOUT, seed = seedPOUT))
			self.timePLIN = ParamGen(Tria(ratePLIN, seed = seedPLIN))

			if simType == 'triag1':
				self.timeSERV = ParamGen(Tri1(rateSERV, seed = seedSERV))

			if simType == 'triag0':
				self.timeSERV = ParamGen(Tri0(rateSERV, seed = seedSERV))

			if simType == 'triag':
				self.timeSERV = ParamGen(Tria(rateSERV, seed = seedSERV))

		if side == 'double':
			self.half_N = N
			self.N = 2 * N
		else:
			assert side == 'single'
			self.N = N

		if angle == 90 and spmatch:
			self.n = spot2blk(N) * CAR_LENGTH
		else:
			self.n = N * LOT_LENGTH

		self.curr = 0.0
		self.start_time = 0.0
		self.eventheap = []
		self.waiting = []
		self.inservice = [None for _ in range(self.N)]
		self.head = None
		self.inCount = 0
		self.outCount = 0	
		self.entry_blocked = self.curr
		self.entry_cleared = self.curr
		self.debug = debug
		self.debug_times = []
		self.debug_speed = []
		self.debug_unit = SIM_UNIT

		self.max_idle = 0.0
		self.max_idle_idx = None
		self.max_pct_idle = 0.0
		self.max_pct_idle_idx = None
		self.prod_time = []
		self.idle_time = []
		self.pct_idle = []	
		self.first_service = None
		self.last_departure = - meanDRIV

		self.wait_out = [{'front_drive': 0.0, 'front_in': 0.0, 'other_drive': 0.0,
						  'back_out': 0.0, 'nearby_out': 0.0, 'oppo_out': 0.0, 
						  'spm_back_out': 0.0, 'total': 0.0, 'veh_count': 0} for _ in range(self.N)]
		
	def add_event(self, event):
		heappush( self.eventheap, event) 

	def state_print(self):

		print (self.curr, self.curr_vehicle.stop_idx, self.curr_typ)
		x = []
		if side == 'single':
			y = [0 for _ in range(spot2blk(self.N))]
		else:
			assert side == 'double'
			y = [0 for _ in range(spot2blk(self.half_N))]
		if not (angle == 0 and mode == 'long'):
			y.append(0)

		for car in self.inservice:
			if car is None:
				x.append(0)
			elif car.status == 3:
				x.append(1)
			elif car.status in [3.5, 4]:
				x.append(2)
			elif car.status == 2:
				x.append(m_in + m_out)
			else:
				assert car.status == 5
				x.append(3)

		car = self.head
		while car is not None:
			if car.status in [1, 6]:
				car.update_loc()
				if np.abs( car.curr_loc % CAR_LENGTH - CAR_LENGTH ) < 1:
					lane_block = ceil(car.curr_loc / CAR_LENGTH)
				else:
					lane_block = floor( car.curr_loc / CAR_LENGTH )
				if lane_block > len(y):
					print ('one vehicle is leaving ...')
				else:
					if y[lane_block - 1] != 0:
						import pdb; pdb.set_trace()
					if car.status == 1:
						y[lane_block - 1] = car.stop_idx
					else:
						y[lane_block - 1] = self.N + 1
			car = car.nex

		print (x)
		print (y)

		if self.N <= 8:
			for idx in range(len(x)):
				if x[idx] == 1:
					print (idx + 1, self.inservice[idx].serv_end)

	def debug_print(self):

		print (self.curr, self.curr_vehicle.stop_idx, self.curr_typ)

		if side == 'single':
			print ('Boarding area ...')
			for car in self.inservice:
				if car is not None:
					print ('Vehicle %s is assigned the spot %s and its current status is %s.' %(car.idx, car.stop_idx, car.status))

		else: 
			assert side == 'double'
			print ('Left boarding area ...')
			for idx in range(self.half_N):
				if self.inservice[2 * idx] is not None:
					car = self.inservice[2 * idx]
					print ('Vehicle %s is assigned the spot %s and its current status is %s.' %(car.idx, car.stop_idx, car.status))

			print ('Right boarding area ...')
			for idx in range(self.half_N):
				if self.inservice[2 * idx + 1] is not None:
					car = self.inservice[2 * idx + 1]
					print ('Vehicle %s is assigned the spot %s and its current status is %s.' %(car.idx, car.stop_idx, car.status))

		print ('Through lane ...')
		car = self.head
		while car != None:
			car.update_loc()
			print ('Vehicle %s is assigned the spot %s; its current status is %s and its current location is %s.' %(car.idx, car.stop_idx, car.status, car.curr_loc) )
			car = car.nex			

	def idle_time_calc(self, curr_vehicle):

		total = self.curr - curr_vehicle.enter_time 
		if spmatch:
			total += meanDRIV
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

			if VEHICLE_IDX is not None and curr_vehicle.idx == VEHICLE_IDX:
				import pdb; pdb.set_trace()

			################################### update system ###################################
			if curr_typ == 'leave_system':
				self.leave_system(curr_vehicle)

			elif curr_typ == 'start_pulling_in':
				self.start_pulling_in(curr_vehicle)

			elif curr_typ == 'start_service':
				self.start_service(curr_vehicle)
			
			elif curr_typ == 'prepare_pulling_out':
				self.prepare_pulling_out(curr_vehicle)

			elif curr_typ == 'finish_pulling_out':
				self.finish_pulling_out(curr_vehicle)

			elif curr_typ == 'add_stop_idx':
				self.add_stop_idx(curr_vehicle)

			else:
				assert curr_typ == 'enter_system'
				self.enter_system(curr_vehicle, debug_idx = VEHICLE_IDX)

			if VEHICLE_IDX is not None and curr_vehicle.idx == VEHICLE_IDX:
				import pdb; pdb.set_trace()

			if self.debug:
				self.state_print()
				if curr_typ == 'enter_system' and curr_vehicle.idx is not None:
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

		self.outCount += 1
		if spmatch and angle == 90:
			curr_vehicle.prod_time += ( (curr_vehicle.dest_to_stop - curr_vehicle.block_idx * CAR_LENGTH) / curr_vehicle.driv )
		elif angle == 90:
			curr_vehicle.prod_time += ( (curr_vehicle.dest_to_stop - (curr_vehicle.stop + dgap) * LOT_LENGTH) / curr_vehicle.driv )
		else:
			curr_vehicle.prod_time += ( (curr_vehicle.dest_to_stop - curr_vehicle.stop * LOT_LENGTH) / curr_vehicle.driv )
		self.idle_time_calc(curr_vehicle)
		self.last_departure = self.curr

		return

	def start_pulling_in(self, curr_vehicle):

		curr_vehicle.update_loc()

		if np.abs( curr_vehicle.curr_loc - curr_vehicle.dest_to_stop ) > 1e-5:
			if curr_vehicle.end_time <= self.curr + SMALL_INTERVAL:
				import pdb; pdb.set_trace() 
			self.add_event( event(curr_vehicle.end_time, curr_vehicle, 'start_pulling_in') )
			return

		heap_len = len(self.eventheap)
		event_holder = []
		next_event = heappop(self.eventheap)
		while self.eventheap != [] and next_event.time - self.curr < 100 * SMALL_INTERVAL:
			assert next_event.time >= self.curr
			if next_event.typ == 'start_pulling_in' and next_event.vehicle.stop > curr_vehicle.stop:
				assert next_event.time + 1e-10 > self.curr
				curr_vehicle.end_time = next_event.time + 1e-10
				curr_vehicle.traj = DLinkedList()
				curr_vehicle.traj.addEnd( changePoint(curr_vehicle.dest_to_stop, self.curr, 0.0) )
				curr_vehicle.traj.addEnd( changePoint(curr_vehicle.dest_to_stop, curr_vehicle.end_time, 'D') )
				for event_visited in event_holder:
					self.add_event( event_visited )
				self.add_event( next_event )
				assert heap_len == len(self.eventheap)
				self.add_event( event(curr_vehicle.end_time, curr_vehicle, 'start_pulling_in') )
				car = curr_vehicle.nex
				while car is not None:
					car.update_loc()
					car.update_traj()
					car = car.nex
				return
			heappush(event_holder, next_event)
			next_event = heappop(self.eventheap)

		for event_visited in event_holder:
			self.add_event( event_visited )
		self.add_event( next_event )
		assert heap_len == len(self.eventheap)

		curr_vehicle.curr_loc = curr_vehicle.dest_to_stop
		delayed = False
		req_time = self.curr

		# if the vehicle has arrived at the desired destination while its enter maneuver blocked by vehicles on the through lane
		# i.e. 1) if its immediate proceeding vehicle is having an exit maneuver at the same spot from the middle lane to the through lane;
		# if (curr_vehicle.prev is not None) and (curr_vehicle.prev.stop_idx == curr_vehicle.stop_idx) and (curr_vehicle.prev.status == 5):
		if (self.inservice[curr_vehicle.stop_idx - 1] is not None):
			assert (self.inservice[curr_vehicle.stop_idx - 1].status == 5)
			assert (curr_vehicle.prev is not None)
			assert (self.inservice[curr_vehicle.stop_idx - 1] == curr_vehicle.prev or side == 'double')
			assert self.inservice[curr_vehicle.stop_idx - 1].pout_end > self.curr
			delayed = True 
			req_time = max( req_time, self.inservice[curr_vehicle.stop_idx - 1].pout_end )

		elif (curr_vehicle.prev is not None) and (curr_vehicle.prev.stop == curr_vehicle.stop) and (curr_vehicle.prev.status == 5):
			assert side == 'double'
			pass

		# or 2) if the immediate proceeding vehicle is stopped.
		elif (curr_vehicle.prev is not None) and (curr_vehicle.prev.status != 2) and (curr_vehicle.prev.curr_loc <= curr_vehicle.dest_to_stop + 1.5 * CAR_LENGTH):
			cp = curr_vehicle.prev.traj.head
			if cp.data.t > self.curr:
				if cp.data.t <= self.curr + 2 * SMALL_INTERVAL:
					delayed = True
					req_time = max( req_time, cp.data.t )
				else:
					import pdb; pdb.set_trace()
			else:
				while (cp.nex is not None) and (cp.nex.data.t <= self.curr):
					cp = cp.nex
				if cp.nex is not None and cp.nex.data.t <= self.curr + SMALL_INTERVAL:
					cp = cp.nex
				assert cp is not None
				if (cp.data.v == 0.0):
					assert cp.nex.data.t > self.curr
					delayed = True
					req_time = max( req_time, cp.nex.data.t )
				elif (spmatch and cp.data.v != 'D' and cp.data.v < rateDRIV and cp.nex.data.t - meanDRIV > self.curr):
					delayed = True
					req_time = max( req_time, cp.nex.data.t - meanDRIV )

		if delayed:
			assert req_time > self.curr
			curr_vehicle.traj = DLinkedList()
			curr_vehicle.traj.addEnd( changePoint(curr_vehicle.dest_to_stop, self.curr, 0.0) )
			curr_vehicle.traj.addEnd( changePoint(curr_vehicle.dest_to_stop, req_time, 'D') )
			curr_vehicle.end_time = req_time
			self.add_event( event(curr_vehicle.end_time, curr_vehicle, 'start_pulling_in') )
			car = curr_vehicle.nex
			while car is not None:
				car.update_loc()
				car.update_traj()
				car = car.nex
			return

		assert self.inservice[curr_vehicle.stop_idx - 1] is None
		self.inservice[curr_vehicle.stop_idx - 1] = curr_vehicle
		curr_vehicle.start_in()

		car = curr_vehicle.nex
		while car != None:
			car.update_loc()
			car.update_traj()
			car = car.nex

		assert curr_vehicle.end_time is not None
		self.add_event( event(curr_vehicle.end_time, curr_vehicle, 'start_service') )
		return

	def start_service(self, curr_vehicle):

		curr_vehicle.start_service()

		if curr_vehicle.prev != None:
			curr_vehicle.prev.nex = curr_vehicle.nex
		else:
			assert self.head == curr_vehicle
			self.head = curr_vehicle.nex

		if curr_vehicle.nex != None:
			car = curr_vehicle.nex
			car.prev = curr_vehicle.prev
			if car.after_plin:
				car.after_plin = False
			while car != None:
				car.update_loc()
				car.update_traj()
				car = car.nex

		curr_vehicle.prev = None
		curr_vehicle.nex = None

		assert curr_vehicle.serv_end > self.curr
		self.add_event( event(curr_vehicle.serv_end, curr_vehicle, 'prepare_pulling_out') )
		return

	def prepare_pulling_out(self, curr_vehicle):

		stop = curr_vehicle.stop
		assert self.inservice[curr_vehicle.stop_idx - 1] == curr_vehicle

		if curr_vehicle.status < 4:
			first_attempt = True
			if curr_vehicle.status == 3:
				assert self.curr == curr_vehicle.serv_end
				self.wait_out[curr_vehicle.stop_idx - 1]['veh_count'] += 1
		else:
			first_attempt = False

		if self.first_service is None:
			self.first_service = self.curr
			first_attempt = False

		if curr_vehicle.idx == VEHICLE_IDX:
			delay_reason = None
			delay_status = None
			delay_speed = None
			delayed, req_time, prev, delay_reason, delay_status, delay_speed = self.check_lane_zero_long(curr_vehicle, first_attempt, delay_reason, delay_status, delay_speed)
		else:
			delayed, req_time, prev = self.check_lane_zero_long(curr_vehicle, first_attempt)

		if delayed:
			if not req_time > self.curr:
				import pdb; pdb.set_trace()
			curr_vehicle.status = 4
			curr_vehicle.end_time = req_time
			self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'prepare_pulling_out') )
			self.wait_out[curr_vehicle.stop_idx - 1]['total'] += (req_time - self.curr)
			if curr_vehicle.idx == VEHICLE_IDX:
				self.debug_times.append( tuple( (delay_reason, delay_status, req_time) ) )
				if delay_reason == 5:
					assert delay_speed is not None
					self.debug_speed.append( delay_speed )
			return

		# vehicles near the entry (i.e. stop <= g_out + 1) might not be able to start an exit maneuver if there is a conflicting entering vehicle 
		# this is a check at the lambda-transitions only (i.e. curr_vehicle.status > 3)
		if spmatch and stop <= g_out + 1 and curr_vehicle.status > 3 and self.eventheap != []:
			broken = False
			car_time = self.curr
			heap_len = len(self.eventheap)
			event_holder = []
			while self.eventheap != []:
				next_event = heappop(self.eventheap)
				heappush(event_holder, next_event)
				if next_event.time <= self.curr + 100 * SMALL_INTERVAL: 
					if next_event.typ == 'enter_system' and self.waiting != []:
						if not free_curb:
							next_spot = - min(self.waiting)
						else:
							next_spot = min(self.waiting)
						if idx2spot(next_spot) > curr_vehicle.stop:
							car_time = next_event.time + 1e-9
							broken = True
							break
						elif idx2spot(next_spot) <= g_out:
							car_time = next_event.time + 1e-9
							broken = True
							break

			for event_visited in event_holder:
				self.add_event( event_visited )
			assert heap_len == len(self.eventheap)

			if broken:
				assert car_time > self.curr
				curr_vehicle.end_time = car_time
				self.add_event( event(curr_vehicle.end_time, curr_vehicle, 'prepare_pulling_out') )
				return

		# this is an optional check and only applies to double-sided design 
		if spmatch and side == 'double' and self.eventheap != []:
			broken = False
			heap_len = len(self.eventheap)
			event_holder = []
			first_service = self.first_service
			if self.first_service is None:
				first_service = self.curr
			while self.eventheap != []:
				next_event = heappop(self.eventheap)
				heappush(event_holder, next_event)
				if curr_vehicle.status == 3 and next_event.time > self.curr + meanDRIV - (self.curr - first_service) % meanDRIV - 0.9 * SMALL_INTERVAL:
					break
				if curr_vehicle.status == 4 and next_event.time > self.curr + meanDRIV - SMALL_INTERVAL:
					break
				if curr_vehicle.stop == next_event.vehicle.stop and curr_vehicle.stop_idx < next_event.vehicle.stop_idx:
					assert curr_vehicle.stop_idx + 1 == next_event.vehicle.stop_idx
					if curr_vehicle.status == 3 and next_event.typ == 'prepare_pulling_out':
						curr_vehicle.status = 4
						curr_vehicle.end_time = self.curr + meanDRIV - (self.curr - first_service) % meanDRIV
						self.add_event( event(curr_vehicle.end_time, curr_vehicle, 'prepare_pulling_out') )
						broken = True
						break	
					if curr_vehicle.status == 4 and (next_event.typ == 'prepare_pulling_out' and (next_event.vehicle.status == 3 or 0 <= next_event.time - self.curr < 100 * SMALL_INTERVAL)):
						curr_vehicle.end_time = self.curr + meanDRIV
						self.add_event( event(curr_vehicle.end_time, curr_vehicle, 'prepare_pulling_out') )
						broken = True
						break
					if curr_vehicle.status == 4 and next_event.typ == 'start_service' and 0 <= next_event.time - self.curr < 100 * SMALL_INTERVAL:
						curr_vehicle.end_time = next_event.time + 1e-9
						self.add_event( event(curr_vehicle.end_time, curr_vehicle, 'prepare_pulling_out') )
						broken = True
						break
				if self.waiting != []:
					if not free_curb:
						next_spot = - min(self.waiting)
					else:
						next_spot = min(self.waiting)
					if control == 'full' and curr_vehicle.stop_idx == 2 and next_event.typ == 'enter_system' and next_spot == 1:
						curr_vehicle.end_time = next_event.time + 1e-9
						self.add_event( event(curr_vehicle.end_time, curr_vehicle, 'prepare_pulling_out') )
						broken = True
						break
			for event_visited in event_holder:
				self.add_event( event_visited )	
			if broken:
				assert heap_len + 1 == len(self.eventheap)
				return
			assert heap_len == len(self.eventheap)

		# Now, the current vehicle is ready to start the exit maneuver	
		curr_vehicle.status = 4
		curr_vehicle.prev = prev
		if prev != None:
			curr_vehicle.nex = prev.nex
			try:
				prev.nex.prev = curr_vehicle
			except:
				pass
			prev.nex = curr_vehicle

		elif self.head != None:
			self.head.prev = curr_vehicle
			curr_vehicle.nex = self.head 
			self.head = curr_vehicle

		else:
			self.head = curr_vehicle

		if spmatch and first_attempt and meanDRIV - (curr_vehicle.serv_end - self.first_service) % meanDRIV > SMALL_INTERVAL:
			assert simType == 'cav'
			assert curr_vehicle.status == 4
			curr_vehicle.status = 5
			curr_vehicle.curr_loc = curr_vehicle.stop * LOT_LENGTH
			curr_vehicle.dest_to_stop = self.n + CAR_LENGTH
			curr_vehicle.pout_start = curr_vehicle.serv_end
			curr_vehicle.pout_time = self.timePOUT.next() - (curr_vehicle.serv_end - self.first_service) % meanDRIV
			curr_vehicle.pout_end = curr_vehicle.pout_time + curr_vehicle.pout_start
			curr_vehicle.prod_time += curr_vehicle.pout_time
			curr_vehicle.update_traj()
		else:
			curr_vehicle.start_out()
		
		assert curr_vehicle.pout_end > self.curr
		self.add_event( event( curr_vehicle.pout_end, curr_vehicle, 'finish_pulling_out') )

		# and we update the trajectories of all upcoming vehicles
		if curr_vehicle.nex is not None:
			car = curr_vehicle.nex
			if car.after_plin:
				if not (prev is not None and prev.status == 2):
					import pdb; pdb.set_trace()
				car.after_plin = False
			while car != None:
				car.update_loc()
				car.update_traj()
				car = car.nex

		# lastly we schedule the replacement vehicles
		new_vehicle = vehicle(self)

		if spmatch:
			assert simType == 'cav'
			enter_time = self.curr + meanDRIV
			if self.curr == curr_vehicle.serv_end or (first_attempt and meanDRIV - (curr_vehicle.serv_end - self.first_service) % meanDRIV > SMALL_INTERVAL):
				enter_time = curr_vehicle.serv_end + meanDRIV - (curr_vehicle.serv_end - self.first_service) % meanDRIV	
			if control == 'full':
				if spot2blk(stop) <= m_out:
					enter_time += meanDRIV * (m_out + 1 - spot2blk(stop))
				enter_time = self.check_enter_zero_long(enter_time, curr_vehicle.stop_idx)
			self.add_event( event(enter_time - 100 * SMALL_INTERVAL, curr_vehicle, 'add_stop_idx') )
			self.add_event( event(enter_time, new_vehicle, 'enter_system') )
		else:
			self.add_event( event( self.curr, new_vehicle, 'enter_system') )
			if not free_curb:
				heappush( self.waiting, (- curr_vehicle.stop_idx) )
			else:
				heappush( self.waiting, curr_vehicle.stop_idx )
		return

	def add_stop_idx(self, curr_vehicle):

		if not free_curb:
			heappush( self.waiting, (- curr_vehicle.stop_idx) )
		else:
			heappush( self.waiting, curr_vehicle.stop_idx )

	def finish_pulling_out(self, curr_vehicle):

		assert curr_vehicle.status == 5
		curr_vehicle.status = 6

		assert self.inservice[curr_vehicle.stop_idx - 1] is not None
		self.inservice[curr_vehicle.stop_idx - 1] = None
		if not (curr_vehicle.end_time is not None and curr_vehicle.end_time > self.curr):
			import pdb; pdb.set_trace()
		self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'leave_system') )				
		return 

	def enter_system(self, curr_vehicle, debug_idx = None):

		assert curr_vehicle.status == 0
		if self.entry_blocked == self.curr:
			self.add_event( event( self.entry_cleared, curr_vehicle, 'enter_system' ) )
			return

		if spmatch:
			assert simType == 'cav'
			if side == 'single' and self.inservice[0] is not None and self.inservice[0].plin_end is not None and self.inservice[0].plin_end + meanDRIV > self.curr + SMALL_INTERVAL:
				self.entry_blocked = self.curr
				self.entry_cleared = self.inservice[0].plin_end + meanDRIV
				self.add_event( event( self.entry_cleared, curr_vehicle, 'enter_system') )
				return	
			if side == 'double' and self.inservice[0] is not None and self.inservice[0].plin_end is not None and self.inservice[0].plin_end + meanDRIV > self.curr + SMALL_INTERVAL:
				self.entry_blocked = self.curr
				self.entry_cleared = self.inservice[0].plin_end + meanDRIV
				self.add_event( event( self.entry_cleared, curr_vehicle, 'enter_system') )
				return
			if side == 'double' and self.inservice[1] is not None and self.inservice[1].plin_end is not None and self.inservice[1].plin_end + meanDRIV > self.curr + SMALL_INTERVAL:
				self.entry_blocked = self.curr
				self.entry_cleared = self.inservice[1].plin_end + meanDRIV
				self.add_event( event( self.entry_cleared, curr_vehicle, 'enter_system') )
				return

			if self.head == None:
				self.inCount += 1
				self.head = curr_vehicle
				heapify(self.waiting)
				if not free_curb:
					curr_vehicle.assign_spot( -heappop(self.waiting) )
				else:
					curr_vehicle.assign_spot( heappop(self.waiting) )
				curr_vehicle.curr_loc = LOT_LENGTH
				if curr_vehicle.stop == 1:						
					curr_vehicle.status = 2
					curr_vehicle.plin_time = self.timePLIN.next()
					curr_vehicle.plin_start = max(0., self.curr - meanDRIV)
					curr_vehicle.plin_end = curr_vehicle.plin_time + curr_vehicle.plin_start
					curr_vehicle.prod_time += curr_vehicle.plin_time
					curr_vehicle.update_traj()
					assert curr_vehicle.end_time != None
					assert self.inservice[curr_vehicle.stop_idx - 1] is None
					self.inservice[curr_vehicle.stop_idx - 1] = curr_vehicle
					self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'start_service') )
					return
				assert curr_vehicle.dest_to_stop >= curr_vehicle.curr_loc
				curr_vehicle.update_traj()
				assert curr_vehicle.end_time != None and curr_vehicle.end_time >= self.curr
				self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'start_pulling_in') )
				return

		# if there is no vehicle on the driving lane
		# the replacement vehicle is assigned to the spot with the largest index
		if self.head == None:
			assert not spmatch
			self.inCount += 1
			self.head = curr_vehicle
			if free_curb and len(self.waiting) >= 3:
				print (self.waiting)
				import pdb; pdb.set_trace()
			heapify(self.waiting)
			if not free_curb:
				curr_vehicle.assign_spot( -heappop(self.waiting) )
			else:
				curr_vehicle.assign_spot( heappop(self.waiting) )
			curr_vehicle.update_traj()
			assert curr_vehicle.end_time != None and curr_vehicle.end_time >= self.curr
			self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'start_pulling_in') )
			return

		# if there are vehicles on the driving lane, find the last one
		req_time = self.curr
		car = self.head
		while car.nex is not None:
			car = car.nex	
		car.update_loc()

		# with 0-degree long spots, 
		# if the last one is making an enter maneuver in to spot 2 (single) or spots 3 & 4 (double) 
		# i.e. the loc of the last vehicle is right at the CAR_LENGTH (i.e. first lane block occupied)
		# the replacement vehicle can enter and be assigned to spot 1 under both access control
		if angle == 0 and mode == 'long' and (not free_curb) and car.curr_loc == CAR_LENGTH and car.status == 2 and (not spmatch or car.plin_start >= self.curr - meanDRIV - 100 * SMALL_INTERVAL):
			stop_idx = None
			if (-2 in self.waiting and side == 'double' and (control == 'full' or -2 == min(self.waiting))):
				stop_idx = 2 
				self.waiting.remove(-2)
			elif (-1 in self.waiting and (control == 'full' or -1 == min(self.waiting))):
				stop_idx = 1
				self.waiting.remove(-1)
			
			if stop_idx is not None:
				if np.abs(car.plin_start - self.curr) < 100 * SMALL_INTERVAL:
					self.entry_blocked = self.curr
					self.entry_cleared = self.curr + meanDRIV
					self.add_event( event( self.entry_cleared, curr_vehicle, 'enter_system') )			
					heappush(self.waiting, (-stop_idx))				
					return

				if free_curb and len(self.waiting) >= 3:
					print (self.waiting)
					import pdb;pdb.set_trace()
				self.inCount += 1
				curr_vehicle.prev = car
				car.nex = curr_vehicle
				curr_vehicle.assign_spot( stop_idx )
				assert curr_vehicle.stop == 1
				assert self.inservice[curr_vehicle.stop_idx - 1] is None
				if not spmatch:
					curr_vehicle.update_traj()
					assert curr_vehicle.end_time != None
					self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'start_pulling_in') )
					return
				else:
					assert simType == 'cav'
					curr_vehicle.status = 2
					curr_vehicle.curr_loc = CAR_LENGTH
					curr_vehicle.plin_time = self.timePLIN.next()
					curr_vehicle.plin_start = max(0., self.curr - meanDRIV)
					curr_vehicle.plin_end = curr_vehicle.plin_time + curr_vehicle.plin_start
					curr_vehicle.prod_time += curr_vehicle.plin_time
					curr_vehicle.update_traj()
					assert curr_vehicle.end_time != None
					self.inservice[curr_vehicle.stop_idx - 1] = curr_vehicle	
					self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'start_service') )
					return

		# if the last one is within CAR_LENGTH (i.e. first lane block occupied)
		# the replacement vehicle cannot enter under either access control
		if car.curr_loc <= CAR_LENGTH + SMALL_INTERVAL:
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

			elif angle == 0 and mode == 'long' and side == 'double' and control == 'full' and (not free_curb) and car.stop == 1 and car.status == 5:
				stop_idx = None
				if car.stop_idx == 1 and -2 in self.waiting:
					stop_idx = 2 
					self.waiting.remove(-2)
				elif car.stop_idx == 2 and -1 in self.waiting:
					stop_idx = 1
					self.waiting.remove(-1)
			
				if stop_idx is not None:
					if free_curb and len(self.waiting) >= 3:
						print (self.waiting)
						import pdb; pdb.set_trace()
					self.inCount += 1
					curr_vehicle.prev = car
					car.nex = curr_vehicle
					curr_vehicle.assign_spot( stop_idx )
					assert curr_vehicle.stop == 1
					assert self.inservice[curr_vehicle.stop_idx - 1] is None
					if not spmatch:
						curr_vehicle.update_traj()
						assert curr_vehicle.end_time != None
						self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'start_pulling_in') )
						return
					else:
						assert simType == 'cav'
						curr_vehicle.status = 2
						curr_vehicle.curr_loc = CAR_LENGTH
						curr_vehicle.plin_time = self.timePLIN.next()
						curr_vehicle.plin_start = max(0., self.curr - meanDRIV)
						curr_vehicle.plin_end = curr_vehicle.plin_time + curr_vehicle.plin_start
						curr_vehicle.prod_time += curr_vehicle.plin_time
						curr_vehicle.update_traj()
						assert curr_vehicle.end_time != None
						self.inservice[curr_vehicle.stop_idx - 1] = curr_vehicle	
						self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'start_service') )
						return
				else:
					assert car_time > self.curr
					req_time = max( req_time, car_time)	 

			else:
				assert car_time > self.curr
				req_time = max( req_time, car_time)

		if spmatch and req_time <= self.curr:
			if car.curr_loc < LOT_LENGTH + CAR_LENGTH - SMALL_INTERVAL:
				if car.status not in [1, 6]:
					print ('!!!!!!!!!!!!!!!!')
					import pdb; pdb.set_trace()
				assert car.dest_to_stop >= LOT_LENGTH + CAR_LENGTH
				car_time = car.calc_time(LOT_LENGTH + CAR_LENGTH)
				req_time = max( req_time, car_time )

		if req_time > self.curr:
			self.entry_blocked = self.curr
			self.entry_cleared = req_time
			self.add_event( event( self.entry_cleared, curr_vehicle, 'enter_system') )
			return

		if debug_idx is not None and curr_vehicle.idx == debug_idx:
			import pdb; pdb.set_trace()

		if control == 'partial':
			if free_curb and len(self.waiting) >= 3:
				print (self.waiting)
				import pdb; pdb.set_trace()
			self.inCount += 1
			curr_vehicle.prev = car
			car.nex = curr_vehicle
			heapify(self.waiting)
			if not free_curb:
				curr_vehicle.assign_spot( -heappop(self.waiting) )
			else:
				curr_vehicle.assign_spot( heappop(self.waiting) )
			if spmatch:
				assert simType == 'cav'
				assert curr_vehicle.prev.curr_loc >= LOT_LENGTH + CAR_LENGTH - SMALL_INTERVAL
				curr_vehicle.curr_loc = CAR_LENGTH
				if spot2blk(curr_vehicle.stop) == 1:
					curr_vehicle.status = 2
					curr_vehicle.plin_time = self.timePLIN.next()
					curr_vehicle.plin_start = max(0., self.curr - meanDRIV)
					curr_vehicle.plin_end = curr_vehicle.plin_time + curr_vehicle.plin_start
					curr_vehicle.prod_time += curr_vehicle.plin_time
					curr_vehicle.update_traj()
					assert curr_vehicle.end_time != None
					assert self.inservice[curr_vehicle.stop_idx - 1] is None
					self.inservice[curr_vehicle.stop_idx - 1] = curr_vehicle
					self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'start_service') )
					return
				else:
					assert curr_vehicle.dest_to_stop >= curr_vehicle.curr_loc
			curr_vehicle.update_traj()
			assert curr_vehicle.end_time != None and curr_vehicle.end_time >= self.curr
			self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'start_pulling_in') )
			return

		assert control == 'full'
		assert curr_vehicle.curr_loc == 0.0
		if spmatch:
			curr_vehicle.curr_loc = CAR_LENGTH

		# car is the last vehicle on the driving lane 
		# and also the prev for the replacement vehicle if the latter can enter
		j_new = []
		last = car

		# j_new will include the spots that the replacement vehicle can head to 
		# without being blocked or delayed by the last vehicle on the lane in expectation
		for j in sorted(self.waiting, reverse = True):
			j = - j
			assert j > 0
			if side == 'double':
				J = idx2spot(j)
			else:
				assert side == 'single'
				J = j

			if (J < car.stop) or (car.status == 6):
				j_new.append(j)

			elif car.status == 2:
				# K_in with K = car.j and J = idx2spot(j)
				assert car.stop_idx != j
				assert car.plin_start <= self.curr
				assert not J == car.stop == 1
				assert car.stop >= 3
				assert car.curr_loc == (car.stop - 1) * LOT_LENGTH >= curr_vehicle.curr_loc
				if ( (car.curr_loc - curr_vehicle.curr_loc) / rateDRIV >= max(0.0, meanPLIN - (self.curr - car.plin_start)) ):
					j_new.append(j)

			elif car.status == 5:
				# K_out with K = car.j and J = idx2spot(j)
				assert car.pout_start <= self.curr <= car.pout_end
				assert (car.stop - 1) * LOT_LENGTH >= curr_vehicle.curr_loc
				if side == 'double' and J == car.stop and j != car.stop_idx:
					j_new.append(j)
				elif spmatch and ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV >= car.pout_end - self.curr - SMALL_INTERVAL:
					j_new.append(j)
				elif (not spmatch) and ( ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV >= max(0.0, meanPOUT - (self.curr - car.pout_start)) ):
					j_new.append(j)

			else:
				assert car.status == 1
				# I_in with K = car.j and J = idx2spot(j)
				assert car.stop_idx != j
				assert car.stop >= 3
				assert car.curr_loc >= curr_vehicle.curr_loc
				if ( (car.curr_loc - curr_vehicle.curr_loc) / rateDRIV >= meanPLIN - 100 * SMALL_INTERVAL ):
					j_new.append(j)

		if j_new != []:
			assert j_new[-1] == max(j_new)
			car_time = 0.0

		else:
			j = - sorted(self.waiting, reverse = True)[0]
			if side == 'double':
				assert (idx2spot(j) >= car.stop)
			else:
				assert side == 'single'
				assert (j >= car.stop)

			if car.status == 2:
				assert car.stop_idx != j
				assert car.plin_start <= self.curr
				assert ( (car.curr_loc - curr_vehicle.curr_loc) / rateDRIV < meanPLIN - (self.curr - car.plin_start) )
				car_time = meanPLIN - (self.curr - car.plin_start) - (car.curr_loc - curr_vehicle.curr_loc) / rateDRIV

			elif car.status == 5:
				assert car.pout_start <= self.curr
				if spmatch:
					assert ( ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV < car.pout_end - self.curr )
					car_time = car.pout_end - self.curr - ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV
				if not spmatch:
					assert ( ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV < meanPOUT - (self.curr - car.pout_start)  )
					car_time = meanPOUT - (self.curr - car.pout_start) - ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV

			else:
				assert car.status == 1
				assert car.stop_idx != j
				assert ( (car.curr_loc - curr_vehicle.curr_loc) / rateDRIV < meanPLIN )
				car_time = meanPLIN - (car.curr_loc - curr_vehicle.curr_loc) / rateDRIV

		while (j_new != []) and (car.prev is not None):

			car = car.prev
			car.update_loc()
			for idx in range(len(j_new)):
				j = j_new[idx]
				if side == 'double':
					J = idx2spot(j)
				else:
					assert side == 'single'
					J = j

				if (J < car.stop) or (car.status == 6):
					pass

				elif car.status == 2:
					# K_in with K = car.j and J = idx2spot(j)
					assert car.stop_idx != j
					assert car.plin_start <= self.curr
					assert car.stop >= 3
					assert car.curr_loc == (car.stop - 1) * LOT_LENGTH >= curr_vehicle.curr_loc
					if ( (car.curr_loc - curr_vehicle.curr_loc) / rateDRIV < max(0.0, meanPLIN - (self.curr - car.plin_start)) ):
						j_new = j_new[:idx]
						if idx == 0:
							car_time = max(car_time, meanPLIN - (self.curr - car.plin_start) - (car.curr_loc - curr_vehicle.curr_loc) / rateDRIV)
						break

				elif car.status == 5:
					# K_out with K = car.j and J = idx2spot(j)
					assert car.pout_start <= self.curr <= car.pout_end
					assert (car.stop - 1) * LOT_LENGTH >= curr_vehicle.curr_loc
					if side == 'double' and J == car.stop and j != car.stop_idx:
						pass
					elif spmatch and ( ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV < car.pout_end - self.curr ):
						j_new = j_new[:idx]
						if idx == 0:
							car_time = max(car_time, car.pout_end - self.curr - ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV) 
						break
					elif not spmatch and ( ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV < max(0.0, meanPOUT - (self.curr - car.pout_start)) ):
						j_new = j_new[:idx]
						if idx == 0:
							car_time = max(car_time, meanPOUT - (self.curr - car.pout_start) - ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV)
						break

				else:
					assert car.status == 1
					# I_in with K = car.j and J = idx2spot(j)
					assert car.stop_idx != j
					assert car.stop >= 3
					assert car.curr_loc >= curr_vehicle.curr_loc
					if ( (car.curr_loc - curr_vehicle.curr_loc) / rateDRIV < meanPLIN ):
						j_new = j_new[:idx]
						if (idx == 0):
							car_time = max(car_time, meanPLIN - (car.curr_loc - curr_vehicle.curr_loc) / rateDRIV)
						break

		if j_new == []:
			car_time = self.curr + car_time
			if car_time == self.curr:
				car_time += SMALL_INTERVAL
			if car_time <= self.curr:
				import pdb; pdb.set_trace()
			curr_vehicle.curr_loc = 0.0
			self.add_event( event( car_time, curr_vehicle, 'enter_system') )
			return

		# car.prev is None
		# i.e. there is at least one spot where a replacement vehicle can head to 
		# without being blocked or delayed by any vehicle already on the lane
		# if there are multiple, choose the largest 
		self.inCount += 1
		assert j_new[-1] == max(j_new)
		assert j_new[0] == min(j_new)
		if not free_curb:
			curr_vehicle.assign_spot( j_new[-1] )
			assert (- j_new[-1]) in self.waiting
			self.waiting.remove( - j_new[-1] )
		else:
			curr_vehicle.assign_spot( j_new[0] )
			assert j_new[0] in self.waiting
			self.waiting.remove( j_new[0] )
			
		curr_vehicle.prev = last
		last.nex = curr_vehicle
		if spmatch:
			assert simType == 'cav'
			assert curr_vehicle.prev.curr_loc >= LOT_LENGTH + CAR_LENGTH - SMALL_INTERVAL
			assert curr_vehicle.curr_loc == CAR_LENGTH
			if spot2blk(curr_vehicle.stop) == 1:
				curr_vehicle.status = 2
				curr_vehicle.plin_time = self.timePLIN.next()
				curr_vehicle.plin_start = max(0., self.curr - meanDRIV)
				curr_vehicle.plin_end = curr_vehicle.plin_time + curr_vehicle.plin_start
				curr_vehicle.prod_time += curr_vehicle.plin_time
				curr_vehicle.update_traj()
				assert curr_vehicle.end_time != None
				assert self.inservice[curr_vehicle.stop_idx - 1] is None
				self.inservice[curr_vehicle.stop_idx - 1] = curr_vehicle
				self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'start_service') )
				return
			else:
				assert curr_vehicle.dest_to_stop >= curr_vehicle.curr_loc
		if not spmatch:
			assert curr_vehicle.curr_loc == 0.0
		curr_vehicle.update_traj()
		assert curr_vehicle.end_time != None and curr_vehicle.end_time >= self.curr
		self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'start_pulling_in') )
		return

	def check_enter_zero_long(self, enter_time, stop_idx):

		assert spmatch and control == 'full'
		car = self.head
		while car is not None:
			if (stop_idx <= car.stop_idx) or (car.status == 6):
				pass
			elif car.status == 2 and car.stop == 1:
				assert car.plin_start <= self.curr
				enter_time = max(enter_time, car.plin_end + meanDRIV)
			elif car.status == 2:
				assert car.plin_start <= self.curr
				if (car.stop - 1) * LOT_LENGTH / rateDRIV < max(0.0, meanPLIN - (self.curr - car.plin_start)):
					enter_time = max(enter_time, car.plin_start + meanPLIN - (car.stop - 1) * LOT_LENGTH / rateDRIV)
			elif car.status == 5 and car.stop == 1:
				assert car.pout_start <= self.curr <= car.pout_end
				enter_time = max(enter_time, car.pout_end + meanDRIV)
			elif car.status == 5:
				assert car.pout_start <= self.curr <= car.pout_end
				if ((car.stop - 1) * LOT_LENGTH - CAR_LENGTH) / rateDRIV < car.pout_end - self.curr:
					enter_time = max(enter_time, car.pout_end - ((car.stop - 1) * LOT_LENGTH - CAR_LENGTH) / rateDRIV)
			else:
				assert car.status == 1
				assert car.stop >= 2 
				if enter_time < self.curr + meanPLIN - car.curr_loc / rateDRIV - 10 * SMALL_INTERVAL:
					enter_time = self.curr + meanPLIN - car.curr_loc / rateDRIV
			car = car.nex
		return enter_time

	def check_lane_zero_long(self, curr_vehicle, first_attempt, delay_reason = None, delay_status = None, delay_speed = None, curr_time = None):

		assert angle == 0 and mode == 'long'

		delayed = False
		if curr_time is None:
			curr_time = self.curr
			if spmatch and first_attempt:
				est_pout_start = curr_vehicle.serv_end
			else:
				est_pout_start = curr_time
		else:
			assert curr_time > self.curr
			est_pout_start = curr_time
		req_time = curr_time

		car = self.head
		prev = None
		stopped = False

		stop = curr_vehicle.stop
		idx = curr_vehicle.idx

		temp_delay = {'front_drive': 0.0, 'front_in': 0.0, 'other_drive': 0.0,
						  'back_out': 0.0, 'nearby_out': 0.0, 'oppo_out': 0.0, 
						  'spm_back_out': 0.0, 'total': 0.0, 'veh_count': 0} 
			
		while car != None:

			car.update_loc()
			if car.idx == idx:
				pass

			elif car.curr_loc >= stop * LOT_LENGTH + CAR_LENGTH - SMALL_INTERVAL:
				prev = car

			elif spmatch and car.stop == stop and car.status == 6:
				prev = car

			elif car.status == 2 and car.stop == stop + 1:
				if not ((side == 'single' and stop <= self.N - 1) or (side == 'double' and stop <= self.half_N - 1)):
					import pdb; pdb.set_trace()
				if not (self.inservice[car.stop_idx - 1] == car):
					import pdb; pdb.set_trace()
				assert car.plin_end >= self.curr > car.plin_start
				if car.plin_end > curr_time:
					delayed = True
					if not car.plin_end > curr_time:
						import pdb; pdb.set_trace()
					if idx == VEHICLE_IDX and car.plin_end > req_time:
						delay_reason = 2
						delay_status = 2
					req_time = max( req_time, car.plin_end)
					temp_delay['front_in'] = max(temp_delay['front_in'], car.plin_end - self.curr)
				else:
					prev = car

			elif car.status == 5: 
				stopped = False
				assert car.stop <= stop
				assert car.dest_to_stop >= stop * LOT_LENGTH + CAR_LENGTH
				if car.stop == stop and car.pout_end - curr_time > SMALL_INTERVAL:
					delayed = True
					assert car.pout_end > curr_time
					if curr_vehicle.idx == VEHICLE_IDX and car.pout_end > req_time:
						delay_reason = 1 
						delay_status = 5
					req_time = max(req_time, car.pout_end)
					temp_delay['oppo_out'] = car.pout_end - self.curr
					break

				if car.stop_idx in out_range(curr_vehicle.stop_idx, self.N):
					if spmatch and first_attempt and np.abs( car.pout_end - (est_pout_start + meanPOUT - (est_pout_start - self.first_service) % meanDRIV) ) < 100 * SMALL_INTERVAL:
						assert simType == 'cav'
						pass
					elif spmatch and 0 <= curr_time - car.pout_start < SMALL_INTERVAL:
						assert simType == 'cav'
						pass
					elif car.pout_end - curr_time > SMALL_INTERVAL:
						delayed = True
						assert car.pout_end > curr_time
						if curr_vehicle.idx == VEHICLE_IDX and car.pout_end > req_time:
							delay_reason = 1 
							delay_status = 5
						req_time = max(req_time, car.pout_end)
						temp_delay['nearby_out'] = max(temp_delay['nearby_out'], car.pout_end - self.curr)

				if spmatch:
					car_time = car.pout_end + ((stop - 1) * LOT_LENGTH - car.curr_loc) / rateDRIV
					if (first_attempt and car_time < est_pout_start + meanPOUT - (est_pout_start - self.first_service) % meanDRIV - 100 * SMALL_INTERVAL) or ((not first_attempt) and car_time < curr_time + meanPOUT - 100 * SMALL_INTERVAL):
						delayed = True
						car_time = car.calc_time( stop * LOT_LENGTH + CAR_LENGTH )
						assert car_time > curr_time
						if idx == VEHICLE_IDX and car_time > req_time:
							delay_reason = 3
							delay_status = 5
						req_time = max( req_time, car_time )
						temp_delay['spm_back_out'] = max(temp_delay['spm_back_out'], car_time - self.curr)

				else:
					car_time = car.pout_start + meanPOUT + ((stop - 1) * LOT_LENGTH - car.curr_loc) / rateDRIV
					if car_time < curr_time + meanPOUT:
						delayed = True
						car_time = car.calc_time( stop * LOT_LENGTH + CAR_LENGTH )
						assert car_time > curr_time
						if idx == VEHICLE_IDX and car_time > req_time:
							delay_reason = 3
							delay_status = 5
						req_time = max( req_time, car_time )
						temp_delay['back_out'] = max(temp_delay['back_out'], car_time - self.curr)

			elif car.status == 2:
				assert (car.stop < stop) or (car.stop == stop and side == 'double')
				if not spmatch and car.stop == stop:
					break
				stopped = True

			elif (car.stop == stop) and (car.status == 1):
				assert (car.curr_loc <= car.dest_to_stop == (stop - 1) * LOT_LENGTH)
				stopped = True			 

			elif (car.stop < stop) and (car.status == 1):
				if stop < 2:
					import pdb; pdb.set_trace()
				assert (car.curr_loc <= car.dest_to_stop <= (stop - 2) * LOT_LENGTH)
				stopped = True

			elif (car.stop == stop + 1) and (car.status == 1):
				assert car.dest_to_stop == stop * LOT_LENGTH
				car_time = self.curr + ((stop - 1) * LOT_LENGTH - car.curr_loc) / rateDRIV

				if stopped and car.prev.status == 1:
					assert car.prev.end_time is not None
					assert car.prev.stop < stop or (car.prev.stop == stop and side == 'double')
					assert (car.prev.stop - 1) * LOT_LENGTH >= car.curr_loc
					car_time = max(meanPLIN + car.prev.end_time, self.curr + ((car.prev.stop - 1) * LOT_LENGTH - car.curr_loc) / rateDRIV) + (stop - car.prev.stop) * LOT_LENGTH / rateDRIV 
				elif stopped:
					assert car.prev.status == 2
					car_time = car.calc_time( (stop - 1) * LOT_LENGTH )

				if ((not spmatch or not first_attempt) and car_time < curr_time + meanPOUT - 100 * SMALL_INTERVAL) or (spmatch and first_attempt and car_time < est_pout_start + meanPOUT - (est_pout_start - self.first_service) % meanDRIV - 100 * SMALL_INTERVAL):
					delayed = True
					if not car.end_time > curr_time - 3 * SMALL_INTERVAL:
						if curr_time > self.curr:
							delayed = True
							break
						import pdb; pdb.set_trace()
					if car.end_time <= curr_time:
						car_time = car.end_time + SMALL_INTERVAL
						if car_time <= curr_time:
							car_time += SMALL_INTERVAL
					else:
						car_time = car.end_time
					if idx == VEHICLE_IDX and car_time > req_time:
						delay_reason = 4
						delay_status = 1
					req_time = max( req_time, car_time )
					temp_delay['front_drive'] = max(temp_delay['front_drive'], car_time - self.curr)

				if stopped:
					break

			else:
				assert (car.status in [1, 6])
				assert car.dest_to_stop >= stop * LOT_LENGTH + CAR_LENGTH
				
				if self.curr < car.calc_time(stop * LOT_LENGTH + CAR_LENGTH) < curr_time:
					prev = car

				else:
					car_time = self.curr + ((stop - 1) * LOT_LENGTH - car.curr_loc) / rateDRIV
					if stopped and car.prev.status == 1:
						assert car.prev.end_time is not None
						assert car.prev.stop < stop or (car.prev.stop == stop and side == 'double')
						assert (car.prev.stop - 1) * LOT_LENGTH >= car.curr_loc
						car_time = max(meanPLIN + car.prev.end_time, self.curr + ((car.prev.stop - 1) * LOT_LENGTH - car.curr_loc) / rateDRIV) + (stop - car.prev.stop) * LOT_LENGTH / rateDRIV 		
					elif stopped:
						assert car.prev.status == 2
						car_time = max(car.prev.plin_end, self.curr + ((car.prev.stop - 1) * LOT_LENGTH - car.curr_loc) / rateDRIV) + (stop - car.prev.stop) * LOT_LENGTH / rateDRIV 

					if ((not spmatch or not first_attempt) and car_time < curr_time + meanPOUT - 100 * SMALL_INTERVAL) or (spmatch and first_attempt and car_time < est_pout_start + meanPOUT - (est_pout_start - self.first_service) % meanDRIV - 100 * SMALL_INTERVAL):
						delayed = True
						car_time = car.calc_time( stop * LOT_LENGTH + CAR_LENGTH )
						if not car_time > curr_time:
							if control == 'full' and curr_time > self.curr and car.dest_to_stop >= stop * LOT_LENGTH + CAR_LENGTH:
								pass
							else:
								import pdb; pdb.set_trace()
						if idx == VEHICLE_IDX and car_time > req_time:
							delay_reason = 5
							delay_status = car.status
							cp = car.traj.head
							if cp.data.t > curr_time:
								import pdb; pdb.set_trace()
							while cp.nex is not None and cp.nex.data.t <= curr_time:
								cp = cp.nex
							if cp.nex is None:
								import pdb; pdb.set_trace()
							assert cp.nex.data.t > curr_time >= cp.data.t 
							delay_speed = cp.data.v
						req_time = max( req_time, car_time )	
						temp_delay['other_drive'] = max(temp_delay['other_drive'], car_time - self.curr)
				
					if stopped:
						break					

			car = car.nex

		for item in temp_delay:
			self.wait_out[curr_vehicle.stop_idx - 1][item] += temp_delay[item] 
		
		if idx == VEHICLE_IDX: 
			return (tuple( (delayed, req_time, prev, delay_reason, delay_status, delay_speed) ))
		else:
			return (tuple( (delayed, req_time, prev) ))

class vehicle():
	
	def __init__(self, sys, getin = False, stop_idx = None):

		self.sys = sys
		self.driv = self.sys.timeDRIV.next() 
		self.status = 0

		self.idx = None
		self.stop_idx = None
		self.stop = None
		self.block_idx = None
		self.curr_loc = None
		self.dest_to_stop = None
		self.enter_time = None
		self.end_time = None
		self.prod_time = 0.0	

		self.after_plin = False
		self.prev = None
		self.nex = None
		self.traj = None

		self.plin_start = None
		self.plin_time = None
		self.plin_end = None

		self.pout_start = None
		self.pout_time = None
		self.pout_end = None

		self.serv_start = None
		self.serv_time = None
		self.serv_end = None

		if getin:
			assert stop_idx is not None
			self.assign_spot(stop_idx)
			self.status = 2
			self.start_service()
		else:
			assert stop_idx is None
			self.curr_loc = 0.0

	def assign_spot(self, stop_idx):

		assert self.status == 0
		self.status = 1
		self.idx = self.sys.inCount
		self.stop_idx = stop_idx
		self.enter_time = self.sys.curr
		self.stop = idx2spot(stop_idx)
		self.block_idx = self.stop
		if (angle == 0) and (mode == 'long'):
			self.dest_to_stop = (self.stop - 1) * LOT_LENGTH
		elif (angle == 0):
			assert (mode == 'short')
			self.dest_to_stop = self.stop * LOT_LENGTH
		else:
			assert (angle == 90)
			self.block_idx = spot2blk(self.stop)
			if spmatch:
				self.dest_to_stop = self.block_idx * CAR_LENGTH
			else:
				self.dest_to_stop = (self.stop + dgap) * LOT_LENGTH

	def start_in(self):

		assert self.curr_loc == self.dest_to_stop
		assert (self.status == 1)
		self.status = 2
		self.plin_time = self.sys.timePLIN.next()
		self.plin_end = self.plin_time + self.sys.curr
		self.plin_start = self.sys.curr
		self.prod_time += (self.dest_to_stop / self.driv + self.plin_time)
		self.update_traj()

	def start_service(self):

		assert (self.status == 2)
		self.status = 3
		if spmatch:
			try:
				self.serv_time = self.sys.service_times[str(self.stop_idx)][0]
				self.sys.service_times[str(self.stop_idx)] = self.sys.service_times[str(self.stop_idx)][1:]
			except:
				self.serv_time = self.sys.timeSERV.next()
		else:
			self.serv_time = self.sys.timeSERV.next()

		self.serv_end = self.serv_time + self.sys.curr
		self.serv_start = self.sys.curr
		self.prod_time += self.serv_time
		self.traj = None

	def start_out(self):

		assert self.status == 4
		self.status = 5
		if spmatch and angle == 90:
			self.curr_loc = self.block_idx * CAR_LENGTH
		elif angle == 90:
			self.curr_loc = (self.stop + dgap) * LOT_LENGTH
		else:
			self.curr_loc = self.stop * LOT_LENGTH
		self.dest_to_stop = self.sys.n + CAR_LENGTH
		self.pout_time = self.sys.timePOUT.next()
		self.pout_end = self.pout_time + self.sys.curr
		self.pout_start = self.sys.curr
		self.prod_time += self.pout_time
		self.update_traj()

	def calc_time(self, loc):

		if self.curr_loc > loc:
			if control == 'full':
				pass
			else:
				import pdb; pdb.set_trace()

		cp = self.traj.head

		if spmatch:
			if self.curr_loc == loc:
				return self.sys.curr
			while cp.nex != None and cp.nex.data.x < loc:
				cp = cp.nex
			if cp.data.v == 0.0:
				if cp.nex.data.x == loc:
					return cp.nex.data.t
				import pdb; pdb.set_trace()
			if cp.nex is None:
				assert cp.data.v == 'D'
				return cp.data.t
			return cp.data.t + (loc - cp.data.x) / cp.data.v

		# if the curr_loc is the target loc,
		# then find the last time that it stays here
		# i.e. the first time that it starts to have nonzero speed.
		# NOTE: nonzero speed can be either positive or 'D'.
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

		# since all trajectories end with speed 'D' i.e. nonzero
		# we can expect that all vehicles with curr_loc == loc already considered
		# except for those starting with cp.data.v > 0.0,
		# which means the vehicle is arrived at the loc
		while cp.nex != None:
			if cp.nex.data.x <= loc:
				cp = cp.nex
			else:
				break

		if cp.data.v == 0.0:
			import pdb;pdb.set_trace()

		# now we start with the changepoint s.t.
		# cp.data.x <= loc < cp.nex.data.x
		# Note: if cp.data.v == 0.0, then cp.data.x = cp.nex.data.x
		# given that cp.data.x <= loc
		# we can expect cp.nex.data.x <= loc thus the iteration will go on 
		# which explains the strict inequality between loc and cp.nex.data.x
		if cp.nex == None:
			assert cp.data.v == 'D'
			return cp.data.t

		return cp.data.t + (loc - cp.data.x) / cp.data.v
	
	def update_loc(self):

		assert (self.status != 3) and (self.status != 4)

		if (self.status == 2):
			if spmatch and spot2blk(self.stop) == 1:
				assert self.curr_loc == CAR_LENGTH
			else:
				assert self.curr_loc == self.dest_to_stop
			return

		if (self.status == 5):
			if spmatch and angle == 90:
				assert self.curr_loc == self.block_idx * CAR_LENGTH
			elif angle == 90:
				assert self.curr_loc == (self.stop + dgap) * LOT_LENGTH
			else:
				assert self.curr_loc == self.stop * LOT_LENGTH
			return

		try:
			cp = self.traj.head
		except:
			import pdb; pdb.set_trace()

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
			if not (self.curr_loc - self.dest_to_stop < 1.5e-04):
				import pdb; pdb.set_trace()
			self.curr_loc = self.dest_to_stop
		return

	def update_traj(self):

		traj = self.traj
		end_time = self.end_time
		assert self.status in [1, 2, 5, 6]

		self.traj = DLinkedList()

		if self.status == 2:
			# i.e. the vehicle has started pulling in but has not started service
			if spmatch and spot2blk(self.stop) == 1:
				assert self.curr_loc == CAR_LENGTH
			else:
				assert (self.curr_loc == self.dest_to_stop)
			if self.plin_end > self.sys.curr:
				self.end_time = self.plin_end
			else:
				assert self.end_time == self.plin_end
			if end_time is not None:
				assert self.end_time >= end_time
			assert self.end_time >= self.sys.curr
			self.traj.addEnd( changePoint(self.curr_loc, self.sys.curr, 0.0) )
			self.traj.addEnd( changePoint(self.curr_loc, self.end_time, 'D') )
			return

		if self.curr_loc == self.dest_to_stop and self.status == 6:
			assert self.end_time is not None
			assert (np.abs(self.end_time - self.sys.curr) < SMALL_INTERVAL)
			if end_time is not None:
				assert self.end_time >= end_time
			assert self.end_time >= self.sys.curr
			self.traj.addEnd( changePoint(self.curr_loc, self.sys.curr, 0.0))
			self.traj.addEnd( changePoint(self.curr_loc, self.end_time, 'D') )
			return

		if self.curr_loc == self.dest_to_stop and self.status == 1:
			if self.end_time is None:
				assert (self.stop == 1) or (spmatch and self.stop == 2) or (angle == 90 and self.block_idx == 1) or (angle == 90 and spmatch and self.block_idx == 2)
				self.end_time = self.sys.curr
			if end_time is not None:
				assert self.end_time >= end_time
			if not self.end_time >= self.sys.curr:
				if not np.abs(self.end_time - self.sys.curr) < SMALL_INTERVAL:
					import pdb; pdb.set_trace()
				self.end_time = self.sys.curr
			self.traj.addEnd( changePoint(self.curr_loc, self.sys.curr, 0.0))
			self.traj.addEnd( changePoint(self.curr_loc, self.end_time, 'D') )
			return

		start_t = self.sys.curr

		if self.status == 5:
			# i.e. the vehicle has started pulling out
			self.update_loc()
			assert self.traj.head == None
			if self.pout_end > self.sys.curr:
				self.traj.addEnd( changePoint(self.curr_loc, self.sys.curr, 0.0) )
				start_t = self.pout_end
			
		if self.prev is None:
			try:
				assert self.sys.head == self
			except:
				print ('Line 541: this car does not have a prev while self.sys.head is not itself!!!')
				import pdb; pdb.set_trace()

			self.traj.addEnd( changePoint(self.curr_loc, start_t, self.driv) )
			self.end_time = start_t + (self.dest_to_stop - self.curr_loc) / self.driv
			if (end_time is not None) and (self.end_time < end_time):
				if spmatch and self.status == 5 and self.idx == self.sys.curr_vehicle.idx:
					import pdb; pdb.set_trace()
					if not ( np.abs(end_time - self.end_time - (self.sys.curr - self.sys.first_service) % meanDRIV) < SMALL_INTERVAL or (self.serv_end == self.pout_start < self.sys.curr and np.abs(end_time - self.end_time) < SMALL_INTERVAL)):
						import pdb; pdb.set_trace()
					pass
				else:
					if not (self.end_time >= end_time - 70 * SMALL_INTERVAL):
						import pdb; pdb.set_trace()
					self.end_time = end_time
			self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D') )
			return

		if self.prev.end_time < start_t + 100 * SMALL_INTERVAL:
			# if self.prev.status = 1, starting the enter maneuver soon => self can travel up to self.prev.dest_to_stop - CAR_LENGTH
			# if self.prev.status = 2, finishing the enter maneuver soon => prev can be ignored
			# if self.prev.status = 5 or 6, leaving the system soon => prev can be ignored 
			if self.prev.status == 1 and self.prev.end_time > start_t:
				if self.curr_loc >= self.prev.dest_to_stop - CAR_LENGTH:
					self.curr_loc = self.prev.dest_to_stop - CAR_LENGTH
					self.traj.addEnd(  changePoint(self.curr_loc, start_t, 0.0) )
					start_t = self.prev.end_time 
			self.traj.addEnd( changePoint(self.curr_loc, start_t, self.driv) )
			self.end_time = start_t + (self.dest_to_stop - self.curr_loc) / self.driv
			if (end_time is not None) and (self.end_time < end_time):
				if spmatch and self.status == 5 and self.idx == self.sys.curr_vehicle.idx:
					import pdb; pdb.set_trace()
					if not ( np.abs(end_time - self.end_time - (self.sys.curr - self.sys.first_service) % meanDRIV) < SMALL_INTERVAL or (self.serv_end == self.pout_start < self.sys.curr and np.abs(end_time - self.end_time) < SMALL_INTERVAL)):
						import pdb; pdb.set_trace()
					pass
				else:
					if not (self.end_time >= end_time - 100 * SMALL_INTERVAL):
						import pdb; pdb.set_trace()
					self.end_time = end_time
			self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D') )
			return

		if self.prev.status == 2 and self.dest_to_stop >= self.prev.curr_loc:

			if angle == 0 and mode == 'long':
				assert self.prev.curr_loc == (self.prev.stop - 1) * LOT_LENGTH

				if not self.after_plin:
					if not self.curr_loc <= self.prev.curr_loc - CAR_LENGTH + 1e-04:
						if self.sys.curr_typ == 'start_pulling_in' and self.sys.curr_vehicle == self.prev:
							print (self.prev.curr_loc - self.curr_loc - CAR_LENGTH)
							import pdb; pdb.set_trace()
							self.curr_loc = self.prev.curr_loc - CAR_LENGTH
						else:
							import pdb; pdb.set_trace()
					self.after_plin = True

					if start_t + (self.prev.curr_loc - self.curr_loc) / self.driv >= self.prev.plin_end:
						self.traj.addEnd( changePoint(self.curr_loc, start_t, self.driv) )
						self.end_time = start_t + (self.dest_to_stop - self.curr_loc) / self.driv
						if end_time is not None and self.end_time < end_time:
							if not self.end_time >= end_time - 62 * SMALL_INTERVAL:
								import pdb; pdb.set_trace()
							self.end_time = end_time
						self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D') )
					else:
						if self.curr_loc < self.prev.curr_loc - CAR_LENGTH:
							self.traj.addEnd( changePoint(self.curr_loc, start_t, self.driv) )
							start_t += (self.prev.curr_loc - CAR_LENGTH - self.curr_loc) / self.driv
						new_speed = CAR_LENGTH / (self.prev.plin_end - start_t)
						self.traj.addEnd( changePoint(self.prev.curr_loc - CAR_LENGTH, start_t, new_speed) )
						self.traj.addEnd( changePoint(self.prev.curr_loc, self.prev.plin_end, self.driv))
						self.end_time = self.prev.plin_end + (self.dest_to_stop - self.prev.curr_loc) / self.driv
						if end_time is not None and self.end_time < end_time:
							if not self.end_time >= end_time - SMALL_INTERVAL:
								import pdb; pdb.set_trace()
							self.end_time = end_time
						self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D') )

				else:
					assert traj is not None
					assert end_time >= self.prev.plin_end
					cp = traj.head
					while cp.nex is not None and cp.nex.data.t <= self.sys.curr:
						cp = cp.nex
					while cp is not None:
						self.traj.addEnd( cp.data )
						cp = cp.nex
					assert self.end_time is not None
				return					

			else:
				if self.status != 5 and np.abs( self.curr_loc + CAR_LENGTH - self.prev.curr_loc ) < SMALL_INTERVAL:
					self.curr_loc = self.prev.curr_loc - CAR_LENGTH

		#####################################################################################
		cp = self.prev.traj.head
		if cp.nex is None:
			import pdb; pdb.set_trace()

		# the assertion above should be fine since cp.nex == None iff either one of the following holds:
		# i) the vehicle is pulling in and thus self.getin is True
		# ii) the vehicle is the head onlane thus self.prev is None
		# it implies that cp.data.v != 'D'
		assert (cp.data.t <= start_t)
		while (cp.nex is not None) and (cp.nex.data.t <= start_t):
			assert cp.data.v != 'D'
			cp = cp.nex

		# the while-loop cannot end with (cp.nex is None)
		# remember that self.prev should already has an accurate trajectory w. a valid end_time
		# if (cp.nex is None), then (self.prev.end_time == cp.nex.data.t)
		# since (self.prev.end_time <= start_t) is already accounted for
		# it is only possible to have (cp.nex.data.t > start_t) and also (cp.nex.data.v = 'D')
		if (cp.nex is None):
			print ('line 575: testing if an option is possible')
			import pdb; pdb.set_trace()
			assert (self.prev.end_time == cp.nex.data.t > start_t)
			assert cp.data.v == 'D'
			# ??? assert self.prev.status == 5
			self.traj.addEnd( changePoint(self.curr_loc, start_t, self.driv) )
			self.end_time = start_t + (self.dest_to_stop - self.curr_loc) / self.driv
			if end_time is not None:
				assert self.end_time >= end_time
			self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D') )
			return
		
		# some sanity check before proceeding
		assert (cp.data.v != 'D') and (cp.nex != None)
		assert (cp.data.t <= start_t < cp.nex.data.t)

		# let us define some counters to keep track of time and location (in the projected traj)
		iter_t = start_t
		if (start_t == self.sys.curr):
			self.prev.update_loc()
			prev_loc = self.prev.curr_loc
		else:
			assert self.status == 5
			prev_loc = cp.data.x + (start_t - cp.data.t) * cp.data.v
			if angle == 90 and 0 <= self.curr_loc + CAR_LENGTH - prev_loc < 1e-12:
				prev_loc = self.curr_loc + CAR_LENGTH
				if cp.data.x < cp.nex.data.x == prev_loc and cp.data.v == 0.0:
					curr = self.prev
					while curr.prev is not None and curr.prev.status != 2:
						curr = curr.prev
					while curr.idx != self.idx:
						curr.update_traj()
						curr = curr.nex
					# import pdb; pdb.set_trace()
					cp = self.prev.traj.head
					assert cp.nex is not None
					assert (cp.data.t <= start_t)
					while (cp.nex is not None) and (cp.nex.data.t <= start_t):
						assert cp.data.v != 'D'
						cp = cp.nex
					assert cp.nex is not None
					assert (cp.data.v != 'D') and (cp.nex != None)
					assert (cp.data.t <= start_t < cp.nex.data.t)	
					if cp.data.x + (start_t - cp.data.t) * cp.data.v != prev_loc:
						import pdb; pdb.set_trace()				

		# the RHS below is the location of prev at start_t
		# check that the distance between the prev and the curr at start_t is at least CAR_LENGTH
		if (self.status != 5) and 0 <= self.curr_loc + CAR_LENGTH - prev_loc < 15e-5:
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
			if self.end_time is None:
				assert (self.stop == 1) or (spmatch and self.stop == 2)
				self.end_time = self.sys.curr
			if end_time is not None:
				assert self.end_time >= end_time
			if not self.end_time >= self.sys.curr:
				if not np.abs(self.end_time - self.sys.curr) < SMALL_INTERVAL:
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
						if spmatch and self.status == 5 and self.idx == self.sys.curr_vehicle.idx:
							import pdb; pdb.set_trace()
							if not ( np.abs(end_time - self.end_time - (self.sys.curr - self.sys.first_service) % meanDRIV) < SMALL_INTERVAL or (self.serv_end == self.pout_start < self.sys.curr and np.abs(end_time - self.end_time) < SMALL_INTERVAL)):
								import pdb; pdb.set_trace()
							pass
						else:
							if not (self.end_time >= end_time - 100 * SMALL_INTERVAL):
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
								if not (self.end_time >= end_time - 21 * SMALL_INTERVAL):
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
										if not (self.end_time >= end_time - SMALL_INTERVAL):
											import pdb; pdb.set_trace()
										self.end_time = end_time
									self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D'))
								else:
									import pdb; pdb.set_trace()
							else:
								self.end_time = iter_t + (self.dest_to_stop - iter_x) / cp.data.v
								if (end_time is not None) and (self.end_time < end_time):
									if not (self.end_time >= end_time - 92 * SMALL_INTERVAL):
										import pdb; pdb.set_trace()
									self.end_time = end_time
								self.traj.addEnd( changePoint(self.dest_to_stop, self.end_time, 'D'))
							return

						if np.abs(iter_x + (cp.nex.data.t - iter_t) * cp.data.v + CAR_LENGTH - cp.nex.data.x) > 1.4e-04:
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
						if spmatch and self.status == 5 and self.idx == self.sys.curr_vehicle.idx:
							import pdb; pdb.set_trace()
							if not ( np.abs(end_time - self.end_time - (self.sys.curr - self.sys.first_service) % meanDRIV) < SMALL_INTERVAL or (self.serv_end == self.pout_start < self.sys.curr and np.abs(end_time - self.end_time) < SMALL_INTERVAL)):
								import pdb; pdb.set_trace()
							pass
						else:
							if not (self.end_time >= end_time - 99 * SMALL_INTERVAL):
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
							if np.abs(iter_x + CAR_LENGTH - prev_loc) < 15e-05:
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
					if spmatch and self.status == 5 and self.idx == self.sys.curr_vehicle.idx:
						import pdb; pdb.set_trace()
						if not ( np.abs(end_time - self.end_time - (self.sys.curr - self.sys.first_service) % meanDRIV) < SMALL_INTERVAL or (self.serv_end == self.pout_start < self.sys.curr and np.abs(end_time - self.end_time) < SMALL_INTERVAL)):
							import pdb; pdb.set_trace()
						pass
					else:
						if not (self.end_time >= end_time - 100 * SMALL_INTERVAL):
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
