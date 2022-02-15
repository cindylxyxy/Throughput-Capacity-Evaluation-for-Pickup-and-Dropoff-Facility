########################################################################################################################################################
# Last checked: Feb 06, 2022 by Xinyu Liu
# This file defines 2 objects: the system object which stores and updates the system state in a trajectory-based simulation,
# the vehicle object which stores and updates the information related to each vehicle in the system;
# and several helper functions.
# This file deals specifically with 0-degree confogurations with short spots.

########################################################################################################################################################
import sys
import json
from math import ceil, floor, sqrt
import numpy as np
from heapq import *
from utils import *
from params import *
from event import *


class systemShort(system):

	def __init__(self, N, seedSERV = None, seedPOUT = None, seedPLIN = None, seedDRIV = None):
		super(systemShort, self).__init__(N, seedSERV, seedPOUT, seedPLIN, seedDRIV)
		self.wait_out = [{'front_drive': 0.0, 'front_in': 0.0, 'other_drive': 0.0,
						  'back_out': 0.0, 'nearby_out': 0.0, 'nearby_in': 0.0, 'oppo_in': 0.0, 
						  'spm_back_out': 0.0, 'total': 0.0, 'veh_count': 0} for _ in range(self.N)]

	def start_pulling_in(self, curr_vehicle):

		curr_vehicle.update_loc()

		if np.abs( curr_vehicle.curr_loc - curr_vehicle.dest_to_stop ) > 1e-5:
			if curr_vehicle.end_time <= self.curr + SMALL_INTERVAL:
				import pdb; pdb.set_trace() 
			self.add_event( event(curr_vehicle.end_time, curr_vehicle, 'start_pulling_in') )
			return

		if self.eventheap != []:
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

		################################################# 0-degree Short Specific #################################################
		# if the vehicle has arrived at the desired destination while its enter maneuver blocked by vehicles on the through lane
		# i.e. 1) if its immediate proceeding vehicle is having an exit maneuver at the spot in front
		if (curr_vehicle.prev is not None) and (curr_vehicle.prev.status == 5):
			assert (curr_vehicle.prev.stop >= curr_vehicle.stop + 1)
			if curr_vehicle.prev.stop == curr_vehicle.stop + 1:
				assert curr_vehicle.prev.pout_end > self.curr
				delayed = True
				req_time = max( req_time, curr_vehicle.prev.pout_end )

		# or 2) if its immediate proceeding vehicle is having an enter maneuver into the spot in front
		elif (curr_vehicle.prev is not None) and (curr_vehicle.prev.status in [1, 2]) and curr_vehicle.prev.stop == curr_vehicle.stop + 1:
			if curr_vehicle.prev.status == 2:
				assert curr_vehicle.prev.plin_end > self.curr
				delayed = True
				req_time = max( req_time, curr_vehicle.prev.plin_end )
			else:
				assert curr_vehicle.prev.status == 1
				assert curr_vehicle.prev.end_time >= self.curr
				delayed = True
				req_time = max( req_time, curr_vehicle.prev.end_time )				

		# or 3) if the immediate proceeding vehicle is stopped.
		elif (curr_vehicle.prev is not None) and (curr_vehicle.prev.curr_loc <= curr_vehicle.dest_to_stop + 1.5 * CAR_LENGTH):
			assert curr_vehicle.status in [1, 6]
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
		###########################################################################################################################
		
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

	def prepare_pulling_out(self, curr_vehicle):

		stop = curr_vehicle.stop
		assert self.inservice[curr_vehicle.stop_idx - 1] == curr_vehicle

		if curr_vehicle.status < 4:
			first_attempt = True
			if curr_vehicle.status == 3:
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
			delayed, req_time, prev, delay_reason, delay_status, delay_speed = self.check_lane_zero_short(curr_vehicle, first_attempt, delay_reason, delay_status, delay_speed)
		else:
			delayed, req_time, prev = self.check_lane_zero_short(curr_vehicle, first_attempt)

		if delayed: 
			if not req_time > self.curr:
				import pdb; pdb.set_trace()
			curr_vehicle.status = 4
			curr_vehicle.end_time = req_time
			self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'prepare_pulling_out') )	
			self.wait_out[curr_vehicle.stop_idx - 1]['total'] += (req_time - self.curr)
			return

		if curr_vehicle.status == 4 and self.eventheap != []:
			broken = False
			heap_len = len(self.eventheap)
			event_holder = []
			while self.eventheap != []:
				next_event = heappop(self.eventheap)
				heappush(event_holder, next_event)
				if next_event.time > self.curr + 1e-05:
					break 
				if next_event.typ == 'prepare_pulling_out' and curr_vehicle.stop_idx < next_event.vehicle.stop_idx:
					if next_event.time <= self.curr:
						import pdb; pdb.set_trace()
					curr_vehicle.end_time = next_event.time
					self.add_event( event(curr_vehicle.end_time, curr_vehicle, 'prepare_pulling_out') )
					self.wait_out[curr_vehicle.stop_idx - 1]['total'] += (curr_vehicle.end_time - self.curr)
					broken = True
					break
				elif spmatch and next_event.typ == 'enter_system' and stop <= g_out + 1: 
					assert self.waiting != []
					if not free_curb:
						next_spot = - min(self.waiting)
					else:
						next_spot = min(self.waiting)
					if idx2spot(next_spot) >= curr_vehicle.stop - 1:
						if (self.inservice[next_spot - 1] is not None) and (self.inservice[next_spot - 1].pout_start > self.curr - 1e-03) and (self.inservice[next_spot - 1].pout_start != self.inservice[next_spot - 1].serv_end):
							import pdb; pdb.set_trace()
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

		if spmatch and self.eventheap != []:
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
				if curr_vehicle.status > 3 and next_event.time > self.curr + meanDRIV - SMALL_INTERVAL:
					break
				if next_event.typ == 'prepare_pulling_out' and next_event.vehicle.stop == stop + 1 and next_event.vehicle.status == 3 and (not self.check_lane_zero_short(next_event.vehicle, True, curr_time = next_event.time)[0]):
					if curr_vehicle.status == 3:
						curr_vehicle.status == 3.5
					curr_vehicle.end_time = next_event.time
					self.add_event( event(curr_vehicle.end_time, curr_vehicle, 'prepare_pulling_out') )
					broken = True
					break
			for event_visited in event_holder:
				self.add_event( event_visited )	
			if broken:
				assert heap_len + 1 == len(self.eventheap)
				return
			assert heap_len == len(self.eventheap)

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
				if curr_vehicle.status > 3 and next_event.time > self.curr + meanDRIV - SMALL_INTERVAL:
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
			for event_visited in event_holder:
				self.add_event( event_visited )	
			if broken:
				assert heap_len + 1 == len(self.eventheap)
				return
			assert heap_len == len(self.eventheap)

		if spmatch and first_attempt:
			est_pout_start = curr_vehicle.serv_end
		elif spmatch and ( (self.curr - self.first_service) % meanDRIV > 1e-03 and meanDRIV - (self.curr - self.first_service) % meanDRIV > 1e-03):
			est_pout_start = self.curr - (self.curr - self.first_service) % meanDRIV
		else:
			est_pout_start = self.curr

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
		elif spmatch and est_pout_start != self.curr:
			assert simType == 'cav'
			assert curr_vehicle.status == 4
			curr_vehicle.status = 5
			curr_vehicle.curr_loc = curr_vehicle.stop * LOT_LENGTH
			curr_vehicle.dest_to_stop = self.n + CAR_LENGTH
			curr_vehicle.pout_start = est_pout_start
			curr_vehicle.pout_time = self.timePOUT.next()
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
				enter_time = self.check_enter_zero_short(enter_time, curr_vehicle.stop)
			self.add_event( event(enter_time - 100 * SMALL_INTERVAL, curr_vehicle, 'add_stop_idx') )
			self.add_event( event(enter_time, new_vehicle, 'enter_system') )
		else:
			self.add_event( event( self.curr, new_vehicle, 'enter_system') )
			if not free_curb:
				heappush( self.waiting, (- curr_vehicle.stop_idx) )
			else:
				heappush( self.waiting, curr_vehicle.stop_idx )
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
				curr_vehicle.curr_loc = CAR_LENGTH
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
			else:
				assert car_time > self.curr
				req_time = max( req_time, car_time )

		if spmatch and req_time <= self.curr:
			if car.curr_loc < 2 * CAR_LENGTH - SMALL_INTERVAL:
				if car.status not in [1, 6]:
					print ('!!!!!!!!!!!!!!!!')
					import pdb; pdb.set_trace()
				assert car.dest_to_stop >= 2 * CAR_LENGTH
				car_time = car.calc_time(2 * CAR_LENGTH)
				req_time = max( req_time, car_time )

		if req_time > self.curr:
			self.entry_blocked = self.curr
			self.entry_cleared = req_time
			self.add_event( event( self.entry_cleared, curr_vehicle, 'enter_system') )
			return

		if debug_idx is not None and curr_vehicle.idx == debug_idx:
			import pdb; pdb.set_trace()

		if control == 'partial':
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
				assert curr_vehicle.prev.curr_loc >= 2 * CAR_LENGTH - SMALL_INTERVAL
				curr_vehicle.curr_loc = LOT_LENGTH
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

			if (J < car.stop - 1) or (car.status == 6):
				j_new.append(j)

			elif car.status == 2:
				# K_in with K = car.j and J = idx2spot(j)
				assert car.stop_idx != j
				assert car.plin_start <= self.curr
				assert car.stop >= 2
				assert (car.stop - 1) * LOT_LENGTH >= curr_vehicle.curr_loc
				if ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV >= max(0.0, meanPLIN - (self.curr - car.plin_start)):
					j_new.append(j)

			elif car.status == 5:
				# K_out with K = car.j and J = idx2spot(j)
				assert car.pout_start <= self.curr <= car.pout_end
				assert (car.stop - 1) * LOT_LENGTH >= curr_vehicle.curr_loc
				if spmatch and ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV >= car.pout_end - self.curr - 1e-05:
					j_new.append(j)
				elif spmatch and ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV >= car.pout_end - self.curr - 200 * SMALL_INTERVAL:
					j_new.append(j)
				elif (not spmatch) and ( (car.stop - 1) * LOT_LENGTH / rateDRIV >= max(0.0, meanPOUT - (self.curr - car.pout_start) - 7e-06) ):
					j_new.append(j)

			else:
				assert car.status == 1
				# I_in with K = car.j and J = idx2spot(j)
				assert car.stop_idx != j
				assert car.stop >= 2
				if not car.curr_loc >= LOT_LENGTH + curr_vehicle.curr_loc - SMALL_INTERVAL:
					import pdb; pdb.set_trace()
				if (car.curr_loc - LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV >= meanPLIN:
					j_new.append(j)

		if j_new != []:
			assert j_new[-1] == max(j_new)
			car_time = 0.0

		else:
			j = - sorted(self.waiting, reverse = True)[0]
			if side == 'double':
				assert (idx2spot(j) >= car.stop - 1)
			else:
				assert side == 'single'
				assert (j >= car.stop - 1)

			if car.status == 2:
				assert car.stop_idx != j
				assert car.plin_start <= self.curr
				assert ( ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV < meanPLIN - (self.curr - car.plin_start) )
				car_time = meanPLIN - (self.curr - car.plin_start) - ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV 

			elif car.status == 5:
				assert car.pout_start <= self.curr <= car.pout_end
				if spmatch:
					assert ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV < car.pout_end - self.curr
					car_time = car.pout_end - self.curr - ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV
				else:
					assert (car.stop - 1) * LOT_LENGTH / rateDRIV < meanPOUT - (self.curr - car.pout_start)
					car_time = meanPOUT - (self.curr - car.pout_start) - (car.stop - 1) * LOT_LENGTH / rateDRIV

			else:
				assert car.status == 1
				assert car.stop_idx != j
				assert (car.curr_loc - LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV < meanPLIN
				car_time = meanPLIN - (car.curr_loc - LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV

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

				if (J < car.stop - 1) or (car.status == 6):
					pass

				elif car.status == 2:
					# K_in with K = car.j and J = idx2spot(j)
					assert car.stop_idx != j
					assert car.plin_start <= self.curr
					assert (car.stop - 1) * LOT_LENGTH >= curr_vehicle.curr_loc
					if ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV < max(0.0, meanPLIN - (self.curr - car.plin_start)):
						j_new = j_new[:idx]
						if idx == 0:
							car_time = max(car_time, meanPLIN - (self.curr - car.plin_start) - ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV)
						break

				elif car.status == 5:
					# K_out with K = car.j and J = idx2spot(j)
					assert car.pout_start <= self.curr <= car.pout_end
					assert (car.stop - 1) * LOT_LENGTH >= curr_vehicle.curr_loc
					if spmatch:
						if ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV < car.pout_end - self.curr - 3e-06:
							if np.abs(((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV - (car.pout_end - self.curr)) > 100 * SMALL_INTERVAL:
								import pdb; pdb.set_trace()
							else:
								j_new = j_new[:idx]
								if idx == 0:
									car_time = max(car_time, car.pout_end - self.curr - ((car.stop - 1) * LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV)
								break
					else:
						if ((car.stop - 1) * LOT_LENGTH / rateDRIV < max(0.0, meanPOUT - (self.curr - car.pout_start) - 3e-06)):
							j_new = j_new[:idx]
							if idx == 0:
								car_time = max(car_time, meanPOUT - (self.curr - car.pout_start) - (car.stop - 1) * LOT_LENGTH / rateDRIV)
							break
	
				else:
					assert car.status == 1
					# I_in with K = car.j and J = idx2spot(j)
					assert car.stop_idx != j
					assert car.curr_loc >= LOT_LENGTH + curr_vehicle.curr_loc
					if (car.curr_loc - LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV < meanPLIN:
						j_new = j_new[:idx]
						if (idx == 0):
							car_time = max(car_time, meanPLIN - (car.curr_loc - LOT_LENGTH - curr_vehicle.curr_loc) / rateDRIV)
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
			assert curr_vehicle.prev.curr_loc >= 2 * CAR_LENGTH - SMALL_INTERVAL
			curr_vehicle.curr_loc = CAR_LENGTH
		assert curr_vehicle.dest_to_stop >= curr_vehicle.curr_loc
		curr_vehicle.update_traj()
		assert curr_vehicle.end_time != None and curr_vehicle.end_time >= self.curr
		self.add_event( event( curr_vehicle.end_time, curr_vehicle, 'start_pulling_in') )
		return

	def check_lane_zero_short(self, curr_vehicle, first_attempt, delay_reason = None, delay_status = None, delay_speed = None, curr_time = None):

		assert angle == 0 and mode == 'short'

		delayed = False
		if curr_time is None:
			curr_time = self.curr
			if spmatch and first_attempt:
				est_pout_start = curr_vehicle.serv_end
				est_pout_end = est_pout_start + meanPOUT - (est_pout_start - self.first_service) % meanDRIV
			elif spmatch and ( (curr_time - self.first_service) % meanDRIV > 1e-03 and meanDRIV - (curr_time - self.first_service) % meanDRIV > 1e-03):
				est_pout_start = curr_time - (curr_time - self.first_service) % meanDRIV
				est_pout_end = est_pout_start + meanPOUT
			else:
				est_pout_start = curr_time
				est_pout_end = est_pout_start + meanPOUT
		else:
			assert curr_time > self.curr
			est_pout_start = curr_time
			est_pout_end = est_pout_start + meanPOUT
		req_time = curr_time

		car = self.head
		prev = None
		stopped = False

		stop = curr_vehicle.stop
		idx = curr_vehicle.idx

		temp_delay = {'front_drive': 0.0, 'front_in': 0.0, 'other_drive': 0.0,
						  'back_out': 0.0, 'nearby_out': 0.0, 'nearby_in': 0.0, 'oppo_in': 0.0, 
						  'spm_back_out': 0.0, 'total': 0.0, 'veh_count': 0} 
		
		while car != None:

			car.update_loc()
			if car.idx == idx:
				pass

			elif car.stop > stop and car.status == 2 and curr_time <= car.plin_end:
				assert car.plin_start <= curr_time <= car.plin_end
				
				if car.stop == stop + 1 and curr_time < car.plin_end - meanPOUT:
					delayed = True
					req_time = max(req_time, car.plin_end - meanPOUT)
					temp_delay['front_in'] = max(temp_delay['front_in'], req_time - self.curr)

				elif (control == 'full' or ptype == 0) and car.stop > stop + 1 and est_pout_end + (car.stop - 1 - stop) * LOT_LENGTH / rateDRIV < max(curr_time, car.plin_start + meanPLIN) - SMALL_INTERVAL:
					assert car.plin_start + meanPLIN > curr_time
					assert car.plin_start + meanPLIN > est_pout_end + (car.stop - 1 - stop) * LOT_LENGTH / rateDRIV
					# assert car.plin_start + meanPLIN > est_pout_start + meanPOUT + (car.stop - 1 - stop) * LOT_LENGTH / rateDRIV + SMALL_INTERVAL
					if meanPLIN + car.plin_start - (meanPOUT + (car.stop - 1 - stop) * LOT_LENGTH / rateDRIV) > req_time:
						delayed = True
						car_time = meanPLIN + car.plin_start - (meanPOUT + (car.stop - 1 - stop) * LOT_LENGTH / rateDRIV)
						req_time = max( req_time, car_time )				
						temp_delay['front_in'] = max(temp_delay['front_in'], car_time - self.curr)
					else:
						prev = car
				else:
					prev = car

			elif car.status == 5 and car.stop == stop + 1:
				if spmatch and car.pout_start == car.serv_end:
					if curr_time == car.pout_start:
						delayed = True
						car_time = car.pout_start + meanDRIV - (car.pout_start - self.first_service) % meanDRIV
					elif first_attempt and est_pout_start - (est_pout_start - self.first_service) % meanDRIV == car.pout_start - (car.pout_start - self.first_service) % meanDRIV:
						delayed = True
						car_time = car.pout_start + meanDRIV - (car.pout_start - self.first_service) % meanDRIV
					elif (not first_attempt) and np.abs( est_pout_start - (car.pout_start - (car.pout_start - self.first_service) % meanDRIV)) <= 100 * SMALL_INTERVAL:
						delayed = True
						car_time = car.pout_start + meanDRIV - (car.pout_start - self.first_service) % meanDRIV
				else:
					if curr_time < meanDRIV + car.pout_start:
						delayed = True
						car_time = min(car.pout_start + meanDRIV, car.pout_end)

				if delayed:
					try:
						if car_time <= curr_time:
							import pdb; pdb.set_trace()
						req_time = max( req_time, car_time )
						temp_delay['nearby_out'] = max(temp_delay['nearby_out'], car_time - self.curr)
					except UnboundLocalError:
						pass
				else:
					prev = car

			elif car.curr_loc >= stop * LOT_LENGTH + CAR_LENGTH - SMALL_INTERVAL:

				if car.stop == stop + 1 and car.status == 1:
					assert car.end_time >= self.curr
					if car.end_time == self.curr:
						delayed = True
						req_time = max( req_time, car.end_time + SMALL_INTERVAL )
						temp_delay['front_drive'] = max(temp_delay['front_drive'], car.end_time + SMALL_INTERVAL - self.curr)
					else:
						delayed = True
						req_time = max( req_time, car.end_time )
						temp_delay['front_drive'] = max(temp_delay['front_drive'], car.end_time - self.curr)

				elif (control == 'full' or ptype == 0) and car.stop > stop and car.status == 1:
					car_time = est_pout_end + (car.stop - stop - 1) * LOT_LENGTH / rateDRIV 
					if car_time < car.end_time + meanPLIN - 2e-05:
						delayed = True
						car_time = car.end_time + meanPLIN - meanPOUT - (car.stop - stop - 1) * LOT_LENGTH / rateDRIV
						req_time = max( req_time, car_time)
						temp_delay['other_drive'] = max(temp_delay['other_drive'], car_time - self.curr)
					else:
						prev = car

				elif control == 'partial' and car.curr_loc < (stop + 1) * LOT_LENGTH + CAR_LENGTH and curr_time < car.calc_time( (stop + 1) * LOT_LENGTH + CAR_LENGTH ) - 1e-05:
					cp = car.traj.head
					if cp.data.t > curr_time:
						if cp.data.t <= curr_time + 2 * SMALL_INTERVAL:
							delayed = True
							req_time = max( req_time, cp.data.t )
							temp_delay['other_drive'] = max(temp_delay['other_drive'], cp.data.t - self.curr)
							import pdb; pdb.set_trace()
						else:
							import pdb; pdb.set_trace()
					else:
						while (cp.nex is not None) and (cp.nex.data.t <= curr_time):
							cp = cp.nex
						if cp.nex is not None and cp.nex.data.t <= curr_time + SMALL_INTERVAL:
							cp = cp.nex
						assert cp is not None
						if (cp.data.v == 0.0):
							assert cp.nex.data.t > curr_time
							delayed = True
							req_time = max( req_time, cp.nex.data.t )
							temp_delay['other_drive'] = max(temp_delay['other_drive'], cp.nex.data.t - self.curr)
						else:
							prev = car			
				else:
					prev = car

			elif car.status == 5:
				stopped = False
				assert car.stop < stop or (car.stop == stop and car.stop_idx != curr_vehicle.stop_idx)
				assert car.dest_to_stop >= stop * LOT_LENGTH + CAR_LENGTH
				assert car.pout_start <= self.curr
				if car.stop_idx in out_range(curr_vehicle.stop_idx, self.N) and car.pout_end - curr_time > SMALL_INTERVAL:
					assert car.pout_end > curr_time
					delayed = True
					req_time = max(req_time, car.pout_end)
					temp_delay['nearby_out'] = max(temp_delay['nearby_out'], car.pout_end - self.curr)
					break

				assert car.stop < stop - 1 or car.pout_end - curr_time <= SMALL_INTERVAL
				assert (stop - 1) * LOT_LENGTH >= car.curr_loc or car.stop == stop
				if spmatch:				
					car_time = car.pout_end + ((stop - 1) * LOT_LENGTH - car.curr_loc) / rateDRIV
					if (first_attempt and car_time < est_pout_start + meanPOUT - (est_pout_start - self.first_service) % meanDRIV - 100 * SMALL_INTERVAL) or ((not first_attempt) and car_time < est_pout_start + meanPOUT - 100 * SMALL_INTERVAL):
						car_time = car.calc_time( stop * LOT_LENGTH + CAR_LENGTH )
						assert car_time > curr_time
						delayed = True
						req_time = max( req_time, car_time )
						temp_delay['spm_back_out'] = max(temp_delay['spm_back_out'], car_time - self.curr)
				else:
					car_time = max(self.curr, meanPOUT + car.pout_start) + ((stop - 1) * LOT_LENGTH - car.curr_loc) / rateDRIV
					if car_time < curr_time + meanPOUT:
						car_time = car.calc_time( stop * LOT_LENGTH + CAR_LENGTH )
						assert car_time > self.curr
						delayed = True
						req_time = max( req_time, car_time )
						temp_delay['back_out'] = max(temp_delay['back_out'], car_time - self.curr)

			elif car.status == 2:
				assert car.plin_start <= self.curr <= car.plin_end
				if car.stop == stop and est_pout_start < car.plin_end - 10 * SMALL_INTERVAL:
					assert (side == 'double' and car.stop_idx != curr_vehicle.stop_idx)
					stopped = False
					delayed = True
					req_time = max( req_time, car.plin_end + SMALL_INTERVAL )
					temp_delay['oppo_in'] = max(temp_delay['oppo_in'], car.plin_end + SMALL_INTERVAL - self.curr)
				elif car.stop >= stop - dgap and est_pout_start + meanDRIV < car.plin_end - 10 * SMALL_INTERVAL and car.plin_end - meanDRIV > req_time:
					stopped = False
					delayed = True
					req_time = max( req_time, car.plin_end - meanDRIV )
					temp_delay['nearby_in'] = max(temp_delay['nearby_in'], car.plin_end - meanDRIV - self.curr)
				elif car.stop < stop - dgap and car.plin_end > curr_time:
					stopped = True

			elif car.status == 1 and car.stop_idx == curr_vehicle.stop_idx:
				assert car.stop == stop
				assert (car.curr_loc < car.dest_to_stop) and (car.curr_loc < (stop - 1) * LOT_LENGTH + SMALL_INTERVAL)
				stopped = True

			elif car.status == 1 and car.stop < stop - 1:
				assert (car.curr_loc <= car.dest_to_stop <= (stop - 1) * LOT_LENGTH)
				stopped = True

			else:
				assert (car.status == 6) or (car.status == 1 and car.stop >= stop - 1 and car.stop_idx != curr_vehicle.stop_idx)
				assert car.dest_to_stop >= stop * LOT_LENGTH + CAR_LENGTH or car.stop == stop or car.stop == stop - 1
				car_time = self.curr + ((stop - 1) * LOT_LENGTH - car.curr_loc) / rateDRIV
				if stopped and car.prev.status == 1 and car.prev.stop_idx == curr_vehicle.stop_idx:
					if car.curr_loc >= (stop - 2) * LOT_LENGTH + SMALL_INTERVAL: 
						import pdb; pdb.set_trace()
					break
				elif stopped and car.prev.status == 1:
					assert car.prev.end_time is not None
					assert car.prev.stop < stop
					if not (car.prev.stop - 1) * LOT_LENGTH >= car.curr_loc - 1.5e-04:
						import pdb; pdb.set_trace()
					car_time = max(meanPLIN + car.prev.end_time, self.curr + max(0.0, ((car.prev.stop - 1) * LOT_LENGTH - car.curr_loc)) / rateDRIV) + (stop - car.prev.stop) * LOT_LENGTH / rateDRIV
				elif stopped:
					assert car.prev.status == 2
					assert car.prev.end_time is not None
					assert car.prev.stop < stop
					if not (car.prev.stop - 1) * LOT_LENGTH >= car.curr_loc - 1.5e-04:
						import pdb; pdb.set_trace()
					car_time = max(car.prev.plin_end, self.curr + max(0.0, ((car.prev.stop - 1) * LOT_LENGTH - car.curr_loc)) / rateDRIV) + (stop - car.prev.stop) * LOT_LENGTH / rateDRIV 

				if ((not spmatch or not first_attempt) and car_time < est_pout_start + meanPOUT - 100 * SMALL_INTERVAL) or (spmatch and first_attempt and car_time < est_pout_start + meanPOUT - (est_pout_start - self.first_service) % meanDRIV - 100 * SMALL_INTERVAL):
					delayed = True
					if car.status == 1 and car.stop <= stop:
						assert car.stop_idx != curr_vehicle.stop_idx
						if car.end_time < self.curr:
							import pdb; pdb.set()
						elif car.end_time == self.curr:
							req_time = max( req_time, car.end_time + SMALL_INTERVAL )
							temp_delay['other_drive'] = max(temp_delay['other_drive'], car.end_time + SMALL_INTERVAL - self.curr)
						else:
							req_time = max( req_time, car.end_time )
							temp_delay['other_drive'] = max(temp_delay['other_drive'], car.end_time - self.curr)
					else:
						assert car.dest_to_stop >= stop * LOT_LENGTH + CAR_LENGTH
						car_time = car.calc_time( stop * LOT_LENGTH + CAR_LENGTH )
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

	def check_enter_zero_short(self, enter_time, stop):

		assert spmatch and control == 'full'
		car = self.head
		while car is not None:
			if (stop < car.stop - 1) or (car.status == 6):
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
				if enter_time < car.end_time + meanPLIN - (car.stop - 1) * LOT_LENGTH / rateDRIV - 10 * SMALL_INTERVAL:
					enter_time = car.end_time + meanPLIN - (car.stop - 1) * LOT_LENGTH / rateDRIV
			car = car.nex
		return enter_time
