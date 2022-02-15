########################################################################################################################################################
# Last checked: Feb 06, 2022 by Xinyu Liu
# This file defines 2 objects: the system object which stores and updates the system state in a trajectory-based simulation,
# the vehicle object which stores and updates the information related to each vehicle in the system;
# and several helper functions.
# This file deals specifically with 90-degree confogurations.

########################################################################################################################################################
import sys
import json
from math import ceil, floor, sqrt
import numpy as np
from heapq import *
from utils import *
from params import *
from event import *


class system90(system):

	def __init__(self, N, seedSERV = None, seedPOUT = None, seedPLIN = None, seedDRIV = None):
		super(system90, self).__init__(N, seedSERV, seedPOUT, seedPLIN, seedDRIV)
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

		if spmatch:

			delayed = False
			req_time = self.curr

			################################################### 90-degree Specific ####################################################
			# if the vehicle has arrived at the desired destination while its enter maneuver blocked by vehicles on the through lane
			# i.e. 1) if its immediate proceeding vehicle is having an exit maneuver at the spot in front
			if curr_vehicle.prev is not None and curr_vehicle.prev.status == 5:
				assert (curr_vehicle.prev.block_idx >= curr_vehicle.block_idx + 1)
				if curr_vehicle.prev.block_idx == curr_vehicle.block_idx + 1:
					assert curr_vehicle.prev.pout_end > self.curr
					delayed = True
					req_time = max( req_time, curr_vehicle.prev.pout_end )

			# or 2) if its immediate proceeding vehicle is having an enter maneuver into the spot in front
			elif curr_vehicle.prev is not None and curr_vehicle.prev.status in [1, 2] and curr_vehicle.prev.block_idx == curr_vehicle.block_idx + 1:
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
	
		if not self.inservice[curr_vehicle.stop_idx - 1] is None:
			import pdb; pdb.set_trace()
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
		block_idx = curr_vehicle.block_idx
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
			delayed, req_time, prev, delay_reason, delay_status, delay_speed = self.check_lane_90(curr_vehicle, first_attempt, delay_reason, delay_status, delay_speed)
		else:
			delayed, req_time, prev = self.check_lane_90(curr_vehicle, first_attempt)

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
				elif spmatch and next_event.typ == 'enter_system' and block_idx <= g_out + 1: 
					assert self.waiting != []
					if not free_curb:
						next_spot = - min(self.waiting)
					else:
						next_spot = min(self.waiting)
					if spot2blk(idx2spot(next_spot)) >= block_idx - 1:
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
				if next_event.typ == 'prepare_pulling_out' and next_event.vehicle.stop <= stop + dgap and next_event.vehicle.stop_idx > curr_vehicle.stop_idx and next_event.vehicle.status == 3 and (not self.check_lane_90(next_event.vehicle, True, curr_time = next_event.time)[0]):
					if curr_vehicle.status == 3:
						curr_vehicle.status = 3.5
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
				if stop == next_event.vehicle.stop and curr_vehicle.stop_idx < next_event.vehicle.stop_idx:
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
			curr_vehicle.curr_loc = curr_vehicle.block_idx * CAR_LENGTH
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
			curr_vehicle.curr_loc = curr_vehicle.block_idx * CAR_LENGTH
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
				if block_idx <= m_out:
					enter_time += meanDRIV * (m_out + 1 - block_idx)
				enter_time = self.check_enter_90(enter_time, curr_vehicle.stop)
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
			if side == 'single':
				affected = [0, 1, 2]
			else:
				assert side == 'double'
				affected = [0, 1, 2, 3, 4, 5]
			for idx in affected:
				if idx < len(self.inservice) and self.inservice[idx] is not None and self.inservice[idx].plin_end is not None and self.inservice[idx].plin_end + meanDRIV > self.curr + SMALL_INTERVAL:
					self.entry_blocked = self.curr
					self.entry_cleared = self.inservice[idx].plin_end + meanDRIV
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
				req_time = max( req_time, car_time)

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
				curr_vehicle.curr_loc = CAR_LENGTH
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

			if (spot2blk(J) < car.block_idx - 1) or (car.status == 6):
				j_new.append(j)

			elif car.status == 2:
				# K_in with k = car.stop_idx and j
				assert car.stop_idx != j
				assert car.plin_start <= self.curr
				if spmatch:
					assert car.block_idx >= 2
					assert (car.block_idx - 1) * CAR_LENGTH >= curr_vehicle.curr_loc
					if ((car.block_idx - 1) * CAR_LENGTH - curr_vehicle.curr_loc) / rateDRIV >= max(0.0, meanPLIN - (self.curr - car.plin_start)):
						j_new.append(j)
				else:
					assert (car.stop - 1) * LOT_LENGTH >= curr_vehicle.curr_loc == 0.0
					if (car.stop - 1) * LOT_LENGTH / rateDRIV >= max(0.0, meanPLIN - (self.curr - car.plin_start)):
						j_new.append(j)

			elif car.status == 5:
				# K_out with k = car.stop_idx and j
				assert car.pout_start <= self.curr <= car.pout_end
				if spmatch:
					assert (car.block_idx - 1) * CAR_LENGTH >= curr_vehicle.curr_loc
					if ((car.block_idx - 1) * CAR_LENGTH - curr_vehicle.curr_loc) / rateDRIV >= car.pout_end - self.curr - 9e-06:
						j_new.append(j)
					elif ((car.block_idx - 1) * CAR_LENGTH - curr_vehicle.curr_loc) / rateDRIV >= car.pout_end - self.curr - 200 * SMALL_INTERVAL:
						j_new.append(j)
				else:
					assert (car.stop - 1) * LOT_LENGTH >= curr_vehicle.curr_loc == 0.0
					if ( (car.stop - 1) * LOT_LENGTH / rateDRIV >= max(0.0, meanPOUT - (self.curr - car.pout_start) - 7e-06) ):
						j_new.append(j)

			else:
				assert car.status == 1
				# I_in with k = car.stop_idx and j
				assert car.stop_idx != j
				assert car.curr_loc >= CAR_LENGTH + curr_vehicle.curr_loc - SMALL_INTERVAL
				if (car.curr_loc - CAR_LENGTH - curr_vehicle.curr_loc) / rateDRIV >= meanPLIN:
					j_new.append(j)

		if j_new != []:
			assert j_new[-1] == max(j_new)
			car_time = 0.0

		else:
			j = - sorted(self.waiting, reverse = True)[0]
			if side == 'double':
				assert (spot2blk(idx2spot(j)) >= car.block_idx - 1)
			else:
				assert side == 'single'
				assert (spot2blk(j) >= car.block_idx - 1)

			if car.status == 2:
				assert car.stop_idx != j
				assert car.plin_start <= self.curr
				if spmatch:
					assert ((car.block_idx - 1) * CAR_LENGTH - curr_vehicle.curr_loc) / rateDRIV < meanPLIN - (self.curr - car.plin_start)
					car_time = meanPLIN - (self.curr - car.plin_start) - ((car.block_idx - 1) * CAR_LENGTH - curr_vehicle.curr_loc) / rateDRIV
				else:
					assert (car.stop - 1) * LOT_LENGTH / rateDRIV < meanPLIN - (self.curr - car.plin_start)
					car_time = meanPLIN - (self.curr - car.plin_start) - (car.stop - 1) * LOT_LENGTH / rateDRIV 

			elif car.status == 5:
				assert car.pout_start <= self.curr <= car.pout_end
				if spmatch:
					assert ((car.block_idx - 1) * CAR_LENGTH - curr_vehicle.curr_loc) / rateDRIV < car.pout_end - self.curr
					car_time = car.pout_end - self.curr - ((car.block_idx - 1) * CAR_LENGTH - curr_vehicle.curr_loc) / rateDRIV
				else:
					assert (car.stop - 1) * LOT_LENGTH / rateDRIV < meanPOUT - (self.curr - car.pout_start)
					car_time = meanPOUT - (self.curr - car.pout_start) - (car.stop - 1) * LOT_LENGTH / rateDRIV

			else:
				assert car.status == 1
				assert car.stop_idx != j
				assert (car.curr_loc - CAR_LENGTH - curr_vehicle.curr_loc) / rateDRIV < meanPLIN
				car_time = meanPLIN - (car.curr_loc - CAR_LENGTH - curr_vehicle.curr_loc) / rateDRIV

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

				if (spot2blk(J) < car.block_idx - 1) or (car.status == 6):
					pass

				elif car.status == 2:
					# K_in with K = car.j and J = idx2spot(j)
					assert car.stop_idx != j
					assert car.plin_start <= self.curr
					if spmatch:
						assert (car.block_idx - 1) * CAR_LENGTH >= curr_vehicle.curr_loc
						if ((car.block_idx - 1) * CAR_LENGTH - curr_vehicle.curr_loc) / rateDRIV < max(0.0, meanPLIN - (self.curr - car.plin_start)):
							j_new = j_new[:idx]
							if idx == 0:
								car_time = max(car_time, meanPLIN - (self.curr - car.plin_start) - ((car.block_idx - 1) * CAR_LENGTH - curr_vehicle.curr_loc) / rateDRIV)
							break
					else:
						assert (car.stop - 1) * LOT_LENGTH >= curr_vehicle.curr_loc == 0.0
						if (car.stop - 1) * LOT_LENGTH / rateDRIV < max(0.0, meanPLIN - (self.curr - car.plin_start)):
							j_new = j_new[:idx]
							if idx == 0:
								car_time = max(car_time, meanPLIN - (self.curr - car.plin_start) - (car.stop - 1) * LOT_LENGTH / rateDRIV)
							break						

				elif car.status == 5:
					# K_out with K = car.j and J = idx2spot(j)
					assert car.pout_start <= self.curr <= car.pout_end
					if spmatch:
						assert (car.block_idx - 1) * CAR_LENGTH >= curr_vehicle.curr_loc
						if ((car.block_idx - 1) * CAR_LENGTH - curr_vehicle.curr_loc) / rateDRIV < car.pout_end - self.curr - 3e-06:
							if ((car.block_idx - 1) * CAR_LENGTH - curr_vehicle.curr_loc) / rateDRIV > car.pout_end - self.curr - 100 * SMALL_INTERVAL:
								import pdb; pdb.set_trace()
							else:
								j_new = j_new[:idx]
								if idx == 0:
									car_time = max(car_time, car.pout_end - self.curr - ((car.block_idx - 1) * CAR_LENGTH - curr_vehicle.curr_loc) / rateDRIV)
								break
					else:
						assert (car.stop - 1) * LOT_LENGTH >= curr_vehicle.curr_loc == 0.0
						if (car.stop - 1) * LOT_LENGTH / rateDRIV < max(0.0, meanPOUT - (self.curr - car.pout_start) - 3e-06):
							j_new = j_new[:idx]
							if idx == 0:
								car_time = max(car_time, meanPOUT - (self.curr - car.pout_start) - (car.stop - 1) * LOT_LENGTH / rateDRIV)
							break

				else:
					assert car.status == 1
					# I_in with K = car.j and J = idx2spot(j)
					assert car.stop_idx != j
					assert car.curr_loc >= CAR_LENGTH + curr_vehicle.curr_loc
					if (car.curr_loc - CAR_LENGTH - curr_vehicle.curr_loc) / rateDRIV < meanPLIN:
						j_new = j_new[:idx]
						if (idx == 0):
							car_time = max(car_time, meanPLIN - (car.curr_loc - CAR_LENGTH - curr_vehicle.curr_loc) / rateDRIV)
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
			assert (j_new[0]) in self.waiting
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

	def check_lane_90(self, curr_vehicle, first_attempt, delay_reason = None, delay_status = None, delay_speed = None, curr_time = None):

		assert angle == 90

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
			if not curr_time > self.curr:
				import pdb; pdb.set_trace()
			est_pout_start = curr_time
			est_pout_end = est_pout_start + meanPOUT
		req_time = curr_time

		car = self.head
		prev = None
		stopped = False
		
		stop = curr_vehicle.stop
		idx = curr_vehicle.idx
		block_idx = curr_vehicle.block_idx

		temp_delay = {'front_drive': 0.0, 'front_in': 0.0, 'other_drive': 0.0,
						  'back_out': 0.0, 'nearby_out': 0.0, 'nearby_in': 0.0, 'oppo_in': 0.0, 
						  'spm_back_out': 0.0, 'total': 0.0, 'veh_count': 0} 

		while car != None:

			car.update_loc()
			if car.idx == idx:
				pass

			elif car.status == 2:
				assert car.plin_start <= self.curr <= car.plin_end
				if car.stop == stop and est_pout_start < car.plin_end - 10 * SMALL_INTERVAL:
					assert (side == 'double' and car.stop_idx != curr_vehicle.stop_idx)
					stopped = False
					delayed = True
					req_time = max( req_time, car.plin_end + SMALL_INTERVAL )
					temp_delay['oppo_in'] = max(temp_delay['oppo_in'], car.plin_end + SMALL_INTERVAL - self.curr)
				elif stop + dgap >= car.stop >= stop - dgap and est_pout_start < car.plin_end - 1e-05:
					stopped = False
					delayed = True
					req_time = max( req_time, car.plin_end + SMALL_INTERVAL )
					temp_delay['nearby_in'] = max(temp_delay['nearby_in'], car.plin_end + SMALL_INTERVAL - self.curr)
				elif car.stop < stop - dgap and car.plin_end > curr_time:
					stopped = True
				elif car.block_idx > block_idx and curr_time <= car.plin_end:
					if car.block_idx == block_idx + 1 and (car.stop > stop + dgap) and curr_time < car.plin_end - meanPOUT:
						delayed = True
						req_time = max( req_time, car.plin_end - meanPOUT )
						temp_delay['front_in'] = max(temp_delay['front_in'], car.plin_end - meanPOUT - self.curr)

					elif control == 'full' or ptype == 0:
						if spmatch and est_pout_end + (car.block_idx - 1 - block_idx) * CAR_LENGTH / rateDRIV < max(curr_time, car.plin_start + meanPLIN) - 2e-05:
							assert car.plin_start + meanPLIN > curr_time
							assert car.plin_start + meanPLIN > est_pout_end + (car.block_idx - 1 - block_idx) * CAR_LENGTH / rateDRIV
							car_time = meanPLIN + car.plin_start - (meanPOUT + (car.block_idx - 1 - block_idx) * CAR_LENGTH / rateDRIV)
							if car_time > req_time:
								delayed = True
								req_time = max( req_time, car_time )
								temp_delay['front_in'] = max(temp_delay['front_in'], car_time - self.curr)
							else:
								prev = car

						elif (not spmatch) and est_pout_end + (car.stop - 1 - stop) * LOT_LENGTH / rateDRIV < max(curr_time, car.plin_start + meanPLIN) - 2e-05:
							assert car.plin_start + meanPLIN > curr_time
							assert car.plin_start + meanPLIN > est_pout_end + (car.stop - 1 - stop) * LOT_LENGTH / rateDRIV
							if not car.plin_start + meanPLIN > est_pout_start + meanPOUT + (car.stop - 1 - stop) * LOT_LENGTH / rateDRIV:
								import pdb; pdb.set_trace()
							else:	
								delayed = True
								car_time = meanPLIN + car.plin_start - (meanPOUT + (car.stop - 1 - stop) * LOT_LENGTH / rateDRIV)
								req_time = max( req_time, car_time )
								temp_delay['front_in'] = max(temp_delay['front_in'], car_time - self.curr)
						else:	
							prev = car
					else:
						prev = car				

			elif car.stop_idx in out_range(curr_vehicle.stop_idx, self.N) and car.status == 5 and car.pout_end - curr_time > SMALL_INTERVAL:
				assert car.pout_end > curr_time
				delayed = True
				req_time = max(req_time, car.pout_end + meanDRIV)
				temp_delay['nearby_out'] = max(temp_delay['nearby_out'], car.pout_end + meanDRIV - self.curr)
				break

			elif car.status == 5 and car.block_idx == block_idx + 1 and spmatch:
				prev = car

			elif spmatch and car.curr_loc >= block_idx * CAR_LENGTH + CAR_LENGTH - SMALL_INTERVAL:

				if car.stop_idx in out_range(curr_vehicle.stop_idx, self.N) and car.status == 6 and car.curr_loc <= (block_idx + 2) * CAR_LENGTH - SMALL_INTERVAL and car.calc_time((block_idx + 2) * CAR_LENGTH) > curr_time:
					delayed = True
					car_time = car.calc_time((block_idx + 2) * CAR_LENGTH)
					req_time = max(req_time, car_time )
					temp_delay['nearby_out'] = max(temp_delay['nearby_out'], car_time - self.curr)
					break

				elif car.block_idx == block_idx + 1 and car.status == 1:
					assert car.end_time >= self.curr
					if car.end_time == self.curr:
						delayed = True
						req_time = max( req_time, car.end_time + SMALL_INTERVAL )
						temp_delay['front_drive'] = max(temp_delay['front_drive'], car.end_time + SMALL_INTERVAL - self.curr)
					else:
						delayed = True
						req_time = max( req_time, car.end_time )
						temp_delay['front_drive'] = max(temp_delay['front_drive'], car.end_time - self.curr)

				elif (control == 'full' or ptype == 0) and car.block_idx > block_idx and car.status == 1:
					if est_pout_end + (car.block_idx - block_idx - 1) * CAR_LENGTH / rateDRIV < car.end_time + meanPLIN - 2e-05:
						car_time = car.end_time + meanPLIN - meanPOUT - (car.block_idx - block_idx - 1) * CAR_LENGTH / rateDRIV 
						if car_time > req_time:
							delayed = True
							req_time = max( req_time, car_time )
							temp_delay['other_drive'] = max(temp_delay['other_drive'], car_time - self.curr)
						else:
							prev = car
					else:
						prev = car

				elif control == 'partial' and car.curr_loc < (block_idx + 1) * CAR_LENGTH + CAR_LENGTH and curr_time < car.calc_time( (block_idx + 1) * CAR_LENGTH + CAR_LENGTH ) - 2e-05:
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

			elif (not spmatch) and car.curr_loc >= (stop + dgap) * LOT_LENGTH + CAR_LENGTH - SMALL_INTERVAL:

				if (control == 'full' or ptype == 0) and car.stop > (stop + dgap) and car.status == 1:
					car_time = est_pout_end + (car.stop - 1 - stop - dgap) * LOT_LENGTH / rateDRIV 
					if car_time < car.end_time + meanPLIN - 2e-05:
						if not car.end_time + meanPLIN > est_pout_start + meanPOUT + (car.stop - 1 - stop - dgap) * LOT_LENGTH / rateDRIV:
							import pdb; pdb.set_trace()
						else:
							delayed = True
							car_time = car.end_time + meanPLIN - meanPOUT - (car.stop - 1 - stop - dgap) * LOT_LENGTH / rateDRIV
							req_time = max( req_time, car_time )
							temp_delay['other_drive'] = max(temp_delay['other_drive'], car_time - self.curr)
					else:
						prev = car
				else:
					prev = car

			elif car.status == 5:
				stopped = False
				assert car.stop < stop - dgap or car.pout_end - curr_time <= SMALL_INTERVAL
				assert car.pout_start <= self.curr	
				if spmatch:	
					assert (block_idx - 1) * CAR_LENGTH >= car.curr_loc	or car.pout_end - curr_time <= SMALL_INTERVAL
					assert car.dest_to_stop >= block_idx * CAR_LENGTH + CAR_LENGTH
					car_time = car.pout_end + ((block_idx - 1) * CAR_LENGTH - car.curr_loc) / rateDRIV
					if (first_attempt and car_time < est_pout_start + meanPOUT - (est_pout_start - self.first_service) % meanDRIV - 100 * SMALL_INTERVAL) or ((not first_attempt) and car_time < est_pout_start + meanPOUT - 100 * SMALL_INTERVAL):
						car_time = car.calc_time( block_idx * CAR_LENGTH + CAR_LENGTH )
						assert car_time > curr_time
						delayed = True
						req_time = max( req_time, car_time )
						temp_delay['spm_back_out'] = max(temp_delay['spm_back_out'], car_time - self.curr)
				else:
					assert (stop - 1) * LOT_LENGTH >= car.curr_loc or car.pout_end - curr_time <= SMALL_INTERVAL		
					if not car.dest_to_stop >= min( (stop + dgap) * LOT_LENGTH + CAR_LENGTH, self.n + CAR_LENGTH):
						import pdb; pdb.set_trace()
					car_time = max(self.curr, meanPOUT + car.pout_start) + ((stop - 1) * LOT_LENGTH - car.curr_loc) / rateDRIV
					if car_time < curr_time + meanPOUT:
						car_time = car.calc_time( (stop + dgap) * LOT_LENGTH + CAR_LENGTH )
						assert car_time > self.curr
						delayed = True
						req_time = max( req_time, car_time )
						temp_delay['back_out'] = max(temp_delay['back_out'], car_time - self.curr)

			elif car.status == 1 and car.stop_idx == curr_vehicle.stop_idx:
				assert car.stop == stop
				assert (car.curr_loc < car.dest_to_stop)
				if spmatch:
					assert car.curr_loc < (block_idx - 1) * CAR_LENGTH + SMALL_INTERVAL
				else:
					assert car.curr_loc < (stop - 1) * LOT_LENGTH + SMALL_INTERVAL
				stopped = True

			elif car.status == 1 and spmatch and car.block_idx < block_idx - 1:
				assert (car.curr_loc <= car.dest_to_stop <= (block_idx - 1) * CAR_LENGTH)
				stopped = True
		
			elif car.status == 1 and (not spmatch) and car.stop < stop - dgap:
				assert (car.curr_loc <= car.dest_to_stop <= (stop - 1) * LOT_LENGTH )
				stopped = True

			elif spmatch:
				assert (car.status == 6) or (car.status == 1 and car.block_idx >= block_idx - 1 and car.stop_idx != curr_vehicle.stop_idx)
				assert car.dest_to_stop >= block_idx * CAR_LENGTH + CAR_LENGTH or car.block_idx == block_idx or car.block_idx == block_idx - 1
				car_time = self.curr + ((block_idx - 1) * CAR_LENGTH - car.curr_loc) / rateDRIV
				if stopped and car.prev.status == 1 and car.prev.stop_idx == curr_vehicle.stop_idx:
					if car.curr_loc >= (block_idx - 2) * CAR_LENGTH + SMALL_INTERVAL: 
						import pdb; pdb.set_trace()
					break
				elif stopped:
					assert car.prev.end_time is not None
					assert car.prev.block_idx < block_idx
					if not (car.prev.block_idx - 1) * CAR_LENGTH >= car.curr_loc - 2e-04:
						import pdb; pdb.set_trace()
					car_time = (block_idx - car.prev.block_idx) * CAR_LENGTH / rateDRIV
					if car.prev.status == 1:		
						car_time += max(meanPLIN + car.prev.end_time, self.curr + max(0.0, ((car.prev.block_idx - 1) * CAR_LENGTH - car.curr_loc)) / rateDRIV)
					else:
						car_time += max(car.prev.plin_end, self.curr + max(0.0, ((car.prev.block_idx - 1) * CAR_LENGTH - car.curr_loc)) / rateDRIV)

				if ((not first_attempt) and car_time < est_pout_start + meanPOUT - 100 * SMALL_INTERVAL) or (first_attempt and car_time < est_pout_start + meanPOUT - (est_pout_start - self.first_service) % meanDRIV - 100 * SMALL_INTERVAL):
					delayed = True
					if car.status == 1 and car.block_idx <= block_idx:
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
						assert car.dest_to_stop >= block_idx * CAR_LENGTH + CAR_LENGTH
						car_time = car.calc_time( block_idx * CAR_LENGTH + CAR_LENGTH )
						req_time = max( req_time, car_time )
						temp_delay['other_drive'] = max(temp_delay['other_drive'], car_time - self.curr)

				if stopped:
					break

			else:
				assert (not spmatch)
				assert (car.status == 6) or (car.status == 1 and car.stop >= stop - dgap and car.stop_idx != curr_vehicle.stop_idx)
				assert car.dest_to_stop >= (stop + dgap) * LOT_LENGTH + CAR_LENGTH or car.stop <= stop + dgap
				car_time = self.curr + ((stop + dgap) * LOT_LENGTH - CAR_LENGTH - car.curr_loc) / rateDRIV
				if stopped and car.prev.status == 1 and car.prev.stop_idx == curr_vehicle.stop_idx:
					if car.curr_loc >= (stop + dgap) * LOT_LENGTH - CAR_LENGTH + SMALL_INTERVAL: 
						import pdb; pdb.set_trace()
					break
				elif stopped:
					assert car.prev.end_time is not None
					assert car.prev.stop < stop - dgap
					if not (car.prev.stop + dgap) * LOT_LENGTH - CAR_LENGTH >= car.curr_loc - 1e-04:
						import pdb; pdb.set_trace()
					car_time = (stop - car.prev.stop) * LOT_LENGTH / rateDRIV
					if car.prev.status == 1:		
						car_time += max(meanPLIN + car.prev.end_time, self.curr + max(0.0, ((car.prev.stop + dgap) * LOT_LENGTH - CAR_LENGTH - car.curr_loc)) / rateDRIV)
					else:
						car_time += max(car.prev.plin_end, self.curr + max(0.0, ((car.prev.stop + dgap) * LOT_LENGTH - CAR_LENGTH - car.curr_loc)) / rateDRIV)

				if car_time < est_pout_start + meanPOUT - 100 * SMALL_INTERVAL:
					delayed = True
					if (car.status == 1 and car.stop <= stop + dgap) or (car.status == 6 and car.stop < stop + dgap):
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
						if not car.dest_to_stop >= (stop + dgap) * LOT_LENGTH + CAR_LENGTH:
							import pdb; pdb.set_trace()
						car_time = car.calc_time( (stop + dgap) * LOT_LENGTH + CAR_LENGTH )
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

	def check_enter_90(self, enter_time, stop):

		assert spmatch and control == 'full'
		car = self.head
		while car is not None:
			if (spot2blk(stop) < car.block_idx - 1) or (car.status == 6):
				pass
			elif car.status == 2 and car.block_idx == 1:
				assert car.plin_start <= self.curr
				enter_time = max(enter_time, car.plin_end + meanDRIV)
			elif car.status == 2:
				assert car.plin_start <= self.curr
				if (car.block_idx - 1) * CAR_LENGTH / rateDRIV < max(0.0, meanPLIN - (self.curr - car.plin_start)):
					enter_time = max(enter_time, car.plin_start + meanPLIN - (car.block_idx - 1) * CAR_LENGTH / rateDRIV)
			elif car.status == 5 and car.block_idx == 1:
				assert car.pout_start <= self.curr <= car.pout_end
				enter_time = max(enter_time, car.pout_end + meanDRIV)
			elif car.status == 5:
				assert car.pout_start <= self.curr <= car.pout_end
				if ((car.block_idx - 1) * CAR_LENGTH - CAR_LENGTH) / rateDRIV < car.pout_end - self.curr:
					enter_time = max(enter_time, car.pout_end - ((car.block_idx - 1) * CAR_LENGTH - CAR_LENGTH) / rateDRIV)
			else:
				assert car.status == 1
				if enter_time < car.end_time + meanPLIN - (car.block_idx - 1) * CAR_LENGTH / rateDRIV - 10 * SMALL_INTERVAL:
					enter_time = car.end_time + meanPLIN - (car.block_idx - 1) * CAR_LENGTH / rateDRIV
			car = car.nex
		return enter_time
