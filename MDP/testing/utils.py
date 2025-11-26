########################################################################################################################################################
# Last updated: Feb 22, 2021 by Xinyu Liu
# This file defines the helper functions for the Markov Chain simulation.
# This file should not be changed unless you know exactly what you are trying to do.
# WARNING: Unintended changes to this file can have unexpected effects to the simulation results.

########################################################################################################################################################
from math import ceil, sqrt
from collections import *
from params import *
import random 


# Definition for the change points in the trajectory, each consisting of x - the location, t - the time, and v - the speed after t. 
changePoint = namedtuple( 'changePoint', ['x', 't', 'v'])

# Definition of the common random number generator class for service times, speeds, enter and exit maneuver times 
class ParamGen():

	def __init__(self, gen):
		self.gen = gen
		self.count = 0

	def next(self):
		self.count += 1
		val = - 0.1
		while val <= 0.0:
			val = next(self.gen)
		return val

# The RNGs of different distributions are adapted from the built-in "random" package in python.
# The tutorial for this package can be access at https://docs.python.org/3/library/random.html.
# Exponential RNG, which takes rate = 1./mean as input and returns times
def Expo(rate, seed = None):
	curr = random.Random(seed)
	while True:
		yield curr.expovariate(rate)

# Uniform RNG, which takes the rate as input and returns times with a mean of 1./rate, a low of 0., and a high of 2./rate
def Unif(rate, seed = None):
	curr = random.Random(seed)
	while True:
		yield curr.uniform(0, 2. / rate)

# Uniform RNG, which takes the rate as input and returns times with a mean of 1./rate, a low of 0., and a high of 2./rate
def UnifSimple(low, high, seed = None):
	curr = random.Random(seed)
	while True:
		yield curr.uniform(low, high)

# Normal RNG, which takes the rate as input and returns times with a mean of 1./rate and a std.dev of sigma
# if sigma is not specified, the default value is 1./rate = mean which is consistent with Exponential distribution
def Norm(rate, sigma = None, seed = None):
	curr = random.Random(seed)
	if sigma == None:
		sigma = 1. / rate
	while True:
		yield curr.gauss(1. / rate, sigma)	

# Triangular RNG variant 1, which takes the rate as input and returns times with a mode of 1./rate, a low of 0., and a high of 2.*mode.
def TriaSimple(mode, low = None, high = None, seed = None):
	if seed is not None:
		curr = random.Random(seed)
	else:
		curr = random.Random()
	if low == None:
		low = 0
	if high == None:
		high = 2. * mode
	while True:
		yield curr.triangular(low, high, mode)

# Triangular RNG variant 1, which takes the rate as input and returns times with a mode of 1./rate, a low of 0., and a high of 2.*mode.
def Tria1(rate, low = None, high = None, seed = None):
	curr = random.Random(seed)
	mode = 1. / rate
	if low == None:
		low = 0
	if high == None:
		high = 2. * mode
	while True:
		yield curr.triangular(low, high, mode)

# Triangular RNG variant 2, which takes the rate as input and returns times with a mode of 1./rate, a low of mode/2., and a high of 3*mode/2.
def Tria2(rate, low = None, high = None, seed = None):
	curr = random.Random(seed)
	mode = 1. / rate
	if low == None:
		low = .5 * mode
	if high == None:
		high = 1.5 * mode
	while True:
		yield curr.triangular(low, high, mode)

# Triangular RNG variant 3, which takes the rate as input and returns times with a mode of 1./rate, a low of mode/3., and a high of 5*mode/3.
def Tria3(rate, low = None, high = None, seed = None):
	curr = random.Random(seed)
	mode = 1. / rate
	if low == None:
		low = mode / 3.
	if high == None:
		high = 5 * mode / 3.
	while True:
		yield curr.triangular(low, high, mode)

# Constant number generator, which takes the rate as input and returns the speed = rate repeatedly 
def Cons(rate, seed = None):
	while True:
		yield rate

# Discrete Uniform RNG, which takes the rate as input and returns the speed with a mean = rate
def Disc(rate, seed = None):
	curr = random.Random(seed)
	while True:
		yield curr.choice([ 0.9*rate, 0.95*rate, rate, 1.05*rate, 1.1*rate])

# Definition for event class, which records the time, vehicle, and type of an event and computes the precedence order between two well-defined events
class event():

	def __init__(self, time, vehicle, typ):
		if time == None:
			import pdb; pdb.set_trace()
		self.time = time
		self.vehicle = vehicle
		self.typ = typ

	def __eq__(self, other):
		return (float(self.time) == float(other.time)) & (self.typ == other.typ) & (self.vehicle.idx == other.vehicle.idx) 

	def __gt__(self, other):
		if float(self.time) == float(other.time):
			if (self.typ == other.typ):
				if (self.vehicle.idx is None) | (other.vehicle.idx is None):
					assert (self.vehicle.idx is None) and (other.vehicle.idx is None)
					assert (self.typ == 'enter_system')
					return False
				return (self.vehicle.idx > other.vehicle.idx)
			return (event_priority[self.typ] > event_priority[other.typ])
		return (float(self.time) > float(other.time))

	def __lt__(self, other):
		if float(self.time) == float(other.time):
			if self.typ == other.typ:
				if (self.vehicle.idx is None) | (other.vehicle.idx is None):
					assert (self.vehicle.idx is None) and (other.vehicle.idx is None)
					assert (self.typ == 'enter_system')
					return False
				return (self.vehicle.idx < other.vehicle.idx)
			return (event_priority[self.typ] < event_priority[other.typ])
		return (float(self.time) < float(other.time))

	def __ne__(self, other):
		return not self == other

	def __ge__(self, other):
		return not self < other

	def __le__(self, other):
		return not self > other 
		
# Helper function to sort a list of nodes efficiently
def sortHelper(headNode, tlen1, tlen2):

	curr = headNode
	print ('starting...')
	while curr.nex:
		print (curr.data)
		curr = curr.nex
	print (curr.data)
	print ('target length is:', tlen1, tlen2 )

	if headNode != None and headNode.nex != None:
		prev = headNode
		curr = headNode.nex
		assert tlen1 >= 1
		assert tlen2 >= 1
		assert tlen1 >= tlen2
		assert tlen1 - tlen2 <= 1
		for _ in range(tlen1 - 1):
			prev = curr
			curr = curr.nex
		prev.nex = None
			
		curr1 = sortHelper(headNode, tlen1 - tlen1//2, tlen1//2)
		curr2 = sortHelper(curr, tlen2 - tlen2//2, tlen2//2)

		if curr1.data <= curr2.data:
			headNode = curr1
			curr1 = curr1.nex
		else:
			headNode = curr2
			curr2 = curr2.nex
		curr = headNode 

		while curr1 is not None and curr2 is not None:
			if curr1.data <= curr2.data:
				curr.nex = curr1
				curr1 = curr1.nex
			else:
				curr.nex = curr2
				curr2 = curr2.nex
			curr = curr.nex

		if curr1 is not None:
			assert curr2 is None
			curr.nex = curr1
		if curr2 is not None:
			assert curr1 is None
			curr.nex = curr2
	return headNode

# Definition for node class
class Node:

	def __init__(self, data = None):
		self.data = data
		self.prev = None
		self.nex = None

# Definition for doubly linked list class
class DLinkedList:

	def __init__(self):
		self.head = None
		self.end = None
		self.length = 0
		
	def __len__(self):
		return self.length

	def print(self):
		curr = self.head
		while (curr):
			print (curr.data)
			curr = curr.nex

	def isEmpty(self):
		assert (self.head is None) == (self.length == 0)
		return (self.head is None)

	def addHead(self, new):
		NewNode = Node(new)
		NewNode.nex = self.head
		self.head = NewNode
		self.end = NewNode
		self.length += 1

	def addEnd(self, new):
		NewNode = Node(new)
		if self.head is None:
			self.head = NewNode
			return
		curr = self.head
		while (curr.nex):
			curr = curr.nex
		curr.nex = NewNode
		NewNode.prev = curr
		self.end = NewNode
		self.length += 1

	def insert(self, loc, new):
		if loc == 0:
			self.addHead(new)
		elif loc >= self.length:
			self.addEnd(new)
		else:
			NewNode = Node(new)
			prev = self.head
			for _ in range(loc - 1):
				prev = prev.nex
			curr = prev.nex
			prev.nex = NewNode
			prev.nex.prev = prev
			prev.nex.nex = curr
			curr.prev = prev.nex
			self.length += 1

	def getVal(self, idx):
		curr = self.head
		while (nstep == idx > 0):
			curr = curr.nex
		return curr.data

	def __contains__(self, data):
		curr = self.head
		while curr != None:
			if curr.data == data:
				return True
			else:
				return False

# Definition for time class
class time():

	def __init__(self, val):
		self.val = round(val, 3)

	def __str__(self):
		return str(self.val)

	def __eq__(self, other):
		return self.val == other.val

	def __gt__(self, other):
		return self.val > other.val

	def __lt__(self, other):
		return self.val < other.val
