########################################################################################################################################################
# Last updated: Oct 10, 2025 by Xinyu Liu
# This file defines the functions to create the input sets and mappings for different configurations. 
# The functions in this file are specific to MDP problems.
# The configuration is specified by the input parameters from the file params.py.
# So far, 6 types of configurations are included:
# == single- and double-sided; 
# == 0-degree long, 0-degree short, and 90-degree configurations.

########################################################################################################################################################
from math import ceil
from params import *

########################################################################################################################################################
# input parameters m_out(j), m_in(j), g_out(j), mu(j); as a function of spot 
def m_out(j):
	return MOUT

def m_in(j):
	return MIIN

def g_out(j):
	return max(m_out(j)-1, 3)

def mu(j):
	return rateSERV

# a helper function which turns the index of a spot in double-sided configuration
# into an equivalent in single-sided configuration
def doub2sing(j):
	if side == 'single':
		return j
	assert side == 'double'
	return ceil(j / 2)

# a helper function which turns the indices of all the spots around the lane block i
# serving as the set function inverse of i(j) s.t., j \in allspot(i) iff i(j) = i
def allspot(i):
	if side == 'single':
		return list( range( nUnit * (i-1) + 1, nUnit * i + 1 ) )
	else:
		assert side == 'double'
		return list( range( 2 * nUnit * (i-1) + 1, 2 * nUnit * i + 1 ) )

########################################################################################################################################################
# the mapping i+ for block i
# which returns the index of next block if the vehicle in block i would move forward
def blk2blk(i):
	return (i + 1)

########################################################################################################################################################
# the mapping i(j) for spot j
# which returns the lane block next to boarding spot j
def spot2blk(j):
	if side == 'double':
		j = doub2sing(j)
	else:
		assert side == 'single'
	
	assert ceil(j / nUnit) >= 1
	return ceil(j / nUnit)

########################################################################################################################################################
# the mapping i_in(j) for spot j
# which returns the block from which a vehicle starts an enter maneuver into spot j
def spot2in(j):
	if in_direction == 'forward':
		# i.e. if the vehicle enters a boarding spot by going forward
		return spot2blk(j) - 1
	else:
		assert in_direction == 'backward'
		# i.e. if the vehicle enters a boarding spot by going backward
		return spot2blk(j) + 1

########################################################################################################################################################
# the mapping i_out(j) for spot j
# which returns the block in which a vehicle arrives after completion of an exit maneuver from spot j
def spot2out(j):
	if out_direction == 'forward':
		# i.e. if the vehicle exits from a boarding spot by going forward
		return spot2blk(j) + 1
	else:
		assert out_direction == 'backward'
		# i.e. if the vehicle exits from a boarding spot by going backward
		return spot2blk(j) - 1 

########################################################################################################################################################
# the mapping i_a(j) for spot j
# which returns the **furthest forward** block in which a vehicle can be when a decision is made to assign it to spot j
def spot2ass(j, control = 'mdp'):
	if control == 'mdp':
		return max(0, spot2blk(j) - 2)
	else:
		return 0

########################################################################################################################################################
# the mapping J(i) for lane block i := {j \in [N]: i_a(j) >= i } 
# such that J(i) is the subset of spots that a vehicle in block i can be assigned
def assignable(i, N, control = 'mdp'):
	i_ = min(doub2sing(N), i+2)
	if control == 'mdp':
		if i > 0:
			return list( range(min(allspot(i_)),min(N, max(allspot(i_)))+1) )
		return list(     range(1,               min(N, max(allspot(i_)))+1) )
	else:
		assert control == 'full' or control == 'partial'
		if i > 0:
			return []
		return list(range(1, N+1))

########################################################################################################################################################
# the mapping i_wi(j) for spot j wwhen a downstream vehicle has completed k \ge 1 enter manuever steps into spot j
# which returns the **furthest forward** block that a forward-moving vehicle can occupy
def spot2wi(j, k):
	if in_direction == 'forward':
		if k <= 0: # k = -1 or k = 0
			return max(0, spot2blk(j) - 2)
		else:
			return spot2blk(j) - 1
	else:
		assert in_direction == 'backward'
		return spot2blk(j) - 1

########################################################################################################################################################
# the mapping i_wo(j) for spot j wwhen a downstream vehicle has completed k \ge 1 exit manuever steps from spot j
# which returns the **furthest forward** block that a forward-moving vehicle can occupy
def spot2wo(j, k):
	if k <= m_out(j) - 2: # k = 0, 1, \cdots, m_{out}(j-1)-2
		return spot2blk(j) - 1
	else:
		assert k == m_out(j) - 1
		return spot2blk(j)

########################################################################################################################################################
# the implementation for K_oo(j,j') with k
def in_Koo(j, j_, k):

	assert (j_ != j)
	if k < 0 or k >= m_out(j_):
		return False 

	if (angle == 0) and (mode == 'long'):
		if spot2blk(j_) == spot2blk(j) - 1:
			return True
		elif spot2blk(j_) == spot2blk(j) and k <= m_out(j_) - 2:
			assert side == 'double'
			return True
		return False

	if (angle == 0) and (mode == 'short'):
		if spot2blk(j_) == spot2blk(j)-1:
			return True
		elif spot2blk(j_) == spot2blk(j):
			assert side == 'double'
			return True
		elif spot2blk(j_) == spot2blk(j)+1 and k == 0:
			return True
		return False

	if (angle == 90):
		if doub2sing(j) <= doub2sing(j_) <= doub2sing(j)+2:
			if doub2sing(j) == doub2sing(j_):
				assert side == 'double'
			return True
		elif doub2sing(j_) < doub2sing(j) and spot2blk(j_) >= spot2blk(j) - 1:
			return True
		return False

########################################################################################################################################################
# alternative implementation for K_oo(j,j')
def range_Koo(j, j_):
	assert (j_ != j)

	if (angle == 0) and (mode == 'long'):
		if spot2blk(j_) == spot2blk(j) - 1:
			return range(m_out(j_))
		elif spot2blk(j_) == spot2blk(j):
			assert side == 'double'
			return range(m_out(j_)-1)
		return range(0)

	if (angle == 0) and (mode == 'short'):
		if spot2blk(j_) == spot2blk(j)-1:
			return range(m_out(j_))
		elif spot2blk(j_) == spot2blk(j):
			assert side == 'double'
			return range(m_out(j_))
		elif spot2blk(j_) == spot2blk(j)+1:
			return range(1)
		return range(0)

	if (angle == 90):
		if doub2sing(j) <= doub2sing(j_) <= doub2sing(j)+2:
			if doub2sing(j) == doub2sing(j_):
				assert side == 'double'
			return range(m_out(j_))
		elif doub2sing(j_) < doub2sing(j) and spot2blk(j_) >= spot2blk(j) - 1:
			return range(m_out(j_))
		return range(0)

########################################################################################################################################################
# the implementation for K_oi(j,j') with k
def in_Koi(j, j_, k):

	assert (j_ != j)
	if k < 0 or k >= m_in(j_):
		return False 

	if (angle == 0) and (mode == 'long'):
		if spot2blk(j_) == spot2blk(j)+1:
			return True
		elif spot2blk(j_) == spot2blk(j)+2 and k == 0:
			return True
		return False

	if (angle == 0) and (mode == 'short'):
		if spot2blk(j_) == spot2blk(j):
			assert side == 'double'
			return True
		if spot2blk(j_) == spot2blk(j)-1 and k <= m_in(j_)-2:
			return True
		elif spot2blk(j_) == spot2blk(j)+1 and k <= m_in(j_) - m_out(j):
			return True
		return False

	if (angle == 90):
		if doub2sing(j)-2 <= doub2sing(j_) <= doub2sing(j)+2:
			return True
		elif doub2sing(j_) < doub2sing(j)-2 and spot2blk(j_) == spot2blk(j) - 1 and k <= m_in(j_) - 2:
			return True
		elif doub2sing(j_) > doub2sing(j)+2 and spot2blk(j_) == spot2blk(j) + 1 and k <= m_in(j_) - m_out(j):
			return True
		return False

########################################################################################################################################################
# alternative implementation for K_oi(j,j')
def range_Koi(j, j_):

	assert (j_ != j)

	if (angle == 0) and (mode == 'long'):
		if spot2blk(j_) == spot2blk(j)+1:
			return range(m_in(j_))
		elif spot2blk(j_) == spot2blk(j)+2:
			return range(1)
		return range(0)

	if (angle == 0) and (mode == 'short'):
		if spot2blk(j_) == spot2blk(j):
			assert side == 'double'
			return range(m_in(j_))
		if spot2blk(j_) == spot2blk(j)-1:
			return range(m_in(j_) - 1)
		elif spot2blk(j_) == spot2blk(j)+1:
			return range(m_in(j_) - m_out(j) + 1)
		return range(0)

	if (angle == 90):
		if doub2sing(j)-2 <= doub2sing(j_) <= doub2sing(j)+2:
			return range(m_in(j_))
		elif doub2sing(j_) < doub2sing(j)-2 and spot2blk(j_) == spot2blk(j) - 1:
			return range(m_in(j_) - 1)
		elif doub2sing(j_) > doub2sing(j)+2 and spot2blk(j_) == spot2blk(j) + 1:
			return range(m_in(j_) - m_out(j) + 1)
		return range(0)

########################################################################################################################################################
# the implementation for K_io(j,j') with k
def in_Kio(j, j_, k):

	if k < 0 or k >= m_out(j_):
		return False

	if (angle == 0) and (mode == 'long'):
		if spot2blk(j_) == spot2blk(j) and k <= m_out(j_) - 2:
			return True
		elif spot2blk(j_) == spot2blk(j)-1 and k >= m_out(j_) - m_in(j):
			return True
		return False 

	if (angle == 0) and (mode == 'short'):
		if spot2blk(j_) == spot2blk(j):
			return True
		elif -1 <= spot2blk(j_) - spot2blk(j) <= 1 and k >= 1:
			return True 
		return False

	if (angle == 90):
		if j_ == j:
			return True
		elif -1 <= spot2blk(j_) - spot2blk(j) <= 1 and k >= 1:
			if spot2blk(j_) == spot2blk(j):
				assert side == 'double'
			return True 
		return False

########################################################################################################################################################
# alternative implementation for K_io(j,j')
def range_Kio(j, j_):

	if (angle == 0) and (mode == 'long'):
		if spot2blk(j_) == spot2blk(j):
			return range(m_out(j_)-1)
		elif spot2blk(j_) == spot2blk(j)-1:
			return range(m_out(j_)-m_in(j), m_out(j_))
		return range(0) 

	if (angle == 0) and (mode == 'short'):
		if spot2blk(j_) == spot2blk(j):
			return range(m_out(j_))
		elif -1 <= spot2blk(j_) - spot2blk(j) <= 1:
			return range(1, m_out(j_))
		return range(0)

	if (angle == 90):
		if j_ == j:
			return range(m_out(j_))
		elif -1 <= spot2blk(j_) - spot2blk(j) <= 1 and k >= 1:
			if spot2blk(j_) == spot2blk(j):
				assert side == 'double'
			return range(1, m_out(j_))
		return range(0)

########################################################################################################################################################
# the implementation for K_ii(j,j') with k
def in_Kii(j, j_, k):

	assert j_ != j
	if k < 0 or k >= m_in(j_):
		return False

	if (angle == 0) and (mode == 'long'):
		if spot2blk(j_) == spot2blk(j) and k >= 1:
			assert side == 'double'
			return True
		return False

	if (angle == 0) and (mode == 'short'):
		if -1 <= spot2blk(j_) - spot2blk(j) <= 1 and k >= 1:
			return True
		return False

	if (angle == 90):
		if -1 <= spot2blk(j_) - spot2blk(j) <= 1 and k >= 1:
			return True
		return False

########################################################################################################################################################
# alternative implementation for K_ii(j,j')
def range_Kii(j, j_):

	assert j_ != j

	if (angle == 0) and (mode == 'long'):
		if spot2blk(j_) == spot2blk(j):
			assert side == 'double'
			return range(1, m_in(j_))
		return range(0)

	if (angle == 0) and (mode == 'short'):
		if -1 <= spot2blk(j_) - spot2blk(j) <= 1:
			return range(1, m_in(j_))
		return range(0)

	if (angle == 90):
		if -1 <= spot2blk(j_) - spot2blk(j) <= 1:
			return range(1, m_in(j_))
		return range(0)

########################################################################################################################################################
# implementation of E_out(j,j_) with k step of completed exit maneuver steps from j_
# checking for conflicts with scheduling an exit maneuver from j

def in_Eout(j, j_, k):

	assert (j != j_)
	if k < 0 or k >= m_out(j_):
		return False 

	if (angle == 0) and (mode == 'long'):
		if spot2blk(j_) == spot2blk(j) and k <= m_out(j_) - 2:
			return True
		if spot2blk(j_) == spot2blk(j) - 1:
			return True
		if spot2blk(j_) < spot2blk(j) - 1 and k >= spot2blk(j) - spot2blk(j_):
			return True
		return False

	if (angle == 0) and (mode == 'short'):
		if spot2blk(j_) == spot2blk(j):
			return True
		if spot2blk(j_) == spot2blk(j) + 1 and k == 0:
			return True
		if spot2blk(j_) == spot2blk(j) - 1:
			return True
		if spot2blk(j_) < spot2blk(j) - 1 and k >= spot2blk(j) - spot2blk(j_):
			return True
		return False

	if (angle == 90):
		if doub2sing(j) - (nUnit-1) <= doub2sing(j_) <= doub2sing(j) + (nUnit-1):
			return True
		elif doub2sing(j_) < doub2sing(j) - (nUnit-1) and k >= spot2blk(j) - spot2blk(j_):
			return True
		return False

########################################################################################################################################################
# alternative implementation of E_out(j,j_)

def range_Eout(j, j_):

	assert (j != j_)

	if (angle == 0) and (mode == 'long'):
		if spot2blk(j_) == spot2blk(j):
			return range( m_out(j_) - 1 )
		if spot2blk(j_) == spot2blk(j) - 1:
			return range( m_out(j_) )
		if spot2blk(j_) < spot2blk(j) - 1:
			return range( spot2blk(j) - spot2blk(j_), m_out(j_) )
		return range(0)

	if (angle == 0) and (mode == 'short'):
		if spot2blk(j_) == spot2blk(j):
			return range( m_out(j_) )
		if spot2blk(j_) == spot2blk(j) + 1:
			return range(1)
		if spot2blk(j_) == spot2blk(j) - 1:
			return range( m_out(j_) )
		if spot2blk(j_) < spot2blk(j) - 1 and k >= spot2blk(j) - spot2blk(j_):
			return range( spot2blk(j) - spot2blk(j_), m_out(j_))
		return range(0)

	if (angle == 90):
		if doub2sing(j) - (nUnit-1) <= doub2sing(j_) <= doub2sing(j) + (nUnit-1):
			return range( m_out(j_) )
		elif doub2sing(j_) < doub2sing(j) - (nUnit-1):
			return range( spot2blk(j) - spot2blk(j_), m_out(j_) )
		return range(0)

########################################################################################################################################################
# implementation of E_in(j, j_) with k steps of completed enter maneuver steps from j_
# checking for conflicts with scheduling an exit maneuver from j

def in_Ein(j, j_, k):

	assert (j != j_)
	if k < 0 or k >= m_in(j_):
		return False

	if (angle == 0) and (mode == 'long'):
		if spot2blk(j_) == spot2blk(j) + 1 and k >= 1:
			return True
		return False

	if (angle == 0) and (mode == 'short'):
		if spot2blk(j_) == spot2blk(j) and k >= 1:
			return True
		if spot2blk(j_) == spot2blk(j) - 1 and 1 <= k <= m_in(j_) - 2:
			return True
		elif spot2blk(j_) == spot2blk(j) + 1 and 1 <= k <= m_in(j_) - m_out(j):
			return True
		elif spot2blk(j_) >= spot2blk(j) + 2 and 1 <= k <= m_in(j_) - m_out(j) - (spot2blk(j_) - spot2blk(j)) + 1:
			return True
		return False

	if (angle == 90):
		if doub2sing(j) - (nUnit-1) <= doub2sing(j_) <= doub2sing(j) + (nUnit-1) and k >= 1:
			return True
		elif doub2sing(j_) > doub2sing(j) + (nUnit-1) and spot2blk(j_) == spot2blk(j) + 1 and 1 <= k <= m_in(j_) - m_out(j):
			return True
		elif spot2blk(j_) >= spot2blk(j) + 2 and 1 <= k <= m_in(j_) - m_out(j) - (spot2blk(j_) - spot2blk(j)) + 1:
			return True
		return False

########################################################################################################################################################
# alternative implementation for E_in(j,j_)

def range_Ein(j, j_):

	assert (j != j_)

	if (angle == 0) and (mode == 'long'):
		if spot2blk(j_) == spot2blk(j) + 1:
			return range( 1, m_in(j_) )
		return range(0)

	if (angle == 0) and (mode == 'short'):
		if spot2blk(j_) == spot2blk(j):
			return range( 1, m_in(j_) )
		if spot2blk(j_) == spot2blk(j) - 1:
			return range( 1, m_in(j_) - 1 )
		elif spot2blk(j_) == spot2blk(j) + 1:
			return range( 1, m_in(j_) - m_out(j) + 1 )
		elif spot2blk(j_) >= spot2blk(j) + 2:
			return range( 1, m_in(j_) - m_out(j) - (spot2blk(j_) - spot2blk(j)) + 2 )
		return range(0)

	if (angle == 90):
		if doub2sing(j) - (nUnit-1) <= doub2sing(j_) <= doub2sing(j) + (nUnit-1) :
			return range( 1, m_in(j_) )
		elif doub2sing(j_) > doub2sing(j) + (nUnit-1) and spot2blk(j_) == spot2blk(j) + 1:
			return range( 1, m_in(j_) - m_out(j) + 1 )
		elif spot2blk(j_) >= spot2blk(j) + 2:
			return range( 1, m_in(j_) - m_out(j) - (spot2blk(j_) - spot2blk(j)) + 2 )
		return range(0)

########################################################################################################################################################
# implementation of I_out(j,j_) with block i in which a vehicle is heading to j_
# checking for conflicts with scheduling an exit maneuver from spot j

def in_Iout(j, j_, i, n):

	assert 0 <= i <= n

	if (angle == 0) and (mode == 'long'):
		if spot2blk(j_) >= spot2blk(j) and j_ != j and spot2blk(j) - 1 - g_out(j) <= i <= spot2blk(j):
			return True
		return False

	if (angle == 0) and (mode == 'short'):
		if (j_ != j) and (spot2blk(j) - 1 <= spot2blk(j_) <= spot2blk(j) + 1) and spot2blk(j) - 1 - g_out(j) <= i <= spot2in(j_):
			return True
		if (spot2blk(j_) >= spot2blk(j) + 2) and spot2blk(j) - 1 - g_out(j) <= i <= spot2in(j_) and i <= spot2blk(j) + m_in(j_) - m_out(j) + 2:
			return True
		return False

	if (angle == 90):
		if (j_ != j) and (doub2sing(j) - (nUnit-1) <= doub2sing(j_) <= doub2sing(j) + (nUnit-1)) and spot2blk(j) - 1 - g_out(j) <= i <= spot2in(j_):
			return True 
		if (j_ != j) and (doub2sing(j_) < doub2sing(j) - (nUnit-1)) and (spot2blk(j) - 1 <= spot2blk(j_)) and spot2blk(j) - 1 - g_out(j) <= i <= spot2in(j_):
			return True
		if doub2sing(j_) > doub2sing(j) + (nUnit-1) and spot2blk(j_) > spot2blk(j) and spot2blk(j) - 1 - g_out(j) <= i <= spot2in(j_) and i <= spot2blk(j) + m_in(j_) - m_out(j) + 2:
			return True
		return False

########################################################################################################################################################
# alternative implementation of I_out(j,j_)

def range_Iout(j, j_, n):

	if (angle == 0) and (mode == 'long'):
		if spot2blk(j_) >= spot2blk(j) and j_ != j:
			return range( max(0, spot2blk(j) - 1 - g_out(j)), spot2blk(j) + 1 )
		return range(0)

	if (angle == 0) and (mode == 'short'):
		if (j_ != j) and (spot2blk(j) - 1 <= spot2blk(j_) <= spot2blk(j) + 1):
			return range( max(0, spot2blk(j) - 1 - g_out(j)), min(n, spot2in(j_)) + 1 )
		if (spot2blk(j_) >= spot2blk(j) + 2):
			return range( max(0, spot2blk(j) - 1 - g_out(j)), min(n, spot2in(j_), spot2blk(j) + m_in(j_) - m_out(j) + 2) + 1 )
		return range(0)

	if (angle == 90):
		if (j_ != j) and (doub2sing(j) - (nUnit-1) <= doub2sing(j_) <= doub2sing(j) + (nUnit-1)):
			return range( max(0, spot2blk(j) - 1 - g_out(j)), min(n, spot2in(j_)) + 1 ) 
		if (j_ != j) and (doub2sing(j_) < doub2sing(j) - (nUnit-1)) and (spot2blk(j) - 1 <= spot2blk(j_)):
			return range( max(0, spot2blk(j) - 1 - g_out(j)), min(n, spot2in(j_)) + 1 )
		if doub2sing(j_) > doub2sing(j) + (nUnit-1) and spot2blk(j_) > spot2blk(j):
			return range( max(0, spot2blk(j) - 1 - g_out(j)), min(n, spot2in(j_), spot2blk(j) + m_in(j_) - m_out(j) + 2) + 1 ) 
		return range(0)

########################################################################################################################################################
# implementation for I_in(j,j_) with block i in which a vehicle is heading to j_
# checking for conflicts with allowing a replacement vehicle enter and be assigned to j

def in_Iin(j, j_, i, n):

	assert j != j_
	assert 0 <= i <= n

	if (angle == 0) and (mode == 'long') :
		if (spot2blk(j_) <= spot2blk(j)):
			return (1 <= i <= m_in(j_) - 1)	
		return False

	if spot2blk(j_) <= spot2blk(j) + 1:
		return (1 <= i <= m_in(j_) + 1)
	
	return False

########################################################################################################################################################
# alternatvie implementation of I_in(j,j_)

def range_Iin(j, j_, n):

	assert j != j_

	if (angle == 0) and (mode == 'long') :
		if (spot2blk(j_) <= spot2blk(j)):
			return range(1, min(n+1, m_in(j_)))	
		return range(0)

	if spot2blk(j_) <= spot2blk(j) + 1:
		return range(1, min(n+1, m_in(j_)+2))
	
	return range(0)

########################################################################################################################################################
# implementation of K_out(j,j_) with k completed exit maneuver steps from spot j_
# when checking for a replacement vehicle assigned to spot j

def in_Kout(j, j_, k):

	if k < 0 or k >= m_out(j_):
		return False

	if j_ == j:
		return (k <= m_out(j_) - spot2blk(j_) - 1)

	if (angle == 0) and (mode == 'long'):
		if spot2blk(j_) <= spot2blk(j):
			return (1 <= k <= m_out(j_) - spot2blk(j_) - 1)
		return False

	if spot2blk(j_) <= spot2blk(j) + 1:
		return (1 <= k <= m_out(j_) - spot2blk(j_) - 1)

	return False

########################################################################################################################################################
# alternative implementation of K_out(j,j_)

def range_Kout(j, j_):

	if j_ == j:
		return range( m_out(j_) - spot2blk(j_) )

	if (angle == 0) and (mode == 'long'):
		if doub2sing(j_) <= doub2sing(j):
			return range(1, m_out(j_) - spot2blk(j_))
		return range(0)

	if spot2blk(j_) <= spot2blk(j) + 1:
		return range(1, m_out(j_) - spot2blk(j_))

	return range(0)

########################################################################################################################################################
# implementation of K_in(j,j_) with k completed enter maneuver steps
# when checking for a replacement vehicle assigned to spot j

def in_Kin(j, j_, k):

	assert j != j_
	if k < 0 or k >= m_in(j_):
		return False
	
	if spot2blk(j_) <= spot2blk(j):
		return (1 <= k <= m_in(j_) - spot2blk(j_))

	if (angle == 0 and mode == 'short') or (angle == 90):
		if spot2blk(j_) <= spot2blk(j)+1:
			return (1 <= k <= m_in(j_) - spot2blk(j_))

	return False

########################################################################################################################################################
# alternative implementation of K_in(j, j_)

def range_Kin(j, j_):

	assert j != j_ 

	if spot2blk(j_) <= spot2blk(j):
		return range(1, m_in(j_) - spot2blk(j_) + 1)

	if (angle == 0 and mode == 'short') or (angle == 90):
		if spot2blk(j_) <= spot2blk(j)+1:
			return range(1, m_in(j_) - spot2blk(j_) + 1)

	return range(0)
