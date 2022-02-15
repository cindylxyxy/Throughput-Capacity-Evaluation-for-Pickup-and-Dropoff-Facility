########################################################################################################################################################
# Last updated: Feb 06, 2022 by Xinyu Liu
# This file defines the functions to create the input sets and mappings for different configurations. 
# The functions in this file are specific to MC problems with full or partial access control modes. 
# The configuration is specified by the input parameters from the file params.py.
# So far, 6 types of configurations are included:
# == single- and double-sided; 
# == 0-degree long, 0-degree short, and 90-degree configurations.

########################################################################################################################################################
from math import ceil
from params import *

########################################################################################################################################################
# [MDP specific] a helper function which turns the indices of all the boarding spots around the lane block i
def allspot(i):
	if angle == 0 and side == 'single':
		return [i]
	elif angle == 0:
		assert side == 'double' 
		return [2 * i - 1, 2 * i]
	else: 
		assert angle == 90 
		if side == 'single':
			return [3 * i - 2, 3 * i - 1, 3 * i]
		else:
			assert side == 'double'
			return list( range(6 * i - 5, 6 * i + 1) )

########################################################################################################################################################
# [MDP specific] a helper function which turns the indices of all the boarding spots that a vehicle in lane block i can be assigned
def enterspot(i):
	if angle == 0 and mode == 'long':
		return allspot(i + g_in + 1)
	else:
		return allspot(i + g_in - 1)

########################################################################################################################################################
# a helper function which turns the index of a spot in double-sided configuration
# into an equivalent in single-sided configuration
def doub2sing(j):
	assert side == 'double'
	return ceil(j / 2)

def idx2spot(j):
	if side == 'double':
		return ceil(j / 2)
	else:
		assert side == 'single'
		return j

########################################################################################################################################################
# the function i(j) for spot j in {1, ..., N} 
# which returns the lane block next to boarding spot j
# a vehicle in any configuration makes an exit maneuver from spot j to block i(j)
def spot2blk(j):
	if side == 'double':
		j = doub2sing(j)
	else:
		assert side == 'single' 
	if nUnit == 1:
		return j
	else:
		return ceil(j / 3)

########################################################################################################################################################
# the function i_in(j) for spot j in {1, ..., N}
# which returns the block from which a vehicle starts an enter maneuver into spot j
def spot2in(j):
	if angle == 0 and mode == 'long':
		# i.e. if the vehicle enters a boarding spot by going forward
		return spot2blk(j) - 1
	else:
		assert (angle == 0 and mode == 'short') or (angle == 90)
		# i.e. if the vehicle enters a boarding spot by going backward
		return spot2blk(j) + 1 

########################################################################################################################################################
# the function i+ for block i in {0, 1, ..., n}
# which returns the index of next block if the vehicle in block i would move forward
def blk2blk(i):
	return (i + 1)

########################################################################################################################################################
# the implementation of E_out for (j, k) in {1, ..., N}^2, j != k, m in {0, 1, ..., m_out - 1}
# i.e. the set of steps of the exit maneuver of a vehicle in spot k 
# s.t. a vehicle in spot j cannot start its exit maneuver without conflict and without its own movement to the exit being delayed
# which returns True if the vehicle in j should not complete the first step of its exit maneuver if a vehicle in k completed m steps of its exit maneuver 
# i.e. True if m in E_out(j, k) and False if m not in E_out(j, k) 
def in_Eout(j, k, m):

	assert (j != k)
	if side == 'double':
		J = doub2sing(j)
		K = doub2sing(k)
	else:
		assert side == 'single'
		J = j
		K = k

	if (K < J - dgap) and (spot2blk(k) > spot2blk(j) - m_out):
		return (spot2blk(j) - spot2blk(k) <= m < m_out)

	if (K > J + dgap) or ( (K < J - dgap) and (spot2blk(k) <= spot2blk(j) - m_out) ):
		return False

	if (angle == 0) and (mode == 'long') and (K == J + 1):
		return False

	if (angle == 0) and (mode == 'long') and (K == J):
		assert side == 'double'
		return (0 <= m < m_out - 1)

	if (angle == 0) and (mode == 'short') and (K == J + 1):
		assert dgap == 1
		return (m == 0)

	return (0 <= m < m_out)

########################################################################################################################################################
# the alternative set implementation of E_out for (j, k) in {1, ..., N}^2, j != k
def range_Eout(j, k):
	
	assert (j != k)
	if side == 'double':
		J = doub2sing(j)
		K = doub2sing(k)
	else:
		assert side == 'single'
		J = j
		K = k

	if (K < J - dgap) and (spot2blk(k) > spot2blk(j) - m_out):
		return range(spot2blk(j) - spot2blk(k), m_out)	

	if (K > J + dgap) or ( (K < J - dgap) and (spot2blk(k) <= spot2blk(j) - m_out) ):
		return range(0)

	if (angle == 0) and (mode == 'long') and (K == J + 1):
		return range(0)

	if (angle == 0) and (mode == 'long') and (K == J):
		assert side == 'double'
		return range(m_out - 1)

	if (angle == 0) and (mode == 'short') and (K == J + 1):
		assert dgap == 1
		return range(1)

	return range(m_out)

########################################################################################################################################################
# the implementation of E_in for (j, k) in {1, ..., N}^2, j != k, m in {1, ..., m_in - 1}
# i.e. the set of steps of the enter maneuver of a vehicle in spot k 
# s.t. a vehicle in spot j cannot start its exit maneuver without conflict and without its own movement to the exit being delayed
# which returns True if the vehicle in j should not complete the first step of its exit maneuver if a vehicle in k completed m steps of its enter maneuver 
# i.e. True if m in E_in(j, k) and False if m not in E_in(j, k) 
def in_Ein(j, k, m):

	assert (j != k)
	if side == 'double':
		J = doub2sing(j)
		K = doub2sing(k)
	else:
		assert side == 'single'
		J = j
		K = k

	if (angle == 0) and (mode == 'long'):
		assert m_in <= m_out
		if K == J + 1:
			return (1 <= m < m_in)
		return False

	if (angle == 0) and (mode == 'short'):
		assert m_in > m_out
		assert m_in > 3
		# we assumed m = m_in - 2
		# thus m_in > 3 is required for the set not to be empty
		if (K == J):
			return (1 <= m < m_in)
		if (K == J - 1):
			return (1 <= m < m_in - 1)
		if (K == J + 1):
			return (1 <= m < m_in - m_out + 1) 
		if (K >= J + 2) and (m_in - m_out >= K - J):
			return (1 <= m < m_in - m_out - (K - J) + 2 )
		return False

	assert angle == 90
	assert m_in > m_out
	assert m_in > 3
	# we assumed m = m_in - 2
	# thus m_in > 3 is required for the set not to be empty
	if (K == J):
		return (1 <= m < m_in)
	elif (J - dgap <= K <= J + dgap):
		return (1 <= m < m_in)
	elif (spot2blk(k) == spot2blk(j) + 1):
		return (1 <= m < m_in - m_out + 1) 
	elif (spot2blk(k) >= spot2blk(j) + 2) and (m_in - m_out >= spot2blk(k) - spot2blk(j)):
		return (1 <= m < m_in - m_out - (spot2blk(k) - spot2blk(j)) + 2 )
	return False

########################################################################################################################################################
# the alternative set implementation of E_in for (j, k) in {1, ..., N}^2, j != k, m in {1, ..., m_in - 1}
def range_Ein(j, k):

	assert (j != k)
	if side == 'double':
		J = doub2sing(j)
		K = doub2sing(k)
	else:
		assert side == 'single'
		J = j
		K = k

	if (angle == 0) and (mode == 'long'):
		assert m_in <= m_out
		if K == J + 1:
			return range(1, m_in)
		return range(0)

	if (angle == 0) and (mode == 'short'):
		assert m_in > m_out
		assert m_in > 3
		# we assumed m = m_in - 2 
		# thus m_in > 3 is required for the set not to be empty
		if (K == J):
			return range(1, m_in)
		if (K == J - 1):
			return range(1, m_in - 1)
		if (K == J + 1):
			return range(1, m_in - m_out + 1)
		if (K >= J + 2) and (m_in - m_out >= K - J):
			return range(1, m_in - m_out - (K - J) + 2 )
		return range(0)

	assert angle == 90
	assert m_in > m_out
	assert m_in > 3
	# we assumed m = m_in - 2
	# thus m_in > 3 is required for the set not to be empty
	if (K == J):
		return range(1, m_in)
	elif (J - dgap <= K <= J + dgap):
		return range(1, m_in)
	elif (spot2blk(k) == spot2blk(j) + 1):
		return range(1, m_in - m_out + 1) 
	elif (spot2blk(k) >= spot2blk(j) + 2) and (m_in - m_out >= spot2blk(k) - spot2blk(j)):
		return range(1, m_in - m_out - (spot2blk(k) - spot2blk(j)) + 2 )
	return range(0)

########################################################################################################################################################
# the implementation of I_out for (j, k) in {1, ..., N} * {1, ..., N + 1}, m in {0, 1, ..., n}, and n is the total number of lane blocks
# i.e. the set of block m s.t. if a vehicle in block m is traveling to spot k (or the exit if k = N + 1), a vehicle in spot j 
# cannot start its exit maneuver without delaying the movement of the vehicle to spot k and without its own movement to the exit being delayed
# which returns True if the vehicle in j should not complete the first step of its exit maneuver if a vehicle in block m is traveling to k
# i.e. True if m in I_out(j, k) and False if m not in I_out(j, k) 
def in_Iout(j, k, m, n):

	if control == 'partial':
		return in_Hout(j, k, m, n)

	if (j == k):
		return False

	if side == 'double':
		J = doub2sing(j)
		K = doub2sing(k)
	else:
		assert side == 'single'
		J = j
		K = k

	if (angle == 0) and (mode == 'long') and (J < K):	
		assert m_in <= m_out
		return ( max(0, J - g_out + 1) <= m < J + 1 )

	if (angle == 0) and (mode == 'short'):
		assert m_in > m_out
		if (J - 1 <= K <= J + 1):
			return ( max(0, J - g_out + 1) <= m < min(n, spot2in(k)) + 1 )
		if (J + 2 <= K):
			return ( max(0, J - g_out + 1) <= m < min(n, spot2in(k), J + m_in - m_out + 2) + 1 )

	if (angle == 90):
		assert m_in > m_out
		if (J - dgap <= K <= J + dgap):
			return ( max(0, spot2blk(j) - g_out + 1) <= m < min(n, spot2in(k)) + 1 )
		if (spot2blk(j) - 1 == spot2blk(k)):
			return ( max(0, spot2blk(j) - g_out + 1) <= m < min(n, spot2in(k)) + 1 )
		if (spot2blk(j) < spot2blk(k)):
			return ( max(0, spot2blk(j) - g_out + 1) <= m < min(n, spot2in(k), spot2blk(j) + m_in - m_out + 2) + 1)
	return False

########################################################################################################################################################
# the alternative set implementation of I_out for (j, k) in {1, ..., N} * {1, ..., N + 1}, m in {0, 1, ..., n}, and n = i(N) is the total number of lane blocks
def range_Iout(j, k, n):

	if control == 'partial':
		return range_Hout(j, k, n)

	if (j == k):
		return range(0)

	if side == 'double':
		J = doub2sing(j)
		K = doub2sing(k)
	else:
		assert side == 'single'
		J = j
		K = k

	if (angle == 0) and (mode == 'long') and (J < K):
		assert m_in <= m_out
		return range( max(0, J - g_out + 1), J + 1)

	if (angle == 0) and (mode == 'short'):
		assert m_in > m_out
		if (J - 1 <= K <= J + 1):
			return range( max(0, J - g_out + 1), min(n, spot2in(k)) + 1 )
		if (J + 2 <= K):
			return range( max(0, J - g_out + 1), min(n, spot2in(k), J + m_in - m_out + 2) + 1 )

	if (angle == 90):
		assert m_in > m_out
		if (J - dgap <= K <= J + dgap):
			return range( max(0, spot2blk(j) - g_out + 1), min(n, spot2in(k)) + 1 )
		if (spot2blk(j) - 1 == spot2blk(k)):
			return range( max(0, spot2blk(j) - g_out + 1), min(n, spot2in(k)) + 1 )
		if (spot2blk(j) < spot2blk(k)):
			return range( max(0, spot2blk(j) - g_out + 1), min(n, spot2in(k), spot2blk(j) + m_in - m_out + 2) + 1)
	
	return range(0)

########################################################################################################################################################
# the implementation of \hat{I}_out for (j, k) in {1, ..., N} * {1, ..., N + 1}, m in {0, 1, ..., n}, and n is the total number of lane blocks
# \hat{I}_out is a variant of I_out under the partial access control
def in_Hout(j, k, m, n):

	if (j == k):
		return False

	if side == 'double':
		J = doub2sing(j)
		K = doub2sing(k)
	else:
		assert side == 'single'
		J = j
		K = k

	if (angle == 0) and (mode == 'long') and (J < K):
		assert m_in <= m_out
		return ( max(1, J - g_out + 1) <= m < J + 1 )

	if (angle == 0) and (mode == 'short'):
		assert (m_in > m_out)
		if (J + 1 >= K >= J - 1):
			return ( max(1, J - g_out + 1) <= m < min(n, spot2in(k)) + 1 )
		if (K > J + 1):
			return range( max(1, J - g_out + 1), min(n, spot2in(k), J + m_in - m_out + 2) + 1 )

	if (angle == 90):
		assert (m_in > m_out)
		if (J + dgap >= K >= J - dgap) or (spot2blk(j) + 1 == spot2blk(k)) or (spot2blk(j) - 1 == spot2blk(k)):
			return ( max(1, spot2blk(j) - g_out + 1) <= m < min(n, spot2in(k)) + 1 )
		if (spot2blk(k) > spot2blk(j) + 1):
			return ( max(1, spot2blk(j) - g_out + 1) <= m < min(n, spot2in(k), spot2blk(j) + m_in - m_out + 2) + 1)	

	return False

########################################################################################################################################################
# the alternative set implementation of \hat{I}_out for (j, k) in {1, ..., N} * {1, ..., N + 1}, m in {0, 1, ..., n}, and n is the total number of lane blocks
def range_Hout(j, k, n):

	if (j == k):
		return range(0)

	if side == 'double':
		J = doub2sing(j)
		K = doub2sing(k)
	else:
		assert side == 'single'
		J = j
		K = k

	if (angle == 0) and (mode == 'long') and (J < K):
		assert m_in <= m_out
		return range(max(1, J - g_out + 1), J + 1)

	if (angle == 0) and (mode == 'short') and (K >= J - 1):
		assert (m_in > m_out)
		if (J + 1 >= K >= J - 1):
			return range( max(1, J - g_out + 1), min(n, spot2in(k)) + 1 )
		if (K > J + 1):
			return range( max(1, J - g_out + 1), min(n, spot2in(k), J + m_in - m_out + 2) + 1 )

	if (angle == 90):
		assert (m_in > m_out)
		if (J + dgap >= K >= J - dgap) or (spot2blk(j) + 1 == spot2blk(k)) or (spot2blk(j) - 1 == spot2blk(k)):
			return range( max(1, spot2blk(j) - g_out + 1), min(n, spot2in(k)) + 1 )
		if (spot2blk(k) > spot2blk(j) + 1):
			return range( max(1, spot2blk(j) - g_out + 1), min(n, spot2in(k), spot2blk(j) + m_in - m_out + 2) + 1)	
	return range(0)

########################################################################################################################################################
# the implementation of \tau_s for (i, j) in \hat{I}_out(j, k) * {1, ..., N} and s = (x, y)
# \tau_s(i,j) returns the time that it will take a vehicle in block i (on its way to spot k) to reach block i(j)
# and it depends on how long the vehicle in block will be delayed.
# This function returns True if \tau_s(i,j) < g_out and False otherwise.
# When there is no existing delay on the vehicle in block i, \tau_s(i,j) = i(j) - i < g_out 
def tau(x, y, i, j):
	for i_ in range(i + 1, spot2blk(j) + 1):
		assert y[i_-1] != j
		if 1 <= spot2blk(y[i_-1]) <= spot2blk(j) and (spot2in(y[i_-1]) - i_ + m_in) + (spot2blk(j) - (spot2blk(y[i_-1]) - 1)) >= g_out:
			return False
	for j_ in range(j - 1, 0, -1):
		if 2 + m_out <= x[j_-1] <= m_in + m_out and spot2blk(j_) > i and (m_in - (x[j_-1] - m_out - 1)) + (spot2blk(j) - (spot2blk(j_) - 1)) >= g_out:
			return False
	return True

########################################################################################################################################################
# the implementation of \gamma_s for j in {1, ..., N} and s = (x, y) when m_in > m_out
# \gamma_s(j) returns whether a vehicle in spot j can start its exit maneuver and arrive at lane block i(j) + 1 upon completion of the maneuver
# and it depends on the state of the system i.e. how long each vehicle between itself and an entering vehicle would be delayed
def gamma(x, y, j):

	assert m_in > m_out
	assert control == 'partial'
	assert (angle == 0 and mode == 'short') or angle == 90

	if y[blk2blk(spot2blk(j)) - 1] == 0:
		# i.e. if the block i(j)+ is empty then vehicle in j can start its exit maneuver
		# if a vehicle is in block i(j)++ making or about to make an enter maneuver, the vehicle in j would have been blocked by E_in(j, k) s.t. i(k) = i(j)+ 
		return True

	i_ = blk2blk(spot2blk(j))
	enter = []
	n = len(y)
	N = len(x)

	assert n == spot2blk(N) + 1

	while blk2blk(i_) <= n and y[blk2blk(i_) - 1] > 0:
		if y[blk2blk(i_) - 1] < N + 1:
			enter.append( tuple((blk2blk(i_), y[blk2blk(i_) - 1])) ) 
		i_ = blk2blk(i_)
	
	if blk2blk(i_) > n:
		# i.e. all lane blocks ahead of i(j) are occupied and the last one in lane block n = i(N) + 1 must be entering a spot k s.t. i(k) = i(N)
		return False

	assert y[blk2blk(i_) - 1] == 0 and y[i_ - 1] > 0

	if blk2blk(blk2blk(i_)) <= n:
		if (spot2in(y[blk2blk(blk2blk(i_)) - 1]) == blk2blk(blk2blk(i_))):
			# i.e. if the vehicle in lane block i_++ is planning to enter a spot beside block i_+ thus the block i_+ is empty
			return False
		if (0 < y[blk2blk(blk2blk(i_)) - 1] <= N) and (2 + m_out <= x[y[blk2blk(blk2blk(i_)) - 1] - 1] <= m_in + 1):
			# i.e. if the vehicle is entering a spot beside block i_+ thus the block i_+ is empty
			# note that if x[y[blk2blk(blk2blk(i_)) - 1] - 1] >= m_in + 2, it is only left with less than (m_in + m_out + 1) - (m_in + 2) = m_out - 1 steps 
			# after m_out - 1 lambda-transitions, this vehicle finishes the enter maneuver
			# after m_out lambda-transitions, a vehicle previously in i_ can move to i_+, thus the vehicle in j can move to i(j) + 1 at its last step
			return False
		if (0 == y[blk2blk(blk2blk(i_)) - 1]):
			for j_ in range(N, 0, -1):
				if spot2blk(j_) <= i_:
					break
				if spot2blk(j_) == blk2blk(i_) and (2 + m_out <= x[j_ - 1] <= m_in + 1):
					return False

	if y[i_ - 1] <= N:
		# i.e. if a vehicle in i_ is entering, it cannot enter any spots k s.t. spot2in(k) == i_ because the block behind itself is occupied
		# it can neither enter any spots k s.t. spot2in(k) < i_
		assert spot2in(y[i_ - 1]) > i_
		if spot2in(y[i_ - 1]) == blk2blk(i_):
			# i.e. i(y_{i_-1}) == i_ the vehicle is right beside its assigned spot, and it will take 1 + m_in - m_out to finish the enter maneuver 
			return False

	if len(enter) >= 1:
		for (block_idx, spot_idx) in enter:
			if spot2blk(spot_idx) == block_idx:
				return False
	
	return True

########################################################################################################################################################
# the implementation of I_in for (j, k) in {1, ..., N}^2, m in {1, ..., n}, and n is the total number of lane blocks
# i.e. the set of block m s.t. if a vehicle in block m is on its way to spot k, 
# then another vehicle cannot enter the facility and be assigned to spot j without being delayed
# which returns True if a vehicle should enter and be assigned to spot j when a vehicle in block m is traveling to spot k
# i.e. True if m in I_in(j, k) and False if m not in I_in(j, k) 
def in_Iin(j, k, m, n):

	if (j == k):
		return False

	if side == 'double':
		J = doub2sing(j)
		K = doub2sing(k)
	else:
		assert side == 'single'
		J = j
		K = k

	if (angle == 0) and (mode == 'long'):
		if (J >= K): 
			return (1 <= m < min(n, m_in - 1) + 1)
		return False

	if spot2blk(j) >= spot2blk(k) - 1:
		return (1 <= m < min(n, m_in + 1) + 1)

	return False

########################################################################################################################################################
# the alternatvie set implementation of I_in for (j, k) in {1, ..., N}^2, m in {1, ..., n, n+1}
def range_Iin(j, k, n):

	if (j == k):
		return range(0)

	if side == 'double':
		J = doub2sing(j)
		K = doub2sing(k)
	else:
		assert side == 'single'
		J = j
		K = k

	if (angle == 0) and (mode == 'long'):
		if (J >= K): 
			return range(1, min(n, m_in - 1) + 1)
		return range(0)

	if spot2blk(j) >= spot2blk(k) - 1:
		return range(1, min(n, m_in + 1) + 1)
	return range(0)

########################################################################################################################################################
# the implementation of K_out for (j, k) in {1, ..., N}^2, j != k, m in {1, ..., m_out - 1}
# i.e. the set of steps of exit maneuver of a vehicle in spot k 
# s.t. another vehicle cannot enter the facility and be assigned to spot j without being delayed
# which returns True if a vehicle should not enter and be assigned to j when a vehicle in spot k completed m steps of its exit maneuver
# i.e. True if m in K_out(j, k) and False if m not in K_out(j, k) 
def in_Kout(j, k, m):

	if (k == j):
		if (spot2blk(j) < m_out):
			return (0 <= m < m_out - spot2blk(j))
		return False

	if side == 'double':
		J = doub2sing(j)
		K = doub2sing(k)
	else:
		assert side == 'single'
		J = j
		K = k

	if (angle == 0) and (mode == 'long'):
		if (J >= K) and (K < m_out):
			return (1 <= m < m_out - K)
		return False

	if (spot2blk(j) >= spot2blk(k) - 1) and (spot2blk(k) < m_out):
		return (1 <= m < m_out - spot2blk(k))
	return False

########################################################################################################################################################
# the alternative set implementation of K_out for (j, k) in {1, ..., N}^2, j != k, m in {1, ..., m_out - 1}
def range_Kout(j, k):

	if (k == j):
		if (spot2blk(j) < m_out):
			return range(0, m_out - spot2blk(j))
		return False

	if side == 'double':
		J = doub2sing(j)
		K = doub2sing(k)
	else:
		assert side == 'single'
		J = j
		K = k

	if (angle == 0) and (mode == 'long'):
		if (J >= K) and (K < m_out):
			return range(1, m_out - K)
		return range(0)

	if (spot2blk(j) >= spot2blk(k) - 1) and (spot2blk(k) < m_out):
		return range(1, m_out - spot2blk(k))
	return range(0)

########################################################################################################################################################
# the implementation of K_in for (j, k) in {1, ..., N}^2, m in {1, ..., m_in - 2}
# i.e. the set of steps of the enter maneuver of a vehicle in spot k 
# s.t. a vehicle cannot enter the system and be assigned to spot j without being delayed
# which returns True if the vehicle assigned to j can enter the system if a vehicle in spot k completed m steps of its enter maneuver
# i.e. True if m in K_in(j, k) and False if m not in K_in(j, k) 
def in_Kin(j, k, m):
	 
	if (j == k):
		return False 

	if side == 'double':
		J = doub2sing(j)
		K = doub2sing(k)
	else:
		assert side == 'single'
		J = j
		K = k

	if (spot2blk(k) <= spot2blk(j)) and (spot2blk(k) < m_in):
		return (1 <= m < m_in - spot2blk(k) + 1)

	if (angle == 0 and mode == 'short') or angle == 90:
		if (spot2blk(k) == spot2blk(j) + 1) and (spot2blk(k) < m_in):
			return (1 <= m < m_in - spot2blk(k) + 1)

	return False

########################################################################################################################################################
# the alternative set implementation of K_in for (j, k) in {1, ..., N}^2, m in {1, ..., m_in - 2}
def range_Kin(j, k):

	if (j == k):
		return range(0) 

	if (spot2blk(k) <= spot2blk(j)) and (spot2blk(k) < m_in):
		return range(1, m_in - spot2blk(k) + 1)

	if (angle == 0 and mode == 'short') or angle == 90:
		if (spot2blk(k) == spot2blk(j) + 1) and (spot2blk(k) < m_in):
			return range(1, m_in - spot2blk(k) + 1)

	return range(0)