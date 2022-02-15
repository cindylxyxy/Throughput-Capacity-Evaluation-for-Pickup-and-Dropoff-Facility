########################################################################################################################################################
# Last updated: Feb 06, 2022 by Xinyu Liu
# This file specifies the input data to all other files required for Markov Decision Process model and simulation, with:
# Configuration -- Single-Sided 90-degree; 
# Input Distributions -- Exponential service times, and deterministic movement times 

########################################################################################################################################################
import math
from math import ceil
import random

# dirname = '../'
dirname = ''
suffix = ''

''' single vs double sided '''
side = 'single'
# side = 'double'

''' 0-degree vs 90-degree '''
# angle = 0
angle = 90

''' only valid for 0-degree '''
''' long vs short spots '''
# mode = 'long'
mode = 'short'

''' mean service time '''
meanSERV = 60. # sec

''' average (desired) speed '''
avgSPEED = 44./3. # ft/sec

''' max number of spots on each side to evaluated '''
if side == 'single':
	LOT_TOTAL = 6
else:
	assert side == 'double'
	LOT_TOTAL = 3

''' discount rate in value iteration for MDP '''
DISC_RATE = 0.005 

########################################################################################################################################################
# The following section of the file should not be changed unless you know exactly what you are trying to do.
# WARNING: Unintended changes to this section can have unexpected effects to the numerical results.

control = 'mdp'

SIM_UNIT = 3600 * SIM_HOUR # sec
SIM_HOUR = 20 # hr
SIM_ITER = 20

if angle == 90:
	# 90-degree configurations
	nUnit = 3
	dgap = 2
	g_in = 3
	meanPOUT = 4.6 # sec
	meanPLIN = 9.7 # sec
	LOT_LENGTH = 10. # ft

if angle == 0:
	# 0-degree configurations
	nUnit = 1
	dgap = 1
	meanPOUT = 5.6 # sec
	if mode == 'long':
		g_in = 1
		meanPLIN = 4.8 # sec
		LOT_LENGTH = 32.5 # ft
	if mode == 'short':
		g_in = 3
		meanPLIN = 16.5 # sec
		LOT_LENGTH = 25.0 # ft

meanDRIV = (LOT_LENGTH * nUnit) / avgSPEED
m_out = int( round(meanPOUT / meanDRIV) ) + 1
m_in = int( round(meanPLIN / meanDRIV) )
g_out = max(m_out, 3)
meanPOUT = m_out * meanDRIV
meanPLIN = m_in * meanDRIV
rateSERV = 1.0 / meanSERV
rateDRIV = 1.0 / meanDRIV 

if angle == 0:
	filename = '%s_%s_%s_%s_%s'%(angle, mode, side, meanSERV, suffix)
elif angle == 90:
	filename ='%s_%s_%s_%s'%(angle, side, meanSERV, suffix)

SEED_SERV = random.Random(2020).sample(range( int(1e12) ), LOT_TOTAL)
SEED_POUT = random.Random(9031).sample(range( int(1e12) ), LOT_TOTAL)
SEED_PLIN = random.Random(4654).sample(range( int(1e12) ), LOT_TOTAL)
SEED_DRIV = random.Random(1203).sample(range( int(1e12) ), LOT_TOTAL)
