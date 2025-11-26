########################################################################################################################################################
# Last updated: Oct 10, 2025 by Xinyu Liu
# This file specifies the input data to all other files required for Markov decision process / Markov chain models and simulations, with:
# Configuration -- Single-Sided 0-degree with long spots; 
# Input Distributions -- Exponential service times, and deterministic movement times 

########################################################################################################################################################
import math
from math import ceil
import random

''' create a directory based on today's date '''
from datetime import date
import os
today = date.today()
if not os.path.exists('%sp/'%today):
	os.makedirs('%sp/'%today)
outdir = '%sp/'%today
outdir = '2025-11-21p/'

debug = False
suffix = ''
spmatch = True

# side = 'single'
side = 'double'

angle = 0
# angle = 90

mode = 'long'
# mode = 'short'

meanSERV = 60.
rateSERV = 1./ 60. # average service time is 60. sec
avgSPEED = 44./ 3.
rateTRVL = 3./ 44. # average speed is 44./3. ft/sec

simType = 'mc'
# simType = 'cav'

SIM_HOUR = 20 # hr
SIM_ITER = 20
DISC_RATE = 0.005 

if spmatch:
	suffix += 'spmatch'
	# assert simType == 'cav'
	SIM_HOUR = 100 # hr
	SIM_ITER = 20

SIM_UNIT = 3600 * SIM_HOUR # sec
if side == 'single':
	LOT_TOTAL = 40
else:
	assert side == 'double'
	LOT_TOTAL = 20

########################################################################################################################################################
# The following section of the file should not be changed unless you know exactly what you are trying to do.
# WARNING: Unintended changes to this section can have unexpected effects to the simulation results.

if angle == 90:
	# 90-degree configurations
	in_direction = 'backward'
	out_direction = 'forward'
	nUnit = 3
	meanPOUT = 4.6 # sec
	meanPLIN = 9.7 # sec
	LOT_LENGTH = 10. # ft

if angle == 0:
	# 0-degree configurations
	nUnit = 1
	meanPOUT = 5.6 # sec
	if mode == 'long':
		in_direction = 'forward'
		out_direction = 'forward'
		meanPLIN = 4.8 # sec
		LOT_LENGTH = 32.5 # ft
	if mode == 'short':
		in_direction = 'backward'
		out_direction = 'forward'
		meanPLIN = 16.5 # sec
		LOT_LENGTH = 25.0 # ft

meanDRIV = (LOT_LENGTH * nUnit) / avgSPEED
MOUT = int( round(meanPOUT / meanDRIV) ) + 1
MIIN = int( round(meanPLIN / meanDRIV) )
meanPOUT = MOUT * meanDRIV
meanPLIN = MIIN * meanDRIV
rateSERV = 1.0 / meanSERV
rateDRIV = 1.0 / meanDRIV 

if angle == 0:
	filename = '%s_%s_%s_%s_%s_%s'%(angle, mode, side, meanSERV, simType, suffix)
elif angle == 90:
	filename ='%s_%s_%s_%s_%s'%(angle, side, meanSERV, simType, suffix)

SEED_SERV = random.Random(2020).sample(range( int(1e12) ), LOT_TOTAL)
SEED_POUT = random.Random(9031).sample(range( int(1e12) ), LOT_TOTAL)
SEED_PLIN = random.Random(4654).sample(range( int(1e12) ), LOT_TOTAL)
SEED_DRIV = random.Random(1203).sample(range( int(1e12) ), LOT_TOTAL)
