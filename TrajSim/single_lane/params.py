########################################################################################################################################################
# Last checked: Feb 06, 2022 by Xinyu Liu
# This file specifies the input data to all other files required for trajectory-based simulation, with:
# Configuration -- Single Lane; 
# Access Control -- Partial 
# Input Distributions -- Exponential service times, and deterministic movement times 

########################################################################################################################################################
import math
from math import ceil
import random


# dirname = '../'
dirname = ''
suffix = ''

''' full vs prtial access control '''
control = 'partial' 
# control = 'full'

''' single vs double sided '''
side = 'single'
# side = 'double'

''' 0-degree vs 90-degree '''
angle = 0
# angle = 90

''' long vs short spots vs single lane '''
# mode = 'long'
# mode = 'short'
mode = 'singlelane'

''' mean service time '''
meanSERV = 90. # sec

''' average (desired) speed '''
rateDRIV = 44./3. # ft/sec

''' types of input distributions '''
# simType = 'det'
# simType = 'cav'
simType = 'exp'
# simType = 'unif'

''' simulation length and warmup period length (if any) '''
WARMUP_HOUR = 2 # hr
SIM_HOUR = 20 # hr
SIM_ITER = 20

''' max number of spots on each side to evaluated '''
if side == 'single':
	LOT_TOTAL = 40
else:
	assert side == 'double'
	LOT_TOTAL = 20
if mode == 'singlelane':
	LOT_TOTAL = 60

########################################################################################################################################################
# The following section of the file should not be changed unless you know exactly what you are trying to do.
# WARNING: Unintended changes to this section can have unexpected effects to the simulation results.

debug = False
spmatch = False
delay = False
free_curb = False
if free_curb:
	suffix = 'free_curb'
VEHICLE_IDX = None
ptype = 0

WARMUP_UNIT = 3600 * WARMUP_HOUR # sec
SIM_UNIT = 3600 * SIM_HOUR # sec

if angle == 90:
	# 90-degree configurations
	nUnit = 3
	dgap = 2
	meanPOUT = 4.6 # sec
	meanPLIN = 9.7 # sec
	LOT_LENGTH = 10. # ft
	LOT_WIDTH = 17.75 # ft
	LANE_WIDTH = 23. # ft

if angle == 0:
	# 0-degree configurations
	nUnit = 1
	dgap = 1
	meanPOUT = 5.6 # 5.6 sec
	LOT_WIDTH = 10. # ft
	LANE_WIDTH = 12. # ft
	if mode == 'long':
		meanPLIN = 4.8 # sec
		LOT_LENGTH = 32.5 # ft
	if mode == 'short':
		meanPLIN = 16.5 # sec
		LOT_LENGTH = 25.0 # ft
	if mode == 'singlelane':
		meanPLIN = 1.
		LOT_LENGTH = 23.0 # ft

if angle == 45:
	# 45-degree configurations
	nUnit = 2 # ?
	dgap = 2 # ?
	meanPOUT = 10.3 # sec
	meanPLIN = 4.2 # sec
	LOT_LENGTH = 13.0 # ft
	LOT_WIDTH = 17.5 # ft
	LANE_WIDTH = 12.0 # ft

CAR_LENGTH = LOT_LENGTH * nUnit # ft
meanDRIV = CAR_LENGTH / rateDRIV # time to drive through one boarding spot in sec
if (angle == 0 and mode == 'short') or (angle == 90 and spmatch):
	meanPLIN += meanDRIV
m_out = int( round(meanPOUT / meanDRIV) )
m_in = int( round(meanPLIN / meanDRIV) ) 
g_out = max(m_out, 3)

if spmatch:
	pass
meanPOUT = m_out * meanDRIV
meanPLIN = m_in * meanDRIV

rateSERV = 1. / meanSERV
ratePOUT = 1. / meanPOUT
ratePLIN = 1. / meanPLIN

if control == 'partial' and (not (angle == 0 and mode == 'long')):
	pcontrol = control + str(ptype)
else:
	pcontrol = control

if angle == 0:
	filename = '%s_%s_%s_%s_%s_%s_%s'%(angle, mode, side, pcontrol, meanSERV, simType, suffix)
elif angle == 90:
	filename ='%s_%s_%s_%s_%s_%s'%(angle, side, pcontrol, meanSERV, simType, suffix)

SEED_SERV = random.Random(2020).sample(range( int(1e12) ), LOT_TOTAL)
SEED_POUT = random.Random(9031).sample(range( int(1e12) ), LOT_TOTAL)
SEED_PLIN = random.Random(4654).sample(range( int(1e12) ), LOT_TOTAL)
SEED_DRIV = random.Random(1203).sample(range( int(1e12) ), LOT_TOTAL)

SMALL_INTERVAL = 1e-7

event_priority = {'leave_system':        1,
				  'finish_pulling_out':  2,
				  'prepare_pulling_out': 3,
				  'start_service':       4,
				  'start_second_enter':  5,
				  'prepare_first_exit':  6,
				  'start_pulling_in':    7,
				  'enter_system':        8,
				  'add_stop_idx':        7.5
				  }
