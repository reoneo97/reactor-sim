from simulator.reaction import ptsa_reaction, ReactionModel
import simulator.const as const
from simulator.reactor import IdealCSTR, IdealPFR
from loguru import logger




temp = 373
pa_feed = 900
M = 5
vol = 500000  # Units of m3
rxn = ptsa_reaction

time_interval = 0.1  # Time Intervals to perform simulation
time_end = 300

reactor_model = IdealCSTR(pa_feed, M, vol, temp, rxn)
reactor_model.run(time_interval, time_end)