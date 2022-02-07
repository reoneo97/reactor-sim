import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

from .. import const
from ..reaction import ReactionModel
from ..reactor import Reactor


class IdealPFR(Reactor):
    def __init__(self,
                 pa_feed: float, M: int, vol: float, temperature: float,
                 rxn_model: ReactionModel):
        # Initial Setup in base model
        super().__init__(pa_feed, M, vol, temperature, rxn_model)

    def run(self, time_interval: float = 0.1, time_end: float = 300):
        # Setup an array to log concentration at every interval

        ts = np.arange(0, time_end, time_interval)
        vol_n = int(self.vol//self.flow_rate) + 1 
        vol_intervals, vol_step = np.linspace(0,self.vol,vol_n,retstep=True)
        vol_intervals = vol_intervals.reshape(-1,1)
        
        # Calculate starting amount in each volume_step
        pa_amt = vol_step/self.flow_rate*self.pa_feed
        ipa_amt = vol_step/self.flow_rate*self.ipa_feed
        ipp_amt = 0
        water_amt = 0

        # Log initial concentration as class variable for conversion calculation
        self.initial_pa_conc = pa_amt/vol_step

        
        # Logging the number of mol
        mol_profile = np.zeros(shape=(vol_n,4))
        mol_profile[:,0] = pa_amt
        mol_profile[:,1] = ipa_amt
        mol_profile[:,2] = ipp_amt
        mol_profile[:,3] = water_amt

        mol_profile = np.concatenate([vol_intervals,mol_profile],axis=1)


        for t in ts[1:]:
            new_mol_profile = np.zeros(shape=(vol_n,6))
            new_mol_profile[:,0] = vol_intervals
            vol_added = time_interval* self.flow_rate
            print(mol_profile[0,1:])
            pa_conc, ipa_conc, ipp_conc, water_conc = [i/self.vol_step for i in mol_profile[0,1:]]
            break
    

    def init_logs(self,ts):
        self.logs = []
    
    def log(self,conc_profile:np.array):
        self.logs.append(conc_profile)
