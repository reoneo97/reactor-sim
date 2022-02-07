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
        self.init_logs()

    def run(self, time_interval: float = 0.1, time_end: float = 300):
        # Setup an array to log concentration at every interval

        ts = np.arange(0, time_end, time_interval)
        vol_n = int(self.vol//self.flow_rate) + 1 
        vol_intervals, vol_step = np.linspace(0,self.vol,vol_n,retstep=True)
        
        # Calculate starting amount in each volume_step
        pa_amt = vol_step/self.flow_rate*self.pa_feed
        ipa_amt = vol_step/self.flow_rate*self.ipa_feed
        ipp_amt = 0
        water_amt = 0

        # Log initial concentration as class variable for conversion calculation
        self.initial_pa_conc = pa_amt/vol_step

        
        # Logging the number of mol

        conc_profile = np.zeros(shape=(vol_n,4))
        conc_profile[:,0] = pa_amt/vol_step
        conc_profile[:,1] = ipa_amt/vol_step
        conc_profile[:,2] = 0
        conc_profile[:,3] = 0

        conc_profile = np.concatenate(
            [vol_intervals.reshape(-1, 1), conc_profile], axis=1)
        self.log(conc_profile)


        for t in tqdm(ts[1:]):
            new_conc_profile = np.zeros(shape=(vol_n,5))
            new_conc_profile[:,0] = vol_intervals
            
            # Concentrations from previous interval 
            pa_conc, ipa_conc, ipp_conc, water_conc = conc_profile[0,1:]
            rate = self.reaction.get_rate(
                ca=pa_conc, cb=ipa_conc, ce=ipp_conc, cw=water_conc,
                T=self.temperature
            )/60
            
            dn = rate*time_interval*vol_step
            pa_in = self.pa_feed*time_interval
            ipa_in = self.ipa_feed*time_interval
            
            pa_amt = pa_in + pa_conc*vol_step -pa_conc*time_interval*self.flow_rate - dn 
            ipa_amt = ipa_in + ipa_conc*vol_step- ipa_conc*time_interval*self.flow_rate - dn
            ipp_amt = ipp_conc*vol_step- ipp_conc*time_interval*self.flow_rate + dn
            water_amt = water_conc*vol_step- water_conc*time_interval*self.flow_rate + dn

            new_conc_profile[0,1] = pa_amt/vol_step
            new_conc_profile[0,2] = ipa_amt/vol_step
            new_conc_profile[0,3] = ipp_amt/vol_step
            new_conc_profile[0,4] = water_amt/vol_step

            for v in range(1,vol_n):
                #Concentration from prev time interval
                pa_conc, ipa_conc, ipp_conc, water_conc = conc_profile[v, 1:]
                rate = self.reaction.get_rate(
                    ca=pa_conc, cb=ipa_conc, ce=ipp_conc, cw=water_conc,
                    T=self.temperature
                )/60

                dn = rate*time_interval*vol_step
                
                #Inflow from Upstream
                pa_in = new_conc_profile[v-1,1]*time_interval*self.flow_rate
                ipa_in = new_conc_profile[v-1,2]*time_interval*self.flow_rate
                ipp_in = new_conc_profile[v-1,3]*time_interval*self.flow_rate
                water_in = new_conc_profile[v-1, 4]*time_interval*self.flow_rate
                #Mass Balance
                pa_amt = pa_in + pa_conc*vol_step -pa_conc*time_interval*self.flow_rate - dn 
                ipa_amt = ipa_in + ipa_conc*vol_step- ipa_conc*time_interval*self.flow_rate - dn
                ipp_amt = ipp_in + ipp_conc*vol_step- ipp_conc*time_interval*self.flow_rate + dn
                water_amt = water_in + water_conc*vol_step- water_conc*time_interval*self.flow_rate + dn

                #Saving Concentration Values
                new_conc_profile[v, 1] = pa_amt/vol_step
                new_conc_profile[v, 2] = ipa_amt/vol_step
                new_conc_profile[v, 3] = ipp_amt/vol_step
                new_conc_profile[v, 4] = water_amt/vol_step
            self.log(new_conc_profile)
            conc_profile = new_conc_profile
        pd.DataFrame(self.logs[1]).to_csv("initial.csv")

        pd.DataFrame(new_conc_profile).to_csv("final.csv")

        data = np.array(self.logs)
        print(data.shape)

    def init_logs(self):
        self.logs = []
    
    def log(self,conc_profile:np.array):
        self.logs.append(conc_profile.copy())
