import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

from .. import const
from ..reaction import ReactionModel
from ..reactor import Reactor


class IdealCSTR(Reactor):
    def __init__(self,
                 pa_feed: float, M: int, vol: float, temperature: float,
                 rxn_model: ReactionModel):
        # Initial Setup in base model
        super().__init__(pa_feed,M, vol,temperature, rxn_model)


    def run(self, time_interval: float = 0.1, time_end: float = 300.):
        """Performs the simulation of the CSTR 

        Args:
            time_interval (float, optional): [description]. Defaults to 0.1.
            time_end (float, optional): [description]. Defaults to 300..
        """
        self.reset_log()
        ts = np.arange(time_interval, time_end, time_interval)
        # Initial Concentration Values
        pa_amt = self.vol/self.flow_rate*self.pa_feed
        ipa_amt = self.vol/self.flow_rate*self.ipa_feed
        ipp_amt = 0
        water_amt = 0

        pa_conc = pa_amt/self.vol
        # Log initial concentration as class variable for conversion calculation
        self.initial_pa_conc = pa_conc
        ipa_conc = ipa_amt/self.vol
        ipp_conc = 0
        water_conc = 0

        self.log(0, pa_conc, ipa_conc, ipp_conc, water_conc)
        for t in tqdm(ts):

            rate = self.reaction.get_rate(
                ca = pa_conc, cb = ipa_conc, ce= ipp_conc, cw = water_conc,
                T=self.temperature
            )/60  # Rate converted from h-1 to min-1
            dc = rate*time_interval
            dn = dc*self.vol

            pa_amt += (self.pa_feed - dn)
            ipa_amt += (self.ipa_feed - dn)
            water_amt += dn
            ipp_amt += dn

            amts = [pa_amt, ipa_amt, water_amt, ipp_amt]
            pa_conc, ipa_conc, ipp_conc, water_conc = [
                i/self.vol for i in amts]

            pa_amt -= self.flow_rate*pa_conc
            ipa_amt -= self.flow_rate*ipa_conc
            ipp_amt -= self.flow_rate*ipp_conc
            water_amt -= self.flow_rate*water_conc

            self.log(t, pa_conc, ipa_conc, ipp_conc, water_conc)

        self.plot_logs()

    def log(self, t, pa_conc, ipa_conc, ipp_conc, water_conc):
        """
        Log Concentration values within the CSTR and calculate the current
        conversion of the reactor 
        """
        conversion = 1-(pa_conc/self.initial_pa_conc)
        self.logs.append(
            [t, pa_conc, ipa_conc, ipp_conc, water_conc, conversion])

    def reset_log(self):
        self.logs = []

    def plot_logs(self):
        col_names = [
            "Time", "Palmitic Acid Concentration", "IPA Concentration",
            "IPP Concentration", "Water Concentration", "Conversion"
        ]
        df = pd.DataFrame(self.logs, columns=col_names)

        df.plot.line(x="Time", y="Conversion")
        plt.show()
