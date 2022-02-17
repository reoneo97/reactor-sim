import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from loguru import logger

from .. import const
from ..reaction import ReactionModel, zna_reaction
from ..reactor import Reactor


class IdealPBR(Reactor):
    def __init__(self,
                 pa_feed: float, M: int, vol: float, weight: float, temperature: float,
                 rxn_model: ReactionModel = zna_reaction, rate_test=1):

        # Initial Setup in base model
        super().__init__(pa_feed, M, temperature, vol, rxn_model=rxn_model)
        # Weight of Catalyst in Grams
        self.weight = weight
        self.rate_test = rate_test
        self.init_logs()

    def run(self, time_interval: float = 0.1, time_end: float = 300):
        # Setup an array to log concentration at every interval

        self.ts = np.arange(0, time_end, time_interval)
        #! Setting number of volume intervals in PBR to a fixed value
        vol_n = 100
        weight_per_vol = self.weight/vol_n
        print("Weight per Vol", weight_per_vol)
        self.vol_intervals, vol_step = np.linspace(
            0, self.vol, vol_n, retstep=True)

        # Calculate starting amount in each volume_step
        pa_amt = vol_step/self.flow_rate*self.pa_feed
        ipa_amt = vol_step/self.flow_rate*self.ipa_feed
        ipp_amt = 0
        water_amt = 0

        # Log initial concentration as class variable for conversion calculation
        self.initial_pa_conc = pa_amt/vol_step

        # Logging the number of mol

        conc_profile = np.zeros(shape=(vol_n, 4))
        conc_profile[:, 0] = pa_amt/vol_step
        conc_profile[:, 1] = ipa_amt/vol_step
        conc_profile[:, 2] = 0
        conc_profile[:, 3] = 0

        conc_profile = np.concatenate(
            [self.vol_intervals.reshape(-1, 1), conc_profile], axis=1)

        self.log(conc_profile)

        for t in self.ts[1:]:
            new_conc_profile = np.zeros(shape=(vol_n, 5))
            new_conc_profile[:, 0] = self.vol_intervals

            # Concentrations from previous interval
            pa_conc, ipa_conc, ipp_conc, water_conc = conc_profile[0, 1:]
            rate = self.reaction.get_rate(
                ca=pa_conc, cb=ipa_conc, ce=ipp_conc, cw=water_conc,
                T=self.temperature
            )*60*self.rate_test  # ! To Remove this rate_test after testing

            # Amount Reacted dN/dt = -r_a' x W => dn  \approx -r_a' x W \Delta_t
            dn = rate*time_interval*weight_per_vol
            pa_in = self.pa_feed*time_interval
            ipa_in = self.ipa_feed*time_interval

            pa_amt = pa_in + pa_conc*vol_step - pa_conc*time_interval*self.flow_rate - dn
            ipa_amt = ipa_in + ipa_conc*vol_step - \
                ipa_conc*time_interval*self.flow_rate - dn
            ipp_amt = ipp_conc*vol_step - ipp_conc*time_interval*self.flow_rate + dn
            water_amt = water_conc*vol_step - water_conc*time_interval*self.flow_rate + dn

            new_conc_profile[0, 1] = pa_amt/vol_step
            new_conc_profile[0, 2] = ipa_amt/vol_step
            new_conc_profile[0, 3] = ipp_amt/vol_step
            new_conc_profile[0, 4] = water_amt/vol_step

            for v in range(1, vol_n):
                # Concentration from prev time interval
                pa_conc, ipa_conc, ipp_conc, water_conc = conc_profile[v, 1:]
                rate = self.reaction.get_rate(
                    ca=pa_conc, cb=ipa_conc, ce=ipp_conc, cw=water_conc,
                    T=self.temperature
                )*60*self.rate_test  # ! To Remove this rate_test after testing

                dn = rate*time_interval*vol_step

                # Inflow from Upstream
                pa_in = new_conc_profile[v-1, 1]*time_interval*self.flow_rate
                ipa_in = new_conc_profile[v-1, 2]*time_interval*self.flow_rate
                ipp_in = new_conc_profile[v-1, 3]*time_interval*self.flow_rate
                water_in = new_conc_profile[v-1, 4] * \
                    time_interval*self.flow_rate
                # Mass Balance
                pa_amt = pa_in + pa_conc*vol_step - pa_conc*time_interval*self.flow_rate - dn
                ipa_amt = ipa_in + ipa_conc*vol_step - \
                    ipa_conc*time_interval*self.flow_rate - dn
                ipp_amt = ipp_in + ipp_conc*vol_step - \
                    ipp_conc*time_interval*self.flow_rate + dn
                water_amt = water_in + water_conc*vol_step - \
                    water_conc*time_interval*self.flow_rate + dn

                # Saving Concentration Values
                new_conc_profile[v, 1] = pa_amt/vol_step
                new_conc_profile[v, 2] = ipa_amt/vol_step
                new_conc_profile[v, 3] = ipp_amt/vol_step
                new_conc_profile[v, 4] = water_amt/vol_step
            self.log(new_conc_profile)
            conc_profile = new_conc_profile

        data = np.array(self.logs)
        logger.info(data.shape)

    def plot_logs(self):
        col_names = [
            "Time", "Palmitic Acid Concentration", "IPA Concentration",
            "IPP Concentration", "Water Concentration", "Conversion"
        ]
        df = pd.DataFrame(self.logs, columns=col_names)

        df.plot.line(x="Time", y="Conversion")
        plt.show()

    def init_logs(self):
        self.logs = []

    def log(self, conc_profile: np.array):
        # print(conc_profile.shape)
        pa_conc = conc_profile[:, 1]
        assert conc_profile.min() >= 0, "There should not be any negative" \
            f" concentrations. Current Length {len(self.logs)}"
        conversion = 1 - pa_conc/self.initial_pa_conc
        conversion = conversion.reshape(-1, 1)
        # print(conc_profile.shape, conversion.shape)

        new_data = np.concatenate([conc_profile, conversion], axis=1)
        self.logs.append(new_data.copy())

    def get_conversion(self):
        data = np.array(self.logs)
        print(data.shape)

        conversion_data = data[:, :, 5]
        return conversion_data

    def get_volume_intervals(self):
        return self.vol_intervals
