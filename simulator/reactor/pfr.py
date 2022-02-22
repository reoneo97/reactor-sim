import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from tqdm import tqdm
from loguru import logger
from typing import Tuple, List

from .. import const
from ..reaction import ReactionModel, ptsa_reaction
from ..reactor import Reactor
from ..energy import mh_cp, MixedHeatCapacity


class IdealPFR(Reactor):
    def __init__(self,
                 pa_feed: float, M: int, vol: float, temperature: float,
                 rxn_model: ReactionModel = ptsa_reaction):
        # Initial Setup in base model
        super().__init__(pa_feed, M, temperature, vol, rxn_model)
        self.init_logs()

    def run(self, time_interval: float = 0.1, time_end: float = 300):
        # Setup an array to log concentration at every interval

        ts = np.arange(0, time_end, time_interval)
        vol_n = int(self.vol//self.flow_rate) + 1
        self.vol_intervals, vol_step = np.linspace(
            0, self.vol, vol_n, retstep=True)

        # Calculate starting amount in each volume_step
        pa_amt = vol_step/self.flow_rate*self.pa_feed
        ipa_amt = vol_step/self.flow_rate*self.ipa_feed
        ipp_amt = 0
        water_amt = 0

        # Log initial concentration as class variable for conversion calculation
        self.initial_pa_conc = pa_amt/vol_step
        print(self.initial_pa_conc, self.initial_pa_conc*self.M)

        # Logging the number of mol

        conc_profile = np.zeros(shape=(vol_n, 4))
        conc_profile[:, 0] = pa_amt/vol_step
        conc_profile[:, 1] = ipa_amt/vol_step
        conc_profile[:, 2] = 0
        conc_profile[:, 3] = 0

        conc_profile = np.concatenate(
            [self.vol_intervals.reshape(-1, 1), conc_profile], axis=1)
        self.log(conc_profile)

        for t in tqdm(ts[1:]):
            new_conc_profile = np.zeros(shape=(vol_n, 5))
            new_conc_profile[:, 0] = self.vol_intervals

            # Initial V_0
            pa_conc, ipa_conc, ipp_conc, water_conc = conc_profile[0, 1:]
            rate = self.reaction.get_rate(
                ca=pa_conc, cb=ipa_conc, ce=ipp_conc, cw=water_conc,
                T=self.temperature
            )*60

            dn = rate*time_interval*vol_step
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
                )*60  # Simulation in min but rate equation in h^-1

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

    def init_logs(self):
        self.logs = []

    def log(self, conc_profile: np.array):
        # print(conc_profile.shape)
        pa_conc = conc_profile[:, 1]
        conversion = 1 - pa_conc/self.initial_pa_conc
        conversion = conversion.reshape(-1, 1)
        # print(conc_profile.shape, conversion.shape)

        new_data = np.concatenate([conc_profile, conversion], axis=1)
        self.logs.append(new_data.copy())

    def get_conversion(self):
        data = np.array(self.logs)

        conversion_data = data[:, :, 5]
        return conversion_data

    def get_volume_intervals(self):
        return self.vol_intervals


class RealPFR(Reactor):
    def __init__(self,
                 pa_feed: float, M: int,
                 L: float, R: float, feed_temp: float, cooler_temp: float,
                 space_interval: float = 100, cp_model: MixedHeatCapacity = mh_cp,
                 rxn_model: ReactionModel = ptsa_reaction,
                 **kwargs
                 ):
        self.pa_feed = pa_feed
        self.M = M
        self.ipa_feed = pa_feed*M
        self.L = L
        self.R = R
        self.feed_temp = feed_temp
        self.cooler_temp = cooler_temp
        self.cross_area = math.pi*R*R
        self.heat_transfer_area = 2*math.pi*R*L
        self.vol = self.cross_area*L
        self.logs = []
        self.space_interval = space_interval
        self.cp_model = cp_model
        super().init_flow_rate()

    def log(self, info):
        """
        Data is stored in a T x L x R x 7 Matrix
            (Time x Radius x Length  x Conc+Conversion + Temp)
        """
        self.logs.append(info)

    def cylinder_area(self, r: float, dr: float, dz: float) -> Tuple[float, float]:
        """
        Helper function to calculate surface area of cylindrical rectangle for 
        heat transfer calculations

        Returns:
        in_area: Heat transfer area for incoming conduction
        out_area: Heat transfer area for outgoing conduction
        """
        return 2*math.pi*(r)*dz, 2*math.pi*(r+dr)*dz

    def cylinder_vol(self, r, dr, dz):
        """
        Helper function to calculate volume of cylindrical rectangle
        Formula = \pi*((r+dr)^2 -r^2)*dl
        """

        return 2*math.pi*r*dr*dz

    def run(self, time_interval: float = 0.1, time_end: float = 300):
        ts = np.arange(0, time_end, time_interval)
        l_step = self.L/self.space_interval
        ls = np.arange(0, self.L, self.space_interval)
        r_step = self.R/self.space_interval
        r_s = np.arange(0, self.L, self.space_interval)

        # Save volumes of cylindrical rectangles for faster calculation
        vol_r_intervals = [self.cylinder_vol(r, r_step, l_step) for r in r_s]
        # Save volumes of flow_rates
        norm_vol = [i/sum(vol_r_intervals) for i in vol_r_intervals]
        flow_rate_intervals = [i*self.flow_rate for i in norm_vol]

        pa_conc_init = self.pa_feed/self.flow_rate
        ipa_conc_init = self.ipa_feed/self.flow_rate

        logger.info(pa_conc_init, ipa_conc_init)
        water_conc_init = 0
        ipp_conc_init = 0
        la_conc_init = 0
        temp = self.feed_temp
        # Units of kmol/m3 -> mol/dm3
        # pa_amt = vol_step/self.flow_rate*self.pa_feed
        # ipa_amt = vol_step/self.flow_rate*self.ipa_feed
        # ipp_amt = 0
        # water_amt = 0

    def calculate_cp(self, amts: List[float], T: float):
        cp = self.cp_model(T, amts)
        return cp
