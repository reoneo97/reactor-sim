from mimetypes import init
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from tqdm import tqdm
from loguru import logger
from typing import Tuple, List, Iterable

from .. import const
from ..reaction import ReactionModel, ptsa_reaction
from ..reactor import Reactor
from ..energy import mh_cp, MixedHeatCapacity

UPPER_BOUND = 1e5


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
                 L: float, R: float, feed_temp: float, heater_temp: float,
                 space_interval: float = 26, cp_model: MixedHeatCapacity = mh_cp,
                 rxn_model: ReactionModel = ptsa_reaction,
                 **kwargs
                 ):
        self.pa_feed = pa_feed
        self.M = M
        self.ipa_feed = pa_feed*M
        self.L = L
        self.R = R
        self.D = R*2
        self.feed_temp = feed_temp
        self.heater_temp = heater_temp
        self.cross_area = math.pi*R*R
        self.heat_transfer_area = 2*math.pi*R*L
        self.vol = self.cross_area*L
        self.logs = []
        self.space_interval = space_interval
        self.cp_model = cp_model
        self.reaction = rxn_model
        super().init_flow_rate()
        self.velocity = self.flow_rate/self.cross_area
        self.get_pressure_drop()

    def log(self, info):
        """
        Data is stored in a T x L x R x 7 Matrix
            (Time x Length x Radius x (Temp + (5 Conc) + Conversion) )
        """
        assert info.shape == (self.space_interval, self.space_interval, 6)
        pa_conc = info[:, :, 1]
        conversion = 1 - pa_conc/self.initial_pa_conc

        new_data = np.concatenate(
            [info, np.expand_dims(conversion, 2)], axis=2)
        self.logs.append(new_data)

    def cylinder_area_r(self, r: float, dr: float, dz: float) -> Tuple[float, float]:
        """
        Helper function to calculate surface area of cylindrical rectangle for
        heat transfer calculations in the r-direction

        Returns:
        in_area: Heat transfer area for incoming conduction
        out_area: Heat transfer area for outgoing conduction
        """
        return 2*math.pi*(r)*dz, 2*math.pi*(r+dr)*dz

    def cylinder_area_z(self, r: float, dr: float) -> Tuple[float, float]:
        """
        Helper function to calculate surface area of cylindrical rectangle for
        heat transfer calculations in the z-direction

        Returns:
        area: Heat transfer area for convection + conduction in z-direction
        """
        return 2*math.pi*r*dr

    def cylinder_vol(self, r, dr, dz):
        """
        Helper function to calculate volume of cylindrical rectangle
        Formula = \pi*((r+dr)^2 -r^2)*dl
        """

        return 2*math.pi*r*dr*dz

    def run(self, time_interval: float = 0.05, time_end: float = 300):
        ts = np.arange(0, time_end, time_interval)

        dz = self.L/self.space_interval
        l_s = np.linspace(0, self.L, self.space_interval)
        dr = self.R/self.space_interval
        r_s = np.linspace(dr, self.R, self.space_interval)

        # Save volumes of cylindrical rectangles for faster calculation
        vol_r_intervals = [self.cylinder_vol(r, dr, dz) for r in r_s]
        area_z_intervals = [self.cylinder_area_z(r, dr) for r in r_s]
        norm_area = [i/sum(area_z_intervals) for i in area_z_intervals]
        # Save volumes of flow_rates
        norm_vol = [i/sum(vol_r_intervals) for i in vol_r_intervals]
        flow_rate_intervals = [i*self.flow_rate for i in norm_area]

        # Initial State of PFR at t=0
        feed_temp = self.feed_temp
        pa_conc_init = self.initial_pa_conc = self.pa_feed/self.flow_rate
        ipa_conc_init = self.ipa_feed/self.flow_rate
        # logger.info(f"Initial PA Concentration: {pa_conc_init:.3f}")
        # logger.info(f"Initial IPA Concentration: {ipa_conc_init:.3f}")
        water_conc_init = 0
        ipp_conc_init = 0
        la_conc_init = 0

        # Storing the amount of heat required from the heating fluid at every time interval
        self.heater_total = np.zeros(shape=(len(ts)))

        # Initial Data Slc is the feed conditions which will be used for the first radial slice
        init_data_slc = np.array([
            feed_temp, pa_conc_init, ipa_conc_init, ipp_conc_init, water_conc_init,
            la_conc_init
        ])
        init_data = np.tile(
            init_data_slc, [self.space_interval, self.space_interval, 1])
        self.log(init_data)
        # Save Previous Data point for calculations
        prev_data_t = init_data
        # Running Simulation from t=1 onwards
        for i_t, t in enumerate(tqdm(ts[1:])):
            new_data = np.zeros(
                shape=(self.space_interval, self.space_interval, 6))
            heater_heat_out = 0
            for i, l in enumerate(l_s):

                # Storing slice of concentrations/ temp from prev radial segment z_{i-1}
                if i == 0:
                    prev_data_l = init_data[0, :, :].copy()
                else:
                    prev_data_l = new_data[i-1, :, :].copy()

                for j, r in enumerate(r_s):
                    # Fetch Values for Calculation
                    slc_vol = vol_r_intervals[j]
                    slc_flow_rate = flow_rate_intervals[j]

                    slc_molar_flow = norm_area[j]*self.flow_rate

                    prev_data_slc_l = prev_data_l[j, :]
                    prev_data_slc_t = prev_data_t[i, j, :]
                    # logger.info(prev_data_slc_l)

                    in_area_r, out_area_r = self.cylinder_area_r(
                        r, dr, dz)
                    area_z = self.cylinder_area_z(r, dr)
                    # Obtaining parameters from same section but at t=i-1
                    init_temp, pa_conc, ipa_conc, ipp_conc, water_conc, la_conc = prev_data_slc_t
                    concs = [pa_conc, ipa_conc, ipp_conc, water_conc]
                    # try:
                    #     assert all([i >= 0 for i in concs]
                    #                ), "Concentrations should be positive"
                    # except:
                    #     logger.warning(concs)
                    #     raise AssertionError(
                    #         "Simulation Error - Negative Concentration")
                    # Mass Balance
                    try:
                        rate = self.reaction.get_rate(
                            ca=pa_conc, cb=ipa_conc, ce=ipp_conc, cw=water_conc,
                            T=init_temp
                        )*60
                        heat = self.reaction.get_heat()
                        dn = rate*slc_vol*time_interval
                        rxn_heat = heat*dn
                        assert abs(dn) < UPPER_BOUND, "Simulation Error"
                        # print(t, i, j, rxn_heat, init_temp)
                    except Exception as e:
                        print()
                        print(e)
                        print(t, i, j, init_temp)
                        quit()
                    # Rate converted from /s to /min

                    # v_o/V  (C_{z-1,t} - C_{z,t-1})
                    conv_term = slc_flow_rate/slc_vol * \
                        (prev_data_slc_l[1:6] - prev_data_slc_t[1:6])
                    dc = conv_term.copy()
                    # PA and IPA, stoichiometric coefficient of -1
                    dc[:2] -= rate
                    # IPP and Water, Stoichiometric cofficient of 1
                    dc[2:4] += rate
                    dcdt = dc*time_interval
                    # logger.debug(prev_data_slc_t[1:6])
                    # logger.info(dcdt)

                    # ?: Should be able to merge the below lines
                    c_t = prev_data_slc_t[1:6] + dcdt

                    try:
                        assert all([i >= 0 for i in c_t]
                                   ), "Concentrations should be positive"
                    except:

                        logger.error(c_t)
                        logger.warning(f"Error at T={t} {i} {j}")
                        logger.warning(f"Temp is {init_temp}")
                        logger.warning(slc_flow_rate/slc_vol)
                        logger.warning(prev_data_slc_l[1:6])
                        logger.warning(prev_data_slc_t[1:6])
                        logger.warning(conv_term)
                        logger.warning(rate)
                        logger.warning(dc)
                        logger.warning(dcdt)
                        raise AssertionError("Simulation Error")
                        # Mass Balance
                    new_data[i, j, 1:6] = c_t
                    # Energy Balance
                    amts = c_t*slc_vol
                    cp = self.calculate_cp(init_temp, amts)
                    # logger.debug(f"R= {r}, C_p = {cp}")

                    # Temperature Values
                    # logger.warning(f"{in_area_r}, {out_area_r}")
                    # Conduction Terms - Radial Direction
                    if j == 0:
                        cond_r_in = (feed_temp - init_temp) / \
                            dr * in_area_r * const.k_e
                    else:
                        cond_r_in = (
                            prev_data_l[j-1, 0] - init_temp)/dr * in_area_r * const.k_e

                    if j == self.space_interval-1:
                        # TODO: Find the Heat Transfer Coefficient for HEater TEmp
                        # TODO: Should remove the dr, since this is not a conduction and there is no KE
                        cond_r_out = (self.heater_temp -
                                      init_temp) * out_area_r * const.h_heater
                        heater_heat_out += cond_r_out
                    else:
                        cond_r_out = (
                            prev_data_t[i, j+1, 0] - init_temp)/dr * out_area_r * const.k_e

                    # Conduction Terms - Z direction - To First find dTdz
                    if i == 0:
                        dtdz = (prev_data_t[i+1, j, 0] - feed_temp)/dz
                    elif i == self.space_interval-1:
                        dtdz = (prev_data_t[i, j, 0] -
                                prev_data_t[i-1, j, 0])/dz
                    else:
                        dtdz = (prev_data_t[i+1, j, 0] -
                                prev_data_t[i-1, j, 0])/dz
                    cond_z = dtdz * area_z * const.k_e
                    # Convection Terms - Z Direction
                    dt = dtdz*dz
                    conv_net = slc_molar_flow * cp * dt
                    # print(cond_z, cond_r_in, cond_r_out, conv_net, rxn_heat)
                    # Cond r_out is positive and cond_r_in is also positive
                    dhdt = cond_z + cond_r_in + cond_r_out - conv_net + rxn_heat
                    # if t == ts[-1] and l == l_s[20]:
                    #     logger.info("Special Log")
                    #     logger.info(dhdt)
                    #     logger.info(conv_net)
                    #     logger.info(cond_z)

                    dh = dhdt*time_interval

                    temp_change = dh/cp
                    # logger.debug(f"Temp Change: {temp_change}")
                    new_data[i, j, 0] = init_temp + temp_change
                    if new_data[i, j, 0] <= 200:
                        logger.error(c_t)
                        logger.warning(
                            f"Error at T={t:3f} Z={l_s[i]:3f} R={r_s[j]:3f}")
                        logger.warning(prev_data_slc_l[0:6])
                        logger.warning(prev_data_slc_t[0:6])
                        logger.warning(
                            "cond_z | cond_r_in | cond_r_out | conv_net | rxn_heat")
                        logger.warning(
                            f"{cond_z:.3f}| {cond_r_in:.3f} | {cond_r_out:.3f} | {conv_net:.3f} | {rxn_heat:.3f}")
                        logger.warning(f"Temp Change: {temp_change}")
                        logger.error(
                            "Terminating Simulation - Temperature Negative")
                        self.log(new_data)
                        prev_data_t = new_data.copy()

                        self.all_data = np.array(self.logs)
                        self.calculate_output_conversion()
                        return self.all_data

            self.log(new_data)
            self.heater_total[i_t+1] = heater_heat_out
            prev_data_t = new_data.copy()

        self.all_data = np.array(self.logs)
        self.calculate_output_conversion()
        return self.all_data

    def calculate_cp(self, T: float, amts):
        # C_p values are in units of
        cp = self.cp_model(T, amts)
        return cp

    def get_dimensions(self,):
        dz = self.L/self.space_interval
        l_s = np.linspace(dz, self.L, self.space_interval)
        dr = self.R/self.space_interval
        r_s = np.linspace(dr, self.R, self.space_interval)
        return l_s, r_s

    def get_radial_volumes(self):
        dz = self.L/self.space_interval
        dr = self.R/self.space_interval
        r_s = np.linspace(dr, self.R, self.space_interval)
        vol_r_intervals = [self.cylinder_vol(r, dr, dz) for r in r_s]
        vol_sum = sum(vol_r_intervals)
        return np.array([i/vol_sum for i in vol_r_intervals])

    def calculate_output_conversion(self):
        output_slc = self.all_data[-1, -1, :, :]
        vol_weights = self.get_radial_volumes()
        self.output_conversion = output_slc[:, -1]@vol_weights
        logger.debug(f"Output Conversion: {self.output_conversion}")
        return self.output_conversion

    def get_heater_total(self):
        return self.heater_total

    def get_pressure_drop(self):
        Q = self.flow_rate
        rho = self.density
        rel_rough = const.rough_steel/self.D
        # Calculating reynolds number

        re = self.velocity*self.D/const.feed_viscosity # This is kinematic viscosity 
        logger.info(f"Reynolds Number: {re}")
        t4 = (7.149/re) ** 0.8981
        t3 = (rel_rough ** 1.1098)/2.857
        t_34 = math.log(t3+t4, 10)
        c_34 = 5.0452/re
        t2 = c_34*t_34
        t1 = -4 * math.log((rel_rough/3.7065 - t2), 10)
        ff = 1/(t1*t1)
        logger.info(f"Friction Factor {ff:.2f}")
        c_ff = 32*ff/(math.pi**2)
        l_d5 = self.L/(self.D**5)
        p_drop = c_ff*l_d5*Q*Q*rho
        self.pressure_drop = p_drop
        logger.info(f"Pressure Drop:{self.pressure_drop} kPa")
