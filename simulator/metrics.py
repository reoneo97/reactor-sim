from simulator.reactor.pfr import RealPFR
from . import const as const
import math
from loguru import logger
from typing import Dict
# Metrics
# Calculation of different optimization metrics

# Fixed cost that will not change based on operating conditions
FIXED_COST = 0


def total_cost(reactor: RealPFR, time_interval: float, n_reactor: int, cepci: int = 773):
    """Function to calculate the operating cost that is required. Note that this
    function only considers the costs which have a large impact on the reactor
    Some costs like pumping cost are excluded from optimization as they have 
    relatively negligible cost and cannot be explicity computed without using Hysys

    Args:
        reactor (RealPFR): PFR Reactor object
    """
    # Pre-Reactor Costing
    results = dict()
    ipa_feed = reactor.ipa_feed*60  # kmol/min -> kmol/hr
    conversion = reactor.output_conversion
    D = reactor.D
    L = reactor.L
    VOL = reactor.vol
    reactor_flow_rate = reactor.flow_rate # m3/min
    # Recovery Values - Based on preliminary hysys separation for base case
    ipp_recovery = 1  # How much IPP can be obtained at 99 % purity after separation
    ptsa_recovery = 1
    # Pumping of Feed
    pump_feed_duty = 13.19  # kW for whole feed
    pump_ipa_duty = 0.0413*ipa_feed  # kW, reactor feed
    pump_duty = pump_feed_duty + pump_ipa_duty
    pump_duty_annual = pump_duty * 3600 * 8000  # kJ/hr
    pump_duty_opex = pump_duty_annual/1e6 * const.electricity_cost
    # Heating of IPA and Feed Mixture

    # Feed_cp here is in terms of kJ/ C instead of kJ/kmol C. Dont need to multiply with anything
    inlet_temp, feed_cp = calculate_props_after_mixing(ipa_feed)
    # Calculating Chemical Cost
    feed_temp = reactor.feed_temp
    preheat_duty = (feed_temp - inlet_temp)*feed_cp
    preheat_duty_kw = preheat_duty/3600  # Duty in KW
    preheat_cost = preheat_duty*8000/1e6*const.mps_cost_gj

    ptsa_loading = reactor.pa_flow*1000*5  # m^3/ min * 5 g/dm3 -> g/min
    ptsa_cost = ptsa_loading*60*8000*const.ptsa_cost_per_g*ptsa_recovery

    ipa_wt_rate = ipa_feed*const.ipa_wt  # kg/hr
    ipa_opex = ipa_wt_rate*const.ipa_cost*8000  # kg/hr * 8000 hr/yr = kg/year

    ipp_out = reactor.pa_feed*conversion*60  # kmol/min -> kmol/hr
    ipp_out_wt = ipp_out * const.ipp_wt  # kmol/hr -> kg/hr
    ipp_profit = ipp_out_wt*const.ipp_cost * 8000 *\
        ipp_recovery  # kg/hr * 8000 hr/yr = kg/year

    # Heating Requirement
    ss_heater_flow = reactor.get_heater_total()[-1]  # kJ/min
    ss_heater_flow = ss_heater_flow/time_interval * 60  # kJ/h
    ss_heater_flow_annual = ss_heater_flow*8000 / 1e6  # GJ
    ss_heater_cost = ss_heater_flow_annual*const.mps_cost_gj

    # Calculation of IPA disposal
    # Estimating 100% Separation
    ipp_vol = (ipp_out_wt/60)/const.ipp_density
    ww_vol = reactor_flow_rate - ipp_vol
    logger.info(f"Wastewater Volume: {ww_vol:.3f} m3/min")
    ww_cost = wastewater_cost(ww_vol)

    pump_capex = get_pump_cc(reactor_flow_rate) 
    # Capital Cost Calculation
    heater_cc = get_heater_cc(inlet_temp, feed_temp, preheat_duty_kw)
    reactor_pv_cc = pressure_vessel(D, L)
    reactor_head_cc = ellipsoidal_head(D)*2  # for 2 head

    reactor_jacket_cc = get_reactor_jacket_cc(VOL)

    results["Pump Operating Cost"] = pump_duty_opex
    results["Preheat Operating Cost"] = preheat_cost
    results["PTSA Reactant Cost"] = ptsa_cost
    results["IPA Reactant Cost"] = ipa_opex
    results["IPP Reactant Profit"] = ipp_profit
    results["Reactor Operating Cost"] = ss_heater_cost * n_reactor
    results["Reactor Wastewater Cost"] = ww_cost
    results["Pump Capital Cost"] = pump_capex * cepci/500
    results["Preheat Capital Cost"] = heater_cc * cepci/500
    results["Reactor Body Capital Cost"] = reactor_pv_cc * \
        n_reactor * cepci/500
    results["Reactor Head Capital Cost"] = reactor_head_cc * \
        n_reactor * cepci/500 *2
    results["Reactor Jacket Capital Cost"] = reactor_jacket_cc * \
        n_reactor * cepci/500
    return results


def pressure_vessel(D, L, P=15, s=const.pv_stress_max, e=1.0):
    # bar = N/mm^2
    num = P*(D*1000)  # m to mm
    denom = (4*s*e) - (0.4*P)

    thickness = num/denom + const.pv_corr_allowance  # Thickness in mm
    logger.debug(f"Pressure Vessel Thickness: {thickness:.3f} mm")
    thickness = min(thickness, const.pv_thick_min)
    OD = D+thickness/500
    shell_weight = math.pi*OD*L*(thickness/1000)*const.steel_density
    logger.debug(f"Body Shell Weight: {shell_weight:.3f} kg")
    return shell_cost(shell_weight)


def ellipsoidal_head(D, P=15, s=const.pv_stress_max, e=1.0):
    # bar = N/mm^2
    num = P*(D*1000)  # m to mm
    denom = (2*s*e) - (0.2*P)

    thickness = num/denom + const.pv_corr_allowance
    logger.debug(f"Pressure Head Thickness: {thickness:.3f} mm")
    thickness = min(thickness, const.pv_thick_min)
    OD = D + (2*thickness/1000)

    shell_weight = (math.pi/24)*(OD**3-D**3)*const.steel_density
    logger.debug(f"Head Shell Weight: {shell_weight:.3f} kg")
    return shell_cost(shell_weight)


def shell_cost(weight: float):
    return const.pv_304_a + const.pv_304_b*(weight**const.pv_304_n)


def get_reactor_jacket_cc(vol):

    return const.jacket_reactor_a + const.jacket_reactor_b * (vol ** const.jacket_reactor_n)


def temp_lm(h_i: float, h_o: float, c_i: float, c_o: float):
    """Calculate the log-mean temperature difference of a heat exchanger given
    the temperature streams

    Args:
        h_i (float): Temperature of hot in
        h_o (float): Temperature of hot out
        c_i (float): Temperature of cold in 
        c_o (float): Temperature of cold out
    """
    dt1 = h_i - c_o
    dt2 = h_o - c_i
    return dt1 - dt2/(math.log(dt1/dt2))

def get_pump_cc(reactor_flow_rate:float):
    """Calculate capital cost for pump
    Args:
        reactor_flow_rate (float): m3/min
    """
    fr_ls = reactor_flow_rate/60*1000
    return const.pump_a + const.pump_b*fr_ls**const.pump_n

def get_heater_cc(inlet_temp: float, feed_temp: float, duty: float):
    """Calculate the capital cost of the heater given the feed temperature 
    that is required. Assuming that the heating fluid is medium pressure steam.


    Args:
        feed_temp (float): _description_
    """
    h_i = const.mps_temp
    h_o = h_i - 0.1
    c_i = inlet_temp
    c_o = feed_temp

    dt_lm = temp_lm(h_i, h_o, c_i, c_o)
    U = 1/(1/4.5 + 1/0.3)  # MPS Heat Transfer Coefficient kW/m2K
    logger.debug(f"Preheat Heat Transfer Coefficient: {U:.3f}")

    A = duty/(dt_lm*U)
    logger.debug(f"Heat Exchanger Area: {A:.3f} m2")
    return 24000 + 46*A**1.2


def calculate_props_after_mixing(ipa_feed,):
    """Helper function to calculate properties and cp of the feed after mixing
    Assumes that C_p can be taken as constant over the interval

    Args:
        ipa_feed (float): Molar flow of IPA into the reactor. kmol/hr
    """
    ipa_feed_wt = ipa_feed*const.ipa_wt

    ipa_cp = 2.978
    ipa_temp = 25
    # Constants Left here for calculation
    feed_1_mass_flow = 2712
    feed_1_cp = 2.204
    feed_1_temp = 45

    feed_5_mass_flow = 18810
    feed_5_cp = 2.132
    feed_5_temp = 30.07

    # total_mass_flow = ipa_feed_wt + feed_1_mass_flow + feed_5_mass_flow

    cp = ipa_feed_wt*ipa_cp + feed_1_cp*feed_1_mass_flow + feed_5_mass_flow*feed_5_cp

    heat_total = ipa_feed_wt*ipa_cp*ipa_temp + feed_1_cp * \
        feed_1_mass_flow*feed_1_temp + feed_5_mass_flow*feed_5_cp*feed_5_temp
    final_temp = heat_total/cp

    return final_temp, cp


def pressure_drop(model:RealPFR):
    L = model.L
    D = model.D
    velocity = model.velocity
    Q = model.flow_rate
    rho = model.density
    rel_rough = const.rough_steel/D
    # Calculating reynolds number

    re = velocity*D/const.feed_viscosity
    t4 = (7.149/re) ** 0.8981
    t3 = (rel_rough ** 1.1098)/2.857
    t_34 = math.log(t3+t4, 10)
    c_34 = 5.0452/re
    t2 = c_34*t_34
    t1 = -4* math.log((rel_rough/3.7065 - t2),10)
    ff = 1/(t1*t1)
    
    c_ff = 32*ff/(math.pi**2)
    l_d5 = L/(D**5)
    p_drop = c_ff*l_d5*Q*Q*rho
    return p_drop



def wastewater_cost(q: float):
    # q must be gal/min
    q_gal = q*264.172
    return 48760*q_gal**0.64


def combine_costs(results: Dict[str, float]):
    op_cost = 0
    cap_cost = 0
    reactant_cost = 0
    for k, v in results.items():
        if "Operating" in k:
            op_cost += v
        elif "Capital" in k:
            cap_cost += v
        elif "Reactant" in k:
            reactant_cost += v
    return op_cost, cap_cost, reactant_cost


def get_tac(reactor: RealPFR, time_interval: float, n_reactor: int, payback: int = 3):
    costs = total_cost(reactor, time_interval, n_reactor)
    op_cost, cap_cost, reactant_cost = combine_costs(costs)
    return op_cost + cap_cost/payback
