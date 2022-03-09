from simulator.reactor.pfr import RealPFR
from . import const as const
import math
# Metrics
# Calculation of different optimization metrics

# Fixed cost that will not change based on operating conditions
FIXED_COST = 0


def operating_cost(reactor: RealPFR, time_interval: float):
    """Function to calculate the operating cost that is required. Note that this
    function only considers the costs which have a large impact on the reactor
    Some costs like pumping cost are excluded from optimization as they have 
    relatively negligible cost and cannot be explicity computed without using Hysys

    Args:
        reactor (RealPFR): PFR Reactor object
    """
    # Pre-Reactor Costing

    # Pumping of Feed
    # This part is largely insignificatn

    # Heating of IPA and Feed Mixture

    # Calculating Chemical Cost
    ipa_feed = reactor.ipa_feed*60  # kmol/min -> kmol/hr
    feed_temp = reactor.feed_temp

    ptsa_loading = reactor.flow_rate/1000*5  # m^3/ min -> 5g dm3/min
    ptsa_cost = ptsa_loading*const.ptsa_cost_per_g

    ipa_wt_rate = ipa_feed*const.ipa_wt  # kg/hr
    ipa_opex = ipa_wt_rate*const.ipa_cost*8000  # kg/hr * 8000 hr/yr = kg/year

    conversion = reactor.output_conversion
    ipp_out = reactor.pa_feed*conversion*60  # kmol/min -> kmol/hr
    ipp_out_wt = ipp_out * const.ipp_wt  # kmol/hr -> kg/hr
    ipp_profit = ipp_out*const.ipp_cost  # kg/hr * 8000 hr/yr = kg/year

    # Heating Requirement
    ss_heater_flow = reactor.get_heater_total()[-1]
    ss_heater_flow = ss_heater_flow/time_interval * 60  # kJ/h
    pass


def find_P(T: float):
    """Find pressure of feed required for the temperature

    Args:
        T (float): Desired temperature of the feed
    """
    pass


def capital_cost(reactor: RealPFR, n_split: int):

    pass


def profit():
    pass


def calculate_props_after_mixing(ipa_feed,):
    """Helper function to calculate properties and cp of the feed after mixing


    Args:
        ipa_feed (float): Molar flow of IPA into the reactor
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


def get_shell_mass(L, D, thickness, rho=8000):
    return math.pi*L*D * thickness*rho


def pressure_drop(L, D, velocity: float, e: float):
    rel_rough = e/D
    # Calculating reynolds number

    re = velocity*D/const.feed_viscosity
    t4 = (7.149/re) ^ 0.8981
    t3 = (rel_rough ^ 1.1098)/2.857
    t_34 = math.log(t3+t4, 10)
