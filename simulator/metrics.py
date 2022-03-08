from simulator.reactor.pfr import RealPFR
from . import const as const
# Metrics
# Calculation of different optimization metrics


def operating_cost(reactor: RealPFR):

    # Calculating Chemical Cost
    ipa_feed = reactor.ipa_feed*60  # kmol/min -> kmol/hr
    feed_temp = reactor.feed_temp

    ptsa_loading = reactor.flow_rate/1000*5  # m^3/ min -> 5g dm3/min
    ptsa_cost = ptsa_loading*const.ptsa_cost_per_g

    ipa_wt_rate = ipa_feed*const.ipa_wt  # kg/hr
    ipa_opex = ipa_wt_rate*const.ipa_cost*8000  # kg/hr * 8000 hr/yr = kg/year

    conversion = reactor.output_conversion
    ipp_out = reactor.pa_feed*conversion*60 #kmol/min -> kmol/hr
    ipp_out_wt = ipp_out * const.ipp_wt # kmol/hr -> kg/hr
    ipp_profit = ipp_out*const.ipp_cost  # kg/hr * 8000 hr/yr = kg/year

    # Calculating Duty Cost
    mix_temp, cp = calculate_props_after_mixing(ipa_feed)
    preheat_duty = (feed_temp - mix_temp)*cp  # /hr cost
    preheat_cost = preheat_duty*8000/const.mps_cost_gj
    # Calculating Capital Cost


def capital_cost():
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
