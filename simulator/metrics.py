from simulator.reactor.pfr import RealPFR
from . import const as const
# Metrics
# Calculation of different optimization metrics


def operating_cost(reactor: RealPFR):

    ptsa_loading = reactor.flow_rate/1000*5  # m^3/ min - 5g dm3/min
    ptsa_cost_min = ptsa_loading*const.ptsa_cost_per_g

    ipa_feed = reactor.ipa_feed  # kmol/min
    ipa_wt_rate = ipa_feed*const.ipa_wt  # kg/min
    ipa_opex = ipa_wt_rate*const.ipa_cost

    conversion = reactor.output_conversion
    ipp_out = reactor.pa_feed*conversion

    pass


def capital_cost():
    pass


def profit():
    pass
