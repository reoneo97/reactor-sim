from simulator.reactor.pfr import RealPFR
from . import const as const
# Metrics
# Calculation of different optimization metrics

# Fixed cost that will not change based on operating conditions
FIXED_COST = 0


def operating_cost(reactor: RealPFR):
    """Function to calculate the operating cost that is required 

    Args:
        reactor (RealPFR): PFR Reactor object
    """
    # Pre-Reactor Costing
    # feed_P =
    # Condense Feed Stream into Liquid

    # Pumping of IPA

    # Pumping of Feed

    # Heating of IPA and Feed Mixture

    ptsa_loading = reactor.flow_rate/1000*5  # m^3/ min - 5g dm3/min
    ptsa_cost = ptsa_loading*const.ptsa_cost_per_g

    ipa_feed = reactor.ipa_feed  # kmol/min
    ipa_wt_rate = ipa_feed*const.ipa_wt  # kg/min
    ipa_opex = ipa_wt_rate*const.ipa_cost

    conversion = reactor.output_conversion
    ipp_out = reactor.pa_feed*conversion
    ipp_profit = ipp_out*const.ipp_cost

    # Separation Calculations

    # Cool Feed Down to Separation Conditions

    # Heating Requirement

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
