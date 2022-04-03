from .. import const
from ..reaction import ReactionModel


class Reactor:
    def __init__(self,
                 pa_feed: float, M: int, temperature: float, vol: float = None,
                 rxn_model: ReactionModel = None):
        """_summary_

        Args:
            pa_feed (float): pa feed in kmol/min
            M (int): _description_
            temperature (float): _description_
            vol (float, optional): _description_. Defaults to None.
            rxn_model (ReactionModel, optional): _description_. Defaults to None.
        """
        self.temperature = temperature
        self.pa_feed = pa_feed
        self.M = M
        self.ipa_feed = M*pa_feed
        self.vol = vol

        self.init_flow_rate()
        # pa_wt = pa_feed*const.palm_acid_wt
        # pa_flow = pa_wt/const.palm_acid_density

        # ipa_wt = self.ipa_feed*const.ipa_wt
        # ipa_flow = ipa_wt/const.ipa_density

        # # To Reach Steady State flow rate in and out must be the same
        # self.flow_rate = ipa_flow + pa_flow
        self.reaction = rxn_model

    def run(self):
        raise NotImplementedError

    def init_flow_rate(self):
        """
        Obtain volumetric flow rate for the reactor. Assumes that the volumetric
        flow rate is uniform in both the axial and radial direction

        All feed and flow rates are based on /min
        """
        pa_wt = self.pa_feed*const.palm_acid_wt
        self.pa_flow = pa_wt/const.palm_acid_density
        
        ipa_wt = self.ipa_feed*const.ipa_wt
        ipa_flow = ipa_wt/const.ipa_density
        self.mass_flow = pa_wt+ipa_wt
        self.flow_rate = ipa_flow + self.pa_flow  # Volumetric Flow

        self.molar_flow_rate = self.pa_feed + self.ipa_feed
        self.density = self.mass_flow/self.flow_rate