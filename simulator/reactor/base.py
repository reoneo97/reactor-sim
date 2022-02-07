from .. import const
from ..reaction import ReactionModel

class Reactor:
    def __init__(self,
                 pa_feed: float, M: int, vol: float, temperature: float,
                 rxn_model: ReactionModel):
        self.temperature = temperature
        self.pa_feed = pa_feed
        self.M = M
        self.ipa_feed = M*pa_feed
        self.vol = vol

        pa_wt = pa_feed*const.palm_acid_wt
        pa_flow = pa_wt/const.palm_acid_density

        ipa_wt = self.ipa_feed*const.ipa_wt
        ipa_flow = ipa_wt/const.ipa_density

        # To Reach Steady State flow rate in and out must be the same
        self.flow_rate = ipa_flow + pa_flow
        self.reaction = rxn_model

    def run(self):
        raise NotImplementedError
