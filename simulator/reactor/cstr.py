from ..reaction import EquilibriumReaction


class CSTR:
    def __init__(self, rxn:EquilibriumReaction):
        self.rxn = rxn

    def volume_sizing(self, feed_rate:float, conversion:float):
        # rxn_rate = self.rxn.get_rate()        
        return 
