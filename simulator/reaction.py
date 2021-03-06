
from typing import List
from math import exp
from .const import const_R


class Reaction:
    """
    Reaction Model class used to store all important information about 
    reaction parameters. 
    """

    def __init__(self, A: float, E_a: float, rxn_order: List[int]):
        """[summary]

        Args:
            A : Pre-exponential Constant 
            E_a (float): Activation Energy in KJ mol
            rxn_order (List[int]): Orders for the reactants
        """

        self.A = A
        self.E_a = E_a
        # self.delta_H = delta_H
        self.rxn_order = rxn_order
        pass

    def arrhenius_eqn(self, T) -> float:
        """Function that returns the rate constant of the reaction given the 
        temperature as a variable.

        Rate Constant that is calculated is in units of mol-1 dm3 min-1
        Args:
            T (float): Temperature
        """
        assert T > 0, "Temperature has to be positive"
        exp_term = -self.E_a/(const_R*T)
        return self.A*exp(exp_term)

    def get_rate(self, conc: List[float], t: float):
        assert len(conc) == len(self.rxn_order), f"There should be {len(self.rxn_order)}" \
            f"concentration parameters supplied, {len(conc)} found"

        # Find Rate Constant k
        k = self.arrhenius_eqn(t)
        prod = 1  # Products
        for order, conc in zip(self.rxn_order, conc):
            prod *= conc**order
        return k*prod


class ArrheniusEquation:

    def __init__(self, A, Ea):
        self.A = A
        self.Ea = Ea

    def arrhenius_eqn(self, T) -> float:
        """Function that returns the rate constant of the reaction given the 
        temperature as a variable.

        Rate Constant that is calculated is in units of mol-1 dm3 min-1
        Args:
            T (float): Temperature
        """
        exp_term = -self.Ea/(const_R*T)
        return self.A*exp(exp_term)

    def __call__(self, T) -> float:
        return self.arrhenius_eqn(T)


class EquilibriumReaction:
    """
    Model to calculate equilibrium constant given parameters
    Model follows general equation of K = Aexp(|G|/RT)        
    """

    def __init__(self, A, G):
        self.A = A
        self.G = G

    def __call__(self, T: float) -> float:
        exp_term = self.G/(const_R*T)
        return self.A*exp(exp_term)


class ReactionModel():

    def __init__(self) -> None:
        super(ReactionModel, self).__init__()


class HomogeneousReaction(ReactionModel):
    """
    Reaction Kinetics model for an equilibrium reaction with homogeneous catalysis
    Equilibrium Reaction that contains both a forward and a backward reaction
    and follows the reaction be
    -r_a = k_1 c_a c_b - k_2 *  c_c * c_d
    """

    def __init__(self, forward_rxn: Reaction, backward_rxn: Reaction, heat: float):
        self.forward = forward_rxn
        self.backward = backward_rxn
        self.heat = heat  # Units of kJ/kmol

    def get_rate(self, ca: float, cb: float, ce: float, cw: float, T: float, heat=False):

        reactant_concs = [ca, cb]
        product_concs = [ce, cw]

        # First convert the reaction rates to s^-1 from h^-1
        forward_rate = self.forward.get_rate(reactant_concs, T)/3600
        backward_rate = self.backward.get_rate(product_concs, T)/3600
        net_rate = forward_rate - backward_rate
        if heat:
            return net_rate, -net_rate*self.heat  # Endothermic Reaction will consume heat
        return net_rate

    def get_heat(self):
        return -self.heat


class LHHWReaction(ReactionModel):
    """
    Reaction Kinetics model for an equilibrium reaction with heterogeneous catalysis
    Follows the LHHW Model where the reaction is surface reaction controlled

    """

    def __init__(self, ks: ArrheniusEquation, KA: EquilibriumReaction,
                 KB: EquilibriumReaction, KE: EquilibriumReaction,
                 KW: EquilibriumReaction, KS: EquilibriumReaction,
                 heat: float):

        self.ks = ks
        self.KA = KA
        self.KB = KB
        self.KE = KE
        self.KW = KW
        self.KS = KS
        self.heat = heat

    def get_rate(self, ca, cb, ce, cw,  T):
        assert all([i >= 0 for i in [ca, cb, ce, cw, T]]
                   ), "Concentration values must be postive"
        # Calculate Surface Rxn Rate Constant
        ks = self.ks(T)
        # Calculate Adsorption Equilibrium Constant
        ka = self.KA(T)
        kb = self.KB(T)
        ke = self.KE(T)
        kw = self.KW(T)
        KS = self.KS(T)
        # Denominator Calculation
        denom = (1 + ka*ca + kb*cb + ke*ce + kw*cw)**2
        # print(denom)
        # Numerator Calculation
        term1 = ka*kb*ca*cb
        term2 = ce*cw*ke*kw/(KS)
        num = ks*(term1 - term2)
        return num/denom


def __get_ptsa_reaction():
    # Reaction Constants are given in units of h^-1
    forward = Reaction(150241, 40770, [1, 1])
    backward = Reaction(7.1166, 13128.3, [1, 1])
    heat = 7984  # kJ mol
    rxn_model = HomogeneousReaction(forward, backward, heat=heat)
    # print(rxn_model.get_rate([1, 100], [1, 1], 500))
    return rxn_model


def __get_zna_reaction():

    KA = EquilibriumReaction(4.71397e-19, 105361.69)
    KB = EquilibriumReaction(4.967627e-17, 87663.523)
    KE = EquilibriumReaction(4.58858e-15, 87893.245)
    KW = EquilibriumReaction(3.4822796e-16, 88018.49)
    ks = ArrheniusEquation(1257.925, 41589.037)
    KS = EquilibriumReaction(75074.872, -32827.8428)
    heat = 7984e3
    reaction_model = LHHWReaction(ks, KA, KB, KE, KW, KS, heat=heat)
    return reaction_model
# print(test(1))


ptsa_reaction = __get_ptsa_reaction()
zna_reaction = __get_zna_reaction()

# print(zna_reaction.get_rate(10, 10, 0, 0, 400)*5)

# print(zna_reaction.ks(373))
