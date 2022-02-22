from .const import const_R
from typing import List


class HeatCapacity():

    def __init__(self, A, B, C):
        self.A = A
        self.B = B
        self.C = C

    def __call__(self, T: float) -> float:
        """Returns Heat Capacity at that temperature


        """
        t1 = self.A + self.B*T
        if self.C:
            t1 += self.C*T*T
        return t1


class MixedHeatCapacity():
    def __init__(self, c_ps: List[HeatCapacity]):
        self.c_ps = c_ps

    def __call__(self, T: float, amts: List[float]) -> float:

        cps = [c_p(T) for c_p in self.c_ps]
        assert len(amts) == len(self.c_ps), \
            "There must be the same number of Cps as concentration. Currently" \
            f" there are {len(amts)} C and {len(self.c_ps)} C_p values"
        final_cp = [conc*cp for conc, cp in zip(amts, cps)]
        return final_cp


pa_cp = HeatCapacity(191.673, 1.23910971, -1.94285714e-04)
ipa_cp = HeatCapacity(515.69904582, -2.51967286, 0.00467143)
water_cp = HeatCapacity(114.76700833, -0.32799223, 0.00064457)
la_cp = HeatCapacity(407.05732346, 0.03792343, 0.00124571)
ipp_cp = HeatCapacity(229.40519469, 1.15042514, 6.85714286e-05)

cps = [pa_cp, ipa_cp, water_cp, la_cp, ipp_cp]
mh_cp = MixedHeatCapacity(cps)
