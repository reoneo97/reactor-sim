from .const import const_R


class HeatCapacity():

    def __init__(self, A, B, C):
        self.A = A
        self.B = B*1e-3
        self.C = C*1e-6

    def __call__(self, T: float) -> float:
        """Returns Heat Capacity at that temperature


        """
        t1 = self.A + self.B*T
        if self.C:
            t1 += self.C*T*T
        return t1*const_R


water_cp = HeatCapacity(4.038, 3.470, 1.450)
 