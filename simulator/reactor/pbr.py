import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

from .. import const
from ..reaction import ReactionModel
from ..reactor import Reactor

class IdealPBR(Reactor):
    def __init__(self,
                 pa_feed: float, M: int, weight: float, temperature: float,
                 rxn_model: ReactionModel):
        # Initial Setup in base model
        super().__init__(pa_feed, M, temperature, rxn_model)
