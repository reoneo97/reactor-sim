from bayes_opt import BayesianOptimization
from bayes_opt.logger import JSONLogger
from bayes_opt.event import Events
import matplotlib.pyplot as plt
from simulator.reactor import RealPFR
from simulator.metrics import total_cost, get_tac, combine_costs
from loguru import logger
from datetime import datetime
import pandas as pd


# Bounded region of parameter space
pbounds = {
    'L': (1., 25.), 'D': (0.25, 6), 'M': (1., 5.), 'feed_temp': (313.15, 433.15),
    'n_reactor': (1, 8.99)
}
info = []

blank = {
    "Pump Operating Cost": 0, "Preheat Operating Cost": 0, "PTSA Reactant Cost": 0,
    "IPA Reactant Cost": 0, "IPP Reactant Profit": 0, "Reactor Operating Cost": 0,
    "Reactor Wastewater Cost": 0, "Preheat Capital Cost": 0, "Reactor Body Capital Cost": 0,
    "Reactor Head Capital Cost": 0, "Reactor Jacket Capital Cost": 0, "conversion": 0,
    "CC": 0, "OC": 0, "TAC": 1e11
}


def run_simulation(L, D, M, feed_temp, n_reactor, time_interval=0.2, time_end=500,
                   payback: int = 3):
    n_reactor = int(n_reactor)
    model = RealPFR(80.31/60/n_reactor, M, L, D/2, feed_temp, 433.15, 26)

    try:
        model.run(time_interval, time_end)
        conversion = model.output_conversion

        cost_info = total_cost(model, time_interval, n_reactor)
        cost_info["conversion"] = conversion
        op_cost, cap_cost, reactant_cost = combine_costs(cost_info)
        tac = op_cost + cap_cost/payback

        cost_info["CC"] = cap_cost
        cost_info["OC"] = op_cost
        cost_info["TAC"] = tac
        info.append(cost_info)
        return LAMBDA*conversion - tac
    except AssertionError:
        info.append(blank)
        return -1e13


LAMBDA = 1e9

optimizer = BayesianOptimization(
    f=run_simulation,
    pbounds=pbounds,
    random_state=1,
)
logger = JSONLogger(path="./logs/logs.json")
optimizer.subscribe(Events.OPTIMIZATION_STEP, logger)
try:
    optimizer.maximize(
        init_points=20,
        n_iter=70,
    )
    info_df = pd.DataFrame(info)
except:
    pass

# for i, res in enumerate(optimizer.res):
#     print("Iteration {}: \n\t{}".format(i, res))
param_info = [i["params"] for i in optimizer.res]
param_df = pd.DataFrame(param_info)
print("Param df:", len(param_df), "info_df", len(info_df))
df = pd.concat([param_df, info_df, ], axis=1)
print(df)
df.to_csv(
    f"./data/{datetime.now().strftime('%Y-%m-%d %H-%M')}_{LAMBDA:.2e}.csv", index=False)
print(optimizer.max)
