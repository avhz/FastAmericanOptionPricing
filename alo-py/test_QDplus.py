import matplotlib.pyplot as plt
import numpy as np
from QDplusAmericanOptionSolver import *
from EuropeanOptionSolver import *

if __name__ == "__main__":
    # unit test one for valuing American option
    K = 10.0  # strike
    r = 0.05  # risk free
    q = 0.02  # dividend yield
    v = 0.20  # volatility
    T = 1.00  # maturity
    option_type = OptionType.Put
    # tau = 0.000870786986

    solver = QDplus(r, q, v, K, option_type)

    print("Price (ATM):", solver.price(T, 10.0), "Bound. = ", solver.exercise_boundary)
    print("Price (ITM):", solver.price(T, 05.0), "Bound. = ", solver.exercise_boundary)
    print("Price (OTM):", solver.price(T, 15.0), "Bound. = ", solver.exercise_boundary)

    # S0 = 10.0  # underlying spot
    # S = np.linspace(1, 4 * S0, 200)
    # plt.plot(S, solver.exercise_boundary_func(S, T), "o-")
    # plt.plot([0, 4 * S0], [0, 0], "r--")
    # plt.ylim([-2 * K, 2 * K])
    # plt.ylabel("target function")
    # plt.xlabel("S*")
    # plt.show()
