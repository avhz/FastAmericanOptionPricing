import matplotlib.pyplot as plt
import numpy as np
from QDplusAmericanOptionSolver import *
from EuropeanOptionSolver import *

if __name__ == '__main__':
    # unit test one for valuing American option
    r = 0.0975729097939295     # risk free
    q = 0.011804520625954162      # dividend yield
    K = 105.61782314803582       # strike
    S0 = 30.543317986992072        # underlying spot
    sigma = 0.2  # volatility
    T = 3  # maturity
    option_type = OptionType.Put
    tau = 0.000870786986

    solver = QDplus(r, q, sigma, K, option_type)
    print("Am price =", solver.price(tau, S0), "exercise boundary = ", solver.exercise_boundary)

    S = np.linspace(1, 4*S0, 200)
    plt.plot(S, solver.exercise_boundary_func(S, tau), 'o-')
    plt.plot([0, 4*S0], [0, 0], 'r--')
    plt.ylim([-2*K, 2 * K])
    plt.ylabel("target function")
    plt.xlabel("S*")
    plt.show()
