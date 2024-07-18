// import numpy as np
// import scipy.stats as stats
// import scipy.optimize
// import EuropeanOptionSolver as europ
// from enum import Enum

// class OptionType(Enum):
//     Call = 1
//     Put = 2

use crate::{european_pricer, TypeFlag};

#[derive(Default)]
pub struct QDPlus {
    pub riskfree: f64,
    pub dividend: f64,
    pub volatility: f64,
    pub strike: f64,
    pub option_type: TypeFlag,

    // backend fields for computation
    backend: QDPlusBackend,
}

#[derive(Default)]
struct QDPlusBackend {
    v_M: f64,
    v_N: f64,
    v_h: f64,
    v_qQD: f64,
    v_qQDdot: f64,
    v_p: f64,
    v_theta: f64,
    v_c: f64,
    v_b: f64,
    v_d1: f64,
    v_d2: f64,
    v_dlogSdh: f64,
    exercise_boundary: f64,
}

impl QDPlus {
    // class QDplus:
    //     """QD+ alogrithm for computing approximated american option price"""
    //     def __init__(self, riskfree, dividend, volatility, strike, option_type):
    //         self.r = riskfree
    //         self.q = dividend
    //         self.sigma = volatility
    //         self.K = strike
    //         self.option_type = option_type
    //         if option_type == OptionType.Call:
    //             self.option_indicator = 1
    //         else:
    //             self.option_indicator = -1
    //         # miscellaneous with tau only
    //         self.v_M = 0
    //         self.v_N = 0
    //         self.v_h = 0
    //         self.v_qQD = 0
    //         self.v_qQDdot = 0
    //         # miscellaneous terms with tau and S
    //         self.v_p = 0
    //         self.v_theta = 0
    //         self.v_c = 0
    //         self.v_b = 0
    //         self.v_d1 = 0
    //         self.v_d2 = 0
    //         self.v_dlogSdh = 0
    //         self.exercise_boundary = 0
    //         self.tolerance = 1e-10
    fn new(
        riskfree: f64,
        dividend: f64,
        volatility: f64,
        strike: f64,
        option_type: TypeFlag,
    ) -> Self {
        Self {
            riskfree,
            dividend,
            volatility,
            strike,
            option_type,

            backend: QDPlusBackend::default(),
        }
    }

    //     def price(self, tau, S):
    //         if tau == 0:
    //             self.exercise_boundary = self.K
    //             return max(S-self.K, 0.0)
    //         self.exercise_boundary = Sb = self.compute_exercise_boundary(tau)
    //         err = self.exercise_boundary_func(Sb, tau)
    //         print("err = ", err)
    //         self.compute_miscellaneous(tau, Sb)
    //         qQD = self.v_qQD
    //         c = self.v_c
    //         b = self.v_b
    //         if self.option_type == OptionType.Put:
    //             pS = europ.EuropeanOption.european_put_value(tau, S, self.r, self.q, self.sigma, self.K)
    //             pSb = europ.EuropeanOption.european_put_value(tau, Sb, self.r, self.q, self.sigma, self.K)
    //         else:
    //             pS = europ.EuropeanOption.european_call_value(tau, S, self.r, self.q, self.sigma, self.K)
    //             pSb = europ.EuropeanOption.european_call_value(tau, Sb, self.r, self.q, self.sigma, self.K)

    //         if self.option_indicator * (Sb - S) <= 0:
    //             return self.option_indicator * (S - self.K)
    //         else:
    //             return pS + (self.K - Sb - pSb)/(1 - b * np.square(np.log(S/Sb)) - c * np.log(S/Sb)) * np.power(S/Sb, qQD)
    fn price(&self, tau: f64, s: f64) -> f64 {
        if tau == 0.0 {
            self.backend.exercise_boundary = self.strike;
            return (s - self.strike).max(0.0);
        }

        self.backend.exercise_boundary = self.compute_exercise_boundary(tau);
        let Sb = self.backend.exercise_boundary;

        let err = self.exercise_boundary_func(Sb, tau);
        println!("err = {}", err);

        self.compute_miscellaneous(tau, Sb);

        let qQD = self.backend.v_qQD;
        let c = self.backend.v_c;
        let b = self.backend.v_b;

        let pS = match self.option_type {
            TypeFlag::Put => european_pricer::european_put_value(
                tau,
                s,
                self.riskfree,
                self.dividend,
                self.volatility,
                self.strike,
            ),
            TypeFlag::Call => european_pricer::european_call_value(
                tau,
                s,
                self.riskfree,
                self.dividend,
                self.volatility,
                self.strike,
            ),
        };
        let pSb = match self.option_type {
            TypeFlag::Put => european_pricer::european_put_value(
                tau,
                Sb,
                self.riskfree,
                self.dividend,
                self.volatility,
                self.strike,
            ),
            TypeFlag::Call => european_pricer::european_call_value(
                tau,
                Sb,
                self.riskfree,
                self.dividend,
                self.volatility,
                self.strike,
            ),
        };

        if self.option_type as i32 * (Sb - s) <= 0 {
            return (s - self.strike).max(0.0);
        } else {
            return pS
                + (self.strike - Sb - pSb) / (1.0 - b * (s / Sb).powi(2) - c * (s / Sb).ln())
                    * (s / Sb).powf(qQD);
        }
    }

    //     def compute_exercise_boundary(self, tau):
    //         if tau == 0:
    //             return self.B_at_zero()
    //         # using x0->0 is critical since there are multiple roots for the target function
    //         res = scipy.optimize.root(self.exercise_boundary_func,x0=self.K, args=(tau,))
    //         #if res.success == False:
    //          #   print("succuess? ", res.success, ", ", res.message, ", res = ", res.x)
    //         return res.x[0]

    //     def B_at_zero(self):
    //         if self.option_type == OptionType.Call:
    //             if self.r <= self.q:
    //                 return self.K
    //             else:
    //                 return self.r/self.q * self.K
    //         else:
    //             if self.r >= self.q:
    //                 return self.K
    //             else:
    //                 return self.r/self.q * self.K
    fn B_at_zero(&self) -> f64 {
        match self.option_type {
            TypeFlag::Call => {
                if self.riskfree <= self.dividend {
                    return self.strike;
                } else {
                    return self.riskfree / self.dividend * self.strike;
                }
            }
            TypeFlag::Put => {
                if self.riskfree >= self.dividend {
                    return self.strike;
                } else {
                    return self.riskfree / self.dividend * self.strike;
                }
            }
        }
    }

    //     def compute_miscellaneous(self, tau, S):
    //         #order cannot be changed
    //         self.v_N = self.N()
    //         self.v_M = self.M()
    //         self.v_h = self.h(tau)
    //         self.v_qQD = self.q_QD(tau)
    //         self.v_qQDdot = self.q_QD_dot()
    //         self.v_d1 = europ.EuropeanOption.d1(tau, S, self.r, self.q, self.sigma, self.K)
    //         self.v_d2 = europ.EuropeanOption.d2(tau, S, self.r, self.q, self.sigma, self.K)
    //         if self.option_type == OptionType.Put:
    //             self.v_p = europ.EuropeanOption.european_put_value(tau, S, self.r, self.q, self.sigma, self.K)
    //         else:
    //             self.v_p = europ.EuropeanOption.european_call_value(tau, S, self.r, self.q, self.sigma, self.K)
    //         self.v_theta = europ.EuropeanOption.european_option_theta(tau, S, self.r, self.q, self.sigma, self.K)
    //         self.v_dlogSdh = self.dlogSdh(tau, S)
    //         self.v_c = self.c(tau, S)
    //         self.v_c0 = self.c0(tau, S)
    //         self.v_b = self.b(tau, S)

    //     def exercise_boundary_func(self, S, tau):
    //         if tau == 0:
    //             if type(S) is float:
    //                 return 0
    //             else:
    //                 return np.ones(S.size) * 0
    //         self.compute_miscellaneous(tau, S)
    //         qQD = self.v_qQD
    //         p = self.v_p
    //         c0 = self.v_c0
    //         d1 = self.v_d1
    //         if self.option_type == OptionType.Call:
    //             ans = (1 - np.exp(-self.q * tau) * stats.norm.cdf(d1)) * S - (qQD) * (S - self.K - p)
    //         else:
    //             ans = (1 - np.exp(-self.q * tau) * stats.norm.cdf(-d1)) * S + (qQD) * (self.K - S - p)
    //         return ans

    //     def q_QD(self, tau):
    //         N = self.v_N
    //         M = self.v_M
    //         h = self.v_h
    //         if self.option_type == OptionType.Call:
    //             return -0.5*(N-1) + 0.5 * np.sqrt((N-1)*(N-1) + 4 * M/h)
    //         else:
    //             return -0.5*(N-1) - 0.5 * np.sqrt((N-1)*(N-1) + 4 * M/h)

    //     def q_QD_dot(self):
    //         N = self.v_N
    //         M = self.v_M
    //         h = self.v_h
    //         return M/(h * h * np.sqrt((N-1)*(N-1) + 4*M/h))

    //     def c0(self, tau, S):
    //         N = self.v_N
    //         M = self.v_M
    //         h = self.v_h
    //         qQD = self.v_qQD
    //         qQDdot = self.v_qQDdot
    //         p = self.v_p
    //         theta = self.v_theta
    //         c = self.v_c
    //         d1 = self.v_d1
    //         d2 = self.v_d2
    //         return - (1-h)*M/(2*qQD + N - 1) * (1/h - (theta*np.exp(self.r * tau))/(self.r*(self.K - S - p)) + qQDdot/(2*qQD+N-1))

    //     def c(self, tau, S):
    //         r = self.r
    //         q = self.q
    //         N = self.v_N
    //         M = self.v_M
    //         h = self.v_h
    //         qQD = self.v_qQD
    //         qQDdot = self.v_qQDdot
    //         p = self.v_p
    //         theta = self.v_theta
    //         c = self.v_c
    //         d1 = self.v_d1
    //         d2 = self.v_d2
    //         dlogSdh = self.v_dlogSdh
    //         c0 = self.c0(tau, S)

    //         return c0 - ((1-h)*M)/(2*qQD + N - 1) \
    //             * ((1 - np.exp(-q * tau)*stats.norm.cdf(-d1))/(self.K - S - p) + qQD/S)\
    //             * dlogSdh

    //     def b(self, tau, S):
    //         N = self.v_N
    //         M = self.v_M
    //         h = self.v_h
    //         qQD = self.v_qQD
    //         qQDdot = self.v_qQDdot
    //         p = self.v_p
    //         theta = self.v_theta
    //         c = self.v_c
    //         d1 = self.v_d1
    //         d2 = self.v_d2
    //         return ((1-h)*M*qQDdot)/(2*(2*qQD + N - 1))

    //     def dlogSdh(self, tau, S):
    //         N = self.v_N
    //         M = self.v_M
    //         h = self.v_h
    //         qQD = self.v_qQD
    //         qQDdot = self.v_qQDdot
    //         p = self.v_p
    //         theta = self.v_theta
    //         c = self.v_c
    //         d1 = self.v_d1
    //         d2 = self.v_d2
    //         r = self.r
    //         q = self.q

    //         dFdh = qQD * theta * np.exp(self.r * tau)/self.r + qQDdot * (self.K - S - p) \
    //             + (S * self.q *np.exp(-self.q*tau) * stats.norm.cdf(-d1))/(r * (1-h)) \
    //             - (S * np.exp(-self.q * tau) * stats.norm.pdf(d1))/(2*r*tau*(1-h))\
    //             * (2*np.log(S/self.K)/(self.sigma * np.sqrt(tau)) - d1)

    //         dFdS = (1 - qQD) * (1 - np.exp(-q * tau) * stats.norm.cdf(-d1)) \
    //                 + (np.exp(-q * tau) * stats.norm.pdf(d1))/(self.sigma * np.sqrt(tau))
    //         return -dFdh/dFdS

    // def h(self, tau):
    //     return 1 - np.exp(-self.r * tau)
    fn h(&self, tau: f64) -> f64 {
        1.0 - f64::exp(-self.riskfree * tau)
    }

    //     def M(self):
    //         return 2 * self.r / (self.sigma * self.sigma)
    fn M(&self) -> f64 {
        2.0 * self.riskfree * self.volatility.powi(2).recip()
    }

    //     def N(self):
    //         return 2 * (self.r - self.q) / (self.sigma * self.sigma)
    fn N(&self) -> f64 {
        2.0 * (self.riskfree - self.dividend) * self.volatility.powi(2).recip()
    }
}
