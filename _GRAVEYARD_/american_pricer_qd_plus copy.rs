// import numpy as np
// import scipy.stats as stats
// import scipy.optimize
// import EuropeanOptionSolver as europ

use crate::{
    european_pricer::{d1, d2, european_call_value, european_option_theta, european_put_value},
    OptionType,
};
use argmin::{
    core::{CostFunction, Executor, State},
    solver::brent::BrentRoot,
};
use statrs::distribution::{Continuous, ContinuousCDF, Normal};
use time::{Date, OffsetDateTime};

#[derive(Debug, Clone, Copy)]
pub struct QDplus {
    pub r: f64,
    pub q: f64,
    pub v: f64,
    pub K: f64,
    pub option_type: OptionType,
    pub expiry_date: Date,

    tolerance: f64,
    // PRIVATE FIELDS
    // These are used in the backend of the QD+ algorithm.
    // option_indicator: i32,
    // v_M: f64,
    // v_N: f64,
    // v_h: f64,
    // v_qQD: f64,
    // v_qQDdot: f64,
    // v_p: f64,
    // v_theta: f64,
    // v_c: f64,
    // v_b: f64,
    // v_d1: f64,
    // v_d2: f64,
    // v_dlogSdh: f64,
}

impl CostFunction for QDplus {
    type Param = f64;
    type Output = f64;

    fn cost(&self, p: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
        Ok(self.exercise_boundary_func(*p))
    }
}

// class QDplus:
//     """QD+ alogrithm for computing approximated american option price"""
impl QDplus {
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

    pub fn new(
        riskfree: f64,
        dividend: f64,
        volatility: f64,
        strike: f64,
        option_type: OptionType,
        expiry_date: Date,
    ) -> Self {
        Self {
            r: riskfree,
            q: dividend,
            v: volatility,
            K: strike,
            option_type,
            expiry_date,

            // // # miscellaneous with tau only
            // v_M: 0.0,
            // v_N: 0.0,
            // v_h: 0.0,
            // v_qQD: 0.0,
            // v_qQDdot: 0.0,
            // // # miscellaneous terms with tau and S
            // v_p: 0.0,
            // v_theta: 0.0,
            // v_c: 0.0,
            // v_b: 0.0,
            // v_d1: 0.0,
            // v_d2: 0.0,
            // v_dlogSdh: 0.0,
            // Sb: 0.0,
            tolerance: f64::EPSILON,
        }
    }

    fn year_fraction(&self) -> f64 {
        (OffsetDateTime::now_utc().date() - self.expiry_date).whole_days() as f64 / 365.25
    }

    // def price(self, tau, S):
    //     if tau == 0:
    //         self.exercise_boundary = self.K
    //         return max(S-self.K, 0.0)

    //     self.exercise_boundary = Sb = self.compute_exercise_boundary(tau)
    //     err = self.exercise_boundary_func(Sb, tau)
    //     print("err = ", err)

    //     self.compute_miscellaneous(tau, Sb)
    //     qQD = self.v_qQD
    //     c = self.v_c
    //     b = self.v_b
    //     if self.option_type == OptionType.Put:
    //         pS = europ.EuropeanOption.european_put_value(tau, S, self.r, self.q, self.sigma, self.K)
    //         pSb = europ.EuropeanOption.european_put_value(tau, Sb, self.r, self.q, self.sigma, self.K)
    //     else:
    //         pS = europ.EuropeanOption.european_call_value(tau, S, self.r, self.q, self.sigma, self.K)
    //         pSb = europ.EuropeanOption.european_call_value(tau, Sb, self.r, self.q, self.sigma, self.K)

    //     if self.option_indicator * (Sb - S) <= 0:
    //         return self.option_indicator * (S - self.K)
    //     else:
    //         return pS + (self.K - Sb - pSb)/(1 - b * np.square(np.log(S/Sb)) - c * np.log(S/Sb)) * np.power(S/Sb, qQD)
    fn price(&self, S: f64) -> f64 {
        let T = self.year_fraction();

        if T == 0.0 {
            // self.exercise_boundary = self.K;
            return f64::max(S - self.K, 0.0);
        }

        let Sb = self.compute_exercise_boundary().unwrap();
        println!("Sb = {}", Sb);
        let err = self.exercise_boundary_func(Sb);
        println!("err = {}", err);

        let (N, M, H, qQD, qQDdot, d1, d2, p, theta, dlogSdh, b, c, c0) = self.precompute(S, T);

        let pS = self.P(T, S);
        let pSb = self.P(T, Sb);

        let (pS, pSb) = match self.option_type {
            OptionType::Put => (
                european_put_value(T, S, self.r, self.q, self.v, self.K),
                european_put_value(T, Sb, self.r, self.q, self.v, self.K),
            ),
            OptionType::Call => (
                european_call_value(T, S, self.r, self.q, self.v, self.K),
                european_call_value(T, Sb, self.r, self.q, self.v, self.K),
            ),
        };

        if self.option_type as i32 as f64 * (Sb - S) <= 0.0 {
            return self.option_type as i32 as f64 * (S - self.K);
        } else {
            return pS
                + (self.K - Sb - pSb)
                    / (1.0 - b * f64::powi(f64::ln(S / Sb), 2) - c * f64::ln(S / Sb))
                    * f64::powf(S / Sb, qQD);
        }
    }

    // def compute_exercise_boundary(self, tau):
    //     if tau == 0:
    //         return self.B_at_zero()
    //     # using x0->0 is critical since there are multiple roots for the target function
    //     res = scipy.optimize.root(self.exercise_boundary_func,x0=self.K, args=(tau,))
    //     #if res.success == False:
    //      #   print("succuess? ", res.success, ", ", res.message, ", res = ", res.x)
    //     return res.x[0]
    fn compute_exercise_boundary(&self) -> Result<f64, argmin::core::Error> {
        let tau = self.year_fraction();

        if tau == 0.0 {
            return Ok(self.B_at_zero());
        }

        let initial_guess = self.K;
        let lower_bound = 0.0;
        let upper_bound = 2.0 * self.K;

        let solver = BrentRoot::new(lower_bound, upper_bound, self.tolerance);
        let executor = Executor::new(*self, solver)
            .configure(|state| state.param(initial_guess).max_iters(10));

        let result = executor.run()?;

        Ok(*result.state().get_best_param().unwrap())
    }

    // def B_at_zero(self):
    //     if self.option_type == OptionType.Call:
    //         if self.r <= self.q:
    //             return self.K
    //         else:
    //             return self.r/self.q * self.K
    //     else:
    //         if self.r >= self.q:
    //             return self.K
    //         else:
    //             return self.r/self.q * self.K
    fn B_at_zero(&self) -> f64 {
        match self.option_type {
            OptionType::Call => match self.r <= self.q {
                true => self.K,
                false => self.r / self.q * self.K,
            },
            OptionType::Put => match self.r >= self.q {
                true => self.K,
                false => self.r / self.q * self.K,
            },
        }
    }

    // def compute_miscellaneous(self, tau, S):
    //     #order cannot be changed
    //     self.v_N = self.N()
    //     self.v_M = self.M()
    //     self.v_h = self.h(tau)
    //     self.v_qQD = self.q_QD(tau)
    //     self.v_qQDdot = self.q_QD_dot()
    //     self.v_d1 = europ.EuropeanOption.d1(tau, S, self.r, self.q, self.sigma, self.K)
    //     self.v_d2 = europ.EuropeanOption.d2(tau, S, self.r, self.q, self.sigma, self.K)
    //     if self.option_type == OptionType.Put:
    //         self.v_p = europ.EuropeanOption.european_put_value(tau, S, self.r, self.q, self.sigma, self.K)
    //     else:
    //         self.v_p = europ.EuropeanOption.european_call_value(tau, S, self.r, self.q, self.sigma, self.K)
    //     self.v_theta = europ.EuropeanOption.european_option_theta(tau, S, self.r, self.q, self.sigma, self.K)
    //     self.v_dlogSdh = self.dlogSdh(tau, S)
    //     self.v_c = self.c(tau, S)
    //     self.v_c0 = self.c0(tau, S)
    //     self.v_b = self.b(tau, S)
    // fn compute_miscellaneous(&mut self, tau: f64, S: f64) {
    //     self.v_N = self.N();
    //     self.v_M = self.M();
    //     self.v_h = self.H(tau);
    //     self.v_qQD = self.q_QD(tau);
    //     self.v_qQDdot = self.q_QD_dot();
    //     self.v_d1 = d1(tau, S, self.r, self.q, self.v, self.K);
    //     self.v_d2 = d2(tau, S, self.r, self.q, self.v, self.K);

    //     self.v_p = match self.option_type {
    //         OptionType::Put => european_put_value(tau, S, self.r, self.q, self.v, self.K),
    //         OptionType::Call => european_call_value(tau, S, self.r, self.q, self.v, self.K),
    //     };

    //     self.v_theta = european_option_theta(tau, S, self.r, self.q, self.v, self.K);
    //     self.v_dlogSdh = self.dlogSdh(tau, S);
    //     self.v_c = self.c(tau, S);
    //     self.v_b = self.b(tau, S);
    // }

    // def exercise_boundary_func(self, S, tau):
    //     if tau == 0:
    //         if type(S) is float:
    //             return 0
    //         else:
    //             return np.ones(S.size) * 0
    //     self.compute_miscellaneous(tau, S)
    //     qQD = self.v_qQD
    //     p = self.v_p
    //     c0 = self.v_c0
    //     d1 = self.v_d1
    //     if self.option_type == OptionType.Call:
    //         ans = (1 - np.exp(-self.q * tau) * stats.norm.cdf(d1)) * S - (qQD) * (S - self.K - p)
    //     else:
    //         ans = (1 - np.exp(-self.q * tau) * stats.norm.cdf(-d1)) * S + (qQD) * (self.K - S - p)
    //     return ans
    fn exercise_boundary_func(&self, S: f64) -> f64 {
        let T = self.year_fraction();

        if T == 0.0 {
            return 0.0;
        }

        let (N, M, H, qQD, qQDdot, d1, d2, p, theta, dlogSdh, b, c, c0) = self.precompute(S, T);

        let gaussian = Normal::standard();

        match self.option_type {
            OptionType::Call => {
                (1. - f64::exp(-self.q * T) * gaussian.cdf(d1)) * S - qQD * (S - self.K - p)
            }
            OptionType::Put => {
                (1. - f64::exp(-self.q * T) * gaussian.cdf(-d1)) * S + qQD * (self.K - S - p)
            }
        }
    }

    // def q_QD(self, tau):
    //     N = self.v_N
    //     M = self.v_M
    //     h = self.v_h
    //     if self.option_type == OptionType.Call:
    //         return -0.5*(N-1) + 0.5 * np.sqrt((N-1)*(N-1) + 4 * M/h)
    //     else:
    //         return -0.5*(N-1) - 0.5 * np.sqrt((N-1)*(N-1) + 4 * M/h)
    // fn q_QD(&self, tau: f64) -> f64 {
    //     let N = self.v_N;
    //     let M = self.v_M;
    //     let h = self.v_h;

    //     match self.option_type {
    //         OptionType::Call => -0.5 * (N - 1.) + 0.5 * ((N - 1.) * (N - 1.) + 4. * M / h).sqrt(),
    //         OptionType::Put => -0.5 * (N - 1.) - 0.5 * ((N - 1.) * (N - 1.) + 4. * M / h).sqrt(),
    //     }
    // }

    // def q_QD_dot(self):
    //     N = self.v_N
    //     M = self.v_M
    //     h = self.v_h
    //     return M/(h * h * np.sqrt((N-1)*(N-1) + 4*M/h))
    // fn q_QD_dot(&self) -> f64 {
    //     let N = self.v_N;
    //     let M = self.v_M;
    //     let h = self.v_h;

    //     M / (h * h * ((N - 1.) * (N - 1.) + 4. * M / h).sqrt())
    // }

    // def c0(self, tau, S):
    //     N = self.v_N
    //     M = self.v_M
    //     h = self.v_h
    //     qQD = self.v_qQD
    //     qQDdot = self.v_qQDdot
    //     p = self.v_p
    //     theta = self.v_theta
    //     c = self.v_c
    //     d1 = self.v_d1
    //     d2 = self.v_d2
    //     return - (1-h)*M/(2*qQD + N - 1) * (1/h - (theta*np.exp(self.r * tau))/(self.r*(self.K - S - p)) + qQDdot/(2*qQD+N-1))
    // fn c0(&self, tau: f64, S: f64) -> f64 {
    //     let N = self.v_N;
    //     let M = self.v_M;
    //     let h = self.v_h;
    //     let qQD = self.v_qQD;
    //     let qQDdot = self.v_qQDdot;
    //     let p = self.v_p;
    //     let theta = self.v_theta;
    //     let c = self.v_c;
    //     let d1 = self.v_d1;
    //     let d2 = self.v_d2;

    //     -(1. - h) * M / (2. * qQD + N - 1.)
    //         * (1. / h - (theta * f64::exp(self.r * tau)) / (self.r * (self.K - S - p))
    //             + qQDdot / (2. * qQD + N - 1.))
    // }

    // def c(self, tau, S):
    //     r = self.r
    //     q = self.q
    //     N = self.v_N
    //     M = self.v_M
    //     h = self.v_h
    //     qQD = self.v_qQD
    //     qQDdot = self.v_qQDdot
    //     p = self.v_p
    //     theta = self.v_theta
    //     c = self.v_c
    //     d1 = self.v_d1
    //     d2 = self.v_d2
    //     dlogSdh = self.v_dlogSdh
    //     c0 = self.c0(tau, S)
    //     return c0 - ((1-h)*M)/(2*qQD + N - 1) \
    //         * ((1 - np.exp(-q * tau)*stats.norm.cdf(-d1))/(self.K - S - p) + qQD/S)\
    //         * dlogSdh
    // fn c(&self, tau: f64, S: f64) -> f64 {
    //     let r = self.r;
    //     let q = self.q;
    //     let N = self.v_N;
    //     let M = self.v_M;
    //     let h = self.v_h;
    //     let qQD = self.v_qQD;
    //     let qQDdot = self.v_qQDdot;
    //     let p = self.v_p;
    //     let theta = self.v_theta;
    //     let c = self.v_c;
    //     let d1 = self.v_d1;
    //     let d2 = self.v_d2;
    //     let dlogSdh = self.v_dlogSdh;
    //     let c0 = self.c0(tau, S);

    //     let mut n = Normal::standard();

    //     c0 - ((1. - h) * M) / (2. * qQD + N - 1.)
    //         * ((1. - f64::exp(-q * tau) * n.cdf(-d1)) / (self.K - S - p) + qQD / S)
    //         * dlogSdh
    // }

    // def b(self, tau, S):
    //     N = self.v_N
    //     M = self.v_M
    //     h = self.v_h
    //     qQD = self.v_qQD
    //     qQDdot = self.v_qQDdot
    //     p = self.v_p
    //     theta = self.v_theta
    //     c = self.v_c
    //     d1 = self.v_d1
    //     d2 = self.v_d2
    //     return ((1-h)*M*qQDdot)/(2*(2*qQD + N - 1))
    // fn b(&self, tau: f64, S: f64) -> f64 {
    //     let N = self.v_N;
    //     let M = self.v_M;
    //     let h = self.v_h;
    //     let qQD = self.v_qQD;
    //     let qQDdot = self.v_qQDdot;
    //     let p = self.v_p;
    //     let theta = self.v_theta;
    //     let c = self.v_c;
    //     let d1 = self.v_d1;
    //     let d2 = self.v_d2;

    //     ((1. - h) * M * qQDdot) / (2. * (2. * qQD + N - 1.))
    // }

    // def dlogSdh(self, tau, S):
    //     N = self.v_N
    //     M = self.v_M
    //     h = self.v_h
    //     qQD = self.v_qQD
    //     qQDdot = self.v_qQDdot
    //     p = self.v_p
    //     theta = self.v_theta
    //     c = self.v_c
    //     d1 = self.v_d1
    //     d2 = self.v_d2
    //     r = self.r
    //     q = self.q
    //     dFdh = qQD * theta * np.exp(self.r * tau)/self.r + qQDdot * (self.K - S - p) \
    //         + (S * self.q *np.exp(-self.q*tau) * stats.norm.cdf(-d1))/(r * (1-h)) \
    //         - (S * np.exp(-self.q * tau) * stats.norm.pdf(d1))/(2*r*tau*(1-h))\
    //         * (2*np.log(S/self.K)/(self.sigma * np.sqrt(tau)) - d1)
    //     dFdS = (1 - qQD) * (1 - np.exp(-q * tau) * stats.norm.cdf(-d1)) \
    //             + (np.exp(-q * tau) * stats.norm.pdf(d1))/(self.sigma * np.sqrt(tau))
    //     return -dFdh/dFdS
    // fn dlogSdh(&self, tau: f64, S: f64) -> f64 {
    //     let N = self.v_N;
    //     let M = self.v_M;
    //     let h = self.v_h;
    //     let qQD = self.v_qQD;
    //     let qQDdot = self.v_qQDdot;
    //     let p = self.v_p;
    //     let theta = self.v_theta;
    //     let c = self.v_c;
    //     let d1 = self.v_d1;
    //     let d2 = self.v_d2;
    //     let r = self.r;
    //     let q = self.q;

    //     let mut n = Normal::standard();

    //     let dFdh = qQD * theta * f64::exp(self.r * tau) / self.r
    //         + qQDdot * (self.K - S - p)
    //         + (S * self.q * f64::exp(-self.q * tau) * n.cdf(-d1)) / (r * (1. - h))
    //         - (S * f64::exp(-self.q * tau) * n.pdf(d1)) / (2. * r * tau * (1. - h))
    //             * (2. * f64::ln(S / self.K) / (self.v * f64::sqrt(tau)) - d1);

    //     let dFdS = (1. - qQD) * (1. - f64::exp(-q * tau) * n.cdf(-d1))
    //         + (f64::exp(-q * tau) * n.pdf(d1)) / (self.v * f64::sqrt(tau));

    //     -dFdh / dFdS
    // }

    // def h(self, tau):
    //     return 1 - np.exp(-self.r * tau)
    fn H(&self, T: f64) -> f64 {
        1.0 - (-self.r * T).exp()
    }

    // def M(self):
    //     return 2 * self.r / (self.sigma * self.sigma)
    fn M(&self) -> f64 {
        2.0 * self.r / (self.v * self.v)
    }

    // def N(self):
    //     return 2 * (self.r - self.q) / (self.sigma * self.sigma)
    fn N(&self) -> f64 {
        2.0 * (self.r - self.q) / (self.v * self.v)
    }

    fn P(&self, S: f64, T: f64) -> f64 {
        match self.option_type {
            OptionType::Put => european_put_value(T, S, self.r, self.q, self.v, self.K),
            OptionType::Call => european_call_value(T, S, self.r, self.q, self.v, self.K),
        }
    }

    /// Compute and pack the intermediate values for the QD+ algorithm.
    fn precompute(
        &self,
        S: f64,
        T: f64,
    ) -> (
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
    ) {
        let mut gaussian = Normal::standard();

        let N = self.N();
        let M = self.M();
        let H = self.H(T);

        let r = self.r;
        let q = self.q;

        let qQD = match self.option_type {
            OptionType::Call => -0.5 * (N - 1.) + 0.5 * ((N - 1.) * (N - 1.) + 4. * M / H).sqrt(),
            OptionType::Put => -0.5 * (N - 1.) - 0.5 * ((N - 1.) * (N - 1.) + 4. * M / H).sqrt(),
        };

        let qQDdot = M / (H * H * ((N - 1.) * (N - 1.) + 4. * M / H).sqrt());

        let d1 = d1(T, S, self.r, self.q, self.v, self.K);
        let d2 = d2(T, S, self.r, self.q, self.v, self.K);

        let p = self.P(T, S);

        let theta = european_option_theta(T, S, self.r, self.q, self.v, self.K);

        let dFdh = qQD * theta * f64::exp(self.r * T) / self.r
            + qQDdot * (self.K - S - p)
            + (S * self.q * f64::exp(-self.q * T) * gaussian.cdf(-d1)) / (r * (1. - H))
            - (S * f64::exp(-self.q * T) * gaussian.pdf(d1)) / (2. * r * T * (1. - H))
                * (2. * f64::ln(S / self.K) / (self.v * f64::sqrt(T)) - d1);

        let dFdS = (1. - qQD) * (1. - f64::exp(-q * T) * gaussian.cdf(-d1))
            + (f64::exp(-q * T) * gaussian.pdf(d1)) / (self.v * f64::sqrt(T));

        let dlogSdh = -dFdh / dFdS;

        let c0 = -(1. - H) * M / (2. * qQD + N - 1.)
            * (1. / H - (theta * f64::exp(self.r * T)) / (self.r * (self.K - S - p))
                + qQDdot / (2. * qQD + N - 1.));

        let c = c0
            - ((1. - H) * M) / (2. * qQD + N - 1.)
                * ((1. - f64::exp(-q * T) * gaussian.cdf(-d1)) / (self.K - S - p) + qQD / S)
                * dlogSdh;

        let b = ((1. - H) * M * qQDdot) / (2. * (2. * qQD + N - 1.));

        (N, M, H, qQD, qQDdot, d1, d2, p, theta, dlogSdh, b, c, c0)
    }
}

#[cfg(test)]
mod test_qdplus {
    use super::*;
    use crate::OptionType;
    use time::macros::date;

    #[test]
    fn test_qdplus() {
        let qdp = QDplus {
            r: 0.05,
            q: 0.0,
            v: 0.2,
            K: 100.0,
            option_type: OptionType::Call,
            expiry_date: date!(2022 - 01 - 01),
            tolerance: f64::EPSILON,
        };

        let S = 100.0;
        let price = qdp.price(S);
        println!("price = {}", price);

        // assert!(1 == 0);
    }
}
