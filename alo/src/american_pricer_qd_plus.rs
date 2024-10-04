use crate::european_pricer::*;
use crate::OptionType;
use argmin::core::{CostFunction, Executor, State};
use argmin::solver::brent::BrentRoot;
use statrs::distribution::{Continuous, ContinuousCDF, Normal};
use time::{Date, OffsetDateTime};

#[derive(Debug, Clone, Copy)]
pub struct QDplus {
    pub strike: f64,
    pub riskfree: f64,
    pub dividend: f64,
    pub volatility: f64,
    pub expiration_date: Date,
    pub evaluation_date: Option<Date>,
    pub option_type: OptionType,
}

const DAYS_IN_YEAR: f64 = 365.25;
const MAX_ITERS: u64 = 10_000;
const LOWER_BOUND: f64 = 0.0;
const UPPER_BOUND_MULTIPLIER: f64 = 100.0;

impl CostFunction for QDplus {
    type Param = f64;
    type Output = f64;

    fn cost(&self, p: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
        Ok(self.exercise_boundary_function(*p))
    }
}

impl QDplus {
    pub fn new(
        strike: f64,
        riskfree: f64,
        dividend: f64,
        volatility: f64,
        option_type: OptionType,
        expiration_date: Date,
        evaluation_date: Option<Date>,
    ) -> Self {
        Self {
            riskfree,
            dividend,
            volatility,
            strike,
            option_type,
            expiration_date,
            evaluation_date,
        }
    }

    pub fn price(&self, underlying: f64) -> f64 {
        let t = self.year_fraction();
        let s = underlying;

        if t == 0.0 {
            return f64::max(s - self.strike, 0.0);
        }

        let s_star = self.compute_exercise_boundary().unwrap();

        let lambda = self.lambda(t);
        let b = self.li_2009_eq30(t);
        let c = self.li_2009_eq31(s, t);

        let p_s = self.value(s, t);
        let p_s_star = self.value(s_star, t);

        match self.option_type as i32 as f64 * (s_star - s) <= 0.0 {
            true => return self.option_type as i32 as f64 * (s - self.strike),
            false => {
                return p_s
                    + (self.strike - s_star - p_s_star)
                        / (1.0 - b * f64::powi(f64::ln(s / s_star), 2) - c * f64::ln(s / s_star))
                        * f64::powf(s / s_star, lambda)
            }
        }
    }

    fn year_fraction(&self) -> f64 {
        let today = OffsetDateTime::now_utc().date();

        match self.evaluation_date {
            Some(date) => {
                if date > self.expiration_date {
                    panic!("The evaluation date must be before the expiration date.");
                }

                (self.expiration_date - date).whole_days() as f64 / DAYS_IN_YEAR
            }
            None => {
                if today > self.expiration_date {
                    panic!(
                        "You must provide an evaluation date if the expiration date is in the past."
                    );
                }

                (self.expiration_date - today).whole_days() as f64 / DAYS_IN_YEAR
            }
        }
    }

    fn compute_exercise_boundary(&self) -> Result<f64, argmin::core::Error> {
        let tau = self.year_fraction();

        if tau == 0.0 {
            return Ok(self.exercise_boundary_at_zero());
        }

        let initial_guess = self.strike;
        let lower_bound = LOWER_BOUND;
        let upper_bound = self.strike * UPPER_BOUND_MULTIPLIER;

        let solver = BrentRoot::new(lower_bound, upper_bound, f64::EPSILON);
        let executor = Executor::new(*self, solver)
            .configure(|state| state.param(initial_guess).max_iters(MAX_ITERS));

        let result = executor.run()?;

        Ok(*result.state().get_best_param().unwrap())
    }

    fn exercise_boundary_at_zero(&self) -> f64 {
        match self.option_type {
            OptionType::Call => match self.riskfree <= self.dividend {
                true => self.strike,
                false => self.riskfree / self.dividend * self.strike,
            },
            OptionType::Put => match self.riskfree >= self.dividend {
                true => self.strike,
                false => self.riskfree / self.dividend * self.strike,
            },
        }
    }

    fn exercise_boundary_function(&self, underlying: f64) -> f64 {
        let t = self.year_fraction();
        let s = underlying;

        if t == 0.0 {
            return 0.0;
        }

        let gaussian = Normal::standard();

        let lambda = self.lambda(t);
        let p = self.value(underlying, t);

        let d1 = d1(
            t,
            s,
            self.riskfree,
            self.dividend,
            self.volatility,
            self.strike,
        );

        match self.option_type {
            OptionType::Call => {
                (1. - f64::exp(-self.dividend * t) * gaussian.cdf(d1)) * underlying
                    - lambda * (underlying - self.strike - p)
            }
            OptionType::Put => {
                (1. - f64::exp(-self.dividend * t) * gaussian.cdf(-d1)) * underlying
                    + lambda * (self.strike - underlying - p)
            }
        }
    }

    fn h(&self, year_fraction: f64) -> f64 {
        1.0 - (-self.riskfree * year_fraction).exp()
    }

    fn scale(&self) -> f64 {
        2.0 * self.riskfree / (self.volatility * self.volatility)
    }

    fn omega(&self) -> f64 {
        2.0 * (self.riskfree - self.dividend) / (self.volatility * self.volatility)
    }

    fn lambda(&self, year_fraction: f64) -> f64 {
        let (w, scale, h) = (self.omega(), self.scale(), self.h(year_fraction));

        match self.option_type {
            OptionType::Call => 0.5 * ((1. - w) + ((w - 1.) * (w - 1.) + 4. * scale / h).sqrt()),
            OptionType::Put => 0.5 * ((1. - w) - ((w - 1.) * (w - 1.) + 4. * scale / h).sqrt()),
        }
    }

    fn lambda_prime(&self, t: f64) -> f64 {
        let (w, scale, h) = (self.omega(), self.scale(), self.h(t));

        scale / (h * h * ((w - 1.) * (w - 1.) + 4. * scale / h).sqrt())
    }

    fn value(&self, underlying: f64, year_fraction: f64) -> f64 {
        let (r, q, v, k) = (self.riskfree, self.dividend, self.volatility, self.strike);

        match self.option_type {
            OptionType::Put => european_put_value(year_fraction, underlying, r, q, v, k),
            OptionType::Call => european_call_value(year_fraction, underlying, r, q, v, k),
        }
    }

    fn c_zero(&self, underlying: f64, year_fraction: f64) -> f64 {
        let s = underlying;
        let t = year_fraction;

        let (r, q, v, k) = (self.riskfree, self.dividend, self.volatility, self.strike);

        let (p, w, scale, h, lambda, lambda_prime, theta) = (
            self.value(s, t),
            self.omega(),
            self.scale(),
            self.h(year_fraction),
            self.lambda(t),
            self.lambda_prime(t),
            european_option_theta(year_fraction, s, r, q, v, k),
        );

        let term_1 = -(1. - h) * scale / (2. * lambda + w - 1.);
        let term_2 = 1. / h;
        let term_3 = f64::exp(r * t) * theta / (r * (self.strike - s - p));
        let term_4 = lambda_prime / (2. * lambda + w - 1.);

        term_1 * (term_2 - term_3 + term_4)
    }

    fn li_2009_eq30(&self, year_fraction: f64) -> f64 {
        let t = year_fraction;

        let w = self.omega();
        let h = self.h(t);
        let scale = self.scale();
        let lambda = self.lambda(t);
        let lambda_prime = self.lambda_prime(t);

        ((1. - h) * scale * lambda_prime) / (2. * (2. * lambda + w - 1.))
    }

    fn li_2009_eq31(&self, s: f64, year_fraction: f64) -> f64 {
        let gaussian = Normal::standard();

        let t = year_fraction;

        let (r, q, v, k) = (self.riskfree, self.dividend, self.volatility, self.strike);

        let d1 = d1(t, s, r, q, v, k);
        let p = self.value(s, t);

        let w = self.omega();
        let c0 = self.c_zero(s, t);
        let h = self.h(t);
        let scale = self.scale();
        let lambda = self.lambda(t);

        c0 - ((1. - h) * scale) / (2. * lambda + w - 1.)
            * ((1. - f64::exp(-q * t) * gaussian.cdf(-d1)) / (self.strike - s - p) + lambda / s)
            * self.li_2009_eq36(s, t)
    }

    fn li_2009_eq36(&self, s: f64, t: f64) -> f64 {
        -self.li_2009_eq37(s, t) / self.li_2009_eq38(s, t)
    }

    fn li_2009_eq37(&self, s: f64, t: f64) -> f64 {
        let gaussian = Normal::standard();
        let r = self.riskfree;
        let q = self.dividend;

        let d1 = d1(t, s, r, q, self.volatility, self.strike);

        let lambda = self.lambda(t);
        let lambda_prime = self.lambda_prime(t);
        let theta = european_option_theta(t, s, r, q, self.volatility, self.strike);
        let p = self.value(t, s);
        let h = self.h(t);

        let term_1 = lambda * theta * f64::exp(r * t) / r;
        let term_2 = lambda_prime * (self.strike - s - p);
        let term_3 = (s * q * f64::exp(-q * t) * gaussian.cdf(-d1)) / (r * (1. - h));
        let term_4 = (s * f64::exp(-q * t) * gaussian.pdf(d1)) / (2. * r * t * (1. - h));
        let term_5 = 2. * f64::ln(s / self.strike) / (self.volatility * f64::sqrt(t)) - d1;

        term_1 + term_2 + term_3 - term_4 * term_5
    }

    fn li_2009_eq38(&self, s: f64, t: f64) -> f64 {
        let gaussian = Normal::standard();

        let r = self.riskfree;
        let q = self.dividend;

        let d1 = d1(t, s, r, q, self.volatility, self.strike);
        let lambda = self.lambda(t);

        let term_1 = 1. - lambda;
        let term_2 = 1. - f64::exp(-q * t) * gaussian.cdf(-d1);
        let term_3 = f64::exp(-q * t) * gaussian.pdf(d1) / (self.volatility * f64::sqrt(t));

        term_1 * term_2 + term_3
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
            riskfree: 0.05,
            dividend: 0.02,
            volatility: 0.2,
            strike: 10.0,
            option_type: OptionType::Call,
            expiration_date: date!(2025 - 07 - 20),
            evaluation_date: Some(date!(2024 - 07 - 20)),
        };

        println!("Year fraction = {:?}", qdp.year_fraction());

        println!(
            "ATM = {}, Boundary = {}",
            qdp.price(10.0),
            qdp.compute_exercise_boundary().unwrap()
        );
        println!(
            "ITM = {}, Boundary = {}",
            qdp.price(05.0),
            qdp.compute_exercise_boundary().unwrap()
        );
        println!(
            "OTM = {}, Boundary = {}",
            qdp.price(15.0),
            qdp.compute_exercise_boundary().unwrap()
        );

        let mut qdp = QDplus {
            riskfree: 0.045,
            dividend: 0.05,
            volatility: 0.25,
            strike: 130.0,
            option_type: OptionType::Call,
            expiration_date: date!(2025 - 07 - 20),
            evaluation_date: Some(date!(2024 - 07 - 20)),
        };

        for i in 1..=20 {
            qdp.expiration_date = date!(2025 - 07 - 20) + time::Duration::days(i * 365);

            println!("{}", qdp.compute_exercise_boundary().unwrap());
        }
    }
}
