use statrs::statistics::Statistics;
use std::f64::consts::PI;

pub struct ChebyshevInterpolation {
    x_values: Vec<f64>,
    y_values: Vec<f64>,
    // x_min: f64,
    // x_max: f64,
    // x_to_cheby: Option<fn(f64, f64, f64) -> f64>,
}

impl ChebyshevInterpolation {
    /// Create a new Chebyshev interpolation object.
    /// Note: The `x_values` need to be sorted in ascending order,
    /// and the `y_values` must be the corresponding
    /// function values at the `x_values`.
    pub fn new(x_values: &[f64], y_values: &[f64]) -> Self {
        assert!(x_values.len() == y_values.len());

        Self {
            x_values: x_values.to_vec(),
            y_values: y_values.to_vec(),
        }
    }

    // def __init__(self, ynodes, x_to_cheby, x_min, x_max):
    //      numpoints = len(ynodes)
    //      self.n = numpoints-1
    //      self.a = [0] * numpoints
    //      self.x_to_cheby = x_to_cheby
    //      self.x_min = x_min
    //      self.x_max = x_max
    //      for k in range(numpoints):
    //          self.a[k] = self.std_coeff(k, ynodes)
    pub fn interpolate(&self, x_values: &[f64]) -> Vec<f64> {
        let coefficients = self.compute_coefficients(&self.y_values);

        x_values
            .iter()
            .map(|x| self.compute_value(*x, &coefficients))
            .collect::<Vec<f64>>()
    }

    /// Perform the Clenshaw algorithm at a single point.
    pub fn clenshaw(&self, z: f64, coefficients: &[f64]) -> f64 {
        let n = self.y_values.len() - 1;

        let mut b0 = coefficients[n] * 0.5;
        let mut b1 = 0.0;
        let mut b2 = 0.0;

        for k in (0..n).rev() {
            (b1, b2) = (b0, b1);
            b0 = coefficients[k] + 2.0 * z * b1 - b2;
        }

        0.5 * (b0 - b2)
    }

    /// Compute a single coefficient of the Chebyshev polynomial.
    pub fn compute_coefficient(&self, k: usize, q: &[f64]) -> f64 {
        let n = q.len();

        let term_1 = (1. / 2. * n as f64) * (q[0] + (-1_f64).powi(n as i32 - 1) * q[n - 1]);

        let term_2 = (2. / n as f64)
            * q.iter()
                .enumerate()
                .map(|(i, q)| q * f64::cos((i * k) as f64 * PI / n as f64))
                .sum::<f64>();

        term_1 + term_2
    }

    /// Compute multiple coefficients of the Chebyshev polynomial.
    pub fn compute_coefficients(&self, y_values: &[f64]) -> Vec<f64> {
        (0..y_values.len())
            .into_iter()
            .map(|i| self.compute_coefficient(i, &self.y_values))
            .collect::<Vec<f64>>()
    }

    /// Compute the Chebyshev polynomial at a single point.
    pub fn compute_value(&self, x: f64, coefficients: &[f64]) -> f64 {
        let (x_min, x_max) = (self.x_values.clone().min(), self.x_values.clone().max());

        self.clenshaw(Self::to_chebyshev_point(x, x_min, x_max), coefficients)
    }

    // def to_cheby_point(self, x, x_min, x_max):
    //      # x in [x_min, x_max] is transformed to [-1, 1]
    //      return np.sqrt(4 * (x - x_min) / (x_max - x_min)) - 1
    pub fn to_chebyshev_point(x: f64, x_min: f64, x_max: f64) -> f64 {
        (4.0 * (x - x_min) / (x_max - x_min)).sqrt() - 1.0
    }

    // def value(self, zv):
    //      ans = []
    //      if zv is float:
    //          to_cheby = [zv]
    //      else:
    //          to_cheby = zv
    //      if self.x_to_cheby is not None:
    //          to_cheby = self.x_to_cheby(zv, self.x_min, self.x_max)
    //      for z in to_cheby:
    //          ans.append(self.std_cheby_single_value(z))
    //      return ans
    // pub fn value(&self, x: f64) -> f64 {
    //     self.std_cheby_single_value(Self::to_cheby_point(
    //         x,
    //         self.x_values.clone().min(),
    //         self.x_values.clone().max(),
    //     ))
    // }

    // def std_cheby_single_value(self, z):
    // """z is the point to be valued between [-1, 1], n_y are the function values at Chebyshev points
    // Iteration using Clenshaw algorithm"""
    // b0 = self.a[self.n] * 0.5
    // b1 = 0
    // b2 = 0
    // for k in range(self.n - 1, -1, -1):
    //     b1, b2 = b0, b1
    //     b0 = self.a[k] + 2 * z * b1 - b2
    // return 0.5 * (b0 - b2)
    // fn std_cheby_single_value(&self, z: f64) -> f64 {
    //     let n = self.a.len() - 1;
    //     let mut b0 = self.a[n] * 0.5;
    //     let mut b1 = 0.0;
    //     let mut b2 = 0.0;
    //     for k in (0..n).rev() {
    //         (b1, b2) = (b0, b1);
    //         b0 = self.a[k] + 2.0 * z * b1 - b2;
    //     }
    //     0.5 * (b0 - b2)
    // }

    // def std_coeff(self, k, n_y):
    //      ans = 0
    //      for i in range(0, self.n+1):
    //          term = n_y[i] * np.cos(i * k * np.pi / self.n)
    //          if i == 0 or i == self.n:
    //              term *= 0.5
    //          ans += term
    //      ans *= (2.0/self.n)
    //      return ans
    // fn std_coeff(k: usize, n_y: &[f64]) -> f64 {
    //     let mut result = 0.0;
    //     let n = n_y.len() - 1;
    //     for i in 0..=n {
    //         let mut term = n_y[i] * f64::cos((i * k) as f64 * PI / n as f64);
    //         if i == 0 || i == n {
    //             term *= 0.5;
    //         }
    //         result += term;
    //     }
    //     result *= 2.0 / n as f64;
    //     result
    // }

    // @staticmethod
    // def get_std_cheby_points(numpoints):
    //     i = np.arange(0, numpoints+1)
    //     return np.cos(i * np.pi/numpoints)
    pub fn get_std_cheby_points(n_points: usize) -> Vec<f64> {
        (0..n_points + 1)
            .into_iter()
            .map(|i| f64::cos(i as f64 * PI / n_points as f64))
            .collect()
    }

    pub fn chebyshev_nodes(n: usize, x_max: f64) -> Vec<f64> {
        (0..n)
            .into_iter()
            .map(|i| {
                let z = -f64::cos(PI * i as f64) / (n as f64);
                (1. + z) * f64::sqrt(x_max) / 2.
            })
            .collect()
    }

    /// Generate Chebyshev nodes of the first kind in the open interval (-1, 1).
    pub fn chebyshev_nodes_first_kind(n: usize) -> Vec<f64> {
        (0..n)
            .into_iter()
            .map(|k| f64::cos(PI * (2. * k as f64 + 1.) / (2. * n as f64)))
            .collect()
    }

    /// Generate Chebyshev nodes of the second kind in the closed interval [-1, 1].
    pub fn chebyshev_nodes_second_kind(n: usize) -> Vec<f64> {
        (0..n)
            .into_iter()
            .map(|k| f64::cos(PI * k as f64 / (n as f64 - 1.)))
            .collect()
    }

    /// Transform Chebyshev nodes to the open interval (a, b).
    pub fn affine_transformation(nodes: &[f64], a: f64, b: f64) -> Vec<f64> {
        let t = |x| 0.5 * ((a + b) + (b - a) * x);

        nodes.into_iter().map(|x| t(x)).collect()
    }

    /// Transform regular points to Chebyshev nodes (reverse affine transformation).
    pub fn affine_transformation_inverse(points: &[f64], a: f64, b: f64) -> Vec<f64> {
        let t = |x| (2. * x - b - a) / (b - a);

        points.into_iter().map(|x| t(x)).collect()
    }

    // // def cheby_to_interval(cheb_points, a, b):
    // //     return a + (b - a) * (cheb_points + 1) * 0.5
    // pub fn chebyshev_to_interval(cheb_points: &[f64], a: f64, b: f64) -> Vec<f64> {
    //     cheb_points
    //         .into_iter()
    //         .map(|point| a + (b - a) * (point + 1.0) * 0.5)
    //         .collect()
    // }

    // // def interval_to_cheby(x, a, b):
    // //     return (2*x - b - a) / (b - a)
    // pub fn interval_to_chebyshev(x: f64, a: f64, b: f64) -> f64 {
    //     (2.0 * x - b - a) / (b - a)
    // }
}

#[cfg(test)]
mod test_chebyshev {
    use std::vec;

    use super::*;
    use rand::distributions::Distribution;
    use rand::{thread_rng, Rng};

    fn trial_function(x: f64) -> f64 {
        x * f64::exp(x)
    }

    fn linspace(start: f64, stop: f64, n: usize) -> Vec<f64> {
        let mut result = Vec::with_capacity(n);
        let step = (stop - start) / (n - 1) as f64;

        for i in 0..n {
            result.push(start + i as f64 * step);
        }

        result
    }

    fn perturb(y: &[f64]) -> Vec<f64> {
        let mut rng = thread_rng();

        y.iter()
            .map(|&x| x + (rng.gen_range(0.0..=1.0) * 2.0 - 1.0) * 0.4 * x)
            .collect::<Vec<f64>>()
    }

    #[test]
    fn test_chebyshev() {
        let n = 15;
        let (a, b) = (0.0, 3.0);

        // let x = ChebyshevInterpolation::get_std_cheby_points(n);
        // let x = ChebyshevInterpolation::chebyshev_nodes_first_kind(n);
        // let x = ChebyshevInterpolation::chebyshev_nodes_second_kind(n);
        // let x = ChebyshevInterpolation::chebyshev_nodes(n);

        // let x = ChebyshevInterpolation::affine_transformation(&x, a, b);

        let x = vec![
            0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0,
        ];
        let y = x.iter().map(|&x| trial_function(x)).collect::<Vec<f64>>();
        let y_perturbed = perturb(&y);

        let mut chebyshev = ChebyshevInterpolation::new(&x, &y_perturbed);

        let xi = linspace(0.0, 1.0, 20);

        let yintrp = chebyshev.interpolate(&xi);

        println!("x = {:?}", x);
        println!("y = {:?}", y);
        println!("y_perturbed = {:?}", y_perturbed);
        println!("xi = {:?}", xi);
        println!("yintrp = {:?}", yintrp);
    }
}
