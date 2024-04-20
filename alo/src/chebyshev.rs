use std::f64::consts::PI;

pub struct ChebyshevInterpolation {
    n: usize,
    a: Vec<f64>,
    x_min: f64,
    x_max: f64,
    // x_to_cheby: Option<fn(f64, f64, f64) -> f64>,
}

impl ChebyshevInterpolation {
    // def __init__(self, ynodes, x_to_cheby, x_min, x_max):
    //      numpoints = len(ynodes)
    //      self.n = numpoints-1
    //      self.a = [0] * numpoints
    //      self.x_to_cheby = x_to_cheby
    //      self.x_min = x_min
    //      self.x_max = x_max
    //      for k in range(numpoints):
    //          self.a[k] = self.std_coeff(k, ynodes)
    pub fn new(ynodes: Vec<f64>, x_min: f64, x_max: f64) -> Self {
        let n_points = ynodes.len();
        let mut a = Vec::with_capacity(n_points);

        for i in 0..n_points {
            a.push(Self::std_coeff(i, &ynodes));
        }

        Self {
            n: n_points - 1,
            a,
            x_min,
            x_max,
        }
    }

    // def to_cheby_point(self, x, x_min, x_max):
    //      # x in [x_min, x_max] is transformed to [-1, 1]
    //      return np.sqrt(4 * (x - x_min) / (x_max - x_min)) - 1
    pub fn to_cheby_point(x: f64, x_min: f64, x_max: f64) -> f64 {
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
    pub fn value(&self, zv: f64) -> f64 {
        let z = Self::to_cheby_point(zv, self.x_min, self.x_max);
        self.std_cheby_single_value(z)
    }

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
    fn std_cheby_single_value(&self, z: f64) -> f64 {
        let mut b0 = self.a[self.n] * 0.5;
        let mut b1 = 0.0;
        let mut b2 = 0.0;
        for k in (0..self.n).rev() {
            b1 = b0;
            b2 = b1;
            b0 = self.a[k] + 2.0 * z * b1 - b2;
        }
        0.5 * (b0 - b2)
    }

    // def std_coeff(self, k, n_y):
    //      ans = 0
    //      for i in range(0, self.n+1):
    //          term = n_y[i] * np.cos(i * k * np.pi / self.n)
    //          if i == 0 or i == self.n:
    //              term *= 0.5
    //          ans += term
    //      ans *= (2.0/self.n)
    //      return ans
    fn std_coeff(k: usize, n_y: &Vec<f64>) -> f64 {
        let mut ans = 0.0;
        let n = n_y.len();

        for i in 0..(n + 1) {
            let mut term = n_y[i] * f64::cos((i * k) as f64 * PI / n as f64);
            if i == 0 || i == n {
                term *= 0.5;
            }
            ans += term;
        }

        ans *= 2.0 / n as f64;
        ans
    }

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
}
