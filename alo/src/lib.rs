pub mod chebyshev;
pub mod european;
pub mod fixed_point_a;
pub mod fixed_point_b;
pub mod qdplus;

#[derive(Debug, Default, Clone, Copy)]
pub enum TypeFlag {
    #[default]
    Call = 1,
    Put = -1,
}
