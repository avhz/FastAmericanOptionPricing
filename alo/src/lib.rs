pub mod american_pricer_a;
pub mod american_pricer_b;
pub mod american_pricer_core;
pub mod american_pricer_qd_plus;
pub mod chebyshev;
pub mod european_pricer;

#[derive(Debug, Default, Clone, Copy)]
pub enum OptionType {
    #[default]
    Call = 1,
    Put = -1,
}

#[derive(Debug, Default, Clone, Copy)]
pub enum OptionStyle {
    #[default]
    European,
    American,
    Bermudan,
}

#[derive(Debug, Default, Clone, Copy)]
pub enum OptionSettlement {
    #[default]
    Cash,
    Physical,
}
