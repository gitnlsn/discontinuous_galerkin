pub mod triangle {
    pub mod element;
    pub mod integrands {
        pub mod dirichlet_constraint;
        pub mod flux_artificial;
        pub mod flux_natural;
        pub mod mass;
        pub mod neumann_constraint;
        pub mod penalty_zero;
        pub mod utils;
    }
    pub mod system_builder {
        pub mod domain;
        pub mod internal {
            pub mod mass;
            pub mod flux_natural;
            pub mod flux_artificial;
        }
        pub mod external {
            pub mod dirichlet;
            pub mod boundary_constraint;
        }
    }
}

pub mod common {
    pub mod interfaces;
    pub mod field;
    pub mod point;
    pub mod edge;
}
