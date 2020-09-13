pub mod triangle {
    pub mod element;
    pub mod boundary_constraint;
    pub mod integrands {
        pub mod dirichlet_constraint;
        pub mod flux_artificial;
        pub mod flux_natural;
        pub mod mass;
        pub mod neumann_constraint;
        pub mod interface_penalty;
        pub mod utils;
    }
    pub mod system_builder {
        pub mod domain;
        pub mod builder;
        pub mod assembler_utils;
        pub mod internal {
            pub mod mass;
            pub mod flux_natural;
            pub mod flux_artificial;
            pub mod jump_penalty;
        }
        pub mod external {
            pub mod dirichlet;
            pub mod neumann;
        }
    }
}

pub mod common {
    pub mod interfaces;
    pub mod field;
    pub mod point;
    pub mod edge;
}
