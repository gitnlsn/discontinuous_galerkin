pub enum DynamicState {
    SteadyState,
    Transient,
}

pub struct DirichletConstraint {
    pub number: usize,
    // value:
}

pub struct NeumannConstraint {
    pub number: usize,
    // value:
}

pub struct ExternalField {
    pub number: usize,
    // value:
}

pub struct HeatProblem {
    pub triangulation: &Triangulation,
    pub conduction: f64,
    pub capacity: f64,
    pub density: f64,
    pub dirichlet_constraints: Vec<DirichletConstraint>,
    pub neumann_constraints: Vec<NeumannConstraint>,
    pub external_field: Vec<ExternalField>,
    pub state: DynamicState,
}
