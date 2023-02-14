mutable struct Environment
    seed::Int64
    instance::String
    max_generations::Int64
    max_time::Int64
    population_size::Int64
    n::Int64
    m::Int64
    s::Int64
    A::Matrix{Float64}
    encoding::String
    initialization_method::String
    initialization_params::Vector{Any}
    selecion_method::String
    parent_selection_method::String
    selection_params::Vector{Any}
    mutation_method::String
    mutation_probability::Float64
    mutation_params::Vector{Any}
    crossover_method::String
    crossover_probability::Float64
    crossover_params::Vector{Any}
    elite_size::Float64
    offspring_size::Float64
end