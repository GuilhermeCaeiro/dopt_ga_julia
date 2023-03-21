using Random
using Plots
using JSON

include("utils.jl")
include("environment.jl")
include("individual.jl")
include("initialization.jl")
include("selection.jl")
include("mutation.jl")
include("crossover.jl")
# include("adaptation.jl")
include("genetic_algorithm.jl")

function read_setup(setup_file)
    experiment_setup = JSON.parsefile(setup_file)
    experiments = []

    _experiments = sort(vec([collect(x) for x in Iterators.product(
        experiment_setup["seed"],
        experiment_setup["instance"],
        experiment_setup["max_generations"],
        experiment_setup["max_time"],
        experiment_setup["population_size"],
        experiment_setup["initialization_method"],
        experiment_setup["selection_method"],
        experiment_setup["parent_selection_method"],
        experiment_setup["mutation_method"],
        experiment_setup["mutation_probability"],
        experiment_setup["crossover_method"],
        experiment_setup["crossover_probability"],
        experiment_setup["elite_size"],
        experiment_setup["offspring_size"],
        experiment_setup["adaptation_method"],
        experiment_setup["generations_until_adaptation"],
        experiment_setup["perform_prlike_crossover"]
    )]))

    for _experiment in _experiments
        experiment = Dict()

        experiment["seed"] = _experiment[1]
        experiment["instance"] = _experiment[2]
        experiment["max_generations"] = _experiment[3]
        experiment["max_time"] = _experiment[4]
        experiment["population_size"] = _experiment[5]
        experiment["initialization_method"] = _experiment[6]
        experiment["selection_method"] = _experiment[7]
        experiment["parent_selection_method"] = _experiment[8]
        experiment["mutation_method"] = _experiment[9]
        experiment["mutation_probability"] = _experiment[10]
        experiment["crossover_method"] = _experiment[11]
        experiment["crossover_probability"] = _experiment[12]
        experiment["elite_size"] = _experiment[13]
        experiment["offspring_size"] = _experiment[14]
        experiment["adaptation_method"] = _experiment[15]
        experiment["generations_until_adaptation"] = _experiment[16]
        experiment["perform_prlike_crossover"] = _experiment[17]
        
        push!(experiments, experiment)
    end

    return experiments
end

function retrieve_additional_params(experiment)
    # seed
    if experiment["seed"] == 0
        experiment["seed"] = rand(1:100000)
    end

    # params
    candidates = ["initialization", "selection", "crossover", "mutation", "adaptation", "parent_selection"]

    #println(experiment)

    for candidate in candidates
        method = experiment[candidate * "_method"][1]
        params = experiment[candidate * "_method"][2]

        experiment[candidate * "_method"] = method
        experiment[candidate * "_params"] = params
    end
end

function write_header(additional_results_fields)
    header_line = ""

    for field in fieldnames(Environment)
        sfield = string(field)
        if sfield == "A"
            continue
        end
        header_line = header_line * string(field) * ";"
    end

    for field in additional_results_fields
        header_line = header_line * field * ";"
    end

    header_line = strip(header_line, ';')

    file = open("results.csv", "w")
    write(file, header_line * "\n")
    close(file)


end

function write_result(environment::Environment, results::Vector{Individual}, total_time::Float64)
    if !isfile("results.csv")
        write_header([
            "total_time_in_seconds",
            "best_cost",
            "solution"
        ])
    end

    result_line = ""

    for field in fieldnames(Environment)
        sfield = string(field)
        #println(field, typeof(field))
        if sfield == "A"
            #println("Ignoring A")
            continue
        end
        result_line = result_line * string(getfield(environment, field)) * ";"
    end

    result_line *= string(total_time) * ";"
    result_line *= string(results[1].fitness) * ";"
    result_line *= string(results[1].chromosome) * ";"

    #println(result_line, "\n")

    result_line = strip(result_line, ';')

    #println(result_line, "\n")

    file = open("results.csv", "a")
    write(file, result_line * "\n")
    close(file)

end

function main()
    experiments = read_setup("experiment_setup.json")

    for experiment in experiments
        retrieve_additional_params(experiment)
        println(experiment)

        A, R, m, n, s = read_instance(experiment["instance"][1], experiment["instance"][2])

        environment = Environment(
            experiment["seed"], # seed
            string(experiment["instance"][1]) * "-" * string(experiment["instance"][2]), # instance
            experiment["max_generations"], # max_generations
            experiment["max_time"], # max_time
            experiment["population_size"], # population_size
            n, # n
            m, # m
            s, # s
            A, # A
            experiment["initialization_method"], # initialization_method
            experiment["initialization_params"], # initialization_params
            experiment["selection_method"], # selecion_method
            experiment["selection_params"], # selection_params
            experiment["parent_selection_method"], # parent_selection_method
            experiment["parent_selection_params"], # parent_selection_params
            experiment["mutation_method"], # mutation_method
            experiment["mutation_probability"], # mutation_probability
            experiment["mutation_params"], # mutation_params
            experiment["crossover_method"], # crossover_method
            experiment["crossover_probability"], # crossover_probability
            experiment["crossover_params"], # crossover_params
            experiment["elite_size"], # elite_size
            experiment["offspring_size"], # offspring_size
            experiment["adaptation_method"], # adaptation_method | accepts "none" and "reset"
            experiment["adaptation_params"], # adaptation_params
            experiment["generations_until_adaptation"], # generations_until_adaptation
            experiment["perform_prlike_crossover"], # perform path-relinking-like crossover
        )

        println(environment.instance)

        start_time = get_time_in_ms()

        Random.seed!(environment.seed)
        
        global execution_statistics = Dict(
            "current_generation" => 0,
            "conventional_of_calc_calls" => Dict(0 => 0),
            "efficient_of_calc_calls" => Dict(0 => 0),
            "parent_zcalc_from_children" => Dict(0 => 0),
        )
        
        ga = GeneticAlgorithm(environment)
        results, solutions = loop(ga)

        total_time = get_time_in_ms() - start_time
        total_time_float = total_time.value/1000.0

        write_result(environment, results, total_time_float)

        #break
    end

    #println(experiments)
    #println([experiment for experiment in experiments])
    #println(length(experiments))

end

main()