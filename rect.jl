include("./src/include.jl")

param_dict = Dict(
    "ipc"       => [10],
    "m"         => collect(range(exp(-9), exp(-4), length=100)),
    "n_pops"    => [20],
    "k"         => [2000],
    "ibd_str"   => [3.0, 8.0, 15.0],
    "mutation_rate" => [10^-7]
)

treatment_df = create_treatments(param_dict)

batch_fits(
    treatment_df,
    num_generations=10000,
    replicates_per_treatment = 50,
    log_frequency = 500,
    base_random_seed = 5,
    dispersal_kernel_type=@rect_diskern
)
