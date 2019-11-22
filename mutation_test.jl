include("./src/include.jl")

param_dict = Dict(
    "ipc"       => [10],
    "m"         => collect(10^-5),
    "n_pops"    => [20],
    "k"         => [2000],
    "ibd_str"   => [8.0],
    "mutation_rate" => [10^-5, 10^-4, 10^-3]
)

treatment_df = create_treatments(param_dict)

batch_fits_multicore(
    treatment_df,
    1,
    num_generations=20000,
    replicates_per_treatment = 100,
    log_frequency = 500,
    base_random_seed = 5,
    metadata_file="mutation_test_metadata.csv",
    data_file="mutation_test.csv",
    dispersal_kernel_type=@ibd_diskern
)
