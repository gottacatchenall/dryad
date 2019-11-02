include("./src/include.jl")

param_dict = Dict(
    "ipc"       => [10],
    "m"         => collect(10^-10:0.0001:10^-2),
    "n_pops"    => [5, 20, 40],
    "k"         => [500, 2000, 4000],
    "ibd_str"   => [3.0, 8.0, 15.0]
)

treatment_df = create_treatments(param_dict)

#===
batch_fits(
        treatment_df,
        num_generations=20000,
        replicates_per_treatment = 100,
        log_frequency = 500,
        base_random_seed = 5,
        metadata_file="more_ibd_metadata.csv",
        data_file="more_ibd.csv",
        dispersal_kernel_type=@ibd_diskern
) ===#

batch_fits_multicore(
    treatment_df,
    2,
    num_generations=20000,
    replicates_per_treatment = 100,
    log_frequency = 500,
    base_random_seed = 5,
    metadata_file="multicore_metadata.csv",
    data_file="multicore.csv",
    dispersal_kernel_type=@ibd_diskern
)
