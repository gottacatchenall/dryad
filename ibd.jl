include("./src/include.jl")

param_dict = Dict(
    "ipc"       => [10],
    "m"         => vcat(collect(0.0001:0.00001:0.003), collect(0.003:0.0001:0.01)),
    "n_pops"    => [20],
    "k"         => [2000],
    "ibd_str"   => [0.5, 1.5, 3.0]
)

treatment_df = create_treatments(param_dict)

batch_fits(
        treatment_df,
        num_generations=20000,
        replicates_per_treatment = 100,
        log_frequency = 500,
        base_random_seed = 5,
        metadata_file="ibd_fits_metadata.csv",
        data_file="ibd_fits.csv",
        dispersal_kernel_type=@ibd_diskern
)
