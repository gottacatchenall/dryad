include("./src/include.jl")

exp_array = collect(range(-9, -3, step=0.5))
exp_array = collect(range(-9, -3, step=0.1))
mig = zeros(Float64, length(exp_array))

for i = 1:length(exp_array)
    mig[i] = exp(exp_array[i])
end

param_dict = Dict(
    "ipc"       => [10],
    "m"         => mig,
    "n_pops"    => [20],
    "k"         => [2000],
    "ibd_str"   => [3.0, 8.0, 15.0],
    "mutation_rate" => [10^-7]
)

treatment_df = create_treatments(param_dict)

batch_fits(
    treatment_df,
    num_generations=10000,
    replicates_per_treatment = 100,
    log_frequency = 500,
    base_random_seed = 5,
    dispersal_kernel_type=@gauss_diskern
)
