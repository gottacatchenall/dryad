include("./src/include.jl")


ibm_metadata = DataFrame(id=[], m=[], s=[], k=[], n_ef=[], n_chromo=[], genome_length=[], init_poly_ct=[])

mp = init_random_metapop(num_populations=5, diskern_type=@uniform_diskern)
set_mp_total_k(mp, 1000)
run_ibm(mp, ibm_metadata, n_gen=500, log_freq=10, init_poly_ct=5.0, selection_strength=5.0)

CSV.write("ibm_metadata.csv", ibm_metadata)
