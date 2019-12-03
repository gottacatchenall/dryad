include("./src/include.jl")


ibm_metadata = DataFrame(id=[], m=[], s=[], k=[], n_ef=[], n_chromo=[], genome_length=[], init_poly_ct=[])

mp = init_random_metapop(num_populations=1, diskern_type=@gauss_diskern)
set_mp_total_k(mp, 1000)
run_ibm(mp, ibm_metadata, n_gen_env_shift = 300, n_gen=1000, log_freq=10, init_poly_ct=10.0, selection_strength=1.5, loci_per_ef = 25, migration_rate=0.01, b = 3.0)

CSV.write("ibm_metadata.csv", ibm_metadata)
