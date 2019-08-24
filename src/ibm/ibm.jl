include("../types.jl")

function run_ibm(mp::metapop, metadata::DataFrame; n_gen::Int64=1000, migration_rate::Float64=0.01, selection_strength::Float64=0.1,  b::Float64=3.0, log_freq::Int64=20, genome_file::String="genome.csv", pop_file::String="pops.csv", id::Int64 = 1, init_poly_ct::Float64 = 3.0, n_ef::Int64=1, n_chromo::Int64=5, genome_length=100.0, k::Int64=2000, rseed::Int64=1)

    g::genome = init_random_genome(n_ef=n_ef, n_chromo=n_chromo, genome_length=genome_length, init_poly_ct=init_poly_ct)

    # Run IBM
    ibm_instance::ibm = init_ibm(mp, g, rs=rseed)
    update_ibm_metadata(metadata, id, migration_rate, selection_strength, k, n_ef, n_chromo, genome_length, init_poly_ct)
    run_n_generations(ibm_instance, n_gen, migration_rate, selection_strength, b, log_freq, genome_file, pop_file, id)

end

function run_n_generations(instance::ibm, n_gen::Int64, m::Float64, s::Float64, b::Float64, log_freq::Int64, genome_file::String, pop_file::String, id::Int64)

    n_log::Int64 = round(n_gen/log_freq)+1
    n_pops = length(instance.mp.populations)
    #n_log_pops::Int64 = n_log*n_pops
    #n_log_global::Int64 =

    genomeDF = DataFrame(gen=[], gst=[], jostd=[], mean_poly_ct_per_pop=[])
    popDF = DataFrame(gen=[], pop=[], w_mean=[], prop_of_k=[], mean_poly_ct=[])

    st = @extant

    for g = 0:n_gen
        if g % log_freq == 0
            log_pt::Int64 = floor(g / log_freq) + 1
            st = run_gen(instance, g, m, s, b, genomeDF, popDF, log_pt)
        else
            st = run_gen(instance, g, m, s, b, genomeDF, popDF, 0)
        end

        if (st == @extinct)
            break
        end

    end

    # fill ID columns
    genomeDF.id = fill(id, length(genomeDF.gen))
    popDF.id = fill(id, length(popDF.gen))

    # append to file if it exists

    ap::Bool = isfile(genome_file)
    CSV.write(genome_file, genomeDF, append=ap)
    CSV.write(pop_file, popDF, append=ap)

end

function run_gen(instance::ibm, gen::Int64, mig_rate::Float64, s::Float64, b::Float64, globalDF::DataFrame, popDF::DataFrame, log_pt::Int64)
    dispersal(instance, migration_rate=mig_rate)
    selection(instance, b=b, s=s)

    if log_pt != 0
        logging(instance, gen, globalDF, popDF, log_pt)
    end

    st = reproduction(instance)

    return st
end
