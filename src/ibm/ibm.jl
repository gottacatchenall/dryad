include("../types.jl")

mutable struct chromosome
    n_loci::Int64
    length::Float64
    selection_str::Array{Float64}
    ef_map::Array{Int64}
    poly_ct::Array{Int64}
end

mutable struct genome
    n_loci::Int64
    n_haplo::Int64
    mutation_rate::Float64
    allele_dict::Array{Array{Float64,1},1}
    chromosomes::Array{chromosome}
end

mutable struct ibm
    mp::metapop
    g::genome
    dk::dispersal_kernel
    genotypes::Array{Float64, 3}
    population_map::Array{Int64, 1}
    fitness_map::Array{Float64,1}
end

function run_n_generations(instance::ibm, n_gen::Int64; migration_rate::Float64=0.01, s::Float64=0.1, log_freq::Int64=20, genome_file::String="genome.csv", pop_file::String="pops.csv", id::Int64 = 1)

    n_log::Int64 = round(n_gen/log_freq)+1
    n_pops = length(instance.mp.populations)
    #n_log_pops::Int64 = n_log*n_pops
    #n_log_global::Int64 =

    genomeDF = DataFrame(gen=[], gst=[], jostd=[], locus=[])
    popDF = DataFrame(gen=[], pop=[], w_mean=[], prop_of_k=[], mean_poly_ct=[])

    for g = 0:n_gen
        if g % log_freq == 0
            log_pt::Int64 = floor(g / log_freq) + 1
            run_gen(instance, g, migration_rate, s, genomeDF, popDF, log_pt)
        else
            run_gen(instance, g, migration_rate, s, genomeDF, popDF, 0)
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

function run_gen(instance::ibm, gen::Int64, mig_rate::Float64, s::Float64, globalDF::DataFrame, popDF::DataFrame, log_pt::Int64)
    dispersal(instance, migration_rate=mig_rate)
    selection(instance, b=3.0, s=s)

    if log_pt != 0
        logging(instance, gen, globalDF, popDF, log_pt)
    end

    reproduction(instance)
end
